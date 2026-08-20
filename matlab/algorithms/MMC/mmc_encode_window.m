function [best_packet, best_xhat, state_next] = mmc_encode_window(x, state, cfg)
%MMC_ENCODE_WINDOW Rate-constrained two-stage MMC for one N-sample window.

N = cfg.N;
x = x(:);
if numel(x) ~= N, error('Window length must equal cfg.N.'); end

% Authors normalize with x_n=x*2^(-kx), kx=ceil(log2(max(abs(x))+eps)).
kx = ceil(log2(max(abs(x)) + 1e-8));
kx = max(0, min(2^cfg.n_kx-1, kx));
xn = x * 2^(-kx);

models = mmc_build_models(xn, state, cfg, kx);
M = numel(models); L = numel(cfg.residual_methods);
nm = max(1, ceil(log2(M)));
nl = max(1, ceil(log2(L)));

best_mse = inf;
best_packet = [];
best_xhat = zeros(N,1);

for im = 1:M
    m = models(im);
    p = numel(m.theta_hat);
    if p == 0
        nx_list = 0;
    else
        nmin = cfg.min_bits_theta*p;
        nmax = min(cfg.max_bits_theta*p, cfg.n_tot);
        nx_list = nmin:cfg.nx_step:nmax;
        if isempty(nx_list), continue; end
    end

    for nx = nx_list
        n_nx = 0;
        if p > 0
            n_nx = max(1,ceil(log2(cfg.max_bits_theta*p + 1)));
        end
        header_stage1 = nm + cfg.n_kx + n_nx;
        if header_stage1 + nx + nl > cfg.n_tot, continue; end

        [theta_q, theta_idx, theta_bits, x_model] = mmc_quantize_model(m, nx, state, cfg);
        r = xn - x_model;

        for il = 1:L
            rname = cfg.residual_methods{il};
            if strcmpi(rname,'none')
                nkr = 0; nnr = 0;
            else
                nkr = cfg.n_kr; nnr = cfg.n_nr;
            end
            nr_budget = cfg.n_tot - (header_stage1 + nx + nl + nkr + nnr);
            if nr_budget < 0, continue; end

            residual = mmc_encode_residual(r, rname, nr_budget, cfg);
            xn_hat = x_model + residual.rhat;
            mse = mean((xn-xn_hat).^2);
            total_bits = header_stage1 + nx + nl + nkr + nnr + residual.nbits;

            if mse < best_mse
                best_mse = mse;
                model_packet = m;
                model_packet.theta_hat = []; % not part of decoder payload
                model_packet.theta_center = m.theta_center;
                model_packet.theta_width = m.theta_width;
                model_packet.theta_idx = theta_idx;
                model_packet.theta_bits = theta_bits;

                best_packet = struct();
                best_packet.model = model_packet;
                best_packet.theta_q = theta_q;
                best_packet.kx = kx;
                best_packet.nx = nx;
                best_packet.nr = residual.nbits;
                best_packet.residual = residual;
                best_packet.header_bits = header_stage1 + nl + nkr + nnr;
                best_packet.total_bits = total_bits;
                best_packet.normalized_mse = mse;
                best_packet.nm = nm; best_packet.nl = nl;
                best_xhat = xn_hat * 2^(kx);
            end
        end
    end
end

if isempty(best_packet)
    error('No feasible MMC configuration for n_tot=%d bits.', cfg.n_tot);
end

state_next = state;
state_next.prev2 = state.prev1;
state_next.prev1 = best_xhat;
state_next.prev_model = best_packet.model;
state_next.prev_model.theta_q = best_packet.theta_q;
if isfield(state,'history_count')
    state_next.history_count = state.history_count + 1;
else
    state_next.history_count = 1;
end
end

function models = mmc_build_models(xn, state, cfg, kx)
% Construct and fit all candidate models for current normalized window.
models = struct('name',{},'family',{},'order',{},'eta',{},'factor',{}, ...
    'theta_hat',{},'theta_center',{},'theta_width',{},'kx',{},'base_family',{});

if cfg.use_none
    models(end+1) = mk('none','none',0,0,0,[],[],[],kx); %#ok<AGROW>
end

if cfg.use_sinusoid
    c = cfg.sin_center; w = cfg.sin_width;
    th = fit_sinusoid(xn,c,w,cfg);
    models(end+1) = mk('sin-1','sin',0,0,0,th,c,w,kx); %#ok<AGROW>
end

for ord = cfg.poly_orders
    [th,c,w] = fit_poly(xn,ord,cfg);
    models(end+1) = mk(sprintf('poly-%d',ord),'poly',ord,0,0,th,c,w,kx); %#ok<AGROW>
end

% Previous reconstructed windows, normalized by current kx.
% A sample-predictive model is undefined before any reconstructed history is
% available. Skipping it on the first window also prevents the all-zero
% predictor matrix from generating rank-deficiency warnings.
xprev = [state.prev2; state.prev1] * 2^(-kx);
hasHistory = isfield(state,'history_count') && state.history_count >= 1;
if hasHistory
    for eta = cfg.sample_pred_eta
        for ord = cfg.sample_pred_orders
            [th,c,w] = fit_sample_pred(xn,xprev,ord,eta);
            models(end+1) = mk(sprintf('samp.-%d-%d',ord,eta),'pred_samples',ord,eta,0,th,c,w,kx); %#ok<AGROW>
        end
    end
end

% Parameter-predictive models: centered on previously selected model parameters.
if ~isempty(state.prev_model) && isfield(state.prev_model,'theta_q') && ~isempty(state.prev_model.theta_q)
    pm = state.prev_model;
    if any(strcmp(pm.family,{'sin','poly','pred_samples'}))
        for factor = cfg.param_pred_factors
            c = pm.theta_q(:).';
            if isfield(pm,'theta_width') && numel(pm.theta_width)==numel(c)
                w = pm.theta_width(:).'/factor;
            else
                w = 0.2*ones(size(c))/factor;
            end
            switch pm.family
                case 'sin'
                    th = fit_sinusoid(xn,c,w,cfg);
                case 'poly'
                    th = fit_poly_bounded(xn,pm.order,c,w,cfg);
                case 'pred_samples'
                    th = fit_sample_pred_bounded(xn,xprev,pm.order,pm.eta,c,w);
            end
            models(end+1) = mk(sprintf('para.-%d-%s',factor,pm.name),'pred_para',pm.order,pm.eta,factor,th,c,w,kx); %#ok<AGROW>
            models(end).base_family = pm.family;
        end
    end
end
end

function s = mk(name,family,order,eta,factor,theta,center,width,kx)
s = struct('name',name,'family',family,'order',order,'eta',eta,'factor',factor, ...
    'theta_hat',theta(:).','theta_center',center(:).','theta_width',width(:).','kx',kx, ...
    'base_family','');
end

function theta = fit_sinusoid(y,center,width,cfg)
% Robust bounded frequency grid + linear amplitude/phase fit.
t = (0:cfg.N-1)'/cfg.fs;
flo = center(2)-width(2)/2; fhi = center(2)+width(2)/2;
if fhi <= flo, fgrid = center(2); else, fgrid = linspace(flo,fhi,81); end
best = inf; theta = center;
for f = fgrid
    A = [cos(2*pi*f*t), sin(2*pi*f*t)];
    ab = A\y;
    amp = hypot(ab(1),ab(2));
    phi = atan2(-ab(2),ab(1)); % a*cos(wt+phi)
    th = [amp f phi];
    lo = center-width/2; hi=center+width/2;
    th = min(max(th,lo),hi);
    yh = th(1)*cos(2*pi*th(2)*t+th(3));
    e = sum((y-yh).^2);
    if e<best, best=e; theta=th; end
end
end

function [theta,c,w] = fit_poly(y,ord,cfg)
c = zeros(1,ord+1); w = cfg.poly_width{ord+1};
tc = linspace(-1,1-1/cfg.N,cfg.N)';
B = cheby_basis(tc,ord);
theta = (B\y).';
theta = min(max(theta,c-w/2),c+w/2);
end

function theta = fit_poly_bounded(y,ord,c,w,cfg)
tc = linspace(-1,1-1/cfg.N,cfg.N)';
B = cheby_basis(tc,ord);
theta = (B\y).';
theta = min(max(theta,c-w/2),c+w/2);
end

function [theta,c,w] = fit_sample_pred(y,xprev,ord,eta)
c = zeros(1,ord); w = 2*ones(1,ord);
X = predictor_matrix(xprev,numel(y),ord,eta);
theta = safe_small_ls(X,y).';
theta = min(max(theta,c-w/2),c+w/2);
end

function theta = fit_sample_pred_bounded(y,xprev,ord,eta,c,w)
X = predictor_matrix(xprev,numel(y),ord,eta);
theta = safe_small_ls(X,y).';
theta = min(max(theta,c-w/2),c+w/2);
end

function theta = safe_small_ls(X,y)
% Robust least-squares for the 1--2 parameter sample predictor.
% For an all-zero or rank-deficient predictor matrix, MATLAB's X\y emits
% repeated rank-deficiency warnings. The minimum-norm pseudoinverse is the
% correct finite least-squares fallback and is inexpensive for <=2 columns.
if isempty(X) || norm(X,'fro') <= 1e-14*max(1,norm(y))
    theta=zeros(size(X,2),1);
    return;
end
G=X.'*X;
scale=max(1,norm(G,'fro'));
if rcond(G) < 1e-12/scale
    theta=pinv(X)*y;
else
    theta=X\y;
end
end

function X = predictor_matrix(xprev,N,ord,eta)
% xprev contains two previous reconstructed windows (2N samples).
X = zeros(N,ord);
for i=1:N
    base = N + i - eta;
    for j=1:ord
        idx = base-j+1;
        if idx>=1 && idx<=numel(xprev), X(i,j)=xprev(idx); end
    end
end
end

function B = cheby_basis(t,ord)
B = zeros(numel(t),ord+1); B(:,1)=1;
if ord>=1, B(:,2)=t; end
for k=2:ord
    B(:,k+1)=2*t.*B(:,k)-B(:,k-1);
end
end

function [theta_q, idx, bits, xhat] = mmc_quantize_model(m,nx,state,cfg)
p = numel(m.theta_hat);
if p==0
    theta_q=[]; idx=[]; bits=[]; xhat=zeros(cfg.N,1); return;
end
bits = allocate_bits_greedy(m,nx,state,cfg);
[theta_q,idx] = quantize_vector(m.theta_hat,bits,m.theta_center,m.theta_width);
xhat = mmc_reconstruct_model(m,theta_q,state,cfg);
end

function bits = allocate_bits_greedy(m,nx,state,cfg)
% Greedy rate-distortion allocation. Starts at minimum bits/parameter and
% adds one bit where it most reduces model-domain MSE.
p=numel(m.theta_hat);
base = cfg.min_bits_theta*ones(1,p);
if sum(base)>nx, base=zeros(1,p); end
bits=base;
while sum(bits)<nx
    bestGain=-inf; bestj=1;
    [q0,~]=quantize_vector(m.theta_hat,bits,m.theta_center,m.theta_width);
    y0=mmc_reconstruct_model(m,q0,state,cfg);
    e0=sum((reconstruct_unquantized(m,state,cfg)-y0).^2);
    for j=1:p
        if bits(j)>=cfg.max_bits_theta, continue; end
        bt=bits; bt(j)=bt(j)+1;
        [q,~]=quantize_vector(m.theta_hat,bt,m.theta_center,m.theta_width);
        yh=mmc_reconstruct_model(m,q,state,cfg);
        gain=e0-sum((reconstruct_unquantized(m,state,cfg)-yh).^2);
        if gain>bestGain, bestGain=gain; bestj=j; end
    end
    if all(bits>=cfg.max_bits_theta), break; end
    bits(bestj)=bits(bestj)+1;
end
end

function y = reconstruct_unquantized(m,state,cfg)
y=mmc_reconstruct_model(m,m.theta_hat,state,cfg);
end

function [q,idx] = quantize_vector(theta,bits,center,width)
p=numel(theta); q=zeros(1,p); idx=zeros(1,p);
for j=1:p
    b=bits(j); c=center(j); w=width(j);
    if b<=0 || w<=0
        idx(j)=0; q(j)=c;
    else
        delta=w/(2^b);
        ind=floor((theta(j)-c)/delta);
        ind=min(ind,2^(b-1)-1); ind=max(ind,-2^(b-1));
        idx(j)=ind; q(j)=delta*(ind+0.5)+c;
    end
end
end
