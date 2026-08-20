function y = mmc_reconstruct_model(m, theta, state, cfg)
%MMC_RECONSTRUCT_MODEL Reconstruct normalized waveform from a model packet.
N=cfg.N; t=(0:N-1)'/cfg.fs;
family=m.family;
if strcmp(family,'pred_para') && isfield(m,'base_family')
    family=m.base_family;
end
switch family
    case 'none'
        y=zeros(N,1);
    case 'sin'
        y=theta(1)*cos(2*pi*theta(2)*t+theta(3));
    case 'poly'
        tc=linspace(-1,1-1/N,N)';
        B=local_cheby(tc,m.order);
        y=B*theta(:);
    case 'pred_samples'
        xprev=[state.prev2;state.prev1]*2^(-m.kx);
        X=local_pred(xprev,N,m.order,m.eta);
        y=X*theta(:);
    otherwise
        error('Unknown model family: %s',family);
end
end

function B=local_cheby(t,ord)
B=zeros(numel(t),ord+1); B(:,1)=1;
if ord>=1,B(:,2)=t;end
for k=2:ord,B(:,k+1)=2*t.*B(:,k)-B(:,k-1);end
end

function X=local_pred(xprev,N,ord,eta)
X=zeros(N,ord);
for i=1:N
    base=N+i-eta;
    for j=1:ord
        idx=base-j+1;
        if idx>=1 && idx<=numel(xprev),X(i,j)=xprev(idx);end
    end
end
end
