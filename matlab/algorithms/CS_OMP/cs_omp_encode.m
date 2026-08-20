function [stream, info] = cs_omp_encode(x, opts)
%CS_OMP_ENCODE  Compressed sensing encoder with OMP reconstruction.
%
%   Implements the compressed-sensing framework for power-system waveform
%   compression/reconstruction (see e.g. "Reconstruction of Harmonic and
%   Transient Electrical Signals Through Compressed Sensing Technique",
%   IEEE Access, 2024, and the dictionary-learning/OMP pipeline of the CN
%   patent CN115827577A by Liu et al., Tianjin Univ., 2023).
%
%   Encoder:  y = Phi * x        (Phi: MxN random Gaussian, M << N)
%   Decoder:  s_hat = argmin ||s||_0  s.t. ||y - Phi*Psi*s||_2 <= eps
%             via orthogonal matching pursuit (OMP), x_hat = Psi*s_hat.
%
%   [stream, info] = CS_OMP_ENCODE(x)
%     opts.M_ratio   M/N measurement ratio (default 0.2)
%     opts.basis     'dft' | 'dct'  (default 'dft')
%     opts.K0        max OMP iterations / assumed sparsity (default 40)
%     opts.tol       relative residual tolerance (default 1e-6)
%     opts.seed      RNG seed for Phi (default 42)
%     opts.ybits     bits per quantized measurement (default 16)
%     opts.RQ        native bits/sample of the raw stream (16)

if nargin<2 || isempty(opts), opts=struct(); end
x=x(:);
if any(~isfinite(x)), error('Input x must be finite.'); end

% merge user options over defaults
D = struct('M_ratio', 0.2, 'basis', 'dft', 'K0', 40, 'tol', 1e-6, ...
           'seed', 42, 'ybits', 16, 'RQ', 16);
fn = fieldnames(D);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = D.(fn{i}); end
end

N = numel(x);
M = max(1, round(opts.M_ratio*N));

rng(opts.seed);
Phi = randn(M, N)/sqrt(M);          % measurement matrix (Gaussian)

switch lower(opts.basis)
    case 'dft'
        Psi = dftmtx(N)/sqrt(N);    % unitary DFT basis
    case 'dct'
        Psi = dctmtx(N);            % orthonormal DCT-II basis
    otherwise
        error('cs_omp:badBasis', 'Unsupported basis "%s".', opts.basis);
end

A  = Phi * Psi;                     % MxN sensing matrix of the sparse signal
y  = Phi * x;                       % measurements
s0 = Psi' * x;                      % true sparse coefficients (reference)

% ---- OMP recovery --------------------------------------------------------
s_hat = omp_recover(A, y, opts.K0, opts.tol);
x_hat = real(Psi * s_hat);

% ---- measurement quantization --------------------------------------------
yscale = max(abs(y));
if yscale == 0, yscale = 1; end
yq = round(y/yscale * 2^(opts.ybits-1));
yq = max(-2^(opts.ybits-1), min(2^(opts.ybits-1)-1, yq));

% ---- bit accounting ------------------------------------------------------
% Quantized measurements at ybits each, plus scale (64 bits), Phi seed
% (32 bits) and codec parameters (64 bits). RQ is the native depth.
bits_meas   = M*opts.ybits;
bits_header = 160;
bits_total  = bits_meas + bits_header;
bits_orig   = N*opts.RQ;

stream.N      = N;
stream.M      = M;
stream.basis  = lower(opts.basis);
stream.seed   = opts.seed;
stream.K0     = opts.K0;
stream.tol    = opts.tol;
stream.ybits  = opts.ybits;
stream.yscale = yscale;
stream.yq     = yq;
stream.RQ     = opts.RQ;
stream.bits_compressed = bits_total;
stream.bits_original   = bits_orig;

info.bits_original   = bits_orig;
info.bits_compressed = bits_total;
info.CR      = bits_orig/bits_total;
info.M       = M;
info.N       = N;
info.sparsity_true  = nnz(abs(s0) > 1e-6*max(abs(s0)));
info.sparsity_rec   = nnz(abs(s_hat) > 1e-6*max(abs(s_hat)));
info.omp_iters      = numel(find(abs(s_hat) > 1e-9));
info.x_hat          = x_hat;        % returned for convenience
end

% ======================================================================
function s_hat = omp_recover(A, y, K0, tol)
% Orthogonal matching pursuit.
%   Solve min ||s||_0 s.t. ||y - A*s||_2 <= tol*||y||_2.
    [M, N] = size(A);
    r  = y;
    supp = zeros(K0, 1);
    As   = zeros(M, K0);            % selected columns of A
    k = 0;
    target = tol*norm(y);
    while k < K0 && norm(r) > target
        k = k + 1;
        c = A' * r;                 % correlation
        [~, idx] = max(abs(c));
        % remove already-selected indices
        for t = 1:k-1
            if idx == supp(t)
                c(idx) = 0;
                [~, idx] = max(abs(c));
            end
        end
        supp(k) = idx;
        As(:,k) = A(:,idx);
        % least-squares on support
        coef = As(:,1:k) \ y;
        r = y - As(:,1:k)*coef;
    end
    s_hat = zeros(N,1);
    if k > 0
        coef = As(:,1:k) \ y;
        s_hat(supp(1:k)) = coef;
    end
end
