function xhat = cs_omp_decode(stream)
%CS_OMP_DECODE  Decoder of the compressed-sensing codec.
%   xhat = CS_OMP_DECODE(stream) re-runs OMP on the stored measurements
%   with the (seeded, reproducible) sensing matrix and returns the
%   reconstructed waveform.

s = stream;
rng(s.seed);
Phi = randn(s.M, s.N)/sqrt(s.M);
switch s.basis
    case 'dft'
        Psi = dftmtx(s.N)/sqrt(s.N);
    case 'dct'
        Psi = dctmtx(s.N);
end
A = Phi * Psi;
y = double(s.yq) * s.yscale / 2^(s.ybits-1);   % dequantize measurements
s_hat = omp_recover(A, y, s.K0, s.tol);
xhat = real(Psi * s_hat);
end

% ======================================================================
function s_hat = omp_recover(A, y, K0, tol)
    [M, N] = size(A);
    r  = y;
    supp = zeros(K0, 1);
    As   = zeros(M, K0);
    k = 0;
    target = tol*norm(y);
    while k < K0 && norm(r) > target
        k = k + 1;
        c = A' * r;
        [~, idx] = max(abs(c));
        for t = 1:k-1
            if idx == supp(t)
                c(idx) = 0;
                [~, idx] = max(abs(c));
            end
        end
        supp(k) = idx;
        As(:,k) = A(:,idx);
        coef = As(:,1:k) \ y;
        r = y - As(:,1:k)*coef;
    end
    s_hat = zeros(N,1);
    if k > 0
        coef = As(:,1:k) \ y;
        s_hat(supp(1:k)) = coef;
    end
end
