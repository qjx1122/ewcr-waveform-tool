function yq = svdcs_lagrange_interp(y, tq, order)
%SVDCS_LAGRANGE_INTERP Local Lagrange polynomial interpolation.
%
% y     : column/row vector sampled at integer MATLAB sample positions 1..N
% tq    : arbitrary query positions in the same 1-based sample coordinate
% order : polynomial order; order=3 uses four neighboring samples
%
% This is a direct, toolbox-free implementation. The paper specifies
% Lagrange interpolation and notes a real-time filter-bank implementation,
% but does not publish the exact filter-bank coefficients.

y = double(y(:));
tq = double(tq(:));
N = numel(y);

m = order + 1;
if m < 2
    error('order must be >= 1');
end
if N < m
    error('Signal is too short for the requested interpolation order.');
end

yq = zeros(size(tq));

for ii = 1:numel(tq)
    t = tq(ii);

    % Center a local stencil around t.
    left = floor(t) - floor((m-1)/2);
    left = max(1, min(left, N-m+1));
    idx = left:(left+m-1);

    nodes = double(idx(:));
    vals = y(idx);

    % Exact sample short-circuit.
    [dmin, jj] = min(abs(nodes - t));
    if dmin < 10*eps(max(1,abs(t)))
        yq(ii) = vals(jj);
        continue;
    end

    L = ones(m,1);
    for a = 1:m
        for b = 1:m
            if a ~= b
                L(a) = L(a) * (t - nodes(b)) / (nodes(a) - nodes(b));
            end
        end
    end
    yq(ii) = sum(vals .* L);
end
end
