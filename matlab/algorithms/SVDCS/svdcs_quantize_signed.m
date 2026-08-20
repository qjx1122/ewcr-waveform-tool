function [q, step] = svdcs_quantize_signed(x, nbits, fullScale)
%SVDCS_QUANTIZE_SIGNED Uniform signed quantizer with explicit full scale.
%
% x is assumed to be in physical or p.u. units. fullScale is the positive
% full-scale magnitude. q is stored as int32 for safe arithmetic.

maxCode = 2^(nbits-1) - 1;
minCode = -2^(nbits-1);
step = fullScale / maxCode;

q = round(double(x) ./ step);
q(q > maxCode) = maxCode;
q(q < minCode) = minCode;
q = int32(q);
end
