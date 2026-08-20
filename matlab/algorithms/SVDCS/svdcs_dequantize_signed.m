function x = svdcs_dequantize_signed(q, step)
%SVDCS_DEQUANTIZE_SIGNED Inverse of svdcs_quantize_signed.
x = double(q) .* step;
end
