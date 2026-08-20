function b = var_bytes(x)
%VAR_BYTES MATLAB in-memory size of one variable, including nested content.
% This is MATLAB representation size, NOT compressed bitstream size and
% NOT target firmware RAM usage.
tmp = x; %#ok<NASGU>
w = whos('tmp');
b = w.bytes;
end
