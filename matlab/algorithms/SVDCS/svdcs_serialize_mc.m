function bytes = svdcs_serialize_mc(mc)
%SVDCS_SERIALIZE_MC Deterministic byte serialization of the 3-row MC matrix.
%
% Paper specifies MC logically as [q; SUM; CNT] and then applies LZW, but it
% does not specify a byte-level file format. This project uses:
%   uint16 q | int64 SUM | uint32 CNT = 14 bytes / MC column
% in little-endian byte order.
%
% MATLAB compatibility note:
%   computer('endian') is not a valid R2025b call. Endianness is obtained
%   from the third output of computer:
%       [~,~,endian] = computer;
%
% This packing choice affects the numerical CR after LZW but does NOT alter
% the spectral-variation algorithm or reconstruction.

n = numel(mc.q);
bytes = zeros(1,14*n,'uint8');
p = 1;

[~,~,endian] = computer;
needSwap = (endian == 'B');

for k = 1:n
    b = to_le_uint16(mc.q(k),needSwap);
    bytes(p:p+1)=b; p=p+2;

    b = to_le_int64(mc.sum(k),needSwap);
    bytes(p:p+7)=b; p=p+8;

    b = to_le_uint32(mc.cnt(k),needSwap);
    bytes(p:p+3)=b; p=p+4;
end
end

function b = to_le_uint16(v,needSwap)
v = uint16(v);
if needSwap, v = swapbytes(v); end
b = typecast(v,'uint8');
end

function b = to_le_uint32(v,needSwap)
v = uint32(v);
if needSwap, v = swapbytes(v); end
b = typecast(v,'uint8');
end

function b = to_le_int64(v,needSwap)
v = int64(v);
if needSwap, v = swapbytes(v); end
b = typecast(v,'uint8');
end
