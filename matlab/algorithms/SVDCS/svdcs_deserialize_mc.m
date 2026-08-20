function mc = svdcs_deserialize_mc(bytes, Q, Nframes)
%SVDCS_DESERIALIZE_MC Inverse of svdcs_serialize_mc.
%
% MATLAB compatibility note:
%   R2025b endianness is queried via [~,~,endian] = computer.

bytes = uint8(bytes(:).');
if mod(numel(bytes),14) ~= 0
    error('Serialized MC length must be a multiple of 14 bytes.');
end
n = numel(bytes)/14;

q   = zeros(1,n,'uint16');
sumv= zeros(1,n,'int64');
cnt = zeros(1,n,'uint32');

[~,~,endian] = computer;
needSwap = (endian == 'B');

p=1;
for k=1:n
    q(k)=from_le_uint16(bytes(p:p+1),needSwap); p=p+2;
    sumv(k)=from_le_int64(bytes(p:p+7),needSwap); p=p+8;
    cnt(k)=from_le_uint32(bytes(p:p+3),needSwap); p=p+4;
end

mc.q=q; mc.sum=sumv; mc.cnt=cnt;
mc.Q=Q; mc.Nframes=Nframes;
end

function v=from_le_uint16(b,needSwap)
v=typecast(uint8(b),'uint16');
if needSwap, v=swapbytes(v); end
end

function v=from_le_uint32(b,needSwap)
v=typecast(uint8(b),'uint32');
if needSwap, v=swapbytes(v); end
end

function v=from_le_int64(b,needSwap)
v=typecast(uint8(b),'int64');
if needSwap, v=swapbytes(v); end
end
