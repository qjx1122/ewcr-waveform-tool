function [packed, nCodes, codes] = lzw_encode12(bytes)
%LZW_ENCODE12 Lossless LZW over uint8 bytes with fixed 12-bit codewords.
%
% Dictionary:
%   initial codes 0..255 = literal bytes
%   new codes 256..4095
%   when full, dictionary simply stops growing (no reset code)
%
% packed is a uint8 byte stream containing 12-bit codes.
% nCodes is required by the decoder because the last packed nibble may pad.

bytes = uint8(bytes(:).');
if isempty(bytes)
    packed = uint8([]);
    nCodes = 0;
    codes = uint16([]);
    return;
end

MAXC = 4095;
% transition(code+1, byte+1) gives a longer-sequence code, zero if absent.
transition = zeros(MAXC+1,256,'uint16');
nextCode = uint16(256);

codes = zeros(1,max(16,ceil(numel(bytes)/2)),'uint16');
nc=0;

w = uint16(bytes(1)); % literal code 0..255

for i=2:numel(bytes)
    b = uint16(bytes(i));
    z = transition(double(w)+1,double(b)+1);

    if z ~= 0
        w = z;
    else
        nc=nc+1;
        if nc > numel(codes)
            codes(end+max(1024,ceil(numel(codes)/2))) = uint16(0);
        end
        codes(nc)=w;

        if nextCode <= MAXC
            transition(double(w)+1,double(b)+1)=nextCode;
            nextCode = nextCode + 1;
        end
        w=b;
    end
end

nc=nc+1;
if nc > numel(codes), codes(end+1)=uint16(0); end
codes(nc)=w;
codes=codes(1:nc);
nCodes=nc;
packed=pack12_codes(codes);
end
