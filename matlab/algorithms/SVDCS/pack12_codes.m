function packed = pack12_codes(codes)
%PACK12_CODES Pack uint16 values 0..4095 into an MSB-first 12-bit stream.

codes=uint16(codes(:).');
if any(codes>4095)
    error('12-bit code exceeds 4095.');
end

n=numel(codes);
packed=zeros(1,ceil(12*n/8),'uint8');
p=1;
i=1;

while i+1<=n
    a=uint16(codes(i));
    b=uint16(codes(i+1));

    packed(p)   = uint8(bitshift(a,-4));
    packed(p+1) = uint8(bitor(bitshift(bitand(a,15),4), bitshift(b,-8)));
    packed(p+2) = uint8(bitand(b,255));
    p=p+3; i=i+2;
end

if i<=n
    a=uint16(codes(i));
    packed(p)   = uint8(bitshift(a,-4));
    packed(p+1) = uint8(bitshift(bitand(a,15),4));
end
end
