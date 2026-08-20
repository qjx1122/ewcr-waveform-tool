function codes = unpack12_codes(packed,nCodes)
%UNPACK12_CODES Inverse of pack12_codes.

packed=uint8(packed(:).');
codes=zeros(1,nCodes,'uint16');
p=1; i=1;

while i+1<=nCodes
    b1=uint16(packed(p));
    b2=uint16(packed(p+1));
    b3=uint16(packed(p+2));

    codes(i)   = bitor(bitshift(b1,4), bitshift(b2,-4));
    codes(i+1) = bitor(bitshift(bitand(b2,15),8), b3);
    p=p+3; i=i+2;
end

if i<=nCodes
    b1=uint16(packed(p));
    b2=uint16(packed(p+1));
    codes(i)=bitor(bitshift(b1,4),bitshift(b2,-4));
end
end
