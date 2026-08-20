function bytes = lzw_decode12(packed, nCodes)
%LZW_DECODE12 Inverse of lzw_encode12.

if nCodes==0
    bytes=uint8([]);
    return;
end

codes = unpack12_codes(packed,nCodes);
MAXC = 4095;

dict = cell(MAXC+1,1);
for k=0:255
    dict{k+1}=uint8(k);
end
nextCode=256;

oldCode = double(codes(1));
if oldCode > 255
    error('Invalid first LZW code.');
end
prev = dict{oldCode+1};

% Dynamic output buffer.
cap=max(1024,4*nCodes);
bytes=zeros(1,cap,'uint8');
nout=0;
[nout,bytes]=append_bytes(nout,bytes,prev);

for i=2:numel(codes)
    c=double(codes(i));

    if c < nextCode && ~isempty(dict{c+1})
        entry=dict{c+1};
    elseif c == nextCode
        entry=[prev, prev(1)];
    else
        error('Invalid LZW code %d at position %d.',c,i);
    end

    [nout,bytes]=append_bytes(nout,bytes,entry);

    if nextCode <= MAXC
        dict{nextCode+1}=[prev, entry(1)];
        nextCode=nextCode+1;
    end
    prev=entry;
end

bytes=bytes(1:nout);
end

function [nout,out]=append_bytes(nout,out,v)
m=numel(v);
if nout+m > numel(out)
    out(end+max(1024,2*m+ceil(numel(out)/2)))=uint8(0);
end
out(nout+1:nout+m)=v;
nout=nout+m;
end
