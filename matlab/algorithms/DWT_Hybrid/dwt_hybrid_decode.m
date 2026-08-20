function xhat = dwt_hybrid_decode(stream)
%DWT_HYBRID_DECODE  Decoder of the multi-stage hybrid DWT codec.
%   xhat = DWT_HYBRID_DECODE(stream) reconstructs the waveform encoded by
%   dwt_hybrid_encode.m (Huffman decode -> delta/approx rebuild -> inverse
%   quantization -> inverse DWT).

s = stream;
N = s.N;
nA = s.Book(1);
M = sum(s.Book(1:end-1)); % total number of wavelet coefficients

% ---- rebuild the approximation (delta) coefficients ----------------------
dA = huff_decode(s.a_dict, s.a_bits, s.delta, nA);
if s.delta
    Aq = cumsum(dA);
else
    Aq = dA;
end

% ---- rebuild the detail coefficients from values + run-length gaps -------
allvals = huff_decode(s.vq_dict, s.vq_bits, false, M - nA);
gaps    = huff_decode(s.gap_dict, s.gap_bits, false, s.nz_count + 1);
Cq = zeros(M,1);
Cq(1:nA) = Aq;
% place nonzero values according to the run-length gaps
vals = allvals(allvals ~= 0);
pos = nA + 1; gi = 1;
for i = 1:numel(vals)
    if gi > numel(gaps), break; end
    pos = pos + gaps(gi); gi = gi + 1;
    if pos > N, break; end
    Cq(pos) = vals(i);
    pos = pos + 1;
end

% ---- inverse quantization + inverse DWT ----------------------------------
C = zeros(size(Cq));
C(1:nA) = Cq(1:nA)*s.scales(1);
pos = nA + 1;
lev_len = s.Book(2:end-1);    % detail lengths, level L..1 order
for j = 1:numel(lev_len)
    seg = pos:pos+lev_len(j)-1;
    C(seg) = Cq(seg)*s.scales(j+1);
    pos = pos + lev_len(j);
end

xhat = waverec(C, s.Book, s.wname);
xhat = xhat(:);
end

% ======================================================================
function out = huff_decode(dict, bits, ~, nOut)
    if isempty(bits)
        out = zeros(nOut,1); return;
    end
    if size(dict,1) < 2
        % single-symbol codec
        out = repmat(dict{1,1}, nOut, 1);
    else
        out = huffmandeco(bits, dict);
    end
    out = out(:);
    if numel(out) < nOut
        out = [out; zeros(nOut-numel(out),1)];
    end
    out = out(1:nOut);
end


