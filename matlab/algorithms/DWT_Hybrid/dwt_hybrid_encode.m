function [stream, info] = dwt_hybrid_encode(x, opts)
%DWT_HYBRID_ENCODE  Multi-stage hybrid DWT compression encoder.
%
%   Follows the modern "multi-stage hybrid" synchro-waveform compression
%   family (e.g. "Lossy Compression: An Online Multi-Stage Technology for
%   High-Fidelity Synchro-Waveform Measurements", IEEE Trans. Industry
%   Applications, 2025):  DWT -> adaptive thresholding -> uniform
%   quantization -> delta coding -> entropy coding.
%
%   Stages implemented here:
%     1. L-level discrete wavelet transform (wavedec, default db4);
%     2. level-adaptive soft thresholding:
%          - 'universal': T_j = sigma_j*sqrt(2*ln(N_j)), sigma_j = MAD/0.6745
%          - 'energy'   : keep the smallest number of coefficients of each
%                         level whose energy reaches thr_frac of the level
%     3. per-level uniform quantization to qbits (signed);
%     4. delta coding of the (slowly varying) approximation coefficients;
%     5. Huffman entropy coding (huffmandict/huffmanenco) of the quantized
%        values and of the run-length gaps between nonzero coefficients.
%
%   [stream, info] = DWT_HYBRID_ENCODE(x)
%     opts.level     decomposition level (default: auto)
%     opts.wname     wavelet name (default 'db4')
%     opts.thr_mode  'universal' | 'energy' (default 'universal')
%     opts.thr_frac  retained energy fraction for 'energy' mode (0.999)
%     opts.qbits     quantization bits per retained coefficient (8)
%     opts.delta     apply delta coding to approximation coefficients (true)
%     opts.RQ        native bits/sample of the raw stream (16)

if nargin<2 || isempty(opts), opts=struct(); end
x=x(:);
if any(~isfinite(x)), error('Input x must be finite.'); end

% merge user options over defaults
D = struct('level', [], 'wname', 'db4', 'thr_mode', 'universal', ...
           'thr_frac', 0.999, 'qbits', 8, 'delta', true, 'RQ', 16);
fn = fieldnames(D);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = D.(fn{i}); end
end

N = numel(x);
if isempty(opts.level)
    L = min(10, max(1, floor(log2(N/16))));
else
    L = min(opts.level, max(1, floor(log2(N/16))));
end

[C, Book] = wavedec(x, L, opts.wname);
% Book = [cA_L, cD_L, ..., cD_1, N];  the last entry is the signal length
nA  = Book(1);
lev_len = Book(2:end-1);  % lengths of detail levels L..1 (this order)
nD = numel(C) - nA;

% ---- thresholding (per detail level) ------------------------------------
Cq  = zeros(size(C));
scales = zeros(L+1,1);        % per-block quantization step
bstart = 1;                   % start of current block in C

% approximation block
blk = C(1:nA);
scales(1) = quant_step(blk, opts.qbits);
Cq(1:nA) = quantize_block(blk, scales(1), opts.qbits);

pos = nA + 1;
for j = L:-1:1
    bj = C(pos:pos+lev_len(L-j+1)-1);    if strcmpi(opts.thr_mode, 'universal')
        sigma = median(abs(bj))/0.6745;
        Tj = sigma*sqrt(2*log(numel(bj)));
        bj = soft_thresh(bj, Tj);
    else  % 'energy'
        bj = energy_retain(bj, opts.thr_frac);
    end
    scales(L-j+2) = quant_step(bj, opts.qbits);
    Cq(pos:pos+numel(bj)-1) = quantize_block(bj, scales(L-j+2), opts.qbits);
    pos = pos + numel(bj);
end

% ---- delta coding of approximation coefficients --------------------------
if opts.delta
    Aq = Cq(1:nA);
    dA = [Aq(1); diff(Aq)];
else
    dA = Cq(1:nA);
end

% ---- entropy coding ------------------------------------------------------
% values of detail levels (already quantized integers)
allvals = Cq(nA+1:end);
[vq_huff, bits_vals] = huff_encode(allvals);

% run-length gaps between nonzero coefficients of the full detail stream
nz = find(allvals ~= 0);
if isempty(nz)
    gaps = [numel(allvals)];
else
    gaps = [nz(1)-1; diff(nz)-1; numel(allvals)-nz(end)];
end
[gaps_huff, bits_gaps] = huff_encode(gaps(:));
nz_count = numel(nz);          % needed by the decoder

% approximation (or delta) coefficients
[a_huff, bits_a] = huff_encode(dA(:));

% ---- bit accounting ------------------------------------------------------
bits_data   = bits_vals + bits_gaps + bits_a;
bits_header = 128 + (L+1)*32 + 64;      % params + per-level scales
bits_total  = bits_data + bits_header;
bits_orig   = N*opts.RQ;

stream.N        = N;
stream.L        = L;
stream.wname    = opts.wname;
stream.thr_mode = opts.thr_mode;
stream.qbits    = opts.qbits;
stream.delta    = opts.delta;
stream.Book     = Book;
stream.scales   = scales;
stream.dA       = dA;
stream.a_dict   = a_huff.dict;
stream.a_bits   = a_huff.bits;
stream.vq_dict  = vq_huff.dict;
stream.vq_bits  = vq_huff.bits;
stream.gap_dict = gaps_huff.dict;
stream.gap_bits = gaps_huff.bits;
stream.nz_count = nz_count;
stream.RQ       = opts.RQ;
stream.bits_compressed = bits_total;
stream.bits_original   = bits_orig;

info.L       = L;
info.bits_original   = bits_orig;
info.bits_compressed = bits_total;
info.CR      = bits_orig/bits_total;
info.bits_data = bits_data;
info.bits_header = bits_header;
info.nA      = nA;
info.nz_frac = nnz(Cq(nA+1:end))/nD;
end

% ======================================================================
function bj = soft_thresh(bj, T)
    bj = sign(bj).*max(abs(bj)-T, 0);
end

function bj = energy_retain(bj, frac)
    [absb, ord] = sort(abs(bj), 'descend');
    cum = cumsum(absb.^2);
    if cum(end) == 0, return; end
    k = find(cum >= frac*cum(end), 1);
    if isempty(k), k = numel(bj); end
    keep = false(size(bj)); keep(ord(1:k)) = true;
    bj(~keep) = 0;
end

function qs = quant_step(blk, qbits)
    mx = max(abs(blk));
    if mx == 0, qs = 1; else, qs = mx/2^(qbits-1); end
end

function vq = quantize_block(blk, qs, qbits)
    vq = round(blk/qs);
    vq = max(-2^(qbits-1), min(2^(qbits-1)-1, vq));
    vq = double(vq);
end

function [out, nbits] = huff_encode(sig)
% Huffman coding of an integer symbol stream; returns dict + coded bits.
    sig = sig(:);
    if isempty(sig)
        out.dict = {}; out.bits = []; nbits = 0; return;
    end
    syms = unique(sig);
    counts = histcounts(sig, [syms; syms(end)+1]);
    prob = counts/sum(counts);
    if numel(syms) < 2
        % degenerate case: single symbol -> fixed 1-bit code
        dict = {syms, {0}};
        enco = zeros(numel(sig), 1);
    else
        dict = huffmandict(syms, prob);
        enco = huffmanenco(sig, dict);
    end
    nbits = numel(enco);
    % dictionary size estimate (symbol value 16 bits + code length)
    nbits = nbits + numel(syms)*48;
    out.dict = dict; out.bits = enco;
end
