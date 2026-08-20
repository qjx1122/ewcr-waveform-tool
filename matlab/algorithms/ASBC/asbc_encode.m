function [stream, info] = asbc_encode(x, fs, opts)
%ASBC_ENCODE  Adaptive subband compression (ASBC) encoder.
%
%   Implements the adaptive subband compression scheme of
%   X. Wang, Y. Liu, L. Tong, "Adaptive Subband Compression for Streaming
%   of Continuous Point-on-Wave and PMU Data", IEEE Trans. Power Systems,
%   vol.36, no.6, pp.5612-5621, 2021.
%
%   Key ideas of the paper (Sec. II-IV):
%     * signal model  x(t) = sum_k ak(t) cos(k*Omega0*t + phi_k(t)) + e(t)
%     * each harmonic k is down-shifted to baseband, low-pass filtered with
%       bandwidth W_k/2, and decimated by S_k = ceil(fs/W_k)   (Eqs. 4,6,7)
%     * block-wise activity detection decides which subbands are encoded
%     * the interharmonic residual e(t) is encoded by an FFT-(k,L) sparse
%       coder when active
%   This implementation realizes the ideal filterbank in the FFT domain;
%   band-limited interpolation (interpft) is used at the decoder. It is a
%   faithful engineering realization of the published architecture.
%
%   [stream, info] = ASBC_ENCODE(x, fs)             uses defaults
%   [stream, info] = ASBC_ENCODE(x, fs, opts)
%     opts.W            harmonic subband bandwidth in Hz   (default 10)
%     opts.Rk           bits/sample baseband quantization  (default 8)
%     opts.block_cycles block length in fundamental cycles (default 5)
%     opts.act_thresh   relative activity threshold        (default 1e-4)
%     opts.ke_max       max FFT coefficients for interharmonic (default 48)
%     opts.energy_frac  retained energy fraction for interharmonic (0.99)
%     opts.RQ           native ADC bits/sample of raw stream (default 16)
%
%   stream: struct with all encoded quantities (see asbc_decode.m).
%   info  : diagnostics (K, subband rates, bit budget, CR...).

arguments
    x (:,1) double {mustBeFinite}
    fs (1,1) double {mustBePositive}
    opts (1,1) struct = struct()
end

% merge user options over defaults
D = struct('W', 10, 'Rk', 8, 'block_cycles', 5, 'act_thresh', 1e-4, ...
           'ke_max', 48, 'energy_frac', 0.99, 'RQ', 16);
fn = fieldnames(D);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = D.(fn{i}); end
end

N    = numel(x);
RQ   = opts.RQ;

% ---- fundamental frequency estimation (peak of spectrum in 45..55 Hz) ---
fft_len = 2^nextpow2(N);
Xf = abs(fft(x, fft_len)).';          % row vector
freqs = (0:fft_len-1)/fft_len*fs;
band = freqs >= 45 & freqs <= 55;
[~, im] = max(Xf .* band);
F0  = freqs(im);
w0  = 2*pi*F0/fs;

% ---- block partition -----------------------------------------------------
blk = round(opts.block_cycles * fs / F0);     % samples per block
nb  = floor(N / blk);
if nb < 2, blk = N; nb = 1; end

% ---- harmonic subband selection K (energy-based, vectorized) ------------
W    = opts.W;
Kmax = floor((fs/2 - W/2)/F0);
bin_spacing = fs/fft_len;
kcenters = round((1:Kmax)*F0/bin_spacing) + 1;
kcenters(kcenters > fft_len) = [];
Ek = zeros(1, Kmax);
Ek(1:numel(kcenters)) = Xf(kcenters).^2;
E1 = Ek(1);
thr = opts.act_thresh * E1;
K = 1;
kk_sel = find(Ek(2:end) > thr);
if ~isempty(kk_sel), K = 1 + max(kk_sel); end

% ---- normalization + quantization scale ----------------------------------
% normalize so max|x| == 1; baseband magnitudes |yk| <= sum(ak)/2 <= 1 for
% the test signals used here; a fixed uniform quantizer of range [-1,1]
% with Rk bits is applied to every active subband (no per-block scales).
xnorm = x / max(abs(x));
scale = max(abs(x));
qstep = 2 / 2^opts.Rk;
Sk = ceil(fs / W);                            % decimation factor (Eq. 6)

% ---- per-block subband processing ---------------------------------------
mask_h = zeros(nb, K);                        % harmonic activity mask
mask_e = zeros(nb, 1);                        % interharmonic activity mask
sub_data = cell(nb, K);
sub_meta = cell(nb, K);
e_fft    = cell(nb, 1);
nn = (0:blk-1)';

keep = (0:blk-1)';
keep = keep <= (W/2)/fs*blk | keep >= blk - (W/2)/fs*blk;  % |f| < W/2

for b = 1:nb
    seg = (b-1)*blk + (1:blk);
    xb  = xnorm(seg);

    % fundamental subband (always active)
    yb = xb .* exp(-1i*w0*nn);
    Yb = fft(yb);
    Yb(~keep) = 0;
    ylp = ifft(Yb);
    yds = ylp(1:Sk:end);
    sub_data{b,1} = quant_uni(yds, opts.Rk, qstep);
    sub_meta{b,1} = numel(yds);
    mask_h(b,1) = 1;

    % higher harmonics: block-wise activity detection (energy detector)
    for k = 2:K
        yb = xb .* exp(-1i*k*w0*nn);
        Yb = fft(yb);
        Yb(~keep) = 0;
        ylp = ifft(Yb);
        if sum(abs(ylp).^2) > opts.act_thresh * sum(abs(xb).^2)
            yds = ylp(1:Sk:end);
            sub_data{b,k} = quant_uni(yds, opts.Rk, qstep);
            sub_meta{b,k} = numel(yds);
            mask_h(b,k) = 1;
        end
    end

    % ---- interharmonic residual  e = x - sum_k sqrt(2)Re(yk e^{jkw0n})
    xrec = zeros(blk,1);
    for k = 1:K
        if mask_h(b,k) == 1
            yds = dequant_uni(sub_data{b,k}, opts.Rk, qstep);
            yup = interp_band(yds, blk);
            xrec = xrec + sqrt(2)*real(yup .* exp(1i*k*w0*nn));
        end
    end
    eb = xb - xrec;
    if sum(eb.^2) > opts.act_thresh * sum(xb.^2)
        mask_e(b) = 1;
        Eb = fft(eb);
        half = floor(blk/2);
        amp = abs(Eb(2:half+1));
        [~, ord] = sort(amp, 'descend');
        cumE = cumsum(amp(ord).^2);
        kk = find(cumE >= opts.energy_frac*sum(amp.^2), 1);
        if isempty(kk), kk = numel(ord); end
        kk = min(kk, opts.ke_max);
        sel = sort(ord(1:kk));
        qe = max(amp) / 2^15;                 % per-block coefficient step
        if qe == 0, qe = 1; end
        e_fft{b}.idx = sel(:);
        e_fft{b}.Re  = round(real(Eb(sel+1))/qe);
        e_fft{b}.Im  = round(imag(Eb(sel+1))/qe);
        e_fft{b}.qe  = qe;
    end
end

% ---- bit accounting (Eqs. 9-14 of the paper) ----------------------------
bits_sub  = 0;
for b = 1:nb
    for k = 1:K
        if mask_h(b,k)
            bits_sub = bits_sub + 2*sub_meta{b,k}*opts.Rk;
        end
    end
end
bits_mask = nb*K + nb;                        % harmonic + interharmonic masks
idx_bits  = max(1, ceil(log2(blk/2)));
bits_ie   = 0;
for b = 1:nb
    if mask_e(b)
        kk = numel(e_fft{b}.idx);
        bits_ie = bits_ie + 32 + kk*(idx_bits + 2*32);   % qe + idx + Re + Im
    end
end
bits_header = 128;                            % codec parameters
bits_scale  = 64;                             % global scale (float64)
bits_total  = bits_header + bits_scale + bits_mask + bits_sub + bits_ie;

info.F0        = F0;
info.W         = W;
info.Sk        = Sk;
info.blk       = blk;
info.nb        = nb;
info.K         = K;
info.RQ_bits   = RQ;
info.bits_original   = N*RQ;
info.bits_compressed = bits_total;
info.CR        = (N*RQ)/bits_total;
info.bits_sub  = bits_sub;
info.bits_interharmonic = bits_ie;
info.bits_mask = bits_mask;
info.active_frac_h = mean(mask_h(:,2:end), 'all');
info.active_frac_e = mean(mask_e);
info.mask_h = mask_h;

% ---- assemble stream -----------------------------------------------------
stream.F0    = F0;
stream.W     = W;
stream.Rk    = opts.Rk;
stream.RQ    = RQ;
stream.fs    = fs;
stream.blk   = blk;
stream.nb    = nb;
stream.Sk    = Sk;
stream.K     = K;
stream.qstep = qstep;
stream.scale = scale;
stream.mask_h = mask_h;
stream.mask_e = mask_e;
stream.sub_data = sub_data;
stream.sub_meta = sub_meta;
stream.e_fft    = e_fft;
stream.N        = N;
stream.bits_compressed = bits_total;
stream.bits_original   = N*RQ;
end

% ======================================================================
function cq = quant_uni(v, R, qstep)
% uniform midtread quantizer of a complex vector, range [-1,1] per axis
    cr = min(2^R-1, max(0, round((real(v)+1)/qstep)));
    ci = min(2^R-1, max(0, round((imag(v)+1)/qstep)));
    cq = cr + 1i*ci;
end

function v = dequant_uni(c, R, qstep) %#ok<INUSD>
    v = real(c)*qstep - 1 + 1i*(imag(c)*qstep - 1);
end

function y = interp_band(v, L)
% band-limited interpolation of a decimated baseband sequence to length L
% (ideal reconstruction by zero-padding the spectrum / sinc interpolation)
    v = v(:);
    if numel(v) == 1
        y = repmat(v, L, 1);       % a single phasor -> constant signal
    else
        y = interpft(v, L);
        y = y(:);
    end
end
