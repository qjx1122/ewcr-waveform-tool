function xhat = asbc_decode(stream)
%ASBC_DECODE  Decoder of the adaptive subband compression (ASBC) codec.
%
%   xhat = ASBC_DECODE(stream) reconstructs the waveform from the stream
%   produced by asbc_encode.m. See asbc_encode.m for the reference and the
%   architecture description.

s = stream;
N = s.N;
blk = s.blk;
nb = s.nb;
w0 = 2*pi*s.F0/s.fs;
nn = (0:blk-1)';

xhat = zeros(N,1);

for b = 1:nb
    seg = (b-1)*blk + (1:blk);
    xrec = zeros(blk,1);

    % harmonic subbands: dequantize -> interpolate -> re-modulate
    for k = 1:s.K
        if s.mask_h(b,k) == 1
            yds = dequant_uni(s.sub_data{b,k}, s.Rk, s.qstep);
            yup = interp_band(yds, blk);
            xrec = xrec + sqrt(2)*real(yup .* exp(1i*k*w0*nn));
        end
    end

    % interharmonic residual from sparse FFT coefficients
    if s.mask_e(b) == 1
        Eb = zeros(blk,1);
        sel = s.e_fft{b}.idx;
        Eb(sel+1) = (s.e_fft{b}.Re + 1i*s.e_fft{b}.Im) * s.e_fft{b}.qe;
        % enforce Hermitian symmetry for the negative frequencies
        % (mirror of 1-based FFT bin (sel+1) is bin (blk+1-sel))
        Eb(blk+1-sel) = conj(Eb(sel+1));
        xrec = xrec + real(ifft(Eb));
    end

    xhat(seg) = xrec;
end

xhat = xhat * s.scale;   % undo normalization
end

% ======================================================================
function v = dequant_uni(c, R, qstep) %#ok<INUSD>
    v = real(c)*qstep - 1 + 1i*(imag(c)*qstep - 1);
end

function y = interp_band(v, L)
    v = v(:);
    if numel(v) == 1
        y = repmat(v, L, 1);
    else
        y = interpft(v, L);
        y = y(:);
    end
end
