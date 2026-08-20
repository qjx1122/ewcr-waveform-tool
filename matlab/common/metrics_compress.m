function m = metrics_compress(x, xhat, bits_original, bits_compressed, t)
%METRICS_COMPRESS  Compute compression performance metrics.
%   m = METRICS_COMPRESS(x, xhat, bits_original, bits_compressed)
%   returns a struct with:
%     CR    : compression ratio (bits_original / bits_compressed)
%     NMSE  : normalized mean-squared error (dB)
%     SNR   : reconstruction signal-to-noise ratio (dB) = -NMSE
%     PSNR  : peak signal-to-noise ratio (dB)
%     PRD   : percent root-mean-square difference (%)
%     MAXE  : maximum absolute error
%     RMSE  : root-mean-squared error
%   x is the original real vector, xhat the reconstruction.
arguments
    x (:,1) double {mustBeFinite}
    xhat (:,1) double {mustBeFinite}
    bits_original (1,1) double {mustBePositive}
    bits_compressed (1,1) double {mustBePositive}
    t (1,1) double = NaN
end

if ~isequal(size(x), size(xhat))
    error('metrics_compress:sizeMismatch', ...
        'x and xhat sizes differ: %s vs %s', mat2str(size(x)), mat2str(size(xhat)));
end

N = numel(x);
e  = x - xhat;
nmse_lin = sum(e.^2) / sum(x.^2);
m.CR    = bits_original / bits_compressed;
m.NMSE  = 10*log10(nmse_lin);
m.SNR   = -m.NMSE;
m.RMSE  = sqrt(mean(e.^2));
m.MAXE  = max(abs(e));
m.PRD   = 100 * sqrt(sum(e.^2)/sum(x.^2));
peak    = max(abs(x));
if peak > 0 && m.RMSE > 0
    m.PSNR = 20*log10(peak / m.RMSE);
else
    m.PSNR = Inf;
end
if ~isnan(t)
    m.elapsed_s = t;
end
end
