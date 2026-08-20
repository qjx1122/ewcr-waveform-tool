function met = svdcs_metrics(x, xhat)
%SVDCS_METRICS Reconstruction metrics used in the paper plus SNR.
%
% Paper definitions:
% NMSE = sum |x-xhat|^2 / sum |x|^2
% COR  = (x^T xhat) / (x^T x)
% RTE  = 100 * sum xhat^2 / sum x^2
%
% COR is implemented exactly as printed in the paper, not as MATLAB corr().

x = double(x(:));
xhat = double(xhat(:));
N = min(numel(x),numel(xhat));
x = x(1:N);
xhat = xhat(1:N);

den = sum(abs(x).^2);
err = x-xhat;

met.N = N;
met.NMSE = sum(abs(err).^2) / max(den,eps);
met.COR = dot(x,xhat) / max(dot(x,x),eps);
met.RTE_percent = 100 * sum(abs(xhat).^2) / max(den,eps);
met.SNR_dB = 10*log10(max(den,eps)/max(sum(abs(err).^2),eps));
met.RMSE = sqrt(mean(abs(err).^2));
met.max_abs_error = max(abs(err));
end
