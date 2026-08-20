function out = mmc_encode_signal(x, cfg)
%MMC_ENCODE_SIGNAL Encode one electrical waveform channel window-by-window.
% x can be a row or column vector. Only complete N-sample windows are used.

if nargin < 2 || isempty(cfg), cfg = mmc_default_config(); end
x = x(:);
N = cfg.N;
W = floor(numel(x)/N);
W = min(W, cfg.max_windows);
if W < 1, error('Signal must contain at least cfg.N samples.'); end
x = x(1:W*N);

xhat = zeros(size(x));
packets = cell(W,1);
metrics = repmat(struct('snr_db',0,'rmse',0,'total_bits',0,'cr',0, ...
    'model','','residual','','nx',0,'nr',0,'kx',0), W, 1);

state = struct();
state.prev2 = zeros(N,1);
state.prev1 = zeros(N,1);
state.prev_model = [];
state.history_count = 0;

for w = 1:W
    idx = (w-1)*N + (1:N);
    xw = x(idx);
    [packet, xw_hat, state_next] = mmc_encode_window(xw, state, cfg);
    xhat(idx) = xw_hat;
    packets{w} = packet;

    e = xw - xw_hat;
    metrics(w).snr_db = 10*log10(sum(xw.^2)/(sum(e.^2)+eps));
    metrics(w).rmse = sqrt(mean(e.^2));
    metrics(w).total_bits = packet.total_bits;
    metrics(w).cr = cfg.raw_bits_per_sample*N / max(packet.total_bits,1);
    metrics(w).model = packet.model.name;
    metrics(w).residual = packet.residual.name;
    metrics(w).nx = packet.nx;
    metrics(w).nr = packet.nr;
    metrics(w).kx = packet.kx;

    if cfg.verbose
        fprintf('w=%3d | bits=%3d | SNR=%7.2f dB | CR=%6.2f | %-12s | %-4s | nx=%2d nr=%3d\n', ...
            w, packet.total_bits, metrics(w).snr_db, metrics(w).cr, ...
            packet.model.name, packet.residual.name, packet.nx, packet.nr);
    end
    state = state_next;
end

out = struct();
out.original = x;
out.reconstructed = xhat;
out.packets = packets;
out.metrics = metrics;
out.cfg = cfg;
out.mean_snr_db = mean([metrics.snr_db]);
out.global_snr_db = 10*log10(sum(x.^2)/(sum((x-xhat).^2)+eps));
out.mean_rmse = mean([metrics.rmse]);
out.mean_bits_per_window = mean([metrics.total_bits]);
out.mean_bps = out.mean_bits_per_window / N;
out.mean_cr = cfg.raw_bits_per_sample*N / out.mean_bits_per_window;
end
