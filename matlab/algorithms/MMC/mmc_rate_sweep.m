function tbl = mmc_rate_sweep(x, bits_per_window, cfg)
%MMC_RATE_SWEEP Evaluate rate-distortion performance on the same waveform.
% Example:
%   cfg=mmc_default_config(); cfg.verbose=false;
%   tbl=mmc_rate_sweep(x,[64 96 128 192 256 384],cfg)
if nargin<3 || isempty(cfg), cfg=mmc_default_config(); end
bits_per_window=bits_per_window(:);
R=zeros(numel(bits_per_window),1); S=R; RM=R; CR=R; BPS=R;
for k=1:numel(bits_per_window)
    c=cfg; c.n_tot=bits_per_window(k); c.verbose=false;
    o=mmc_encode_signal(x,c);
    R(k)=bits_per_window(k); S(k)=o.global_snr_db; RM(k)=o.mean_rmse;
    CR(k)=o.mean_cr; BPS(k)=o.mean_bps;
    fprintf('%4d bits/window | %6.3f bit/sample | SNR=%7.2f dB | CR=%7.2f:1\n',R(k),BPS(k),S(k),CR(k));
end
tbl=table(R,BPS,S,RM,CR,'VariableNames',{'BitsPerWindow','BitsPerSample','SNR_dB','RMSE','CompressionRatio'});
end
