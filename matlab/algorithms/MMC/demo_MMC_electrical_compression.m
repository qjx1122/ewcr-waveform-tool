clear; clc; close all;
% Demo: MATLAB engineering implementation of
% "Multiple-Model Coding Scheme for Electrical Signal Compression"

cfg = mmc_default_config();
cfg.verbose = true;
cfg.n_tot = 128;          % 1 bit/sample for N=128
cfg.max_windows = 25;     % quick demo; set inf for full record

fs=cfg.fs; t=(0:1/fs:1-1/fs)';
% Synthetic electrical waveform: fundamental + harmonics + amplitude step + transient
x = 230*sqrt(2)*sin(2*pi*50*t+0.15) ...
  + 0.04*230*sqrt(2)*sin(2*pi*150*t+0.7) ...
  + 0.025*230*sqrt(2)*sin(2*pi*250*t-0.4);
x(t>=0.18 & t<0.28)=0.82*x(t>=0.18 & t<0.28);          % voltage sag
x = x + 35*exp(-220*(t-0.32)).*sin(2*pi*850*(t-0.32)).*(t>=0.32); % transient
rng(1); x=x+0.6*randn(size(x));

out = mmc_encode_signal(x,cfg);
xdec = mmc_decode_signal(out.packets,cfg);
roundtrip_err = max(abs(xdec-out.reconstructed));

fprintf('\n===== MMC summary =====\n');
fprintf('Windows            : %d\n',numel(out.packets));
fprintf('Global SNR         : %.2f dB\n',out.global_snr_db);
fprintf('Mean RMSE          : %.4f\n',out.mean_rmse);
fprintf('Mean bits/sample   : %.4f\n',out.mean_bps);
fprintf('Mean CR (16-bit)   : %.2f : 1\n',out.mean_cr);
fprintf('Decoder roundtrip  : max |enc-dec| = %.3g\n',roundtrip_err);

n=numel(out.original); tt=(0:n-1)'/fs;
figure('Name','MMC reconstruction');
plot(tt,out.original,'LineWidth',1); hold on;
plot(tt,out.reconstructed,'--','LineWidth',1);
xlabel('Time (s)'); ylabel('Amplitude'); grid on;
legend('Original','MMC reconstructed');
title(sprintf('MMC: SNR %.2f dB, mean CR %.2f:1',out.global_snr_db,out.mean_cr));

figure('Name','MMC error');
plot(tt,out.original-out.reconstructed,'LineWidth',1);
xlabel('Time (s)'); ylabel('Reconstruction error'); grid on;

figure('Name','Per-window performance');
plot([out.metrics.snr_db],'-o','LineWidth',1);
xlabel('Window'); ylabel('SNR (dB)'); grid on;
