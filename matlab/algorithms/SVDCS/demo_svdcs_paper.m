%% DEMO_SVDCS_PAPER
% Paper-oriented demonstration of:
% E. B. Kapisch et al.,
% "Spectral Variation-Based Signal Compression Technique for Gapless
% Power Quality Waveform Recording in Smart Grids,"
% IEEE Transactions on Industrial Informatics, 18(7), 2022.
%
% Paper synthetic setup:
%   f0 = 60 Hz, Nppc = 128, fs = 7.680 kHz
%   G = 0.03, beta = 0.1
%   16-bit samples and FFT real/imag parts
%
% This demo is intentionally much shorter than the paper's 10--13 minute
% cases so it runs quickly.

clear; clc; close all;

cfg=svdcs_default_config();
cfg.fs=7680;
cfg.f0=60;
cfg.Nppc=128;
cfg.sync_mode='nominal';
cfg.G=0.03;
cfg.beta=0.1;
cfg.signal_scale=1; % generated signal is already p.u.
cfg.verbose=true;

fs=cfg.fs;
f0=cfg.f0;
dur=20; % seconds
t=(0:1/fs:dur-1/fs).';

% Base waveform.
amp=ones(size(t));
phase=zeros(size(t));
x=amp.*sin(2*pi*f0*t+phase);

% Sag.
m=t>=2 & t<3.5;
x(m)=0.70*sin(2*pi*f0*t(m));

% Swell.
m=t>=5 & t<6.2;
x(m)=1.30*sin(2*pi*f0*t(m));

% Odd/even harmonic-rich interval.
m=t>=8 & t<10;
x(m)=sin(2*pi*f0*t(m)) ...
    +0.10*sin(2*pi*3*f0*t(m)+0.2) ...
    +0.06*sin(2*pi*5*f0*t(m)-0.4) ...
    +0.04*sin(2*pi*2*f0*t(m)+0.1);

% Time-varying interharmonic.
m=t>=11 & t<13;
tau=t(m)-11;
fi=190+40*(tau/max(tau));
phi_i=2*pi*cumsum(fi)/fs;
x(m)=sin(2*pi*f0*t(m))+0.06*sin(phi_i);

% Damped oscillatory transient.
m=t>=14 & t<14.15;
tau=t(m)-14;
x(m)=x(m)+0.25*exp(-tau/0.025).*sin(2*pi*900*tau);

% Phase jump.
m=t>=16;
x(m)=sin(2*pi*f0*t(m)+15*pi/180);

% Small measurement noise.
rng(2);
sigp=mean(x.^2);
noisePower=sigp/10^(50/10);
x=x+sqrt(noisePower)*randn(size(x));

out=svdcs_encode(x,cfg);

% Verify reconstruction from the actual LZW byte payload.
xhat2=svdcs_decode_from_payload(out.packet);
fprintf('Payload-decoder max difference: %.3g\n',max(abs(xhat2-out.xhat)));

% Plot a representative interval.
tt=(0:numel(out.x_sync)-1).'/fs;
figure('Name','SVDCS reconstruction');
plot(tt,out.x_sync,'-'); hold on;
plot(tt,out.xhat,'--');
xlabel('Time (s)'); ylabel('Amplitude (p.u.)');
legend('Original synchronized','SVDCS reconstructed');
title(sprintf('SVDCS: CR=%.1f:1, NMSE=%.3g', ...
    out.stats.CR_payload_only,out.metrics.NMSE));
grid on;
xlim([7.8 10.2]);

figure('Name','Reconstruction error');
plot(tt,out.x_sync-out.xhat);
xlabel('Time (s)'); ylabel('Error (p.u.)');
title('SVDCS reconstruction error');
grid on;

% Rate-quality sweep through G.
Gvalues=[0.005 0.01 0.02 0.03 0.05 0.08];
T=svdcs_threshold_sweep(x,cfg,Gvalues);
disp(T);

figure('Name','CR-NMSE trade-off');
semilogy(T.CR,T.NMSE,'o-');
xlabel('Compression ratio, payload-only');
ylabel('NMSE');
title('SVDCS threshold trade-off');
grid on;
