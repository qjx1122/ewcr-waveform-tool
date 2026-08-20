%% DEMO_SVDCS_6400HZ
% Direct setup for 6.4 kHz / 50 Hz electrical waveform data.
%
% Because 6400/50 = 128 exactly, each nominal fundamental cycle naturally
% contains Nppc=128 samples and is especially convenient for SVDCS.

clear; clc; close all;

cfg=svdcs_default_config();
cfg.fs=6400;
cfg.f0=50;
cfg.Nppc=128;
cfg.sync_mode='nominal';
cfg.G=0.03;
cfg.beta=0.10;
cfg.signal_scale=1;     % set to nominal peak if your x is in volts/amps
cfg.verbose=true;

fs=cfg.fs;
t=(0:1/fs:30-1/fs).';

x=sin(2*pi*50*t);
x=x+0.04*sin(2*pi*150*t)+0.025*sin(2*pi*250*t);

% Load-like change.
m=t>=10 & t<18;
x(m)=0.82*sin(2*pi*50*t(m))+0.10*sin(2*pi*250*t(m));

% Short switching transient.
m=t>=22 & t<22.08;
tau=t(m)-22;
x(m)=x(m)+0.18*exp(-tau/0.018).*sin(2*pi*700*tau);

out=svdcs_encode(x,cfg);

fprintf('\nFor your own 6.4 kHz waveform:\n');
fprintf('  replace x above by your measured vector.\n');
fprintf('  if x is in physical units, set cfg.signal_scale to nominal peak.\n');
