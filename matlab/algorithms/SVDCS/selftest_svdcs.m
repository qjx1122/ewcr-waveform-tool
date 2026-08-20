%% SELFTEST_SVDCS
% Lightweight consistency tests. No toolbox is required.

clear; clc;

% LZW test
rng(1);
b=uint8([repmat(uint8('ABRACADABRA'),1,50), randi([0 255],1,500,'uint8')]);
[p,nc]=lzw_encode12(b);
br=lzw_decode12(p,nc);
assert(isequal(b,br),'LZW round-trip failed.');

% SVDCS test on exactly periodic 50 Hz waveform.
cfg=svdcs_default_config();
cfg.fs=6400;
cfg.f0=50;
cfg.Nppc=128;
cfg.sync_mode='nominal';
cfg.G=0.03;
cfg.verbose=false;

t=(0:1/cfg.fs:2-1/cfg.fs).';
x=sin(2*pi*50*t)+0.05*sin(2*pi*250*t);

o=svdcs_encode(x,cfg);
xr=svdcs_decode_from_payload(o.packet);

assert(numel(xr)==numel(o.x_sync),'Decoded length mismatch.');
assert(max(abs(xr-o.xhat))<1e-10,'Payload decoding differs from direct MC decoding.');
assert(o.metrics.COR>0.99,'Unexpectedly poor reconstruction.');
assert(o.metrics.RTE_percent>95 && o.metrics.RTE_percent<105,'Unexpected RTE.');

fprintf('All SVDCS self-tests passed.\n');
fprintf('CR payload-only = %.2f:1, NMSE = %.3g, SNR = %.2f dB\n', ...
    o.stats.CR_payload_only,o.metrics.NMSE,o.metrics.SNR_dB);
