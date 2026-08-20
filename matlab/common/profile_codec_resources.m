function P = profile_codec_resources(algo,x,fs,cfg,opts)
%PROFILE_CODEC_RESOURCES Repeated encode/decode timing + resource accounting.
%
% Timing protocol:
% - one cold run is recorded separately;
% - cfg.profile.warmup warm-up runs are discarded;
% - cfg.profile.repeats steady-state runs are summarized by median/p90;
% - encode and decode are timed separately by run_codec_unified.
%
% The timings are MATLAB/host-specific. Record MATLAB version, computer and
% CPU in exp5_runtime_memory.m before comparing across machines.
if nargin<5,opts=struct();end

% Cold run includes first-call/JIT/cache effects.
r0=run_codec_unified(algo,x,fs,cfg,opts);

for k=1:cfg.profile.warmup
    run_codec_unified(algo,x,fs,cfg,opts);
end

ne=cfg.profile.repeats;
te=zeros(ne,1); td=zeros(ne,1);
r=r0;
for k=1:ne
    r=run_codec_unified(algo,x,fs,cfg,opts);
    te(k)=r.enc_time_s; td(k)=r.dec_time_s;
end

P=struct();
P.algo=char(algo);
P.cold_enc_s=r0.enc_time_s;
P.cold_dec_s=r0.dec_time_s;
P.enc_median_s=median(te);
P.enc_mean_s=mean(te);
P.enc_p90_s=percentile_simple(te,90);
P.dec_median_s=median(td);
P.dec_mean_s=mean(td);
P.dec_p90_s=percentile_simple(td,90);
P.signal_duration_s=numel(x)/fs;
P.enc_realtime_factor=P.enc_median_s/P.signal_duration_s;
P.dec_realtime_factor=P.dec_median_s/P.signal_duration_s;
P.enc_throughput_Msps=numel(x)/P.enc_median_s/1e6;
P.dec_throughput_Msps=numel(x)/P.dec_median_s/1e6;
P.logical_payload_bytes=r.resource.logical_payload_bytes;
P.matlab_packet_bytes=r.resource.matlab_packet_bytes;
P.encoder_working_est_bytes=r.resource.encoder_working_est_bytes;
P.decoder_working_est_bytes=r.resource.decoder_working_est_bytes;
P.edge_state_est_bytes=r.resource.edge_state_est_bytes;
P.source_mfile_bytes=r.resource.source_mfile_bytes;
P.CR=r.metrics.CR;
P.SNR_dB=r.metrics.SNR_dB;
P.N=r.metrics.N;
P.resource_note=r.resource.note;
end
