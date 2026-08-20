function [T,Tscale,Env] = exp5_runtime_memory(cfg)
%EXP5_RUNTIME_MEMORY Encode/decode time and edge-resource accounting.
%
% IMPORTANT: MATLAB timing is machine-specific. This experiment writes the
% environment metadata alongside results. Use median warm timings for
% cross-algorithm comparison on one machine. Memory columns are deterministic
% estimates of dominant algorithm arrays and exclude MATLAB/JIT runtime.
algos={'ASBC','DWT-Hybrid','CS-OMP','MMC','SVDCS'};
rows=struct([]);

Env=struct();
Env.matlab_version=version;
Env.computer=computer;
Env.arch=computer('arch');
Env.timestamp=datestr(now,'yyyy-mm-dd HH:MM:SS');
Env.fs=cfg.fs; Env.f0=cfg.f0; Env.raw_bits=cfg.raw_bits_per_sample;
Env.processor_identifier=getenv('PROCESSOR_IDENTIFIER');
Env.processor_architecture=getenv('PROCESSOR_ARCHITECTURE');
try
    Env.cpu=char(java.lang.System.getProperty('os.arch'));
catch
    Env.cpu='unavailable';
end

for ii=1:numel(cfg.profile.scenarios)
    sig=cfg.profile.scenarios{ii};
    x=gen_ieee1159_signals(sig,cfg.fs,cfg.cycles);
    for j=1:numel(algos)
        fprintf('[Profile] %-14s %-10s\n',sig,algos{j});
        P=profile_codec_resources(algos{j},x,cfg.fs,cfg,struct());
        rows=[rows; p_row(P,'profile',sig,sprintf('%d cycles',cfg.cycles))]; %#ok<AGROW>
    end
end
T=structrows_to_table(rows);
writetable(T,fullfile(cfg.results_dir,'exp5_runtime_memory.csv'));
save(fullfile(cfg.results_dir,'exp5_runtime_memory.mat'),'T','Env');
fid=fopen(fullfile(cfg.results_dir,'exp5_environment.txt'),'w');
if fid>0
    fprintf(fid,'MATLAB version: %s\n',Env.matlab_version);
    fprintf(fid,'Computer: %s\n',Env.computer);
    fprintf(fid,'Arch: %s\n',Env.arch);
    fprintf(fid,'Timestamp: %s\n',Env.timestamp);
    fprintf(fid,'PROCESSOR_IDENTIFIER: %s\n',Env.processor_identifier);
    fprintf(fid,'PROCESSOR_ARCHITECTURE: %s\n',Env.processor_architecture);
    fprintf(fid,'fs: %.12g Hz\n',Env.fs);
    fprintf(fid,'f0: %.12g Hz\n',Env.f0);
    fprintf(fid,'raw bits/sample: %d\n',Env.raw_bits);
    fclose(fid);
end

% Scaling experiment: how time/RAM grow with record length.
rows2=struct([]);
if cfg.profile.run_scaling
    for cyc=cfg.profile.scaling_cycles
        x=gen_ieee1159_signals('complex',cfg.fs,cyc);
        for j=1:numel(algos)
            fprintf('[Profile-scale] %2d cyc %-10s\n',cyc,algos{j});
            P=profile_codec_resources(algos{j},x,cfg.fs,cfg,struct());
            rows2=[rows2; p_row(P,'scaling','complex',sprintf('%d cycles',cyc))]; %#ok<AGROW>
        end
    end
end
Tscale=structrows_to_table(rows2);
if ~isempty(Tscale)
    writetable(Tscale,fullfile(cfg.results_dir,'exp5_runtime_memory_scaling.csv'));
end
save(fullfile(cfg.results_dir,'exp5_runtime_memory_scaling.mat'),'Tscale');
end

function r=p_row(P,dataset,signal,param)
% Build a scalar row struct using field-by-field assignment for maximum
% compatibility across MATLAB releases. Text fields are plain char.
r=struct();
r.dataset=char(dataset);
r.signal=char(signal);
r.algo=char(P.algo);
r.param=char(param);

r.N=P.N;
r.CR=P.CR;
r.SNR_dB=P.SNR_dB;
r.signal_duration_s=P.signal_duration_s;
r.cold_enc_s=P.cold_enc_s;
r.cold_dec_s=P.cold_dec_s;
r.enc_median_s=P.enc_median_s;
r.enc_mean_s=P.enc_mean_s;
r.enc_p90_s=P.enc_p90_s;
r.dec_median_s=P.dec_median_s;
r.dec_mean_s=P.dec_mean_s;
r.dec_p90_s=P.dec_p90_s;
r.enc_realtime_factor=P.enc_realtime_factor;
r.dec_realtime_factor=P.dec_realtime_factor;
r.enc_throughput_Msps=P.enc_throughput_Msps;
r.dec_throughput_Msps=P.dec_throughput_Msps;
r.logical_payload_bytes=P.logical_payload_bytes;
r.matlab_packet_bytes=P.matlab_packet_bytes;
r.encoder_working_est_bytes=P.encoder_working_est_bytes;
r.decoder_working_est_bytes=P.decoder_working_est_bytes;
r.edge_state_est_bytes=P.edge_state_est_bytes;
r.source_mfile_bytes=P.source_mfile_bytes;
r.resource_note=char(P.resource_note);
end
