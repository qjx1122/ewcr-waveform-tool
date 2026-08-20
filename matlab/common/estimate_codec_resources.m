function R = estimate_codec_resources(algo,x,fs,cfg,opts,encoded,extra,bits)
%ESTIMATE_CODEC_RESOURCES Portable resource accounting for the five codecs.
%
% Outputs deliberately separate four different quantities:
%  1) logical_payload_bytes: actual codec bit-accounting / 8;
%  2) matlab_packet_bytes: in-memory MATLAB struct/cell representation;
%  3) encoder/decoder_working_est_bytes: deterministic estimate of dominant
%     numeric working arrays used by CURRENT MATLAB implementation;
%  4) edge_state_est_bytes: persistent algorithm state/buffer estimate.
%
% Working-memory estimates are lower-bound engineering estimates and exclude
% MATLAB runtime/JIT/toolbox overhead, allocator fragmentation, stack, OS,
% plotting and unrelated workspace variables. They are intended to compare
% algorithms, not to claim exact MCU/FPGA RAM.

if nargin<5, opts=struct(); end
N = numel(x);
a = upper(char(algo));
R = struct();
R.logical_payload_bytes = bits/8;
R.matlab_packet_bytes = var_bytes(encoded);
R.encoder_working_est_bytes = NaN;
R.decoder_working_est_bytes = NaN;
R.edge_state_est_bytes = NaN;
R.source_mfile_bytes = NaN;
R.note = '';

root = fileparts(fileparts(mfilename('fullpath')));
switch a
    case 'ASBC'
        folder=fullfile(root,'algorithms','ASBC');
        fft_len=2^nextpow2(N);
        blk=getfield_def(extra,'blk',min(N,round(5*fs/cfg.f0))); %#ok<GFLD>
        K=getfield_def(extra,'K',max(1,floor((fs/2)/cfg.f0))); %#ok<GFLD>
        % Major arrays in current encoder: x/xnorm, Xf/freqs, one complex
        % FFT-domain block chain yb/Yb/ylp plus block recon/residual.
        R.encoder_working_est_bytes = ...
            3*8*N + 2*8*fft_len + 4*16*blk + 4*8*blk + 8*max(1,K)*2;
        R.decoder_working_est_bytes = ...
            8*N + 2*16*blk + 3*8*blk + R.logical_payload_bytes;
        % FFT-block implementation needs one block plus per-subband activity
        % and small filter/modulation state between blocks.
        R.edge_state_est_bytes = 8*blk + 64*max(1,K) + 1024;
        R.source_mfile_bytes=dir_mfile_bytes(folder);
        R.note='Current MATLAB realization uses FFT-domain ideal filtering; a polyphase FIR FPGA/DSP implementation has different RAM/compute.';

    case {'DWT','DWT-HYBRID'}
        folder=fullfile(root,'algorithms','DWT_Hybrid');
        % wavedec/C/Cq/temp/symbol/gap work arrays: O(N). Huffman dynamic
        % objects vary; include packet RAM separately.
        R.encoder_working_est_bytes = 8*(8*N) + R.matlab_packet_bytes;
        R.decoder_working_est_bytes = 8*(6*N) + R.matlab_packet_bytes;
        % Current code is record/block DWT, so one transform block must be resident.
        R.edge_state_est_bytes = 8*(3*N) + 4096;
        R.source_mfile_bytes=dir_mfile_bytes(folder);
        R.note='Current implementation is batch/block DWT with MATLAB Huffman dictionaries; lifting/streaming DWT can materially reduce target RAM.';

    case {'CS','CS-OMP'}
        folder=fullfile(root,'algorithms','CS_OMP');
        p=merge_local(cfg.CS,opts);
        M=max(1,round(p.M_ratio*N));
        K0=p.K0;
        if strcmpi(p.basis,'dft')
            bPsi=16*N*N; % complex double
            bA=16*M*N;   % complex double
            bAs=16*M*K0;
            vec=16*(5*N+3*M);
        else
            bPsi=8*N*N;
            bA=8*M*N;
            bAs=8*M*K0;
            vec=8*(5*N+3*M);
        end
        bPhi=8*M*N;
        % Encoder currently performs OMP for diagnostics, so it carries nearly
        % the same dominant matrices as the decoder.
        R.encoder_working_est_bytes=bPhi+bPsi+bA+bAs+vec;
        R.decoder_working_est_bytes=bPhi+bPsi+bA+bAs+vec;
        % Persistent/transmitted state can be small (seed+measurements), but
        % current decoder reconstructs dense matrices in working RAM.
        R.edge_state_est_bytes=8*M + 8*K0 + 4096;
        R.source_mfile_bytes=dir_mfile_bytes(folder);
        R.note='This is the current dense-matrix MATLAB implementation. Structured/implicit sensing + FFT can reduce RAM by orders of magnitude; current encoder also runs OMP only for diagnostics.';

    case 'MMC'
        folder=fullfile(root,'algorithms','MMC');
        p=cfg.MMC; f=fieldnames(opts); for i=1:numel(f),p.(f{i})=opts.(f{i});end
        Nw=round(fs/cfg.f0);
        % Explicit orthonormal DCT matrix cached in current residual coder.
        bDct=8*Nw*Nw;
        bModels=8*Nw*30 + 128*1024; % conservative candidate/search arrays
        R.encoder_working_est_bytes=bDct+bModels+R.matlab_packet_bytes;
        R.decoder_working_est_bytes=bDct+8*Nw*8+R.logical_payload_bytes;
        R.edge_state_est_bytes=2*8*Nw + 16*1024; % two previous windows + prior model
        R.source_mfile_bytes=dir_mfile_bytes(folder);
        R.note='Encoder cost is dominated by candidate-model fitting and rate-distortion search; decoder evaluates only the selected model/residual and is much lighter.';

    case 'SVDCS'
        folder=fullfile(root,'algorithms','SVDCS');
        Nppc=round(fs/cfg.f0); nF=max(1,floor(N/Nppc)); Q=2*Nppc;
        % Current LZW encoder allocates transition(4096,256,uint16)=2 MiB.
        bLzw=4096*256*2;
        bFrames=8*nF*Nppc;
        bFFT=16*nF*Nppc;
        bAfloat=8*nF*Q;
        bAcode=4*nF*Q;
        bState=Q*(8+8+8+4); % MAX/MIN/SUM/CNT approximate integer state
        R.encoder_working_est_bytes=bLzw+bFrames+bFFT+bAfloat+bAcode+bState+R.matlab_packet_bytes;
        % Decoder's cell LZW dictionary is variable; use the same 2 MiB order
        % plus MR/FFT arrays as a reproducible comparison estimate.
        R.decoder_working_est_bytes=bLzw+8*nF*Q+16*nF*Nppc+R.matlab_packet_bytes;
        R.edge_state_est_bytes=bState+8*Nppc+64*1024;
        R.source_mfile_bytes=dir_mfile_bytes(folder);
        R.note='The 2 MiB term comes from this MATLAB fixed-table LZW implementation, not from the SVDCS paper FPGA. Paper reports logic utilization but does not give exact RAM used.';
    otherwise
        error('Unknown codec %s',algo);
end
end

function v=getfield_def(s,f,d)
if isstruct(s)&&isfield(s,f),v=s.(f);else,v=d;end
end
function a=merge_local(a,b)
f=fieldnames(b);for i=1:numel(f),a.(f{i})=b.(f{i});end
end
