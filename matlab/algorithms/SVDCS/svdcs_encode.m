function out = svdcs_encode(x, cfg)
%SVDCS_ENCODE End-to-end MATLAB implementation of the SVDCS paper.
%
% out.packet contains everything required for synchronous-domain decoding.
% out.metrics evaluates the reconstructed synchronized waveform.
%
% Paper pipeline:
% input -> frequency/synchronization -> one-cycle frames -> FFT ->
% spectral variation (MAX/MIN/SUM/CNT) -> MC -> LZW

if nargin<2 || isempty(cfg)
    cfg=svdcs_default_config();
end

x=double(x(:));
xpu=x/cfg.signal_scale;

% 16-bit ADC emulation used in paper experiments.
if cfg.input_quantize
    [xcode, inStep]=svdcs_quantize_signed(xpu,cfg.input_bits,cfg.input_full_scale_pu);
    xpu=svdcs_dequantize_signed(xcode,inStep);
else
    inStep=NaN;
end

% Synchronous one-cycle segmentation.
[frames, prep]=svdcs_prepare_cycles(xpu,cfg);
[nFrames,Nppc]=size(frames);

% FFT per cycle. Division by Nppc gives coefficients in signal-amplitude
% units, making G=0.03 interpretable as a p.u.-scale tolerance.
X=fft(frames,[],2)/Nppc;

% Interleave Re/Im exactly as Eq. (2): [Re(bin1),Im(bin1),...,Re(binN),Im(binN)].
Afloat=zeros(nFrames,2*Nppc);
Afloat(:,1:2:end)=real(X);
Afloat(:,2:2:end)=imag(X);

if cfg.fft_quantize
    [Acode,fftStep]=svdcs_quantize_signed(Afloat,cfg.fft_bits,cfg.fft_full_scale_pu);
else
    % Use a very fine integer grid if quantization is disabled so the same
    % integer SUM/CNT code path remains available.
    fftStep=1e-12;
    Acode=int32(round(Afloat/fftStep));
end

Q=2*Nppc;
q=1:Q;
gamma=cfg.G*(q.^(-cfg.beta));
gammaCode=gamma/fftStep;

mc=svdcs_spectral_variation_encode(Acode,gammaCode);

% Paper's second stage: LZW applied to MC.
mcRaw=svdcs_serialize_mc(mc);
if cfg.use_lzw
    [payload,nCodes]=lzw_encode12(mcRaw);
    mcCheckBytes=lzw_decode12(payload,nCodes);
    if ~isequal(mcRaw,mcCheckBytes)
        error('Internal LZW round-trip failed.');
    end
else
    payload=mcRaw;
    nCodes=0;
end

packet.mc=mc;
packet.Nppc=Nppc;
packet.Nframes=nFrames;
packet.fftStep=fftStep;
packet.signal_scale=cfg.signal_scale;
packet.cfg=cfg;
packet.prep=prep;
packet.payload=payload;
packet.lzw_ncodes=nCodes;
packet.mc_raw_nbytes=numel(mcRaw);
packet.payload_nbytes=numel(payload);

% Independent logical decoder.
xhat=svdcs_decode(packet);

% Reference is the synchronized, ADC-quantized waveform.
xref=reshape(frames.',[],1)*cfg.signal_scale;
met=svdcs_metrics(xref,xhat);

origBits=numel(xref)*cfg.original_bits_per_sample;

% Since the paper does not specify byte-level MC packing metadata, report
% two CRs:
% 1) payload-only: closest to "LZW(MC)" concept in the paper.
% 2) engineering-total: payload + explicit compact metadata allowance.
metaBytes=64; % documented accounting placeholder for deployment metadata
payloadBits=max(1,8*numel(payload));
totalBits=max(1,8*(numel(payload)+metaBytes));

stats.original_bits=origBits;
stats.original_bytes=origBits/8;
stats.nFrames=nFrames;
stats.Nppc=Nppc;
stats.nMC=numel(mc.q);
stats.mc_raw_bytes=numel(mcRaw);
stats.lzw_payload_bytes=numel(payload);
stats.metadata_bytes_assumed=metaBytes;
stats.CR_preLZW=origBits/max(1,8*numel(mcRaw));
stats.CR_payload_only=origBits/payloadBits;
stats.CR_engineering_total=origBits/totalBits;
stats.novelty_density=(numel(mc.q)-Q)/(max(1,(nFrames-1)*Q));

out.packet=packet;
out.x_sync=xref;
out.xhat=xhat;
out.metrics=met;
out.stats=stats;
out.gamma=gamma;
out.gammaCode=gammaCode;

if cfg.verbose
    fprintf('\n===== SVDCS spectral-variation compression =====\n');
    fprintf('Frames              : %d\n',nFrames);
    fprintf('Samples/frame       : %d\n',Nppc);
    fprintf('MC columns          : %d\n',stats.nMC);
    fprintf('Novelty density     : %.6g\n',stats.novelty_density);
    fprintf('MC raw bytes        : %d\n',stats.mc_raw_bytes);
    fprintf('LZW payload bytes   : %d\n',stats.lzw_payload_bytes);
    fprintf('CR before LZW       : %.2f : 1\n',stats.CR_preLZW);
    fprintf('CR payload-only     : %.2f : 1\n',stats.CR_payload_only);
    fprintf('CR engineering      : %.2f : 1\n',stats.CR_engineering_total);
    fprintf('NMSE                : %.6g\n',met.NMSE);
    fprintf('COR (paper eq.)     : %.8f\n',met.COR);
    fprintf('RTE                 : %.6f %%\n',met.RTE_percent);
    fprintf('SNR                 : %.2f dB\n',met.SNR_dB);
end
end
