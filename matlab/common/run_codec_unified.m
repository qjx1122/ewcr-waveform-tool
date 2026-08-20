function r = run_codec_unified(algo,x,fs,cfg,opts)
%RUN_CODEC_UNIFIED Common adapter for ASBC, DWT, CS-OMP, MMC and SVDCS.
% Returns unified distortion/rate metrics, one-shot encode/decode timings,
% the encoded MATLAB object, and deterministic resource estimates.
if nargin<5, opts=struct(); end
x=x(:);
tEnc=tic;

switch upper(algo)
    case 'ASBC'
        p=merge_struct(cfg.ASBC,opts);
        [s,info]=asbc_encode(x,fs,p);
        t1=toc(tEnc); tDec=tic; xhat=asbc_decode(s); t2=toc(tDec);
        bits=s.bits_compressed; extra=info; encoded=s;

    case {'DWT','DWT-HYBRID'}
        p=merge_struct(cfg.DWT,opts);
        [s,info]=dwt_hybrid_encode(x,p);
        t1=toc(tEnc); tDec=tic; xhat=dwt_hybrid_decode(s); t2=toc(tDec);
        bits=s.bits_compressed; extra=info; encoded=s;

    case {'CS','CS-OMP'}
        p=merge_struct(cfg.CS,opts);
        [s,info]=cs_omp_encode(x,p);
        t1=toc(tEnc); tDec=tic; xhat=cs_omp_decode(s); t2=toc(tDec);
        bits=s.bits_compressed; extra=info; encoded=s;

    case 'MMC'
        p=cfg.MMC;
        f=fieldnames(opts); for i=1:numel(f), p.(f{i})=opts.(f{i}); end
        p.fs=fs; p.fn=cfg.f0; p.N=round(fs/cfg.f0); p.verbose=false;
        out=mmc_encode_signal(x,p);
        t1=toc(tEnc); tDec=tic; xhat=mmc_decode_signal(out.packets,p); t2=toc(tDec);
        bits=sum(cellfun(@(z) z.total_bits,out.packets));
        extra=struct('mean_bits_per_window',out.mean_bits_per_window,'mean_bps',out.mean_bps);
        encoded=out.packets;

    case 'SVDCS'
        p=cfg.SVDCS;
        f=fieldnames(opts); for i=1:numel(f), p.(f{i})=opts.(f{i}); end
        p.fs=fs; p.f0=cfg.f0; p.Nppc=round(fs/cfg.f0); p.verbose=false;
        out=svdcs_encode(x,p);
        t1=toc(tEnc); tDec=tic; xhat=svdcs_decode_from_payload(out.packet); t2=toc(tDec);
        switch lower(cfg.svdcs_rate_mode)
            case 'payload'
                bits=8*out.stats.lzw_payload_bytes;
            case 'engineering'
                bits=8*(out.stats.lzw_payload_bytes+out.stats.metadata_bytes_assumed);
            otherwise
                error('Unknown cfg.svdcs_rate_mode.');
        end
        extra=out.stats; encoded=out.packet;

    otherwise
        error('Unknown codec %s.',algo);
end

N=min(numel(x),numel(xhat));
xref=x(1:N); xhat=xhat(1:N);
m=unified_metrics(xref,xhat,bits,fs,cfg.f0,cfg.raw_bits_per_sample);
res=estimate_codec_resources(algo,x,fs,cfg,opts,encoded,extra,bits);
% IMPORTANT: construct the unified result by field assignment.
% Do NOT use struct('encoded',encoded,...) here because MMC encoded data are
% a cell array (one packet per window). MATLAB's struct(name,value,...) treats
% a cell value as a request to build a struct ARRAY. That silently expands r
% into one struct per MMC packet and breaks all downstream metadata/table code.
fprintf('    %-12s | Encode: %9.4f ms | Decode: %9.4f ms | Total: %9.4f ms\n', ...
    char(algo), 1000*t1, 1000*t2, 1000*(t1+t2));

r = struct();
r.algo = char(algo);
r.xhat = xhat;
r.metrics = m;
r.enc_time_s = t1;
r.dec_time_s = t2;
r.extra = extra;
r.encoded = encoded;
r.bits_compressed = bits;
r.resource = res;
end

function a=merge_struct(a,b)
f=fieldnames(b); for i=1:numel(f), a.(f{i})=b.(f{i}); end
end
