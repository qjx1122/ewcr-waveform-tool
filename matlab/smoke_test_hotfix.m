%% SMOKE_TEST_HOTFIX
% Fast five-codec regression test for v2.6.
%
% It intentionally exercises the complete encode -> decode adapter path for:
%   ASBC, DWT-Hybrid, CS-OMP, MMC, SVDCS
%
% It also checks:
%   - MMC unified result remains a scalar struct although encoded packets are cells
%   - MMC does not emit the previous rank-deficiency warning on first window
%   - SVDCS MC byte serialization/deserialization and LZW payload round-trip
%   - unified result row / table conversion

clear; clc;

root=fileparts(mfilename('fullpath'));
addpath(genpath(root));

cfg=benchmark_config();
cfg.cycles=2;                 % deliberately short smoke signal
cfg.profile.run_scaling=false;

fprintf('Running v2.6 five-codec smoke test...\n');

x=gen_ieee1159_signals('pure',cfg.fs,cfg.cycles);
algos={'ASBC','DWT-Hybrid','CS-OMP','MMC','SVDCS'};

% Do not initialize with an empty-field struct: MATLAB does not allow
% assigning a populated struct into a structurally different empty struct.
% Instead, use the first real result row to establish the schema.
rows=[];

for k=1:numel(algos)
    a=algos{k};
    fprintf('  %-12s ... ',a);

    lastwarn('');
    r=run_codec_unified(a,x,cfg.fs,cfg,struct());

    assert(isstruct(r) && isscalar(r), ...
        'run_codec_unified(%s) must return exactly one scalar struct.',a);
    assert(ischar(r.algo) && strcmp(r.algo,a), ...
        'Unified result algo metadata mismatch for %s.',a);
    assert(numel(r.xhat)>0 && all(isfinite(r.xhat)), ...
        '%s decoder produced empty or nonfinite output.',a);
    assert(isfinite(r.metrics.CR) && r.metrics.CR>0, ...
        '%s produced invalid compression ratio.',a);
    assert(isfinite(r.metrics.NMSE_dB), ...
        '%s produced invalid NMSE_dB.',a);

    if strcmp(a,'MMC')
        assert(iscell(r.encoded), ...
            'MMC r.encoded should contain the complete packet cell array.');
        assert(numel(r.encoded)==cfg.cycles, ...
            'MMC packet count should equal the number of complete test cycles.');
        [msg,~]=lastwarn;
        if ~isempty(msg) && ~isempty(strfind(lower(msg),'rank')) %#ok<STREMP>
            error('FiveCodecBenchmark:RankWarning', ...
                'MMC still emitted a rank-related warning: %s',msg);
        end
    end

    if strcmp(a,'SVDCS')
        pkt=r.encoded;
        assert(isfield(pkt,'payload') && isa(pkt.payload,'uint8'), ...
            'SVDCS payload is missing or is not uint8.');
        assert(pkt.payload_nbytes==numel(pkt.payload), ...
            'SVDCS payload byte count is inconsistent.');
        % Force the actual payload decoder path again.
        xr=svdcs_decode_from_payload(pkt);
        assert(numel(xr)==numel(r.xhat), ...
            'SVDCS payload decoder length mismatch.');
        assert(max(abs(xr-r.xhat)) < 1e-10, ...
            'SVDCS payload decoder disagrees with unified decoder.');
    end

    rowk=result_row(r,'IEEE1159','pure','smoke');
    if k==1
        rows=repmat(rowk,numel(algos),1);   % schema-safe preallocation
    else
        rows(k,1)=rowk;
    end
    fprintf('PASS (CR=%.3g, SNR=%.2f dB)\n',r.metrics.CR,r.metrics.SNR_dB);
end

T=structrows_to_table(rows);
assert(height(T)==5,'Five-codec smoke table must contain five rows.');
assert(strcmp(T.algo{1},'ASBC'));
assert(strcmp(T.algo{4},'MMC'));
assert(strcmp(T.algo{5},'SVDCS'));

disp(T(:,{'algo','CR','SNR_dB','enc_time_s','dec_time_s', ...
    'logical_payload_bytes','encoder_working_est_bytes','decoder_working_est_bytes'}));

fprintf('v2.6 five-codec smoke test PASSED.\n');
