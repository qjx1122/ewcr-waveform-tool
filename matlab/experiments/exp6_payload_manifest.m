function T = exp6_payload_manifest(cfg) %#ok<INUSD>
%EXP6_PAYLOAD_MANIFEST Export human-readable compressed-payload audit.
T=codec_payload_manifest();
writetable(T,fullfile(cfg.results_dir,'exp6_payload_manifest.csv'));
end
