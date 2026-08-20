function T = exp4_noise_robustness(cfg)
%EXP4_NOISE_ROBUSTNESS All five codecs under measurement-noise levels.
base=gen_ieee1159_signals('complex',cfg.fs,cfg.cycles);
algos={'ASBC','DWT-Hybrid','CS-OMP','MMC','SVDCS'};
rows=struct([]);
for s=cfg.noise_snr_db
    x=add_awgn_snr(base,s,cfg.seed);
    for j=1:numel(algos)
        fprintf('[Noise] SNR=%g %-10s\n',s,algos{j});
        r=run_codec_unified(algos{j},x,cfg.fs,cfg,struct());
        if isinf(s), p='inputSNR=Inf'; else, p=sprintf('inputSNR=%g',s); end
        rows=[rows;result_row(r,'noise','complex',p)]; %#ok<AGROW>
    end
end
T=structrows_to_table(rows);
writetable(T,fullfile(cfg.results_dir,'exp4_noise_robustness.csv'));
save(fullfile(cfg.results_dir,'exp4_noise_robustness.mat'),'T');
end
