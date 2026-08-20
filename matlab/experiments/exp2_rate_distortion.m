function T = exp2_rate_distortion(cfg)
%EXP2_RATE_DISTORTION Rate-distortion curves on the same complex waveform.
x=gen_ieee1159_signals('complex',cfg.fs,cfg.cycles);
rows=struct([]);

% ASBC: W sweep.
for W=cfg.rd.asbc_W
    r=run_codec_unified('ASBC',x,cfg.fs,cfg,struct('W',W,'Rk',8));
    rows=[rows;result_row(r,'IEEE1159','complex',sprintf('W=%g',W))]; %#ok<AGROW>
end
% ASBC: quantizer sweep at W=20 Hz.
for Rk=cfg.rd.asbc_Rk
    r=run_codec_unified('ASBC',x,cfg.fs,cfg,struct('W',20,'Rk',Rk));
    rows=[rows;result_row(r,'IEEE1159','complex',sprintf('Rk=%d',Rk))]; %#ok<AGROW>
end

for e=cfg.rd.dwt_energy
    r=run_codec_unified('DWT-Hybrid',x,cfg.fs,cfg,struct('qbits',8,'thr_mode','energy','thr_frac',e));
    rows=[rows;result_row(r,'IEEE1159','complex',sprintf('E=%.4f',e))]; %#ok<AGROW>
end

for mr=cfg.rd.cs_ratio
    r=run_codec_unified('CS-OMP',x,cfg.fs,cfg,struct('M_ratio',mr,'basis','dft','K0',48));
    rows=[rows;result_row(r,'IEEE1159','complex',sprintf('M/N=%.2f',mr))]; %#ok<AGROW>
end

for nt=cfg.rd.mmc_bits
    r=run_codec_unified('MMC',x,cfg.fs,cfg,struct('n_tot',nt));
    rows=[rows;result_row(r,'IEEE1159','complex',sprintf('n_tot=%d',nt))]; %#ok<AGROW>
end

for G=cfg.rd.svdcs_G
    r=run_codec_unified('SVDCS',x,cfg.fs,cfg,struct('G',G));
    rows=[rows;result_row(r,'IEEE1159','complex',sprintf('G=%.3g',G))]; %#ok<AGROW>
end

T=structrows_to_table(rows);
writetable(T,fullfile(cfg.results_dir,'exp2_rate_distortion.csv'));
save(fullfile(cfg.results_dir,'exp2_rate_distortion.mat'),'T');
end
