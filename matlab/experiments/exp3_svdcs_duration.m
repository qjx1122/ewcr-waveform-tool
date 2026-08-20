function T = exp3_svdcs_duration(cfg)
%EXP3_SVDCS_DURATION Demonstrate SVDCS duration/event-density dependence.
rows=struct([]);
for cyc=cfg.duration_cycles
    x=make_long_pattern(cfg.fs,cyc,'sparse');
    r=run_codec_unified('SVDCS',x,cfg.fs,cfg,struct('G',0.03));
    row=result_row(r,'long','sparse',sprintf('cycles=%d',cyc));
    rows=[rows;row]; %#ok<AGROW>
end
for cyc=cfg.duration_cycles(1:end-1)
    x=make_long_pattern(cfg.fs,cyc,'periodic_complex');
    r=run_codec_unified('SVDCS',x,cfg.fs,cfg,struct('G',0.03));
    row=result_row(r,'long','periodic_complex',sprintf('cycles=%d',cyc));
    rows=[rows;row]; %#ok<AGROW>
end

% Long-record rate-distortion curve for SVDCS on 200-cycle periodic complex data.
x=make_long_pattern(cfg.fs,200,'periodic_complex');
for G=cfg.rd.svdcs_G
    r=run_codec_unified('SVDCS',x,cfg.fs,cfg,struct('G',G));
    row=result_row(r,'long','periodic_complex_200',sprintf('G=%.3g',G));
    rows=[rows;row]; %#ok<AGROW>
end
T=structrows_to_table(rows);
writetable(T,fullfile(cfg.results_dir,'exp3_svdcs_duration.csv'));
save(fullfile(cfg.results_dir,'exp3_svdcs_duration.mat'),'T');
end
