function T = exp1_fixed_scenarios(cfg)
%EXP1_FIXED_SCENARIOS Ten IEEE-1159-style conditions, five codecs.
kinds={'pure','sag','swell','interruption','harmonics','osc_transient','notch','flicker','spike','complex'};
algos={'ASBC','DWT-Hybrid','CS-OMP','MMC','SVDCS'};
rows=struct([]);

for i=1:numel(kinds)
    x=gen_ieee1159_signals(kinds{i},cfg.fs,cfg.cycles);
    for j=1:numel(algos)
        fprintf('[Fixed] %-14s %-10s\n',kinds{i},algos{j});
        xx=x;
        % Quick mode is useful for development only; standard mode keeps
        % the full 10-cycle record for every codec.
        if cfg.quick_mode && strcmp(algos{j},'CS-OMP')
            xx=x(1:min(numel(x),round(4*cfg.fs/cfg.f0)));
        end
        r=run_codec_unified(algos{j},xx,cfg.fs,cfg,struct());
        row=result_row(r,'IEEE1159',kinds{i},'fixed');
        rows=[rows;row]; %#ok<AGROW>
    end
end
T=structrows_to_table(rows);
writetable(T,fullfile(cfg.results_dir,'exp1_fixed_scenarios.csv'));
save(fullfile(cfg.results_dir,'exp1_fixed_scenarios.mat'),'T');
end
