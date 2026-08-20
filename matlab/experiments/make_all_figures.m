function make_all_figures(cfg,T1,T2,T3,T4)
%MAKE_ALL_FIGURES Publication-oriented benchmark figures.
if ~exist(cfg.figures_dir,'dir'),mkdir(cfg.figures_dir);end
algos=unique(T2.algo,'stable');
marks={'o','s','^','d','v'};

% --- Rate-distortion
f=figure('Visible','off','Position',[80 80 850 600]); hold on; grid on;
for i=1:numel(algos)
    z=T2(strcmp(T2.algo,algos{i}),:);
    [xr,ord]=sort(z.rate_bps);
    plot(xr,z.NMSE_dB(ord),['-' marks{i}],'LineWidth',1.3,'MarkerSize',6,'DisplayName',algos{i});
end
xlabel('Rate R (compressed bits/sample)'); ylabel('NMSE (dB, lower is better)');
title('Five-codec rate-distortion comparison: IEEE1159 complex signal'); legend('Location','best');
save_figure_compat(f,fullfile(cfg.figures_dir,'fig_rd_five_codecs.png'),220); close(f);

% --- Fixed condition
f=figure('Visible','off','Position',[50 50 1050 620]);
if exist('tiledlayout','file')==2
    tiledlayout(2,1,'TileSpacing','compact'); nexttile;
else
    subplot(2,1,1);
end
hold on; grid on;
for i=1:numel(algos)
    z=T1(strcmp(T1.algo,algos{i}),:);
    plot(1:height(z),z.CR,['-' marks{i}],'LineWidth',1.2,'DisplayName',algos{i});
end
set(gca,'YScale','log'); ylabel('Compression ratio');
ix=strcmp(T1.algo,algos{1}); xticks(1:sum(ix)); xticklabels(T1.signal(ix)); xtickangle(25); legend('Location','bestoutside');
title('Fixed operating points across ten disturbance types');
if exist('tiledlayout','file')==2, nexttile; else, subplot(2,1,2); end
hold on; grid on;
for i=1:numel(algos)
    z=T1(strcmp(T1.algo,algos{i}),:);
    plot(1:height(z),min(z.SNR_dB,100),['-' marks{i}],'LineWidth',1.2,'DisplayName',algos{i});
end
ylabel('Reconstruction SNR (dB, capped at 100)');
xticks(1:sum(ix)); xticklabels(T1.signal(ix)); xtickangle(25);
save_figure_compat(f,fullfile(cfg.figures_dir,'fig_fixed_scenarios.png'),220); close(f);

% --- SVDCS duration scaling
sel=strcmp(T3.signal,'sparse'); z=T3(sel,:);
cyc=cellfun(@(s) sscanf(s,'cycles=%d'),z.param);
f=figure('Visible','off','Position',[100 100 760 520]); loglog(cyc,z.CR,'o-','LineWidth',1.4); grid on;
xlabel('Record length (fundamental cycles)'); ylabel('Compression ratio');
title('SVDCS duration scaling under sparse events');
save_figure_compat(f,fullfile(cfg.figures_dir,'fig_svdcs_duration.png'),220); close(f);

% --- SVDCS long record rate-distortion
z=T3(strcmp(T3.signal,'periodic_complex_200'),:);
f=figure('Visible','off','Position',[100 100 760 520]); plot(z.rate_bps,z.NMSE_dB,'o-','LineWidth',1.4); grid on;
xlabel('Rate R (compressed bits/sample)'); ylabel('NMSE (dB)');
title('SVDCS rate-distortion on 200-cycle repeated complex record');
save_figure_compat(f,fullfile(cfg.figures_dir,'fig_svdcs_long_rd.png'),220); close(f);

if nargin>=5 && ~isempty(T4)
    f=figure('Visible','off','Position',[100 100 820 540]); hold on; grid on;
    for i=1:numel(algos)
        z=T4(strcmp(T4.algo,algos{i}),:);
        labels=cellfun(@parse_snr_label,z.param);
        plot(labels,z.SNR_dB,['-' marks{i}],'LineWidth',1.3,'DisplayName',algos{i});
    end
    xlabel('Input SNR (dB; Inf displayed at 60)'); ylabel('Reconstruction SNR (dB)');
    title('Noise robustness'); legend('Location','best');
    save_figure_compat(f,fullfile(cfg.figures_dir,'fig_noise_robustness.png'),220); close(f);
end
end

function v=parse_snr_label(s)
s=strrep(s,'inputSNR=','');
v=str2double(s);
if isinf(v)||isnan(v),v=60;end
end
