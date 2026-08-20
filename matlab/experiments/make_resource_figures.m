function make_resource_figures(cfg,T)
%MAKE_RESOURCE_FIGURES Timing/resource figures from exp5 results.
if isempty(T),return;end
sig=char(cfg.profile.primary_scenario);
S=T(strcmp(T.signal,sig),:);
if isempty(S),return;end
wanted={'ASBC','DWT-Hybrid','CS-OMP','MMC','SVDCS'};
[~,ord]=ismember(wanted,S.algo);
ord=ord(ord>0); S=S(ord,:);

figure('Visible','off');
bar([S.enc_median_s S.dec_median_s]);
set(gca,'XTickLabel',S.algo); ylabel('Time (s)');
legend('Encode median','Decode median','Location','best');
title(sprintf('Encode/decode time: %s',sig)); grid on;
save_figure_compat(gcf,fullfile(cfg.figures_dir,'exp5_runtime.png'),180); close(gcf);

figure('Visible','off');
bar([S.encoder_working_est_bytes S.decoder_working_est_bytes]/1024^2);
set(gca,'XTickLabel',S.algo); ylabel('Estimated working memory (MiB)');
legend('Encoder','Decoder','Location','best');
title(sprintf('Current MATLAB implementation working-set estimate: %s',sig)); grid on;
save_figure_compat(gcf,fullfile(cfg.figures_dir,'exp5_working_memory.png'),180); close(gcf);

figure('Visible','off');
bar([S.logical_payload_bytes S.matlab_packet_bytes]/1024);
set(gca,'XTickLabel',S.algo); ylabel('Size (KiB)');
legend('Logical compressed payload','MATLAB in-memory packet','Location','best');
title(sprintf('Payload vs MATLAB container: %s',sig)); grid on;
save_figure_compat(gcf,fullfile(cfg.figures_dir,'exp5_payload_memory.png'),180); close(gcf);
end
