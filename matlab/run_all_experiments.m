function run_all_experiments()
%RUN_ALL_EXPERIMENTS Entry point for five-codec benchmark + resource audit.
ROOT=fileparts(mfilename('fullpath'));
addpath(genpath(ROOT));
if ~preflight_check()
    error('FiveCodecBenchmark:Preflight','Preflight check failed.');
end
cfg=benchmark_config();
if ~exist(cfg.results_dir,'dir'),mkdir(cfg.results_dir);end
if ~exist(cfg.figures_dir,'dir'),mkdir(cfg.figures_dir);end
rng(cfg.seed);

fprintf('\n=== Experiment 1: fixed scenarios ===\n');
T1=exp1_fixed_scenarios(cfg);
fprintf('\n=== Experiment 2: rate-distortion ===\n');
T2=exp2_rate_distortion(cfg);
fprintf('\n=== Experiment 3: SVDCS duration / event density ===\n');
T3=exp3_svdcs_duration(cfg);
T4=[];
if cfg.run_noise
    fprintf('\n=== Experiment 4: noise robustness ===\n');
    T4=exp4_noise_robustness(cfg);
end
fprintf('\n=== Experiment 5: encode/decode timing + memory/resource profiling ===\n');
[T5,T5scale,Env]=exp5_runtime_memory(cfg);
fprintf('\n=== Experiment 6: compressed-payload manifest ===\n');
T6=exp6_payload_manifest(cfg);

make_all_figures(cfg,T1,T2,T3,T4);
make_resource_figures(cfg,T5);
save(fullfile(cfg.results_dir,'all_experiments.mat'),'cfg','T1','T2','T3','T4','T5','T5scale','T6','Env');

fprintf('\nFinished. CSV/MAT results: %s\n',cfg.results_dir);
fprintf('Figures: %s\n',cfg.figures_dir);
fprintf('Resource timing is host-specific; compare medians only within the same MATLAB machine.\n');
end
