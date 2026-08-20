function cfg = benchmark_config()
%BENCHMARK_CONFIG Unified experiment settings for the five-codec benchmark.

cfg.fs = 12800;
cfg.f0 = 50;
cfg.cycles = 10;
cfg.raw_bits_per_sample = 16;
cfg.seed = 42;
cfg.results_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');
cfg.figures_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'figures');
cfg.run_noise = true;
cfg.quick_mode = false;

% Fixed operating points.
cfg.ASBC = struct('W',20,'Rk',8,'block_cycles',5,'act_thresh',1e-4, ...
    'ke_max',48,'energy_frac',0.99,'RQ',16);
cfg.DWT = struct('qbits',8,'thr_mode','energy','thr_frac',0.99,'delta',true,'RQ',16);
cfg.CS = struct('M_ratio',0.20,'basis','dft','K0',48,'tol',1e-6,'seed',42,'ybits',16,'RQ',16);

cfg.MMC = mmc_default_config();
cfg.MMC.fs = cfg.fs;
cfg.MMC.fn = cfg.f0;
cfg.MMC.N = round(cfg.fs/cfg.f0);     % one synchronized cycle
cfg.MMC.n_tot = cfg.MMC.N;            % 1 bit/sample budget
cfg.MMC.raw_bits_per_sample = 16;
cfg.MMC.verbose = false;

cfg.SVDCS = svdcs_default_config();
cfg.SVDCS.fs = cfg.fs;
cfg.SVDCS.f0 = cfg.f0;
cfg.SVDCS.Nppc = round(cfg.fs/cfg.f0);
cfg.SVDCS.sync_mode = 'nominal';
cfg.SVDCS.G = 0.03;
cfg.SVDCS.beta = 0.10;
cfg.SVDCS.signal_scale = 1;
cfg.SVDCS.verbose = false;
% Which SVDCS rate is used as the unified CR: 'payload' or 'engineering'.
cfg.svdcs_rate_mode = 'payload';

% Rate-distortion sweeps.
cfg.rd.asbc_W = [5 20 100];
cfg.rd.asbc_Rk = [4 6 8 10 12];
cfg.rd.dwt_energy = [0.95 0.97 0.99 0.999 0.9999];
cfg.rd.cs_ratio = [0.05 0.10 0.20 0.30 0.40];
cfg.rd.mmc_bits = [128 192 256 384 512];
cfg.rd.svdcs_G = [0.005 0.01 0.02 0.03 0.05 0.08 0.12];

cfg.duration_cycles = [10 50 200 1000 5000];
cfg.noise_snr_db = [Inf 45 35 25];

% Runtime / edge-resource profiling.
cfg.profile.warmup = 1;       % discarded warm-up repetitions
cfg.profile.repeats = 3;      % steady-state repetitions; increase to 10+ for publication
cfg.profile.scenarios = {'pure','complex','osc_transient'};
cfg.profile.primary_scenario = 'complex';
cfg.profile.run_scaling = true;
cfg.profile.scaling_cycles = [2 4 10]; % keep CS dense-DFT RAM manageable

end
