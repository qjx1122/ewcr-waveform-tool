function cfg = mmc_default_config()
%MMC_DEFAULT_CONFIG Configuration for the MATLAB engineering reproduction
% of Presvots et al., Multiple-Model Coding Scheme for Electrical Signal Compression.
%
% The architecture follows the paper/released Python implementation:
%   Stage 1: none / sinusoid / Chebyshev polynomial / sample prediction /
%            parameter prediction
%   Stage 2: none / DCT residual / DWT residual
%   Search : rate-constrained model + residual selection per window.
%
% Note: the public Python code uses Antonini/Khan embedded arithmetic coding
% for residual bit-plane coding. This MATLAB v1 keeps the same transform-
% residual competition but uses a deterministic sparse fixed-budget coder.

cfg.fs = 6400;
cfg.fn = 50;
cfg.N  = 128;          % 20 ms at 6.4 kHz
cfg.n_tot = 128;      % bits/window = 1 bit/sample
cfg.raw_bits_per_sample = 16;

% Parameter quantization
cfg.min_bits_theta = 2;
cfg.max_bits_theta = 10;
cfg.nx_step = 2;      % same idea as released exhaustive code (skip every other bit)

% Header accounting (close to public Python implementation)
cfg.n_kx = 5;
cfg.n_kr = 5;
cfg.n_nr = 10;        % residual payload-length / symbol-count field

% Model dictionary
cfg.use_none = true;
cfg.use_sinusoid = true;
cfg.poly_orders = 0:9;
cfg.sample_pred_orders = 1:2;
cfg.sample_pred_eta = 0:1;
cfg.param_pred_factors = [2 10];

% Sinusoidal prior on normalized waveform
cfg.sin_center = [0.75, cfg.fn, 0];
cfg.sin_width  = [0.5, 0.2, 2*pi];

% Polynomial coefficient standard deviations from authors' public main.py
sigma_poly = {
    [0.215], ...
    [0.2150 0.3812], ...
    [0.2409 0.3812 0.4191], ...
    [0.2409 0.2812 0.4192 0.2539], ...
    [0.2258 0.2812 0.3572 0.2539 0.1363], ...
    [0.2258 0.2877 0.3572 0.2310 0.1364 0.0654], ...
    [0.2268 0.2877 0.3611 0.2310 0.1264 0.0655 0.0426], ...
    [0.2268 0.2870 0.3611 0.2316 0.1264 0.0601 0.0427 0.0326], ...
    [0.2268 0.2871 0.3609 0.2316 0.1262 0.0601 0.0383 0.0327 0.0280], ...
    [0.2268 0.2871 0.3609 0.2317 0.1261 0.0593 0.0384 0.0289 0.0281 0.0244]};
cfg.poly_width = cell(size(sigma_poly));
for k = 1:numel(sigma_poly)
    cfg.poly_width{k} = min(1.55, 7*sigma_poly{k});
end

% Residual methods
cfg.residual_methods = {'none','DCT','DWT'};
cfg.residual_value_bits = 10;  % signed magnitude quantizer in sparse coder
cfg.residual_min_coeffs = 1;

% Search controls
cfg.verbose = true;
cfg.max_windows = inf;
end
