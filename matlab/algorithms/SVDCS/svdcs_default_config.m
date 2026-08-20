function cfg = svdcs_default_config()
%SVDCS_DEFAULT_CONFIG Default configuration for the spectral-variation
% compression system (SVDCS) described by Kapisch et al., IEEE TII 2022.
%
% Core paper parameters:
%   Gamma(q) = G * q^(-beta)
%   G = 0.03, beta = 0.1
%   one fundamental cycle per FFT frame
%   real and imaginary FFT-bin parts processed independently
%   final lossless stage: LZW
%
% IMPORTANT:
% The paper does not publish the exact source code of the zero-crossing
% frequency estimator [26], the FPGA Lagrange filter-bank coefficients, or
% the binary packing convention used before LZW. Those pieces are therefore
% implemented here as transparent engineering equivalents. The spectral
% variation core follows Fig. 2 / Sec. II-E directly.

cfg.fs = 7680;                 % Hz; paper synthetic default: 60*128
cfg.f0 = 60;                   % nominal fundamental frequency, Hz
cfg.Nppc = 128;                % samples per synchronized fundamental cycle

% Synchronization / interpolation
% 'nominal'    : frame by nominal-cycle duration; best for fixed-frequency data
% 'zero_cross' : positive-going zero crossings + local Lagrange interpolation
cfg.sync_mode = 'nominal';
cfg.lagrange_order = 3;        % cubic local Lagrange interpolation

% Signal scaling.
% Supply x in per-unit and leave signal_scale=1, OR set signal_scale to the
% nominal peak value (e.g. 311 V for 220 Vrms) for physical-unit data.
cfg.signal_scale = 1.0;

% ADC-like input quantization (paper used 16-bit)
cfg.input_quantize = true;
cfg.input_bits = 16;
cfg.input_full_scale_pu = 2.0; % signed full-scale magnitude in p.u.

% FFT-bin quantization (paper says real/imag FFT parts were quantized to 16 bit)
cfg.fft_quantize = true;
cfg.fft_bits = 16;
cfg.fft_full_scale_pu = 2.0;   % full-scale for normalized FFT coefficients

% Spectral-variation tolerance: gamma_q = G*q^(-beta)
cfg.G = 0.03;
cfg.beta = 0.10;

% Paper stores full FFT (not only the non-redundant positive-frequency half)
cfg.store_full_fft = true;

% Bit accounting
cfg.original_bits_per_sample = 16;

% LZW implementation in this project:
% fixed 12-bit LZW codewords, dictionary grows from 256 to 4096 entries,
% then stops growing (no reset code). This is lossless and deterministic.
cfg.use_lzw = true;
cfg.lzw_code_bits = 12;

% Diagnostics
cfg.verbose = true;
end
