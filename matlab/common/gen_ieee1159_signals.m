function [x, info] = gen_ieee1159_signals(kind, fs, cycles, params)
%GEN_IEEE1159_SIGNALS  Generate synthetic power-quality disturbance
%waveforms following the disturbance models of IEEE Std 1159 and the
%power-quality compression literature (Tcheou et al., IEEE TSG 2014).
%
%   x = GEN_IEEE1159_SIGNALS(kind, fs, cycles) generates an Nx1 vector.
%   kind: 'pure' | 'sag' | 'swell' | 'interruption' | 'harmonics' |
%         'osc_transient' | 'notch' | 'flicker' | 'spike' | 'complex'
%   fs   : sampling frequency in Hz (e.g. 12800)
%   cycles: total number of 50 Hz cycles in the record
%   params: optional struct overriding the default disturbance parameters
%           (see below). Returns struct info with parameter values used.
%
%   The fundamental amplitude is 1 p.u.; nominal frequency F0 = 50 Hz.
if nargin < 3 || isempty(cycles), cycles = 10; end
if nargin < 4 || isempty(params), params = struct(); end
if iscell(kind), kind = kind{1}; end
kind = char(kind);
if ~(isscalar(fs) && isnumeric(fs) && isfinite(fs) && fs>0)
    error('gen_ieee1159_signals:badFs','fs must be a positive finite scalar.');
end
if ~(isscalar(cycles) && isnumeric(cycles) && isfinite(cycles) && cycles>0)
    error('gen_ieee1159_signals:badCycles','cycles must be a positive finite scalar.');
end
if ~isstruct(params), error('params must be a struct.'); end

F0   = 50;                       % system frequency (Hz)
N    = round(fs * cycles / F0);  % total samples
n    = (0:N-1).';
w    = 2*pi*F0/fs;

% ---- defaults (IEEE 1159 ranges) --------------------------------------
d = struct('sag_mag',     0.7,    ...  % sag magnitude (p.u.)
           'sag_start',   0.2,    ...  % start fraction of record
           'sag_dur',     0.3,    ...  % duration fraction of record
           'swell_mag',   1.3,    ...
           'swell_start', 0.2,    ...
           'swell_dur',   0.3,    ...
           'harm_a',      [0.05 0.03 0.02 0.01], ... % harmonic amps
           'harm_h',      [3 5 7 11],               ... % harmonic orders
           'trans_amp',   0.6,    ...  % oscillatory transient amplitude
           'trans_f',     550,    ...  % oscillatory transient freq (Hz)
           'trans_start', 0.35,   ...
           'notch_depth', 0.25,   ...
           'notch_frac',  0.05,   ...  % notch duty cycle
           'flicker_m',   0.05,   ...  % flicker modulation depth
           'flicker_f',   8,      ...  % flicker modulation freq (Hz)
           'spike_amp',   1.2,    ...
           'spike_n',     3,      ...
           'noise_snr',   inf);
f = fieldnames(d);
for i = 1:numel(f)
    if isfield(params, f{i})
        d.(f{i}) = params.(f{i});
    end
end
info = d; info.fs = fs; info.F0 = F0; info.N = N;

base = sin(w*n);   % pure sinusoid

switch kind
    case 'pure'
        x = base;

    case 'sag'
        seg = make_seg(N, d.sag_start, d.sag_dur);
        x = base .* (1 - (1-d.sag_mag)*seg);

    case 'swell'
        seg = make_seg(N, d.swell_start, d.swell_dur);
        x = base .* (1 + (d.swell_mag-1)*seg);

    case 'interruption'
        seg = make_seg(N, 0.2, 0.3);
        x = base .* (1 - seg);           % magnitude drops to ~0

    case 'harmonics'
        x = base;
        for i = 1:numel(d.harm_a)
            x = x + d.harm_a(i)*sin(d.harm_h(i)*w*n + 0.3*i);
        end
        x = x / (1 + sum(d.harm_a));     % keep RMS ~1 p.u.

    case 'osc_transient'
        x = base;
        seg = make_seg(N, d.trans_start, 0.05);
        env = exp(-(n - d.trans_start*N) / (0.02*N));   % exponential decay
        x = x + d.trans_amp*sin(2*pi*d.trans_f/fs*n).*env.*seg;

    case 'notch'
        x = base;
        period = round(fs/F0);                 % samples per cycle
        notch_w = max(1, round(d.notch_frac*period));
        idx = false(N,1);
        peak_pos = round(period*(0.25:1:cycles-0.75));  % peaks of sin at 1/4 cycle
        for i = 1:numel(peak_pos)
            p = peak_pos(i);
            if p > 0 && p <= N
                lo = max(1, p-notch_w); hi = min(N, p+notch_w);
                idx(lo:hi) = true;
            end
        end
        x(idx) = (1 - d.notch_depth)*x(idx);

    case 'flicker'
        x = base .* (1 + d.flicker_m*sin(2*pi*d.flicker_f/fs*n));

    case 'spike'
        x = base;
        rng(7);  % deterministic
        spike_pos = round(linspace(0.1*N, 0.9*N, d.spike_n)) + randi(5, d.spike_n, 1);
        for i = 1:d.spike_n
            p = min(max(spike_pos(i), 2), N-1);
            x(p) = d.spike_amp * sign(x(p));
        end

    case 'complex'   % sag + harmonics + oscillatory transient
        x = base;
        seg = make_seg(N, 0.15, 0.4);
        x = x .* (1 - (1-0.75)*seg);
        for i = 1:numel(d.harm_a)
            x = x + d.harm_a(i)*sin(d.harm_h(i)*w*n + 0.3*i);
        end
        x = x/(1 + sum(d.harm_a));
        env = exp(-(n - 0.55*N) / (0.01*N));
        x = x + 0.5*sin(2*pi*700/fs*n).*env.*make_seg(N, 0.55, 0.02);

    otherwise
        error('gen_ieee1159_signals:unknownKind', 'Unknown kind "%s".', kind);
end

% noise (optional)
if isfinite(d.noise_snr)
    rng(1);
    pn = sum(x.^2)/N;
    x = x + sqrt(pn/(10^(d.noise_snr/10)))*randn(N,1);
end

x = x(:);
info.kind = kind;
end

% ----------------------------------------------------------------------
function seg = make_seg(N, start_frac, dur_frac)
    a = max(1, round(start_frac*N));
    b = min(N, a + round(dur_frac*N));
    seg = zeros(N,1);
    seg(a:b) = 1;
end
