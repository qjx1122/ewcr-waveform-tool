function xhat = mmc_decode_signal(packets, cfg)
%MMC_DECODE_SIGNAL Independent decoder for packet structures produced by
% mmc_encode_signal. This reconstructs only from quantized model parameters,
% residual payload, and decoder state.

if nargin < 2 || isempty(cfg), cfg = mmc_default_config(); end
W = numel(packets); N = cfg.N;
xhat = zeros(W*N,1);
state.prev2 = zeros(N,1);
state.prev1 = zeros(N,1);
state.prev_model = [];
state.history_count = 0;

for w = 1:W
    p = packets{w};
    xn_model = mmc_reconstruct_model(p.model, p.theta_q, state, cfg);
    rn = mmc_decode_residual(p.residual, cfg);
    xw = (xn_model + rn) * 2^(p.kx);
    idx = (w-1)*N + (1:N);
    xhat(idx) = xw;

    state.prev2 = state.prev1;
    state.prev1 = xw;
    state.prev_model = p.model;
    state.prev_model.theta_q = p.theta_q;
    state.history_count = state.history_count + 1;
end
end
