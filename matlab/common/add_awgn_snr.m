function y = add_awgn_snr(x,snrdb,seed)
%ADD_AWGN_SNR Toolbox-free deterministic AWGN injection.
if isinf(snrdb), y=x; return; end
rng(seed); x=x(:); p=mean(x.^2); np=p/10^(snrdb/10);
y=x+sqrt(np)*randn(size(x));
end
