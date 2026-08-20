function m = unified_metrics(x,xhat,bitsCompressed,fs,f0,rawBits)
%UNIFIED_METRICS Rate, waveform distortion and power-signal feature errors.
if nargin<6 || isempty(rawBits), rawBits=16; end
x=x(:); xhat=xhat(:);
N=min(numel(x),numel(xhat)); x=x(1:N); xhat=xhat(1:N);
e=x-xhat;
Es=sum(x.^2); Ee=sum(e.^2);
m.N=N;
m.bits_original=N*rawBits;
m.bits_compressed=bitsCompressed;
m.CR=m.bits_original/max(bitsCompressed,1);
m.rate_bps=bitsCompressed/max(N,1);
m.NMSE_lin=Ee/max(Es,eps);
m.NMSE_dB=10*log10(max(m.NMSE_lin,realmin));
m.SNR_dB=-m.NMSE_dB;
m.RMSE=sqrt(mean(e.^2));
m.PRD=100*sqrt(Ee/max(Es,eps));
m.MAXE=max(abs(e));
peak=max(abs(x));
m.PSNR_dB=20*log10(max(peak,eps)/max(m.RMSE,eps));

% Fundamental amplitude/phase by direct least-squares projection.
t=(0:N-1)'/fs;
A=[cos(2*pi*f0*t), sin(2*pi*f0*t)];
co=A\x; ch=A\xhat;
amp=@(c) hypot(c(1),c(2));
ph=@(c) atan2(-c(2),c(1));
m.fund_amp=amp(co); m.fund_amp_hat=amp(ch);
m.fund_amp_relerr=abs(m.fund_amp_hat-m.fund_amp)/max(m.fund_amp,eps);
dp=angle(exp(1i*(ph(ch)-ph(co))));
m.fund_phase_err_deg=abs(dp)*180/pi;

% Harmonic distortion estimate from integer-order projections.
Hmax=min(25,floor((fs/2)/f0));
Ah=zeros(Hmax,1); Ahh=zeros(Hmax,1);
for h=1:Hmax
    B=[cos(2*pi*h*f0*t),sin(2*pi*h*f0*t)];
    c=B\x; d=B\xhat;
    Ah(h)=hypot(c(1),c(2)); Ahh(h)=hypot(d(1),d(2));
end
m.THD=sqrt(sum(Ah(2:end).^2))/max(Ah(1),eps);
m.THD_hat=sqrt(sum(Ahh(2:end).^2))/max(Ahh(1),eps);
m.THD_abs_error=abs(m.THD_hat-m.THD);
end
