function T = svdcs_threshold_sweep(x, cfg, Gvalues)
%SVDCS_THRESHOLD_SWEEP Sweep G to obtain CR-quality trade-off.

if nargin<3 || isempty(Gvalues)
    Gvalues=[0 0.005 0.01 0.02 0.03 0.05 0.08 0.12];
end

oldVerbose=cfg.verbose;
cfg.verbose=false;

n=numel(Gvalues);
G=zeros(n,1);
CR=zeros(n,1);
CRtotal=zeros(n,1);
NMSE=zeros(n,1);
COR=zeros(n,1);
RTE=zeros(n,1);
SNR=zeros(n,1);
MCcols=zeros(n,1);

for k=1:n
    cfg.G=Gvalues(k);
    o=svdcs_encode(x,cfg);
    G(k)=cfg.G;
    CR(k)=o.stats.CR_payload_only;
    CRtotal(k)=o.stats.CR_engineering_total;
    NMSE(k)=o.metrics.NMSE;
    COR(k)=o.metrics.COR;
    RTE(k)=o.metrics.RTE_percent;
    SNR(k)=o.metrics.SNR_dB;
    MCcols(k)=o.stats.nMC;
end

cfg.verbose=oldVerbose; %#ok<NASGU>
T=table(G,CR,CRtotal,NMSE,COR,RTE,SNR,MCcols);
end
