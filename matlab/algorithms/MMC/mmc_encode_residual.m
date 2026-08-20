function residual = mmc_encode_residual(r, method, budget, cfg)
%MMC_ENCODE_RESIDUAL Fixed-budget transform residual codec.
% Paper-faithful transform competition; engineering simplification of the
% Antonini/Khan arithmetic bit-plane coder in the released Python code.

r=r(:); method=upper(method);
residual=struct('name',method,'nbits',0,'kr',0,'payload',[],'rhat',zeros(size(r)), ...
    'transform','','value_bits',cfg.residual_value_bits);
if strcmp(method,'NONE') || budget<=0
    residual.name='none'; return;
end

switch method
    case 'DCT'
        T=dct_matrix(cfg.N); c=T*r; invfun=@(z) T'*z;
    case 'DWT'
        c=haar_forward(r); invfun=@(z) haar_inverse(z);
    otherwise
        error('Unknown residual method %s',method);
end

% Power-of-two coefficient normalization, analogous to public implementation.
kr=max(0,ceil(log2(max(abs(c))+1e-12)));
kr=min(2^cfg.n_kr-1,kr);
cn=c*2^(-kr);

% Sparse fixed-budget payload: [index, signed quantized coefficient].
% It is deterministic and has explicit real bit accounting.
index_bits=ceil(log2(cfg.N));
value_bits=cfg.residual_value_bits;
bits_per_coeff=index_bits+value_bits;
K=floor(budget/bits_per_coeff);
K=min(cfg.N,max(0,K));
if K<cfg.residual_min_coeffs
    residual.name=method; residual.transform=method; residual.kr=kr; return;
end

[~,ord]=sort(abs(cn),'descend'); sel=sort(ord(1:K));
qmax=2^(value_bits-1)-1;
q=round(cn(sel)*qmax);
q=max(-qmax,min(qmax,q));
chat=zeros(cfg.N,1); chat(sel)=q/qmax;
c_hat=chat*2^(kr);
rhat=invfun(c_hat);

residual.name=method;
residual.transform=method;
residual.kr=kr;
residual.nbits=K*bits_per_coeff;
residual.payload.indices=sel;
residual.payload.q=q;
residual.payload.index_bits=index_bits;
residual.payload.value_bits=value_bits;
residual.rhat=rhat;
end

function T=dct_matrix(N)
persistent cacheN cacheT
if isempty(cacheN)||cacheN~=N
    T=zeros(N,N);
    for k=0:N-1
        a=sqrt(2/N); if k==0,a=sqrt(1/N);end
        T(k+1,:)=a*cos(pi*(0.5+(0:N-1))*k/N);
    end
    cacheN=N; cacheT=T;
else
    T=cacheT;
end
end

function c=haar_forward(x)
x=x(:); N=numel(x); c=zeros(N,1); temp=x; len=N; pos=N;
while len>1
    a=(temp(1:2:len)+temp(2:2:len))/sqrt(2);
    d=(temp(1:2:len)-temp(2:2:len))/sqrt(2);
    c(pos-numel(d)+1:pos)=d;
    pos=pos-numel(d); temp(1:numel(a))=a; len=numel(a);
end
c(1)=temp(1);
end

function x=haar_inverse(c)
c=c(:); N=numel(c); a=c(1); pos=2; len=1;
while len<N
    d=c(pos:pos+len-1); pos=pos+len;
    temp=zeros(2*len,1);
    temp(1:2:end)=(a+d)/sqrt(2);
    temp(2:2:end)=(a-d)/sqrt(2);
    a=temp; len=2*len;
end
x=a;
end
