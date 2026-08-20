function xhat = svdcs_decode(packet)
%SVDCS_DECODE Reconstruct synchronized waveform from an SVDCS packet.

mc=packet.mc;

% Reconstruct FFT_ARRAY matrix in quantized-code units.
MRcode=svdcs_reconstruct_MR(mc);
MR=MRcode*packet.fftStep;

N=packet.Nppc;
nF=packet.Nframes;

X=zeros(nF,N);
X(:,:)=MR(:,1:2:end) + 1j*MR(:,2:2:end);

% Encoder used X = FFT(frame)/N, therefore inverse is IFFT(N*X).
framesHat=real(ifft(X*N,[],2));
xhat=reshape(framesHat.',[],1)*packet.signal_scale;
end
