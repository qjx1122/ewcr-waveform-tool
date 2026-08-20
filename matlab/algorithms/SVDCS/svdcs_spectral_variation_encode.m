function mc = svdcs_spectral_variation_encode(A, gammaCode)
%SVDCS_SPECTRAL_VARIATION_ENCODE Core spectral-variation compressor.
%
% Direct implementation of Fig. 2 / Sec. II-E:
%   MAX(q), MIN(q), SUM(q), CNT(q)
%   novelty if |MAX(q)-MIN(q)| > gamma_q
%   on novelty, output [q, SUM(q), CNT(q)] and reset that q-position.
%
% Inputs
%   A         : Nframes x Q integer FFT_ARRAY values
%   gammaCode : 1 x Q tolerances in the same integer-code units
%
% Output
%   mc.q   : uint16 q-position index (1-based MATLAB)
%   mc.sum : int64 accumulated FFT-bin-part code
%   mc.cnt : uint32 number of frames represented by this run
%
% The final active run for every q is flushed at end-of-stream. This flush
% is required for finite-record reconstruction.

[nFrames,Q] = size(A);
if nFrames < 1
    error('A must contain at least one frame.');
end
if numel(gammaCode) ~= Q
    error('gammaCode length must match A columns.');
end
if Q > double(intmax('uint16'))
    error('Q exceeds uint16 q-index capacity.');
end

A = int64(A);
gammaCode = double(gammaCode(:).');

MAXV = A(1,:);
MINV = A(1,:);
SUMV = A(1,:);
CNTV = ones(1,Q,'uint32');

% Moderate initial allocation. Long recordings can still grow as required.
cap = max(4*Q, min(100000, max(Q, round(0.001*nFrames*Q))));
qOut   = zeros(1,cap,'uint16');
sumOut = zeros(1,cap,'int64');
cntOut = zeros(1,cap,'uint32');
nOut = 0;

for i = 2:nFrames
    NEW = A(i,:);

    newMax = max(MAXV,NEW);
    newMin = min(MINV,NEW);
    novelty = double(abs(newMax-newMin)) > gammaCode;

    idx = find(novelty);
    nNew = numel(idx);

    if nNew > 0
        needed = nOut+nNew;
        if needed > numel(qOut)
            grow = max(needed-numel(qOut), max(Q,ceil(0.5*numel(qOut))));
            qOut(end+grow) = uint16(0);
            sumOut(end+grow) = int64(0);
            cntOut(end+grow) = uint32(0);
        end

        rr = nOut+(1:nNew);
        qOut(rr) = uint16(idx);
        sumOut(rr) = SUMV(idx);
        cntOut(rr) = CNTV(idx);
        nOut = needed;
    end

    % No novelty -> update running range, sum and count.
    keep = ~novelty;
    MAXV(keep) = newMax(keep);
    MINV(keep) = newMin(keep);
    SUMV(keep) = SUMV(keep) + NEW(keep);
    CNTV(keep) = CNTV(keep) + uint32(1);

    % Novelty -> start a new run from current frame.
    MAXV(novelty) = NEW(novelty);
    MINV(novelty) = NEW(novelty);
    SUMV(novelty) = NEW(novelty);
    CNTV(novelty) = uint32(1);
end

% End-of-stream flush for all q-positions.
needed=nOut+Q;
if needed > numel(qOut)
    grow=needed-numel(qOut);
    qOut(end+grow)=uint16(0);
    sumOut(end+grow)=int64(0);
    cntOut(end+grow)=uint32(0);
end
rr=nOut+(1:Q);
qOut(rr)=uint16(1:Q);
sumOut(rr)=SUMV;
cntOut(rr)=CNTV;
nOut=needed;

mc.q   = qOut(1:nOut);
mc.sum = sumOut(1:nOut);
mc.cnt = cntOut(1:nOut);
mc.Q = Q;
mc.Nframes = nFrames;
end
