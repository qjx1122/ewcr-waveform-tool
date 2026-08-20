function MR = svdcs_reconstruct_MR(mc)
%SVDCS_RECONSTRUCT_MR Decode MC into the reconstruction matrix MR.
%
% For each q-position, every MC column contains SUM and CNT for one run.
% The reconstructed value for that run is SUM/CNT and is repeated CNT times.
% Each q-position is decoded independently until Nframes values are obtained.

Q = mc.Q;
Nframes = mc.Nframes;
MR = zeros(Nframes,Q);

pos = ones(1,Q); % next frame position to fill for each q

for e = 1:numel(mc.q)
    q = double(mc.q(e));
    c = double(mc.cnt(e));
    if c < 1
        error('Invalid CNT=0 in MC.');
    end

    a = pos(q);
    b = a + c - 1;
    if b > Nframes
        error('MC decoding overflow at event %d, q=%d.', e, q);
    end

    meanValue = double(mc.sum(e)) / c;
    MR(a:b,q) = meanValue;
    pos(q) = b + 1;
end

if any(pos ~= Nframes+1)
    bad = find(pos ~= Nframes+1,1);
    error('MC does not fully cover q=%d: filled %d of %d frames.', ...
        bad, pos(bad)-1, Nframes);
end
end
