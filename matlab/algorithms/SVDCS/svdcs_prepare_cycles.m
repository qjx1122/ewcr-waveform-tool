function [frames, info] = svdcs_prepare_cycles(x, cfg)
%SVDCS_PREPARE_CYCLES Synchronous one-cycle segmentation for SVDCS.
%
% Output:
%   frames : Nframes x Nppc, each row is exactly one fundamental cycle.
%
% Modes:
%   nominal    - nominal-cycle timing, locally Lagrange-resampled
%   zero_cross - positive-going zero-crossing cycle boundaries, then
%                locally Lagrange-resampled
%
% The paper uses a zero-crossing frequency estimator and Lagrange
% interpolation. Because the exact estimator from cited Ref. [26] is not
% specified in this paper, zero_cross mode is a transparent approximation.

x = double(x(:));
fs = cfg.fs;
f0 = cfg.f0;
M  = cfg.Nppc;

switch lower(cfg.sync_mode)
    case 'nominal'
        T = fs / f0; % samples per nominal fundamental cycle
        % Keep every complete cycle. The original v1 reproduction used
        % floor((N-1)/T), which drops one valid frame when N/T is integer.
        nFrames = floor(numel(x) / T);
        while nFrames > 0
            lastStart = 1 + (nFrames-1)*T;
            lastQuery = lastStart + (M-1)*(T/M);
            if lastQuery <= numel(x) + 10*eps(numel(x)), break; end
            nFrames = nFrames - 1;
        end
        if nFrames < 1
            error('Input is shorter than one complete nominal cycle.');
        end

        frames = zeros(nFrames, M);
        starts = 1 + (0:nFrames-1) * T;

        for k = 1:nFrames
            tq = starts(k) + (0:M-1) * (T/M);
            frames(k,:) = svdcs_lagrange_interp(x, tq, cfg.lagrange_order).';
        end

        info.boundaries = [starts, starts(end)+T];
        info.mode = 'nominal';
        info.nFrames = nFrames;

    case 'zero_cross'
        % Positive-going zero crossings, linearly localized.
        idx = find(x(1:end-1) <= 0 & x(2:end) > 0);
        if isempty(idx)
            error('No positive-going zero crossings found.');
        end

        zc = zeros(size(idx));
        for k = 1:numel(idx)
            i = idx(k);
            den = x(i+1)-x(i);
            if abs(den) < eps
                frac = 0;
            else
                frac = -x(i)/den;
            end
            zc(k) = i + frac;
        end

        % Remove implausibly close crossings caused by notches/noise.
        Tnom = fs/f0;
        keep = true(size(zc));
        last = zc(1);
        for k = 2:numel(zc)
            if zc(k)-last < 0.55*Tnom
                keep(k) = false;
            else
                last = zc(k);
            end
        end
        zc = zc(keep);

        % Reject grossly implausible cycle lengths.
        d = diff(zc);
        good = (d > 0.55*Tnom) & (d < 1.60*Tnom);
        pairIdx = find(good);
        if isempty(pairIdx)
            error('No plausible zero-crossing cycle intervals remain.');
        end

        frames = zeros(numel(pairIdx), M);
        bounds = zeros(numel(pairIdx),2);
        for kk = 1:numel(pairIdx)
            k = pairIdx(kk);
            a = zc(k);
            b = zc(k+1);
            tq = a + (0:M-1)*(b-a)/M;
            frames(kk,:) = svdcs_lagrange_interp(x, tq, cfg.lagrange_order).';
            bounds(kk,:) = [a b];
        end

        info.boundaries = bounds;
        info.mode = 'zero_cross';
        info.nFrames = size(frames,1);

    otherwise
        error('Unknown cfg.sync_mode: %s', cfg.sync_mode);
end
end
