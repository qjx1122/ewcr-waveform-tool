function x = make_long_pattern(fs,cyclesTotal,eventDensity)
%MAKE_LONG_PATTERN Long 50-Hz signal with controllable event density.
% eventDensity: 'sparse'|'medium'|'dense'|'periodic_complex'.
F0=50; one=round(fs/F0);
base=sin(2*pi*F0/fs*(0:one-1)').';
base=base(:);
x=repmat(base,cyclesTotal,1);

switch lower(eventDensity)
    case 'periodic_complex'
        pat=gen_ieee1159_signals('complex',fs,10);
        rep=ceil(numel(x)/numel(pat)); x=repmat(pat,rep,1); x=x(1:cyclesTotal*one);
    otherwise
        switch lower(eventDensity)
            case 'sparse', gap=100;
            case 'medium', gap=30;
            case 'dense', gap=10;
            otherwise, error('Unknown event density.');
        end
        p=gap;
        while p+10 <= cyclesTotal
            seg=gen_ieee1159_signals('complex',fs,10);
            idx=(p*one+1):((p+10)*one);
            x(idx)=seg;
            p=p+gap;
        end
end
end
