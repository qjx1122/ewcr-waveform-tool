function ok = preflight_check()
%PREFLIGHT_CHECK Check paths/toolbox functions before running benchmark.
fprintf('FiveCodecBenchmark preflight\n');
fprintf('MATLAB: %s\n',version);
fprintf('Computer: %s (%s)\n',computer,computer('arch'));
req = {'wavedec','waverec','huffmandict','huffmanenco','dftmtx','dctmtx'};
ok=true;
for k=1:numel(req)
    p=which(req{k});
    if isempty(p)
        fprintf('  MISSING: %-14s\n',req{k}); ok=false;
    else
        fprintf('  OK:      %-14s %s\n',req{k},p);
    end
end
sp=which('string','-all');
if ischar(sp), sp={sp}; end
if isempty(sp)
    fprintf('  Note: string() not found; package v2.6 does not require it.\n');
elseif numel(sp)>1
    fprintf('  Note: multiple string implementations found. v2.6 avoids string() calls.\n');
    for k=1:numel(sp), fprintf('         %s\n',sp{k}); end
end
if ok
    fprintf('Preflight PASSED. Run: run_all_experiments\n');
else
    fprintf('Preflight FAILED: install/enable the missing required toolbox functions.\n');
end
end
