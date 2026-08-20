function b = dir_mfile_bytes(folder)
%DIR_MFILE_BYTES Total bytes of .m source files below folder.
% This is useful for source-package comparison only. It is not compiled
% edge-device flash usage.
b = 0;
if ~isfolder(folder), return; end
D = dir(fullfile(folder,'**','*.m'));
if ~isempty(D), b = sum([D.bytes]); end
end
