function writeText(p, txt)
% WRITETEXT  Write text to a file byte-for-byte (fwrite: no newline translation).
%   Shared by the codegen driver and both generators so cache comparisons
%   (needsCompile) see identical bytes on every platform.
fid = fopen(p, 'w');
if fid < 0
    error('gdsge:codegen:cannotWrite', 'Cannot write %s', p);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, txt);
end
