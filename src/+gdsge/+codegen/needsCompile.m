function tf = needsCompile(cppText, cacheFile)
% NEEDSCOMPILE  True when cppText differs from the cached copy (or no cache).
%   Replicates the old gdsge_codegen cache semantics: a text-identical
%   mex_<model>.cache means the MEX binary is already up to date.
tf = true;
if exist(cacheFile, 'file')
    tf = ~strcmp(fileread(cacheFile), cppText);
end
end
