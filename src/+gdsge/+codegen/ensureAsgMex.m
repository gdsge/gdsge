function ensureAsgMex()
% ENSUREASGMEX  Compile src/kernels/asg_mex.cpp when needed (hash cache
%   asg_mex.cache beside it, same mechanism as mex_<model>.cache). Needed
%   before ASG codegen: emitCompile reads asg.get_mex_constants(). Deliberate
%   improvement over the old toolbox's ship-prebuilt + manual compile_asg.m;
%   invisible to callers. Compile flags mirror compile_asg.m (minus -v).
here      = fileparts(mfilename('fullpath'));   % src/+gdsge/+codegen
srcRoot   = fileparts(fileparts(here));         % src
repoRoot  = fileparts(srcRoot);
kernels   = fullfile(srcRoot, 'kernels');
includeDir = fullfile(repoRoot, 'include');
cppFile   = fullfile(kernels, 'asg_mex.cpp');
cacheFile = fullfile(kernels, 'asg_mex.cache');
mexFile   = fullfile(kernels, ['asg_mex.' mexext]);
% Cache key = the .cpp PLUS the asg.h header it includes, so edits to either
% trigger a rebuild (needsCompile otherwise watches only the .cpp text).
cacheText = [fileread(cppFile) newline fileread(fullfile(includeDir, 'asg.h'))];
if exist(mexFile, 'file') == 3 && ~gdsge.codegen.needsCompile(cacheText, cacheFile)
    return;
end
fprintf('Compiling asg_mex (cache-gated):\n');
oldCd = pwd; restore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(kernels);
% -DUSE_OMP activates the #ifdef USE_OMP parallel-for loops in asg.h (the eval
% inside push_to_info/at_valid + the surplus loops); /openmp alone is inert
% because those pragmas are guarded. Disjoint-write loops -> results unchanged.
mexArgs = gdsge.codegen.mexOmpArgs();   % per-platform OpenMP compile+link flags
mex('-DUSE_OMP', ['-I' includeDir], 'asg_mex.cpp', mexArgs{:});
gdsge.codegen.writeText(cacheFile, cacheText);
end
