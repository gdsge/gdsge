function ensureInterpEvalMex()
% ENSUREINTERPEVALMEX  Compile src/kernels/interp_eval_mex.cpp when needed
%   (hash cache interp_eval_mex.cache beside it, same mechanism as
%   ensureSplineConstructMex). Generic uniform-order evaluator over the
%   double-only core include/interp_eval_double.h — the myppual_mex replacement
%   for MATLAB-side interpolation. The cache key includes interp_eval_double.h
%   so header edits trigger a rebuild. No BLAS link (pure double math).
here       = fileparts(mfilename('fullpath'));   % src/+gdsge/+codegen
srcRoot    = fileparts(fileparts(here));          % src
repoRoot   = fileparts(srcRoot);
kernels    = fullfile(srcRoot, 'kernels');
includeDir = fullfile(repoRoot, 'include');
cppFile    = fullfile(kernels, 'interp_eval_mex.cpp');
cacheFile  = fullfile(kernels, 'interp_eval_mex.cache');
mexFile    = fullfile(kernels, ['interp_eval_mex.' mexext]);
cppText    = [fileread(cppFile), fileread(fullfile(includeDir, 'interp_eval_double.h'))];
if exist(mexFile, 'file') == 3 && ~gdsge.codegen.needsCompile(cppText, cacheFile)
    return;
end
fprintf('Compiling interp_eval_mex (cache-gated):\n');
oldCd = pwd; restore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(kernels);
mexArgs = gdsge.codegen.mexOmpArgs();   % per-platform OpenMP compile+link flags
mex(['-I' includeDir], '-DUSE_OMP', 'interp_eval_mex.cpp', mexArgs{:});
gdsge.codegen.writeText(cacheFile, cppText);
end
