function ensureSplineConstructMex()
% ENSURESPLINECONSTRUCTMEX  Compile src/kernels/interp_construct_mex.cpp when
%   needed (hash cache interp_construct_mex.cache beside it, same mechanism as
%   ensureAsgMex / mex_<model>.cache). The kernel is the fused spline/linear
%   constructor consumed by gdsge.runtime.constructSplines. Flags: /openmp,
%   -DUSE_OMP, /O2; /fp:precise is the MSVC default and MUST stay (it preserved
%   bit-exactness with the original vendored constructor the kernel replaced).
here      = fileparts(mfilename('fullpath'));   % src/+gdsge/+codegen
srcRoot   = fileparts(fileparts(here));         % src
repoRoot  = fileparts(srcRoot);
kernels   = fullfile(srcRoot, 'kernels');
cppFile   = fullfile(kernels, 'interp_construct_mex.cpp');
cacheFile = fullfile(kernels, 'interp_construct_mex.cache');
mexFile   = fullfile(kernels, ['interp_construct_mex.' mexext]);
cppText   = fileread(cppFile);
if exist(mexFile, 'file') == 3 && ~gdsge.codegen.needsCompile(cppText, cacheFile)
    return;
end
fprintf('Compiling interp_construct_mex (cache-gated):\n');
includeDir = fullfile(repoRoot, 'include');
oldCd = pwd; restore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(kernels);
mexArgs = gdsge.codegen.mexOmpArgs();   % per-platform OpenMP compile+link flags
mex(['-I' includeDir], '-DUSE_OMP', 'interp_construct_mex.cpp', mexArgs{:});
gdsge.codegen.writeText(cacheFile, cppText);
end
