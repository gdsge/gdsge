function args = mexOmpArgs()
% MEXOMPARGS  Platform-specific mex() flag arguments to compile AND link the
%   OpenMP flat kernels (asg_mex, interp_construct_mex, interp_eval_mex).
%   Splat into the mex call: mex(..., args{:}).
%
%   The kernels call omp_set_num_threads()/omp pragmas (guarded by -DUSE_OMP),
%   so the link step must pull in the OpenMP runtime. COMPFLAGS is Windows-only,
%   so passing only "/openmp" left Linux/Mac builds with an undefined reference
%   to omp_set_num_threads at link time. Mirror the base toolbox's per-platform
%   kernel recipes (source/compile_asg{,_linux,_macOS}.m):
%     Windows (MSVC) : /openmp via COMPFLAGS; MSVC links the OpenMP runtime itself.
%     Linux   (gcc)  : -fopenmp on both compile (CXXFLAGS) and link (LDFLAGS).
%     macOS          : -fopenmp to compile, libomp (-lomp) to link.
if ispc
    args = {'OPTIMFLAGS=/O2 /DNDEBUG', ...
            'COMPFLAGS=$COMPFLAGS /wd4267 /wd4068 /wd4091 /openmp'};
elseif ismac
    % macOS mex compiles with exceptions OFF by default, but the flat_hash_map
    % kernel (asg.h) calls .at(), which can throw -> need -fexceptions.
    % -fpermissive mirrors the base toolbox's compile_asg_macOS.m.
    args = {'CXXFLAGS=$CXXFLAGS -fopenmp -fexceptions -fpermissive', ...
            'CFLAGS=$CFLAGS -fopenmp -fexceptions -fpermissive', ...
            'LDFLAGS=$LDFLAGS -lomp'};
else  % isunix && ~ismac  (Linux)
    args = {'CXXFLAGS=$CXXFLAGS -fopenmp -fexceptions', ...
            'LDFLAGS=$LDFLAGS -fopenmp'};
end
end
