function ensurePath()
% ENSUREPATH  Make the flat runtime kernels in src/kernels loadable.
%   (1) Register src/kernels on the MATLAB path so the flat kernel names
%       (interp_construct_mex, interp_eval_mex, asg_mex, gen_discrete_markov_rn)
%       resolve. The folder is located relative to THIS file, so a single
%       `addpath src` is enough for callers — this hook adds the rest.
%   (2) On Windows, append that directory to the system PATH so the MEX
%       kernels can load essential_blas.dll.
%   Idempotent. Called first by gdsge_codegen and by every generated
%   iter_<model>/simulate_<model>, ahead of any kernel call.
kernelsDir = fullfile( ...
    fileparts(fileparts(fileparts(mfilename('fullpath')))), 'kernels');

if ~any(strcmp(strsplit(path, pathsep), kernelsDir))
    addpath(kernelsDir);
end

if ~ispc; return; end
if ~any(strcmp(strsplit(getenv('PATH'), ';'), kernelsDir))
    setenv('PATH', [getenv('PATH') ';' kernelsDir]);
end
end
