% +GDSGE/+CODEGEN  Code generators consuming the IR.
%   codegen        - Public driver: <model>.gmod in cwd -> IR json + generated files + MEX
%   generateMatlab - IR -> iter_<model>.m / simulate_<model>.m
%   generateCxx    - IR -> mex_<model>.cpp / compile_<model>.m (adept autodiff)
%   needsCompile   - True when cpp text differs from mex_<model>.cache (or no cache)
%   writeText      - Shared byte-exact text-file writer
%   codeWriter     - Indented-line buffer shared by all codegen backends
%   +mat           - MATLAB-emitting backend (section emitters)
%   +cxx           - C++ adept-autodiff backend (section emitters)
