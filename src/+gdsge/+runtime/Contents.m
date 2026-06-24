% +GDSGE/+RUNTIME  Hand-written runtime library used by generated code.
%   unpackOptions     - whitelisted GDSGE_OPTIONS unpacking (errors on unknown)
%   ensurePath        - essential_blas.dll system-PATH setup
%   solveProblems     - MEX driver + resolve cascade (owns the MEX caller contract)
%   applyWarmUp       - WarmUp SOL/LB/UB transfer
%   constructSplines  - interp_construct_mex interpolants + GDSGE_SPLINE_VEC
%   computeMetric     - sup-norm interp-update metric (NaN on error)
%   printIterProgress - structured progress line
%   reportUnconverged - failed-point report with state coordinates
