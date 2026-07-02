% +GDSGE/+RUNTIME  Hand-written runtime library used by generated code.
%   unpackOptions        - whitelisted GDSGE_OPTIONS unpacking (errors on unknown)
%   ensurePath           - src/kernels path + DLL-loadability setup (self-register hook)
%   solveProblems        - MEX driver + resolve cascade (owns the MEX caller contract)
%   solveProblemsAsg     - ASG-variant MEX driver (randomized-restart resolve)
%   applyWarmUp          - WarmUp SOL/LB/UB transfer
%   constructSplines     - interp_construct_mex interpolants + GDSGE_SPLINE_VEC
%   computeMetric        - sup-norm interp-update metric (NaN on error)
%   printIterProgress    - checkpoint progress line + compact retry summary
%   printResolveProgress - throttled resolve round heartbeat + cap-exhaustion line
%   diagnoseConvergence  - tiered convergence diagnostics (summary/full)
%   reportUnconverged    - failed-point report with state coordinates
%   retryStatsInit       - zeroed retry-stats accumulator
%   retryStatsAdd        - add one solve call's retry rounds to the current bucket
%   retryStatsEndIter    - close the per-major-iteration retry bucket
