gdsge_codegen('safe_assets');
% Multi-root grid point: pin the RNG for a deterministic root and raise the
% retry cap so the nonlinear solve converges.
rng('default');
options.MaxMinorIter = 10000;
options.DiagnoseMinorIter = inf;
IterRslt = iter_safe_assets(options);
