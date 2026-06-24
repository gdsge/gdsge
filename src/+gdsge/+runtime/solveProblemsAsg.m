function [sol, f, aux, eqVal, optInfo, diag] = solveProblemsAsg(mexFn, sol, lb, ub, data, f, aux, eqVal, cfg)
% SOLVEPROBLEMSASG  Solve one batch of ASG-proposed grid problems. Old
%   iter_solve_problem_asg_template.m parity:
%     1) solve everything;
%     2) warm unconverged from the CURRENT-iteration sol interp
%        (cfg.solInterpNew), tighten if cfg.useAdaptiveBound, re-solve;
%     3) restore remaining unconverged from the PREVIOUS-iteration sol interp
%        (cfg.solInterpOld), tighten — NO solve here (the randomize loop is
%        the next step);
%     4) randomized restarts until tolSol or maxMinorIter, with the
%        UseAdaptiveBoundInSol hitting-bounds adjustment (cfg.adaptInSol).
%   The init problem passes solInterpNew/solInterpOld = [], reducing the
%   cascade to 1+4 (old iter_solve_init_problem_asg_template.m).
%
%   CONTRACT: the MEX reads named variables from THIS function's workspace
%   via mexGetVariable('caller',...). Keep every mexFn(...) call in this
%   body. ASG models read GDSGE_ASG_HANDLE (the var-interp asg object of the
%   previous iteration); GDSGE_USE_BROYDEN_NOW is 0 on the ASG path, always
%   (old params_template.m:27 sets it once and no ASG template flips it).
%
%   cfg: tolSol tolFun solMaxIter numThreads debugEvalOnly useBroyden
%        finiteDiffDelta useBroydenNow taskName asgHandle solInterpNew
%        solInterpOld evalArrayIdx evalGrids useAdaptiveBound adaptTight
%        adaptInSol maxMinorIter verboseRetry
%   diag: .solved (logical row) .needResolved .minorIters
%   diag differs from solveProblems: adds .solved/.minorIters, omits .worstF.

% ---- MEX caller-workspace contract -----------------------------------------
TolFun = cfg.tolFun;                          %#ok<NASGU>
TolSol = cfg.tolSol;                          %#ok<NASGU>
SolMaxIter = cfg.solMaxIter;                  %#ok<NASGU>
NumThreads = cfg.numThreads;                  %#ok<NASGU>
GDSGE_DEBUG_EVAL_ONLY = cfg.debugEvalOnly;    %#ok<NASGU>
UseBroyden = cfg.useBroyden;                  %#ok<NASGU>
FiniteDiffDelta = cfg.finiteDiffDelta;        %#ok<NASGU>
GDSGE_USE_BROYDEN_NOW = cfg.useBroydenNow;    %#ok<NASGU>
MEX_TASK_NAME = cfg.taskName;                 %#ok<NASGU>
MEX_TASK_INIT = 0;                            %#ok<NASGU>
MEX_TASK_INF_HORIZON = 1;                     %#ok<NASGU>
GDSGE_SPLINE_VEC = [];                        %#ok<NASGU>  unused on ASG path; MEX contract compliance
GDSGE_ASG_HANDLE = cfg.asgHandle;             %#ok<NASGU>
% ----------------------------------------------------------------------------

% 1) Initial solve — attempt all problems from the supplied initial guess
f(:) = 1e20;
skip = zeros(1, size(sol, 2));
[sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);

% 2) Warm unconverged from the current-iteration sol interp, re-solve.
%    (Old: iter_solve_problem_asg_template.m warm-start block.)
if ~isempty(cfg.solInterpNew) && cfg.solInterpNew.get_current_level() >= 0
    needResolved = (f > cfg.tolSol) | isnan(f);
    if any(needResolved)
        sol0 = cfg.solInterpNew.eval_vec(cfg.evalArrayIdx, cfg.evalGrids);
        sol(:, needResolved) = sol0(:, needResolved);
        skip(:) = 0;
        skip(~needResolved) = 1;
        if cfg.useAdaptiveBound == 1
            [lb, ub] = cfg.adaptTight(sol, lb, ub);
        end
        [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);
    end
end

% 3) Restore remaining unconverged from the previous-iteration sol interp.
%    No re-solve here — the random-restart loop handles them.
%    (Old: the "restore from old" block before the minor-iter loop.)
if ~isempty(cfg.solInterpOld) && cfg.solInterpOld.get_current_level() >= 0
    needResolved = (f > cfg.tolSol) | isnan(f);
    if any(needResolved)  % skip eval_vec when everything already converged
        sol0 = cfg.solInterpOld.eval_vec(cfg.evalArrayIdx, cfg.evalGrids);
        sol(:, needResolved) = sol0(:, needResolved);
        if cfg.useAdaptiveBound == 1
            [lb, ub] = cfg.adaptTight(sol, lb, ub);
        end
    end
end

% 4) Randomized restarts — bounded by maxMinorIter.
%    (Old: the minor-iter loop in iter_solve_problem_asg_template.m,
%     including the UseAdaptiveBoundInSol hitting-bounds adjustment.)
minorIter = 0;
% Bare-cfg callers (unit tests) omit resolvePrintFreq/majorIter; default them.
if isfield(cfg, 'resolvePrintFreq'); resolveFreq = cfg.resolvePrintFreq; else; resolveFreq = 100; end
if isfield(cfg, 'majorIter'); majorIter = cfg.majorIter; else; majorIter = 0; end
while ((max(isnan(f)) || max(f(:)) > cfg.tolSol) && minorIter < cfg.maxMinorIter)
    x0Rand = rand(size(sol)) .* (ub - lb) + lb;
    needResolved = (f > cfg.tolSol) | isnan(f);
    sol(:, needResolved) = x0Rand(:, needResolved);
    skip(:) = 0;
    skip(~needResolved) = 1;
    [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);
    if ~isempty(cfg.adaptInSol)
        lbOld = lb; ubOld = ub;
        [lb, ub] = cfg.adaptInSol(sol, lb, ub);
        hitLower = abs(sol - lbOld) < 1e-8;
        hitUpper = abs(sol - ubOld) < 1e-8;
        lb(~hitLower) = lbOld(~hitLower);
        ub(~hitUpper) = ubOld(~hitUpper);
    end
    % Deliberate deviation: the old ASG template increments only under
    % UseAdaptiveBoundInSol (template line 66), leaving the loop unbounded when
    % that flag is 0. MaxMinorIter defaults to inf, so behavior is identical at
    % defaults; with a finite cap the new code honors it (improvement).
    minorIter = minorIter + 1;
    gdsge.runtime.printResolveProgress('round', 'asg ', majorIter, minorIter, 1, ...
        nnz((f > cfg.tolSol) | isnan(f)), max(f(:)), resolveFreq, cfg.verboseRetry, false);
end

gdsge.runtime.printResolveProgress('summary', 'asg ', majorIter, minorIter, 0, ...
    nnz((f > cfg.tolSol) | isnan(f)), max(f(:)), resolveFreq, cfg.verboseRetry, ...
    (max(isnan(f)) || max(f(:)) > cfg.tolSol));

diag = struct();
diag.needResolved = (f > cfg.tolSol) | isnan(f);
diag.solved = ~diag.needResolved;
diag.minorIters = minorIter;
end
