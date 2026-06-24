function [sol, f, aux, eqVal, optInfo, diag] = solveProblems(mexFn, sol, lb, ub, data, f, aux, eqVal, cfg)
% SOLVEPROBLEMS  Drive the model MEX over all grid-point problems until every
%   residual is below tolSol (or maxMinorIter is exhausted). Behavior-parity
%   port of the old generated resolve cascade (see scratch/old_generated/).
%
%   CONTRACT: the MEX reads named variables from ITS CALLER's workspace via
%   mexGetVariable('caller',...). Every mexFn(...) call below must therefore
%   stay in THIS function's body (never a subfunction), and the contract
%   variables below must exist here as locals.
%
%   cfg: tolSol tolFun solMaxIter numThreads debugEvalOnly useBroyden
%        finiteDiffDelta useBroydenNow taskName splineVec ppNames ppCell
%        maxMinorIter probSize useNearestNeighbor verboseRetry adaptInSol asgHandle
%        useMexResolve
%   ppNames/ppCell: per-interp pp structs the MEX pulls by model-specific
%        name (GET_STRUCT(GDSGE_PP_<interp>)); created as locals below.
%   adaptInSol: [] or @(sol,lb,ub) -> [lb,ub] (UseAdaptiveBoundInSol hook)
%   diag: .needResolved (logical row) .worstF .minorIters

% ---- MEX caller-workspace contract (read via mexGetVariable) --------------
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
GDSGE_SPLINE_VEC = cfg.splineVec;             %#ok<NASGU>
GDSGE_PROBLEM_STRIDES = [];                   %#ok<NASGU>  % MEX in-resolve gate (set per sweep call)
useMexResolve = isfield(cfg,'useMexResolve') && cfg.useMexResolve;
if useMexResolve
    GDSGE_RESOLVE_STRIDES = cumprod([1, cfg.probSize(1:end-1)]);  % prod(probSize(1:d-1)), d=1..D
end
if isfield(cfg, 'asgHandle')
    GDSGE_ASG_HANDLE = cfg.asgHandle;         %#ok<NASGU>
end
for iPP = 1:numel(cfg.ppNames)
    % model-named interp structs, e.g. GDSGE_PP_ps_future (identifiers come
    % from the parser-validated IR; eval only creates workspace locals)
    assert(isvarname(cfg.ppNames{iPP}));
    eval([cfg.ppNames{iPP} ' = cfg.ppCell{iPP};']);
end
% ---------------------------------------------------------------------------

f(:) = 1e20;
skip = zeros(1, size(sol, 2));
[sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);
GDSGE_USE_BROYDEN_NOW = 0;   %#ok<NASGU> % parity: Broyden only on the first call

minorIter = 0;
numNeedResolvedAfter = inf;

% Bare-cfg callers (unit tests) omit resolvePrintFreq/majorIter; default them.
if isfield(cfg, 'resolvePrintFreq'); resolveFreq = cfg.resolvePrintFreq; else; resolveFreq = 100; end
if isfield(cfg, 'majorIter'); majorIter = cfg.majorIter; else; majorIter = 0; end

useMexRandomize = isfield(cfg,'useMexRandomize') && cfg.useMexRandomize ...
    && isempty(cfg.adaptInSol) && useMexResolve;
if useMexRandomize
    % In-MEX randomized restart: each MEX call runs up to `batch` full
    % minor-iterations (neighbor sweep + random restart) internally. MATLAB
    % keeps diagnostics, the cap, and the exit decision at batch boundaries.
    GDSGE_RANDOM_SEED = cfg.randomSeed;            %#ok<NASGU> caller-workspace contract
    GDSGE_RANDOMIZE_SALT = cfg.randomizeSalt;      %#ok<NASGU>
    GDSGE_RANDOMIZE_TRIAL_OFFSET = 0;             %#ok<NASGU>
    GDSGE_RANDOMIZE_BATCH = 0;                    %#ok<NASGU>
    trialOffset = 0;
    while ((max(isnan(f)) || max(f(:)) > cfg.tolSol) && minorIter < cfg.maxMinorIter)
        batch = min(cfg.randomizeBatch, cfg.maxMinorIter - minorIter);
        GDSGE_RANDOMIZE_BATCH = batch;                       %#ok<NASGU>
        GDSGE_RANDOMIZE_TRIAL_OFFSET = trialOffset;          %#ok<NASGU>
        GDSGE_PROBLEM_STRIDES = GDSGE_RESOLVE_STRIDES;       %#ok<NASGU> enable in-MEX neighbor sweep
        % The MEX breaks its internal trial loop early once everything converges,
        % so it reports how many restart rounds it actually ran via this caller
        % var (mexPutVariable). Pre-seed = batch so a MATLAB fake MEX that does
        % not report leaves the legacy full-batch count untouched.
        GDSGE_RANDOMIZE_TRIALS_USED = batch;                 %#ok<NASGU>
        skip(:) = 1;                                         % Pass-0 no-op; restarts gate on f
        [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);
        GDSGE_RANDOMIZE_BATCH = 0;                  %#ok<NASGU> off for any later call
        GDSGE_PROBLEM_STRIDES = [];                 %#ok<NASGU>
        actualTrials = GDSGE_RANDOMIZE_TRIALS_USED; % real restart rounds this call
        minorIter = minorIter + actualTrials;
        trialOffset = trialOffset + actualTrials;
        gdsge.runtime.printResolveProgress('round', '', majorIter, minorIter, actualTrials, ...
            nnz((f > cfg.tolSol) | isnan(f)), max(f(:)), resolveFreq, cfg.verboseRetry, false);

        % Periodic diagnostics heartbeat (now at batch granularity).
        if isfield(cfg,'diagnoseAt') && cfg.diagnoseAt > 0 ...
                && mod(minorIter, cfg.diagnoseAt) == 0 && minorIter < cfg.maxMinorIter ...
                && (max(isnan(f)) || max(f(:)) > cfg.tolSol)
            fprintf('[gdsge] convergence diagnostics (retry %d, continuing):\n%s\n', ...
                minorIter, gdsge.runtime.diagnoseConvergence(sol, lb, ub, f, eqVal, cfg, 'summary'));
        end
    end
else
while ((max(isnan(f)) || max(f(:)) > cfg.tolSol) && minorIter < cfg.maxMinorIter)
    if useMexResolve
        % In-MEX neighbor sweep. Replicate the OLD inner-while ENTRY gate exactly
        % (numNeedResolvedAfter persists across outer iterations, init = inf), then
        % let the MEX run the entire progress loop internally in ONE call.
        needResolved = (f > cfg.tolSol) | isnan(f);
        numNeedResolved = sum(needResolved);
        if numNeedResolvedAfter ~= numNeedResolved
            skip(:) = 1;                                   % Pass 0 no-ops; strides drive the sweep
            GDSGE_PROBLEM_STRIDES = GDSGE_RESOLVE_STRIDES; %#ok<NASGU> % enable in-MEX resolve
            [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);
            GDSGE_PROBLEM_STRIDES = [];                    %#ok<NASGU> % off for the restart call below
            numNeedResolvedAfter = sum((f > cfg.tolSol) | isnan(f));
        end
    elseif cfg.useNearestNeighbor
        needResolved = (f > cfg.tolSol) | isnan(f);
        numNeedResolved = sum(needResolved);
        while numNeedResolvedAfter ~= numNeedResolved
            needResolved = (f > cfg.tolSol) | isnan(f);
            numNeedResolved = sum(needResolved);
            for iDim = 1:length(cfg.probSize)
                stride = prod(cfg.probSize(1:iDim-1));

                % copy from the lower neighbor (parity: stale needResolved
                % during the sweep, exactly like the old generated loop)
                needResolved = (f > cfg.tolSol) | isnan(f);
                skip(:) = 1;
                for i = 1:numel(f)
                    if needResolved(i) && i-stride >= 1 && ~needResolved(i-stride)
                        sol(:, i) = sol(:, i-stride);
                        skip(i) = 0;
                    end
                end
                [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);

                % copy from the upper neighbor
                needResolved = (f > cfg.tolSol) | isnan(f);
                skip(:) = 1;
                for i = 1:numel(f)
                    if needResolved(i) && i+stride <= numel(f) && ~needResolved(i+stride)
                        sol(:, i) = sol(:, i+stride);
                        skip(i) = 0;
                    end
                end
                [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);
            end
            numNeedResolvedAfter = sum((f > cfg.tolSol) | isnan(f));
        end
    end

    % randomize restarts for whatever is left (parity: old lines 376-385)
    x0Rand = rand(size(sol)) .* (ub - lb) + lb;
    needResolved = (f > cfg.tolSol) | isnan(f);
    sol(:, needResolved) = x0Rand(:, needResolved);
    skip(:) = 0;
    skip(~needResolved) = 1;
    [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);

    if ~isempty(cfg.adaptInSol)
        % UseAdaptiveBoundInSol hook (parity: old lines 387-408)
        lbOld = lb; ubOld = ub;
        [lb, ub] = cfg.adaptInSol(sol, lb, ub);
        hitLower = abs(sol - lbOld) < 1e-8;
        hitUpper = abs(sol - ubOld) < 1e-8;
        lb(~hitLower) = lbOld(~hitLower);
        ub(~hitUpper) = ubOld(~hitUpper);
        minorIter = minorIter + 1;   % parity: old counts an extra minor iter here
    end
    minorIter = minorIter + 1;

    roundStep = 1 + double(~isempty(cfg.adaptInSol));
    gdsge.runtime.printResolveProgress('round', '', majorIter, minorIter, roundStep, ...
        nnz((f > cfg.tolSol) | isnan(f)), max(f(:)), resolveFreq, cfg.verboseRetry, false);

    % Periodic diagnostics: a summary heartbeat every cfg.diagnoseAt retries
    % while still stuck. The final block prints after the loop, so skip the cap
    % iteration here to avoid a duplicate at minorIter == maxMinorIter.
    if isfield(cfg,'diagnoseAt') && cfg.diagnoseAt > 0 ...
            && mod(minorIter, cfg.diagnoseAt) == 0 && minorIter < cfg.maxMinorIter ...
            && (max(isnan(f)) || max(f(:)) > cfg.tolSol)
        fprintf('[gdsge] convergence diagnostics (retry %d, continuing):\n%s\n', ...
            minorIter, gdsge.runtime.diagnoseConvergence(sol, lb, ub, f, eqVal, cfg, 'summary'));
    end
end
end   % closes the `else` wrapping the legacy resolve loop (useMexRandomize branch)

gdsge.runtime.printResolveProgress('summary', '', majorIter, minorIter, 0, ...
    nnz((f > cfg.tolSol) | isnan(f)), max(f(:)), resolveFreq, cfg.verboseRetry, ...
    (max(isnan(f)) || max(f(:)) > cfg.tolSol));

% Exhaustion: print the full block; error on the iter path, warn-and-continue
% on simulate (gated on the diagnostic fields, so bare-cfg callers are unchanged).
if isfield(cfg,'diagnoseAt') && (max(isnan(f)) || max(f(:)) > cfg.tolSol)
    fprintf('[gdsge] FINAL convergence diagnostics (retry %d):\n%s\n', ...
        minorIter, gdsge.runtime.diagnoseConvergence(sol, lb, ub, f, eqVal, cfg, 'full'));
    if isfield(cfg,'errorOnNonconvergence') && cfg.errorOnNonconvergence
        error('gdsge:runtime:solveNotConverged', ...
            'Solve did not converge: %d/%d points after %d retries. See diagnostics above.', ...
            nnz((f > cfg.tolSol) | isnan(f)), numel(f), minorIter);
    end
end

diag = struct();
diag.needResolved = (f > cfg.tolSol) | isnan(f);
diag.worstF = max(f(:));
diag.minorIters = minorIter;
end
