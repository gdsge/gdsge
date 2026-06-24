# Bounded Retries + Tiered Convergence Diagnostics — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the silent, unbounded random-restart loop in the cartesian solve path with a bounded retry budget plus named, tiered diagnostics — so a non-converging model stops with an actionable report (`iter` errors; `simulate` warns and continues) instead of looping forever.

**Architecture:** All behavior is centralized in `gdsge.runtime.solveProblems` (shared by `iter` and `simulate`, both Jacobian backends). A new pure formatter `gdsge.runtime.diagnoseConvergence` builds the report; codegen threads named labels + the schedule into `GDSGE_CFG`. The new behavior is gated on `isfield(cfg,'diagnoseAt')`, so existing bare-cfg unit tests keep their exact current behavior.

**Tech Stack:** MATLAB R2025b, `matlab.unittest`. Source under `src/+gdsge/`, tests auto-discovered under `tests/` (`TestSuite.fromFolder(..., 'IncludingSubfolders', true)`).

**Spec:** `docs/superpowers/specs/2026-06-14-convergence-diagnostics-design.md`

---

## Conventions for every command

All commands run from the repo root `D:\refactor_gdsge`. MATLAB is
`C:\Program Files\MATLAB\R2025b\bin\matlab.exe` (assume `matlab` is on PATH;
otherwise use the full path).

**Run one test class** (template — substitute the path):

```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/runtime/tDiagnoseConvergence.m'); assert(~any([r.Failed]))"
```

A failing assertion / test makes `assert` throw, so the process exits non-zero.

**Run the full suite** (exit 0 = all pass; report at `tests/results/junit.xml`):

```
matlab -batch "cd('tests'); run_tests"
```

**Golden rule:** never add a toolbox source to a persistent MATLAB path; each
command above `addpath`s exactly `src/` for one process and exits.

---

## Task 0: Branch

- [ ] **Step 1: Create the feature branch**

```
git checkout -b convergence-diagnostics
```

---

## Task 1: `solComponentNames` codegen helper

Maps each `GDSGE_SOL`/`LB`/`UB` row to a named policy component (`kp`, `w1n(3)`),
length = total solution rows. Shared by `emitIter` and `emitSimulate`.

**Files:**
- Create: `src/+gdsge/+codegen/+mat/solComponentNames.m`
- Test: `tests/codegen/tSolComponentNames.m`

- [ ] **Step 1: Write the failing test**

Create `tests/codegen/tSolComponentNames.m`:

```matlab
classdef tSolComponentNames < matlab.unittest.TestCase
    methods (Test)
        function scalarAndArrayComponents(tc)
            ir.variables.policy = { ...
                struct('name','c',  'length',1, 'slot',[1 1]), ...
                struct('name','kp', 'length',3, 'slot',[2 4])};
            names = gdsge.codegen.mat.solComponentNames(ir);
            tc.verifyEqual(names, {'c','kp(1)','kp(2)','kp(3)'});
        end
        function orderedBySlotNotDeclaration(tc)
            % declaration order shuffled; output must follow slot order
            ir.variables.policy = { ...
                struct('name','b','length',2,'slot',[3 4]), ...
                struct('name','a','length',2,'slot',[1 2])};
            names = gdsge.codegen.mat.solComponentNames(ir);
            tc.verifyEqual(names, {'a(1)','a(2)','b(1)','b(2)'});
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tSolComponentNames.m'); assert(~any([r.Failed]))"
```

Expected: FAIL — `Unrecognized function 'gdsge.codegen.mat.solComponentNames'`.

- [ ] **Step 3: Write the implementation**

Create `src/+gdsge/+codegen/+mat/solComponentNames.m`:

```matlab
function names = solComponentNames(ir)
% SOLCOMPONENTNAMES  Per-row labels for the GDSGE_SOL/LB/UB solution stack.
%   Returns a 1×maxDim cell aligned with the policy slot layout: a scalar
%   policy var contributes its name; an array var of length L contributes
%   name(1)..name(L). Ordered by slot so it matches GDSGE_SOL's row order.
pol = ir.variables.policy;
maxDim = 0;
for i = 1:numel(pol); maxDim = max(maxDim, pol{i}.slot(2)); end
names = repmat({''}, 1, maxDim);
for i = 1:numel(pol)
    p = pol{i};
    lo = p.slot(1); hi = p.slot(2);
    if p.length == 1
        names{lo} = p.name;
    else
        for k = lo:hi
            names{k} = sprintf('%s(%d)', p.name, k - lo + 1);
        end
    end
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tSolComponentNames.m'); assert(~any([r.Failed]))"
```

Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```
git add src/+gdsge/+codegen/+mat/solComponentNames.m tests/codegen/tSolComponentNames.m
git commit -m "feat(codegen): solComponentNames maps SOL rows to named policy components" -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: `diagnoseConvergence` runtime formatter

Pure formatter producing the diagnostic report (NaN, bounds-pinned among
unconverged, near-bound among converged, dominant equation, worst-K points).

**Files:**
- Create: `src/+gdsge/+runtime/diagnoseConvergence.m`
- Test: `tests/runtime/tDiagnoseConvergence.m`

- [ ] **Step 1: Write the failing test**

Create `tests/runtime/tDiagnoseConvergence.m`:

```matlab
classdef tDiagnoseConvergence < matlab.unittest.TestCase
    methods (Static)
        function [sol,lb,ub,f,eqVal,cfg] = fixture()
            % 2 components × 4 problems. Pt1 converged (c near LB); pts 2-4
            % unconverged with kp pinned at UB; pt4 is a NaN residual.
            sol   = [0.5   1    1    1;        % kp
                     0.005 0.5  0.5  0.5];     % c
            lb    = zeros(2,4);
            ub    = ones(2,4);
            f     = [0 10 10 NaN];
            eqVal = [0  50  60  NaN;           % eq#1 inflated
                     0  0.5 0.6 NaN];          % eq#2 small
            cfg = struct('tolSol',1e-8,'probSize',[1 4],'solNames',{{'kp','c'}});
        end
    end
    methods (Test)
        function fullReportNamesAllSignals(tc)
            [sol,lb,ub,f,eqVal,cfg] = tDiagnoseConvergence.fixture();
            msg = gdsge.runtime.diagnoseConvergence(sol,lb,ub,f,eqVal,cfg,'full');
            tc.verifyTrue(contains(msg, '3 of 4 problems unconverged'));
            tc.verifyTrue(contains(msg, 'NaN residuals: 1'));
            tc.verifyTrue(contains(msg, 'kp: 3 at UB'));         % pinned, named
            tc.verifyTrue(contains(msg, 'c: 1 near LB'));        % near-bound, converged
            tc.verifyTrue(contains(msg, 'eq#1'));                % dominant equation
            tc.verifyTrue(contains(msg, 'dominant'));
        end
        function summaryIsHeadlineSubset(tc)
            [sol,lb,ub,f,eqVal,cfg] = tDiagnoseConvergence.fixture();
            msg = gdsge.runtime.diagnoseConvergence(sol,lb,ub,f,eqVal,cfg,'summary');
            tc.verifyTrue(contains(msg, 'kp: 3 at UB'));         % headline kept
            tc.verifyTrue(contains(msg, 'eq#1'));
            tc.verifyFalse(contains(msg, 'near LB'));            % near-bound is full-only
        end
        function fallsBackToIndexNamesWhenAbsent(tc)
            [sol,lb,ub,f,eqVal,cfg] = tDiagnoseConvergence.fixture();
            cfg = rmfield(cfg, 'solNames');
            msg = gdsge.runtime.diagnoseConvergence(sol,lb,ub,f,eqVal,cfg,'summary');
            tc.verifyTrue(contains(msg, 'x1: 3 at UB'));         % row 1 -> x1
        end
        function pointTableAppendedWhenGridsPresent(tc)
            sol = [1 1 1]; lb = zeros(1,3); ub = ones(1,3);
            f = [0 7 3]; eqVal = [0 7 3];
            cfg = struct('tolSol',1e-8,'probSize',[1 3],'solNames',{{'kp'}}, ...
                'stateNames',{{'w1'}},'stateGrids',{{[0.1 0.5 0.9]}});
            msg = gdsge.runtime.diagnoseConvergence(sol,lb,ub,f,eqVal,cfg,'full');
            tc.verifyTrue(contains(msg, 'w1='));                 % reportUnconverged table
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/runtime/tDiagnoseConvergence.m'); assert(~any([r.Failed]))"
```

Expected: FAIL — `Unrecognized function 'gdsge.runtime.diagnoseConvergence'`.

- [ ] **Step 3: Write the implementation**

Create `src/+gdsge/+runtime/diagnoseConvergence.m`:

```matlab
function msg = diagnoseConvergence(sol, lb, ub, f, eqVal, cfg, level)
% DIAGNOSECONVERGENCE  Human-readable diagnostics for a stuck cartesian solve.
%   Pure formatter — returns the report text; the caller decides fprintf/error.
%   level : 'summary' (tier-1 headlines) | 'full' (all signals + point table)
%   cfg   : .tolSol .probSize .solNames (cell, numel == size(sol,1));
%           optional .stateNames/.stateGrids (point table appended iff present).
%   Signals: NaN count; bounds pinned among unconverged (per named component);
%   near-bound among converged (full only); dominant equation residual vs the
%   other equations; worst-K grid points (full only — iter path).
tolSol = cfg.tolSol;
[nSol, nProb] = size(sol);

if isfield(cfg,'solNames') && numel(cfg.solNames) == nSol
    names = cfg.solNames;
else
    names = arrayfun(@(r) sprintf('x%d', r), 1:nSol, 'UniformOutput', false);
end

unc  = (f > tolSol) | isnan(f);          % 1 × nProb
conv = ~unc;
nUnc = nnz(unc);
span = max(1, ub - lb);                  % nSol × nProb; guards tiny ranges

L = {};
L{end+1} = sprintf('%d of %d problems unconverged (TolSol=%g).', nUnc, nProb, tolSol);

% --- NaN ---
nNan = nnz(isnan(f));
if nNan > 0
    L{end+1} = sprintf('  NaN residuals: %d (suspect log/sqrt of a negative or divide-by-zero in the model body).', nNan);
end

% --- bounds pinned among unconverged (per component) ---
atLB = abs(sol - lb) <= 1e-6 .* span;
atUB = abs(ub - sol) <= 1e-6 .* span;
pinLB = sum(atLB(:, unc), 2);
pinUB = sum(atUB(:, unc), 2);
[~, ord] = sort(pinLB + pinUB, 'descend');
shown = false;
for k = 1:min(5, nSol)
    r = ord(k);
    if pinLB(r) == 0 && pinUB(r) == 0; break; end
    if ~shown
        L{end+1} = '  Bounds (unconverged points pinned within 1e-6 of span):';
        shown = true;
    end
    p = {};
    if pinUB(r) > 0; p{end+1} = sprintf('%d at UB', pinUB(r)); end %#ok<AGROW>
    if pinLB(r) > 0; p{end+1} = sprintf('%d at LB', pinLB(r)); end %#ok<AGROW>
    L{end+1} = sprintf('    %s: %s (of %d unconverged) -- bound may be too tight', ...
        names{r}, strjoin(p, ', '), nUnc);
end

% --- dominant equation residual (vs the other equations) ---
if ~isempty(eqVal) && nUnc > 0
    rowMax = max(abs(eqVal(:, unc)), [], 2);   % nEq × 1 (max ignores NaN)
    rowMax(isnan(rowMax)) = 0;
    [sorted, eord] = sort(rowMax, 'descend');
    if numel(sorted) > 1; medOthers = median(sorted(2:end)); else; medOthers = sorted(1); end
    if medOthers == 0; medOthers = eps; end
    L{end+1} = sprintf('  Largest equation residuals over unconverged points (vs median of the others = %.3g):', medOthers);
    for k = 1:min(3, numel(sorted))
        e = eord(k);
        tag = '';
        if k == 1 && sorted(1) > 5 * medOthers
            tag = '  <-- dominant; likely a bug or needs normalization';
        end
        L{end+1} = sprintf('    eq#%d: %.3g (%.1fx others)%s', e, rowMax(e), rowMax(e) / medOthers, tag);
    end
end

if strcmp(level, 'summary')
    msg = strjoin(L, newline);
    return;
end

% --- full-only: near-bound among converged (per component) ---
nearBand = 0.01;
nearLB = sum(((sol - lb) <= nearBand .* span) & conv, 2);
nearUB = sum(((ub - sol) <= nearBand .* span) & conv, 2);
[~, nord] = sort(nearLB + nearUB, 'descend');
shown = false;
for k = 1:min(5, nSol)
    r = nord(k);
    if nearLB(r) == 0 && nearUB(r) == 0; break; end
    if ~shown
        L{end+1} = '  Bounds (converged points within 1% of span):';
        shown = true;
    end
    p = {};
    if nearUB(r) > 0; p{end+1} = sprintf('%d near UB', nearUB(r)); end %#ok<AGROW>
    if nearLB(r) > 0; p{end+1} = sprintf('%d near LB', nearLB(r)); end %#ok<AGROW>
    L{end+1} = sprintf('    %s: %s', names{r}, strjoin(p, ', '));
end

% --- full-only: worst-K grid points (iter path supplies stateGrids) ---
if isfield(cfg,'stateGrids') && ~isempty(cfg.stateGrids)
    L{end+1} = gdsge.runtime.reportUnconverged(unc, f, cfg.probSize, ...
        cfg.stateNames, cfg.stateGrids, 5);
end

msg = strjoin(L, newline);
end
```

- [ ] **Step 4: Run the test to verify it passes**

```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/runtime/tDiagnoseConvergence.m'); assert(~any([r.Failed]))"
```

Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```
git add src/+gdsge/+runtime/diagnoseConvergence.m tests/runtime/tDiagnoseConvergence.m
git commit -m "feat(runtime): diagnoseConvergence formatter (bounds/NaN/dominant-eq/points)" -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Wire the bounded retry + tiered diagnostics into `solveProblems`

Add a once-only tier-1 print at `cfg.diagnoseAt`, and a post-loop exhaustion
handler: print the full block, then `error(...)` iff
`cfg.errorOnNonconvergence`. All gated on `isfield(cfg,'diagnoseAt')` so the
existing bare-cfg unit tests are unchanged.

**Files:**
- Modify: `src/+gdsge/+runtime/solveProblems.m`
- Test: `tests/runtime/tSolveProblemsCap.m` (create)
- Regression guard: `tests/runtime/tSolveProblems.m` (must stay green, no edits)

- [ ] **Step 1: Write the failing test**

Create `tests/runtime/tSolveProblemsCap.m`:

```matlab
classdef tSolveProblemsCap < matlab.unittest.TestCase
    % Bounded-retry + tiered-diagnostics behavior of solveProblems when the
    % full diagnostic cfg is supplied (the generated iter/simulate path).
    methods (Test)
        function iterPathThrowsAfterCap(tc)
            cfg = capCfg();
            cfg.errorOnNonconvergence = true;        % iter path
            tc.verifyError(@() runNever(cfg), 'gdsge:runtime:solveNotConverged');
        end
        function simulatePathPrintsAndReturns(tc)
            cfg = capCfg();
            cfg.errorOnNonconvergence = false;       % simulate path
            out = evalc('[fOut, dg] = runNever(cfg);');
            tc.verifyTrue(all(fOut > cfg.tolSol));
            tc.verifyTrue(all(dg.needResolved));
            tc.verifyEqual(dg.minorIters, 4);
            tc.verifyTrue(contains(out, 'convergence diagnostics (retry 2'));  % tier-1
            tc.verifyTrue(contains(out, 'FINAL convergence diagnostics'));     % tier-2
            tc.verifyTrue(contains(out, 'eq#1'));
        end
        function disabledWithoutDiagnoseAt(tc)
            % No diagnoseAt -> exact legacy behavior: cap, return, no print/throw.
            cfg = capCfg();
            cfg = rmfield(cfg, 'diagnoseAt');
            out = evalc('[fOut, dg] = runNever(cfg);');
            tc.verifyEqual(out, '');
            tc.verifyEqual(dg.minorIters, 4);
            tc.verifyTrue(all(dg.needResolved));
        end
    end
end

function [fOut, dg] = runNever(cfg)
n = 2;
[~, fOut, ~, ~, ~, dg] = gdsge.runtime.solveProblems(@fakeNever, ...
    zeros(2,n), zeros(2,n), ones(2,n), zeros(1,n), 1e20*ones(1,n), ...
    zeros(1,n), zeros(2,n), cfg);
end

function [sol, f, aux, eqVal, optInfo] = fakeNever(sol, ~, ~, ~, skip, f, aux, eqVal)
f(~logical(skip)) = 1e20;            % never below tolSol
eqVal(1,:) = 100; eqVal(2,:) = 1;    % eq#1 dominant
optInfo = zeros(1, numel(f));
end

function cfg = capCfg()
cfg = struct('tolSol',1e-8,'tolFun',1e-8,'solMaxIter',200,'numThreads',1, ...
    'debugEvalOnly',0,'useBroyden',0,'finiteDiffDelta',1e-6,'useBroydenNow',0, ...
    'taskName',1,'splineVec',struct(),'ppNames',{{}},'ppCell',{{}}, ...
    'maxMinorIter',4,'probSize',[1 2],'useNearestNeighbor',false, ...
    'verboseRetry',false,'adaptInSol',[], ...
    'diagnoseAt',2,'solNames',{{'a','b'}});
end
```

- [ ] **Step 2: Run the test to verify it fails**

```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/runtime/tSolveProblemsCap.m'); assert(~any([r.Failed]))"
```

Expected: FAIL — `iterPathThrowsAfterCap` does not throw, and the printed-text
assertions fail (no diagnostics emitted yet).

- [ ] **Step 3: Modify `solveProblems`**

In `src/+gdsge/+runtime/solveProblems.m`, just before the `while` loop add the
tier-1 flag. Change:

```matlab
minorIter = 0;
numNeedResolvedAfter = inf;
while ((max(isnan(f)) || max(f(:)) > cfg.tolSol) && minorIter < cfg.maxMinorIter)
```

to:

```matlab
minorIter = 0;
numNeedResolvedAfter = inf;
diagnosedTier1 = false;
while ((max(isnan(f)) || max(f(:)) > cfg.tolSol) && minorIter < cfg.maxMinorIter)
```

Then, inside the loop, immediately after the existing `verboseRetry` block:

```matlab
    if cfg.verboseRetry
        fprintf('  resolve round %d: %d unconverged, worst residual %g\n', ...
            minorIter, nnz((f > cfg.tolSol) | isnan(f)), max(f(:)));
    end
```

insert the tier-1 block (still inside the `while`, before its closing `end`):

```matlab
    if isfield(cfg,'diagnoseAt') && ~diagnosedTier1 && minorIter >= cfg.diagnoseAt ...
            && (max(isnan(f)) || max(f(:)) > cfg.tolSol)
        fprintf('[gdsge] convergence diagnostics (retry %d, continuing):\n%s\n', ...
            minorIter, gdsge.runtime.diagnoseConvergence(sol, lb, ub, f, eqVal, cfg, 'summary'));
        diagnosedTier1 = true;
    end
```

Finally, after the `while ... end` and before the `diag = struct();` line, insert
the exhaustion handler:

```matlab
if isfield(cfg,'diagnoseAt') && (max(isnan(f)) || max(f(:)) > cfg.tolSol)
    fprintf('[gdsge] FINAL convergence diagnostics (retry %d):\n%s\n', ...
        minorIter, gdsge.runtime.diagnoseConvergence(sol, lb, ub, f, eqVal, cfg, 'full'));
    if isfield(cfg,'errorOnNonconvergence') && cfg.errorOnNonconvergence
        error('gdsge:runtime:solveNotConverged', ...
            'Solve did not converge: %d/%d points after %d retries. See diagnostics above.', ...
            nnz((f > cfg.tolSol) | isnan(f)), numel(f), minorIter);
    end
end
```

(The existing `diag` struct return stays exactly as-is.)

- [ ] **Step 4: Run the new test + the regression guard to verify they pass**

```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests({'tests/runtime/tSolveProblemsCap.m','tests/runtime/tSolveProblems.m'}); assert(~any([r.Failed]))"
```

Expected: PASS — new cap tests pass; the existing `tSolveProblems` (incl.
`maxMinorIterCapsAndReportsDiag`, which supplies no `diagnoseAt`) is unchanged.

- [ ] **Step 5: Commit**

```
git add src/+gdsge/+runtime/solveProblems.m tests/runtime/tSolveProblemsCap.m
git commit -m "feat(runtime): bound retries + tiered diagnostics in solveProblems" -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Options — `MaxMinorIter` default 20 + new `DiagnoseMinorIter`

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitSetup.m:22`
- Modify: `src/+gdsge/+codegen/+mat/optionsWhitelist.m:6-12`
- Test: `tests/codegen/tEmitSetup.m:31` (update) + add a case
- Test: `tests/codegen/tOptionsWhitelist.m` (add a case)

- [ ] **Step 1: Update the failing tests first**

In `tests/codegen/tEmitSetup.m`, change line 31 from:

```matlab
            tc.verifyEqual(ws.MaxMinorIter, inf);
```

to:

```matlab
            tc.verifyEqual(ws.MaxMinorIter, 20);
            tc.verifyEqual(ws.DiagnoseMinorIter, 10);
```

In `tests/codegen/tOptionsWhitelist.m`, add this test inside `methods (Test)`:

```matlab
        function includesDiagnoseMinorIter(tc)
            ir = gdsgefix.minimalIR();
            lit = gdsge.codegen.mat.optionsWhitelist(ir, {});
            tc.verifyTrue(contains(lit, '''DiagnoseMinorIter'''));
        end
```

- [ ] **Step 2: Run the tests to verify they fail**

```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests({'tests/codegen/tEmitSetup.m','tests/codegen/tOptionsWhitelist.m'}); assert(~any([r.Failed]))"
```

Expected: FAIL — `MaxMinorIter` is still `inf`, `DiagnoseMinorIter` undefined,
whitelist lacks `DiagnoseMinorIter`.

- [ ] **Step 3: Make the source changes**

In `src/+gdsge/+codegen/+mat/emitSetup.m`, change line 22 from:

```matlab
w.add('MaxMinorIter = inf;');
```

to:

```matlab
w.add('MaxMinorIter = 20;');         % bounded retries; was inf (could loop forever)
w.add('DiagnoseMinorIter = 10;');    % retry count at which tiered diagnostics first print
```

In `src/+gdsge/+codegen/+mat/optionsWhitelist.m`, add `'DiagnoseMinorIter'` to
the `base` cell — change the last line of the literal from:

```matlab
    'CONSTRUCT_OUTPUT','NumThreads','WarmUp','UseMexResolve'};
```

to:

```matlab
    'CONSTRUCT_OUTPUT','NumThreads','WarmUp','UseMexResolve','DiagnoseMinorIter'};
```

- [ ] **Step 4: Run the tests to verify they pass**

```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests({'tests/codegen/tEmitSetup.m','tests/codegen/tOptionsWhitelist.m'}); assert(~any([r.Failed]))"
```

Expected: PASS.

- [ ] **Step 5: Commit**

```
git add src/+gdsge/+codegen/+mat/emitSetup.m src/+gdsge/+codegen/+mat/optionsWhitelist.m tests/codegen/tEmitSetup.m tests/codegen/tOptionsWhitelist.m
git commit -m "feat(codegen): MaxMinorIter default 20 + DiagnoseMinorIter option" -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Thread the diagnostic cfg into `emitIter` (and drop the dead warning block)

`emitIter` builds `GDSGE_CFG` once before the loop — add the four diagnostic
fields there, and remove the now-dead post-return `if any(NeedResolved)
warning(...)` block (the variable `NeedResolved` and `IterRslt.NeedResolved`
**stay**).

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitIter.m`
- Snapshot: `tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt` (regen)

- [ ] **Step 1: Add the solNames literal near the other `*Cell` literals**

In `src/+gdsge/+codegen/+mat/emitIter.m`, after the line:

```matlab
stateGridCell = ['{' stateList '}'];                                % {w1}
```

add:

```matlab
solComp     = gdsge.codegen.mat.solComponentNames(ir);
solNamesLit = ['{''' strjoin(solComp, ''',''') '''}'];              % {'c',...,'w1n(8)'}
```

- [ ] **Step 2: Add the diagnostic cfg fields to the built-once cfg block**

Still in `emitIter.m`, in the `GDSGE_CFG` built-once section, after:

```matlab
w.add('GDSGE_CFG.verboseRetry = ~NoPrint;');
```

add:

```matlab
w.add('GDSGE_CFG.solNames = %s;', solNamesLit);
w.add('GDSGE_CFG.stateNames = %s;', stateNameCell);
w.add('GDSGE_CFG.stateGrids = %s;', stateGridCell);
w.add('GDSGE_CFG.diagnoseAt = DiagnoseMinorIter;');
w.add('GDSGE_CFG.errorOnNonconvergence = true;');
```

- [ ] **Step 3: Remove the dead post-return warning block**

Still in `emitIter.m`, delete these three `w.add` lines (keep the preceding
`NeedResolved = GDSGE_DIAG.needResolved;` line):

```matlab
w.add('if any(NeedResolved)');
w.add('    warning(''gdsge:runtime:unconverged'', ''%%s'', gdsge.runtime.reportUnconverged(NeedResolved, GDSGE_F, GDSGE_SIZE, %s, %s, 5));', ...
    stateNameCell, stateGridCell);
w.add('end');
```

- [ ] **Step 4: Regenerate the HL1996 snapshot and review the diff**

```
matlab -batch "run('tests/HeatonLucas1996/codegen/regen_snapshots.m')"
git --no-pager diff -- tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt
```

Expected diff, confined to: `MaxMinorIter = 20;` + new `DiagnoseMinorIter = 10;`
in the setup block; five new `GDSGE_CFG.*` lines (solNames/stateNames/stateGrids/
diagnoseAt/errorOnNonconvergence); the removed three-line `if any(NeedResolved)
warning(...)` block. `NeedResolved = GDSGE_DIAG.needResolved;`,
`IterRslt.NeedResolved = NeedResolved;`, and the `printIterProgress(... nnz(NeedResolved) ...)`
line must all still be present. The C++/compile goldens must be unchanged.

- [ ] **Step 5: Run the HL1996 codegen tests to verify the snapshot matches**

```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests({'tests/HeatonLucas1996/codegen/tSnapshotHL1996.m','tests/codegen/tEmitResultIter.m'}); assert(~any([r.Failed]))"
```

Expected: PASS (snapshot matches the regenerated golden; `IterRslt.NeedResolved`
guard still green).

- [ ] **Step 6: Commit**

```
git add src/+gdsge/+codegen/+mat/emitIter.m tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt
git commit -m "feat(codegen): thread diagnostic cfg into emitIter; drop dead warning block" -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Thread the diagnostic cfg into `emitSimulate`

Give simulate the `solNames` + `diagnoseAt` + `errorOnNonconvergence = false`
fields (no `stateNames`/`stateGrids`, so the formatter prints the per-component
signals but skips the grid-point table; the existing per-period warning is the
period locator).

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitSimulate.m`
- Snapshot: `tests/HeatonLucas1996/codegen/golden/simulate_HL1996_golden.txt` (regen)

- [ ] **Step 1: Add the solNames literal**

In `src/+gdsge/+codegen/+mat/emitSimulate.m`, after:

```matlab
stateNameCell = ['{''' strjoin(ir.states.names, ''',''') '''}'];
```

add:

```matlab
solComp     = gdsge.codegen.mat.solComponentNames(ir);
solNamesLit = ['{''' strjoin(solComp, ''',''') '''}'];
```

- [ ] **Step 2: Add the diagnostic cfg fields**

Still in `emitSimulate.m`, after:

```matlab
w.add('GDSGE_CFG.verboseRetry = false;');
```

add:

```matlab
w.add('GDSGE_CFG.solNames = %s;', solNamesLit);
w.add('GDSGE_CFG.diagnoseAt = DiagnoseMinorIter;');
w.add('GDSGE_CFG.errorOnNonconvergence = false;');
```

(Leave the existing per-period `if any(GDSGE_DIAG.needResolved) warning(...)`
block untouched.)

- [ ] **Step 3: Regenerate the snapshot and review the diff**

```
matlab -batch "run('tests/HeatonLucas1996/codegen/regen_snapshots.m')"
git --no-pager diff -- tests/HeatonLucas1996/codegen/golden/simulate_HL1996_golden.txt
```

Expected diff, confined to: `MaxMinorIter = 20;` + `DiagnoseMinorIter = 10;` in
the setup block, and three new `GDSGE_CFG.*` lines
(solNames/diagnoseAt/errorOnNonconvergence). The per-period warning line stays.

- [ ] **Step 4: Run the simulate snapshot test**

```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/HeatonLucas1996/codegen/tSnapshotHL1996.m'); assert(~any([r.Failed]))"
```

Expected: PASS.

- [ ] **Step 5: Commit**

```
git add src/+gdsge/+codegen/+mat/emitSimulate.m tests/HeatonLucas1996/codegen/golden/simulate_HL1996_golden.txt
git commit -m "feat(codegen): thread diagnostic cfg into emitSimulate (warn-and-continue)" -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: Full-suite backward-compat gate (the cap=20 risk check)

The end-to-end gates re-run `iter` and `simulate` for every corpus model under
the new `MaxMinorIter = 20`. The risk is **safe_assets** (multi-root grid point;
convergence depends on the RNG state at the restart — see the project memory on
safe_assets RNG-restart fragility): if any outer iteration needs more than 20
random restarts on a single point, the new hard error breaks a passing gate.

**Files:** none (validation only).

- [ ] **Step 1: Run the full suite**

```
matlab -batch "cd('tests'); run_tests"
```

Expected: exit 0; all tests pass (the prior count plus the new
`tSolComponentNames`, `tDiagnoseConvergence`, `tSolveProblemsCap`). Functional
`IterRslt`/`SimuRslt` goldens are bit-exact (converging models never enter the
new code paths).

- [ ] **Step 2: If (and only if) an end-to-end gate now fails with `gdsge:runtime:solveNotConverged`**

Diagnose before changing the default:

1. Identify the failing model and read the printed diagnostic block (it names
   the stuck component / dominant equation).
2. Confirm the gate pins `rng` before `iter` (the safe_assets/HL1996 gates do).
   An unpinned restart sequence is a test artifact, not a real >20-retry need.
3. Re-run that single gate a few times pinned. If it converges within 20 every
   time, the failure was RNG-order noise — re-pin and move on.
4. Only if a point genuinely needs >20 restarts every run: raise the default in
   `emitSetup.m` (`MaxMinorIter = 30;` say), update `tEmitSetup.m` to match,
   regen the HL1996 snapshots (`run('tests/HeatonLucas1996/codegen/regen_snapshots.m')`),
   and document the reason in this plan + the PROGRESS changelog. Do **not**
   silently widen the cap without recording why.

- [ ] **Step 3: Commit (only if Step 2 required a change)**

```
git add -A
git commit -m "fix(codegen): raise MaxMinorIter default to N (safe_assets needs >20 restarts)" -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: Update PROGRESS.md changelog

**Files:**
- Modify: `PROGRESS.md` (add a dated changelog entry under `## Changelog`)

- [ ] **Step 1: Prepend a changelog entry**

Add, as the newest entry under `## Changelog` in `PROGRESS.md`:

```markdown
- 2026-06-14: **Bounded retries + tiered convergence diagnostics (cartesian path).**
  The inner resolve loop in `gdsge.runtime.solveProblems` no longer retries random
  restarts forever: `MaxMinorIter` default `inf → 20`, a new `DiagnoseMinorIter`
  (default 10) tier-1 threshold, and a new pure formatter
  `gdsge.runtime.diagnoseConvergence` (NaN count; bounds pinned among unconverged
  per named policy component; near-bound among converged; dominant equation residual
  vs the others; worst-K grid points). On exhaustion the **iter** path prints the full
  block and throws `gdsge:runtime:solveNotConverged` (`cfg.errorOnNonconvergence=true`);
  the **simulate** path prints the same block but warns and continues (now bounded).
  All new behavior gated on `isfield(cfg,'diagnoseAt')`, so converging models are
  byte-identical (functional goldens unchanged; only the two HL1996 text snapshots
  regenerated). New helper `gdsge.codegen.mat.solComponentNames` maps SOL rows to
  named components, threaded into emitIter/emitSimulate cfg. ASG path + outer-loop
  `MaxIter` stalls remain out of scope. Spec
  `docs/superpowers/specs/2026-06-14-convergence-diagnostics-design.md`, plan
  `docs/superpowers/plans/2026-06-14-convergence-diagnostics.md`. Branch
  `convergence-diagnostics`.
```

- [ ] **Step 2: Run the full suite once more (PROGRESS.md is not code, but confirm nothing drifted)**

```
matlab -batch "cd('tests'); run_tests"
```

Expected: exit 0.

- [ ] **Step 3: Commit**

```
git add PROGRESS.md
git commit -m "docs: log convergence-diagnostics work in PROGRESS.md" -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-review (performed against the spec)

**Spec coverage:**
- Bounded retry (`MaxMinorIter inf→20`) → Task 4.
- Tier-1 at `DiagnoseMinorIter=10`, tier-2 at cap → Task 3 (solveProblems) + Task 4 (option).
- Hard error on iter / warn-continue on simulate (`cfg.errorOnNonconvergence`) → Task 3 + Tasks 5/6.
- Named per-component bounds (pinned + near) → Task 1 (names) + Task 2 (formatter).
- NaN signal → Task 2.
- Dominant equation vs others → Task 2.
- Worst-K point table (iter only) → Task 2 (gated on stateGrids) + Task 5 (iter supplies grids) / Task 6 (simulate omits them).
- `isfield(cfg,'diagnoseAt')` gate preserves bare-cfg behavior → Task 3 (`disabledWithoutDiagnoseAt`) + the untouched `tSolveProblems`.
- Snapshot regen confined to HL1996 text goldens; functional goldens bit-exact → Tasks 5/6/7.
- safe_assets cap risk → Task 7.

**Placeholder scan:** none — every code/test step contains complete code; every command has expected output.

**Type/name consistency:** `cfg.diagnoseAt`, `cfg.errorOnNonconvergence`,
`cfg.solNames`, `cfg.stateNames`, `cfg.stateGrids` are used identically across
solveProblems (Task 3), the formatter (Task 2), and the codegen threading
(Tasks 5/6). `gdsge.codegen.mat.solComponentNames` and
`gdsge.runtime.diagnoseConvergence` are defined (Tasks 1/2) before first use
(Tasks 3/5/6). `NeedResolved` is preserved everywhere it is read
(`IterRslt.NeedResolved`, `printIterProgress`); only the dead warning block is
removed (Task 5).
