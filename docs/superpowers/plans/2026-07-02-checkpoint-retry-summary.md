# Checkpoint Retry Summary Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the per-major-iteration `[iter N] resolve converged after M rounds` console lines with one compact retry summary line printed at each `PrintFreq` checkpoint; non-convergence diagnostics/errors unchanged.

**Architecture:** A pure stats struct (`retryStatsInit` / `retryStatsAdd` / `retryStatsEndIter` in `+gdsge/+runtime`) accumulates retry rounds per major iteration; `printIterProgress` gains an optional stats argument that prints the summary line and returns the stats reset. `printResolveProgress`'s `'summary'` mode goes silent unless the `MaxMinorIter` cap was hit. The generated iter loops (cartesian + ASG) get three emitted lines each to thread the struct through.

**Tech Stack:** MATLAB R2025b (`C:\Program Files\MATLAB\R2025b\bin\matlab.exe`), `matlab.unittest`, existing codegen golden snapshots.

**Spec:** `docs/superpowers/specs/2026-07-02-checkpoint-retry-summary-design.md`

## Global Constraints

- Run MATLAB headless: `matlab -batch "<expr>"` from the repo root `D:\refactor_gdsge`. Exit code reflects errors.
- **One MATLAB process at a time** — never run two `matlab -batch` concurrently (OpenMP saturates all cores).
- **Never add a toolbox source to a persistent MATLAB path.** Tests add `src/` per process (`run_tests.m` does this; single-class commands below do it inline).
- Public API and `IterRslt`/`SimuRslt` struct shapes are **frozen**. No new gmod option, no IR/schema change, no model-JSON snapshot change.
- Simulate path untouched (it already runs `verboseRetry=false`).
- The non-convergence path is unchanged: `stopped at MaxMinorIter` line + FINAL diagnostics + `error` (iter) / warn-and-continue (simulate).
- TDD: failing test first. Small, frequent commits.
- Test report of record: `tests/results/junit.xml` + exit code (`results.tap` accumulates stale documents — ignore it).

## File Map

| File | Change |
|---|---|
| `src/+gdsge/+runtime/retryStatsInit.m` | create — zeroed accumulator |
| `src/+gdsge/+runtime/retryStatsAdd.m` | create — add one solve call's rounds to the current bucket |
| `src/+gdsge/+runtime/retryStatsEndIter.m` | create — close the bucket once per major iteration |
| `src/+gdsge/+runtime/printIterProgress.m` | modify — optional 9th arg prints/reset the summary line |
| `src/+gdsge/+runtime/printResolveProgress.m` | modify — `'summary'` mode silent unless `hitCap` |
| `src/+gdsge/+codegen/+mat/emitIter.m` | modify — thread `GDSGE_RETRY_STATS` through the cartesian loop |
| `src/+gdsge/+codegen/+mat/emitIterAsg.m` | modify — same for the ASG main refinement loop |
| `tests/runtime/tRetryStats.m` | create — unit tests for the three helpers |
| `tests/runtime/tReporting.m` | modify — printIterProgress stats tests; converged-summary now silent |
| `tests/runtime/tSolveProblems.m` | modify — verbose converged run prints nothing; verbose cap still prints |
| `tests/codegen/tEmitAsg.m` | modify — contains-assertions for the ASG wiring |
| `tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt` | regenerate via `regen_snapshots.m` |

---

### Task 1: Retry-stats accumulator helpers

**Files:**
- Create: `src/+gdsge/+runtime/retryStatsInit.m`
- Create: `src/+gdsge/+runtime/retryStatsAdd.m`
- Create: `src/+gdsge/+runtime/retryStatsEndIter.m`
- Test: `tests/runtime/tRetryStats.m` (new)

**Interfaces:**
- Consumes: nothing.
- Produces (used by Tasks 2, 4, 5):
  - `stats = gdsge.runtime.retryStatsInit()` → struct with double fields `itersSeen, itersRetried, totalRounds, worstIter, worstRounds, curRounds`, all `0`.
  - `stats = gdsge.runtime.retryStatsAdd(stats, minorIters)` — `[]` for `stats` initializes; adds `minorIters` (a solve call's `diag.minorIters`) into `curRounds`.
  - `stats = gdsge.runtime.retryStatsEndIter(stats, iter)` — closes the per-iteration bucket: `itersSeen+1`; if `curRounds > 0` also `itersRetried+1`, `totalRounds+curRounds`, worst-tracking (`worstIter`/`worstRounds` = the major iteration with the most rounds); clears `curRounds`.

- [ ] **Step 1: Write the failing test**

Create `tests/runtime/tRetryStats.m`:

```matlab
classdef tRetryStats < matlab.unittest.TestCase
    % Pure accumulator for per-major-iteration retry activity; summarized at
    % the PrintFreq checkpoint by printIterProgress (see the 2026-07-02
    % checkpoint-retry-summary spec).
    methods (Test)
        function initIsAllZero(tc)
            s = gdsge.runtime.retryStatsInit();
            tc.verifyEqual(s.itersSeen, 0);
            tc.verifyEqual(s.itersRetried, 0);
            tc.verifyEqual(s.totalRounds, 0);
            tc.verifyEqual(s.worstIter, 0);
            tc.verifyEqual(s.worstRounds, 0);
            tc.verifyEqual(s.curRounds, 0);
        end
        function addFromEmptyInitializes(tc)
            s = gdsge.runtime.retryStatsAdd([], 3);
            tc.verifyEqual(s.curRounds, 3);
            tc.verifyEqual(s.itersSeen, 0);   % bucket not closed yet
        end
        function multipleAddsAccumulateOneBucket(tc)
            % ASG calls Add once per refinement level within one major iter
            s = gdsge.runtime.retryStatsAdd([], 2);
            s = gdsge.runtime.retryStatsAdd(s, 5);
            s = gdsge.runtime.retryStatsEndIter(s, 7);
            tc.verifyEqual(s.itersSeen, 1);
            tc.verifyEqual(s.itersRetried, 1);
            tc.verifyEqual(s.totalRounds, 7);
            tc.verifyEqual(s.worstIter, 7);
            tc.verifyEqual(s.worstRounds, 7);
            tc.verifyEqual(s.curRounds, 0);   % bucket cleared
        end
        function zeroRetryIterCountsSeenOnly(tc)
            s = gdsge.runtime.retryStatsAdd([], 0);
            s = gdsge.runtime.retryStatsEndIter(s, 3);
            tc.verifyEqual(s.itersSeen, 1);
            tc.verifyEqual(s.itersRetried, 0);
            tc.verifyEqual(s.totalRounds, 0);
        end
        function worstTracksMaxAcrossIters(tc)
            s = gdsge.runtime.retryStatsInit();
            s = gdsge.runtime.retryStatsAdd(s, 4); s = gdsge.runtime.retryStatsEndIter(s, 1);
            s = gdsge.runtime.retryStatsAdd(s, 9); s = gdsge.runtime.retryStatsEndIter(s, 2);
            s = gdsge.runtime.retryStatsAdd(s, 6); s = gdsge.runtime.retryStatsEndIter(s, 3);
            tc.verifyEqual(s.itersSeen, 3);
            tc.verifyEqual(s.itersRetried, 3);
            tc.verifyEqual(s.totalRounds, 19);
            tc.verifyEqual(s.worstIter, 2);
            tc.verifyEqual(s.worstRounds, 9);
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (from repo root):
```
matlab -batch "addpath('src'); results = runtests('tests/runtime/tRetryStats.m'); disp(table(results)); exit(any([results.Failed]))"
```
Expected: FAIL / error — `Unrecognized ... gdsge.runtime.retryStatsInit`.

- [ ] **Step 3: Write the implementation**

Create `src/+gdsge/+runtime/retryStatsInit.m`:

```matlab
function stats = retryStatsInit()
% RETRYSTATSINIT  Zeroed retry-stats accumulator.
%   Tracks resolve-retry activity per major iteration between PrintFreq
%   checkpoints; see retryStatsAdd / retryStatsEndIter / printIterProgress.
stats = struct('itersSeen', 0, 'itersRetried', 0, 'totalRounds', 0, ...
    'worstIter', 0, 'worstRounds', 0, 'curRounds', 0);
end
```

Create `src/+gdsge/+runtime/retryStatsAdd.m`:

```matlab
function stats = retryStatsAdd(stats, minorIters)
% RETRYSTATSADD  Add one solve call's retry rounds to the current bucket.
%   Pass [] to start a fresh accumulator. Call once per solveProblems /
%   solveProblemsAsg call (the ASG refinement loop calls several times per
%   major iteration); retryStatsEndIter closes the bucket.
if isempty(stats); stats = gdsge.runtime.retryStatsInit(); end
stats.curRounds = stats.curRounds + minorIters;
end
```

Create `src/+gdsge/+runtime/retryStatsEndIter.m`:

```matlab
function stats = retryStatsEndIter(stats, iter)
% RETRYSTATSENDITER  Close the current major iteration's retry bucket.
%   Call exactly once per major iteration, after all solves and before the
%   printIterProgress checkpoint, so itersSeen counts iterations (not solve
%   calls) and the N/M checkpoint denominator is right.
if isempty(stats); stats = gdsge.runtime.retryStatsInit(); end
stats.itersSeen = stats.itersSeen + 1;
if stats.curRounds > 0
    stats.itersRetried = stats.itersRetried + 1;
    stats.totalRounds = stats.totalRounds + stats.curRounds;
    if stats.curRounds > stats.worstRounds
        stats.worstRounds = stats.curRounds;
        stats.worstIter = iter;
    end
end
stats.curRounds = 0;
end
```

- [ ] **Step 4: Run test to verify it passes**

Same command as Step 2. Expected: all 5 pass, exit 0.

- [ ] **Step 5: Commit**

```
git add src/+gdsge/+runtime/retryStatsInit.m src/+gdsge/+runtime/retryStatsAdd.m src/+gdsge/+runtime/retryStatsEndIter.m tests/runtime/tRetryStats.m
git commit -m "feat(runtime): retry-stats accumulator for checkpoint summaries"
```

---

### Task 2: `printIterProgress` optional stats argument

**Files:**
- Modify: `src/+gdsge/+runtime/printIterProgress.m` (whole file, currently 11 lines)
- Test: `tests/runtime/tReporting.m` (add methods after `returnsWhetherPrinted`, ~line 36)

**Interfaces:**
- Consumes: `gdsge.runtime.retryStatsInit()` and the stats struct fields from Task 1.
- Produces (used by Tasks 4, 5): `[printed, retryStats] = gdsge.runtime.printIterProgress(iter, metric, maxF, nUnconverged, elapsedSec, printFreq, noPrint, stopFlag, retryStats)`. The 9th argument is optional; all existing 8-arg call sites keep working. When the checkpoint line prints and `retryStats.itersRetried > 0`, a second line prints: `  resolve: %d/%d iters retried, %d rounds total, worst iter %d (%d rounds)`. Whenever the checkpoint prints, the returned stats are reset via `retryStatsInit()`.

- [ ] **Step 1: Write the failing tests**

In `tests/runtime/tReporting.m`, insert after the `returnsWhetherPrinted` method (line 36):

```matlab
        function iterStatsLinePrintsAndResets(tc)
            s = [];
            for it = 71:80
                s = gdsge.runtime.retryStatsAdd(s, double(it == 74)*58);
                s = gdsge.runtime.retryStatsEndIter(s, it);
            end
            [out, ~, s2] = evalc('gdsge.runtime.printIterProgress(80, 0.5, 1e-9, 0, 1, 10, 0, false, s)');
            tc.verifyTrue(contains(out, 'Iter:80'));
            tc.verifyTrue(contains(out, ...
                'resolve: 1/10 iters retried, 58 rounds total, worst iter 74 (58 rounds)'));
            tc.verifyEqual(s2.itersSeen, 0);        % reset after printing
            tc.verifyEqual(s2.itersRetried, 0);
        end
        function iterStatsSilentWhenNoRetries(tc)
            s = [];
            for it = 1:10
                s = gdsge.runtime.retryStatsAdd(s, 0);
                s = gdsge.runtime.retryStatsEndIter(s, it);
            end
            [out, ~, s2] = evalc('gdsge.runtime.printIterProgress(10, 0.5, 1e-9, 0, 1, 10, 0, false, s)');
            tc.verifyTrue(contains(out, 'Iter:10'));
            tc.verifyFalse(contains(out, 'resolve:'));
            tc.verifyEqual(s2.itersSeen, 0);        % still reset on print
        end
        function iterStatsCarriedWhenNotPrinting(tc)
            s = gdsge.runtime.retryStatsAdd([], 5);
            s = gdsge.runtime.retryStatsEndIter(s, 7);
            [out, ~, s2] = evalc('gdsge.runtime.printIterProgress(7, 0.5, 1e-9, 0, 1, 10, 0, false, s)');
            tc.verifyEqual(out, '');                % iter 7, freq 10: no print
            tc.verifyEqual(s2.itersRetried, 1);     % accumulates across window
            tc.verifyEqual(s2.totalRounds, 5);
        end
```

(The existing `printsOnFreqAndStop` / `silentOffFreqOrNoPrint` / `returnsWhetherPrinted` methods double as the 8-arg back-compat tests — leave them untouched.)

- [ ] **Step 2: Run tests to verify the new ones fail**

```
matlab -batch "addpath('src'); results = runtests('tests/runtime/tReporting.m'); disp(table(results)); exit(any([results.Failed]))"
```
Expected: the three new methods FAIL (`Too many input arguments` / `Too many output arguments`); all pre-existing methods pass.

- [ ] **Step 3: Write the implementation**

Replace the entire body of `src/+gdsge/+runtime/printIterProgress.m` with:

```matlab
function [printed, retryStats] = printIterProgress(iter, metric, maxF, nUnconverged, elapsedSec, printFreq, noPrint, stopFlag, retryStats)
% PRINTITERPROGRESS  One structured progress line, PrintFreq/NoPrint gated.
%   Returns true if it printed (caller resets its elapsed timer on true).
%   Optional retryStats (see retryStatsInit/Add/EndIter): when the checkpoint
%   prints and any covered iteration retried, a second indented "resolve:"
%   line summarizes the window; the returned stats are reset whenever the
%   checkpoint prints, so counts cover exactly the window since the last
%   printed checkpoint.
if nargin < 9; retryStats = []; end
printed = false;
if noPrint; return; end
if mod(iter, printFreq) == 0 || stopFlag
    fprintf('Iter:%d, Metric:%g, maxF:%g, unconverged:%d, elapsed:%.1fs\n', ...
        iter, metric, maxF, nUnconverged, elapsedSec);
    if ~isempty(retryStats)
        if retryStats.itersRetried > 0
            fprintf('  resolve: %d/%d iters retried, %d rounds total, worst iter %d (%d rounds)\n', ...
                retryStats.itersRetried, retryStats.itersSeen, retryStats.totalRounds, ...
                retryStats.worstIter, retryStats.worstRounds);
        end
        retryStats = gdsge.runtime.retryStatsInit();
    end
    printed = true;
end
end
```

- [ ] **Step 4: Run tests to verify they pass**

Same command as Step 2. Expected: all of `tReporting` passes, exit 0.

- [ ] **Step 5: Commit**

```
git add src/+gdsge/+runtime/printIterProgress.m tests/runtime/tReporting.m
git commit -m "feat(runtime): printIterProgress prints compact retry summary at checkpoints"
```

---

### Task 3: Silence the converged resolve summary

**Files:**
- Modify: `src/+gdsge/+runtime/printResolveProgress.m:25-36` (the `'summary'` case) and its header comment (lines 10-13)
- Test: `tests/runtime/tReporting.m` (`resolveSummaryConverged`, lines 71-76)
- Test: `tests/runtime/tSolveProblems.m` (two new methods; uses existing local fakes `fakeMexNarrow`, `fakeMexNever`, and `baseCfg`)

`solveProblems.m` and `solveProblemsAsg.m` need **no edits**: they already pass `hitCap` to the `'summary'` call; gating inside `printResolveProgress` makes the converged call a no-op (this is how the spec's "fires only when the cap was hit" is realized).

**Interfaces:**
- Consumes: nothing new.
- Produces: `printResolveProgress('summary', ...)` prints ONLY when `hitCap` is true (wording unchanged: `resolve stopped at MaxMinorIter=%d: %d points still unconverged, worst residual %g`); prints nothing when converged. `'round'` mode unchanged.

- [ ] **Step 1: Write the failing tests**

In `tests/runtime/tReporting.m`, replace the `resolveSummaryConverged` method (lines 71-76) with:

```matlab
        function resolveSummaryConvergedIsSilent(tc)
            % converged retries are summarized at the PrintFreq checkpoint
            % (printIterProgress + retryStats*), not per iteration
            out = evalc(['gdsge.runtime.printResolveProgress(' ...
                '''summary'', '''', 37, 7, 0, 0, 0, 100, true, false);']);
            tc.verifyEqual(out, '');
        end
```

In `tests/runtime/tSolveProblems.m`, add two methods inside `methods (Test)` (after `convergedPointsAreSkippedOnRestart`, line 85):

```matlab
        function verboseConvergedRetriesPrintNothing(tc)
            % retries happen (restarts needed) but console stays clean
            rng(42);
            n = 4;
            cfg = baseCfg([1 n]);
            cfg.useNearestNeighbor = false;
            cfg.verboseRetry = true;
            data = 0.5 * ones(1, n);
            sol = zeros(1, n);
            [out, ~, f] = evalc(['gdsge.runtime.solveProblems(@fakeMexNarrow, ' ...
                'sol, zeros(1,n), ones(1,n), data, 1e20*ones(1,n), ' ...
                'zeros(1,n), zeros(1,n), cfg)']);
            tc.verifyTrue(all(f <= cfg.tolSol));    % it did converge via restarts
            tc.verifyEqual(out, '');
        end
        function verboseCapStillPrintsStoppedLine(tc)
            n = 2;
            cfg = baseCfg([1 n]);
            cfg.useNearestNeighbor = false;
            cfg.maxMinorIter = 3;
            cfg.verboseRetry = true;
            out = evalc(['gdsge.runtime.solveProblems(@fakeMexNever, ' ...
                'zeros(1,n), zeros(1,n), ones(1,n), zeros(1,n), ' ...
                '1e20*ones(1,n), zeros(1,n), zeros(1,n), cfg);']);
            tc.verifyTrue(contains(out, 'resolve stopped at MaxMinorIter=3'));
            tc.verifyFalse(contains(out, 'converged after'));
        end
```

- [ ] **Step 2: Run tests to verify the new ones fail**

```
matlab -batch "addpath('src'); results = [runtests('tests/runtime/tReporting.m'), runtests('tests/runtime/tSolveProblems.m')]; disp(table(results)); exit(any([results.Failed]))"
```
Expected: `resolveSummaryConvergedIsSilent` and `verboseConvergedRetriesPrintNothing` FAIL (output contains `resolve converged after N rounds`); everything else passes.

- [ ] **Step 3: Write the implementation**

In `src/+gdsge/+runtime/printResolveProgress.m`, replace the `'summary'` case (lines 25-36) with:

```matlab
    case 'summary'
        % Converged retries print nothing here: they are summarized at the
        % PrintFreq checkpoint (printIterProgress + retryStats*). Only the
        % MaxMinorIter-exhaustion line remains a per-solve print.
        if ~hitCap
            return;
        end
        fprintf(['  [iter %d] %sresolve stopped at MaxMinorIter=%d: ' ...
            '%d points still unconverged, worst residual %g\n'], ...
            majorIter, label, minorIter, nUnconverged, worstF);
```

Also update the header comment lines 12-13 from

```
%   mode = 'summary' : one end-of-loop line. hitCap selects the wording:
%                      false -> "converged after N rounds"
%                      true  -> "stopped at MaxMinorIter=N: M still unconverged".
```

to

```
%   mode = 'summary' : end-of-loop line, printed ONLY when hitCap is true
%                      ("stopped at MaxMinorIter=N: M still unconverged");
%                      converged solves are summarized at the PrintFreq
%                      checkpoint instead (printIterProgress + retryStats*).
```

(The `if minorIter <= 0; return; end` guard becomes unreachable-but-harmless; delete it.)

- [ ] **Step 4: Run tests to verify they pass**

Same command as Step 2, plus the cap-path suite:
```
matlab -batch "addpath('src'); results = runtests('tests/runtime/tSolveProblemsCap.m'); disp(table(results)); exit(any([results.Failed]))"
```
Expected: all pass (tSolveProblemsCap uses `verboseRetry=false`, unaffected).

- [ ] **Step 5: Commit**

```
git add src/+gdsge/+runtime/printResolveProgress.m tests/runtime/tReporting.m tests/runtime/tSolveProblems.m
git commit -m "feat(runtime): per-iteration resolve summary prints only on MaxMinorIter exhaustion"
```

---

### Task 4: Wire the cartesian iter loop (`emitIter.m`) + regen golden

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitIter.m:204, 220-221, 250-252`
- Regenerate: `tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt` (via `tests/HeatonLucas1996/codegen/regen_snapshots.m`)
- Gate: `tests/HeatonLucas1996/codegen/tSnapshotHL1996.m` (no edits — equality against the golden)

**Interfaces:**
- Consumes: `retryStatsInit/Add/EndIter` (Task 1), 9-arg `printIterProgress` (Task 2).
- Produces: generated `iter_<model>.m` contains, verbatim:
  - `GDSGE_RETRY_STATS = gdsge.runtime.retryStatsInit();` (before the main `while(~stopFlag)`)
  - `GDSGE_RETRY_STATS = gdsge.runtime.retryStatsAdd(GDSGE_RETRY_STATS, GDSGE_DIAG.minorIters);` (after the solve)
  - `GDSGE_RETRY_STATS = gdsge.runtime.retryStatsEndIter(GDSGE_RETRY_STATS, GDSGE_Iter);` and a two-output `printIterProgress(..., GDSGE_RETRY_STATS)` call guarded by `if GDSGE_PRINTED` (at the checkpoint).

- [ ] **Step 1: Make the emitter changes**

In `src/+gdsge/+codegen/+mat/emitIter.m`:

(a) Before `w.add('stopFlag = false;');` (line 204), insert:

```matlab
w.add('GDSGE_RETRY_STATS = gdsge.runtime.retryStatsInit();');
```

(b) After `w.add('NeedResolved = GDSGE_DIAG.needResolved;');` (line 220), insert:

```matlab
w.add('GDSGE_RETRY_STATS = gdsge.runtime.retryStatsAdd(GDSGE_RETRY_STATS, GDSGE_DIAG.minorIters);');
```

(c) Replace the checkpoint block (lines 250-252)

```matlab
w.add('if gdsge.runtime.printIterProgress(GDSGE_Iter, GDSGE_Metric, max(GDSGE_F), nnz(NeedResolved), toc(GDSGE_timer), PrintFreq, NoPrint, stopFlag)');
w.add('    GDSGE_timer = tic;');
w.add('end');
```

with

```matlab
w.add('GDSGE_RETRY_STATS = gdsge.runtime.retryStatsEndIter(GDSGE_RETRY_STATS, GDSGE_Iter);');
w.add('[GDSGE_PRINTED, GDSGE_RETRY_STATS] = gdsge.runtime.printIterProgress(GDSGE_Iter, GDSGE_Metric, max(GDSGE_F), nnz(NeedResolved), toc(GDSGE_timer), PrintFreq, NoPrint, stopFlag, GDSGE_RETRY_STATS);');
w.add('if GDSGE_PRINTED');
w.add('    GDSGE_timer = tic;');
w.add('end');
```

Do NOT touch `emitIterInit.m` or the WarmUp segment — init solves have no checkpoint (spec: their converged chatter just disappears, which Task 3 already accomplished).

- [ ] **Step 2: Run the snapshot test to see the expected mismatch**

```
matlab -batch "addpath('src'); results = runtests('tests/HeatonLucas1996/codegen/tSnapshotHL1996.m'); disp(table(results)); exit(any([results.Failed]))"
```
Expected: FAIL — emitted text no longer equals the golden (this is the "failing test" for a snapshot-speced emitter).

- [ ] **Step 3: Regenerate the golden and review the diff**

```
matlab -batch "addpath('src'); cd('tests/HeatonLucas1996/codegen'); regen_snapshots"
git diff tests/HeatonLucas1996/codegen/golden/
```
Expected diff: ONLY `iter_HL1996_golden.txt` changes, and only the three insertions/replacement from Step 1 (plus no change to `simulate_HL1996_golden.txt`, `mex_HL1996_golden_cpp.txt`, `compile_HL1996_golden.txt`). If anything else changed, stop and investigate before committing.

- [ ] **Step 4: Run the snapshot + functional gates**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); results = [runtests('tests/HeatonLucas1996/codegen/tSnapshotHL1996.m'), runtests('tests/HeatonLucas1996/codegen/tFunctionalHL1996.m')]; disp(table(results)); exit(any([results.Failed]))"
```
Expected: PASS (the functional test executes the generated iter file, exercising the new emitted lines end to end).

- [ ] **Step 5: Commit**

```
git add src/+gdsge/+codegen/+mat/emitIter.m tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt
git commit -m "feat(codegen): cartesian iter loop records retries, summarizes at checkpoint"
```

---

### Task 5: Wire the ASG iter loop (`emitIterAsg.m`)

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitIterAsg.m:93, 117-118, 145-147`
- Test: `tests/codegen/tEmitAsg.m` (extend `iterAsgHasRefinementLoop`, lines 37-59)

**Interfaces:**
- Consumes: `retryStatsInit/Add/EndIter` (Task 1), 9-arg `printIterProgress` (Task 2).
- Produces: generated `iter_<model>.m` (ASG) with the same three wiring lines as Task 4, except the checkpoint call uses `nnz(~GDSGE_solved)` for the unconverged count.

**Placement rules (important):**
- The `retryStatsAdd` line goes in the MAIN refinement loop only — immediately after the `addProposeAndSolve(...)` call at `emitIterAsg.m:117-118`. Do NOT put it inside the `addProposeAndSolve` subfunction: that subfunction is also called at line 167 for the CONSTRUCT_OUTPUT pass, which runs AFTER the checkpoint (inside the SaveFreq block) and must not leak rounds into the next window.
- The ASG init segment (`addInitSegment`, solve at line 330) stays untouched.

- [ ] **Step 1: Write the failing test**

In `tests/codegen/tEmitAsg.m`, add to the `iterAsgHasRefinementLoop` method (after line 58's `GDSGE_OPTIONS_VALID` assertion):

```matlab
            % checkpoint retry summary wiring (2026-07-02 spec): init before
            % the loop, one Add per refinement-level solve, EndIter + 9-arg
            % printIterProgress at the checkpoint
            tc.verifyTrue(contains(txt, 'GDSGE_RETRY_STATS = gdsge.runtime.retryStatsInit();'));
            tc.verifyTrue(contains(txt, 'GDSGE_RETRY_STATS = gdsge.runtime.retryStatsAdd(GDSGE_RETRY_STATS, GDSGE_DIAG.minorIters);'));
            tc.verifyTrue(contains(txt, 'GDSGE_RETRY_STATS = gdsge.runtime.retryStatsEndIter(GDSGE_RETRY_STATS, GDSGE_Iter);'));
            tc.verifyTrue(contains(txt, '[GDSGE_PRINTED, GDSGE_RETRY_STATS] = gdsge.runtime.printIterProgress(GDSGE_Iter, GDSGE_Metric, max(GDSGE_F), nnz(~GDSGE_solved), toc(GDSGE_timer), PrintFreq, NoPrint, stopFlag, GDSGE_RETRY_STATS);'));
            % exactly ONE Add call site: the CONSTRUCT_OUTPUT re-solve (after
            % the checkpoint) must not accumulate into the next window
            tc.verifyEqual(numel(strfind(txt, 'retryStatsAdd')), 1);
```

- [ ] **Step 2: Run test to verify it fails**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); results = runtests('tests/codegen/tEmitAsg.m'); disp(table(results)); exit(any([results.Failed]))"
```
Expected: `iterAsgHasRefinementLoop` FAILS on the new assertions; other methods pass.

- [ ] **Step 3: Make the emitter changes**

In `src/+gdsge/+codegen/+mat/emitIterAsg.m`:

(a) Before `w.add('stopFlag = false;');` (line 93), insert:

```matlab
w.add('GDSGE_RETRY_STATS = gdsge.runtime.retryStatsInit();');
```

(b) After the main-loop `addProposeAndSolve` call (lines 117-118):

```matlab
addProposeAndSolve(w, main, pack, m, ...
    'GDSGE_SOL_ASG_INTERP_NEW', 'GDSGE_SOL_ASG_INTERP', 8);
```

insert:

```matlab
w.add('GDSGE_RETRY_STATS = gdsge.runtime.retryStatsAdd(GDSGE_RETRY_STATS, GDSGE_DIAG.minorIters);');
```

(c) Replace the checkpoint block (lines 145-147)

```matlab
w.add('if gdsge.runtime.printIterProgress(GDSGE_Iter, GDSGE_Metric, max(GDSGE_F), nnz(~GDSGE_solved), toc(GDSGE_timer), PrintFreq, NoPrint, stopFlag)');
w.add('    GDSGE_timer = tic;');
w.add('end');
```

with

```matlab
w.add('GDSGE_RETRY_STATS = gdsge.runtime.retryStatsEndIter(GDSGE_RETRY_STATS, GDSGE_Iter);');
w.add('[GDSGE_PRINTED, GDSGE_RETRY_STATS] = gdsge.runtime.printIterProgress(GDSGE_Iter, GDSGE_Metric, max(GDSGE_F), nnz(~GDSGE_solved), toc(GDSGE_timer), PrintFreq, NoPrint, stopFlag, GDSGE_RETRY_STATS);');
w.add('if GDSGE_PRINTED');
w.add('    GDSGE_timer = tic;');
w.add('end');
```

- [ ] **Step 4: Run test to verify it passes**

Same command as Step 2. Expected: all of `tEmitAsg` passes, exit 0.

- [ ] **Step 5: Commit**

```
git add src/+gdsge/+codegen/+mat/emitIterAsg.m tests/codegen/tEmitAsg.m
git commit -m "feat(codegen): ASG iter loop records retries, summarizes at checkpoint"
```

---

### Task 6: Full-suite verification + Mendoza output eyeball

**Files:**
- No source changes. Gate: `tests/results/junit.xml` + exit code; visual check of the Mendoza console output.

**Interfaces:**
- Consumes: everything above.
- Produces: green suite; confirmation that the Mendoza example now prints compactly.

- [ ] **Step 1: Run the full suite**

```
matlab -batch "cd('tests'); run_tests"
```
Expected: exit 0. Takes a few minutes (~300 tests; the end-to-end gates are 5-18s each). Check `tests/results/junit.xml` for failures — NOT `results.tap` (it accumulates stale documents across runs).

- [ ] **Step 2: Eyeball the Mendoza end-to-end output**

The Mendoza model is the case that motivated this change; its end-to-end test runs `iter_mendoza2010` with default printing.

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); results = runtests('tests/Mendoza2010/codegen/tEndToEndMendoza2010.m'); disp(table(results)); exit(any([results.Failed]))"
```
Expected: PASS, and in the captured console output:
- NO `resolve converged after` lines anywhere;
- checkpoint lines (`Iter:10, Metric:...`) optionally followed by a single indented `  resolve: N/M iters retried, R rounds total, worst iter K (W rounds)` line;
- round-heartbeat lines (`resolve round 100:`) only if some single solve ran ≥100 rounds (unlikely; absence is fine).

If iterations converge without retries in this fixture, the `resolve:` line legitimately never appears — the `Iter:` lines with no interleaved resolve chatter are the pass signal.

- [ ] **Step 3: Commit (only if anything changed)**

Nothing should need committing here; if a test surfaced a fix, commit it with a `fix:` message referencing the failing test.

---

## Self-Review Notes

- Spec coverage: behavior (Tasks 2-5), components 1-5 (Tasks 1-5), goldens (Task 4), edge cases (denominator → Task 1 EndIter tests; skipped checkpoints → Task 2 `iterStatsCarriedWhenNotPrinting`; non-convergence unchanged → Task 3 cap tests + untouched `tSolveProblemsCap`), testing section (Tasks 1-3 unit, Task 3 evalc fake-MEX, Task 4 golden regen, Task 6 suite + eyeball). Spec's "solveProblems/solveProblemsAsg call fires only on hitCap" is realized by gating inside `printResolveProgress` with zero caller churn — behavior identical.
- The spec names `retryStatsAdd([] , ...)` as the initializer; `retryStatsInit` is the shared implementation of that initializer (Add still accepts `[]`), used directly by the emitters and by `printIterProgress` for the reset.
