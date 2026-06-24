# Throttled Resolve-Round Printing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Throttle the resolve-round progress line to print every `ResolvePrintFreq` rounds (default 100) and print a single end-of-loop resolve summary, uniformly across all three resolve sites.

**Architecture:** A new pure helper `gdsge.runtime.printResolveProgress` owns the window-crossing throttle logic and the summary text. A new `ResolvePrintFreq` option (default 100) is threaded through the parser/IR/codegen chain exactly like the existing `PrintFreq`. The two solver runtimes (`solveProblems.m`, `solveProblemsAsg.m`) call the helper instead of inline `fprintf`, gated by the existing `cfg.verboseRetry` master switch.

**Tech Stack:** MATLAB R2025b, `+gdsge` namespaced package, MATLAB unit-test framework (`matlab.unittest`).

Spec: `docs/superpowers/specs/2026-06-16-resolve-print-throttle-design.md`

Run tests: `matlab -batch "cd('tests'); run_tests"` (PowerShell 7 unavailable; junit.xml + exit code authoritative). Per project memory, **never run two MATLAB processes concurrently**.

---

## File Structure

- **Create** `src/+gdsge/+runtime/printResolveProgress.m` — the throttle+summary helper (one responsibility: format/gate resolve output).
- **Modify** `tests/runtime/tReporting.m` — add unit tests for the helper (sits beside the existing `printIterProgress` tests).
- **Modify** option-plumbing files (one line each):
  - `src/+gdsge/+parser/defaultSetupCode.m`
  - `src/+gdsge/+parser/resolveOptions.m`
  - `src/+gdsge/+ir/schema.m`
  - `src/+gdsge/+codegen/+mat/emitSetup.m`
  - `src/+gdsge/+codegen/+mat/optionsWhitelist.m`
  - `src/+gdsge/+codegen/+mat/emitIterInit.m`
  - `src/+gdsge/+codegen/+mat/emitIter.m`
  - `src/+gdsge/+codegen/+mat/emitIterAsg.m`
- **Modify** call sites:
  - `src/+gdsge/+runtime/solveProblems.m`
  - `src/+gdsge/+runtime/solveProblemsAsg.m`

---

## Task 1: The `printResolveProgress` helper (TDD)

**Files:**
- Create: `src/+gdsge/+runtime/printResolveProgress.m`
- Test: `tests/runtime/tReporting.m`

- [ ] **Step 1: Write the failing tests**

Add these methods inside the `methods (Test)` block of `tests/runtime/tReporting.m`, after the `returnsWhetherPrinted` method (before the `reportUnconverged` section):

```matlab
% ---- printResolveProgress -----------------------------------------
function resolveRoundSilentBetweenBoundaries(tc)
    % minorIter 50, step 1, freq 100 -> no boundary crossed
    out = evalc(['gdsge.runtime.printResolveProgress(' ...
        '''round'', '''', 50, 1, 3, 0.2, 100, true, false);']);
    tc.verifyEqual(out, '');
end
function resolveRoundPrintsOncePerWindow(tc)
    % minorIter 100, step 1, freq 100 -> crosses one boundary
    out = evalc(['gdsge.runtime.printResolveProgress(' ...
        '''round'', '''', 100, 1, 3, 0.2, 100, true, false);']);
    tc.verifyTrue(contains(out, 'resolve round 100'));
    tc.verifyTrue(contains(out, '3 unconverged'));
end
function resolveRoundBatchPrintsOnce(tc)
    % a single batch jumps 0 -> 250 (step 250), freq 100: print exactly once
    out = evalc(['gdsge.runtime.printResolveProgress(' ...
        '''round'', '''', 250, 250, 5, 0.9, 100, true, false);']);
    nlines = numel(strsplit(strtrim(out), char(10)));
    tc.verifyEqual(nlines, 1);
    tc.verifyTrue(contains(out, 'resolve round 250'));
end
function resolveMasterSwitchOff(tc)
    out = evalc(['gdsge.runtime.printResolveProgress(' ...
        '''round'', '''', 100, 1, 3, 0.2, 100, false, false);']);
    tc.verifyEqual(out, '');
    out2 = evalc(['gdsge.runtime.printResolveProgress(' ...
        '''summary'', '''', 5, 0, 0, 0, 100, false, false);']);
    tc.verifyEqual(out2, '');
end
function resolveSummaryConverged(tc)
    out = evalc(['gdsge.runtime.printResolveProgress(' ...
        '''summary'', '''', 7, 0, 0, 0, 100, true, false);']);
    tc.verifyTrue(contains(out, 'resolve converged after 7 rounds'));
end
function resolveSummaryNoneWhenZeroRounds(tc)
    out = evalc(['gdsge.runtime.printResolveProgress(' ...
        '''summary'', '''', 0, 0, 0, 0, 100, true, false);']);
    tc.verifyEqual(out, '');
end
function resolveSummaryCapHit(tc)
    out = evalc(['gdsge.runtime.printResolveProgress(' ...
        '''summary'', '''', 200, 0, 4, 1.5, 100, true, true);']);
    tc.verifyTrue(contains(out, 'resolve stopped at MaxMinorIter=200'));
    tc.verifyTrue(contains(out, '4 points still unconverged'));
    tc.verifyTrue(contains(out, 'worst residual 1.5'));
end
function resolveAsgLabelPrefix(tc)
    out = evalc(['gdsge.runtime.printResolveProgress(' ...
        '''round'', ''asg '', 100, 1, 3, 0.2, 100, true, false);']);
    tc.verifyTrue(contains(out, 'asg resolve round 100'));
end
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: FAIL — `Unrecognized function or variable 'gdsge.runtime.printResolveProgress'` (or "Undefined function") on the new `tReporting` methods.

- [ ] **Step 3: Write the helper**

Create `src/+gdsge/+runtime/printResolveProgress.m`:

```matlab
function printResolveProgress(mode, label, minorIter, step, nUnconverged, worstF, ...
        freq, verboseRetry, hitCap)
% PRINTRESOLVEPROGRESS  Throttled resolve-round line + end-of-loop summary.
%   mode = 'round'   : per-round line, printed only when a `freq` window
%                      boundary is crossed this round (step = the increment to
%                      minorIter this round; 1 for the +1 paths, batch for the
%                      in-MEX randomize path). Prints at most once per window
%                      even when a single batch spans several windows.
%   mode = 'summary' : one end-of-loop line. hitCap selects the wording:
%                      false -> "converged after N rounds"
%                      true  -> "stopped at MaxMinorIter=N: M still unconverged".
%   label  : '' for cartesian, 'asg ' for ASG (prefixes the message).
%   Master gate: prints nothing when verboseRetry is false.
if ~verboseRetry
    return;
end
switch mode
    case 'round'
        if floor(minorIter / freq) ~= floor((minorIter - step) / freq)
            fprintf('  %sresolve round %d: %d unconverged, worst residual %g\n', ...
                label, minorIter, nUnconverged, worstF);
        end
    case 'summary'
        if minorIter <= 0
            return;
        end
        if hitCap
            fprintf(['  %sresolve stopped at MaxMinorIter=%d: ' ...
                '%d points still unconverged, worst residual %g\n'], ...
                label, minorIter, nUnconverged, worstF);
        else
            fprintf('  %sresolve converged after %d rounds\n', label, minorIter);
        end
    otherwise
        error('gdsge:runtime:printResolveProgress:badMode', ...
            'mode must be ''round'' or ''summary'', got ''%s''.', mode);
end
end
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: PASS (all `tReporting` methods, including the 8 new ones; exit code 0).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+runtime/printResolveProgress.m tests/runtime/tReporting.m
git commit -m "feat(runtime): printResolveProgress helper — throttled round line + summary"
```

---

## Task 2: Thread the `ResolvePrintFreq` option (default 100)

No new test of its own — Task 3/4 exercise it end-to-end via the existing solver tests, and the option-default is covered by the full suite's codegen goldens. This task is pure plumbing mirroring `PrintFreq`.

**Files:**
- Modify: `src/+gdsge/+parser/defaultSetupCode.m`
- Modify: `src/+gdsge/+parser/resolveOptions.m:35`
- Modify: `src/+gdsge/+ir/schema.m:70`
- Modify: `src/+gdsge/+codegen/+mat/emitSetup.m:15`
- Modify: `src/+gdsge/+codegen/+mat/optionsWhitelist.m:6`
- Modify: `src/+gdsge/+codegen/+mat/emitIterInit.m:42`
- Modify: `src/+gdsge/+codegen/+mat/emitIter.m:189`
- Modify: `src/+gdsge/+codegen/+mat/emitIterAsg.m:265` and `:325`

- [ ] **Step 1: Add the parser default**

In `src/+gdsge/+parser/defaultSetupCode.m`, find the line:

```matlab
    'PrintFreq = 10;', ...
```

Add immediately after it:

```matlab
    'ResolvePrintFreq = 100;', ...
```

- [ ] **Step 2: Resolve the option value**

In `src/+gdsge/+parser/resolveOptions.m`, after line 35:

```matlab
o.printFreq     = getf(ws, 'PrintFreq', 10);
```

add:

```matlab
o.resolvePrintFreq = getf(ws, 'ResolvePrintFreq', 100);
```

- [ ] **Step 3: Register in the IR schema**

In `src/+gdsge/+ir/schema.m`, after line 70:

```matlab
    'printFreq',    opt(fScalar()), ...
```

add:

```matlab
    'resolvePrintFreq', opt(fScalar()), ...
```

- [ ] **Step 4: Emit the setup assignment**

In `src/+gdsge/+codegen/+mat/emitSetup.m`, after line 15:

```matlab
w.add('PrintFreq = %s;', mat2str(getOpt(opt, 'printFreq', 10), 17));
```

add:

```matlab
w.add('ResolvePrintFreq = %s;', mat2str(getOpt(opt, 'resolvePrintFreq', 100), 17));
```

- [ ] **Step 5: Whitelist the user-facing option**

In `src/+gdsge/+codegen/+mat/optionsWhitelist.m`, the `base` cell array (line 6) lists `'PrintFreq'`. Add `'ResolvePrintFreq'` to that same cell, e.g. change:

```matlab
base = {'TolEq','TolSol','TolFun','PrintFreq','NoPrint','SaveFreq','NoSave', ...
```

to:

```matlab
base = {'TolEq','TolSol','TolFun','PrintFreq','ResolvePrintFreq','NoPrint','SaveFreq','NoSave', ...
```

- [ ] **Step 6: Set `GDSGE_CFG.resolvePrintFreq` in the three iter emitters**

In each of `emitIterInit.m` (line 42), `emitIter.m` (line 189), and `emitIterAsg.m` (lines 265 **and** 325), find the line:

```matlab
w.add('GDSGE_CFG.verboseRetry = ~NoPrint;');
```

and add immediately after each occurrence:

```matlab
w.add('GDSGE_CFG.resolvePrintFreq = ResolvePrintFreq;');
```

(That is one insertion in `emitIterInit.m`, one in `emitIter.m`, and two in `emitIterAsg.m`. Do **not** touch `emitSimulate.m` / `emitSimulateAsg.m`, which set `verboseRetry = false` and stay silent.)

- [ ] **Step 7: Run the full suite to verify nothing regressed**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: PASS (exit 0). The codegen goldens now include the new `ResolvePrintFreq = 100;` setup line and `GDSGE_CFG.resolvePrintFreq` assignments. If golden-comparison tests fail purely on these added lines, update the affected golden files (they are regenerated expectations, not behavior changes) and re-run.

- [ ] **Step 8: Commit**

```bash
git add src/+gdsge/+parser/defaultSetupCode.m src/+gdsge/+parser/resolveOptions.m src/+gdsge/+ir/schema.m src/+gdsge/+codegen/+mat/emitSetup.m src/+gdsge/+codegen/+mat/optionsWhitelist.m src/+gdsge/+codegen/+mat/emitIterInit.m src/+gdsge/+codegen/+mat/emitIter.m src/+gdsge/+codegen/+mat/emitIterAsg.m tests/
git commit -m "feat(codegen): thread ResolvePrintFreq option (default 100)"
```

---

## Task 3: Wire the helper into `solveProblems.m` (cartesian)

**Files:**
- Modify: `src/+gdsge/+runtime/solveProblems.m`
- Test: `tests/runtime/tSolveProblems.m` (verify still green; bare-cfg path)

This file has three insertion points: the `useMexRandomize` branch (currently silent), the legacy `else` branch (currently prints every round), and one end-of-loop summary covering both.

- [ ] **Step 1: Add a guarded freq read near the top of the function**

In `src/+gdsge/+runtime/solveProblems.m`, just after `minorIter = 0;` (line 54), add:

```matlab
% Bare-cfg callers (unit tests) omit resolvePrintFreq; default to 100.
if isfield(cfg, 'resolvePrintFreq'); resolveFreq = cfg.resolvePrintFreq; else; resolveFreq = 100; end
```

- [ ] **Step 2: Add the throttled round line to the `useMexRandomize` branch**

In the `useMexRandomize` while-loop, after `trialOffset = trialOffset + batch;` (line 78) add:

```matlab
        gdsge.runtime.printResolveProgress('round', '', minorIter, batch, ...
            nnz((f > cfg.tolSol) | isnan(f)), max(f(:)), resolveFreq, cfg.verboseRetry, false);
```

- [ ] **Step 3: Replace the legacy per-round print**

In the legacy `else` branch, replace the existing block (lines 159-162):

```matlab
    if cfg.verboseRetry
        fprintf('  resolve round %d: %d unconverged, worst residual %g\n', ...
            minorIter, nnz((f > cfg.tolSol) | isnan(f)), max(f(:)));
    end
```

with a step-aware helper call. The increment this round is 1, or 2 when `adaptInSol` added an extra (line 155). Compute it from the values: replace the block with:

```matlab
    roundStep = 1 + double(~isempty(cfg.adaptInSol));
    gdsge.runtime.printResolveProgress('round', '', minorIter, roundStep, ...
        nnz((f > cfg.tolSol) | isnan(f)), max(f(:)), resolveFreq, cfg.verboseRetry, false);
```

- [ ] **Step 4: Add the end-of-loop summary**

After the line `end   % closes the `else` wrapping the legacy resolve loop ...` (line 174) and before the exhaustion diagnostics block (line 176 `if isfield(cfg,'diagnoseAt') ...`), add:

```matlab
gdsge.runtime.printResolveProgress('summary', '', minorIter, 0, ...
    nnz((f > cfg.tolSol) | isnan(f)), max(f(:)), resolveFreq, cfg.verboseRetry, ...
    (max(isnan(f)) || max(f(:)) > cfg.tolSol));
```

- [ ] **Step 5: Run the cartesian solver tests**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: PASS. `tSolveProblems`, `tInMexRandomizeDriver`, `tInMexResolveDriver`, `tSolveProblemsCap` stay green (they set `verboseRetry=false`, so the helper is silent and the bare-cfg `resolveFreq` default is exercised).

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+runtime/solveProblems.m
git commit -m "feat(runtime): throttled resolve round + summary in solveProblems"
```

---

## Task 4: Wire the helper into `solveProblemsAsg.m` (ASG)

**Files:**
- Modify: `src/+gdsge/+runtime/solveProblemsAsg.m`
- Test: `tests/runtime/tSolveProblemsAsg.m` (verify still green)

- [ ] **Step 1: Add a guarded freq read near the top**

In `src/+gdsge/+runtime/solveProblemsAsg.m`, just after `minorIter = 0;` (line 82), add:

```matlab
% Bare-cfg callers (unit tests) omit resolvePrintFreq; default to 100.
if isfield(cfg, 'resolvePrintFreq'); resolveFreq = cfg.resolvePrintFreq; else; resolveFreq = 100; end
```

- [ ] **Step 2: Replace the per-round print**

Replace the existing block (lines 103-106):

```matlab
    if cfg.verboseRetry
        fprintf('  asg resolve round %d: %d unconverged, worst residual %g\n', ...
            minorIter, nnz((f > cfg.tolSol) | isnan(f)), max(f(:)));
    end
```

with (ASG advances `minorIter` by 1 per round, so step = 1):

```matlab
    gdsge.runtime.printResolveProgress('round', 'asg ', minorIter, 1, ...
        nnz((f > cfg.tolSol) | isnan(f)), max(f(:)), resolveFreq, cfg.verboseRetry, false);
```

- [ ] **Step 3: Add the end-of-loop summary**

After the `while` loop ends (line 107 `end`) and before `diag = struct();` (line 109), add:

```matlab
gdsge.runtime.printResolveProgress('summary', 'asg ', minorIter, 0, ...
    nnz((f > cfg.tolSol) | isnan(f)), max(f(:)), resolveFreq, cfg.verboseRetry, ...
    (max(isnan(f)) || max(f(:)) > cfg.tolSol));
```

- [ ] **Step 4: Run the ASG solver tests**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: PASS. `tSolveProblemsAsg` stays green.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+runtime/solveProblemsAsg.m
git commit -m "feat(runtime): throttled resolve round + summary in solveProblemsAsg"
```

---

## Task 5: Full-suite verification

- [ ] **Step 1: Run the entire suite**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: PASS, exit 0. Confirm via `tests/results/junit.xml` (TAP file may carry stale lines — per project memory, junit.xml + exit code are authoritative).

- [ ] **Step 2: Manual smoke (optional but recommended)**

Pick one end-to-end model gate that drives many resolve rounds and confirm by eye that the resolve output is now throttled (a line roughly every 100 rounds) and ends with one summary line. If none is convenient, the unit tests in Task 1 plus the green suite are sufficient evidence.

- [ ] **Step 3: Final commit (if any golden/test tweaks remain uncommitted)**

```bash
git add -A
git commit -m "test: resolve-print throttle — full suite green"
```
