# TechBubbles Retry Cap Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `TechBubblesV4nblogKRAv3.gmod` run end-to-end on the new toolbox and match the old toolbox's results, by removing the resolve-retry cap and closing the three defects that removing it exposes.

**Architecture:** Four independent source changes plus a new corpus gate. `emitSetup` stops emitting a finite `MaxMinorIter`; `solveProblems` replaces its `mod`-based diagnostics heartbeat with geometric backoff; the simulate path gets a sample-shaped unconverged reporter and a per-checkpoint retry aggregate; `ensurePath` gains a once-per-session source-tree integrity check. Nothing touches the IR schema, the solver, or the C++ backends.

**Tech Stack:** MATLAB R2025b, `matlab.unittest`, the `+gdsge` package (`gdsge.codegen.mat.*` emitters over a unit-tested `gdsge.runtime.*` library).

**Spec:** `docs/superpowers/specs/2026-07-31-techbubbles-retry-cap-design.md`

## Global Constraints

- **Branch:** all work on `fix/techbubbles-retry-cap` (already created; the spec commit is `babaf85`).
- **MATLAB path rule:** never add a toolbox source to a persistent path. Each run `addpath`s **exactly one** source — old `base_package/gdsge/source` **or** new `src/` — in its **own `matlab -batch` process**, restoring the path on cleanup. See `CLAUDE.md`.
- **One MATLAB process at a time.** Each saturates all cores via OpenMP; concurrent runs contend and distort timings.
- **Run the suite with:** `matlab -batch "cd('tests'); run_tests"` (PowerShell 7 is not installed, so `tests/run.ps1` will not run). Exit 0 = all pass. **`tests/results/junit.xml` and the exit code are authoritative — `results.tap` accumulates stale documents across runs.**
- **No IR schema change.** `MaxMinorIter` and `DiagnoseMinorIter` are generated-file defaults, not IR options. No existing model's `ir/*.gdsge.json` snapshot changes.
- **TDD:** write the failing test first, watch it fail, then implement. Small, frequent commits.
- **Commit trailers:** end every commit message with

  ```
  Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_011HSDe48haPWHSrWuTrGanS
  ```

---

### Task 1: Install self-check (`checkInstall`)

Turns `MATLAB:scriptNotAFunction` on a truncated source file into one actionable error. Independent of every other task.

**Files:**
- Create: `src/+gdsge/+runtime/checkInstall.m`
- Modify: `src/+gdsge/+runtime/ensurePath.m` (add the call at the top of the function body)
- Modify: `src/+gdsge/+runtime/Contents.m` (add the one-line index entry)
- Test: `tests/runtime/tCheckInstall.m` (create)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `gdsge.runtime.checkInstall(srcDir)` — `srcDir` optional. With no argument it derives the tree from its own location and short-circuits after the first call in a MATLAB session; with an argument it always scans (this is how tests point it at a fixture). Returns nothing; raises `gdsge:install:corrupt`.

- [ ] **Step 1: Write the failing test**

Create `tests/runtime/tCheckInstall.m`:

```matlab
classdef tCheckInstall < matlab.unittest.TestCase
    % Guards against the failure mode reported from a macOS install: a
    % truncated .m file makes MATLAB classify it as a script, and calling it
    % dies with MATLAB:scriptNotAFunction naming an internal file.
    methods (Test)
        function healthyTreePasses(tc)
            root = makeTree(tc, struct('good', 'function y = good(x)\ny = x;\nend\n'));
            tc.verifyWarningFree(@() gdsge.runtime.checkInstall(root));
        end
        function emptyFileIsReported(tc)
            root = makeTree(tc, struct('good', 'function good()\nend\n', 'broken', ''));
            err = tc.verifyError(@() gdsge.runtime.checkInstall(root), 'gdsge:install:corrupt');
            tc.verifyTrue(contains(err.message, 'broken.m'));
            tc.verifyFalse(contains(err.message, 'good.m'));
        end
        function commentOnlyFileIsReported(tc)
            root = makeTree(tc, struct('broken', '% just a comment\n'));
            tc.verifyError(@() gdsge.runtime.checkInstall(root), 'gdsge:install:corrupt');
        end
        function contentsMIsExempt(tc)
            % Contents.m files are comment-only by design.
            root = makeTree(tc, struct('good', 'function good()\nend\n'));
            fid = fopen(fullfile(root, '+gdsge', 'Contents.m'), 'w');
            fprintf(fid, '%% package index\n'); fclose(fid);
            tc.verifyWarningFree(@() gdsge.runtime.checkInstall(root));
        end
        function classdefIsAccepted(tc)
            root = makeTree(tc, struct('good', 'classdef good < handle\nend\n'));
            tc.verifyWarningFree(@() gdsge.runtime.checkInstall(root));
        end
        function allBadFilesReportedInOneError(tc)
            root = makeTree(tc, struct('b1', '', 'b2', ''));
            err = tc.verifyError(@() gdsge.runtime.checkInstall(root), 'gdsge:install:corrupt');
            tc.verifyTrue(contains(err.message, 'b1.m'));
            tc.verifyTrue(contains(err.message, 'b2.m'));
        end
        function scansNestedPackageFolders(tc)
            root = makeTree(tc, struct());
            nested = fullfile(root, '+gdsge', '+runtime');
            mkdir(nested);
            fid = fopen(fullfile(nested, 'deep.m'), 'w'); fclose(fid);   % empty
            err = tc.verifyError(@() gdsge.runtime.checkInstall(root), 'gdsge:install:corrupt');
            tc.verifyTrue(contains(err.message, 'deep.m'));
        end
    end
end

function root = makeTree(tc, files)
% Build a throwaway src-shaped tree: <root>/+gdsge/<name>.m for each field.
root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
mkdir(fullfile(root, '+gdsge'));
names = fieldnames(files);
for i = 1:numel(names)
    fid = fopen(fullfile(root, '+gdsge', [names{i} '.m']), 'w');
    fprintf(fid, '%s', sprintf(files.(names{i})));
    fclose(fid);
end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "addpath('src'); r = runtests('tests/runtime/tCheckInstall.m'); disp(table(r)); exit(any([r.Failed]))"`

Expected: FAIL — `Unable to resolve the name 'gdsge.runtime.checkInstall'`.

- [ ] **Step 3: Write the implementation**

Create `src/+gdsge/+runtime/checkInstall.m`:

```matlab
function checkInstall(srcDir)
% CHECKINSTALL  Verify every toolbox source file actually declares a function.
%   A truncated or empty .m file makes MATLAB classify it as a script; calling
%   it then fails with MATLAB:scriptNotAFunction naming an internal file, which
%   tells the user nothing about what is wrong or how to fix it. This turns
%   that into one actionable error naming every damaged file.
%
%   checkInstall()        scan the installed tree, once per MATLAB session.
%   checkInstall(srcDir)  scan srcDir unconditionally (tests).
%
%   Contents.m is exempt: those files are comment-only by design. The scan
%   costs a few milliseconds over ~70 small files and runs from ensurePath,
%   ahead of any generated iter_/simulate_ work.
persistent checked
if nargin < 1
    if ~isempty(checked); return; end
    checked = true;    % set before scanning: a failure must not re-scan forever
    srcDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % -> src/
end

files = [ dir(fullfile(srcDir, '*.m'));
          dir(fullfile(srcDir, 'kernels', '*.m'));
          dir(fullfile(srcDir, '+gdsge', '**', '*.m')) ];

bad = {};
for i = 1:numel(files)
    if strcmp(files(i).name, 'Contents.m'); continue; end
    p = fullfile(files(i).folder, files(i).name);
    txt = fileread(p);
    if isempty(regexp(txt, '^\s*(function|classdef)\b', 'once', 'lineanchors'))
        bad{end+1} = p; %#ok<AGROW>
    end
end

if ~isempty(bad)
    error('gdsge:install:corrupt', ...
        ['GDSGE install looks incomplete or corrupted - these files declare no ' ...
         'function:\n  %s\n' ...
         'Restore them (e.g. `git checkout -- src/`) or re-clone the toolbox.'], ...
        strjoin(bad, sprintf('\n  ')));
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "addpath('src'); r = runtests('tests/runtime/tCheckInstall.m'); disp(table(r)); exit(any([r.Failed]))"`

Expected: PASS, 7 tests.

- [ ] **Step 5: Wire it into `ensurePath`**

In `src/+gdsge/+runtime/ensurePath.m`, insert as the first executable line of the function body — immediately after the header comment block ending `%   Idempotent. Called first by gdsge_codegen and by every generated`/`%   iter_<model>/simulate_<model>, ahead of any kernel call.` and before `kernelsDir = fullfile( ...`:

```matlab
gdsge.runtime.checkInstall();   % once per session; clear error on a damaged tree
```

- [ ] **Step 6: Add the `Contents.m` index entry**

In `src/+gdsge/+runtime/Contents.m`, add a line alongside the existing entries (keep the column alignment used by its neighbours):

```matlab
%   checkInstall         - verify the source tree declares functions (install guard)
```

- [ ] **Step 7: Run the runtime test folder to confirm nothing regressed**

Run: `matlab -batch "addpath('src'); r = runtests('tests/runtime'); disp(table(r)); exit(any([r.Failed]))"`

Expected: PASS. `tEnsurePath` must still pass — `checkInstall` scans the real tree and finds it healthy.

- [ ] **Step 8: Commit**

```bash
git add src/+gdsge/+runtime/checkInstall.m src/+gdsge/+runtime/ensurePath.m \
        src/+gdsge/+runtime/Contents.m tests/runtime/tCheckInstall.m
git commit -m "feat(runtime): install self-check for truncated source files

A truncated .m file makes MATLAB classify it as a script; calling it dies
with MATLAB:scriptNotAFunction naming an internal file. checkInstall runs
once per session from ensurePath and reports every damaged file with the
remedy.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_011HSDe48haPWHSrWuTrGanS"
```

---

### Task 2: Uncap `MaxMinorIter`; silence the simulate heartbeat by default

The root-cause fix. This model needs up to ~5,300 restarts per simulate period; the 200 cap silently truncated results.

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitSetup.m:23-24`
- Modify: `tests/codegen/tEmitSetup.m:31-32` and the `simulateTargetUsesIrSimulate` test
- Regenerate: `tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt`, `.../simulate_HL1996_golden.txt`

**Interfaces:**
- Consumes: nothing.
- Produces: generated `iter_*.m`/`simulate_*.m` carry `MaxMinorIter = inf;`; `simulate_*.m` additionally carries `DiagnoseMinorIter = 0;` while `iter_*.m` keeps `DiagnoseMinorIter = 10;`.

- [ ] **Step 1: Write the failing test**

In `tests/codegen/tEmitSetup.m`, change the two assertions in `frozenDefaultsPresent` (currently lines 31-32) to:

```matlab
            tc.verifyEqual(ws.MaxMinorIter, inf);
            tc.verifyEqual(ws.DiagnoseMinorIter, 10);      % iter target
```

and add this test method after `simulateTargetUsesIrSimulate`:

```matlab
        function simulateTargetSilencesRetryHeartbeat(tc)
            % Heavy retrying is normal in simulate (TechBubblesV4nblogKRAv3
            % needs ~5300 restarts/period), so the per-solve diagnostics
            % heartbeat is off by default there; the SimuPrintFreq aggregate
            % replaces it. The field is still emitted, so the FINAL
            % cap-exhaustion block still prints.
            ir = buildHL1996IR();
            ws = evalSetupText(gdsge.codegen.mat.emitSetup(ir, 'simulate'));
            tc.verifyEqual(ws.DiagnoseMinorIter, 0);
            tc.verifyEqual(ws.MaxMinorIter, inf);
        end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "addpath('src'); r = runtests('tests/codegen/tEmitSetup.m'); disp(table(r)); exit(any([r.Failed]))"`

Expected: FAIL — `frozenDefaultsPresent` reports 200 ~= Inf, and `simulateTargetSilencesRetryHeartbeat` reports 10 ~= 0.

- [ ] **Step 3: Write the implementation**

In `src/+gdsge/+codegen/+mat/emitSetup.m`, replace lines 23-24:

```matlab
w.add('MaxMinorIter = 200;');        % bounded retries; was inf (could loop forever)
w.add('DiagnoseMinorIter = 10;');    % print a diagnostic heartbeat every N retries while stuck
```

with:

```matlab
% Old-toolbox parity. A finite cap does not fail loudly: the simulate path
% writes the unconverged solution into SimuRslt and feeds it forward as the
% next period's state, so results are silently truncated. TechBubblesV4nblogKRAv3
% legitimately needs ~5300 randomized restarts per period. Users who want a
% runaway guard set MaxMinorIter themselves, which also arms
% errorOnNonconvergence on the iter path.
w.add('MaxMinorIter = inf;');
if strcmp(target, 'simulate')
    % Retrying in bulk is normal here, so a per-solve heartbeat would fire
    % every period. printSimuRetryStats aggregates at the SimuPrintFreq
    % checkpoint instead. The field is still emitted, so the FINAL
    % cap-exhaustion diagnostics block still prints.
    w.add('DiagnoseMinorIter = 0;');
else
    w.add('DiagnoseMinorIter = 10;');   % heartbeat every N retries while stuck
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "addpath('src'); r = runtests('tests/codegen/tEmitSetup.m'); disp(table(r)); exit(any([r.Failed]))"`

Expected: PASS.

- [ ] **Step 5: Regenerate the HL1996 codegen snapshots**

Run: `matlab -batch "addpath('src'); addpath('tests/HeatonLucas1996/codegen'); regen_snapshots"`

Then review the diff — the snapshot **is** the spec of the output:

```bash
git diff tests/HeatonLucas1996/codegen/golden/
```

Expected: exactly two files changed. `iter_HL1996_golden.txt` line 17 `MaxMinorIter = 200;` → `inf;`. `simulate_HL1996_golden.txt` line 17 the same, plus line 18 `DiagnoseMinorIter = 10;` → `0;`. **Nothing else.** If any other line moved, stop and investigate before committing.

- [ ] **Step 6: Run the codegen test folder**

Run: `matlab -batch "addpath('src'); r = runtests('tests/codegen'); disp(table(r)); exit(any([r.Failed]))"`

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/+gdsge/+codegen/+mat/emitSetup.m tests/codegen/tEmitSetup.m \
        tests/HeatonLucas1996/codegen/golden/
git commit -m "fix(codegen): MaxMinorIter defaults to inf again (old-toolbox parity)

The 200 cap did not fail loudly - the simulate path wrote the unconverged
solution into SimuRslt and fed it forward, silently truncating results.
TechBubblesV4nblogKRAv3 needs ~5300 randomized restarts per period and
produced garbage (q reached 987 against a true max of 2.5).

Simulate also defaults DiagnoseMinorIter to 0: bulk retrying is normal
there and the per-solve heartbeat would fire every period.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_011HSDe48haPWHSrWuTrGanS"
```

---

### Task 3: Geometric backoff for the diagnostics heartbeat

With the cap gone, `mod(minorIter, diagnoseAt) == 0` would print ~530 blocks per period. The `mod` test is also wrong for the in-MEX branch, which advances `minorIter` in batches and can jump across a multiple without landing on it.

**Files:**
- Modify: `src/+gdsge/+runtime/solveProblems.m` (lines 54-59 area, 93-98, 178-183, and a new local subfunction at the end)
- Modify: `tests/runtime/tSolveProblemsCap.m` (`recurringDiagnosticsEveryInterval`)
- Test: `tests/runtime/tSolveProblemsCap.m` (add two tests)

**Interfaces:**
- Consumes: nothing.
- Produces: heartbeat fires at `diagnoseAt`, then at each doubling (10, 20, 40, 80, …). `diagnoseAt <= 0` disables it entirely. Behaviour on cap exhaustion (the FINAL block, `errorOnNonconvergence`) is unchanged.

- [ ] **Step 1: Write the failing tests**

In `tests/runtime/tSolveProblemsCap.m`, **replace** the `recurringDiagnosticsEveryInterval` method with:

```matlab
        function diagnosticsBackOffGeometrically(tc)
            % diagnoseAt=2, cap=7: fires at 2 then 4 (next threshold 8 is never
            % reached). The cap iteration is covered by the FINAL block, not a
            % heartbeat. Backoff keeps an uncapped solve to ~log2(N) blocks.
            cfg = capCfg();
            cfg.maxMinorIter = 7;
            cfg.errorOnNonconvergence = false;
            out = evalc('runNever(cfg);');
            tc.verifyEqual(numel(strfind(out, 'continuing):')), 2);
            tc.verifyTrue(contains(out, '(retry 2,'));
            tc.verifyTrue(contains(out, '(retry 4,'));
            tc.verifyTrue(contains(out, 'FINAL convergence diagnostics'));
        end
        function diagnosticsDoubleAcrossAWideRun(tc)
            % 2, 4, 8, 16 over a 20-retry run: four blocks, not ten.
            cfg = capCfg();
            cfg.maxMinorIter = 20;
            cfg.errorOnNonconvergence = false;
            out = evalc('runNever(cfg);');
            tc.verifyEqual(numel(strfind(out, 'continuing):')), 4);
            tc.verifyTrue(contains(out, '(retry 16,'));
        end
        function diagnosticsOffWhenDiagnoseAtIsZero(tc)
            % Characterization: emitSetup emits DiagnoseMinorIter = 0 for
            % simulate, and the field is still present, so the FINAL block must
            % still print. Passes before and after the backoff change - it pins
            % behaviour Task 2 now depends on but nothing asserted before.
            cfg = capCfg();
            cfg.diagnoseAt = 0;
            cfg.errorOnNonconvergence = false;
            out = evalc('runNever(cfg);');
            tc.verifyEqual(numel(strfind(out, 'continuing):')), 0);
            tc.verifyTrue(contains(out, 'FINAL convergence diagnostics'));
        end
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `matlab -batch "addpath('src'); r = runtests('tests/runtime/tSolveProblemsCap.m'); disp(table(r)); exit(any([r.Failed]))"`

Expected: 2 of the 3 FAIL. `diagnosticsBackOffGeometrically` sees 3 blocks (the current cadence prints at 2, 4, 6) and `diagnosticsDoubleAcrossAWideRun` sees 9 (2, 4, …, 18). `diagnosticsOffWhenDiagnoseAtIsZero` already passes — it is a characterization test pinning behaviour Task 2 now relies on.

- [ ] **Step 3: Add the backoff state**

In `src/+gdsge/+runtime/solveProblems.m`, after line 59 (`if isfield(cfg, 'majorIter'); majorIter = cfg.majorIter; else; majorIter = 0; end`) add:

```matlab
% Diagnostics heartbeat: geometric backoff. Bare-cfg callers (unit tests) omit
% diagnoseAt entirely, which disables it.
if isfield(cfg, 'diagnoseAt'); diagnoseAt = cfg.diagnoseAt; else; diagnoseAt = 0; end
nextDiagnose = diagnoseAt;
```

- [ ] **Step 4: Replace both heartbeat blocks with the helper call**

In the `useMexRandomize` branch, replace lines 93-98:

```matlab
        % Periodic diagnostics heartbeat (now at batch granularity).
        if isfield(cfg,'diagnoseAt') && cfg.diagnoseAt > 0 ...
                && mod(minorIter, cfg.diagnoseAt) == 0 && minorIter < cfg.maxMinorIter ...
                && (max(isnan(f)) || max(f(:)) > cfg.tolSol)
            fprintf('[gdsge] convergence diagnostics (retry %d, continuing):\n%s\n', ...
                minorIter, gdsge.runtime.diagnoseConvergence(sol, lb, ub, f, eqVal, cfg, 'summary'));
        end
```

with:

```matlab
        nextDiagnose = heartbeat(nextDiagnose, diagnoseAt, minorIter, cfg, sol, lb, ub, f, eqVal);
```

In the legacy branch, replace lines 175-183 (the comment block plus the `if`):

```matlab
    % Periodic diagnostics: a summary heartbeat every cfg.diagnoseAt retries
    % while still stuck. The final block prints after the loop, so skip the cap
    % iteration here to avoid a duplicate at minorIter == maxMinorIter.
    if isfield(cfg,'diagnoseAt') && cfg.diagnoseAt > 0 ...
            && mod(minorIter, cfg.diagnoseAt) == 0 && minorIter < cfg.maxMinorIter ...
            && (max(isnan(f)) || max(f(:)) > cfg.tolSol)
        fprintf('[gdsge] convergence diagnostics (retry %d, continuing):\n%s\n', ...
            minorIter, gdsge.runtime.diagnoseConvergence(sol, lb, ub, f, eqVal, cfg, 'summary'));
    end
```

with:

```matlab
    nextDiagnose = heartbeat(nextDiagnose, diagnoseAt, minorIter, cfg, sol, lb, ub, f, eqVal);
```

- [ ] **Step 5: Add the local subfunction**

At the end of `src/+gdsge/+runtime/solveProblems.m`, after the final `end` of the main function, append:

```matlab
function nextDiagnose = heartbeat(nextDiagnose, diagnoseAt, minorIter, cfg, sol, lb, ub, f, eqVal)
% Diagnostics heartbeat with geometric backoff: fires at diagnoseAt, then at
% each doubling. Threshold-crossing rather than mod(), because the in-MEX
% randomize branch advances minorIter by a whole batch and would otherwise
% jump clean over a trigger. The cap iteration is left to the FINAL block, so
% there is no duplicate at minorIter == maxMinorIter.
%
% NOTE: this is a print-only helper. It must never call the model MEX - the
% MEX reads its inputs from ITS CALLER's workspace, so every mexFn call has to
% stay in solveProblems' own body (see the CONTRACT comment at the top).
if diagnoseAt <= 0 || minorIter < nextDiagnose || minorIter >= cfg.maxMinorIter
    return;
end
if ~(max(isnan(f)) || max(f(:)) > cfg.tolSol)
    return;                       % everything converged; nothing to report
end
fprintf('[gdsge] convergence diagnostics (retry %d, continuing):\n%s\n', ...
    minorIter, gdsge.runtime.diagnoseConvergence(sol, lb, ub, f, eqVal, cfg, 'summary'));
while nextDiagnose <= minorIter
    nextDiagnose = nextDiagnose * 2;
end
end
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `matlab -batch "addpath('src'); r = runtests('tests/runtime'); disp(table(r)); exit(any([r.Failed]))"`

Expected: PASS, including `tSolveProblems`, `tSolveProblemsAsg`, `tInMexResolveDriver`, and `tInMexRandomizeDriver`.

- [ ] **Step 7: Commit**

```bash
git add src/+gdsge/+runtime/solveProblems.m tests/runtime/tSolveProblemsCap.m
git commit -m "fix(runtime): geometric backoff for the diagnostics heartbeat

mod(minorIter, diagnoseAt) printed ~530 blocks per period once the retry
cap was lifted, and the in-MEX randomize branch advances minorIter by a
whole batch, so a step could jump clean over a trigger. Threshold-crossing
with doubling fires ~log2(N) times and cannot be skipped.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_011HSDe48haPWHSrWuTrGanS"
```

---

### Task 4: Sample-shaped unconverged reporter

`reportUnconverged` maps a linear index through `ind2sub(probSize, ...)` onto grid subscripts. Simulate passes `probSize = [1 GDSGE_NPROB]` with two state names, so the loop asks for `subs{3}` and dies with `MATLAB:badsubscript`. The simulate warn-and-continue path crashes whenever it fires.

**Files:**
- Create: `src/+gdsge/+runtime/reportUnconvergedSamples.m`
- Modify: `src/+gdsge/+runtime/Contents.m` (index entry)
- Modify: `src/+gdsge/+codegen/+mat/emitSimulate.m:156-159`
- Test: `tests/runtime/tReporting.m` (add tests after `reportAllConverged`)

**Interfaces:**
- Consumes: nothing.
- Produces: `msg = gdsge.runtime.reportUnconvergedSamples(needResolved, f, shockVec, stateNames, stateVals, topK)`. `needResolved` is a logical row over samples; `shockVec` is the per-sample shock index; `stateVals` is a cell of per-sample state vectors aligned with `stateNames`; all indexed by the linear problem (= sample) index. Returns `'All problems converged.'` when nothing needs resolving. `reportUnconverged` is unchanged.

- [ ] **Step 1: Write the failing test**

In `tests/runtime/tReporting.m`, add after the `reportAllConverged` method:

```matlab
        % ---- reportUnconvergedSamples --------------------------------------
        function reportSamplesNamesWorstPoints(tc)
            % Simulate solves one problem per sample; coordinates are the
            % sample's own shock/state values, not grid subscripts.
            need  = logical([0 1 0 1]);
            f     = [0 3 0 9];
            shock = [1 2 1 2];
            logK  = [0.1 0.2 0.3 0.4];
            b     = [-1 -2 -3 -4];
            msg = gdsge.runtime.reportUnconvergedSamples(need, f, shock, ...
                {'logK','b'}, {logK, b}, 5);
            tc.verifyTrue(contains(msg, '2 of 4'));
            tc.verifyTrue(contains(msg, 'sample=4'));      % worst residual first
            tc.verifyTrue(contains(msg, 'shock=2'));
            tc.verifyTrue(contains(msg, 'logK=0.4'));
            tc.verifyTrue(contains(msg, 'b=-4'));
            tc.verifyTrue(contains(msg, 'residual=9'));
            % worst first: sample 4 must appear before sample 2
            tc.verifyLessThan(strfind(msg, 'sample=4'), strfind(msg, 'sample=2'));
        end
        function reportSamplesTwoStatesDoesNotThrow(tc)
            % Regression: the old call reused reportUnconverged with
            % probSize=[1 N], which indexed subs{3} and threw
            % MATLAB:badsubscript for any two-state model.
            need = true(1,3); f = [1 2 3];
            tc.verifyWarningFree(@() gdsge.runtime.reportUnconvergedSamples( ...
                need, f, [1 1 1], {'a','b'}, {[1 2 3], [4 5 6]}, 5));
        end
        function reportSamplesHonoursTopK(tc)
            need = true(1,4); f = [4 3 2 1];
            msg = gdsge.runtime.reportUnconvergedSamples(need, f, ones(1,4), ...
                {'x'}, {1:4}, 2);
            tc.verifyTrue(contains(msg, '4 of 4 problems unconverged; worst 2:'));
            tc.verifyEqual(numel(strfind(msg, 'residual=')), 2);
        end
        function reportSamplesAllConverged(tc)
            msg = gdsge.runtime.reportUnconvergedSamples(false(1,4), zeros(1,4), ...
                ones(1,4), {'x'}, {1:4}, 5);
            tc.verifyEqual(msg, 'All problems converged.');
        end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "addpath('src'); r = runtests('tests/runtime/tReporting.m'); disp(table(r)); exit(any([r.Failed]))"`

Expected: FAIL — `Unable to resolve the name 'gdsge.runtime.reportUnconvergedSamples'`.

- [ ] **Step 3: Write the implementation**

Create `src/+gdsge/+runtime/reportUnconvergedSamples.m`:

```matlab
function msg = reportUnconvergedSamples(needResolved, f, shockVec, stateNames, stateVals, topK)
% REPORTUNCONVERGEDSAMPLES  Human-readable report of unconverged simulate samples.
%   The simulate path solves one problem per sample, so a problem's coordinates
%   are that sample's own current shock and state values - not grid subscripts,
%   which is reportUnconverged's job. Everything here is indexed by the linear
%   problem (= sample) index.
%
%   needResolved : logical row over samples
%   f            : residual per sample
%   shockVec     : per-sample shock index
%   stateNames   : cell of state names
%   stateVals    : cell of per-sample state vectors, aligned with stateNames
%   Returns the report text; the caller decides warning vs fprintf.
idx = find(needResolved);
if isempty(idx)
    msg = 'All problems converged.';
    return;
end
[~, order] = sort(f(idx), 'descend');
idx = idx(order);
k = min(topK, numel(idx));
lines = cell(1, k+1);
lines{1} = sprintf('%d of %d problems unconverged; worst %d:', ...
    nnz(needResolved), numel(needResolved), k);
for i = 1:k
    j = idx(i);
    parts = sprintf('sample=%d, shock=%d', j, shockVec(j));
    for s = 1:numel(stateNames)
        parts = sprintf('%s, %s=%g', parts, stateNames{s}, stateVals{s}(j));
    end
    lines{i+1} = sprintf('  [%s] residual=%g', parts, f(j));
end
msg = strjoin(lines, newline);
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "addpath('src'); r = runtests('tests/runtime/tReporting.m'); disp(table(r)); exit(any([r.Failed]))"`

Expected: PASS.

- [ ] **Step 5: Switch the simulate emitter to the new reporter**

In `src/+gdsge/+codegen/+mat/emitSimulate.m`, replace lines 156-159:

```matlab
w.add('if any(GDSGE_DIAG.needResolved)');
w.add('    warning(''gdsge:runtime:unconverged'', ''period %%d: %%s'', GDSGE_t, gdsge.runtime.reportUnconverged(GDSGE_DIAG.needResolved, GDSGE_F, [1 GDSGE_NPROB], %s, {%s}, 5));', ...
    stateNameCell, stateList);
w.add('end');
```

with:

```matlab
% Sample-shaped report: problems are samples here, so coordinates come from the
% per-sample `shock`/state vectors assigned at the top of the period loop, not
% from ind2sub over a grid (which threw MATLAB:badsubscript for >1 state).
w.add('if any(GDSGE_DIAG.needResolved)');
w.add('    warning(''gdsge:runtime:unconverged'', ''period %%d: %%s'', GDSGE_t, gdsge.runtime.reportUnconvergedSamples(GDSGE_DIAG.needResolved, GDSGE_F, shock, %s, {%s}, 5));', ...
    stateNameCell, stateList);
w.add('end');
```

- [ ] **Step 6: Add the `Contents.m` index entry**

In `src/+gdsge/+runtime/Contents.m`, next to the existing `reportUnconverged` line:

```matlab
%   reportUnconvergedSamples - same, for the per-sample simulate path
```

- [ ] **Step 7: Regenerate the HL1996 simulate snapshot and review**

Run: `matlab -batch "addpath('src'); addpath('tests/HeatonLucas1996/codegen'); regen_snapshots"`

```bash
git diff tests/HeatonLucas1996/codegen/golden/
```

Expected: only `simulate_HL1996_golden.txt`, only the `reportUnconverged(... [1 GDSGE_NPROB] ...)` call becoming `reportUnconvergedSamples(... shock ...)`.

- [ ] **Step 8: Run the runtime and codegen folders**

Run: `matlab -batch "addpath('src'); r = runtests({'tests/runtime','tests/codegen'}); disp(table(r)); exit(any([r.Failed]))"`

Expected: PASS.

- [ ] **Step 9: Commit**

```bash
git add src/+gdsge/+runtime/reportUnconvergedSamples.m src/+gdsge/+runtime/Contents.m \
        src/+gdsge/+codegen/+mat/emitSimulate.m tests/runtime/tReporting.m \
        tests/HeatonLucas1996/codegen/golden/
git commit -m "fix(runtime): sample-shaped unconverged report for simulate

emitSimulate passed probSize=[1 NPROB] with the model's state names, so
reportUnconverged indexed subs{3} and threw MATLAB:badsubscript for any
model with more than one state - the warn-and-continue path crashed
whenever it fired.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_011HSDe48haPWHSrWuTrGanS"
```

---

### Task 5: Per-checkpoint retry aggregate for simulate

Replaces the per-period heartbeat that Task 2 turned off with one aggregate line at the existing `SimuPrintFreq` checkpoint. `retryStatsInit`/`Add`/`EndIter` are already generic counters — "iter" is only a label — so they are reused verbatim.

**Files:**
- Create: `src/+gdsge/+runtime/printSimuRetryStats.m`
- Modify: `src/+gdsge/+runtime/Contents.m` (index entry)
- Modify: `src/+gdsge/+codegen/+mat/emitSimulate.m` (four insertions; anchors below)
- Test: `tests/runtime/tRetryStats.m` (add tests)

**Interfaces:**
- Consumes: `gdsge.runtime.retryStatsInit()`, `retryStatsAdd(stats, minorIters)`, `retryStatsEndIter(stats, iter)` — all existing, unchanged.
- Produces: `stats = gdsge.runtime.printSimuRetryStats(stats)` — prints one line when `stats.itersRetried > 0`, and **always** returns a fresh accumulator so counts cover exactly the window since the last checkpoint (same contract as `printIterProgress`'s stats handling).

- [ ] **Step 1: Write the failing test**

In `tests/runtime/tRetryStats.m`, add these methods inside the existing `methods (Test)` block:

```matlab
        function simuStatsLinePrintsAndResets(tc)
            s = [];
            for period = 1:10
                s = gdsge.runtime.retryStatsAdd(s, double(period == 3)*5253);
                s = gdsge.runtime.retryStatsEndIter(s, period);
            end
            [out, s2] = evalc('gdsge.runtime.printSimuRetryStats(s)');
            tc.verifyTrue(contains(out, '1/10 periods retried'));
            tc.verifyTrue(contains(out, '5253 rounds total'));
            tc.verifyTrue(contains(out, 'worst period 3 (5253 rounds)'));
            tc.verifyEqual(s2.itersSeen, 0);          % reset for the next window
            tc.verifyEqual(s2.itersRetried, 0);
        end
        function simuStatsSilentWhenNoRetries(tc)
            s = [];
            for period = 1:10
                s = gdsge.runtime.retryStatsAdd(s, 0);
                s = gdsge.runtime.retryStatsEndIter(s, period);
            end
            [out, s2] = evalc('gdsge.runtime.printSimuRetryStats(s)');
            tc.verifyEqual(out, '');
            tc.verifyEqual(s2.itersSeen, 0);
        end
        function simuStatsHandlesEmptyAccumulator(tc)
            [out, s2] = evalc('gdsge.runtime.printSimuRetryStats([])');
            tc.verifyEqual(out, '');
            tc.verifyEqual(s2.itersSeen, 0);
        end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "addpath('src'); r = runtests('tests/runtime/tRetryStats.m'); disp(table(r)); exit(any([r.Failed]))"`

Expected: FAIL — `Unable to resolve the name 'gdsge.runtime.printSimuRetryStats'`.

- [ ] **Step 3: Write the implementation**

Create `src/+gdsge/+runtime/printSimuRetryStats.m`:

```matlab
function stats = printSimuRetryStats(stats)
% PRINTSIMURETRYSTATS  One aggregate resolve line at a SimuPrintFreq checkpoint.
%   The simulate counterpart of the "resolve:" line printIterProgress emits for
%   the iteration loop. Models where bulk retrying is normal (TechBubblesV4nblogKRAv3
%   needs ~5300 restarts on its hardest periods) would otherwise report either
%   nothing at all or one diagnostics block per period.
%
%   Prints only when some period in the window retried, and ALWAYS returns a
%   fresh accumulator, so counts cover exactly the window since the last
%   checkpoint. Feed it the stats built by retryStatsInit/Add/EndIter.
if isempty(stats)
    stats = gdsge.runtime.retryStatsInit();
    return;
end
if stats.itersRetried > 0
    fprintf('  resolve: %d/%d periods retried, %d rounds total, worst period %d (%d rounds)\n', ...
        stats.itersRetried, stats.itersSeen, stats.totalRounds, ...
        stats.worstIter, stats.worstRounds);
end
stats = gdsge.runtime.retryStatsInit();
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "addpath('src'); r = runtests('tests/runtime/tRetryStats.m'); disp(table(r)); exit(any([r.Failed]))"`

Expected: PASS.

- [ ] **Step 5: Wire the accumulator into the simulate emitter**

Four insertions in `src/+gdsge/+codegen/+mat/emitSimulate.m`. Apply them **bottom-up** so earlier line numbers stay valid.

**(d)** Inside the `SimuPrintFreq` block, between line 178 (`w.add('    fprintf(''elapsed:%%.1fs\\n'', toc(GDSGE_timer));');`) and line 179 (`w.add('    GDSGE_timer = tic;');`), insert:

```matlab
w.add('    GDSGE_RETRY_STATS = gdsge.runtime.printSimuRetryStats(GDSGE_RETRY_STATS);');
```

**(c)** Immediately before line 167 (`w.add('if mod(GDSGE_t,SimuPrintFreq)==0');`), insert:

```matlab
w.add('GDSGE_RETRY_STATS = gdsge.runtime.retryStatsEndIter(GDSGE_RETRY_STATS, GDSGE_t);');
```

**(b)** Immediately after line 159 (the `w.add('end');` that closes the unconverged-warning block from Task 4), insert:

```matlab
w.add('GDSGE_RETRY_STATS = gdsge.runtime.retryStatsAdd(GDSGE_RETRY_STATS, GDSGE_DIAG.minorIters);');
```

**(a)** Between line 96 (`w.add('GDSGE_SHOCK_VAR_INDEX_BASE = ...');`) and line 97 (`w.add('GDSGE_timer = tic;');`), insert:

```matlab
w.add('GDSGE_RETRY_STATS = gdsge.runtime.retryStatsInit();');
```

- [ ] **Step 6: Regenerate the HL1996 simulate snapshot and review**

Run: `matlab -batch "addpath('src'); addpath('tests/HeatonLucas1996/codegen'); regen_snapshots"`

```bash
git diff tests/HeatonLucas1996/codegen/golden/
```

Expected: only `simulate_HL1996_golden.txt`, exactly four added lines, correctly indented (the two inside the loop body at 4 spaces, the checkpoint one at 8).

- [ ] **Step 7: Run the runtime and codegen folders**

Run: `matlab -batch "addpath('src'); r = runtests({'tests/runtime','tests/codegen'}); disp(table(r)); exit(any([r.Failed]))"`

Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git add src/+gdsge/+runtime/printSimuRetryStats.m src/+gdsge/+runtime/Contents.m \
        src/+gdsge/+codegen/+mat/emitSimulate.m tests/runtime/tRetryStats.m \
        tests/HeatonLucas1996/codegen/golden/
git commit -m "feat(runtime): aggregate simulate retry stats at the SimuPrintFreq checkpoint

Replaces the per-period diagnostics heartbeat (off by default since the
retry cap was lifted) with one line per checkpoint, reusing the generic
retryStats accumulator the iteration loop already uses.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_011HSDe48haPWHSrWuTrGanS"
```

---

### Task 6: Capture the TechBubbles golden from the old toolbox

Goldens come from the **old** toolbox, in its own MATLAB process with only `base_package/gdsge/source` on the path.

**Files:**
- Create: `tests/TechBubblesV4nblogKRAv3/TechBubblesV4nblogKRAv3.gmod` (copy)
- Create: `tests/golden/capture_TechBubblesV4nblogKRAv3.m`
- Create: `tests/TechBubblesV4nblogKRAv3/golden/IterRslt.mat`, `.../golden/SimuRslt.mat`

**Interfaces:**
- Consumes: nothing from earlier tasks (the old toolbox is independent of all of them).
- Produces: `golden/IterRslt.mat` (`IterRslt.Iter == 355`, `Metric < 1e-6`) and `golden/SimuRslt.mat` (fields `shock logK b q qd ER r qn IK`; `shock` is `6 x 101`, the rest `6 x 100`), consumed by Task 7.

- [ ] **Step 1: Copy the model into the test corpus**

```bash
mkdir -p tests/TechBubblesV4nblogKRAv3/golden tests/TechBubblesV4nblogKRAv3/ir \
         tests/TechBubblesV4nblogKRAv3/parser tests/TechBubblesV4nblogKRAv3/codegen
cp base_package/gdsge/tests/TechBubblesV4nblogKRAv3/TechBubblesV4nblogKRAv3.gmod \
   tests/TechBubblesV4nblogKRAv3/
```

- [ ] **Step 2: Write the capture script**

Create `tests/golden/capture_TechBubblesV4nblogKRAv3.m`:

```matlab
function capture_TechBubblesV4nblogKRAv3()
% Golden capture for TechBubblesV4nblogKRAv3 from the OLD toolbox.
% One-source-per-process: addpath base_package/gdsge/source only, controlled
% cd, restore on cleanup.
%
% This model is the corpus's only >200-retry case: its simulate needs up to
% ~5300 randomized restarts per period (the old toolbox's retry loop is
% unbounded). It also has multi-root grid points, so the randomized restart is
% order-sensitive - pin rng before iter AND before simulate, with the same seed
% the end-to-end gate uses.
%
% The gmod asks for 100 samples x 10000 periods (~7h). Reduce to 6 x 100 for a
% gate: measured 21s, versus 17.7s for 6x20 - cost is dominated by the hard
% early periods, so the longer run is nearly free.
here     = fileparts(mfilename('fullpath'));            % tests/golden
repoRoot = fileparts(fileparts(here));
oldSrc   = fullfile(repoRoot, 'base_package', 'gdsge', 'source');
modelDir = fullfile(repoRoot, 'tests', 'TechBubblesV4nblogKRAv3');
goldenDir = fullfile(modelDir, 'golden');

work = tempname; mkdir(work);
copyfile(fullfile(modelDir, 'TechBubblesV4nblogKRAv3.gmod'), work);

oldPath = path; restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
oldCd = pwd; cdRestore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
addpath(oldSrc);
cd(work);

t0 = tic;
gdsge_codegen('TechBubblesV4nblogKRAv3');
fprintf('codegen+compile wall-clock: %.1fs\n', toc(t0));

t1 = tic;
rng(0823);
IterRslt = iter_TechBubblesV4nblogKRAv3();
fprintf('iter wall-clock: %.1fs\n', toc(t1));
fprintf('Iter=%d  Metric=%g\n', IterRslt.Iter, IterRslt.Metric);

t2 = tic;
rng(0823);
SimuRslt = simulate_TechBubblesV4nblogKRAv3(IterRslt, ...
    struct('num_periods', 100, 'num_samples', 6));
fprintf('simulate wall-clock: %.1fs\n', toc(t2));

save(fullfile(goldenDir, 'IterRslt.mat'), 'IterRslt', '-v7');
save(fullfile(goldenDir, 'SimuRslt.mat'), 'SimuRslt', '-v7');
fprintf('SimuRslt fields: %s\n', strjoin(fieldnames(SimuRslt)', ', '));
fprintf('GOLDEN SAVED to %s\n', goldenDir);
end
```

- [ ] **Step 3: Run the capture**

Run: `matlab -batch "addpath('tests/golden'); capture_TechBubblesV4nblogKRAv3"`

Expected output includes `Iter=355  Metric=9.74867e-07` and `GOLDEN SAVED`. Iter takes ~8 s, simulate ~20-45 s.

**Note:** the old toolbox's MEX may raise `0xc0000005` (access violation) at process teardown when clearing the MEX — this is the known `vcomp` unload crash in the *old* toolbox and is harmless **provided `GOLDEN SAVED` printed first**. If it printed, the goldens are on disk and valid.

- [ ] **Step 4: Verify the goldens loaded and have the expected shape**

Run:

```bash
matlab -batch "g=load('tests/TechBubblesV4nblogKRAv3/golden/IterRslt.mat'); \
  fprintf('Iter=%d Metric=%g states=%d,%d\n', g.IterRslt.Iter, g.IterRslt.Metric, \
  numel(g.IterRslt.var_state.logK), numel(g.IterRslt.var_state.b)); \
  s=load('tests/TechBubblesV4nblogKRAv3/golden/SimuRslt.mat'); \
  disp(size(s.SimuRslt.shock)); disp(size(s.SimuRslt.q));"
```

Expected: `Iter=355 Metric=9.74867e-07 states=10,10`, then `6 101` and `6 100`.

- [ ] **Step 5: Check the golden file sizes**

```bash
ls -lh tests/TechBubblesV4nblogKRAv3/golden/
```

Expected: both well under 1 MB (the grid is 10x10x2). If either exceeds ~10 MB, stop — something captured more than intended.

- [ ] **Step 6: Commit**

```bash
git add tests/TechBubblesV4nblogKRAv3/TechBubblesV4nblogKRAv3.gmod \
        tests/TechBubblesV4nblogKRAv3/golden/ \
        tests/golden/capture_TechBubblesV4nblogKRAv3.m
git commit -m "test(golden): capture TechBubblesV4nblogKRAv3 from the old toolbox

Iter=355, Metric=9.75e-07; simulate reduced to 6x100 (the gmod asks for
100x10000, ~7h). The corpus's only >200-retry model.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_011HSDe48haPWHSrWuTrGanS"
```

---

### Task 7: Front-end test, IR snapshot, and the end-to-end gate

The regression lock. Without this, Task 2 can silently revert.

**Files:**
- Create: `tests/TechBubblesV4nblogKRAv3/parser/tFrontEndTechBubbles.m`
- Create: `tests/TechBubblesV4nblogKRAv3/ir/TechBubblesV4nblogKRAv3.gdsge.json`
- Create: `tests/TechBubblesV4nblogKRAv3/codegen/tEndToEndTechBubbles.m`
- Create: `tests/TechBubblesV4nblogKRAv3/tGoldenTechBubbles.m`

**Interfaces:**
- Consumes: Tasks 2-5 (the uncapped default and the two simulate-path fixes must be in place or the gate fails); Task 6's `golden/IterRslt.mat` and `golden/SimuRslt.mat`.
- Produces: nothing consumed downstream. `tests/ir/tRoundtrip.m` auto-discovers `**/ir/*.gdsge.json`, so the new snapshot is round-tripped automatically with no change to that test.

- [ ] **Step 1: Write the golden-shape test**

Create `tests/TechBubblesV4nblogKRAv3/tGoldenTechBubbles.m`:

```matlab
classdef tGoldenTechBubbles < matlab.unittest.TestCase
    % Golden files exist, load, and have the expected shape.
    % Captured from the OLD toolbox: Iter=355, Metric=9.74867e-07.
    methods (Test)
        function goldenLoadsAndHasShape(tc)
            here = fileparts(mfilename('fullpath'));
            g = load(fullfile(here, 'golden', 'IterRslt.mat'));
            R = g.IterRslt;
            tc.verifyEqual(R.shock_num, 2);
            tc.verifyEqual(numel(R.var_state.logK), 10);
            tc.verifyEqual(numel(R.var_state.b), 10);
            tc.verifyEqual(R.Iter, 355);
            tc.verifyLessThan(R.Metric, 1e-6);
            tc.verifyTrue(isfield(R.var_policy, 'bp'));       % array policy [2]
            tc.verifyTrue(isfield(R.var_aux, 'ER'));
            s = load(fullfile(here, 'golden', 'SimuRslt.mat'));
            tc.verifyEqual(size(s.SimuRslt.shock), [6 101]);  % samples x (periods+1)
            tc.verifyEqual(size(s.SimuRslt.q), [6 100]);
        end
    end
end
```

- [ ] **Step 2: Write the front-end test**

Create `tests/TechBubblesV4nblogKRAv3/parser/tFrontEndTechBubbles.m`:

```matlab
classdef tFrontEndTechBubbles < matlab.unittest.TestCase
    % Front-end gate: TechBubblesV4nblogKRAv3.gmod -> schema-valid IR.
    % The model declares an array policy (bp[2]), a model_init block with
    % var_policy_init/var_aux_init, adaptive(2) inbounds, and SIMU_RESOLVE=1.
    methods (Test)
        function irIsValidAndHasExpectedShape(tc)
            here = fileparts(mfilename('fullpath'));    % tests/.../parser
            modelDir = fileparts(here);
            gmod = fileread(fullfile(modelDir, 'TechBubblesV4nblogKRAv3.gmod'));
            % The gmod assigns params after later declaration keywords, which
            % warns by design; the warnings are not what this test is about.
            ws = warning('off', 'all');
            restore = onCleanup(@() warning(ws));
            ir = gdsge.parser.parseFrontEnd(gmod, 'TechBubblesV4nblogKRAv3', modelDir);

            r = gdsge.ir.validate(ir);
            tc.verifyTrue(r.pass, strjoin(r.errors, ' | '));

            tc.verifyEqual(ir.shocks.count, 2);
            tc.verifyEqual(ir.states.names, {'logK','b'});
            % 7 unknowns: q qn IK c1 logKp + bp[2]
            tc.verifyEqual(sum(cellfun(@(p) p.length, ir.variables.policy)), 7);
            tc.verifyFalse(isempty(ir.modelInit), 'model_init block was dropped');
            tc.verifyEqual(ir.options.simuResolve, 1);
            tc.verifyEqual(ir.options.simuInterp, 0);

            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));
        end
    end
end
```

- [ ] **Step 3: Run it to verify it passes**

Run: `matlab -batch "addpath('src'); r = runtests('tests/TechBubblesV4nblogKRAv3/parser'); disp(table(r)); exit(any([r.Failed]))"`

Expected: PASS. This is a characterization test — the front end already parses this model correctly, and every asserted field name was verified against the real IR while writing this plan (`shocks.count`, `states.names = {'logK','b'}`, `options.simuResolve`/`simuInterp`, `modelInit` a non-empty struct, policy lengths summing to 7 over `q qn IK c1 logKp bp`). If one fails, print the real IR and correct the assertion — do **not** weaken it to a tautology:

```bash
matlab -batch "addpath('src'); \
  d = 'tests/TechBubblesV4nblogKRAv3'; \
  gmod = fileread(fullfile(d,'TechBubblesV4nblogKRAv3.gmod')); warning('off','all'); \
  ir = gdsge.parser.parseFrontEnd(gmod,'TechBubblesV4nblogKRAv3',d); \
  disp(ir.options); disp(ir.states); disp(fieldnames(ir.shocks)');"
```

- [ ] **Step 4: Write the IR snapshot**

Run:

```bash
matlab -batch "addpath('src'); \
  d = 'tests/TechBubblesV4nblogKRAv3'; \
  gmod = fileread(fullfile(d,'TechBubblesV4nblogKRAv3.gmod')); \
  warning('off','all'); \
  ir = gdsge.parser.parseFrontEnd(gmod,'TechBubblesV4nblogKRAv3',d); \
  fid = fopen(fullfile(d,'ir','TechBubblesV4nblogKRAv3.gdsge.json'),'w'); \
  fwrite(fid, gdsge.ir.encode(ir)); fclose(fid); disp('IR snapshot written');"
```

- [ ] **Step 5: Confirm the auto-discovered round-trip test picks it up**

Run: `matlab -batch "addpath('src'); r = runtests('tests/ir'); disp(table(r)); exit(any([r.Failed]))"`

Expected: PASS. `tests/ir/tRoundtrip.m:101` globs `**/ir/*.gdsge.json`, so the new snapshot is covered with no edit to that file.

- [ ] **Step 6: Write the end-to-end gate**

Create `tests/TechBubblesV4nblogKRAv3/codegen/tEndToEndTechBubbles.m`:

```matlab
classdef tEndToEndTechBubbles < matlab.unittest.TestCase
    % GATE (TechBubblesV4nblogKRAv3, autodiff): public API end-to-end vs golden.
    %
    % This is the corpus's ONLY >200-retry model: its simulate needs up to ~5300
    % randomized restarts per period. It exists to lock in MaxMinorIter = inf -
    % with a finite cap the simulate silently truncates and q reaches 987
    % against a true maximum of 2.5. Do not "fix" a failure here by capping.
    %
    % Multi-root grid points make the randomized restart order-sensitive (the
    % safe_assets fragility), so rng is pinned before iter AND before simulate,
    % and comparisons are tolerance-based rather than bit-exact.
    % Slow: MEX compile + 355 iterations + a retry-heavy simulate (~30s total).
    properties (Constant)
        RelTol = 1e-4; AbsTol = 1e-4;
    end
    methods (Test, TestTags = {'Slow'})
        function publicApiMatchesGolden(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'TechBubblesV4nblogKRAv3.gmod'), work);

            ir = gdsge_codegen('TechBubblesV4nblogKRAv3');

            mexFile = fullfile(work, ['mex_TechBubblesV4nblogKRAv3.' mexext]);
            artifacts = {'iter_TechBubblesV4nblogKRAv3.m', ...
                'simulate_TechBubblesV4nblogKRAv3.m', ...
                'mex_TechBubblesV4nblogKRAv3.cpp', ...
                'compile_TechBubblesV4nblogKRAv3.m', ...
                'TechBubblesV4nblogKRAv3.gdsge.json', ...
                'mex_TechBubblesV4nblogKRAv3.cache'};
            for i = 1:numel(artifacts)
                tc.assertTrue(exist(fullfile(work, artifacts{i}), 'file') == 2, artifacts{i});
            end
            tc.assertTrue(exist(mexFile, 'file') == 3, 'MEX did not compile');

            % the uncapped default is the whole point of this gate
            simuTxt = fileread(fullfile(work, 'simulate_TechBubblesV4nblogKRAv3.m'));
            tc.verifyTrue(contains(simuTxt, 'MaxMinorIter = inf;'), ...
                'simulate must not cap resolve retries: this model needs ~5300');
            tc.verifyTrue(contains(simuTxt, 'reportUnconvergedSamples'), ...
                'simulate must use the sample-shaped unconverged reporter');

            % the written IR is the contract: it must round-trip
            decoded = gdsge.ir.decode(fileread( ...
                fullfile(work, 'TechBubblesV4nblogKRAv3.gdsge.json')));
            tc.verifyTrue(gdsge.ir.isequalIR(decoded, ir), 'IR JSON does not round-trip');

            rng(0823);                     % multi-root restart determinism
            IterRslt = iter_TechBubblesV4nblogKRAv3(struct('SaveFreq', inf, 'NoSave', 1));

            golden = load(fullfile(modelDir, 'golden', 'IterRslt.mat'));
            G = golden.IterRslt;
            tc.verifyLessThan(IterRslt.Metric, 1e-6);
            tc.verifyEqual(IterRslt.Iter, G.Iter, 'iteration count diverged from golden');
            r = gdsgetest.compareNumericClose(IterRslt.var_policy, G.var_policy, ...
                tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_aux, G.var_aux, ...
                tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));

            rng(0823);
            SimuRslt = simulate_TechBubblesV4nblogKRAv3(IterRslt, ...
                struct('num_periods', 100, 'num_samples', 6));
            gs = load(fullfile(modelDir, 'golden', 'SimuRslt.mat')); GS = gs.SimuRslt;
            % same vintage + same RNG (gen_discrete_markov_rn) -> shock path matches
            tc.verifyEqual(SimuRslt.shock, GS.shock, 'shock path diverged');
            % every simulated field: a capped run diverges by orders of magnitude
            fields = {'logK','b','q','qd','ER','r','qn','IK'};
            for i = 1:numel(fields)
                r = gdsgetest.compareNumericClose(SimuRslt.(fields{i}), ...
                    GS.(fields{i}), 1e-3, 1e-3);
                tc.verifyTrue(r.pass, sprintf('%s: %s', fields{i}, ...
                    strjoin(r.failures, newline)));
            end
        end
    end
end
```

- [ ] **Step 7: Run the gate**

Run: `matlab -batch "addpath('src'); r = runtests('tests/TechBubblesV4nblogKRAv3','IncludingSubfolders',true); disp(table(r)); exit(any([r.Failed]))"`

Expected: PASS, ~60 s (compile + iter + simulate).

If the simulate comparison fails while `Iter` matches, check the residual scale before touching tolerances: a **capped** run is wrong by orders of magnitude (`q` ~987 vs ~2.5), whereas restart-ordering noise is ~1e-10. Only the latter justifies a tolerance discussion.

- [ ] **Step 8: Commit**

```bash
git add tests/TechBubblesV4nblogKRAv3/
git commit -m "test: TechBubblesV4nblogKRAv3 end-to-end gate

The corpus's only >200-retry model, locking in MaxMinorIter = inf. Front-end
IR check, auto-discovered IR round-trip snapshot, golden-shape test, and a
full public-API gate asserting every simulated field against the old
toolbox.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_011HSDe48haPWHSrWuTrGanS"
```

---

### Task 8: Full suite, docs, and changelog

**Files:**
- Modify: `PROGRESS.md` (new entry at the top of the Changelog section)
- Modify: `Changelog.md` (note: the dev changelog is `Changelog.md`, **not** `CHANGELOG.md`)

**Interfaces:**
- Consumes: all previous tasks.
- Produces: nothing.

- [ ] **Step 1: Run the complete suite**

Run: `matlab -batch "cd('tests'); run_tests"`

Expected: exit 0. This takes several minutes (the corpus gates dominate).

- [ ] **Step 2: Confirm the result from the authoritative source**

```bash
grep -o 'tests="[0-9]*"\|failures="[0-9]*"\|errors="[0-9]*"' tests/results/junit.xml | head
```

Expected: `failures="0"` and `errors="0"`, with a test count of roughly 615 (600 before this work, plus ~15 new). **Do not read `tests/results/results.tap`** — it accumulates documents across runs and shows stale `not ok` lines.

- [ ] **Step 3: Verify no unintended golden churn**

```bash
git status --short
git diff --stat origin/master...HEAD -- tests/ | tail -20
```

Expected: the only modified pre-existing goldens are `tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt` and `simulate_HL1996_golden.txt`. Every other change under `tests/` is a new file. If any other model's golden or IR snapshot changed, stop and investigate.

- [ ] **Step 4: Update `PROGRESS.md`**

Add at the top of the `## Changelog` section (immediately above the `2026-07-03: **Preamble MATLAB passthrough...**` entry):

```markdown
- 2026-07-31: **Uncapped resolve retries + simulate diagnostics + install self-check.**
  `TechBubblesV4nblogKRAv3` ran on the old toolbox but its simulate produced garbage on the
  new one: `emitSetup` emitted `MaxMinorIter = 200` where the old toolbox uses `inf`, and the
  model needs up to ~5300 randomized restarts per period. The cap did not fail loudly — the
  simulate path wrote the unconverged solution into `SimuRslt` and fed it forward, so `q`
  reached 987 against a true maximum of 2.5. Restoring `inf` reproduces the old toolbox to
  ≤9.1e-10 on every field with the shock path bit-exact, at parity cost (44.7s old vs 52.5s
  new for 20 periods × 100 samples). Three defects the cap had been masking are also fixed:
  `reportUnconverged` threw `MATLAB:badsubscript` on the simulate path for any multi-state
  model (new `reportUnconvergedSamples`); the diagnostics heartbeat used `mod(minorIter,
  diagnoseAt)`, which is unbounded uncapped and skippable by the in-MEX batch step (now
  geometric backoff); and simulate now aggregates retry stats at the `SimuPrintFreq`
  checkpoint instead of a per-period block. Separately, `checkInstall` (run once per session
  from `ensurePath`) turns a truncated source file — reported from a macOS install as
  `MATLAB:scriptNotAFunction` naming an internal file — into one actionable error. New
  corpus gate `tests/TechBubblesV4nblogKRAv3` is the only >200-retry model in the suite.
  No IR schema change. Spec:
  `docs/superpowers/specs/2026-07-31-techbubbles-retry-cap-design.md`.
```

- [ ] **Step 5: Update `Changelog.md`**

The dev changelog is `Changelog.md` (capital C, lowercase rest) — **not** `CHANGELOG.md`. It follows Keep a Changelog with `## [version] - date` sections. This work is unreleased, so add a new section directly under the intro paragraph and above `## [0.2.1] - 2026-07-03`:

```markdown
## [Unreleased]

Driven by `TechBubblesV4nblogKRAv3`, an asset-pricing model with growth shocks whose
simulation needs thousands of randomized solver restarts per period.

### Fixed

- Resolve retries are unbounded again by default (`MaxMinorIter = inf`, as in the
  pre-refactor toolbox). The previous default of 200 did not fail loudly: `simulate`
  wrote the unconverged solution into `SimuRslt` and fed it forward as the next
  period's state, so a model needing more restarts produced plausible-looking but
  wholly wrong paths. Set `MaxMinorIter` yourself to restore a bound — doing so also
  arms the existing hard error on the iteration path.
- `simulate` no longer crashes with `MATLAB:badsubscript` when it reports unconverged
  samples. The report was grid-shaped and mis-indexed for any model with more than one
  state; it now names the sample and its own shock and state values.

### Changed

- `simulate` reports resolve activity once per `SimuPrintFreq` checkpoint
  (`resolve: N/M periods retried, ... worst period P`) instead of a full diagnostics
  block every few retries. Set `DiagnoseMinorIter` to re-enable the per-solve block.
- The diagnostics heartbeat backs off geometrically (at `DiagnoseMinorIter`, then each
  doubling) rather than firing at a fixed interval, so a long solve reports about
  `log2(N)` times instead of `N/interval`.

### Added

- An install self-check runs once per session and reports a truncated or empty toolbox
  source file as `gdsge:install:corrupt`, naming every damaged file and the remedy,
  instead of MATLAB's `MATLAB:scriptNotAFunction` pointing at an internal file.
```

If a release number has already been chosen for this work, use it in place of `[Unreleased]` — **never stamp a version number without asking.**

- [ ] **Step 6: Commit**

```bash
git add PROGRESS.md Changelog.md
git commit -m "docs: log the uncapped-retry fix, simulate diagnostics, install self-check

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_011HSDe48haPWHSrWuTrGanS"
```

---

## Verification summary

| Spec section | Task |
|---|---|
| §3 Fix A — `MaxMinorIter = inf` | 2 |
| §4 Fix B — `reportUnconvergedSamples` | 4 |
| §5.1 backoff | 3 |
| §5.2 simulate `DiagnoseMinorIter = 0` | 2 |
| §5.3 simulate retry aggregate | 5 |
| §6 Fix D — `checkInstall` | 1 |
| §7 corpus gate | 6, 7 |
| §8 verification | every task's test steps, plus 8 |

§9 (in-MEX randomize for simulate; the four parser warnings on this gmod) is explicitly out of scope and has no task.
