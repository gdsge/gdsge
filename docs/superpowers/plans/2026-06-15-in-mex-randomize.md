# In-MEX randomized restart (`UseMexRandomize`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the random-restart trial loop (and the outer-iter-1 initial guess) out of MATLAB into the C++ MEX, driven by a deterministic, fixed-seedable per-grid-point RNG, batching up to `MexRandomizeBatch` (default 100) minor-iterations per MEX call.

**Architecture:** A counter-based RNG helper in `mex.tpl.cpp` makes every random draw a pure function of `(seed, salt, point, trial, component)` — thread-independent and batch-resumable. `task.tpl.cpp` gains a batched restart loop (reusing the existing neighbor-resolve sweep, refactored into a lambda). `gdsge.runtime.solveProblems` drives it via the existing MEX caller-workspace contract, keeping diagnostics / cap / exit decisions in MATLAB at batch boundaries. Codegen adds three public options; the iter initial guess becomes a deterministic midpoint when the feature is on. Cartesian iter path only; ASG and simulate/init keep today's MATLAB path. Default-on; iter goldens re-baselined after a tolerance-equivalence check.

**Tech Stack:** MATLAB R2025b (`matlab.unittest`), C++ MEX (MSVC, OpenMP), the project's `gdsge.codegen`/`gdsge.runtime` packages.

**Spec:** `docs/superpowers/specs/2026-06-15-in-mex-randomize-design.md`
**Branch:** `inmex-randomize`

**How to run tests:** `matlab -batch "cd('tests'); run_tests"` (exit 0 = all pass; report at `tests/results/junit.xml`). A single class: `matlab -batch "cd('tests'); runtests('tInMexRandomizeDriver')"`. Run MATLAB **one process at a time** (OpenMP saturates all cores).

---

## File structure

| File | Responsibility | Action |
|------|----------------|--------|
| `templates/cxx/mex.tpl.cpp` | Top-level MEX file (emitted once). Home of the counter-based RNG helper. | Modify |
| `templates/cxx/task.tpl.cpp` | Per-task solve body. Neighbor sweep → lambda; new batched randomize loop. | Modify |
| `src/+gdsge/+runtime/solveProblems.m` | MATLAB resolve driver. New `useMexRandomize` branch. | Modify |
| `src/+gdsge/+codegen/+mat/optionsWhitelist.m` | Frozen option surface. | Modify |
| `src/+gdsge/+codegen/+mat/emitSetup.m` | Generated-file option defaults. | Modify |
| `src/+gdsge/+codegen/+mat/emitIter.m` | Iter cfg wiring + initial-guess branch. | Modify |
| `tests/runtime/tInMexRandomizeDriver.m` | Fast unit test of the MATLAB driver (fake MEX). | Create |
| `tests/codegen/tEmitRandomizeOptions.m` | Whitelist + emitSetup defaults for the 3 new options. | Create |
| `tests/codegen/tEmitIterRandomize.m` | emitIter cfg wiring + midpoint init branch. | Create |
| `tests/HeatonLucas1996/codegen/tInMexRandomizeDeterminism.m` | Thread-independence gate (Slow). | Create |
| `tests/HeatonLucas1996/codegen/tInMexRandomizeAB_HL1996.m` | `UseMexRandomize=0` == old golden; `=1` within tolerance (Slow). | Create |
| `tests/Barro_et_al_2017/codegen/tInMexRandomizeSafeAssets.m` | Decisive restart-heavy A/B + tolerance gate (Slow). | Create |
| `tests/HeatonLucas1996/codegen/golden/*.txt` | Regenerated iter + cpp snapshots. | Regenerate |
| `tests/*/golden/IterRslt.mat` | Re-baselined iter goldens where output changed. | Regenerate (Task 9) |
| `PROGRESS.md` | Changelog entry. | Modify |

---

## Task 1: C++ counter-based RNG helper

Add a deterministic, thread-independent RNG to the once-emitted `mex.tpl.cpp`. (`task.tpl.cpp` is instantiated **twice** — `task_inf_horizon` and `task_init` — so a free function there would be an ODR redefinition; it must live in `mex.tpl.cpp`.)

**Files:**
- Modify: `templates/cxx/mex.tpl.cpp` (add `#include <cstdint>` and the helper before the task functions)

- [ ] **Step 1: Add the include**

In `templates/cxx/mex.tpl.cpp`, after the existing `#include "cmath"` (line 28), add:

```cpp
#include "cmath"
#include <cstdint>
```

- [ ] **Step 2: Add the RNG helper before the forward declarations**

In `templates/cxx/mex.tpl.cpp`, immediately after `using adept::Stack;` (line 47) and before the `#ifdef HAS_INIT` forward declaration block, insert:

```cpp
// ---------------------------------------------------------------------------
// Counter-based RNG for in-MEX randomized restart (UseMexRandomize).
// Every draw is a pure function of (seed, salt, point, trial, component), so
// results are independent of OpenMP thread count/scheduling and resumable
// across batches. splitmix64 mixing -> uniform double in [0,1).
// ---------------------------------------------------------------------------
static inline uint64_t gdsge_splitmix64(uint64_t x)
{
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

static inline double gdsge_rng_uniform(uint64_t seed, uint64_t salt, uint64_t point,
                                       uint64_t trial, uint64_t comp)
{
  uint64_t s = seed;
  s = gdsge_splitmix64(s ^ (salt  * 0x9E3779B97F4A7C15ULL + 0x1ULL));
  s = gdsge_splitmix64(s ^ (point * 0xD1B54A32D192ED03ULL + 0x2ULL));
  s = gdsge_splitmix64(s ^ (trial * 0x9E3779B97F4A7C15ULL + 0x5ULL));
  s = gdsge_splitmix64(s ^ (comp  * 0xD1B54A32D192ED03ULL + 0x7ULL));
  // 53-bit mantissa -> [0,1)
  return (double)(s >> 11) * (1.0 / 9007199254740992.0);
}
```

- [ ] **Step 3: Verify it compiles as part of a model build (deferred to Task 3)**

No standalone C++ test harness exists; this helper is exercised by the Slow determinism gate (Task 7) and validated to compile when the HL1996 MEX is built there. For now, confirm the helper text is present:

Run: `matlab -batch "t=fileread(fullfile('templates','cxx','mex.tpl.cpp')); assert(contains(t,'gdsge_rng_uniform')); assert(contains(t,'#include <cstdint>')); disp('ok')"`
Expected: prints `ok`.

- [ ] **Step 4: Commit**

```bash
git add templates/cxx/mex.tpl.cpp
git commit -m "feat(codegen): add counter-based RNG helper to mex template (UseMexRandomize)"
```

---

## Task 2: C++ batched randomize loop in `task.tpl.cpp`

Refactor the existing in-MEX neighbor-resolve block into a reusable lambda, then add the batched restart loop guarded by `GDSGE_RANDOMIZE_BATCH`.

**Files:**
- Modify: `templates/cxx/task.tpl.cpp` (lines 56–171)

- [ ] **Step 1: Read the new randomize controls from the caller workspace**

In `templates/cxx/task.tpl.cpp`, immediately after the existing strides block (the closing `}` of the `GDSGE_STRIDES_MX` block, currently line 67), insert:

```cpp
  // Optional in-MEX randomized restart (UseMexRandomize). Absent/<=0 => today's
  // behavior (no in-MEX randomize; MATLAB drives restarts). Read from caller.
  int      GDSGE_RAND_BATCH = 0;
  uint64_t GDSGE_RAND_SEED = 0, GDSGE_RAND_SALT = 0, GDSGE_RAND_TRIAL_OFFSET = 0;
  {
    const mxArray* p = mexGetVariablePtr("caller", "GDSGE_RANDOMIZE_BATCH");
    if (p != 0 && mxIsNumeric(p) && mxGetNumberOfElements(p) > 0)
      GDSGE_RAND_BATCH = (int) mxGetScalar(p);
    if (GDSGE_RAND_BATCH > 0) {
      const mxArray* ps = mexGetVariablePtr("caller", "GDSGE_RANDOM_SEED");
      if (ps != 0 && mxGetNumberOfElements(ps) > 0) GDSGE_RAND_SEED = (uint64_t) mxGetScalar(ps);
      const mxArray* pa = mexGetVariablePtr("caller", "GDSGE_RANDOMIZE_SALT");
      if (pa != 0 && mxGetNumberOfElements(pa) > 0) GDSGE_RAND_SALT = (uint64_t) mxGetScalar(pa);
      const mxArray* pt = mexGetVariablePtr("caller", "GDSGE_RANDOMIZE_TRIAL_OFFSET");
      if (pt != 0 && mxGetNumberOfElements(pt) > 0) GDSGE_RAND_TRIAL_OFFSET = (uint64_t) mxGetScalar(pt);
    }
  }
```

- [ ] **Step 2: Wrap the existing neighbor-resolve block in a lambda**

In `templates/cxx/task.tpl.cpp`, replace the current `if (GDSGE_NUM_DIMS > 0)` header (line 124) so the body becomes a reusable lambda. Change:

```cpp
  // In-MEX nearby-warmup resolve (deterministic snapshot semantics; bit-exact
  // with solveProblems.m's MATLAB sweep). Runs only when strides were supplied.
  if (GDSGE_NUM_DIMS > 0)
  {
    vector<char> GDSGE_SOLVED0(GDSGE_NPROB + 1);
```

to:

```cpp
  // In-MEX nearby-warmup resolve (deterministic snapshot semantics; bit-exact
  // with solveProblems.m's MATLAB sweep). Reusable: called by the standalone
  // resolve sweep and by each batched randomize trial below.
  auto GDSGE_neighbor_fixpoint = [&]()
  {
    if (GDSGE_NUM_DIMS <= 0) return;
    vector<char> GDSGE_SOLVED0(GDSGE_NPROB + 1);
```

Then change the **closing** brace of that block (currently line 171, `  }`) to close the lambda with a semicolon:

```cpp
      GDSGE_AFTER = GDSGE_count_unconv();
    }
  };
```

- [ ] **Step 3: Add the batched randomize loop + standalone-sweep dispatch**

In `templates/cxx/task.tpl.cpp`, **insert the following immediately before the function's final closing brace `}`** (the original last line of the file — do NOT add another brace; the existing `}` stays and closes the function). The block goes right after the lambda's closing `};` from Step 2:

```cpp

  // Dispatch: batched in-MEX randomize (UseMexRandomize) takes precedence; else
  // the standalone neighbor sweep (UseMexResolve); else nothing (today's Pass-0-only).
  if (GDSGE_RAND_BATCH > 0)
  {
    GDSGE_neighbor_fixpoint();   // propagate warm-started / Pass-0 solutions first
    for (int GDSGE_TRIAL = 0; GDSGE_TRIAL < GDSGE_RAND_BATCH; GDSGE_TRIAL++)
    {
      int GDSGE_NLEFT = 0;
      for (int i = 1; i <= GDSGE_NPROB; i++) if (!(GDSGE_F(i) <= TolSol)) GDSGE_NLEFT++;
      if (GDSGE_NLEFT == 0) break;

      uint64_t GDSGE_T = GDSGE_RAND_TRIAL_OFFSET + (uint64_t) GDSGE_TRIAL;
#pragma omp parallel for schedule(dynamic)
      for (int GDSGE_I = 1; GDSGE_I <= GDSGE_NPROB; GDSGE_I++)
      {
        if (GDSGE_F(GDSGE_I) <= TolSol) continue;   // keep converged points
        for (int GDSGE_C = 1; GDSGE_C <= NUM_EQUATIONS; GDSGE_C++)
        {
          double u = gdsge_rng_uniform(GDSGE_RAND_SEED, GDSGE_RAND_SALT,
                                       (uint64_t)(GDSGE_I - 1), GDSGE_T, (uint64_t)(GDSGE_C - 1));
          GDSGE_SOL(GDSGE_C, GDSGE_I) =
              GDSGE_LB(GDSGE_C, GDSGE_I) + u * (GDSGE_UB(GDSGE_C, GDSGE_I) - GDSGE_LB(GDSGE_C, GDSGE_I));
        }
        GDSGE_solve_one(GDSGE_I);
      }
      GDSGE_neighbor_fixpoint();
    }
  }
  else if (GDSGE_NUM_DIMS > 0)
  {
    GDSGE_neighbor_fixpoint();   // standalone UseMexResolve sweep (unchanged behavior)
  }
```

The existing in-MEX resolve block (now the lambda) no longer runs itself — the dispatch above calls it. Make sure the old `if (GDSGE_NUM_DIMS > 0) { ... }` invocation is fully converted to the lambda definition in Step 2 (it must not also execute inline).

- [ ] **Step 4: Verify the template parses (sanity)**

Run: `matlab -batch "t=fileread(fullfile('templates','cxx','task.tpl.cpp')); assert(count(t,'GDSGE_neighbor_fixpoint')>=3); assert(count(t,'gdsge_rng_uniform')==1); disp('ok')"`
Expected: prints `ok` (lambda defined + called twice in dispatch; RNG used once).

- [ ] **Step 5: Commit**

```bash
git add templates/cxx/task.tpl.cpp
git commit -m "feat(codegen): batched in-MEX randomize loop + neighbor-fixpoint lambda"
```

---

## Task 3: MATLAB driver — `solveProblems` `useMexRandomize` branch

Add the batched-call loop. The `useMexRandomize=false` path stays byte-for-byte as today.

**Files:**
- Modify: `src/+gdsge/+runtime/solveProblems.m`
- Test: `tests/runtime/tInMexRandomizeDriver.m` (Create)

- [ ] **Step 1: Write the failing driver unit test**

Create `tests/runtime/tInMexRandomizeDriver.m`:

```matlab
classdef tInMexRandomizeDriver < matlab.unittest.TestCase
    % Fast unit test of solveProblems' MATLAB driver under cfg.useMexRandomize.
    % The fake MEX is a NESTED function passed directly so solveProblems is its
    % caller; evalin('caller',...) reads solveProblems' workspace, exactly how
    % the real MEX reads GDSGE_RANDOMIZE_* via mexGetVariable('caller',...).
    methods (Static)
        function [calls, fOut, diag] = runWithFake(maxMinorIter, batch, nNeeded)
            calls = {};
            n = 0;
            cfg = struct('tolFun',1e-8,'tolSol',1e-8,'solMaxIter',200, ...
                'numThreads',1,'debugEvalOnly',0,'useBroyden',0, ...
                'finiteDiffDelta',1e-6,'useBroydenNow',0,'taskName',1, ...
                'splineVec',[],'ppNames',{{}},'ppCell',{{}}, ...
                'maxMinorIter',maxMinorIter,'probSize',[2 3], ...
                'useNearestNeighbor',true,'verboseRetry',false,'adaptInSol',[], ...
                'useMexResolve',true,'useMexRandomize',true, ...
                'randomizeBatch',batch,'randomSeed',7,'randomizeSalt',3);
            np = 6;                                  % prod([2 3])
            s = zeros(1,np); l = zeros(1,np); u = ones(1,np);
            d = zeros(1,np); f = 1e20*ones(1,np); a = zeros(1,np); e = zeros(1,np);
            [~, fOut, ~, ~, ~, diag] = gdsge.runtime.solveProblems( ...
                @fake, s,l,u,d,f,a,e, cfg);

            function [sol, f, aux, eqVal, optInfo] = fake(sol, lb, ub, data, skip, f, aux, eqVal) %#ok<INUSD>
                n = n + 1;
                rec = struct('n',n,'skip',skip(:)');
                if n == 1
                    rec.kind = 'pass0';
                    f(:) = 1e20; f(1) = 0;             % Pass 0 converges only problem 1
                else
                    rec.kind = 'batch';
                    rec.batch  = evalin('caller','GDSGE_RANDOMIZE_BATCH');
                    rec.seed   = evalin('caller','GDSGE_RANDOM_SEED');
                    rec.salt   = evalin('caller','GDSGE_RANDOMIZE_SALT');
                    rec.offset = evalin('caller','GDSGE_RANDOMIZE_TRIAL_OFFSET');
                    rec.strides = evalin('caller','GDSGE_PROBLEM_STRIDES');
                    % "converge nNeeded more points per batch call"
                    nowConv = nnz(f <= 1e-8);
                    target = min(np, nowConv + nNeeded);
                    f(1:target) = 0; f(target+1:end) = 1e20;
                end
                calls{end+1} = rec;
                optInfo = zeros(1, numel(f));
            end
        end
    end
    methods (Test)
        function batchCallsCarryContractAndConverge(tc)
            % batch=2, converge 2 more each batch call: needs problems 2..6 (5)
            % after Pass-0 gives 1 => 3 batch calls (3,5,6 -> but capped at 6).
            [calls, fOut, diag] = tInMexRandomizeDriver.runWithFake(inf, 2, 2);
            tc.verifyEqual(fOut, zeros(1,6), 'all problems should end converged');
            batchCalls = calls(strcmp(cellfun(@(c)c.kind,calls,'uni',0), 'batch'));
            tc.verifyGreaterThanOrEqual(numel(batchCalls), 1);
            b1 = batchCalls{1};
            tc.verifyEqual(b1.batch, 2, 'batch size passed to MEX');
            tc.verifyEqual(b1.seed, 7, 'seed passed to MEX');
            tc.verifyEqual(b1.salt, 3, 'salt passed to MEX');
            tc.verifyEqual(b1.offset, 0, 'first batch trialOffset = 0');
            tc.verifyEqual(b1.strides(:)', [1 2], 'strides = cumprod([1 probSize(1:end-1)])');
            tc.verifyEqual(b1.skip, ones(1,6), 'batch calls pass skip(:)=1');
        end
        function capStopsAndReportsDiag(tc)
            % nNeeded=0 => never fully converges; cap=10, batch=4 => offsets 0,4,8.
            [calls, fOut, diag] = tInMexRandomizeDriver.runWithFake(10, 4, 0);
            batchCalls = calls(strcmp(cellfun(@(c)c.kind,calls,'uni',0), 'batch'));
            offsets = cellfun(@(c) c.offset, batchCalls);
            tc.verifyEqual(offsets, [0 4 8], 'trialOffset advances by batch each call');
            batches = cellfun(@(c) c.batch, batchCalls);
            tc.verifyEqual(batches, [4 4 2], 'batch is cap-aware: min(batch, cap-minorIter)');
            tc.verifyGreaterThan(max(fOut), 1e-8, 'still unconverged at the cap');
            tc.verifyEqual(diag.minorIters, 10, 'minorIter reaches the cap');
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); r=runtests('tInMexRandomizeDriver'); disp(table(r))"`
Expected: FAIL — `solveProblems` ignores `cfg.useMexRandomize` (still runs the MATLAB restart path), so batch contract variables aren't set / offsets wrong.

- [ ] **Step 3: Implement the `useMexRandomize` branch in `solveProblems`**

In `src/+gdsge/+runtime/solveProblems.m`, after the Pass-0 call (currently line 51) and the two lines initializing `minorIter`/`numNeedResolvedAfter` (lines 54–55), wrap the existing `while` loop so the randomize branch is taken first. Replace:

```matlab
minorIter = 0;
numNeedResolvedAfter = inf;
while ((max(isnan(f)) || max(f(:)) > cfg.tolSol) && minorIter < cfg.maxMinorIter)
```

with:

```matlab
minorIter = 0;
numNeedResolvedAfter = inf;

useMexRandomize = isfield(cfg,'useMexRandomize') && cfg.useMexRandomize ...
    && isempty(cfg.adaptInSol) && useMexResolve;
if useMexRandomize
    % In-MEX randomized restart: each MEX call runs up to `batch` full
    % minor-iterations (neighbor sweep + random restart) internally. MATLAB
    % keeps diagnostics, the cap, and the exit decision at batch boundaries.
    GDSGE_RANDOM_SEED = cfg.randomSeed;            %#ok<NASGU> caller-workspace contract
    GDSGE_RANDOMIZE_SALT = cfg.randomizeSalt;      %#ok<NASGU>
    GDSGE_RANDOMIZE_TRIAL_OFFSET = 0;              %#ok<NASGU>
    GDSGE_RANDOMIZE_BATCH = 0;                     %#ok<NASGU>
    trialOffset = 0;
    while ((max(isnan(f)) || max(f(:)) > cfg.tolSol) && minorIter < cfg.maxMinorIter)
        batch = min(cfg.randomizeBatch, cfg.maxMinorIter - minorIter);
        GDSGE_RANDOMIZE_BATCH = batch;                       %#ok<NASGU>
        GDSGE_RANDOMIZE_TRIAL_OFFSET = trialOffset;          %#ok<NASGU>
        GDSGE_PROBLEM_STRIDES = GDSGE_RESOLVE_STRIDES;       %#ok<NASGU> enable neighbor sweep
        skip(:) = 1;                                         % Pass-0 no-op; restarts gate on f
        [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);
        GDSGE_RANDOMIZE_BATCH = 0;                  %#ok<NASGU> off for any later call
        GDSGE_PROBLEM_STRIDES = [];                 %#ok<NASGU>
        minorIter = minorIter + batch;
        trialOffset = trialOffset + batch;

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
```

Then add a closing `end` for the new `else` block. The existing `while` body (lines 57–140) is unchanged; add `end` right **before** the "Exhaustion" comment block (currently line 142), i.e. after the existing `while`'s `end` (line 140):

```matlab
    end
end   % <-- closes the `else` wrapping the legacy while-loop
```

(There is now: `if useMexRandomize ... else <legacy while ... end> end`.)

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "cd('tests'); r=runtests('tInMexRandomizeDriver'); assert(all([r.Passed])); disp('PASS')"`
Expected: PASS.

- [ ] **Step 5: Run the existing solveProblems tests (no regression)**

Run: `matlab -batch "cd('tests'); r=runtests({'tSolveProblems','tInMexResolveDriver'}); assert(all([r.Passed])); disp('PASS')"`
Expected: PASS — the `useMexRandomize=false` / absent path is unchanged.

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+runtime/solveProblems.m tests/runtime/tInMexRandomizeDriver.m
git commit -m "feat(runtime): in-MEX randomize driver branch in solveProblems"
```

---

## Task 4: Codegen — options whitelist + defaults

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/optionsWhitelist.m`
- Modify: `src/+gdsge/+codegen/+mat/emitSetup.m`
- Test: `tests/codegen/tEmitRandomizeOptions.m` (Create — self-contained, uses `buildHL1996IR`)

- [ ] **Step 1: Write the failing test**

Create `tests/codegen/tEmitRandomizeOptions.m` (self-contained; `buildHL1996IR` is a full IR fixture already in the repo):

```matlab
classdef tEmitRandomizeOptions < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addIrPath(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', 'HeatonLucas1996', 'ir')));
        end
    end
    methods (Test)
        function optionsWhitelisted(tc)
            lit = gdsge.codegen.mat.optionsWhitelist(buildHL1996IR(), {});
            for opt = {'UseMexRandomize','MexRandomizeBatch','MexRandomSeed'}
                tc.verifyTrue(contains(lit, ['''' opt{1} '''']), ...
                    sprintf('%s must be whitelisted', opt{1}));
            end
        end
        function defaultsEmitted(tc)
            txt = gdsge.codegen.mat.emitSetup(buildHL1996IR(), 'iter');
            tc.verifyTrue(contains(txt, 'UseMexRandomize = 1;'));
            tc.verifyTrue(contains(txt, 'MexRandomizeBatch = 100;'));
            tc.verifyTrue(contains(txt, 'MexRandomSeed = 0;'));
        end
    end
end
```

- [ ] **Step 2: Run to verify failure**

Run: `matlab -batch "cd('tests'); r=runtests('tEmitRandomizeOptions'); disp(table(r))"`
Expected: FAIL — options/defaults not yet present.

- [ ] **Step 3: Add the options to the whitelist**

In `src/+gdsge/+codegen/+mat/optionsWhitelist.m`, extend the `base` cell (the last line ends `'UseMexResolve','DiagnoseMinorIter'};`). Change it to:

```matlab
    'CONSTRUCT_OUTPUT','NumThreads','WarmUp','UseMexResolve','DiagnoseMinorIter', ...
    'UseMexRandomize','MexRandomizeBatch','MexRandomSeed'};
```

- [ ] **Step 4: Add the defaults to emitSetup**

In `src/+gdsge/+codegen/+mat/emitSetup.m`, after the `UseMexResolve = 1;` line (line 44), add:

```matlab
w.add('UseMexResolve = 1;');
w.add('UseMexRandomize = 1;');     % restart loop runs in the MEX (C++ RNG); 0 = MATLAB path
w.add('MexRandomizeBatch = 100;'); % minor-iterations per MEX call before returning to MATLAB
w.add('MexRandomSeed = 0;');       % fixed seed; reproducible out of the box
```

- [ ] **Step 5: Run to verify pass**

Run: `matlab -batch "cd('tests'); r=runtests('tEmitRandomizeOptions'); assert(all([r.Passed])); disp('PASS')"`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+codegen/+mat/optionsWhitelist.m src/+gdsge/+codegen/+mat/emitSetup.m tests/codegen/tEmitRandomizeOptions.m
git commit -m "feat(codegen): whitelist + defaults for UseMexRandomize/MexRandomizeBatch/MexRandomSeed"
```

---

## Task 5: Codegen — `emitIter` cfg wiring + initial-guess branch

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitIter.m` (initial guess ~line 94; cfg block ~lines 178–186; per-iter cfg updates ~lines 206–208)
- Test: `tests/HeatonLucas1996/codegen/tSnapshotHL1996.m` exercises this via the golden text; we add a focused contains-test first.

- [ ] **Step 1: Write the failing focused test**

Create `tests/codegen/tEmitIterRandomize.m`:

```matlab
classdef tEmitIterRandomize < matlab.unittest.TestCase
    methods (Test)
        function emitsRandomizeCfgAndInitBranch(tc)
            here = fileparts(mfilename('fullpath'));
            addpath(fullfile(here, '..', 'HeatonLucas1996', 'ir'));
            ir = buildHL1996IR();
            txt = gdsge.codegen.mat.emitIter(ir);
            tc.verifyTrue(contains(txt, 'GDSGE_CFG.useMexRandomize = UseMexRandomize;'));
            tc.verifyTrue(contains(txt, 'GDSGE_CFG.randomizeBatch = MexRandomizeBatch;'));
            tc.verifyTrue(contains(txt, 'GDSGE_CFG.randomSeed = MexRandomSeed;'));
            tc.verifyTrue(contains(txt, 'GDSGE_CFG.randomizeSalt = GDSGE_Iter;'));
            tc.verifyTrue(contains(txt, 'if UseMexRandomize'));
            tc.verifyTrue(contains(txt, '(GDSGE_LB + GDSGE_UB) / 2'));
        end
    end
end
```

- [ ] **Step 2: Run to verify failure**

Run: `matlab -batch "cd('tests'); r=runtests('tEmitIterRandomize'); disp(table(r))"`
Expected: FAIL — those lines aren't emitted yet.

- [ ] **Step 3: Replace the initial-guess line with a runtime branch**

In `src/+gdsge/+codegen/+mat/emitIter.m`, replace (line 94–95):

```matlab
w.add('GDSGE_X0 = rand(size(GDSGE_SOL)) .* (GDSGE_UB-GDSGE_LB) + GDSGE_LB;');
w.add('GDSGE_SOL(:) = GDSGE_X0;');
```

with:

```matlab
w.add('if UseMexRandomize');
w.add('    GDSGE_SOL(:) = (GDSGE_LB + GDSGE_UB) / 2;   %% deterministic midpoint; restarts drawn in the MEX');
w.add('else');
w.add('    GDSGE_X0 = rand(size(GDSGE_SOL)) .* (GDSGE_UB-GDSGE_LB) + GDSGE_LB;');
w.add('    GDSGE_SOL(:) = GDSGE_X0;');
w.add('end');
```

- [ ] **Step 4: Add randomize cfg fields to the once-built cfg block**

In `src/+gdsge/+codegen/+mat/emitIter.m`, after the `GDSGE_CFG.useMexResolve = UseMexResolve;` line (line 181), add:

```matlab
w.add('GDSGE_CFG.useMexResolve = UseMexResolve;');
w.add('GDSGE_CFG.useMexRandomize = UseMexRandomize;');
w.add('GDSGE_CFG.randomizeBatch = MexRandomizeBatch;');
w.add('GDSGE_CFG.randomSeed = MexRandomSeed;');
```

- [ ] **Step 5: Set the per-iteration salt inside the loop**

In `src/+gdsge/+codegen/+mat/emitIter.m`, in the per-iteration cfg-update block (after `GDSGE_CFG.useBroydenNow = ...;`, line 206), add:

```matlab
w.add('GDSGE_CFG.useBroydenNow = double(GDSGE_Iter>1)*GDSGE_USE_BROYDEN;');
w.add('GDSGE_CFG.randomizeSalt = GDSGE_Iter;');
```

- [ ] **Step 6: Run to verify pass**

Run: `matlab -batch "cd('tests'); r=runtests('tEmitIterRandomize'); assert(all([r.Passed])); disp('PASS')"`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/+gdsge/+codegen/+mat/emitIter.m tests/codegen/tEmitIterRandomize.m
git commit -m "feat(codegen): emitIter threads randomize cfg + midpoint init under UseMexRandomize"
```

---

## Task 6: Regenerate HL1996 text snapshots (mechanical)

The generated iter `.m` and the `.cpp` changed, so the committed text snapshots must be regenerated. **Review the diff** — the snapshot IS the spec of the output.

**Files:**
- Regenerate: `tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt`, `mex_HL1996_golden_cpp.txt`

- [ ] **Step 1: Regenerate snapshots**

Run: `matlab -batch "cd('tests/HeatonLucas1996/codegen'); regen_snapshots"`
Expected: prints `Snapshots written to ...`.

- [ ] **Step 2: Review the diff**

Run: `git --no-pager diff -- tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt tests/HeatonLucas1996/codegen/golden/mex_HL1996_golden_cpp.txt`
Expected: only the new `if UseMexRandomize` init branch + the `GDSGE_CFG.*randomize*` lines in iter; the RNG helper + batched loop + `GDSGE_neighbor_fixpoint` lambda + `GDSGE_RANDOMIZE_*` reads in cpp. No unrelated churn.

- [ ] **Step 3: Run the snapshot tests**

Run: `matlab -batch "cd('tests'); r=runtests({'tSnapshotHL1996','tSnapshotCxxHL1996'}); assert(all([r.Passed])); disp('PASS')"`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt tests/HeatonLucas1996/codegen/golden/mex_HL1996_golden_cpp.txt
git commit -m "test(codegen): regenerate HL1996 iter+cpp snapshots for UseMexRandomize"
```

---

## Task 7: Determinism gate (thread-independence, Slow)

Prove the C++ RNG is independent of `NumThreads` — the core "fix seed" guarantee. Compiles the HL1996 MEX (~minutes).

**Files:**
- Create: `tests/HeatonLucas1996/codegen/tInMexRandomizeDeterminism.m`

- [ ] **Step 1: Write the determinism test**

Create `tests/HeatonLucas1996/codegen/tInMexRandomizeDeterminism.m`:

```matlab
classdef tInMexRandomizeDeterminism < matlab.unittest.TestCase
    % The in-MEX counter RNG must be independent of NumThreads: same seed/salt,
    % a fixed batch, a non-converged start -> bit-identical sol/f at 1 vs many
    % threads. This is what makes MexRandomSeed actually reproducible.
    methods (Test, TestTags = {'Slow'})
        function threadCountDoesNotChangeResult(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(modelDir, 'ir')));
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            oldCd = pwd; cd(work); tc.addTeardown(@() cd(oldCd));
            gdsge.runtime.ensurePath();
            gdsge.codegen.generateCxx(buildHL1996IR(), work);
            compile_HL1996();
            tc.assertTrue(exist(fullfile(work, ['mex_HL1996.' mexext]), 'file') == 3, ...
                'new MEX did not compile');

            S = gdsgetest.buildHL1996MexInputs();
            % Force restarts: start from the bound midpoint (Pass 0 won't fully converge).
            sol0 = (S.lb + S.ub) / 2;
            S.cfg.useMexResolve = true;
            S.cfg.useMexRandomize = true;
            S.cfg.randomizeBatch = 5;
            S.cfg.randomSeed = 0;
            S.cfg.randomizeSalt = 1;
            S.cfg.maxMinorIter = 5;          % exactly one batch
            S.cfg.useNearestNeighbor = true;

            cfg1 = S.cfg; cfg1.numThreads = 1;
            cfgN = S.cfg; cfgN.numThreads = max(2, feature('numcores'));
            [solA, fA] = gdsge.runtime.solveProblems(@mex_HL1996, sol0, S.lb, S.ub, ...
                S.data, S.f0, S.aux0, S.eqval0, cfg1);
            [solB, fB] = gdsge.runtime.solveProblems(@mex_HL1996, sol0, S.lb, S.ub, ...
                S.data, S.f0, S.aux0, S.eqval0, cfgN);

            tc.verifyTrue(isequal(solA, solB), 'sol differs across NumThreads (RNG not deterministic)');
            tc.verifyTrue(isequal(fA, fB), 'f differs across NumThreads');
            % Sanity: restarts actually fired (midpoint start did not solve all in Pass 0).
            tc.verifyGreaterThan(nnz(fA > S.cfg.tolSol) + 1, 1, ...
                'expected some points to need restarts from the midpoint start');
        end
    end
end
```

- [ ] **Step 2: Run it**

Run: `matlab -batch "cd('tests'); r=runtests('tInMexRandomizeDeterminism'); assert(all([r.Passed])); disp('PASS')"`
Expected: PASS — `solA`/`solB` bit-identical. (If it fails, the RNG depends on thread scheduling — check that the draw is keyed only by `(seed,salt,point,trial,comp)` and the neighbor sweep uses snapshot semantics.)

- [ ] **Step 3: Commit**

```bash
git add tests/HeatonLucas1996/codegen/tInMexRandomizeDeterminism.m
git commit -m "test: in-MEX randomize is independent of NumThreads (determinism gate)"
```

---

## Task 8: A/B + tolerance gates (Slow)

Two gates: (a) `UseMexRandomize=0` still reproduces the OLD goldens bit-exactly (the MATLAB path is unchanged at runtime); (b) `UseMexRandomize=1` lands within tolerance of the OLD golden. safe_assets is the decisive restart-heavy case.

**Files:**
- Create: `tests/HeatonLucas1996/codegen/tInMexRandomizeAB_HL1996.m`
- Create: `tests/Barro_et_al_2017/codegen/tInMexRandomizeSafeAssets.m`

- [ ] **Step 1: Write the HL1996 A/B + tolerance gate**

Create `tests/HeatonLucas1996/codegen/tInMexRandomizeAB_HL1996.m`:

```matlab
classdef tInMexRandomizeAB_HL1996 < matlab.unittest.TestCase
    % (a) UseMexRandomize=0 reproduces the committed OLD golden bit-exactly
    %     (MATLAB randomize path unchanged at runtime).
    % (b) UseMexRandomize=1 lands within tolerance of the same golden.
    methods (Test, TestTags = {'Slow'})
        function offMatchesGoldenBitExact_onMatchesTolerance(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            g = load(fullfile(modelDir, 'golden', 'IterRslt.mat')); G = g.IterRslt;
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'HL1996.gmod'), work);
            oldCd = pwd; cd(work); tc.addTeardown(@() cd(oldCd));
            gdsge_codegen('HL1996');

            % (a) MATLAB path, RNG pinned exactly as the golden was captured.
            rng(0823);
            R0 = iter_HL1996(struct('SaveFreq',inf,'NoSave',1,'UseMexRandomize',0));
            tc.verifyEqual(R0.Iter, G.Iter, 'UseMexRandomize=0 iter count must match golden');
            flds = {'var_policy','var_aux','var_interp'};
            for i = 1:numel(flds)
                fn = fieldnames(G.(flds{i}));
                for k = 1:numel(fn)
                    tc.verifyTrue(isequal(R0.(flds{i}).(fn{k}), G.(flds{i}).(fn{k})), ...
                        sprintf('UseMexRandomize=0 %s.%s must be bit-identical to golden', ...
                        flds{i}, fn{k}));
                end
            end

            % (b) in-MEX randomize, fixed seed -> within tolerance of the golden.
            R1 = iter_HL1996(struct('SaveFreq',inf,'NoSave',1,'UseMexRandomize',1,'MexRandomSeed',0));
            for i = 1:numel(flds)
                fn = fieldnames(G.(flds{i}));
                for k = 1:numel(fn)
                    r = gdsgetest.compareNumericClose(R1.(flds{i}).(fn{k}), G.(flds{i}).(fn{k}), 1e-4, 1e-6);
                    tc.verifyTrue(r.pass, sprintf('UseMexRandomize=1 %s.%s out of tolerance:\n%s', ...
                        flds{i}, fn{k}, strjoin(r.failures, newline)));
                end
            end
        end
    end
end
```

- [ ] **Step 2: Run the HL1996 gate**

Run: `matlab -batch "cd('tests'); r=runtests('tInMexRandomizeAB_HL1996'); assert(all([r.Passed])); disp('PASS')"`
Expected: PASS. (HL1996 has unique roots; the midpoint/seeded path should converge to the same fixed point well within `1e-4`.)

- [ ] **Step 3: Write the safe_assets gate (decisive — actually exercises restarts)**

Create `tests/Barro_et_al_2017/codegen/tInMexRandomizeSafeAssets.m`:

```matlab
classdef tInMexRandomizeSafeAssets < matlab.unittest.TestCase
    % safe_assets is restart-heavy and multi-root: the decisive case.
    % (a) UseMexRandomize=0 reproduces the OLD golden bit-exactly (rng pinned).
    % (b) UseMexRandomize=1 (fixed seed) converges to a valid equilibrium that
    %     matches the golden within tolerance where the root is unique; multi-root
    %     points may differ pointwise but the simulated moments agree (and the
    %     new output becomes the re-baselined golden in Task 9).
    methods (Test, TestTags = {'Slow'})
        function offBitExact_onValidEquilibrium(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            g = load(fullfile(modelDir, 'golden', 'IterRslt.mat')); G = g.IterRslt;
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'safe_assets.gmod'), work);
            oldCd = pwd; cd(work); tc.addTeardown(@() cd(oldCd));
            gdsge_codegen('safe_assets');

            % (a) MATLAB path, unbounded retries (restart-heavy), rng pinned.
            rng('default');
            R0 = iter_safe_assets(struct('SaveFreq',inf,'NoSave',1,'UseMexRandomize',0, ...
                'MaxMinorIter',inf,'DiagnoseMinorIter',inf));
            tc.verifyEqual(R0.Iter, G.Iter, 'UseMexRandomize=0 iter count must match golden');
            tc.verifyEqual(R0.Metric, G.Metric, 'UseMexRandomize=0 metric must match golden');

            % (b) in-MEX randomize, fixed seed, unbounded retries.
            R1 = iter_safe_assets(struct('SaveFreq',inf,'NoSave',1,'UseMexRandomize',1, ...
                'MexRandomSeed',0,'MaxMinorIter',inf,'DiagnoseMinorIter',inf));
            tc.verifyLessThan(R1.Metric, 1e-6, 'UseMexRandomize=1 must converge');
            % Equilibrium closeness: median pointwise agreement is tight even if a
            % few multi-root points differ. (Re-baseline in Task 9 makes it exact.)
            fn = fieldnames(G.var_interp);
            for k = 1:numel(fn)
                d = abs(R1.var_interp.(fn{k})(:) - G.var_interp.(fn{k})(:));
                s = max(1, abs(G.var_interp.(fn{k})(:)));
                tc.verifyLessThan(median(d ./ s), 1e-3, ...
                    sprintf('UseMexRandomize=1 var_interp.%s median rel-diff too large', fn{k}));
            end
        end
    end
end
```

- [ ] **Step 4: Run the safe_assets gate**

Run: `matlab -batch "cd('tests'); r=runtests('tInMexRandomizeSafeAssets'); disp(table(r))"`
Expected: PASS on (a) bit-exact and (b) converges + median rel-diff < 1e-3. If (b)'s median check is too strict because the seeded path lands on a materially different (but valid) equilibrium, **stop and report to the owner** with the observed diffs — per spec §7 this is the multi-root re-baseline judgment call, not an auto-fix.

- [ ] **Step 5: Commit**

```bash
git add tests/HeatonLucas1996/codegen/tInMexRandomizeAB_HL1996.m tests/Barro_et_al_2017/codegen/tInMexRandomizeSafeAssets.m
git commit -m "test: A/B (off=golden) + tolerance (on) gates for in-MEX randomize"
```

---

## Task 9: Re-baseline the iter goldens + full-suite green

With the default now `UseMexRandomize=1`, every cartesian end-to-end gate that compares `IterRslt` bit-exactly will fail where the seeded path changed the output. For each, confirm tolerance-equivalence (Task 8 pattern), then re-capture the golden under the fixed default seed. **Run MATLAB one process at a time.**

**Files:**
- Regenerate (only where output changed): `tests/HeatonLucas1996/golden/IterRslt.mat`, `tests/Barro_et_al_2017/golden/IterRslt.mat`, `tests/Mendoza2010/golden/IterRslt.mat`, `tests/GLSW*/golden/IterRslt.mat`, `tests/Cao2011EZ/golden/IterRslt.mat`, `tests/CaoNie2016/golden/IterRslt.mat`
- Modify the affected end-to-end gates' RNG pinning to also pin `MexRandomSeed` (replace `rng(...)` reliance for the solve restart with `MexRandomSeed=0`; shock-draw `rng(0823)` in simulate stays).

- [ ] **Step 1: Identify which gates actually changed**

Run the full Slow suite and capture failures:

Run: `matlab -batch "cd('tests'); r=runtests('IncludingSubfolders',true); t=table(r); disp(t(~[r.Passed],:))"`
Expected: a list of failing end-to-end gates (clean unique-root models like HL1996 may already pass bit-exact and need no re-baseline; multi-root/restart-heavy ones fail).

- [ ] **Step 2: For each failing end-to-end gate, confirm tolerance-equivalence**

For each failing model, run the new default vs the old golden through `gdsgetest.compareNumericClose` (mirror Task 8). If a model lands on a materially different (but valid, converged) equilibrium at multi-root points, **stop and report to the owner** for sign-off before re-baselining (spec §7).

- [ ] **Step 3: Re-capture the changed goldens under the fixed default seed**

For each confirmed model, regenerate its golden by running the new default and saving `IterRslt`. Use each model's existing capture/test entry point in its `golden/` or `codegen/` folder (read the model's `tEndToEnd*` to see how it loads/saves). Pattern (adapt per model):

```matlab
% from the model's working dir, after gdsge_codegen('<model>')
IterRslt = iter_<model>(struct('SaveFreq',inf,'NoSave',1));   % default => UseMexRandomize=1, seed 0
save(fullfile(modelGoldenDir,'IterRslt.mat'),'IterRslt');
```

- [ ] **Step 4: Pin `MexRandomSeed` in the re-baselined gates**

In each affected `tEndToEnd*`/`tFusedConstruct*`/`tInMexResolve*` gate that previously relied on `rng(...)` before iter to reproduce the golden, add `'MexRandomSeed',0` to the options struct (and keep `UseMexRandomize` default). The solve restart no longer reads MATLAB `rng`; only simulate's shock draws do (`rng(0823)` stays). Note: `tInMexResolve*` A/B gates compare `UseMexResolve` on/off — set `UseMexRandomize=0` in BOTH arms there so they stay a clean neighbor-sweep A/B (the randomize re-baseline is covered by the Task 8 gates).

- [ ] **Step 5: Run the full suite to green**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: exit 0. Inspect `tests/results/junit.xml` for the count.

- [ ] **Step 6: Commit**

```bash
git add tests
git commit -m "test: re-baseline iter goldens + pin MexRandomSeed for UseMexRandomize default"
```

---

## Task 10: Docs + PROGRESS changelog

**Files:**
- Modify: `docs/user-guide.md` (options table — add the three options)
- Modify: `PROGRESS.md` (changelog entry)

- [ ] **Step 1: Document the options in the user guide**

In `docs/user-guide.md`, in the options reference, add rows for `UseMexRandomize` (default 1; in-MEX C++ randomized restart; 0 = MATLAB path), `MexRandomizeBatch` (default 100; minor-iterations per MEX call), `MexRandomSeed` (default 0; fixed RNG seed for reproducibility). Note the iter goldens were re-baselined and that `rng` no longer affects iter restarts (only simulate shock draws).

- [ ] **Step 2: Add the PROGRESS changelog entry**

In `PROGRESS.md`, add a dated entry under Changelog summarizing: in-MEX randomized restart (cartesian iter path), counter-based per-point RNG (thread-independent, fixed seed), batched (default 100) returning to MATLAB for diagnostics/cap/exit, default-on, iter goldens re-baselined after tolerance check, ASG + simulate/init deferred, measured speedup on safe_assets, full suite count, branch `inmex-randomize`.

- [ ] **Step 3: Commit**

```bash
git add docs/user-guide.md PROGRESS.md
git commit -m "docs: document UseMexRandomize options + PROGRESS changelog"
```

---

## Task 11: Perf check (optional, recommended)

Quantify the win on the restart-heavy model the feature targets.

**Files:**
- Create: `tests/perf/inmex_randomize_bench.m` (mirror `tests/perf/inmex_resolve_bench.m`)

- [ ] **Step 1: Write the benchmark**

Read `tests/perf/inmex_resolve_bench.m` and mirror it: time `iter_safe_assets` with `UseMexRandomize=1` (seed 0) vs `UseMexRandomize=0` (rng pinned), assert identical convergence (Iter/Metric) is NOT required (paths differ), report wall-clock of each. Print both timings.

- [ ] **Step 2: Run and record**

Run: `matlab -batch "cd('tests/perf'); inmex_randomize_bench"`
Expected: prints both wall-clocks; the `=1` path faster (fewer MATLAB↔MEX round-trips). Record the numbers in the PROGRESS changelog entry from Task 10.

- [ ] **Step 3: Commit**

```bash
git add tests/perf/inmex_randomize_bench.m
git commit -m "test(perf): in-MEX randomize benchmark on safe_assets"
```

---

## Self-review notes (for the implementer)

- **Spec coverage:** RNG design → Tasks 1,7; batched loop → Task 2; MATLAB driver → Task 3; options/defaults → Task 4; emitIter init+cfg → Task 5; snapshots → Task 6; determinism → Task 7; A/B + tolerance → Task 8; re-baseline → Task 9; docs → Task 10; perf → Task 11. ASG/simulate/init are explicitly untouched (no task).
- **Caller-workspace contract:** every `GDSGE_RANDOMIZE_*` and `GDSGE_RANDOM_SEED` must be a **local** in `solveProblems` (set with `%#ok<NASGU>`), never in a subfunction — same rule as `GDSGE_PROBLEM_STRIDES`.
- **Naming consistency:** cfg fields `useMexRandomize`, `randomizeBatch`, `randomSeed`, `randomizeSalt`; caller vars `GDSGE_RANDOMIZE_BATCH`, `GDSGE_RANDOM_SEED`, `GDSGE_RANDOMIZE_SALT`, `GDSGE_RANDOMIZE_TRIAL_OFFSET`; options `UseMexRandomize`, `MexRandomizeBatch`, `MexRandomSeed`. Use these exact spellings everywhere.
- **Do NOT** set `useMexRandomize` in `emitSimulate`/`emitIterInit`/`emitIterAsg` — those stay on the MATLAB path.
- **Multi-root judgment (Task 8/9):** if a seeded run lands on a materially different equilibrium, stop and get owner sign-off; do not loosen tolerances to force green.
