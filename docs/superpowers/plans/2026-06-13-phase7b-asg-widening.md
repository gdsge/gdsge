# Phase 7b — ASG Widening: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** CaoKS2016 and Bianchi2011_asg (`USE_ASG=1`) run end-to-end through the public API and match freshly captured old-toolbox goldens.

**Architecture:** Capture both goldens first; vendor the public `asg` class + `asg_mex` (cache-gated auto-compile); settle the small front-end/IR delta (Asg options, setup-name whitelist, transition-ref pool); add the ASG C++ interp path (class-handle interop through the shared layout) and the ASG MATLAB emitters (refinement-loop iter, ASG resolve-simulate) over a new `gdsge.runtime.solveProblemsAsg`; close vertically — CaoKS2016 (vanilla driver, core loop) then Bianchi (WarmUp/SkipModelInit/MaxIter/shock-override surface).

**Tech Stack:** MATLAB R2025b (`matlab -batch`, native `matlab.unittest`), existing `gdsge.parser` / `gdsge.ir` / `gdsge.codegen` / `gdsge.runtime` packages, MSVC for MEX.

**Spec:** `docs/superpowers/specs/2026-06-13-phase7b-asg-widening-design.md`

---

## Conventions for every task

- Work from the repo root `D:\refactor_gdsge`.
- `matlab` means `C:\Program Files\MATLAB\R2025b\bin\matlab.exe`. In PowerShell invoke it as
  `& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "<expr>"`.
  (PowerShell 7 / `pwsh` is **not** installed — never use `tests/run.ps1`.)
- Single-test run pattern (exit code 0 = pass):

  ```
  matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','<folder>','<TestFile>.m')); disp(table(results)); exit(any([results.Failed]))"
  ```

- Full suite: `matlab -batch "cd('tests'); run_tests"` — exit 0 = all pass;
  `tests/results/junit.xml` is authoritative (ignore `results.tap`, it appends).
- The MATLAB path policy is sacred: never `addpath` persistently; one toolbox source per
  `matlab -batch` process. Golden capture adds ONLY `base_package/gdsge/source`; everything
  else adds ONLY `src/` (+ `src/kernels` + `tests`).
- If a golden gate fails: use superpowers:systematic-debugging. **Never loosen tolerances.**

## Verified facts this plan builds on (do not re-derive)

- **IR groundwork exists:** `interpMethod ∈ {spline,asg,pchip}` + `asgMaxLevel`/`asgThreshold`
  in the schema (`src/+gdsge/+ir/schema.m:41-53`); validator enforces asg ⊥ `var_tensor` and
  requires `asgMaxLevel` (`validate.m:72-79`); `resolveOptions.m:16-19` already extracts the
  two Asg options; `generateCxx.m:10-14` refuses non-spline with a "Phase 7b" error.
- **The C++ POP sequence is method-agnostic.** Old ASG and spline pop the identical
  `GDSGE_DATA` column: shock_num → params → transitions → shock-value grids → shockIdx →
  state values (`gdsge_parser.m:1279-1317`). Only the MATLAB-side per-problem rows differ:
  ASG packs `GDSGE_evalArrayIdx;GDSGE_evalGrids` instead of expanded tensors
  (`asg_propose_grids_and_solve_template.m:25`). `gdsge.codegen.dataLayout` needs **no
  change**; `emitDataPack` only gains an `asgPack` string.
- **New templates already carry the ASG hooks:** `templates/cxx/mex.tpl.cpp:38-40` has the
  `#ifdef USE_ASG` include of `class_handle.hpp`/`asg_adouble.h`; `task.tpl.cpp:63` declares
  `int GDSGE_EVAL_GRAD_FLAG = 0;`, `INTERP_GET_CODE` is at task scope (line 45),
  `INTERP_GET_THREAD_CODE` at thread scope (line 65). `GET_MX_ARRAY` is defined in
  `include/mm_lite.h:82` and `include/` already holds `asg.h`, `asg_adouble.h`,
  `asg_matlab.h`, `class_handle.hpp` (vendored in Phase 5).
- **Old ASG runs with Broyden-now off, always:** `params_template.m:27` sets
  `GDSGE_USE_BROYDEN_NOW = 0;` once; no ASG template ever flips it. So the ASG solve passes
  `cfg.useBroydenNow = 0` everywhere — exact parity.
- **Old ASG solve cascade** (`iter_solve_problem_asg_template.m`): solve all → warm
  unconverged from the **current**-iteration sol interp (+ adaptive tighten if
  `UseAdaptiveBound`), solve → restore remaining from the **previous**-iteration sol interp
  (+ tighten, **no solve**) → randomized restarts with the `UseAdaptiveBoundInSol`
  hitting-bounds logic, bounded by `MaxMinorIter`. The **init** cascade
  (`iter_solve_init_problem_asg_template.m`) is just: solve all → randomized restarts.
  The **simulate** cascade (`simulate_resolve_asg_template.m:92-116`) is: solve →
  randomized restarts (+ `UseAdaptiveBoundInSol` block) — exactly what
  `gdsge.runtime.solveProblems` does with `useNearestNeighbor=false`.
- **Old ASG IterRslt fields** (`iter_inf_horizon_asg_template.m:97-107` +
  `iter_inf_output_construct_asg_template.m:21-23`): `Metric, MetricVec, Iter, shock_num,
  shock_trans, var_shock, var_state, params, asg_interp_struct, sol_asg_interp_struct,
  var_others` (+ `output_var_index`, `asg_output_struct` when outputs exist). **No**
  `var_policy/var_aux/var_interp/pp/GDSGE_PROB/Metric0/output_interp` — the ASG result
  shape differs from spline by design and is frozen as-is.
- **`asg` class API** (`base_package/gdsge/source/asg.m`): `asg(inputGridsCell, numVec,
  numArray)`; `get_eval_grids(threshold)` → `[evalArrayIdx, evalGrids, evalGridsLength,
  evalGridsUnscaled]`; `push_eval_results_at_valid(value, valid)`;
  `push_eval_results_at_grids(arrayIdx, gridsUnscaled, value, level)`;
  `eval_vec(arrayIdx, sites)`; `get_current_level()` (fresh object → -1);
  `get_grids_info()`, `get_grids_info_at_level(level)`; `convert_to_struct()` / static
  `construct_from_struct(s)`; static `get_mex_constants()` → `[ASG_MAX_DIM, ASG_MAX_NVEC,
  ASG_MAX_LEVEL]`; static `compute_inf_metric(a,b)` → `[metric, metricVec]`; property
  `objectHandle` (uint64 for `GDSGE_ASG_HANDLE`).
- **Old ASG option defaults** (`default_mod.nmod`): `AsgMinLevel=4`, `AsgMaxLevel=10`,
  `AsgThreshold=1e-2`, `AsgOutputMaxLevel=10`, `AsgOutputThreshold=1e-2`.
  `GDSGE_ASG_FIX_GRID` (default 0) and its fixed-grid branch are deliberately **not**
  reproduced (spec §2 — no driver exercises it; the option is not parsed, so a gmod/driver
  setting it errors informatively via the whitelist).
- **Options whitelist** (`emitIter.m:29-39`) already contains params + shock names +
  transition names + state names — Bianchi's `options.yT/yN/shock_trans/b` already pass.
  Only **setup-assigned names** (`bMin`, `bMax`, …) are missing; the eval sandbox `ws` in
  `parseDeclarations` has them all as fieldnames.
- **`parseDeclarations.m:41`** builds `ir.shocks.transitions` as
  `struct('shock_trans', …)` only. CaoKS's `GDSGE_EXPECT{ … | shock_trans2 }` pipes on a
  **parameter** (`shock_trans2`, a 4×4 copy). The transRef checks live at
  `analyzeModel.m:35,124`, `validate.m:199-202` (`buildPools`), and
  `emitModelBody.m:11,50`; the C++ EXPECT lowering (`emitModelBody.m:66-67`) emits
  `<transRef>((shock)+N*(GDSGE_iter-1))` and the param pop already defines
  `#define shock_trans2(idx) …` (`emitPop.m:18-19`) — so widening the three ref pools to
  square-matrix params is the **whole** fix.
- **Binary precedent:** `src/kernels/myppual_mex.mexw64` is committed; `*.cache` is
  git-ignored; `asg_mex.mexw64` does NOT match the `mex_*.mexw64` ignore glob → commit it.
- **Gmod facts:** CaoKS2016: `TolEq=1e-4`, `AsgMaxLevel=10`, `AsgThreshold=1e-4`,
  states `K,X` (21 pts each), 4 shocks, 8 policies, `var_output Kp Xp kp1 kp2`,
  `model_init` with `var_policy_init c1 c2`, interp updates `c1Future = c1; c2Future = c2;`.
  Bianchi2011_asg (`bianchi2011.gmod`): no TolEq (→ default 1e-6), `AsgMaxLevel=14`,
  `AsgThreshold=1e-5`, state `b` (101 pts), 16 shocks (gmod holds **placeholder**
  `shock_trans=zeros(16)`, `yT=yN=ones(1,16)`; the driver overrides from
  `shock_process.mat`), 4 policies, `var_aux c lambda bNext`, `var_policy_init dummy` +
  `var_aux_init c lambda`, init equation literal `0;`, `var_output bNext pN c`.
  Both `var_simu` lists ⊆ `var_output` (no promotion happens). Both use `SIMU_RESOLVE`
  (default). Simulate transitions are unprimed (`K'=Kp`, `X'=Xp`, `b'=bNext`) — the 7a
  `primed` flag handles them. All `inbound`/`inbound_init` bounds in both models are
  setup-constant expressions (no state/shock-dependent bounds on the ASG path).
- **`codeWriter`**: `w.add(fmt,…)` auto-indents per `w.in()`/`w.out()`; `w.addRaw(txt)`
  does NOT — wrap raw fragments with `gdsge.codegen.mat.indentBy(txt, '<spaces>')` (see
  `emitIter.m:248-253`). Verify by reading `src/+gdsge/+codegen/codeWriter.m` before
  writing the emitters.

## File structure

| File | Action | Responsibility |
|---|---|---|
| `tests/CaoKS2016/CaoKS2016.gmod` | Create (copy) | Model source under test |
| `tests/Bianchi2011_asg/bianchi2011.gmod`, `shock_process.mat` | Create (copy) | Model source + shock fixture |
| `tests/golden/capture_CaoKS2016.m`, `capture_bianchi2011.m` | Create | Old-toolbox golden captures |
| `tests/CaoKS2016/tGoldenCaoKS2016.m`, `tests/Bianchi2011_asg/tGoldenBianchi2011.m` | Create | Golden-integrity tests |
| `src/asg.m` | Create (copy) | Frozen public `asg` class (flat) |
| `src/kernels/asg_mex.cpp` | Create (copy) | ASG MEX source (vendored) |
| `src/kernels/asg_mex.mexw64` | Create (compiled) | Committed binary (myppual_mex precedent) |
| `src/+gdsge/+codegen/ensureAsgMex.m` | Create | Cache-gated asg_mex auto-compile |
| `tests/runtime/tAsgClass.m`, `tests/codegen/tEnsureAsgMex.m` | Create | asg smoke + compile-cache tests |
| `src/+gdsge/+parser/resolveOptions.m` | Modify | AsgMinLevel/AsgOutput* options; ASG⊥SIMU_INTERP error |
| `src/+gdsge/+parser/defaultSetupCode.m` | Modify | New Asg option defaults |
| `src/+gdsge/+parser/parseDeclarations.m`, `assemblePartialIR.m` | Modify | `setupNames` into the IR |
| `src/+gdsge/+parser/analyzeModel.m` | Modify | transRef accepts square-matrix params |
| `src/+gdsge/+ir/schema.m` | Modify | `setupNames`; 3 new Asg option fields |
| `src/+gdsge/+ir/validate.m` | Modify | transitions pool += square-matrix params |
| `src/+gdsge/+codegen/+cxx/emitModelBody.m` | Modify | transNames += square-matrix params |
| `src/+gdsge/+codegen/+mat/optionsWhitelist.m` | Create | Shared whitelist literal (+ setupNames) |
| `src/+gdsge/+codegen/+mat/emitIter.m`, `emitSimulate.m` | Modify | Use optionsWhitelist |
| `tests/HeatonLucas1996/ir/buildHL1996IR.m` + JSON/codegen snapshots | Modify/Regen | New IR fields ripple |
| `tests/<7a models>/ir/*.gdsge.json` | Regenerate | setupNames ripple |
| `tests/CaoKS2016/parser/tFrontEndCaoKS2016.m`, `tests/Bianchi2011_asg/parser/tFrontEndBianchi2011.m` | Create | Front-end gates |
| `templates/cxx/interp_asg_construct.tpl.cpp`, `interp_asg_get.tpl.cpp`, `interp_asg_prepare_space.tpl.cpp` | Create (copy) | ASG C++ interp fragments |
| `src/+gdsge/+codegen/+cxx/emitInterp.m` | Modify | asg branch |
| `src/+gdsge/+codegen/+cxx/emitCompile.m` | Modify | `-DUSE_ASG -DASG_MAX_*` defines |
| `src/+gdsge/+codegen/generateCxx.m` | Modify | Unlock asg; ensureAsgMex; pchip/tensor → "7c" |
| `src/+gdsge/+runtime/solveProblems.m` | Modify | `GDSGE_ASG_HANDLE` contract local |
| `src/+gdsge/+runtime/solveProblemsAsg.m` | Create | ASG solve cascade |
| `tests/runtime/tSolveProblemsAsg.m` | Create | Cascade unit test (mock MEX) |
| `src/+gdsge/+codegen/+mat/emitBounds.m` | Modify | Optional ASG name-rewrite mode |
| `src/+gdsge/+codegen/+mat/emitDataPack.m` | Modify | `frag.asgPack` |
| `src/+gdsge/+codegen/+mat/emitIterAsg.m` | Create | ASG iter (init + refinement + output + IterRslt) |
| `src/+gdsge/+codegen/+mat/emitSimulateAsg.m` | Create | ASG resolve-simulate |
| `src/+gdsge/+codegen/generateMatlab.m` | Modify | asg branch |
| `tests/+gdsgefix/minimalIRAsg.m`, `tests/codegen/tEmitAsg.m` | Create | Emitter unit tests |
| `tests/CaoKS2016/codegen/tEndToEndCaoKS2016.m`, `tests/Bianchi2011_asg/codegen/tEndToEndBianchi2011.m` | Create | End-to-end golden gates (Slow) |
| `PROGRESS.md` | Modify | Check off 7b |

---

### Task 1: Branch setup

**Files:** none (git only)

- [ ] **Step 1: Create the phase branch**

```bash
git checkout -b phase7b-asg-widening
```

Expected: `Switched to a new branch 'phase7b-asg-widening'`. (`git status` must be clean
first; if not, stop and report.)

---

### Task 2: Golden capture — CaoKS2016 + Bianchi2011_asg

**Files:**
- Create: `tests/CaoKS2016/CaoKS2016.gmod`, `tests/Bianchi2011_asg/bianchi2011.gmod`, `tests/Bianchi2011_asg/shock_process.mat` (verbatim copies)
- Create: `tests/golden/capture_CaoKS2016.m`, `tests/golden/capture_bianchi2011.m`
- Create: `tests/CaoKS2016/tGoldenCaoKS2016.m`, `tests/Bianchi2011_asg/tGoldenBianchi2011.m`

- [ ] **Step 1: Copy the model sources into the test tree**

```bash
mkdir -p tests/CaoKS2016/golden tests/Bianchi2011_asg/golden
cp base_package/gdsge/tests/CaoKS2016/CaoKS2016.gmod tests/CaoKS2016/
cp base_package/gdsge/tests/Bianchi2011_asg/bianchi2011.gmod tests/Bianchi2011_asg/
cp base_package/gdsge/tests/Bianchi2011_asg/shock_process.mat tests/Bianchi2011_asg/
```

- [ ] **Step 2: Sanity-check the OLD toolbox's asg_mex** (ships prebuilt; must load in R2025b)

```
matlab -batch "addpath(fullfile('base_package','gdsge','source')); [d,v,l]=asg.get_mex_constants(); fprintf('ASG_MAX_DIM=%d NVEC=%d LEVEL=%d\n',d,v,l)"
```

Expected: three positive integers. If this errors (binary incompatible), recompile it in
place: `matlab -batch "cd(fullfile('base_package','gdsge','source')); compile_asg"` and
re-run the check. (base_package is git-ignored; modifying it is fine.)

- [ ] **Step 3: Write the capture scripts**

Create `tests/golden/capture_CaoKS2016.m`:

```matlab
function capture_CaoKS2016()
% Golden capture for CaoKS2016 from the OLD toolbox. Protocol per the Phase 7b
% spec: full-convergence iter exactly as test.m (PrintFreq=10, SaveFreq=100;
% gmod TolEq=1e-4), reduced seeded simulate (test.m runs the gmod default
% 100x10000; golden uses 6x1000 — SIMU_RESOLVE re-solves every period).
here     = fileparts(mfilename('fullpath'));          % tests/golden
repoRoot = fileparts(fileparts(here));
oldSrc   = fullfile(repoRoot, 'base_package', 'gdsge', 'source');
modelDir = fullfile(repoRoot, 'tests', 'CaoKS2016');
goldenDir = fullfile(modelDir, 'golden');

work = tempname; mkdir(work);
copyfile(fullfile(modelDir, 'CaoKS2016.gmod'), work);

oldPath = path; restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(oldSrc);
oldCd = pwd; cdRestore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(work);

t0 = tic;
gdsge_codegen('CaoKS2016');
IterOptions.PrintFreq = 10;
IterOptions.SaveFreq = 100;
IterRslt = iter_CaoKS2016(IterOptions);
fprintf('iter wall-clock: %.1fs (Iter=%d)\n', toc(t0), IterRslt.Iter);

t0 = tic;
rng(0823);
SimuRslt = simulate_CaoKS2016(IterRslt, struct('num_samples', 6, 'num_periods', 1000));
fprintf('simulate wall-clock: %.1fs\n', toc(t0));

save(fullfile(goldenDir, 'IterRslt.mat'), 'IterRslt', '-v7');
save(fullfile(goldenDir, 'SimuRslt.mat'), 'SimuRslt', '-v7');
fprintf('Golden captured: Iter=%d Metric=%g\n', IterRslt.Iter, IterRslt.Metric);
end
```

Create `tests/golden/capture_bianchi2011.m`:

```matlab
function capture_bianchi2011()
% Golden capture for Bianchi2011_asg from the OLD toolbox. The FULL two-stage
% driver sequence is captured as-is (spec 7b: the staging itself is surface
% under test): stage 1 = MaxIter=50 with runtime shock overrides from
% shock_process.mat; stage 2 = WarmUp re-solve to convergence with widened
% state bounds and SkipModelInit=1. Reduced seeded simulate 6x1000.
here     = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(here));
oldSrc   = fullfile(repoRoot, 'base_package', 'gdsge', 'source');
modelDir = fullfile(repoRoot, 'tests', 'Bianchi2011_asg');
goldenDir = fullfile(modelDir, 'golden');

work = tempname; mkdir(work);
copyfile(fullfile(modelDir, 'bianchi2011.gmod'), work);
copyfile(fullfile(modelDir, 'shock_process.mat'), work);

oldPath = path; restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(oldSrc);
oldCd = pwd; cdRestore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(work);

t0 = tic;
gdsge_codegen('bianchi2011');
options = struct;
shock_process = load('shock_process.mat');
options.shock_trans = shock_process.shock_trans;
options.yT = shock_process.yT;
options.yN = shock_process.yN;
options.MaxIter = 50;
IterRslt1 = iter_bianchi2011(options);
fprintf('stage-1 iter wall-clock: %.1fs (Iter=%d, Metric=%g)\n', ...
    toc(t0), IterRslt1.Iter, IterRslt1.Metric);

t0 = tic;
options.MaxIter = inf;
options.WarmUp = IterRslt1;
options.bMin = -1.1;
options.bMax = 0.0;
options.b = [options.bMin, options.bMax];
options.SkipModelInit = 1;
IterRslt = iter_bianchi2011(options);
fprintf('stage-2 iter wall-clock: %.1fs (Iter=%d)\n', toc(t0), IterRslt.Iter);

t0 = tic;
rng(0823);
SimuRslt = simulate_bianchi2011(IterRslt, struct('num_samples', 6, 'num_periods', 1000));
fprintf('simulate wall-clock: %.1fs\n', toc(t0));

save(fullfile(goldenDir, 'IterRslt1.mat'), 'IterRslt1', '-v7');
save(fullfile(goldenDir, 'IterRslt.mat'),  'IterRslt',  '-v7');
save(fullfile(goldenDir, 'SimuRslt.mat'),  'SimuRslt',  '-v7');
fprintf('Golden captured: Iter=%d Metric=%g\n', IterRslt.Iter, IterRslt.Metric);
end
```

- [ ] **Step 4: Run the two captures (old toolbox; wall-clock unknown — log it)**

```
matlab -batch "cd(fullfile('tests','golden')); capture_CaoKS2016"
matlab -batch "cd(fullfile('tests','golden')); capture_bianchi2011"
```

Expected: each prints `Golden captured: Iter=<n> Metric=<m>` with the final Metric below
TolEq (1e-4 CaoKS, 1e-6 Bianchi). **Record the printed Iter values, the stage-1 Bianchi
Metric, and wall-clock times** — gates and PROGRESS.md need them. Check golden sizes
(`ls -la tests/*/golden/`): ASG structs should be a few MB at most; if any file exceeds
~50 MB, stop and reconsider before committing. If a run is prohibitively long (hours),
STOP and escalate per spec §10.

- [ ] **Step 5: Write the golden-integrity tests**

Create `tests/CaoKS2016/tGoldenCaoKS2016.m`:

```matlab
classdef tGoldenCaoKS2016 < matlab.unittest.TestCase
    % Phase 7b: golden files exist, load, and have the expected ASG shape.
    methods (Test)
        function goldenLoadsAndHasShape(tc)
            here = fileparts(mfilename('fullpath'));
            g = load(fullfile(here, 'golden', 'IterRslt.mat'));
            R = g.IterRslt;
            tc.verifyEqual(R.shock_num, 4);
            tc.verifyLessThan(R.Metric, 1e-4);
            tc.verifyTrue(isfield(R, 'MetricVec'));
            tc.verifyEqual(numel(R.var_state.K), 21);
            tc.verifyEqual(numel(R.var_state.X), 21);
            tc.verifyEqual(R.asg_interp_struct.numDim, 2);
            tc.verifyEqual(R.asg_interp_struct.numVec, 2);     % c1Future,c2Future
            tc.verifyEqual(R.asg_interp_struct.numArray, 4);
            tc.verifyEqual(R.sol_asg_interp_struct.numVec, 8); % 8 policies
            tc.verifyEqual(R.asg_output_struct.numVec, 4);     % Kp,Xp,kp1,kp2
            tc.verifyEqual(R.output_var_index.Kp, 1);
            % the public asg class round-trips the struct
            A = asg.construct_from_struct(R.asg_interp_struct);
            grids = A.get_grids_info();
            tc.verifyEqual(numel(grids), 4);
            s = load(fullfile(here, 'golden', 'SimuRslt.mat'));
            tc.verifyEqual(size(s.SimuRslt.shock), [6, 1001]);
            tc.verifyTrue(isfield(s.SimuRslt, 'kp1'));
            tc.verifyTrue(isfield(s.SimuRslt, 'K'));
        end
    end
end
```

Create `tests/Bianchi2011_asg/tGoldenBianchi2011.m`:

```matlab
classdef tGoldenBianchi2011 < matlab.unittest.TestCase
    methods (Test)
        function goldenLoadsAndHasShape(tc)
            here = fileparts(mfilename('fullpath'));
            g1 = load(fullfile(here, 'golden', 'IterRslt1.mat'));
            tc.verifyEqual(g1.IterRslt1.Iter, 50);             % MaxIter=50 stage
            g = load(fullfile(here, 'golden', 'IterRslt.mat'));
            R = g.IterRslt;
            tc.verifyEqual(R.shock_num, 16);
            tc.verifyLessThan(R.Metric, 1e-6);
            tc.verifyEqual(numel(R.var_state.b), 2);           % options.b override
            tc.verifyEqual(R.asg_interp_struct.numDim, 1);
            tc.verifyEqual(R.asg_interp_struct.numVec, 1);     % lambda_interp
            tc.verifyEqual(R.asg_interp_struct.numArray, 16);
            tc.verifyEqual(R.sol_asg_interp_struct.numVec, 4); % nbNext,mu,cT,pN
            tc.verifyEqual(R.asg_output_struct.numVec, 3);     % bNext,pN,c
            % the overridden (non-placeholder) shock process is in the result
            tc.verifyGreaterThan(max(abs(R.shock_trans(:))), 0);
            tc.verifyGreaterThan(max(abs(R.var_shock.yT - 1)), 1e-6);
            s = load(fullfile(here, 'golden', 'SimuRslt.mat'));
            tc.verifyEqual(size(s.SimuRslt.shock), [6, 1001]);
            tc.verifyTrue(isfield(s.SimuRslt, 'c'));
            tc.verifyTrue(isfield(s.SimuRslt, 'pN'));
            tc.verifyTrue(isfield(s.SimuRslt, 'b'));
        end
    end
end
```

- [ ] **Step 6: Run the integrity tests** (need the OLD source on the path for `asg`
  until Task 3 vendors it — run them with the old source, this once):

```
matlab -batch "addpath(fullfile('base_package','gdsge','source')); results = [runtests(fullfile('tests','CaoKS2016','tGoldenCaoKS2016.m')), runtests(fullfile('tests','Bianchi2011_asg','tGoldenBianchi2011.m'))]; disp(table(results)); exit(any([results.Failed]))"
```

Expected: PASS (2/2). If a shape assertion fails, the capture diverged from this plan's
assumptions — inspect the golden, fix the assertion or the capture, do not proceed blind.
(After Task 3 these tests run under `src/` like everything else; nothing to change in the
test files themselves.)

- [ ] **Step 7: Commit**

```bash
git add tests/CaoKS2016 tests/Bianchi2011_asg tests/golden
git commit -m "test(golden): capture CaoKS2016 + Bianchi2011_asg ASG goldens from the old toolbox"
```

---

### Task 3: Vendor the public ASG runtime + cache-gated asg_mex compile

**Files:**
- Create: `src/asg.m`, `src/kernels/asg_mex.cpp` (verbatim copies), `src/+gdsge/+codegen/ensureAsgMex.m`
- Test: `tests/runtime/tAsgClass.m`, `tests/codegen/tEnsureAsgMex.m`

- [ ] **Step 1: Copy the vendored sources** (as-is; parent spec: kernels move, not rewritten)

```bash
cp base_package/gdsge/source/asg.m src/asg.m
cp base_package/gdsge/source/asg_mex.cpp src/kernels/asg_mex.cpp
```

- [ ] **Step 2: Write the failing tests**

Create `tests/codegen/tEnsureAsgMex.m`:

```matlab
classdef tEnsureAsgMex < matlab.unittest.TestCase
    % ensureAsgMex compiles src/kernels/asg_mex.cpp once (hash cache beside
    % it, asg_mex.cache) and skips when source+cache are unchanged.
    methods (Test)
        function compilesOnceThenSkips(tc)
            gdsge.codegen.ensureAsgMex();
            mexFile = fullfile(kernelsDir(), ['asg_mex.' mexext]);
            tc.assertTrue(exist(mexFile, 'file') == 3);
            before = dir(mexFile);
            gdsge.codegen.ensureAsgMex();        % second call: skip path
            after = dir(mexFile);
            tc.verifyEqual(after.datenum, before.datenum, ...
                'second ensureAsgMex call recompiled despite unchanged source');
            [d, v, l] = asg.get_mex_constants();
            tc.verifyGreaterThan(d, 0);
            tc.verifyGreaterThan(v, 0);
            tc.verifyGreaterThan(l, 0);
        end
    end
end

function k = kernelsDir()
k = fileparts(which('asg_mex.cpp'));
end
```

Create `tests/runtime/tAsgClass.m`:

```matlab
classdef tAsgClass < matlab.unittest.TestCase
    % Smoke tests for the vendored public asg class over asg_mex: refinement
    % round-trip on a known smooth function, struct serialization, metric.
    methods (TestClassSetup)
        function compileMex(tc) %#ok<MANU>
            gdsge.codegen.ensureAsgMex();
        end
    end
    methods (Test)
        function refinesAndEvaluatesKnownFunction(tc)
            % f(x) = x.^2 on [0,1], 1 vector, 2 shock arrays
            A = asg({linspace(0, 1, 11)}, 1, 2);
            tc.verifyEqual(A.get_current_level(), -1);
            while A.get_current_level() < 6
                [idx, grids] = A.get_eval_grids(0.0);
                if isempty(grids); break; end
                A.push_eval_results(grids.^2);   % same f for both arrays
            end
            sites = [0.05 0.3 0.62 0.9];
            for arrayIdx = 1:2
                v = A.eval_vec(arrayIdx*ones(1, numel(sites)), sites);
                tc.verifyLessThan(max(abs(v - sites.^2)), 1e-3);
            end
        end
        function structRoundTripPreservesEvaluation(tc)
            A = asg({linspace(-1, 1, 5)}, 1, 1);
            while A.get_current_level() < 5
                [~, grids] = A.get_eval_grids(0.0);
                if isempty(grids); break; end
                A.push_eval_results(sin(grids));
            end
            s = A.convert_to_struct();
            tc.verifyEqual(s.numDim, 1);
            tc.verifyEqual(s.numVec, 1);
            B = asg.construct_from_struct(s);
            sites = linspace(-0.95, 0.95, 7);
            va = A.eval_vec(ones(1, 7), sites);
            vb = B.eval_vec(ones(1, 7), sites);
            tc.verifyEqual(va, vb);
            [metric, metricVec] = asg.compute_inf_metric(A, B);
            tc.verifyEqual(metric, 0);
            tc.verifyEqual(numel(metricVec), 1);
        end
    end
end
```

- [ ] **Step 3: Run them — verify they fail** (`ensureAsgMex` undefined):

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','codegen','tEnsureAsgMex.m')); disp(table(results)); exit(any([results.Failed]))"
```

- [ ] **Step 4: Implement `ensureAsgMex`**

Create `src/+gdsge/+codegen/ensureAsgMex.m`:

```matlab
function ensureAsgMex()
% ENSUREASGMEX  Compile src/kernels/asg_mex.cpp when needed (hash cache
%   asg_mex.cache beside it, same mechanism as mex_<model>.cache). Needed
%   before ASG codegen: emitCompile reads asg.get_mex_constants(). Deliberate
%   improvement over the old toolbox's ship-prebuilt + manual compile_asg.m;
%   invisible to callers. Compile flags mirror compile_asg.m (minus -v).
here      = fileparts(mfilename('fullpath'));   % src/+gdsge/+codegen
srcRoot   = fileparts(fileparts(here));         % src
repoRoot  = fileparts(srcRoot);
kernels   = fullfile(srcRoot, 'kernels');
cppFile   = fullfile(kernels, 'asg_mex.cpp');
cacheFile = fullfile(kernels, 'asg_mex.cache');
mexFile   = fullfile(kernels, ['asg_mex.' mexext]);
cppText   = fileread(cppFile);
if exist(mexFile, 'file') == 3 && ~gdsge.codegen.needsCompile(cppText, cacheFile)
    return;
end
fprintf('Compiling asg_mex (cache-gated):\n');
includeDir = fullfile(repoRoot, 'include');
oldCd = pwd; restore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(kernels);
mex(['-I' includeDir], 'asg_mex.cpp', ...
    'OPTIMFLAGS=/O2 /DNDEBUG', ...
    'COMPFLAGS=$COMPFLAGS /wd4267 /wd4068 /wd4091 /openmp');
gdsge.codegen.writeText(cacheFile, cppText);
end
```

- [ ] **Step 5: Run both test files — verify PASS**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','codegen','tEnsureAsgMex.m')), runtests(fullfile('tests','runtime','tAsgClass.m'))]; disp(table(results)); exit(any([results.Failed]))"
```

If `push_eval_results` (used by tAsgClass) does not exist with that exact name on the
vendored class, check `src/asg.m` for the actual method (the explorer found
`push_eval_results(value)` for the grids from the last `get_eval_grids`) and fix the TEST
— never edit the vendored class.

- [ ] **Step 6: Commit** (include the compiled `asg_mex.mexw64` — committed-binary
  precedent `src/kernels/myppual_mex.mexw64`; `asg_mex.cache` is ignored by `*.cache`)

```bash
git add src/asg.m src/kernels/asg_mex.cpp src/kernels/asg_mex.mexw64 src/+gdsge/+codegen/ensureAsgMex.m tests/runtime/tAsgClass.m tests/codegen/tEnsureAsgMex.m
git commit -m "feat(runtime): vendor public asg class + cache-gated asg_mex compile"
```

---

### Task 4: Front-end/IR — Asg options, setupNames, transRef pool, front-end gates

**Files:**
- Modify: `src/+gdsge/+parser/resolveOptions.m`, `defaultSetupCode.m`, `parseDeclarations.m`, `assemblePartialIR.m`, `analyzeModel.m`
- Modify: `src/+gdsge/+ir/schema.m`, `src/+gdsge/+ir/validate.m`, `src/+gdsge/+codegen/+cxx/emitModelBody.m`
- Create: `src/+gdsge/+codegen/+mat/optionsWhitelist.m`
- Modify: `src/+gdsge/+codegen/+mat/emitIter.m`, `emitSimulate.m`
- Modify/Regen: `tests/HeatonLucas1996/ir/buildHL1996IR.m`, all committed `*.gdsge.json` snapshots, HL1996 codegen snapshots, `docs/ir-schema.md`
- Create: `tests/CaoKS2016/parser/tFrontEndCaoKS2016.m`, `tests/Bianchi2011_asg/parser/tFrontEndBianchi2011.m`
- Test: additions to `tests/parser/tResolveOptions.m`; new `tests/parser/tTransRefParam.m`, `tests/codegen/tOptionsWhitelist.m`

- [ ] **Step 1: Write the failing unit tests**

Create `tests/codegen/tOptionsWhitelist.m`:

```matlab
classdef tOptionsWhitelist < matlab.unittest.TestCase
    % The shared GDSGE_OPTIONS_VALID literal: frozen options + model names +
    % gmod setup-assigned names (old semantics: options override anything the
    % declaration region defined — Bianchi passes bMin/bMax).
    methods (Test)
        function includesSetupAndModelNames(tc)
            ir = gdsgefix.minimalIR();
            ir.setupNames = {'beta','kMin','kMax'};
            lit = gdsge.codegen.mat.optionsWhitelist(ir, {'init'});
            tc.verifyTrue(contains(lit, '''kMin'''));
            tc.verifyTrue(contains(lit, '''kMax'''));
            tc.verifyTrue(contains(lit, '''WarmUp'''));
            tc.verifyTrue(contains(lit, '''SkipModelInit'''));
            tc.verifyTrue(contains(lit, '''init'''));
            % no duplicates even when a param is also a setup name
            tc.verifyEqual(numel(strfind(lit, '''beta''')), 1);
        end
        function setupNamesOptional(tc)
            ir = gdsgefix.minimalIR();
            if isfield(ir, 'setupNames'); ir = rmfield(ir, 'setupNames'); end
            lit = gdsge.codegen.mat.optionsWhitelist(ir, {});
            tc.verifyTrue(contains(lit, '''TolEq'''));
        end
    end
end
```

(`minimalIR`'s param set may not literally contain `beta` — read the fixture and pin the
duplicate-check to one of its actual param names.)

Add to `tests/parser/tResolveOptions.m` inside `methods (Test)`:

```matlab
        function asgOptionsFlowThrough(tc)
            ws = struct('USE_ASG', 1, 'USE_SPLINE', 0, ...
                'AsgMaxLevel', 14, 'AsgThreshold', 1e-5);
            o = gdsge.parser.resolveOptions(ws);
            tc.verifyEqual(o.interpMethod, 'asg');
            tc.verifyEqual(o.asgMinLevel, 4);            % old default_mod.nmod default
            tc.verifyEqual(o.asgMaxLevel, 14);
            tc.verifyEqual(o.asgThreshold, 1e-5);
            tc.verifyEqual(o.asgOutputMaxLevel, 10);
            tc.verifyEqual(o.asgOutputThreshold, 1e-2);
        end
        function asgRejectsSimuInterp(tc)
            ws = struct('USE_ASG', 1, 'USE_SPLINE', 0, ...
                'SIMU_INTERP', 1, 'SIMU_RESOLVE', 0);
            tc.verifyError(@() gdsge.parser.resolveOptions(ws), ...
                'gdsge:parser:asgSimuInterp');
        end
```

Create `tests/parser/tTransRefParam.m`:

```matlab
classdef tTransRefParam < matlab.unittest.TestCase
    % CaoKS2016 surface: a reduction may pipe on a PARAMETER whose value is a
    % shock_num x shock_num matrix (GDSGE_EXPECT{ ... | shock_trans2 }). The
    % transRef pools in analyzeModel / validate / emitModelBody accept it; a
    % non-square parameter still errors.
    methods (Test)
        function squareParamPipeAccepted(tc)
            ir = gdsge.parser.parseFrontEnd(gmodText('trans2'), 'tiny');
            tc.verifyTrue(gdsge.ir.validate(ir));
            red = ir.model.statements{cellfun(@(s) strcmp(s.type,'reduction'), ...
                ir.model.statements)};
            tc.verifyEqual(red.transRef, 'trans2');
        end
        function nonSquareParamPipeErrors(tc)
            tc.verifyError(@() gdsge.parser.parseFrontEnd(gmodText('beta'), 'tiny'), ...
                'gdsge:parser:badTransRef');
        end
    end
end

function txt = gmodText(pipeName)
txt = sprintf([ ...
    'parameters beta trans2;\n' ...
    'beta = 0.95;\n' ...
    'var_shock z;\n' ...
    'shock_num = 2;\n' ...
    'shock_trans = [0.9 0.1; 0.1 0.9];\n' ...
    'trans2 = shock_trans;\n' ...
    'z = [0.9 1.1];\n' ...
    'var_state k;\n' ...
    'k = linspace(0.5,1.5,3);\n' ...
    'var_policy c;\n' ...
    'inbound c 0 2;\n' ...
    'var_interp vf;\n' ...
    'initial vf k;\n' ...
    'vf = c;\n' ...
    'model;\n' ...
    '  cf'' = vf''(k);\n' ...
    '  ev = GDSGE_EXPECT{ cf'' | ' pipeName ' };\n' ...
    '  equations;\n' ...
    '    c - beta*ev;\n' ...
    '  end;\n' ...
    'end;\n' ...
    'simulate;\n' ...
    '  num_periods = 5;\n' ...
    '  num_samples = 1;\n' ...
    '  initial k 1;\n' ...
    '  initial shock 1;\n' ...
    '  var_simu c;\n' ...
    '  k'' = c;\n' ...
    'end;\n']);
end
```

(If the existing parser raises a different error ID for an unknown pipe, read the actual
ID from `analyzeModel.m` and align the test to it — the ID there today is
`gdsge:parser:badTransRef`.)

- [ ] **Step 2: Run all three — verify the new tests fail**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','parser','tResolveOptions.m')), runtests(fullfile('tests','parser','tTransRefParam.m')), runtests(fullfile('tests','codegen','tOptionsWhitelist.m'))]; disp(table(results)); exit(any([results.Failed]))"
```

- [ ] **Step 3: Implement the option + pool changes**

`src/+gdsge/+parser/resolveOptions.m` — replace the asg block (lines 16-19) with:

```matlab
if strcmp(method, 'asg')
    o.asgMinLevel        = getf(ws, 'AsgMinLevel', 4);
    o.asgMaxLevel        = getf(ws, 'AsgMaxLevel', 10);
    o.asgThreshold       = getf(ws, 'AsgThreshold', 1e-2);
    o.asgOutputMaxLevel  = getf(ws, 'AsgOutputMaxLevel', 10);
    o.asgOutputThreshold = getf(ws, 'AsgOutputThreshold', 1e-2);
end
```

and after the existing simuResolve/simuInterp exclusivity check add:

```matlab
if strcmp(method, 'asg') && o.simuInterp ~= 0
    error('gdsge:parser:asgSimuInterp', ...
        'SIMU_INTERP with USE_ASG is not supported yet (deferred to Phase 7c).');
end
```

`src/+gdsge/+parser/defaultSetupCode.m` — next to the existing `AsgMaxLevel`/`AsgThreshold`
defaults add:

```matlab
    'AsgMinLevel = 4;', ...
    'AsgOutputMaxLevel = 10;', ...
    'AsgOutputThreshold = 1e-2;', ...
```

`src/+gdsge/+ir/schema.m` — in the `options` spec, after `'asgThreshold', opt(fScalar()), ...`:

```matlab
    'asgMinLevel',        opt(fScalar()), ...
    'asgOutputMaxLevel',  opt(fScalar()), ...
    'asgOutputThreshold', opt(fScalar()), ...
```

and in `s.root`, after `'setup', opt(fText()), ...`:

```matlab
    'setupNames', opt(fList(fText())), ...
```

`src/+gdsge/+parser/parseDeclarations.m` — add to the `out = struct(...)` literal (keep
the `{}` wrapping for cell-valued fields):

```matlab
    'setupNames', {reshape(setdiff(fieldnames(ws), {'ans'}, 'stable'), 1, [])}, ...
```

`src/+gdsge/+parser/assemblePartialIR.m` — after `ir.setup = decl.setupText;`:

```matlab
ir.setupNames = decl.setupNames;
```

`src/+gdsge/+ir/validate.m` — in `buildPools` (line ~200), extend the transitions pool:

```matlab
pools.transitions = {};
if isfield(ir,'shocks') && isfield(ir.shocks,'transitions') && isstruct(ir.shocks.transitions)
    pools.transitions = fieldnames(ir.shocks.transitions);
end
% CaoKS2016 surface: square shock_num x shock_num parameters are legal
% reduction-pipe targets (the C++ pop defines <name>(idx) for array params).
if isfield(ir,'params') && isfield(ir,'shocks') && isfield(ir.shocks,'count')
    n = ir.shocks.count;
    for iP = 1:numel(ir.params)
        if numel(ir.params{iP}.value) == n*n
            pools.transitions{end+1} = ir.params{iP}.name; %#ok<AGROW>
        end
    end
end
```

`src/+gdsge/+parser/analyzeModel.m` — both transRef checks (lines 35 and 124) become:

```matlab
            if ~(isfield(ir.shocks.transitions, s.transRef) || isSquareParam(ir, s.transRef))
```

and append this local function at the end of the file:

```matlab
function tf = isSquareParam(ir, name)
% A parameter whose value is shock_num x shock_num is a legal reduction
% transition (CaoKS2016: GDSGE_EXPECT{ ... | shock_trans2 }).
tf = false;
n = ir.shocks.count;
for i = 1:numel(ir.params)
    if strcmp(ir.params{i}.name, name) && numel(ir.params{i}.value) == n*n
        tf = true;
        return;
    end
end
end
```

`src/+gdsge/+codegen/+cxx/emitModelBody.m` — line 11 becomes:

```matlab
transNames = fieldnames(ir.shocks.transitions);
for iP = 1:numel(ir.params)
    if numel(ir.params{iP}.value) == ir.shocks.count^2
        transNames{end+1} = ir.params{iP}.name; %#ok<AGROW>
    end
end
```

(No change to the EXPECT emission itself — `<transRef>((shock)+N*(GDSGE_iter-1))` works
for params because `emitPop` already defines `#define <name>(idx)` for array params.)

- [ ] **Step 4: Extract the shared options whitelist (+ setupNames)**

Create `src/+gdsge/+codegen/+mat/optionsWhitelist.m`:

```matlab
function lit = optionsWhitelist(ir, extra)
% OPTIONSWHITELIST  The GDSGE_OPTIONS_VALID cell literal shared by all
%   generated entry points: frozen option surface + model names + gmod
%   setup-assigned names (old semantics: GDSGE_OPTIONS may override anything
%   the declaration-region code defined, e.g. Bianchi's bMin/bMax).
base = {'TolEq','TolSol','TolFun','PrintFreq','NoPrint','SaveFreq','NoSave', ...
    'SimuPrintFreq','SimuSaveFreq','MaxIter','MaxMinorIter','num_samples', ...
    'num_periods','SolMaxIter','UseBroyden','FiniteDiffDelta','GDSGE_USE_BROYDEN', ...
    'GDSGE_DEBUG_EVAL_ONLY','INTERP_ORDER','EXTRAP_ORDER','OutputInterpOrder', ...
    'IterSaveAll','SkipModelInit','UseAdaptiveBound','UseAdaptiveBoundInSol', ...
    'EnforceSimuStateInbound','REUSE_WARMUP_SOL','INTERP_WARMUP_SOL', ...
    'CONSTRUCT_OUTPUT','NumThreads','WarmUp'};
modelNames = {};
for i = 1:numel(ir.params); modelNames{end+1} = ir.params{i}.name; end %#ok<AGROW>
modelNames = [modelNames, ir.shocks.names, fieldnames(ir.shocks.transitions)', ir.states.names];
setupNames = {};
if isfield(ir, 'setupNames'); setupNames = reshape(ir.setupNames, 1, []); end
wl = unique([base, modelNames, setupNames, reshape(extra, 1, [])], 'stable');
lit = ['{''' strjoin(wl, ''',''') '''}'];
end
```

In `src/+gdsge/+codegen/+mat/emitIter.m`, delete the inline `base`/`modelNames`/`wl`
block (lines 29-40) and replace the `GDSGE_OPTIONS_VALID` line with:

```matlab
w.add('GDSGE_OPTIONS_VALID = %s;', gdsge.codegen.mat.optionsWhitelist(ir, {}));
```

In `src/+gdsge/+codegen/+mat/emitSimulate.m`, same replacement (lines 44-54) with:

```matlab
w.add('GDSGE_OPTIONS_VALID = %s;', gdsge.codegen.mat.optionsWhitelist(ir, {'init','GEN_SHOCK_START_PERIOD'}));
```

- [ ] **Step 5: Ripple the new IR fields through the reference artifacts**

1. `tests/HeatonLucas1996/ir/buildHL1996IR.m`: add `ir.setupNames = {...}` right after
   `ir.setup`. Obtain the exact list from the parser:

```
matlab -batch "addpath('src'); ir = gdsge.parser.parseFrontEnd(fileread(fullfile('tests','HeatonLucas1996','HL1996.gmod')),'HL1996'); fprintf('''%s'',', ir.setupNames{:}); fprintf('\n')"
```

   Paste the printed list as the cell literal (drop the trailing comma).

2. Regenerate the HL1996 reference JSON + the 7a model IR snapshots + the schema doc:

```
matlab -batch "addpath('src'); addpath(fullfile('tests','HeatonLucas1996','ir')); ir = buildHL1996IR(); fid=fopen(fullfile('tests','HeatonLucas1996','ir','HL1996.gdsge.json'),'w'); fwrite(fid, gdsge.ir.encode(ir)); fclose(fid);"
matlab -batch "addpath('src'); for m = {{'Barro_et_al_2017','safe_assets'},{'Mendoza2010','mendoza2010'},{'GLSW2020','GLSW_interp'}}; d=m{1}; ir = gdsge.parser.parseFrontEnd(fileread(fullfile('tests',d{1},[d{2} '.gmod'])), d{2}); fid=fopen(fullfile('tests',d{1},'ir',[d{2} '.gdsge.json']),'w'); fwrite(fid, gdsge.ir.encode(ir)); fclose(fid); end"
matlab -batch "addpath('src'); gdsge.ir.gendoc(fullfile('docs','ir-schema.md'));"
```

   (First check where each 7a IR snapshot actually lives — `ls tests/*/ir/` — and match
   the exact paths/filenames; the loop above assumes `tests/<Folder>/ir/<model>.gdsge.json`.)

3. Regenerate the HL1996 codegen snapshots (the whitelist literal changed in every
   generated file):

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); addpath(fullfile('tests','HeatonLucas1996','ir')); run(fullfile('tests','HeatonLucas1996','codegen','regen_snapshots.m'))"
```

   Review `git diff` on all regenerated files: the JSONs gain exactly `setupNames` (and
   nothing else); the codegen snapshots change only in the `GDSGE_OPTIONS_VALID` line.
   Unexplained churn = stop and debug.

- [ ] **Step 6: Write the two front-end gates**

Create `tests/CaoKS2016/parser/tFrontEndCaoKS2016.m`:

```matlab
classdef tFrontEndCaoKS2016 < matlab.unittest.TestCase
    % Phase 7b front-end gate: parse -> validate -> JSON round-trip ->
    % structural assertions, per the 7a gate pattern.
    methods (Test)
        function frontEndProducesValidIR(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            ir = gdsge.parser.parseFrontEnd( ...
                fileread(fullfile(modelDir, 'CaoKS2016.gmod')), 'CaoKS2016');
            tc.verifyTrue(gdsge.ir.validate(ir));
            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));

            tc.verifyEqual(ir.options.interpMethod, 'asg');
            tc.verifyEqual(ir.options.asgMaxLevel, 10);
            tc.verifyEqual(ir.options.asgThreshold, 1e-4);
            tc.verifyEqual(ir.options.asgMinLevel, 4);
            tc.verifyEqual(ir.options.tolEq, 1e-4);
            tc.verifyEqual(ir.options.simuResolve, 1);

            tc.verifyEqual(ir.states.names, {'K','X'});
            tc.verifyEqual(ir.shocks.count, 4);
            tc.verifyEqual(numel(ir.variables.policy), 8);
            tc.verifyEqual(numel(ir.variables.interp), 2);
            tc.verifyEqual(ir.variables.output, {'Kp','Xp','kp1','kp2'});
            tc.verifyEqual(numel(ir.model.equations), 8);
            tc.verifyEmpty(ir.variables.tensor);   % the var_tensor block is commented out

            % model_init: square 2x2 system
            tc.verifyTrue(isfield(ir, 'modelInit'));
            tc.verifyEqual(numel(ir.modelInit.variables.policyInit), 2);
            tc.verifyEqual(numel(ir.modelInit.equations), 2);

            % the named-transition pipe on a square parameter
            reds = ir.model.statements(cellfun(@(s) strcmp(s.type,'reduction'), ...
                ir.model.statements));
            refs = cellfun(@(s) s.transRef, reds, 'UniformOutput', false);
            tc.verifyTrue(ismember('shock_trans2', refs));
            tc.verifyTrue(ismember('shock_trans', refs));

            % setup-assigned names feed the runtime override whitelist
            tc.verifyTrue(all(ismember({'KMin','KMax','XPts','EZTrans'}, ir.setupNames)));
        end
    end
end
```

Create `tests/Bianchi2011_asg/parser/tFrontEndBianchi2011.m`:

```matlab
classdef tFrontEndBianchi2011 < matlab.unittest.TestCase
    methods (Test)
        function frontEndProducesValidIR(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            ir = gdsge.parser.parseFrontEnd( ...
                fileread(fullfile(modelDir, 'bianchi2011.gmod')), 'bianchi2011');
            tc.verifyTrue(gdsge.ir.validate(ir));
            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));

            tc.verifyEqual(ir.options.interpMethod, 'asg');
            tc.verifyEqual(ir.options.asgMaxLevel, 14);
            tc.verifyEqual(ir.options.asgThreshold, 1e-5);
            tc.verifyEqual(ir.options.tolEq, 1e-6);        % gmod sets none -> default

            tc.verifyEqual(ir.states.names, {'b'});
            tc.verifyEqual(ir.shocks.count, 16);
            tc.verifyEqual(numel(ir.variables.policy), 4);
            auxNames = cellfun(@(a) a.name, ir.variables.aux, 'UniformOutput', false);
            tc.verifyEqual(auxNames, {'c','lambda','bNext'});
            tc.verifyEqual(ir.variables.interp, {'lambda_interp'});
            tc.verifyEqual(ir.variables.output, {'bNext','pN','c'});
            tc.verifyEqual(numel(ir.model.equations), 4);

            % placeholder shock process (the driver overrides at runtime)
            tc.verifyEqual(ir.shocks.transitions.shock_trans, zeros(16));

            % degenerate init system: dummy policy vs literal-0 equation
            tc.verifyTrue(isfield(ir, 'modelInit'));
            tc.verifyEqual(numel(ir.modelInit.variables.policyInit), 1);
            tc.verifyEqual(ir.modelInit.variables.policyInit{1}.name, 'dummy');
            tc.verifyEqual(numel(ir.modelInit.variables.auxInit), 2);
            tc.verifyEqual(numel(ir.modelInit.equations), 1);

            % bMin/bMax (driver overrides) are setup-assigned names
            tc.verifyTrue(all(ismember({'bMin','bMax','bPts'}, ir.setupNames)));
        end
    end
end
```

- [ ] **Step 7: Run the gates — fix what actually breaks**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','CaoKS2016','parser')), runtests(fullfile('tests','Bianchi2011_asg','parser'))]; disp(table(results)); exit(any([results.Failed]))"
```

Likely-green areas (7a built them): `model_init`, `*_init` declarations, named interp
calls. Plausible breakage: the literal `0;` init equation (if `parseExpr`/`parseModel`
rejects a bare numeric statement, fix `parseModel`'s equation branch to accept any
expression — the AST `num` node exists); validator complaints about the placeholder
all-zeros `shock_trans` (if a normalization check exists, confirm what the validator
actually requires and relax ONLY for values, never shape). Use systematic-debugging; add
a focused unit test for any parser fix.

- [ ] **Step 8: Run the parser/ir suites + HL1996 + 7a front-end suites — verify green**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','parser'),'IncludingSubfolders',true), runtests(fullfile('tests','ir'),'IncludingSubfolders',true), runtests(fullfile('tests','codegen'),'IncludingSubfolders',true), runtests(fullfile('tests','HeatonLucas1996'),'IncludingSubfolders',true), runtests(fullfile('tests','Barro_et_al_2017'),'IncludingSubfolders',true), runtests(fullfile('tests','Mendoza2010','parser')), runtests(fullfile('tests','GLSW2020','parser'))]; disp(table(results)); exit(any([results.Failed]))"
```

Expected: PASS, including the HL1996 Slow end-to-end gate (the whitelist line changed in
generated files — the runtime behavior must not).

- [ ] **Step 9: Commit**

```bash
git add src/+gdsge tests docs/ir-schema.md
git commit -m "feat(parser,ir): ASG options, setupNames whitelist, square-param transition refs"
```

---

### Task 5: C++ backend — ASG interp templates, compile defines, driver unlock

**Files:**
- Create: `templates/cxx/interp_asg_construct.tpl.cpp`, `interp_asg_get.tpl.cpp`, `interp_asg_prepare_space.tpl.cpp`
- Modify: `src/+gdsge/+codegen/+cxx/emitInterp.m`, `emitCompile.m`, `src/+gdsge/+codegen/generateCxx.m`, `src/+gdsge/+runtime/solveProblems.m`
- Create: `tests/+gdsgefix/minimalIRAsg.m`
- Test: `tests/codegen/tEmitAsg.m` (cxx half)

- [ ] **Step 1: Copy the three template fragments** (verbatim from the old toolbox; the
  new pipeline substitutes the same holes)

```bash
cp base_package/gdsge/source/code_template/interp_asg_construct_template.cpp templates/cxx/interp_asg_construct.tpl.cpp
cp base_package/gdsge/source/code_template/interp_asg_get_template.cpp templates/cxx/interp_asg_get.tpl.cpp
cp base_package/gdsge/source/code_template/interp_asg_prepare_space_template.cpp templates/cxx/interp_asg_prepare_space.tpl.cpp
```

- [ ] **Step 2: Create the ASG IR fixture**

Create `tests/+gdsgefix/minimalIRAsg.m`:

```matlab
function ir = minimalIRAsg()
% MINIMALIRASG  minimalIR with interpMethod='asg' + the Asg option set.
%   Schema-valid; used by the emitter unit tests (no MEX compile).
ir = gdsgefix.minimalIR();
ir.options.interpMethod = 'asg';
ir.options.asgMinLevel = 2;
ir.options.asgMaxLevel = 4;
ir.options.asgThreshold = 1e-2;
ir.options.asgOutputMaxLevel = 4;
ir.options.asgOutputThreshold = 1e-2;
ir.options.simuResolve = 1;
ir.options.simuInterp = 0;
assert(gdsge.ir.validate(ir));
end
```

(Read `tests/+gdsgefix/minimalIR.m` first; if its options struct lacks `simuResolve`
fields or the validator wants more for asg, satisfy the validator — the `assert` is the
spec here.)

- [ ] **Step 3: Write the failing cxx tests**

Create `tests/codegen/tEmitAsg.m` (the iter/simulate methods arrive in Tasks 7-8; start
with the cxx half):

```matlab
classdef tEmitAsg < matlab.unittest.TestCase
    % ASG emitter unit tests on the minimal fixture: key contract lines must
    % appear in the generated text. The real artifacts are pinned by the
    % CaoKS2016/Bianchi end-to-end gates; these tests catch emitter wiring
    % regressions cheaply.
    methods (Test)
        function cxxInterpSectionsAreAsg(tc)
            ir = gdsgefix.minimalIRAsg();
            frag = gdsge.codegen.cxx.emitInterp(ir);
            tc.verifyTrue(contains(frag.getCode, 'GET_MX_ARRAY(GDSGE_ASG_HANDLE)'));
            tc.verifyTrue(contains(frag.getCode, 'AsgInterpArrayAdoubleEvaluator'));
            tc.verifyFalse(contains(frag.getCode, 'GDSGE_SPLINE_VEC'));
            % in-place numeric substitution of DIM/NVEC; LEVEL stays a define
            tc.verifyTrue(contains(frag.threadCode, 'GDSGE_ASG_CELL[(ASG_MAX_LEVEL+2)*1]'));
            tc.verifyTrue(contains(frag.threadCode, 'GDSGE_INTERP_RSLT_adouble[1]'));
            % per-interp evaluator lambda for the fixture's interp var
            tc.verifyTrue(contains(frag.threadCode, [ir.variables.interp{1} '_adouble']));
        end
        function compileGainsAsgDefines(tc)
            gdsge.codegen.ensureAsgMex();
            ir = gdsgefix.minimalIRAsg();
            txt = gdsge.codegen.cxx.emitCompile(ir, 'include');
            tc.verifyTrue(contains(txt, '-DUSE_ASG'));
            tc.verifyTrue(contains(txt, '-DASG_MAX_DIM='));
            tc.verifyTrue(contains(txt, '-DASG_MAX_NVEC='));
            tc.verifyTrue(contains(txt, '-DASG_MAX_LEVEL='));
        end
        function generateCxxAcceptsAsg(tc)
            gdsge.codegen.ensureAsgMex();
            ir = gdsgefix.minimalIRAsg();
            out = tempname; mkdir(out);
            files = gdsge.codegen.generateCxx(ir, out);
            cpp = fileread(files.cppFile);
            tc.verifyTrue(contains(cpp, 'GET_MX_ARRAY(GDSGE_ASG_HANDLE)'));
        end
    end
end
```

(The numeric literals in `cxxInterpSectionsAreAsg` — `*1]`, `_adouble[1]` — assume
minimalIR has exactly 1 state and 1 interp var; read `tests/+gdsgefix/minimalIR.m` first
and adjust the expected numbers to the fixture's actual counts.)

- [ ] **Step 4: Run it — verify it fails** (spline fragments / "Phase 7b" refusal):

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','codegen','tEmitAsg.m')); disp(table(results)); exit(any([results.Failed]))"
```

- [ ] **Step 5: Implement**

`src/+gdsge/+codegen/+cxx/emitInterp.m` — at the top, after the `numInterp == 0` early
return, branch:

```matlab
if strcmp(ir.options.interpMethod, 'asg')
    frag = emitInterpAsg(ir);
    return;
end
```

and append the local function (old parser parity, `gdsge_parser.m:1170-1241`; substring
substitution order matters — ADOUBLE before DOUBLE before VAR):

```matlab
function frag = emitInterpAsg(ir)
% ASG variant: one shared evaluator over the class handle. getCode (task
% scope) converts GDSGE_ASG_HANDLE; threadCode declares per-thread scratch
% (DIM/NVEC substituted in place; ASG_MAX_LEVEL stays a compile-time define)
% plus the GDSGE_INTERP_VEC lambdas and one named lambda per interp var. The
% evaluator sites are always the full state vector (the parser forbids
% position-arg interps under ASG via the old rule; both corpus models comply).
import gdsge.codegen.cxx.readTemplate
states = ir.states.names;
numInterp = numel(ir.interp);
getCode = readTemplate('interp_asg_construct.tpl.cpp');
threadCode = readTemplate('interp_asg_prepare_space.tpl.cpp');
threadCode = strrep(threadCode, 'ASG_MAX_DIM',  num2str(numel(states)));
threadCode = strrep(threadCode, 'ASG_MAX_NVEC', num2str(numInterp));
getTpl = readTemplate('interp_asg_get.tpl.cpp');
for j = 1:numInterp
    filled = strrep(getTpl, 'INTERP_NAME', ir.interp{j}.name);
    filled = strrep(filled, 'INTERP_IDX', num2str(j - 1));
    threadCode = [threadCode newline filled]; %#ok<AGROW>
end
adoubleArgs = strjoin(cellfun(@(v) ['adouble ' v], states, 'UniformOutput', false), ',');
doubleArgs  = strjoin(cellfun(@(v) ['double '  v], states, 'UniformOutput', false), ',');
threadCode = strrep(threadCode, 'ADOUBLE_VAR_NAME', adoubleArgs);
threadCode = strrep(threadCode, 'DOUBLE_VAR_NAME',  doubleArgs);
threadCode = strrep(threadCode, 'VAR_NAME', strjoin(states, ','));
frag = struct('getCode', getCode, 'threadCode', threadCode);
end
```

Note the prepare-space template substitutes `ASG_MAX_DIM`/`ASG_MAX_NVEC` **before** the
per-interp blocks are appended, so the appended get-blocks keep nothing to substitute but
names — same as the old parser. `ASG_MAX_LEVEL` must survive into the emitted C++
(compile-time define) — the strrep of `ASG_MAX_DIM` does not touch it because the
template spells them differently; the unit test pins this.

`src/+gdsge/+codegen/+cxx/emitCompile.m` — after the `extraDef`/HAS_INIT block add:

```matlab
if strcmp(ir.options.interpMethod, 'asg')
    [asgMaxDim, asgMaxNvec, asgMaxLevel] = asg.get_mex_constants();
    if numel(ir.states.names) > asgMaxDim
        error('gdsge:codegen:asgDim', ...
            'model has %d states but asg_mex was built with ASG_MAX_DIM=%d', ...
            numel(ir.states.names), asgMaxDim);
    end
    if ir.options.asgMaxLevel > asgMaxLevel
        error('gdsge:codegen:asgLevel', ...
            'AsgMaxLevel=%d exceeds asg_mex ASG_MAX_LEVEL=%d', ...
            ir.options.asgMaxLevel, asgMaxLevel);
    end
    extraDef = strtrim([extraDef ' -DUSE_ASG' ...
        ' -DASG_MAX_DIM='   num2str(asgMaxDim) ...
        ' -DASG_MAX_NVEC='  num2str(asgMaxNvec) ...
        ' -DASG_MAX_LEVEL=' num2str(asgMaxLevel)]);
end
```

Also update its header comment (drop "ASG … not yet wired").

`src/+gdsge/+codegen/generateCxx.m` — replace the refusal block (lines 10-22) with:

```matlab
if strcmp(ir.options.interpMethod, 'pchip')
    error('gdsge:codegen:unsupported', 'interpMethod pchip: Phase 7c');
end
if strcmp(ir.options.interpMethod, 'asg')
    gdsge.codegen.ensureAsgMex();   % emitCompile reads asg.get_mex_constants()
end
if ~isempty(ir.hooks.cxx)
    error('gdsge:codegen:unsupported', 'cxx hook blocks: Phase 7c');
end
if ~isempty(ir.variables.tensor)
    error('gdsge:codegen:unsupported', 'var_tensor: Phase 7c');
end
```

`src/+gdsge/+runtime/solveProblems.m` — in the contract block, after the
`GDSGE_SPLINE_VEC = cfg.splineVec;` line add:

```matlab
if isfield(cfg, 'asgHandle')
    GDSGE_ASG_HANDLE = cfg.asgHandle;         %#ok<NASGU>  ASG models: the MEX
end                                           % reads it via mexGetVariable
```

- [ ] **Step 6: Run the tests — verify PASS**, then the codegen + HL1996 fast suites
  (spline emitInterp untouched paths must stay byte-stable):

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','codegen','tEmitAsg.m')), runtests(fullfile('tests','codegen'),'IncludingSubfolders',false)]; disp(table(results)); exit(any([results.Failed]))"
```

- [ ] **Step 7: Commit**

```bash
git add templates/cxx src/+gdsge tests/+gdsgefix/minimalIRAsg.m tests/codegen/tEmitAsg.m
git commit -m "feat(codegen-cxx): ASG interp path via class handle + ASG compile defines"
```

---

### Task 6: `gdsge.runtime.solveProblemsAsg` — the ASG solve cascade

**Files:**
- Create: `src/+gdsge/+runtime/solveProblemsAsg.m`
- Test: `tests/runtime/tSolveProblemsAsg.m`

- [ ] **Step 1: Write the failing test** (mock MEX handle — no compilation; a fake `asg`
  stand-in struct provides `get_current_level`/`eval_vec` via a tiny local class)

Create `tests/runtime/tSolveProblemsAsg.m`:

```matlab
classdef tSolveProblemsAsg < matlab.unittest.TestCase
    % The ASG cascade, old iter_solve_problem_asg_template.m parity:
    %   solve -> warm from NEW sol interp + solve -> restore from OLD (no
    %   solve) -> randomized restarts. Verified with a mock mexFn whose
    %   "model" converges only from a good initial guess.
    methods (Test)
        function warmStartFromNewInterpRescues(tc)
            % mock problem: 3 problems, solution x*=2; mexFn "converges" a
            % column only if |x0 - 2| < 0.5, else leaves f huge.
            calls = {};
            function [sol,f,aux,eqval,opt] = mockMex(sol,lb,ub,data,skip,f,aux,eqval) %#ok<INUSL>
                calls{end+1} = struct('sol',sol,'skip',skip); %#ok<AGROW>
                for i = 1:size(sol,2)
                    if skip(i); continue; end
                    if abs(sol(1,i) - 2) < 0.5
                        sol(1,i) = 2; f(i) = 0;
                    else
                        f(i) = 1e10;
                    end
                end
                opt = [];
            end
            cfg = baseCfg();
            cfg.solInterpNew = fakeInterp(2.1);   % good warm start
            cfg.solInterpOld = fakeInterp(-50);   % bad (never reached)
            sol0 = zeros(1,3); lb = -100*ones(1,3); ub = 100*ones(1,3);
            [sol, f, ~, ~, ~, diag] = gdsge.runtime.solveProblemsAsg(@mockMex, ...
                sol0, lb, ub, zeros(1,3), 1e20*ones(1,3), zeros(1,3), 1e20*ones(1,3), cfg);
            tc.verifyEqual(sol, 2*ones(1,3));
            tc.verifyEqual(f, zeros(1,3));
            tc.verifyTrue(all(diag.solved));
            tc.verifyEqual(numel(calls), 2);      % initial call + warm-start call
            tc.verifyEqual(calls{2}.sol(1,:), 2.1*ones(1,3));   % warm values used
        end
        function emptyInterpsFallThroughToRandom(tc)
            function [sol,f,aux,eqval,opt] = mockMex(sol,lb,ub,data,skip,f,aux,eqval) %#ok<INUSL>
                for i = 1:size(sol,2)
                    if skip(i); continue; end
                    if abs(sol(1,i) - 2) < 0.5; sol(1,i) = 2; f(i) = 0; else; f(i) = 1e10; end
                end
                opt = [];
            end
            cfg = baseCfg();
            cfg.solInterpNew = [];                % init-problem shape of the cascade
            cfg.solInterpOld = [];
            cfg.maxMinorIter = 500;
            rng(0);                                % bounds [1.6,2.4]: random restart hits
            [sol, ~, ~, ~, ~, diag] = gdsge.runtime.solveProblemsAsg(@mockMex, ...
                zeros(1,2), 1.6*ones(1,2), 2.4*ones(1,2), zeros(1,2), ...
                1e20*ones(1,2), zeros(1,2), 1e20*ones(1,2), cfg);
            tc.verifyEqual(sol, 2*ones(1,2));
            tc.verifyTrue(all(diag.solved));
        end
    end
end

function cfg = baseCfg()
cfg = struct('tolSol', 1e-8, 'tolFun', 1e-8, 'solMaxIter', 100, 'numThreads', 1, ...
    'debugEvalOnly', 0, 'useBroyden', 0, 'finiteDiffDelta', 1e-6, 'useBroydenNow', 0, ...
    'taskName', 1, 'asgHandle', uint64(0), 'evalArrayIdx', ones(1, 3), ...
    'evalGrids', zeros(1, 3), 'useAdaptiveBound', 0, 'adaptTight', [], ...
    'adaptInSol', [], 'maxMinorIter', 5, 'verboseRetry', false);
end

function h = fakeInterp(warmValue)
h = FakeAsgInterp(warmValue);
end
```

And create the tiny stand-in `tests/runtime/FakeAsgInterp.m`:

```matlab
classdef FakeAsgInterp < handle
    % Duck-typed stand-in for the asg class in solveProblemsAsg unit tests.
    properties; warmValue; end
    methods
        function obj = FakeAsgInterp(v); obj.warmValue = v; end
        function l = get_current_level(~); l = 3; end
        function v = eval_vec(obj, arrayIdx, grids) %#ok<INUSD>
            v = obj.warmValue * ones(1, size(grids, 2));
        end
    end
end
```

- [ ] **Step 2: Run it — verify it fails** (`solveProblemsAsg` undefined).

- [ ] **Step 3: Implement** — create `src/+gdsge/+runtime/solveProblemsAsg.m`:

```matlab
function [sol, f, aux, eqVal, optInfo, diag] = solveProblemsAsg(mexFn, sol, lb, ub, data, f, aux, eqVal, cfg)
% SOLVEPROBLEMSASG  Solve one batch of ASG-proposed grid problems. Old
%   iter_solve_problem_asg_template.m parity:
%     1) solve everything;
%     2) warm unconverged from the CURRENT-iteration sol interp
%        (cfg.solInterpNew), tighten if cfg.useAdaptiveBound, re-solve;
%     3) restore remaining unconverged from the PREVIOUS-iteration sol interp
%        (cfg.solInterpOld), tighten — NO solve here (the randomize loop is
%        the next call);
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
GDSGE_SPLINE_VEC = [];                        %#ok<NASGU>
GDSGE_ASG_HANDLE = cfg.asgHandle;             %#ok<NASGU>
% ----------------------------------------------------------------------------

f(:) = 1e20;
skip = zeros(1, size(sol, 2));
[sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);

% 2) warm unconverged from the current-iteration sol interp, re-solve
if ~isempty(cfg.solInterpNew) && cfg.solInterpNew.get_current_level() >= 0
    needResolved = (f > cfg.tolSol) | isnan(f);
    sol0 = cfg.solInterpNew.eval_vec(cfg.evalArrayIdx, cfg.evalGrids);
    sol(:, needResolved) = sol0(:, needResolved);
    skip(:) = 0;
    skip(~needResolved) = 1;
    if cfg.useAdaptiveBound == 1
        [lb, ub] = cfg.adaptTight(sol, lb, ub);
    end
    [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);
end

% 3) restore remaining unconverged from the previous-iteration sol interp
if ~isempty(cfg.solInterpOld) && cfg.solInterpOld.get_current_level() >= 0
    needResolved = (f > cfg.tolSol) | isnan(f);
    sol0 = cfg.solInterpOld.eval_vec(cfg.evalArrayIdx, cfg.evalGrids);
    sol(:, needResolved) = sol0(:, needResolved);
    if cfg.useAdaptiveBound == 1
        [lb, ub] = cfg.adaptTight(sol, lb, ub);
    end
end

% 4) randomized restarts
minorIter = 0;
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
    minorIter = minorIter + 1;
    if cfg.verboseRetry
        fprintf('  asg resolve round %d: %d unconverged, worst residual %g\n', ...
            minorIter, nnz((f > cfg.tolSol) | isnan(f)), max(f(:)));
    end
end

diag = struct();
diag.needResolved = (f > cfg.tolSol) | isnan(f);
diag.solved = ~diag.needResolved;
diag.minorIters = minorIter;
end
```

- [ ] **Step 4: Run the test — verify PASS.**

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+runtime/solveProblemsAsg.m tests/runtime/tSolveProblemsAsg.m tests/runtime/FakeAsgInterp.m
git commit -m "feat(runtime): solveProblemsAsg — old ASG solve-cascade parity"
```

---

### Task 7: MATLAB codegen — `emitIterAsg`

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitBounds.m` (ASG rewrite mode), `emitDataPack.m` (`asgPack`)
- Create: `src/+gdsge/+codegen/+mat/emitIterAsg.m`
- Modify: `src/+gdsge/+codegen/generateMatlab.m`
- Test: additions to `tests/codegen/tEmitAsg.m`

- [ ] **Step 1: Write the failing tests** — add to `tests/codegen/tEmitAsg.m`:

```matlab
        function iterAsgHasRefinementLoop(tc)
            ir = gdsgefix.minimalIRAsg();
            txt = gdsge.codegen.mat.emitIterAsg(ir);
            tc.verifyTrue(contains(txt, 'GDSGE_ASG_HANDLE = GDSGE_ASG_INTERP.objectHandle;'));
            tc.verifyTrue(contains(txt, 'get_eval_grids(AsgThreshold)'));
            tc.verifyTrue(contains(txt, 'get_eval_grids(0.0)'));
            tc.verifyTrue(contains(txt, 'push_eval_results_at_valid(GDSGE_evalRslts, GDSGE_solved)'));
            tc.verifyTrue(contains(txt, 'push_eval_results_at_grids(GDSGE_evalArrayIdx(GDSGE_solved)'));
            tc.verifyTrue(contains(txt, '[GDSGE_Metric,GDSGE_MetricVec] = asg.compute_inf_metric(GDSGE_ASG_INTERP_NEW, GDSGE_ASG_INTERP);'));
            tc.verifyTrue(contains(txt, 'gdsge.runtime.solveProblemsAsg(@mex_'));
            tc.verifyTrue(contains(txt, 'IterRslt.asg_interp_struct = GDSGE_ASG_INTERP.convert_to_struct();'));
            tc.verifyTrue(contains(txt, 'IterRslt.sol_asg_interp_struct = GDSGE_SOL_ASG_INTERP.convert_to_struct();'));
            tc.verifyTrue(contains(txt, 'GDSGE_DATA(:) = '));
            tc.verifyTrue(contains(txt, 'GDSGE_evalArrayIdx;GDSGE_evalGrids]'));
            % ASG WarmUp surface
            tc.verifyTrue(contains(txt, 'asg.construct_from_struct(GDSGE_OPTIONS.WarmUp.asg_interp_struct)'));
            tc.verifyTrue(contains(txt, 'asg.construct_from_struct(GDSGE_OPTIONS.WarmUp.sol_asg_interp_struct)'));
            % spline machinery must be absent
            tc.verifyFalse(contains(txt, 'constructSplines'));
            tc.verifyFalse(contains(txt, 'GDSGE_TENSOR_shockIdx'));
            % setup-assigned names reach the runtime override whitelist
            tc.verifyTrue(contains(txt, 'GDSGE_OPTIONS_VALID'));
        end
```

- [ ] **Step 2: Run — verify the new test fails** (`emitIterAsg` undefined).

- [ ] **Step 3: Implement the two small fragment changes**

`src/+gdsge/+codegen/+mat/emitDataPack.m` — add before `frag.maxData`:

```matlab
frag.asgPack = sprintf('GDSGE_DATA(:) = [repmat([%s],1,GDSGE_NPROB);GDSGE_evalArrayIdx;GDSGE_evalGrids];', ...
    innerStr);
```

`src/+gdsge/+codegen/+mat/emitBounds.m` — signature becomes
`function frag = emitBounds(ir, mode)` with:

```matlab
if nargin < 2; mode = 'tensor'; end
```

and the reps line (line 17) becomes:

```matlab
if strcmp(mode, 'asg')
    % ASG: per-problem columns are the proposed eval grids, not tensors.
    reps = {};
    for k = 1:numel(ir.states.names)
        reps{end+1} = sprintf('GDSGE_evalGrids(%d,:)', k); %#ok<AGROW>
    end
    for k = 1:numel(ir.shocks.names)
        reps{end+1} = sprintf('%s(GDSGE_evalArrayIdx)', ir.shocks.names{k}); %#ok<AGROW>
    end
else
    reps = cellfun(@(n) ['GDSGE_TENSOR_' n '(:)'''], names, 'UniformOutput', false);
end
```

(`names` keeps its current `[states, shocks]` order — the reps lists above must follow
the same order. Both corpus models have constant bounds, so this rewrite is exercised
only by unit tests; it exists so state-dependent bounds fail loudly correct, not
silently wrong.)

- [ ] **Step 4: Implement `emitIterAsg`**

Create `src/+gdsge/+codegen/+mat/emitIterAsg.m`:

```matlab
function txt = emitIterAsg(ir)
% EMITITERASG  iter_<model>.m for interpMethod='asg'.
%   Old-template parity: iter_inf_horizon_asg_template.m with
%   asg_propose_grids_and_solve_template.m inlined (interp loop + output
%   construct) and iter_init_asg_template.m as the init segment. Differences
%   from the spline emitIter: no state-space tensors — each fixed-point
%   iteration rebuilds the interpolant by level-by-level adaptive refinement;
%   per-problem GDSGE_DATA rows are the proposed eval grids; the MEX
%   evaluates the PREVIOUS iteration's interpolant through GDSGE_ASG_HANDLE.
%   Deliberately not emitted (spec 7b §2 — no corpus driver exercises them):
%   the GDSGE_ASG_FIX_GRID branch, PRE_ITER/POST_ITER hooks. Progress goes
%   through gdsge.runtime.printIterProgress (7a precedent: console format is
%   a deliberate simplification; results are unaffected).
m           = ir.modelName;
bounds      = gdsge.codegen.mat.emitBounds(ir, 'asg');
pack        = gdsge.codegen.mat.emitDataPack(ir);
unpk        = gdsge.codegen.mat.emitSolUnpack(ir);
interpNames = ir.variables.interp;
stateList   = strjoin(ir.states.names, ',');
stateGridCell = ['{' stateList '}'];
interpSemi  = strjoin(interpNames, ';');
nAux = 0;
for i = 1:numel(ir.variables.aux); nAux = max(nAux, ir.variables.aux{i}.slot(2)); end
numPolicyTotal = sum(cellfun(@(p) p.length, ir.variables.policy));
outputs = ir.variables.output;
lenOf = struct();
for i = 1:numel(ir.variables.policy); lenOf.(ir.variables.policy{i}.name) = ir.variables.policy{i}.length; end
for i = 1:numel(ir.variables.aux);    lenOf.(ir.variables.aux{i}.name)    = ir.variables.aux{i}.length;    end
numOutputTotal = 0;
for i = 1:numel(outputs); numOutputTotal = numOutputTotal + lenOf.(outputs{i}); end
outputSemi = strjoin(outputs, ';');
main = struct('bounds', bounds, 'unpack', unpk.unpack, 'unpackAux', unpk.unpackAux, ...
    'maxDim', bounds.maxDim, 'nAux', max(nAux, 1));

w = gdsge.codegen.codeWriter();
w.add('%% Generated by gdsge.codegen.generateMatlab — do not edit.');
w.add('function [IterRslt,IterFlag] = iter_%s(GDSGE_OPTIONS)', m);
w.add('gdsge.runtime.ensurePath();');
w.blank();
w.addRaw(gdsge.codegen.mat.emitSetup(ir, 'iter'));
w.blank();
w.add('GDSGE_OPTIONS_VALID = %s;', gdsge.codegen.mat.optionsWhitelist(ir, {}));
w.add('if nargin < 1; GDSGE_OPTIONS = struct(); end');
w.add('gdsge.runtime.unpackOptions(GDSGE_OPTIONS, GDSGE_OPTIONS_VALID);');
w.blank();
w.add('assert(exist(''shock_num'',''var'')==1);');
for i = 1:numel(ir.shocks.names)
    w.add('assert(length(%s)==%d);', ir.shocks.names{i}, ir.shocks.count);
end
tn = fieldnames(ir.shocks.transitions);
for i = 1:numel(tn)
    w.add('assert(size(%s,1)==%d);', tn{i}, ir.shocks.count);
    w.add('assert(size(%s,2)==%d);', tn{i}, ir.shocks.count);
end
w.blank();

% ---- model-init segment (ASG refinement on the init problem) --------------
if isfield(ir, 'modelInit')
    addInitSegment(w, ir, pack, m);
    w.blank();
end

w.add('GDSGE_SOL_ASG_INTERP = asg(%s,%d,shock_num);', stateGridCell, numPolicyTotal);
w.add('GDSGE_Metric = 1;');
w.add('GDSGE_MetricVec = [];');
w.add('GDSGE_Iter = 0;');
w.blank();
w.add('if nargin>=1 && isfield(GDSGE_OPTIONS,''WarmUp'')');
w.add('    if isfield(GDSGE_OPTIONS.WarmUp,''asg_interp_struct'')');
w.add('        GDSGE_ASG_INTERP = asg.construct_from_struct(GDSGE_OPTIONS.WarmUp.asg_interp_struct);');
w.add('    end');
w.add('    if isfield(GDSGE_OPTIONS.WarmUp,''Iter'')');
w.add('        GDSGE_Iter = GDSGE_OPTIONS.WarmUp.Iter;');
w.add('    end');
w.add('    if isfield(GDSGE_OPTIONS.WarmUp,''sol_asg_interp_struct'')');
w.add('        GDSGE_SOL_ASG_INTERP = asg.construct_from_struct(GDSGE_OPTIONS.WarmUp.sol_asg_interp_struct);');
w.add('    end');
w.add('end');
w.blank();
w.add('stopFlag = false;');
w.add('GDSGE_timer = tic;');
w.add('while(~stopFlag)');
w.in();
w.add('GDSGE_Iter = GDSGE_Iter+1;');
w.blank();
w.add('GDSGE_ASG_HANDLE = GDSGE_ASG_INTERP.objectHandle;');
w.add('GDSGE_ASG_INTERP_NEW = asg(%s,%d,shock_num);', stateGridCell, numel(interpNames));
w.add('GDSGE_SOL_ASG_INTERP_NEW = asg(%s,%d,shock_num);', stateGridCell, numPolicyTotal);
w.add('GDSGE_ASG_STORE_evalArrayIdx = cell(0);');
w.add('GDSGE_ASG_STORE_evalGridsUnscaled = cell(0);');
w.add('GDSGE_ASG_STORE_output = cell(0);');
w.blank();
w.add('while GDSGE_ASG_INTERP_NEW.get_current_level < AsgMaxLevel');
w.in();
w.add('if GDSGE_ASG_INTERP_NEW.get_current_level<AsgMinLevel');
w.add('    [GDSGE_evalArrayIdx,GDSGE_evalGrids,GDSGE_evalGridsLength,GDSGE_evalGridsUnscaled] = GDSGE_ASG_INTERP_NEW.get_eval_grids(0.0); %%#ok<ASGLU>');
w.add('else');
w.add('    [GDSGE_evalArrayIdx,GDSGE_evalGrids,GDSGE_evalGridsLength,GDSGE_evalGridsUnscaled] = GDSGE_ASG_INTERP_NEW.get_eval_grids(AsgThreshold); %%#ok<ASGLU>');
w.add('end');
w.add('if isempty(GDSGE_evalGrids); break; end');
w.blank();
addProposeAndSolve(w, ir, main, pack, m, ...
    'GDSGE_SOL_ASG_INTERP_NEW', 'GDSGE_SOL_ASG_INTERP', 8);
w.blank();
for i = 1:numel(ir.interp)
    w.add('%s = %s;', ir.interp{i}.name, ir.interp{i}.updateExpr);
end
w.add('GDSGE_evalRslts = [%s];', interpSemi);
w.add('GDSGE_SOL_ASG_INTERP_NEW.push_eval_results_at_grids(GDSGE_evalArrayIdx(GDSGE_solved), GDSGE_evalGridsUnscaled(:, GDSGE_solved), GDSGE_SOL(:, GDSGE_solved), GDSGE_SOL_ASG_INTERP_NEW.get_current_level);');
w.add('GDSGE_ASG_INTERP_NEW.push_eval_results_at_valid(GDSGE_evalRslts, GDSGE_solved);');
w.add('GDSGE_ASG_STORE_evalArrayIdx = [GDSGE_ASG_STORE_evalArrayIdx,GDSGE_evalArrayIdx];');
w.add('GDSGE_ASG_STORE_evalGridsUnscaled = [GDSGE_ASG_STORE_evalGridsUnscaled,GDSGE_evalGridsUnscaled];');
if ~isempty(outputs)
    w.add('GDSGE_ASG_STORE_output = [GDSGE_ASG_STORE_output,[%s]];', outputSemi);
end
w.out();
w.add('end');
w.blank();
w.add('[GDSGE_Metric,GDSGE_MetricVec] = asg.compute_inf_metric(GDSGE_ASG_INTERP_NEW, GDSGE_ASG_INTERP);');
w.add('GDSGE_ASG_INTERP = GDSGE_ASG_INTERP_NEW;');
w.add('GDSGE_SOL_ASG_INTERP = GDSGE_SOL_ASG_INTERP_NEW;');
w.blank();
w.add('stopFlag = GDSGE_Metric<TolEq || GDSGE_Iter>=MaxIter;');
w.blank();
w.add('if gdsge.runtime.printIterProgress(GDSGE_Iter, GDSGE_Metric, max(GDSGE_F), nnz(~GDSGE_solved), toc(GDSGE_timer), PrintFreq, NoPrint, stopFlag)');
w.add('    GDSGE_timer = tic;');
w.add('end');
w.blank();
w.add('if ( mod(GDSGE_Iter,SaveFreq)==0 || stopFlag == true )');
w.in();
if ~isempty(outputs)
    w.add('if CONSTRUCT_OUTPUT==1');
    w.in();
    w.add('GDSGE_ASG_HANDLE = GDSGE_ASG_INTERP.objectHandle;');
    w.add('GDSGE_ASG_INTERP_OUTPUT = asg(%s,%d,shock_num);', stateGridCell, numOutputTotal);
    w.add('while GDSGE_ASG_INTERP_OUTPUT.get_current_level < AsgOutputMaxLevel');
    w.in();
    w.add('[GDSGE_TEMP_grids, GDSGE_TEMP_surplus, GDSGE_TEMP_levels, GDSGE_TEMP_unscaledGrids] = GDSGE_ASG_INTERP.get_grids_info_at_level(GDSGE_ASG_INTERP_OUTPUT.get_current_level+1); %%#ok<ASGLU>');
    w.add('GDSGE_evalArrayIdx = [];');
    w.add('for GDSGE_I_ARRAY=1:shock_num');
    w.add('    GDSGE_evalArrayIdx = [GDSGE_evalArrayIdx,GDSGE_I_ARRAY*ones(1,size(GDSGE_TEMP_grids{GDSGE_I_ARRAY},2))]; %%#ok<AGROW>');
    w.add('end');
    w.add('GDSGE_evalGrids = cat(2,GDSGE_TEMP_grids{:});');
    w.add('GDSGE_evalGridsUnscaled = cat(2,GDSGE_TEMP_unscaledGrids{:});');
    w.add('if isempty(GDSGE_evalGrids); break; end');
    addProposeAndSolve(w, ir, main, pack, m, ...
        'GDSGE_SOL_ASG_INTERP', 'GDSGE_SOL_ASG_INTERP', 16);
    w.add('GDSGE_evalRslts = [%s];', outputSemi);
    w.add('GDSGE_ASG_INTERP_OUTPUT.push_eval_results_at_grids(GDSGE_evalArrayIdx(GDSGE_solved), GDSGE_evalGridsUnscaled(:, GDSGE_solved), GDSGE_evalRslts(:, GDSGE_solved), GDSGE_ASG_INTERP_OUTPUT.get_current_level);');
    w.out();
    w.add('end');
    w.add('output_var_index=struct();');
    idx = 1;
    for i = 1:numel(outputs)
        w.add('output_var_index.%s=%d:%d;', outputs{i}, idx, idx + lenOf.(outputs{i}) - 1);
        idx = idx + lenOf.(outputs{i});
    end
    w.add('IterRslt.output_var_index = output_var_index;');
    w.add('IterRslt.asg_output_struct = GDSGE_ASG_INTERP_OUTPUT.convert_to_struct();');
    w.out();
    w.add('end');
end
w.add('IterRslt.Metric = GDSGE_Metric;');
w.add('IterRslt.MetricVec = GDSGE_MetricVec;');
w.add('IterRslt.Iter = GDSGE_Iter;');
w.add('IterRslt.shock_num = shock_num;');
w.add('IterRslt.shock_trans = shock_trans;');
for i = 1:numel(ir.shocks.names)
    w.add('IterRslt.var_shock.%s = %s;', ir.shocks.names{i}, ir.shocks.names{i});
end
for i = 1:numel(ir.states.names)
    w.add('IterRslt.var_state.%s = %s;', ir.states.names{i}, ir.states.names{i});
end
for i = 1:numel(ir.params)
    w.add('IterRslt.params.%s = %s;', ir.params{i}.name, ir.params{i}.name);
end
w.add('IterRslt.params.GDSGE_EMPTY = GDSGE_EMPTY;');
w.add('IterRslt.asg_interp_struct = GDSGE_ASG_INTERP.convert_to_struct();');
w.add('IterRslt.sol_asg_interp_struct = GDSGE_SOL_ASG_INTERP.convert_to_struct();');
w.add('IterRslt.var_others.GDSGE_EMPTY = GDSGE_EMPTY;');
w.add('if ~NoSave');
w.add('    if IterSaveAll');
w.add('        save([''IterRslt_%s_'' num2str(GDSGE_Iter) ''.mat'']);', m);
w.add('    else');
w.add('        save([''IterRslt_%s_'' num2str(GDSGE_Iter) ''.mat''],''IterRslt'');', m);
w.add('    end');
w.add('end');
w.out();
w.add('end');
w.out();
w.add('end');
w.blank();
w.add('IterFlag = 0;');
w.add('end');
w.blank();
w.add('function [GDSGE_LB,GDSGE_UB] = GDSGE_ADAPT_TIGHT(GDSGE_SOL,GDSGE_LB,GDSGE_UB)');
if isempty(bounds.adaptiveTight)
    w.add('%% no adaptive bounds in this model');
else
    w.addRaw(bounds.adaptiveTight);
end
w.add('end');
txt = w.str();
end

% =========================================================================
function addProposeAndSolve(w, ir, main, pack, m, newSolName, oldSolName, depth)
% One batch at the proposed grids (old asg_propose_grids_and_solve_template):
% bounds, pre-warm from the PREVIOUS-iteration sol interp, allocate, pack,
% solve through gdsge.runtime.solveProblemsAsg, unpack rows + GDSGE_solved.
ind = repmat(' ', 1, depth);
w.add('GDSGE_NPROB = size(GDSGE_evalGrids,2);');
w.add('GDSGE_SIZE = [1,GDSGE_NPROB];');
w.addRaw(gdsge.codegen.mat.indentBy(main.bounds.init, ind));
w.add('GDSGE_SOL0 = (GDSGE_LB + GDSGE_UB)/2;');
w.add('GDSGE_SOL = GDSGE_SOL0;');
w.add('if %s.get_current_level>=0', oldSolName);
w.add('    GDSGE_SOL0 = %s.eval_vec(GDSGE_evalArrayIdx,GDSGE_evalGrids);', oldSolName);
w.add('    GDSGE_SOL = GDSGE_SOL0;');
w.add('    if UseAdaptiveBound==1');
w.add('        [GDSGE_LB,GDSGE_UB] = GDSGE_ADAPT_TIGHT(GDSGE_SOL,GDSGE_LB,GDSGE_UB);');
w.add('    end');
w.add('end');
w.add('GDSGE_EQVAL = 1e20*ones(%d,GDSGE_NPROB);', main.maxDim);
w.add('GDSGE_F = 1e20*ones(1,GDSGE_NPROB);');
w.add('GDSGE_AUX = zeros(%d,GDSGE_NPROB);', main.nAux);
w.add('GDSGE_DATA = zeros(%d,GDSGE_NPROB);', pack.maxData);
w.add('%s', pack.asgPack);
for cfgLine = gdsge.codegen.mat.emitCfgCommon(); w.add('%s', cfgLine{1}); end
w.add('GDSGE_CFG.useBroydenNow = 0;');
w.add('GDSGE_CFG.taskName = MEX_TASK_INF_HORIZON;');
w.add('GDSGE_CFG.asgHandle = GDSGE_ASG_HANDLE;');
w.add('GDSGE_CFG.solInterpNew = %s;', newSolName);
w.add('GDSGE_CFG.solInterpOld = %s;', oldSolName);
w.add('GDSGE_CFG.evalArrayIdx = GDSGE_evalArrayIdx;');
w.add('GDSGE_CFG.evalGrids = GDSGE_evalGrids;');
w.add('GDSGE_CFG.useAdaptiveBound = UseAdaptiveBound;');
w.add('GDSGE_CFG.adaptTight = @GDSGE_ADAPT_TIGHT;');
w.add('GDSGE_CFG.maxMinorIter = MaxMinorIter;');
w.add('GDSGE_CFG.verboseRetry = ~NoPrint;');
w.add('if UseAdaptiveBoundInSol==1');
w.add('    GDSGE_CFG.adaptInSol = @GDSGE_ADAPT_TIGHT;');
w.add('else');
w.add('    GDSGE_CFG.adaptInSol = [];');
w.add('end');
w.add('[GDSGE_SOL,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL,GDSGE_OPT_INFO,GDSGE_DIAG] = gdsge.runtime.solveProblemsAsg(@mex_%s, GDSGE_SOL, GDSGE_LB, GDSGE_UB, GDSGE_DATA, GDSGE_F, GDSGE_AUX, GDSGE_EQVAL, GDSGE_CFG); %%#ok<ASGLU>', m);
w.add('GDSGE_solved = GDSGE_DIAG.solved;');
w.addRaw(gdsge.codegen.mat.indentBy(main.unpack, ind));
w.addRaw(gdsge.codegen.mat.indentBy(main.unpackAux, ind));
end

% =========================================================================
function addInitSegment(w, ir, pack, m)
% Old iter_init_asg_template.m parity: ASG refinement on the init problem
% seeds the var_interp interpolant. Gated ONLY by SkipModelInit. The init
% cascade has no warm/restore steps (solInterpNew/Old = []) and no adaptive
% tighten. asgHandle is unused by task_init (its INTERP_GET_CODE is empty).
iv = gdsge.codegen.initView(ir);
boundsInit = gdsge.codegen.mat.emitBounds(iv, 'asg');
unpkInit = gdsge.codegen.mat.emitSolUnpack(iv);
nAuxInit = 0;
for i = 1:numel(iv.variables.aux); nAuxInit = max(nAuxInit, iv.variables.aux{i}.slot(2)); end
stateList = strjoin(ir.states.names, ',');
interpNames = ir.variables.interp;
interpSemi = strjoin(interpNames, ';');
w.add('%%%% Model-init solve (ASG; old semantics: gated by SkipModelInit only)');
w.add('if ~SkipModelInit');
w.in();
w.add('GDSGE_ASG_INTERP = asg({%s},%d,shock_num);', stateList, numel(interpNames));
w.add('while GDSGE_ASG_INTERP.get_current_level < AsgMaxLevel');
w.in();
w.add('if GDSGE_ASG_INTERP.get_current_level<AsgMinLevel');
w.add('    [GDSGE_evalArrayIdx,GDSGE_evalGrids,GDSGE_evalGridsLength] = GDSGE_ASG_INTERP.get_eval_grids(0.0); %%#ok<ASGLU>');
w.add('else');
w.add('    [GDSGE_evalArrayIdx,GDSGE_evalGrids,GDSGE_evalGridsLength] = GDSGE_ASG_INTERP.get_eval_grids(AsgThreshold); %%#ok<ASGLU>');
w.add('end');
w.add('if isempty(GDSGE_evalGrids); break; end');
w.add('GDSGE_NPROB = size(GDSGE_evalGrids,2);');
w.add('GDSGE_SIZE = [1,GDSGE_NPROB];');
w.addRaw(gdsge.codegen.mat.indentBy(boundsInit.init, '        '));
w.add('GDSGE_SOL = (GDSGE_LB + GDSGE_UB)/2;');
w.add('GDSGE_EQVAL = 1e20*ones(%d,GDSGE_NPROB);', boundsInit.maxDim);
w.add('GDSGE_F = 1e20*ones(1,GDSGE_NPROB);');
w.add('GDSGE_AUX = zeros(%d,GDSGE_NPROB);', max(nAuxInit, 1));
w.add('GDSGE_DATA = zeros(%d,GDSGE_NPROB);', pack.maxData);
w.add('%s', pack.asgPack);
for cfgLine = gdsge.codegen.mat.emitCfgCommon(); w.add('%s', cfgLine{1}); end
w.add('GDSGE_CFG.useBroydenNow = 0;');
w.add('GDSGE_CFG.taskName = MEX_TASK_INIT;');
w.add('GDSGE_CFG.asgHandle = [];');
w.add('GDSGE_CFG.solInterpNew = [];');
w.add('GDSGE_CFG.solInterpOld = [];');
w.add('GDSGE_CFG.evalArrayIdx = GDSGE_evalArrayIdx;');
w.add('GDSGE_CFG.evalGrids = GDSGE_evalGrids;');
w.add('GDSGE_CFG.useAdaptiveBound = 0;');
w.add('GDSGE_CFG.adaptTight = [];');
w.add('GDSGE_CFG.adaptInSol = [];');
w.add('GDSGE_CFG.maxMinorIter = MaxMinorIter;');
w.add('GDSGE_CFG.verboseRetry = ~NoPrint;');
w.add('[GDSGE_SOL,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL,GDSGE_OPT_INFO,GDSGE_DIAG] = gdsge.runtime.solveProblemsAsg(@mex_%s, GDSGE_SOL, GDSGE_LB, GDSGE_UB, GDSGE_DATA, GDSGE_F, GDSGE_AUX, GDSGE_EQVAL, GDSGE_CFG); %%#ok<ASGLU>', m);
w.add('GDSGE_solved = GDSGE_DIAG.solved;');
w.addRaw(gdsge.codegen.mat.indentBy(unpkInit.unpack, '        '));
w.addRaw(gdsge.codegen.mat.indentBy(unpkInit.unpackAux, '        '));
for i = 1:numel(ir.interp)
    w.add('%s = %s;', ir.interp{i}.name, ir.interp{i}.initialExpr);
end
w.add('GDSGE_evalRslts = [%s];', interpSemi);
w.add('GDSGE_ASG_INTERP.push_eval_results_at_valid(GDSGE_evalRslts, GDSGE_solved);');
w.out();
w.add('end');
w.out();
w.add('end');
end
```

Before writing, read `src/+gdsge/+codegen/codeWriter.m` and `emitCfgCommon.m`: confirm
`w.add` auto-indents under `w.in()` and what indent string a given nesting depth needs
for the `indentBy` calls (the `depth` argument of `addProposeAndSolve`: 8 spaces at
interp-loop depth, 16 at the output-construct depth inside the save block — adjust to
the writer's actual indent width if it differs from 4 per level).

(`generateMatlab.m` stays unchanged in this task — its ASG branch needs `emitSimulateAsg`
and lands as part of Task 8.)

- [ ] **Step 5: Run the tests — verify the emitIterAsg assertions PASS**, plus the full
  codegen + 7a suites (the `emitBounds` signature change and `asgPack` must not disturb
  spline output):

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','codegen'),'IncludingSubfolders',true), runtests(fullfile('tests','HeatonLucas1996','codegen'))]; disp(table(results)); exit(any([results.Failed]))"
```

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+codegen tests/codegen/tEmitAsg.m
git commit -m "feat(codegen-mat): emitIterAsg — ASG refinement-loop iter with init + output construct"
```

---

### Task 8: MATLAB codegen — `emitSimulateAsg` + generateMatlab branch

**Files:**
- Create: `src/+gdsge/+codegen/+mat/emitSimulateAsg.m`
- Modify: `src/+gdsge/+codegen/generateMatlab.m`
- Test: additions to `tests/codegen/tEmitAsg.m`

- [ ] **Step 1: Write the failing tests** — add to `tests/codegen/tEmitAsg.m`:

```matlab
        function simulateAsgWarmStartsFromSolInterp(tc)
            ir = gdsgefix.minimalIRAsg();
            txt = gdsge.codegen.mat.emitSimulateAsg(ir);
            tc.verifyTrue(contains(txt, 'GDSGE_ASG_INTERP = asg.construct_from_struct(IterRslt.asg_interp_struct);'));
            tc.verifyTrue(contains(txt, 'GDSGE_ASG_HANDLE = GDSGE_ASG_INTERP.objectHandle;'));
            tc.verifyTrue(contains(txt, 'GDSGE_SOL_ASG_INTERP = asg.construct_from_struct(IterRslt.sol_asg_interp_struct);'));
            tc.verifyTrue(contains(txt, 'GDSGE_SOL = GDSGE_SOL_ASG_INTERP.eval_vec(SimuRslt.shock(:,GDSGE_t)'''));
            tc.verifyTrue(contains(txt, 'gdsge.runtime.solveProblems(@mex_'));
            tc.verifyTrue(contains(txt, 'GDSGE_CFG.asgHandle = GDSGE_ASG_HANDLE;'));
            tc.verifyTrue(contains(txt, 'gen_discrete_markov_rn'));
            tc.verifyFalse(contains(txt, 'myppual'));
            tc.verifyFalse(contains(txt, 'GDSGE_SPLINE_VEC = IterRslt.pp'));
        end
        function generateMatlabBranchesOnAsg(tc)
            ir = gdsgefix.minimalIRAsg();
            out = tempname; mkdir(out);
            files = gdsge.codegen.generateMatlab(ir, out);
            tc.verifyTrue(contains(fileread(files.iterFile), 'solveProblemsAsg'));
            tc.verifyTrue(contains(fileread(files.simulateFile), 'sol_asg_interp_struct'));
        end
```

- [ ] **Step 2: Run — verify they fail.**

- [ ] **Step 3: Implement `emitSimulateAsg`**

Create `src/+gdsge/+codegen/+mat/emitSimulateAsg.m`:

```matlab
function txt = emitSimulateAsg(ir)
% EMITSIMULATEASG  simulate_<model>.m for interpMethod='asg' (SIMU_RESOLVE;
%   old simulate_resolve_asg_template.m parity). Per period: warm-start the
%   solver from the solution interpolant (sol_asg_interp_struct), re-solve
%   through the MEX with the converged var-interp ASG (asg_interp_struct)
%   evaluated inside the equations via GDSGE_ASG_HANDLE. The retry cascade is
%   first-call + randomized restarts — exactly gdsge.runtime.solveProblems
%   with useNearestNeighbor=false. Old parity notes: bounds come from the
%   gmod inbound code only (no GDSGE_PROB widening — the ASG IterRslt has no
%   GDSGE_PROB); both corpus models have setup-constant bounds.
bounds = gdsge.codegen.mat.emitBounds(ir);
pack   = gdsge.codegen.mat.emitDataPack(ir);
unpk   = gdsge.codegen.mat.emitSolUnpack(ir);
simu   = gdsge.codegen.mat.emitResultSimu(ir);
m      = ir.modelName;
nAux   = 0;
for i = 1:numel(ir.variables.aux); nAux = max(nAux, ir.variables.aux{i}.slot(2)); end
stateNameCell = ['{''' strjoin(ir.states.names, ''',''') '''}'];
stateList = strjoin(ir.states.names, ',');
simuStateRows = strjoin(cellfun(@(s) ['SimuRslt.' s '(:,GDSGE_t)'''], ...
    ir.states.names, 'UniformOutput', false), ';');

w = gdsge.codegen.codeWriter();
w.add('%% Generated by gdsge.codegen.generateMatlab — do not edit.');
w.add('function SimuRslt = simulate_%s(IterRslt,GDSGE_OPTIONS)', m);
w.add('gdsge.runtime.ensurePath();');
w.blank();
w.addRaw(gdsge.codegen.mat.emitSetup(ir, 'simulate'));
w.blank();
w.add('GEN_SHOCK_START_PERIOD = 1;');
w.blank();
w.add('%% ---- pull the solved model from IterRslt (explicit assignments)');
tn = fieldnames(ir.shocks.transitions);
for i = 1:numel(tn)
    w.add('%s = IterRslt.shock_trans;', tn{i});
end
for i = 1:numel(ir.params)
    w.add('%s = IterRslt.params.%s;', ir.params{i}.name, ir.params{i}.name);
end
for i = 1:numel(ir.shocks.names)
    w.add('%s = IterRslt.var_shock.%s;', ir.shocks.names{i}, ir.shocks.names{i});
end
for i = 1:numel(ir.states.names)
    w.add('%s = IterRslt.var_state.%s;', ir.states.names{i}, ir.states.names{i});
end
w.add('shock_num = IterRslt.shock_num;');
w.blank();
w.add('GDSGE_OPTIONS_VALID = %s;', gdsge.codegen.mat.optionsWhitelist(ir, {'init','GEN_SHOCK_START_PERIOD'}));
w.add('if nargin < 2; GDSGE_OPTIONS = struct(); end');
w.add('gdsge.runtime.unpackOptions(GDSGE_OPTIONS, GDSGE_OPTIONS_VALID);');
w.blank();
w.add('%% ---- reconstruct the converged interpolants');
w.add('GDSGE_ASG_INTERP = asg.construct_from_struct(IterRslt.asg_interp_struct);');
w.add('GDSGE_ASG_HANDLE = GDSGE_ASG_INTERP.objectHandle;');
w.add('GDSGE_SOL_ASG_INTERP = asg.construct_from_struct(IterRslt.sol_asg_interp_struct);');
w.blank();
w.add('GDSGE_NPROB = num_samples;');
w.add('SimuRslt.shock = ones(num_samples,num_periods+1);');
w.addRaw(simu.prealloc);
w.blank();
w.addRaw(simu.init);
w.blank();
w.add('if nargin>1 && isfield(GDSGE_OPTIONS,''init'')');
w.addRaw(indent4(simu.initOverwrite));
w.add('end');
w.add('if any(SimuRslt.shock(:,1)>shock_num)');
w.add('    error(''initial shock exceeds shock_num'');');
w.add('end');
w.blank();
w.addRaw(bounds.init);
w.blank();
w.add('GDSGE_EQVAL = 1e20*ones(%d,GDSGE_NPROB);', bounds.maxDim);
w.add('GDSGE_F = 1e20*ones(1,GDSGE_NPROB);');
w.add('GDSGE_SOL = zeros(%d,GDSGE_NPROB);', bounds.maxDim);
w.add('GDSGE_AUX = zeros(%d,GDSGE_NPROB);', max(nAux, 1));
w.blank();
w.add('SimuRslt.shock(:,GEN_SHOCK_START_PERIOD:end) = gen_discrete_markov_rn(shock_trans,num_samples,length(GEN_SHOCK_START_PERIOD:num_periods+1),...');
w.add('    SimuRslt.shock(:,GEN_SHOCK_START_PERIOD));');
w.blank();
w.add('%s', pack.simuData0);
w.add('GDSGE_SHOCK_VAR_INDEX_BASE = ([0:num_samples-1]'')*shock_num;');
w.add('GDSGE_timer = tic;');
w.add('for GDSGE_t=1:num_periods');
w.in();
w.add('shock = SimuRslt.shock(:,GDSGE_t);');
for i = 1:numel(ir.states.names)
    w.add('%s=SimuRslt.%s(:,GDSGE_t);', ir.states.names{i}, ir.states.names{i});
end
w.blank();
w.add('%% solution-interp warm start (old line 78)');
w.add('GDSGE_SOL = GDSGE_SOL_ASG_INTERP.eval_vec(SimuRslt.shock(:,GDSGE_t)'',[%s]);', simuStateRows);
w.blank();
w.add('if UseAdaptiveBoundInSol==1');
w.add('    GDSGE_LB_OLD = GDSGE_LB; GDSGE_UB_OLD = GDSGE_UB;');
w.add('    [GDSGE_LB,GDSGE_UB] = GDSGE_ADAPT_TIGHT(GDSGE_SOL,GDSGE_LB,GDSGE_UB);');
w.add('    GDSGE_SOL_hitting_lower_bound = abs(GDSGE_SOL - GDSGE_LB_OLD) < 1e-8;');
w.add('    GDSGE_SOL_hitting_upper_bound = abs(GDSGE_SOL - GDSGE_UB_OLD) < 1e-8;');
w.add('    GDSGE_LB(~GDSGE_SOL_hitting_lower_bound) = GDSGE_LB_OLD(~GDSGE_SOL_hitting_lower_bound);');
w.add('    GDSGE_UB(~GDSGE_SOL_hitting_upper_bound) = GDSGE_UB_OLD(~GDSGE_SOL_hitting_upper_bound);');
w.add('end');
w.blank();
w.add('%s', pack.simuPack);
w.blank();
for cfgLine = gdsge.codegen.mat.emitCfgCommon(); w.add('%s', cfgLine{1}); end
w.add('GDSGE_CFG.useBroydenNow = 0;');
w.add('GDSGE_CFG.taskName = MEX_TASK_INF_HORIZON;');
w.add('GDSGE_CFG.splineVec = [];');
w.add('GDSGE_CFG.ppNames = {}; GDSGE_CFG.ppCell = {};');
w.add('GDSGE_CFG.asgHandle = GDSGE_ASG_HANDLE;');
w.add('GDSGE_CFG.maxMinorIter = MaxMinorIter;');
w.add('GDSGE_CFG.probSize = [1 GDSGE_NPROB];');
w.add('GDSGE_CFG.useNearestNeighbor = false;');
w.add('GDSGE_CFG.verboseRetry = false;');
w.add('if UseAdaptiveBoundInSol==1');
w.add('    GDSGE_CFG.adaptInSol = @GDSGE_ADAPT_TIGHT;');
w.add('else');
w.add('    GDSGE_CFG.adaptInSol = [];');
w.add('end');
w.add('[GDSGE_SOL,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL,GDSGE_OPT_INFO,GDSGE_DIAG] = gdsge.runtime.solveProblems(@mex_%s, GDSGE_SOL, GDSGE_LB, GDSGE_UB, GDSGE_DATA, GDSGE_F, GDSGE_AUX, GDSGE_EQVAL, GDSGE_CFG); %%#ok<ASGLU>', m);
w.add('if any(GDSGE_DIAG.needResolved)');
w.add('    warning(''gdsge:runtime:unconverged'', ''period %%d: %%s'', GDSGE_t, gdsge.runtime.reportUnconverged(GDSGE_DIAG.needResolved, GDSGE_F, [1 GDSGE_NPROB], %s, {%s}, 5));', ...
    stateNameCell, stateList);
w.add('end');
w.blank();
w.addRaw(indent4(unpk.unpack));
w.addRaw(indent4(unpk.unpackAux));
w.blank();
w.add('GDSGE_SHOCK_VAR_LINEAR_INDEX = SimuRslt.shock(:,GDSGE_t+1) + GDSGE_SHOCK_VAR_INDEX_BASE;');
w.addRaw(indent4(simu.assign));
w.blank();
w.add('if mod(GDSGE_t,SimuPrintFreq)==0');
w.add('    fprintf(''Periods: %%d\\n'', GDSGE_t);');
w.add('    SimuRsltNames = fieldnames(SimuRslt);');
w.add('    for GDSGE_field = 1:length(SimuRsltNames)');
w.add('        fprintf(''%%8s'', SimuRsltNames{GDSGE_field});');
w.add('    end');
w.add('    fprintf(''\\n'');');
w.add('    for GDSGE_field = 1:length(SimuRsltNames)');
w.add('        fprintf(''%%8.4g'', SimuRslt.(SimuRsltNames{GDSGE_field})(1,GDSGE_t));');
w.add('    end');
w.add('    fprintf(''\\n'');');
w.add('    fprintf(''elapsed:%%.1fs\\n'', toc(GDSGE_timer));');
w.add('    GDSGE_timer = tic;');
w.add('end');
w.add('if mod(GDSGE_t,SimuSaveFreq)==0');
w.add('    save([''SimuRslt_%s_'' num2str(GDSGE_t) ''.mat''], ''SimuRslt'');', m);
w.add('end');
w.out();
w.add('end');
w.add('end');
w.blank();
w.add('function [GDSGE_LB,GDSGE_UB] = GDSGE_ADAPT_TIGHT(GDSGE_SOL,GDSGE_LB,GDSGE_UB)');
if isempty(bounds.adaptiveTight)
    w.add('%% no adaptive bounds in this model');
else
    w.addRaw(bounds.adaptiveTight);
end
w.add('end');
txt = w.str();
end

function s = indent4(s)
s = gdsge.codegen.mat.indentBy(s, '    ');
end
```

Then `src/+gdsge/+codegen/generateMatlab.m` — replace the body with:

```matlab
files = struct();
files.iterFile = fullfile(outDir, ['iter_' ir.modelName '.m']);
if strcmp(ir.options.interpMethod, 'asg')
    gdsge.codegen.writeText(files.iterFile, gdsge.codegen.mat.emitIterAsg(ir));
    simTxt = gdsge.codegen.mat.emitSimulateAsg(ir);   % parser guarantees SIMU_RESOLVE
else
    gdsge.codegen.writeText(files.iterFile, gdsge.codegen.mat.emitIter(ir));
    if isfield(ir.options, 'simuInterp') && ir.options.simuInterp == 1
        simTxt = gdsge.codegen.mat.emitSimulateInterp(ir);
    else
        simTxt = gdsge.codegen.mat.emitSimulate(ir);
    end
end
files.simulateFile = fullfile(outDir, ['simulate_' ir.modelName '.m']);
gdsge.codegen.writeText(files.simulateFile, simTxt);
```

- [ ] **Step 4: Run the full codegen suite + 7a/HL1996 codegen — verify green:**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','codegen'),'IncludingSubfolders',true), runtests(fullfile('tests','HeatonLucas1996','codegen'))]; disp(table(results)); exit(any([results.Failed]))"
```

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen tests/codegen/tEmitAsg.m
git commit -m "feat(codegen-mat): emitSimulateAsg + ASG branch in generateMatlab"
```

---

### Task 9: CaoKS2016 end-to-end gate

**Files:**
- Create: `tests/CaoKS2016/codegen/tEndToEndCaoKS2016.m`

- [ ] **Step 1: Write the gate**

Create `tests/CaoKS2016/codegen/tEndToEndCaoKS2016.m`:

```matlab
classdef tEndToEndCaoKS2016 < matlab.unittest.TestCase
    % PHASE 7b GATE (CaoKS2016): public API end-to-end vs goldens, driven
    % exactly like the old test.m (iter with PrintFreq/SaveFreq overrides,
    % default SIMU_RESOLVE simulate), plus the driver's asg post-processing.
    % ASG comparisons are evaluation-based: tiny numeric drift can change
    % WHICH grids refine, so interpolants are compared by evaluating both on
    % the golden's grid set, never surplus-by-surplus. Slow.
    properties (Constant)
        RelTol = 1e-3;
        AbsTol = 1e-3;
    end
    methods (Test, TestTags = {'Slow'})
        function publicApiMatchesGolden(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'CaoKS2016.gmod'), work);

            ir = gdsge_codegen('CaoKS2016');

            mexFile = fullfile(work, ['mex_CaoKS2016.' mexext]);
            artifacts = {'iter_CaoKS2016.m', 'simulate_CaoKS2016.m', ...
                'mex_CaoKS2016.cpp', 'compile_CaoKS2016.m', ...
                'CaoKS2016.gdsge.json', 'mex_CaoKS2016.cache'};
            for i = 1:numel(artifacts)
                tc.assertTrue(exist(fullfile(work, artifacts{i}), 'file') == 2, ...
                    artifacts{i});
            end
            tc.assertTrue(exist(mexFile, 'file') == 3, 'MEX did not compile');
            decoded = gdsge.ir.decode(fileread(fullfile(work, 'CaoKS2016.gdsge.json')));
            tc.verifyTrue(gdsge.ir.isequalIR(decoded, ir), 'IR JSON does not round-trip');

            before = dir(mexFile);
            gdsge_codegen('CaoKS2016');
            after = dir(mexFile);
            tc.verifyEqual(after.datenum, before.datenum, ...
                'second gdsge_codegen run recompiled despite unchanged C++');

            % ---- iter (test.m parity)
            IterOptions.PrintFreq = 10;
            IterOptions.SaveFreq = 100;
            IterRslt = iter_CaoKS2016(IterOptions);

            golden = load(fullfile(modelDir, 'golden', 'IterRslt.mat'));
            G = golden.IterRslt;
            tc.verifyLessThan(IterRslt.Metric, 1e-4);
            tc.verifyLessThan(abs(IterRslt.Iter - G.Iter), 0.2*G.Iter + 20, ...
                sprintf('Iter=%d vs golden %d', IterRslt.Iter, G.Iter));
            r = gdsgetest.compareNumericClose(IterRslt.shock_trans, G.shock_trans, 1e-12, 1e-15);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_state, G.var_state, 1e-12, 1e-15);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));

            % ---- interpolants compared by evaluation on the golden grids
            tc.verifyTrue(isfield(IterRslt, 'asg_interp_struct'));
            tc.verifyTrue(isfield(IterRslt, 'sol_asg_interp_struct'));
            tc.verifyTrue(isfield(IterRslt, 'asg_output_struct'));
            assertInterpClose(tc, IterRslt.asg_interp_struct, G.asg_interp_struct, ...
                tc.RelTol, tc.AbsTol);
            assertInterpClose(tc, IterRslt.asg_output_struct, G.asg_output_struct, ...
                tc.RelTol, tc.AbsTol);
            tc.verifyEqual(IterRslt.output_var_index, G.output_var_index);

            % ---- the old driver's post-processing runs against new results
            A = asg.construct_from_struct(IterRslt.asg_interp_struct);
            [grids, surplus] = A.get_grids_info(); %#ok<ASGLU>
            tc.verifyEqual(numel(grids), 4);
            v = A.eval_vec(1*ones(1, size(grids{1}, 2)), grids{1});
            tc.verifyTrue(all(isfinite(v(:))));

            % ---- seeded reduced simulate (capture protocol)
            rng(0823);
            SimuRslt = simulate_CaoKS2016(IterRslt, ...
                struct('num_samples', 6, 'num_periods', 1000));
            gs = load(fullfile(modelDir, 'golden', 'SimuRslt.mat'));
            GS = gs.SimuRslt;
            tc.verifyEqual(SimuRslt.shock, GS.shock, ...
                'shock paths differ — check shock_trans bit-identity first');
            flds = {'K','X','kp1','kp2'};
            T0 = 100;
            for i = 1:numel(flds)
                a = SimuRslt.(flds{i}); b = GS.(flds{i});
                tc.verifyEqual(size(a), size(b), flds{i});
                r = gdsgetest.compareNumericClose(a(:,1:T0), b(:,1:T0), 1e-2, 1e-2);
                tc.verifyTrue(r.pass, sprintf('%s early path: %s', flds{i}, ...
                    strjoin(r.failures, newline)));
            end
        end
    end
end

function assertInterpClose(tc, newStruct, goldStruct, relTol, absTol)
% Evaluate both interpolants on the GOLDEN's stored grid set, per shock array.
A = asg.construct_from_struct(newStruct);
B = asg.construct_from_struct(goldStruct);
grids = B.get_grids_info();
for j = 1:numel(grids)
    g = grids{j};
    if isempty(g); continue; end
    idx = j*ones(1, size(g, 2));
    va = A.eval_vec(idx, g);
    vb = B.eval_vec(idx, g);
    r = gdsgetest.compareNumericClose(va, vb, relTol, absTol);
    tc.verifyTrue(r.pass, sprintf('interp eval mismatch, shock %d: %s', j, ...
        strjoin(r.failures, newline)));
end
end
```

- [ ] **Step 2: Run the gate** (budget: compile + the recorded capture wall-clock; this
  is the first time the whole ASG pipeline runs)

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','CaoKS2016','codegen','tEndToEndCaoKS2016.m')); disp(table(results)); exit(any([results.Failed]))"
```

Expected: PASS. When it does not pass first try, debug bottom-up with
superpowers:systematic-debugging — typical layers: C++ compile errors (template fill),
`GDSGE_DATA` shape mismatch (asgPack vs POP), interp evaluation inside equations
(handle/buffer wiring), refinement-loop logic (compare the generated `iter_CaoKS2016.m`
side-by-side with an old-toolbox-generated one: run `gdsge_codegen('CaoKS2016')` in a
scratch dir with the OLD source on the path and diff structurally). Never loosen the
gate's tolerances; fix causes.

- [ ] **Step 3: Run the full suite — everything stays green**

```
matlab -batch "cd('tests'); run_tests"
```

- [ ] **Step 4: Commit**

```bash
git add tests/CaoKS2016
git commit -m "test(gate): CaoKS2016 green end-to-end vs golden (ASG core loop)"
```

---

### Task 10: Bianchi2011_asg end-to-end gate (WarmUp / SkipModelInit / MaxIter / shock overrides)

**Files:**
- Create: `tests/Bianchi2011_asg/codegen/tEndToEndBianchi2011.m`

- [ ] **Step 1: Write the gate**

Create `tests/Bianchi2011_asg/codegen/tEndToEndBianchi2011.m`:

```matlab
classdef tEndToEndBianchi2011 < matlab.unittest.TestCase
    % PHASE 7b GATE (Bianchi2011_asg): the full two-stage driver — stage-1
    % iter with runtime shock overrides (options.shock_trans/yT/yN from
    % shock_process.mat) and MaxIter=50; stage-2 WarmUp re-solve with widened
    % bounds (options.b/bMin/bMax) and SkipModelInit=1 — then the seeded
    % reduced simulate and the driver's asg_output_struct post-processing.
    % Slow.
    properties (Constant)
        RelTol = 1e-3;
        AbsTol = 1e-3;
    end
    methods (Test, TestTags = {'Slow'})
        function fullDriverSequenceMatchesGolden(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'bianchi2011.gmod'), work);
            copyfile(fullfile(modelDir, 'shock_process.mat'), work);

            gdsge_codegen('bianchi2011');

            % ---- stage 1: shock overrides + MaxIter
            options = struct;
            shock_process = load('shock_process.mat');
            options.shock_trans = shock_process.shock_trans;
            options.yT = shock_process.yT;
            options.yN = shock_process.yN;
            options.MaxIter = 50;
            IterRslt1 = iter_bianchi2011(options);
            tc.verifyEqual(IterRslt1.Iter, 50);
            % the overrides reached the result (gmod values are placeholders)
            r = gdsgetest.compareNumericClose(IterRslt1.shock_trans, ...
                shock_process.shock_trans, 1e-12, 1e-15);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt1.var_shock.yT, ...
                shock_process.yT, 1e-12, 1e-15);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            g1 = load(fullfile(modelDir, 'golden', 'IterRslt1.mat'));
            assertInterpClose(tc, IterRslt1.asg_interp_struct, ...
                g1.IterRslt1.asg_interp_struct, tc.RelTol, tc.AbsTol);

            % ---- stage 2: WarmUp + widened bounds + SkipModelInit
            options.MaxIter = inf;
            options.WarmUp = IterRslt1;
            options.bMin = -1.1;
            options.bMax = 0.0;
            options.b = [options.bMin, options.bMax];
            options.SkipModelInit = 1;
            IterRslt = iter_bianchi2011(options);

            golden = load(fullfile(modelDir, 'golden', 'IterRslt.mat'));
            G = golden.IterRslt;
            tc.verifyLessThan(IterRslt.Metric, 1e-6);
            tc.verifyLessThan(abs(IterRslt.Iter - G.Iter), 0.2*G.Iter + 20, ...
                sprintf('Iter=%d vs golden %d', IterRslt.Iter, G.Iter));
            tc.verifyEqual(IterRslt.var_state.b, [-1.1, 0]);
            assertInterpClose(tc, IterRslt.asg_interp_struct, G.asg_interp_struct, ...
                tc.RelTol, tc.AbsTol);
            assertInterpClose(tc, IterRslt.asg_output_struct, G.asg_output_struct, ...
                tc.RelTol, tc.AbsTol);
            tc.verifyEqual(IterRslt.output_var_index, G.output_var_index);

            % ---- seeded reduced simulate
            rng(0823);
            SimuRslt = simulate_bianchi2011(IterRslt, ...
                struct('num_samples', 6, 'num_periods', 1000));
            gs = load(fullfile(modelDir, 'golden', 'SimuRslt.mat'));
            GS = gs.SimuRslt;
            tc.verifyEqual(SimuRslt.shock, GS.shock, ...
                'shock paths differ — check overridden shock_trans bit-identity');
            flds = {'b','c','pN'};
            T0 = 100;
            for i = 1:numel(flds)
                a = SimuRslt.(flds{i}); b = GS.(flds{i});
                tc.verifyEqual(size(a), size(b), flds{i});
                r = gdsgetest.compareNumericClose(a(:,1:T0), b(:,1:T0), 1e-2, 1e-2);
                tc.verifyTrue(r.pass, sprintf('%s early path: %s', flds{i}, ...
                    strjoin(r.failures, newline)));
            end

            % ---- the old driver's asg_output_struct post-processing
            A = asg.construct_from_struct(IterRslt.asg_output_struct);
            grids = A.get_grids_info();
            tc.verifyEqual(numel(grids), 16);
            for j = [1 4 16]
                g = grids{j};
                if isempty(g); continue; end
                lenGrid = size(g, 2);
                bNext_fval = A.eval(j*ones(1, lenGrid), 1*ones(1, lenGrid), g);
                tc.verifyTrue(all(isfinite(bNext_fval(:))));
            end
        end
    end
end

function assertInterpClose(tc, newStruct, goldStruct, relTol, absTol)
A = asg.construct_from_struct(newStruct);
B = asg.construct_from_struct(goldStruct);
grids = B.get_grids_info();
for j = 1:numel(grids)
    g = grids{j};
    if isempty(g); continue; end
    idx = j*ones(1, size(g, 2));
    va = A.eval_vec(idx, g);
    vb = B.eval_vec(idx, g);
    r = gdsgetest.compareNumericClose(va, vb, relTol, absTol);
    tc.verifyTrue(r.pass, sprintf('interp eval mismatch, shock %d: %s', j, ...
        strjoin(r.failures, newline)));
end
end
```

(Check `src/asg.m` for `eval`'s exact argument order — the old driver calls
`GDSGE_ASG_INTERP.eval(j*ones(1,lenGrid), bNext_idx*ones(1,lenGrid), grid)`, i.e.
(arrayIdx, vecIdx, sites) — and align the post-processing block to it.)

- [ ] **Step 2: Run the gate**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','Bianchi2011_asg','codegen','tEndToEndBianchi2011.m')); disp(table(results)); exit(any([results.Failed]))"
```

Expected: PASS. The new surface here relative to Task 9 is exactly: whitelist acceptance
of `bMin`/`bMax` (setup names), runtime shock overrides reaching the packed data, the
ASG WarmUp reconstruct + `Iter` resume, `SkipModelInit=1` skipping the init segment, and
`MaxIter` staging. If stage-1 errors on an option name, the whitelist is missing a name —
check `ir.setupNames` for the Bianchi IR first. Debug with systematic-debugging; never
loosen tolerances.

- [ ] **Step 3: Run the full suite**

```
matlab -batch "cd('tests'); run_tests"
```

Expected: exit 0, all suites green (the 7a/HL1996 gates prove no regression).

- [ ] **Step 4: Commit**

```bash
git add tests/Bianchi2011_asg
git commit -m "test(gate): Bianchi2011_asg green end-to-end vs golden (WarmUp/SkipModelInit/shock overrides)"
```

---

### Task 11: Close out the phase

**Files:**
- Modify: `PROGRESS.md`

- [ ] **Step 1: Full suite, one more time, fresh process**

```
matlab -batch "cd('tests'); run_tests"
```

Expected: exit 0. Record the total test count.

- [ ] **Step 2: Update PROGRESS.md** — flip Phase 7b to ☑ with a sub-list mirroring the
  7a entry (goldens with recorded Iter/Metric numbers and wall-clocks, vendored asg
  runtime + ensureAsgMex, front-end delta, emitters, runtime cascade, both gates), set
  Phase 7c to NEXT, and add a dated changelog entry summarizing: both models green
  through the public API vs fresh goldens; the asg public class vendored with cache-gated
  auto-compile (deliberate improvement over ship-prebuilt); setupNames whitelist widening;
  square-param transition refs; the ASG IterRslt field set frozen; ASG console print via
  printIterProgress (deliberate simplification, results unaffected); golden sizes.

- [ ] **Step 3: Commit**

```bash
git add PROGRESS.md
git commit -m "docs(progress): Phase 7b complete — ASG widening green end-to-end"
```

- [ ] **Step 4: Hand off** — implementation complete on branch `phase7b-asg-widening`;
  use superpowers:finishing-a-development-branch to decide merge/PR/cleanup with the user.
