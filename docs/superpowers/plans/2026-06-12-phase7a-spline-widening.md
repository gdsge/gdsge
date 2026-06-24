# Phase 7a — Spline-Path Widening: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Barro_et_al_2017 (`safe_assets`), Mendoza2010 (`mendoza2010`) and GLSW2020 (`GLSW_interp`) run end-to-end through the public API and match freshly captured goldens.

**Architecture:** Batch-capture three goldens from the old toolbox first; settle the shared parser/IR front-end (gmod `setup` replay, `model_init` + `*_init` declarations, var_simu→output→aux promotion, transition `primed` flag); then go vertical per model — Barro (near-free), Mendoza (init task + SIMU_INTERP + 2-D states), GLSW (simulate overrides) — each closed by an end-to-end golden gate.

**Tech Stack:** MATLAB R2025b (`matlab -batch`, native `matlab.unittest`), existing `gdsge.parser` / `gdsge.ir` / `gdsge.codegen` / `gdsge.runtime` packages, MSVC for MEX.

**Spec:** `docs/superpowers/specs/2026-06-12-phase7a-spline-widening-design.md`

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
  else adds ONLY `src/`.
- If a golden gate fails: use superpowers:systematic-debugging. **Never loosen tolerances.**

## Verified facts this plan builds on (do not re-derive)

- `splitBlocks` already recognizes `model_init` (`src/+gdsge/+parser/splitBlocks.m:11`);
  `splitStatements` splits on newline at depth 0, so Barro's semicolon-less
  `var_policy shr_x1 Rf omega1n[2]` line parses as its own statement.
- The IR schema already defines `modelInit` (optional) and options fields
  `simuPrintFreq`/`simuSaveFreq` (`src/+gdsge/+ir/schema.m:46-47,71-75`); `resolveOptions`
  does **not** yet extract the two freqs.
- All four reduction kinds (EXPECT/MIN/MAX/PROD) already lower in
  `src/+gdsge/+codegen/+cxx/emitModelBody.m:36-54`.
- `templates/cxx/mex.tpl.cpp` already carries the `#ifdef HAS_INIT` dispatch and the
  `TASK_INIT_CODE` placeholder; `templates/cxx/compile.tpl.m` carries `EXTRA_DEF`.
  `MODEL_CODE` sits **inside** the task function (`templates/cxx/task.tpl.cpp:80`), so the
  init task's `MODEL_1` lambdas cannot collide with the main task's.
- `gdsge.codegen.dataLayout` does not depend on policy/aux — the init task pops the
  **identical** GDSGE_DATA layout as the main task.
- HL1996 sets no `INTERP_ORDER` → order 4 is already the exercised default; GLSW's
  `INTERP_ORDER=4` needs no new spline work.
- The old generator replays the gmod declaration-region code verbatim in BOTH generated
  files (`setParamsCode`/`code3`); the new `emitSetup` instead emits numeric params — which
  breaks Barro (`Re_n` in bounds) and GLSW (`Ngrid/a1_lb/a1_ub` in grid text). Task 4
  switches to replay.
- The old parser promotes `var_simu` names into `var_output` (gdsge_parser.m:619-635) and
  missing output names into `var_aux` (lines 643-681), printing
  `The following var_simulate are added to var_output: ...` /
  `The following var_output are added to var_aux: ...`. Mendoza/Barro depend on this.
- Old transition lowering: primed RHS (`cTilde' = cTildeNext'`) indexes with
  `GDSGE_SHOCK_VAR_LINEAR_INDEX`; **unprimed** RHS (`k' = kNext`) gets `(:)`. The new
  `emitResultSimu.m:32` always emits the linear index — wrong for Mendoza/GLSW. Task 3/4
  adds the `primed` flag.
- Old `iter_template.m` gates the init solve **only** on `SkipModelInit`; it runs even
  under WarmUp (the warm-up interpolants then supersede it). Replicate exactly.
- Old `SIMU_INTERP` simulate: `GDSGE_PP = IterRslt.output_interp;` + per-period
  `myppual_mex` eval at `[shock; states]`, output names assigned from
  `output_var_index` rows; state clamp via `EnforceSimuStateInbound`.
  (`gdsge_parser.m:1684-1737`, `code_template/simulate_interp_template.m`.)
- `solveProblems` cfg contract: `src/+gdsge/+runtime/solveProblems.m:11-16`; with
  `cfg.splineVec=[] , ppNames={}, ppCell={}` no spline state is needed — `task_init` has
  empty `INTERP_GET_CODE` and never reads `GDSGE_SPLINE_VEC`.

## File structure

| File | Action | Responsibility |
|---|---|---|
| `tests/golden/capture_safe_assets.m` | Create | Old-toolbox golden capture, Barro |
| `tests/golden/capture_mendoza2010.m` | Create | Old-toolbox golden capture, Mendoza (incl. fine-grid re-solve) |
| `tests/golden/capture_GLSW_interp.m` | Create | Old-toolbox golden capture, GLSW (incl. init-override batch) |
| `tests/Barro_et_al_2017/safe_assets.gmod` | Create (copy) | Model source under test |
| `tests/Mendoza2010/mendoza2010.gmod` | Create (copy) | Model source under test |
| `tests/GLSW2020/GLSW_interp.gmod` | Create (copy) | Model source under test |
| `tests/<Model>/tGolden<Model>.m` (×3) | Create | Golden-integrity tests |
| `tests/<Model>/parser/tFrontEnd<Model>.m` (×3) | Create | Front-end gates: parse → validate → round-trip → structure |
| `tests/<Model>/codegen/tEndToEnd<Model>.m` (×3) | Create | End-to-end golden gates (Slow) |
| `src/+gdsge/+ir/schema.m` | Modify | `setup` section; `transItem.primed`; init-bound item with text name |
| `src/+gdsge/+parser/parseDeclarations.m` | Modify | Expose `setupText`; init variable tables + `boundsInit` |
| `src/+gdsge/+parser/parseVarDecls.m` | Modify | `var_policy_init`/`var_aux_init`/`inbound_init` parsing |
| `src/+gdsge/+parser/parseSimulate.m` | Modify | Transition `primed` flag |
| `src/+gdsge/+parser/resolveOptions.m` | Modify | Extract `SimuPrintFreq`/`SimuSaveFreq` |
| `src/+gdsge/+parser/resolveOutputs.m` | Create | var_simu→output→aux promotion (old parity) |
| `src/+gdsge/+parser/parseFrontEnd.m` | Modify | model_init block; promotion call |
| `src/+gdsge/+parser/assemblePartialIR.m` | Modify | `ir.setup`, `ir.modelInit`, init bounds gate |
| `src/+gdsge/+parser/analyzeModel.m` | Modify | model_init semantic pass |
| `src/+gdsge/+codegen/initView.m` | Create | modelInit → main-shaped IR view for emitter reuse |
| `src/+gdsge/+codegen/generateCxx.m` | Modify | Fill `TASK_INIT_CODE`; drop modelInit refusal |
| `src/+gdsge/+codegen/generateMatlab.m` | Modify | SIMU_INTERP branch |
| `src/+gdsge/+codegen/+cxx/emitTask.m` | Modify | `taskName` parameter |
| `src/+gdsge/+codegen/+cxx/emitCompile.m` | Modify | `max(main,init)` MAXDIM; `-DHAS_INIT` |
| `src/+gdsge/+codegen/+mat/emitSetup.m` | Modify | Replay `ir.setup` verbatim (replaces numeric emission) |
| `src/+gdsge/+codegen/+mat/emitIter.m` | Modify | Insert init segment |
| `src/+gdsge/+codegen/+mat/emitIterInit.m` | Create | The `if ~SkipModelInit` init-solve segment |
| `src/+gdsge/+codegen/+mat/emitSimulateInterp.m` | Create | SIMU_INTERP simulate variant |
| `src/+gdsge/+codegen/+mat/emitResultSimu.m` | Modify | primed-aware transitions |
| `tests/+gdsgefix/minimalIR.m` | Modify | `setup` + transition `primed` |
| `tests/+gdsgefix/minimalIRWithInit.m` | Create | minimal IR with a modelInit section |
| `tests/HeatonLucas1996/ir/buildHL1996IR.m` | Modify | `setup` text + transition primed |
| `tests/HeatonLucas1996/ir/HL1996.gdsge.json` | Regenerate | Reference JSON golden |
| `tests/HeatonLucas1996/codegen/golden/*` | Regenerate | Codegen snapshots (setup replay) |
| `docs/ir-schema.md` | Regenerate | gendoc no-drift |
| `tests/parser/*`, `tests/codegen/*` | Modify/Create | Unit tests per task below |
| `PROGRESS.md` | Modify | Check off 7a |

---

### Task 1: Branch setup

**Files:** none (git only)

- [ ] **Step 1: Create the phase branch**

```bash
git checkout -b phase7a-spline-widening
```

Expected: `Switched to a new branch 'phase7a-spline-widening'`. (`git status` must be clean first; if not, stop and report.)

---

### Task 2: Golden capture — all three models

**Files:**
- Create: `tests/Barro_et_al_2017/safe_assets.gmod`, `tests/Mendoza2010/mendoza2010.gmod`, `tests/GLSW2020/GLSW_interp.gmod` (verbatim copies)
- Create: `tests/golden/capture_safe_assets.m`, `tests/golden/capture_mendoza2010.m`, `tests/golden/capture_GLSW_interp.m`
- Create: `tests/Barro_et_al_2017/tGoldenSafeAssets.m`, `tests/Mendoza2010/tGoldenMendoza2010.m`, `tests/GLSW2020/tGoldenGLSW.m`

- [ ] **Step 1: Copy the three gmod files into the test tree**

```bash
mkdir -p tests/Barro_et_al_2017/golden tests/Mendoza2010/golden tests/GLSW2020/golden
cp base_package/gdsge/tests/Barro_et_al_2017/safe_assets.gmod tests/Barro_et_al_2017/
cp base_package/gdsge/tests/Mendoza2010/mendoza2010.gmod tests/Mendoza2010/
cp base_package/gdsge/tests/GLSW2020/GLSW_interp.gmod tests/GLSW2020/
```

- [ ] **Step 2: Write the capture scripts** (pattern: `tests/golden/capture_HL1996.m` — old source on path, temp cwd, full-convergence iter exactly as the model's `test.m`, reduced seeded simulate)

Create `tests/golden/capture_safe_assets.m`:

```matlab
function capture_safe_assets()
% Golden capture for Barro_et_al_2017 from the OLD toolbox. Protocol per the
% Phase 7a spec: full-convergence iter (test.m parity: no options beyond
% SaveFreq), reduced seeded simulate (test.m never simulates; the generated
% simulate_safe_assets is public surface, so we pin it too).
here     = fileparts(mfilename('fullpath'));          % tests/golden
repoRoot = fileparts(fileparts(here));
oldSrc   = fullfile(repoRoot, 'base_package', 'gdsge', 'source');
modelDir = fullfile(repoRoot, 'tests', 'Barro_et_al_2017');
goldenDir = fullfile(modelDir, 'golden');

work = tempname; mkdir(work);
copyfile(fullfile(modelDir, 'safe_assets.gmod'), work);

oldPath = path; restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(oldSrc);
oldCd = pwd; cdRestore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(work);

t0 = tic;
gdsge_codegen('safe_assets');
options = struct; options.SaveFreq = inf;
IterRslt = iter_safe_assets(options);
fprintf('iter wall-clock: %.1fs\n', toc(t0));

t0 = tic;
rng(0823);
SimuRslt = simulate_safe_assets(IterRslt, struct('num_samples', 6, 'num_periods', 1000));
fprintf('simulate wall-clock: %.1fs\n', toc(t0));

save(fullfile(goldenDir, 'IterRslt.mat'), 'IterRslt', '-v7');
save(fullfile(goldenDir, 'SimuRslt.mat'), 'SimuRslt', '-v7');
fprintf('Golden captured: Iter=%d Metric=%g\n', IterRslt.Iter, IterRslt.Metric);
end
```

Create `tests/golden/capture_mendoza2010.m`:

```matlab
function capture_mendoza2010()
% Golden capture for Mendoza2010 from the OLD toolbox. test.m parity for the
% iter chain: no-arg solve on default grids, then the 80x80 fine-grid WarmUp
% re-solve (that re-solve is itself a feature under test). Reduced seeded
% simulate (test.m runs 24x50000; golden uses 6x1000 per the spec bound).
here     = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(here));
oldSrc   = fullfile(repoRoot, 'base_package', 'gdsge', 'source');
modelDir = fullfile(repoRoot, 'tests', 'Mendoza2010');
goldenDir = fullfile(modelDir, 'golden');

work = tempname; mkdir(work);
copyfile(fullfile(modelDir, 'mendoza2010.gmod'), work);

oldPath = path; restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(oldSrc);
oldCd = pwd; cdRestore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(work);

t0 = tic;
gdsge_codegen('mendoza2010');
IterRslt1 = iter_mendoza2010;
fprintf('coarse iter wall-clock: %.1fs (Iter=%d)\n', toc(t0), IterRslt1.Iter);

t0 = tic;
options = struct;
options.cTilde = linspace(IterRslt1.var_state.cTilde(1), IterRslt1.var_state.cTilde(end), 80);
options.k      = linspace(IterRslt1.var_state.k(1),      IterRslt1.var_state.k(end),      80);
options.WarmUp = IterRslt1;
IterRslt = iter_mendoza2010(options);
fprintf('fine-grid iter wall-clock: %.1fs (Iter=%d)\n', toc(t0), IterRslt.Iter);

t0 = tic;
rng(0823);
SimuRslt = simulate_mendoza2010(IterRslt, struct('num_samples', 6, 'num_periods', 1000));
fprintf('simulate wall-clock: %.1fs\n', toc(t0));

save(fullfile(goldenDir, 'IterRslt1.mat'), 'IterRslt1', '-v7');
save(fullfile(goldenDir, 'IterRslt.mat'),  'IterRslt',  '-v7');
save(fullfile(goldenDir, 'SimuRslt.mat'),  'SimuRslt',  '-v7');
fprintf('Golden captured: Iter=%d Metric=%g\n', IterRslt.Iter, IterRslt.Metric);
end
```

Create `tests/golden/capture_GLSW_interp.m`:

```matlab
function capture_GLSW_interp()
% Golden capture for GLSW2020 from the OLD toolbox. test.m parity for iter
% (no options; gmod sets TolEq=1e-8, SIMU_INTERP=1). Two simulate goldens:
% the reduced seeded baseline and a small deterministic init-override batch
% (pins the options.init/num_samples/num_periods feature the IRF uses).
here     = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(here));
oldSrc   = fullfile(repoRoot, 'base_package', 'gdsge', 'source');
modelDir = fullfile(repoRoot, 'tests', 'GLSW2020');
goldenDir = fullfile(modelDir, 'golden');

work = tempname; mkdir(work);
copyfile(fullfile(modelDir, 'GLSW_interp.gmod'), work);

oldPath = path; restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(oldSrc);
oldCd = pwd; cdRestore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(work);

t0 = tic;
gdsge_codegen('GLSW_interp');
IterRslt = iter_GLSW_interp;
fprintf('iter wall-clock: %.1fs (Iter=%d)\n', toc(t0), IterRslt.Iter);

t0 = tic;
rng(0823);
SimuRslt = simulate_GLSW_interp(IterRslt, struct('num_samples', 6, 'num_periods', 1000));
fprintf('simulate wall-clock: %.1fs\n', toc(t0));

t0 = tic;
irfOpts = struct;
irfOpts.init = struct('shock', [ones(3,1); 2*ones(3,1)], 'a1', zeros(6,1));
irfOpts.num_samples = 6;
irfOpts.num_periods = 100;
rng(0823);
SimuRsltIrf = simulate_GLSW_interp(IterRslt, irfOpts);
fprintf('irf simulate wall-clock: %.1fs\n', toc(t0));

save(fullfile(goldenDir, 'IterRslt.mat'),    'IterRslt',    '-v7');
save(fullfile(goldenDir, 'SimuRslt.mat'),    'SimuRslt',    '-v7');
save(fullfile(goldenDir, 'SimuRsltIrf.mat'), 'SimuRsltIrf', '-v7');
fprintf('Golden captured: Iter=%d Metric=%g\n', IterRslt.Iter, IterRslt.Metric);
end
```

- [ ] **Step 3: Run the three captures (old toolbox; wall-clock unknown — log it)**

```
matlab -batch "cd(fullfile('tests','golden')); capture_safe_assets"
matlab -batch "cd(fullfile('tests','golden')); capture_mendoza2010"
matlab -batch "cd(fullfile('tests','golden')); capture_GLSW_interp"
```

Expected: each prints `Golden captured: Iter=<n> Metric=<m>` with Metric below the model's
TolEq (1e-6, 1e-6, 1e-8 respectively) and `.mat` files appear under each
`tests/<Model>/golden/`. **Record the printed Iter values and wall-clock times** — the
gates and PROGRESS.md need them. If a run is prohibitively long (hours), STOP and
escalate per spec §6 before building gates.

- [ ] **Step 4: Write the golden-integrity tests**

Create `tests/Barro_et_al_2017/tGoldenSafeAssets.m`:

```matlab
classdef tGoldenSafeAssets < matlab.unittest.TestCase
    % Phase 7a: golden files exist, load, and have the expected shape.
    methods (Test)
        function goldenLoadsAndHasShape(tc)
            here = fileparts(mfilename('fullpath'));
            g = load(fullfile(here, 'golden', 'IterRslt.mat'));
            R = g.IterRslt;
            tc.verifyEqual(R.shock_num, 2);
            tc.verifyTrue(isfield(R.var_policy, 'shr_x1'));
            tc.verifyTrue(isfield(R.var_policy, 'omega1n'));
            tc.verifyEqual(numel(R.var_state.omega1), 501);
            tc.verifyLessThan(R.Metric, 1e-6);
            s = load(fullfile(here, 'golden', 'SimuRslt.mat'));
            tc.verifyEqual(size(s.SimuRslt.shock), [6, 1001]);
            tc.verifyTrue(isfield(s.SimuRslt, 'Rf'));
        end
    end
end
```

Create `tests/Mendoza2010/tGoldenMendoza2010.m`:

```matlab
classdef tGoldenMendoza2010 < matlab.unittest.TestCase
    methods (Test)
        function goldenLoadsAndHasShape(tc)
            here = fileparts(mfilename('fullpath'));
            g = load(fullfile(here, 'golden', 'IterRslt.mat'));
            R = g.IterRslt;
            tc.verifyEqual(R.shock_num, 8);
            tc.verifyTrue(isfield(R.var_policy, 'cTildeNext'));
            tc.verifyEqual(numel(R.var_state.cTilde), 80);   % fine-grid re-solve
            tc.verifyEqual(numel(R.var_state.k), 80);
            tc.verifyLessThan(R.Metric, 1e-6);
            % the promoted var_simu names landed in var_aux (old parser behavior)
            tc.verifyTrue(isfield(R.var_aux, 'gdp'));
            g1 = load(fullfile(here, 'golden', 'IterRslt1.mat'));
            tc.verifyLessThan(g1.IterRslt1.Metric, 1e-6);
            s = load(fullfile(here, 'golden', 'SimuRslt.mat'));
            tc.verifyEqual(size(s.SimuRslt.shock), [6, 1001]);
            tc.verifyTrue(isfield(s.SimuRslt, 'gdp'));
        end
    end
end
```

Create `tests/GLSW2020/tGoldenGLSW.m`:

```matlab
classdef tGoldenGLSW < matlab.unittest.TestCase
    methods (Test)
        function goldenLoadsAndHasShape(tc)
            here = fileparts(mfilename('fullpath'));
            g = load(fullfile(here, 'golden', 'IterRslt.mat'));
            R = g.IterRslt;
            tc.verifyEqual(R.shock_num, 2);
            tc.verifyTrue(isfield(R.var_policy, 'a1n'));
            tc.verifyEqual(numel(R.var_state.a1), 301);
            tc.verifyLessThan(R.Metric, 1e-8);
            s = load(fullfile(here, 'golden', 'SimuRslt.mat'));
            tc.verifyEqual(size(s.SimuRslt.shock), [6, 1001]);
            irf = load(fullfile(here, 'golden', 'SimuRsltIrf.mat'));
            tc.verifyEqual(size(irf.SimuRsltIrf.shock), [6, 101]);
            tc.verifyEqual(irf.SimuRsltIrf.shock(:,1), [1;1;1;2;2;2]);
        end
    end
end
```

- [ ] **Step 5: Run the integrity tests**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','Barro_et_al_2017','tGoldenSafeAssets.m')), runtests(fullfile('tests','Mendoza2010','tGoldenMendoza2010.m')), runtests(fullfile('tests','GLSW2020','tGoldenGLSW.m'))]; disp(table(results)); exit(any([results.Failed]))"
```

Expected: PASS (3/3). If a shape assertion fails, the capture diverged from this plan's
assumptions — inspect the golden, fix the assertion or the capture, do not proceed blind.

- [ ] **Step 6: Commit** (gmod copies, capture scripts, goldens, integrity tests)

```bash
git add tests/Barro_et_al_2017 tests/Mendoza2010 tests/GLSW2020 tests/golden
git commit -m "test(golden): capture Barro/Mendoza/GLSW goldens from the old toolbox"
```

---

### Task 3: IR schema — `setup` section + transition `primed` + init-bound name

One schema commit so the HL1996 reference IR is regenerated once.

**Files:**
- Modify: `src/+gdsge/+ir/schema.m`
- Modify: `src/+gdsge/+parser/parseDeclarations.m` (expose `setupText`)
- Modify: `src/+gdsge/+parser/parseSimulate.m` (`primed` flag)
- Modify: `src/+gdsge/+parser/assemblePartialIR.m` (`ir.setup`)
- Modify: `tests/+gdsgefix/minimalIR.m`, `tests/HeatonLucas1996/ir/buildHL1996IR.m`
- Regenerate: `tests/HeatonLucas1996/ir/HL1996.gdsge.json`, `docs/ir-schema.md`
- Test: `tests/parser/tParseSimulate.m`, `tests/ir/` suite (existing)

- [ ] **Step 1: Write the failing tests**

Add to `tests/parser/tParseSimulate.m` inside `methods (Test)`:

```matlab
        function transitionPrimedFlag(tc)
            sim = gdsge.parser.parseSimulate(sprintf( ...
                'num_periods = 10;\nnum_samples = 2;\ninitial k 0.5;\ninitial shock 1;\nvar_simu c;\nk'' = kNext;\nw'' = wn'';\n'));
            tc.verifyEqual(sim.transitions{1}, ...
                struct('state','k','expr','kNext','primed',0));
            tc.verifyEqual(sim.transitions{2}, ...
                struct('state','w','expr','wn','primed',1));
        end
```

- [ ] **Step 2: Run it — verify it fails** (missing `primed` field)

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','parser','tParseSimulate.m')); disp(table(results)); exit(any([results.Failed]))"
```

- [ ] **Step 3: Implement the schema + parser changes**

In `src/+gdsge/+ir/schema.m`:

1. transItem (line 30) becomes:

```matlab
transItem = fStruct(structOf('state', fRef('states'), 'expr', fText(), 'primed', fScalar()));
```

2. Add an init-bound item next to `boundItem` (the init bounds reference `policyInit`
   names, which are outside the `policyAux` ref domain):

```matlab
initBoundItem = fStruct(structOf( ...
    'name', fText(), ...
    'lower', fText(), ...
    'upper', fText(), ...
    'adaptiveFactor', opt(fScalar())));
```

and in the `modelInit` spec replace `'bounds', fList(boundItem)` with
`'bounds', fList(initBoundItem)`.

3. In `s.root` add, immediately after `'params', fList(paramItem), ...`:

```matlab
    'setup',     opt(fText()), ...
```

In `src/+gdsge/+parser/parseSimulate.m` replace the transition branch (lines 21-25) with:

```matlab
        m = regexp(st, '^([A-Za-z_]\w*)''?\s*=\s*(.+)$', 'tokens', 'once');
        if ~isempty(m)
            rhs = strtrim(m{2});
            sim.transitions{end+1} = struct('state', m{1}, ...
                'expr', strtrim(regexprep(rhs, '''', '')), ...
                'primed', double(contains(rhs, ''''))); %#ok<AGROW>
        end
```

In `src/+gdsge/+parser/parseDeclarations.m`, change the `out = struct(...)` (line 63) to
also carry the setup text (the same re-terminated body that was eval'd, WITHOUT the
parser defaults):

```matlab
out = struct('params', {params}, 'shocks', shocks, 'states', states, ...
    'variables', variables, 'bounds', {decls.bounds}, 'interp', {interp}, ...
    'ws', ws, 'setupText', setupBody);
```

In `src/+gdsge/+parser/assemblePartialIR.m`, after `ir.params = decl.params;` add:

```matlab
ir.setup     = decl.setupText;
```

In `tests/+gdsgefix/minimalIR.m`: add after the `ir.params` line:

```matlab
ir.setup = sprintf('beta = 0.95;\nz = 1.0;\nshock_num = 1;\nshock_trans = 1.0;\nk = linspace(0,1,3);');
```

and change its transition to `struct('state','k','expr','c','primed',0)`.

- [ ] **Step 4: Update the HL1996 reference IR**

In `tests/HeatonLucas1996/ir/buildHL1996IR.m`: add `'primed', 1` to the `w1`
transition struct (HL1996: `w1' = w1n'`), and add the `setup` field. Obtain the exact
text the parser produces (after Step 3 the parser emits it):

```
matlab -batch "addpath('src'); ir = gdsge.parser.parseFrontEnd(fileread(fullfile('tests','HeatonLucas1996','HL1996.gmod')),'HL1996'); fid=fopen('setup_dump.txt','w'); fwrite(fid, ir.setup); fclose(fid);"
```

Paste `setup_dump.txt` verbatim into `buildHL1996IR.m` as the `ir.setup` literal
(`sprintf`-escape newlines and quotes), then delete `setup_dump.txt`.

- [ ] **Step 5: Regenerate the JSON golden + schema doc**

```
matlab -batch "addpath('src'); addpath(fullfile('tests','HeatonLucas1996','ir')); ir = buildHL1996IR(); fid=fopen(fullfile('tests','HeatonLucas1996','ir','HL1996.gdsge.json'),'w'); fwrite(fid, gdsge.ir.encode(ir)); fclose(fid);"
matlab -batch "addpath('src'); gdsge.ir.gendoc(fullfile('docs','ir-schema.md'));"
```

Review `git diff` on both files: the JSON should gain exactly `setup` and the `primed`
flags; `docs/ir-schema.md` should gain the new fields. Unexplained churn = stop and debug.

- [ ] **Step 6: Run the parser + ir suites — verify green**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','parser'),'IncludingSubfolders',true), runtests(fullfile('tests','ir'),'IncludingSubfolders',true), runtests(fullfile('tests','HeatonLucas1996','parser')), runtests(fullfile('tests','HeatonLucas1996','ir'))]; disp(table(results)); exit(any([results.Failed]))"
```

Expected: PASS. `tFullIRHL1996` (parser output == reference) is the canary here.

- [ ] **Step 7: Commit**

```bash
git add src/+gdsge/+ir/schema.m src/+gdsge/+parser tests/+gdsgefix tests/HeatonLucas1996/ir tests/parser docs/ir-schema.md
git commit -m "feat(ir): setup section + transition primed flag + init-bound item"
```

---

### Task 4: `emitSetup` replays the gmod setup; primed-aware transitions

This fixes the two latent generators bugs Barro/Mendoza/GLSW would hit: undefined setup
intermediates (`Re_n`, `Ngrid`, `a1_lb`…) and wrongly-indexed unprimed transitions.

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitSetup.m`, `src/+gdsge/+codegen/+mat/emitResultSimu.m`
- Regenerate: `tests/HeatonLucas1996/codegen/golden/*` snapshots
- Test: existing snapshot + functional suites

- [ ] **Step 1: Modify `emitSetup`** — replace the three numeric sections (lines 56-79:
  `%% ---- model parameters`, `%% ---- shocks`, `%% ---- state grids`) with a verbatim
  replay of `ir.setup` (old-generator `setParamsCode` semantics — the setup code defines
  params, shock values, `shock_num`, `shock_trans`, state grids, and any intermediates
  like `Re_n` in原 original order, exactly once):

```matlab
w.add('%% ---- gmod setup (declaration-region code replayed verbatim, old setParamsCode parity)');
if isfield(ir, 'setup') && ~isempty(ir.setup)
    w.addRaw(ir.setup);
end
```

Keep everything above (options block, num_samples/num_periods) unchanged. Update the
file's header comment to describe the replay.

- [ ] **Step 2: Modify `emitResultSimu`** — transitions branch on `primed`
  (old parity: primed RHS → linear index; unprimed → `(:)`):

```matlab
for i = 1:numel(ir.simulate.transitions)
    t = ir.simulate.transitions{i};
    if t.primed
        w4.add('SimuRslt.%s(:,GDSGE_t+1) = %s(GDSGE_SHOCK_VAR_LINEAR_INDEX);', t.state, t.expr);
    else
        w4.add('SimuRslt.%s(:,GDSGE_t+1) = %s(:);', t.state, t.expr);
    end
end
```

- [ ] **Step 3: Regenerate the HL1996 codegen snapshots and review**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); addpath(fullfile('tests','HeatonLucas1996','ir')); run(fullfile('tests','HeatonLucas1996','codegen','regen_snapshots.m'))"
```

`git diff` the four golden txt files. Expected change: the params/shocks/grids numeric
lines replaced by the HL1996 gmod setup text, and nothing else in iter; simulate
additionally keeps its transition line identical (HL1996's only transition is primed).
The C++ snapshots must be **byte-identical** (no C++ change in this task) — if they
differ, stop and debug.

- [ ] **Step 4: Run the fast suites + the HL1996 Slow end-to-end gate**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','HeatonLucas1996'),'IncludingSubfolders',true); disp(table(results)); exit(any([results.Failed]))"
```

Expected: PASS including `tEndToEndHL1996` (budget ~10 min: compile + 209 iterations +
simulation). The replayed setup must reproduce the exact same numbers.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/+mat tests/HeatonLucas1996/codegen/golden
git commit -m "feat(codegen-mat): replay gmod setup verbatim; primed-aware simulate transitions"
```

---

### Task 5: `resolveOptions` extracts `SimuPrintFreq` / `SimuSaveFreq`

**Files:**
- Modify: `src/+gdsge/+parser/resolveOptions.m`
- Test: `tests/parser/tResolveOptions.m`

- [ ] **Step 1: Write the failing test** — add to `tests/parser/tResolveOptions.m`:

```matlab
        function simuFreqsFlowThrough(tc)
            ws = struct('SimuPrintFreq', 10000, 'SimuSaveFreq', inf);
            o = gdsge.parser.resolveOptions(ws);
            tc.verifyEqual(o.simuPrintFreq, 10000);
            tc.verifyEqual(o.simuSaveFreq, inf);
        end
```

- [ ] **Step 2: Run it — verify it fails** (no such fields)

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','parser','tResolveOptions.m')); disp(table(results)); exit(any([results.Failed]))"
```

- [ ] **Step 3: Implement** — in `resolveOptions.m` after the `o.saveFreq` line add:

```matlab
o.simuPrintFreq = getf(ws, 'SimuPrintFreq', 1000);
o.simuSaveFreq  = getf(ws, 'SimuSaveFreq', Inf);
```

(`emitSetup` already reads `opt.simuPrintFreq`/`opt.simuSaveFreq` with these defaults;
the schema already declares both fields.)

- [ ] **Step 4: Run the test — verify PASS.** Note: `tFullIRHL1996` compares parser
  output to the reference IR — the reference's options struct must gain the two fields.
  Add to `buildHL1996IR.m`'s options struct: `'simuPrintFreq', 1000, 'simuSaveFreq', Inf`
  and regenerate `HL1996.gdsge.json` (same command as Task 3 Step 5), then re-run the
  HL1996 parser/ir suites.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/resolveOptions.m tests/parser/tResolveOptions.m tests/HeatonLucas1996/ir
git commit -m "feat(parser): SimuPrintFreq/SimuSaveFreq flow into IR options"
```

---

### Task 6: `resolveOutputs` — var_simu → var_output → var_aux promotion

**Files:**
- Create: `src/+gdsge/+parser/resolveOutputs.m`
- Modify: `src/+gdsge/+parser/parseFrontEnd.m`
- Test: `tests/parser/tResolveOutputs.m`

- [ ] **Step 1: Write the failing tests**

Create `tests/parser/tResolveOutputs.m`:

```matlab
classdef tResolveOutputs < matlab.unittest.TestCase
    % Old-generator parity: var_simu names promote into var_output (skipping
    % states and existing outputs); output names missing from policy/aux are
    % appended to aux with length 1, extending the slot layout.
    methods (Test)
        function promotesSimuToOutputAndAux(tc)
            [v2, out] = runPromote(vars(), {'c','gdp','mu'}, {'k'});
            tc.verifyEqual(v2.output, {'c','gdp','mu'});
            % c is aux already; gdp and mu are new: mu is policy (no aux add),
            % gdp is neither -> appended to aux after slot 2
            names = cellfun(@(a) a.name, v2.aux, 'UniformOutput', false);
            tc.verifyEqual(names, {'c','gdp'});
            tc.verifyEqual(v2.aux{2}.slot, [3 3]);   % c spans [1 2]
            tc.verifyTrue(contains(out, 'added to var_output'));
            tc.verifyTrue(contains(out, 'added to var_aux'));
        end
        function skipsStatesAndExisting(tc)
            v = vars();
            v.output = {'c'};
            v2 = runPromote(v, {'k','c'}, {'k'});
            tc.verifyEqual(v2.output, {'c'});        % k is a state, c existing
            tc.verifyEqual(numel(v2.aux), 1);        % nothing promoted to aux
        end
        function errorsOnDeclaredStateOutput(tc)
            v = vars();
            v.output = {'k'};
            tc.verifyError(@() gdsge.parser.resolveOutputs(v, {}, {'k'}), ...
                'gdsge:parser:stateInOutput');
        end
    end
end

function v = vars()
v = struct( ...
    'policy', {{ struct('name','mu','length',1,'slot',[1 1]) }}, ...
    'aux',    {{ struct('name','c','length',2,'slot',[1 2]) }}, ...
    'interp', {{}}, 'tensor', {{}}, 'output', {{}}, 'others', {{}});
end

function [v2, out] = runPromote(v, varSimu, states)
out = evalc('v2 = gdsge.parser.resolveOutputs(v, varSimu, states);');
end
```

- [ ] **Step 2: Run it — verify it fails** (`resolveOutputs` undefined).

- [ ] **Step 3: Implement** — create `src/+gdsge/+parser/resolveOutputs.m`:

```matlab
function variables = resolveOutputs(variables, varSimu, stateNames)
% RESOLVEOUTPUTS  Old-generator output promotion (gdsge_parser.m:619-685):
%   1) error if a declared var_output is a state;
%   2) var_simu names that are neither states nor existing outputs append to
%      variables.output, in var_simu order;
%   3) output names not found in policy or aux append to variables.aux with
%      length 1 (model-body locals the MEX must export), extending the slot
%      layout. Prints the same informational lines as the old parser.
for i = 1:numel(variables.output)
    if ismember(variables.output{i}, stateNames)
        error('gdsge:parser:stateInOutput', ...
            'state %s should not be in var_output', variables.output{i});
    end
end
added = {};
for i = 1:numel(varSimu)
    n = varSimu{i};
    if ismember(n, stateNames) || ismember(n, variables.output); continue; end
    variables.output{end+1} = n;
    added{end+1} = n; %#ok<AGROW>
end
if ~isempty(added)
    fprintf('The following var_simulate are added to var_output: %s\n', strjoin(added, ' '));
end

polNames = cellfun(@(p) p.name, variables.policy, 'UniformOutput', false);
auxNames = cellfun(@(a) a.name, variables.aux, 'UniformOutput', false);
pos = 1;
for i = 1:numel(variables.aux); pos = max(pos, variables.aux{i}.slot(2) + 1); end
added = {};
for i = 1:numel(variables.output)
    n = variables.output{i};
    if ismember(n, polNames) || ismember(n, auxNames); continue; end
    variables.aux{end+1} = struct('name', n, 'length', 1, 'slot', [pos, pos]);
    auxNames{end+1} = n; %#ok<AGROW>
    pos = pos + 1;
    added{end+1} = n; %#ok<AGROW>
end
if ~isempty(added)
    fprintf('The following var_output are added to var_aux: %s\n', strjoin(added, ' '));
end
end
```

In `parseFrontEnd.m`, after the `sim = gdsge.parser.parseSimulate(...)` line add:

```matlab
decl.variables = gdsge.parser.resolveOutputs(decl.variables, sim.varSimu, decl.states.names);
```

- [ ] **Step 4: Run tests — verify PASS**, then run the full HL1996 folder (HL1996's
  var_simu ⊆ var_output, so this must be a **no-op** there — `tFullIRHL1996` and the
  snapshots confirm):

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','parser','tResolveOutputs.m')), runtests(fullfile('tests','HeatonLucas1996','parser'))]; disp(table(results)); exit(any([results.Failed]))"
```

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/resolveOutputs.m src/+gdsge/+parser/parseFrontEnd.m tests/parser/tResolveOutputs.m
git commit -m "feat(parser): var_simu -> var_output -> var_aux promotion (old parity)"
```

### Task 7: Parser front-end — `*_init` declarations + `model_init` block

**Files:**
- Modify: `src/+gdsge/+parser/parseVarDecls.m`, `parseDeclarations.m`, `parseFrontEnd.m`, `assemblePartialIR.m`, `analyzeModel.m`
- Test: `tests/parser/tParseVarDecls.m`, `tests/parser/tParseDeclarations.m`, `tests/parser/tModelInitFrontEnd.m` (new)

- [ ] **Step 1: Write the failing tests**

Add to `tests/parser/tParseVarDecls.m`:

```matlab
        function initDeclarationsParsed(tc)
            d = gdsge.parser.parseVarDecls(sprintf( ...
                'var_policy_init L v;\nvar_aux_init flow q d;\ninbound_init L 0 1000;\ninbound_init v 0 1000;'));
            tc.verifyEqual(numel(d.policyInit), 2);
            tc.verifyEqual(d.policyInit{1}.name, 'L');
            tc.verifyEqual(numel(d.auxInit), 3);
            tc.verifyEqual(numel(d.boundsInit), 2);
            tc.verifyEqual(d.boundsInit{2}, struct('name','v','lower','0','upper','1000'));
            tc.verifyEmpty(d.bounds);   % inbound_init no longer leaks into bounds
        end
```

Create `tests/parser/tModelInitFrontEnd.m` (a tiny synthetic model with an init block;
1 shock so the system stays small):

```matlab
classdef tModelInitFrontEnd < matlab.unittest.TestCase
    % model_init flows gmod -> ir.modelInit -> validate -> JSON round-trip.
    methods (Test)
        function modelInitPopulatesIR(tc)
            ir = gdsge.parser.parseFrontEnd(gmodText(), 'tiny');
            tc.verifyTrue(isfield(ir, 'modelInit'));
            mi = ir.modelInit;
            tc.verifyEqual(mi.variables.policyInit{1}, ...
                struct('name','c0','length',1,'slot',[1 1]));
            tc.verifyEqual(mi.variables.auxInit{1}.name, 'q0');
            tc.verifyEqual(numel(mi.bounds), 1);
            tc.verifyEqual(numel(mi.equations), 1);
            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));
        end
        function initBlockWithoutVarsErrors(tc)
            txt = strrep(gmodText(), sprintf('var_policy_init c0;\n'), '');
            txt = strrep(txt, sprintf('inbound_init c0 0 2;\n'), '');
            tc.verifyError(@() gdsge.parser.parseFrontEnd(txt, 'tiny'), ...
                'gdsge:parser:modelInitNoVars');
        end
        function initVarsWithoutBlockErrors(tc)
            txt = regexprep(gmodText(), 'model_init;.*?end;\s*model;', 'model;');
            tc.verifyError(@() gdsge.parser.parseFrontEnd(txt, 'tiny'), ...
                'gdsge:parser:initVarsNoBlock');
        end
        function initPolicyNeedsInitBound(tc)
            txt = strrep(gmodText(), sprintf('inbound_init c0 0 2;\n'), '');
            tc.verifyError(@() gdsge.parser.parseFrontEnd(txt, 'tiny'), ...
                'gdsge:parser:missingBound');
        end
        function initSystemMustBeSquare(tc)
            txt = strrep(gmodText(), 'var_policy_init c0;', 'var_policy_init c0 d0;');
            txt = strrep(txt, sprintf('inbound_init c0 0 2;\n'), ...
                sprintf('inbound_init c0 0 2;\ninbound_init d0 0 2;\n'));
            tc.verifyError(@() gdsge.parser.parseFrontEnd(txt, 'tiny'), ...
                'gdsge:parser:notSquare');
        end
    end
end

function txt = gmodText()
txt = sprintf([ ...
    'parameters beta;\n' ...
    'beta = 0.95;\n' ...
    'var_shock z;\n' ...
    'shock_num = 1;\n' ...
    'shock_trans = 1;\n' ...
    'z = 1;\n' ...
    'var_state k;\n' ...
    'k = linspace(0.5,1.5,3);\n' ...
    'var_policy_init c0;\n' ...
    'var_aux_init q0;\n' ...
    'inbound_init c0 0 2;\n' ...
    'var_policy c;\n' ...
    'inbound c 0 2;\n' ...
    'var_interp vf;\n' ...
    'initial vf q0;\n' ...
    'vf = q;\n' ...
    'model_init;\n' ...
    '  q0 = beta*k;\n' ...
    '  equations;\n' ...
    '    c0 - q0;\n' ...
    '  end;\n' ...
    'end;\n' ...
    'model;\n' ...
    '  q = beta*c;\n' ...
    '  equations;\n' ...
    '    c - k;\n' ...
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

(In this fixture `var_simu c` promotes `c` into var_output; `c` is policy, so no aux
promotion happens. `q` is a model-body local consumed by the interp update — interp
updates are opaque text evaluated in the generated iter, so only parsing/validation
matters here.)

- [ ] **Step 2: Run both test files — verify the new tests fail** with missing-field /
  missing-error-id failures.

- [ ] **Step 3: Implement the parser changes**

`parseVarDecls.m`:
- Extend the result template (line 8-11) with `'policyInit',{{}}, 'auxInit',{{}}, 'boundsInit',{{}}`.
- Replace the `case {'inbound','inbound_init'}` and `case {'var_policy_init','var_aux_init','model_init'}` arms with:

```matlab
        case 'inbound'
            decls.bounds{end+1} = parseBound(rest(st,key)); %#ok<AGROW>
        case 'inbound_init'
            decls.boundsInit{end+1} = parseBound(rest(st,key)); %#ok<AGROW>
        case 'var_policy_init'
            decls.policyInit = [decls.policyInit, parseSizedNames(rest(st,key))];
        case 'var_aux_init'
            decls.auxInit = [decls.auxInit, parseSizedNames(rest(st,key))];
```

(The defensive `model_init` keyword arm is removed — `splitBlocks` extracts the block
before declText is built, so the token can no longer appear here.)

`parseDeclarations.m` — extend the output struct (now with init tables; note the `{}`
wrapping for cell-valued fields inside `struct()`):

```matlab
out = struct('params', {params}, 'shocks', shocks, 'states', states, ...
    'variables', variables, 'bounds', {decls.bounds}, 'interp', {interp}, ...
    'ws', ws, 'setupText', setupBody, ...
    'boundsInit', {decls.boundsInit}, ...
    'initVariables', struct( ...
        'policyInit', {assignSlots(decls.policyInit)}, ...
        'auxInit',    {assignSlots(decls.auxInit)}));
```

`parseFrontEnd.m` — after the `model = gdsge.parser.parseModel(...)` line:

```matlab
modelInit = [];
hasInitVars = ~isempty(decl.initVariables.policyInit);
if isfield(sb.blocks, 'model_init')
    if ~hasInitVars
        error('gdsge:parser:modelInitNoVars', ...
            'model_init block present but no var_policy_init declared');
    end
    initBody = gdsge.parser.parseModel(sb.blocks.model_init);
    modelInit = struct( ...
        'variables',  decl.initVariables, ...
        'bounds',     {decl.boundsInit}, ...
        'statements', {initBody.statements}, ...
        'equations',  {initBody.equations});
elseif hasInitVars
    error('gdsge:parser:initVarsNoBlock', ...
        'var_policy_init declared but no model_init block found');
end
```

and pass it through: `ir = gdsge.parser.assemblePartialIR(modelName, decl, sim, options, hooks, model, modelInit);`

`assemblePartialIR.m` — signature gains the optional 7th arg:

```matlab
function ir = assemblePartialIR(modelName, decl, sim, options, hooks, model, modelInit)
...
if nargin < 7; modelInit = []; end
```

after `ir.model = model;` add:

```matlab
if ~isempty(modelInit)
    ir.modelInit = modelInit;
end
```

and after the existing bounds-completeness gate add the init twin:

```matlab
if ~isempty(modelInit)
    bn = cellfun(@(b) b.name, modelInit.bounds, 'UniformOutput', false);
    for i = 1:numel(modelInit.variables.policyInit)
        pn = modelInit.variables.policyInit{i}.name;
        if ~ismember(pn, bn)
            error('gdsge:parser:missingBound', ...
                'Init policy variable "%s" has no inbound_init', pn);
        end
    end
end
```

`analyzeModel.m` — at the end of the main function body add:

```matlab
if isfield(ir, 'modelInit')
    analyzeInit(ir);
end
```

and append the local function (reuses `checkNames`):

```matlab
function analyzeInit(ir)
% Same discipline as the main pass, with the init variable tables. The init
% problem runs before any interpolant exists, so interp calls are an error;
% only shocks are primeable (no policy priming in a static problem).
known = {'GDSGE_Iter','TASK','shock'};
for i = 1:numel(ir.params); known{end+1} = ir.params{i}.name; end %#ok<AGROW>
known = [known, reshape(ir.states.names, 1, []), reshape(ir.shocks.names, 1, [])];
polNames = cellfun(@(p) p.name, ir.modelInit.variables.policyInit, 'UniformOutput', false);
auxNames = cellfun(@(a) a.name, ir.modelInit.variables.auxInit, 'UniformOutput', false);
known = [known, reshape(polNames, 1, []), reshape(auxNames, 1, [])];
primeable = reshape(ir.shocks.names, 1, []);

stmts = ir.modelInit.statements;
for i = 1:numel(stmts)
    s = stmts{i};
    switch s.type
        case 'assign'
            checkNames(s.expr, known, primeable, ...
                sprintf('model_init statement %d (assign %s)', i, s.target));
            known{end+1} = s.target; %#ok<AGROW>
        case 'reduction'
            if ~isfield(ir.shocks.transitions, s.transRef)
                error('gdsge:parser:badTransRef', ...
                    'model_init statement %d: reduction transition "%s" is not declared', ...
                    i, s.transRef);
            end
            checkNames(s.body, known, primeable, ...
                sprintf('model_init statement %d (reduction %s)', i, s.target));
            known{end+1} = s.target; %#ok<AGROW>
        case 'interpCall'
            error('gdsge:parser:interpInInit', ...
                'model_init statement %d: GDSGE_INTERP_VEC is not available in model_init', i);
    end
end
for i = 1:numel(ir.modelInit.equations)
    checkNames(ir.modelInit.equations{i}.expr, known, primeable, ...
        sprintf('model_init equation %d', i));
end

nUnknown = 0;
for i = 1:numel(ir.modelInit.variables.policyInit)
    nUnknown = nUnknown + ir.modelInit.variables.policyInit{i}.length;
end
nEq = 0;
for i = 1:numel(ir.modelInit.equations)
    if ir.modelInit.equations{i}.primed
        nEq = nEq + ir.shocks.count;
    else
        nEq = nEq + 1;
    end
end
if nUnknown ~= nEq
    error('gdsge:parser:notSquare', ...
        'model_init system is not square: %d unknowns vs %d equations', nUnknown, nEq);
end
end
```

- [ ] **Step 4: Run the parser suite — verify PASS** (including the five new
  tModelInitFrontEnd tests and the untouched HL1996 gates):

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','parser'),'IncludingSubfolders',true), runtests(fullfile('tests','HeatonLucas1996','parser'))]; disp(table(results)); exit(any([results.Failed]))"
```

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser tests/parser
git commit -m "feat(parser): model_init block + var_policy_init/var_aux_init/inbound_init -> ir.modelInit"
```

---

### Task 8: Front-end gates for the three models

**Files:**
- Create: `tests/Barro_et_al_2017/parser/tFrontEndSafeAssets.m`
- Create: `tests/Mendoza2010/parser/tFrontEndMendoza2010.m`
- Create: `tests/GLSW2020/parser/tFrontEndGLSW.m`

- [ ] **Step 1: Write the three gates** (these are expected to pass if Tasks 3-7 are
  complete — they are the front-end acceptance tests; any failure here is a real
  front-end gap to fix now, before codegen)

Create `tests/Barro_et_al_2017/parser/tFrontEndSafeAssets.m`:

```matlab
classdef tFrontEndSafeAssets < matlab.unittest.TestCase
    % Phase 7a front-end gate: safe_assets.gmod -> schema-valid IR.
    methods (Test)
        function parsesValidatesRoundTrips(tc)
            ir = parseIt();
            tc.verifyEqual(ir.shocks.count, 2);
            tc.verifyEqual(ir.states.names, {'omega1'});
            % policy: shr_x1, Rf, omega1n[2] -> 4 unknowns
            tc.verifyEqual(cellfun(@(p) p.length, ir.variables.policy), [1 1 2]);
            % equations: eq1, eq2, omega_future_consis' (primed)
            tc.verifyEqual(numel(ir.model.equations), 3);
            tc.verifyEqual(ir.model.equations{3}.primed, 1);
            tc.verifyFalse(isfield(ir, 'modelInit'));
            % bound text references the setup intermediate Re_n
            tc.verifyTrue(contains(ir.setup, 'Re_n'));
            bounds = ir.bounds; rf = bounds{2};
            tc.verifyEqual(rf.name, 'Rf');
            tc.verifyEqual(rf.lower, 'Re_n(2)');
            % GDSGE_MIN parsed as a reduction statement
            kinds = {};
            for i = 1:numel(ir.model.statements)
                s = ir.model.statements{i};
                if strcmp(s.type, 'reduction'); kinds{end+1} = s.kind; end %#ok<AGROW>
            end
            tc.verifyTrue(ismember('MIN', kinds));
            % no var_output declared: all four var_simu promoted
            tc.verifyEqual(ir.variables.output, {'Rf','K1','b1','expectedRe'});
            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));
        end
    end
end

function ir = parseIt()
here = fileparts(mfilename('fullpath'));
gmodPath = fullfile(fileparts(here), 'safe_assets.gmod');
ir = gdsge.parser.parseFrontEnd(fileread(gmodPath), 'safe_assets');
end
```

Create `tests/Mendoza2010/parser/tFrontEndMendoza2010.m`:

```matlab
classdef tFrontEndMendoza2010 < matlab.unittest.TestCase
    % Phase 7a front-end gate: mendoza2010.gmod -> schema-valid IR with
    % modelInit, promotion, SIMU_INTERP options, adaptive bounds, 2-D states.
    methods (Test)
        function parsesValidatesRoundTrips(tc)
            ir = parseIt();
            tc.verifyEqual(ir.shocks.count, 8);
            tc.verifyEqual(ir.states.names, {'cTilde','k'});
            % policy: L v mu nbNext kNext cTildeNext[8] -> 13 unknowns
            tc.verifyEqual(sum(cellfun(@(p) p.length, ir.variables.policy)), 13);
            tc.verifyEqual(numel(ir.model.equations), 6);
            tc.verifyEqual(ir.model.equations{6}.primed, 1);

            % ---- modelInit
            tc.verifyTrue(isfield(ir, 'modelInit'));
            mi = ir.modelInit;
            pn = cellfun(@(p) p.name, mi.variables.policyInit, 'UniformOutput', false);
            tc.verifyEqual(pn, {'L','v'});
            an = cellfun(@(a) a.name, mi.variables.auxInit, 'UniformOutput', false);
            tc.verifyEqual(an, {'flow','q','d'});
            tc.verifyEqual(numel(mi.bounds), 2);
            tc.verifyEqual(numel(mi.equations), 2);

            % ---- options
            tc.verifyEqual(ir.options.simuInterp, 1);
            tc.verifyEqual(ir.options.simuResolve, 0);
            tc.verifyEqual(ir.options.saveFreq, inf);

            % ---- adaptive bounds present (5 of them, factor 1.5)
            nAdaptive = 0;
            for i = 1:numel(ir.bounds)
                if isfield(ir.bounds{i}, 'adaptiveFactor'); nAdaptive = nAdaptive + 1; end
            end
            tc.verifyEqual(nAdaptive, 5);

            % ---- promotion: declared output {cTildeNext,kNext} + 14 var_simu
            tc.verifyEqual(ir.variables.output(1:2), {'cTildeNext','kNext'});
            tc.verifyEqual(numel(ir.variables.output), 16);
            auxNames = cellfun(@(a) a.name, ir.variables.aux, 'UniformOutput', false);
            tc.verifyTrue(all(ismember({'c','Y','inv','b','bNext','nx','gdp', ...
                'b_gdp','nx_gdp','lev','wkcptl'}, auxNames)));

            % ---- transitions: cTilde' = cTildeNext' (primed), k' = kNext (not)
            tc.verifyEqual(ir.simulate.transitions{1}.primed, 1);
            tc.verifyEqual(ir.simulate.transitions{2}.primed, 0);

            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));
        end
    end
end

function ir = parseIt()
here = fileparts(mfilename('fullpath'));
gmodPath = fullfile(fileparts(here), 'mendoza2010.gmod');
ir = gdsge.parser.parseFrontEnd(fileread(gmodPath), 'mendoza2010');
end
```

Create `tests/GLSW2020/parser/tFrontEndGLSW.m`:

```matlab
classdef tFrontEndGLSW < matlab.unittest.TestCase
    % Phase 7a front-end gate: GLSW_interp.gmod -> schema-valid IR.
    methods (Test)
        function parsesValidatesRoundTrips(tc)
            ir = parseIt();
            tc.verifyEqual(ir.shocks.count, 2);
            tc.verifyEqual(ir.states.names, {'a1'});
            tc.verifyEqual(sum(cellfun(@(p) p.length, ir.variables.policy)), 5);
            tc.verifyEqual(numel(ir.model.equations), 5);

            % ---- modelInit: c1_shr / {P1, log_lambda1, log_lambda2} / 1 eq
            tc.verifyTrue(isfield(ir, 'modelInit'));
            mi = ir.modelInit;
            tc.verifyEqual(mi.variables.policyInit{1}.name, 'c1_shr');
            tc.verifyEqual(numel(mi.variables.auxInit), 3);
            tc.verifyEqual(numel(mi.equations), 1);

            % ---- options: the GLSW flag block
            tc.verifyEqual(ir.options.tolEq, 1e-8);
            tc.verifyEqual(ir.options.interpOrder, 4);
            tc.verifyEqual(ir.options.simuInterp, 1);
            tc.verifyEqual(ir.options.simuPrintFreq, 10000);
            tc.verifyEqual(ir.options.simuSaveFreq, inf);
            tc.verifyEqual(ir.options.saveFreq, inf);

            % ---- grid text references setup intermediates; replay carries them
            tc.verifyTrue(contains(ir.setup, 'Ngrid'));
            tc.verifyTrue(contains(ir.setup, 'a1_lb'));

            % ---- var_simu {a2,P1,r,c1_shr} all already declared outputs: no growth
            tc.verifyEqual(ir.variables.output, {'a1n','a2','P1','r','c1_shr'});

            % ---- transition a1' = a1n (unprimed RHS)
            tc.verifyEqual(ir.simulate.transitions{1}, ...
                struct('state','a1','expr','a1n','primed',0));

            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));
        end
    end
end

function ir = parseIt()
here = fileparts(mfilename('fullpath'));
gmodPath = fullfile(fileparts(here), 'GLSW_interp.gmod');
ir = gdsge.parser.parseFrontEnd(fileread(gmodPath), 'GLSW_interp');
end
```

- [ ] **Step 2: Run all three — fix any real front-end gap they expose**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','Barro_et_al_2017','parser')), runtests(fullfile('tests','Mendoza2010','parser')), runtests(fullfile('tests','GLSW2020','parser'))]; disp(table(results)); exit(any([results.Failed]))"
```

Expected: PASS. Likely first-run failures to investigate (these are known-unknowns, fix
with a failing unit test in `tests/parser/` first, then the fix):
- Barro `inbound Rf Re_n(2) Re_n(1)` — `parseBound` splits on whitespace; `Re_n(2)` is
  one token, fine. But Barro's `omega1n` bound names a **sized** policy (`omega1n[2]`):
  bound name lookup must match the base name.
- Mendoza's `[flow_future',q_plus_d_future'] = GDSGE_INTERP_VEC'(cTildeNext',kNext);`
  — a 2-target primed interp call with a non-primed second argument; Phase 3 built this
  for HL1996's 4-target form, so it should parse.
- assertion text mismatches (e.g. exact output ordering) — verify against the IR by
  printing, then correct the assertion ONLY if the IR is genuinely right (old-parity
  promotion order: declared outputs first, then var_simu order).

- [ ] **Step 3: Commit**

```bash
git add tests/Barro_et_al_2017/parser tests/Mendoza2010/parser tests/GLSW2020/parser
git commit -m "test(parser): Phase 7a front-end gates - safe_assets, mendoza2010, GLSW_interp"
```

---

### Task 9: C++ init task (`task_init`) + compile defines

**Files:**
- Create: `src/+gdsge/+codegen/initView.m`
- Modify: `src/+gdsge/+codegen/+cxx/emitTask.m`, `src/+gdsge/+codegen/+cxx/emitCompile.m`, `src/+gdsge/+codegen/generateCxx.m`
- Create: `tests/+gdsgefix/minimalIRWithInit.m`
- Test: `tests/codegen/tInitView.m` (new), `tests/codegen/tGenerateCxxInit.m` (new)

- [ ] **Step 1: Write the failing tests**

Create `tests/+gdsgefix/minimalIRWithInit.m`:

```matlab
function ir = minimalIRWithInit()
% minimalIR + a one-unknown model_init section (c0 with bound, q0 aux,
% one equation). Used by init-task codegen tests.
ir = gdsgefix.minimalIR();
B = @gdsge.ir.node.binop; NM = @gdsge.ir.node.name;
ir.modelInit = struct( ...
    'variables', struct( ...
        'policyInit', {{ struct('name','c0','length',1,'slot',[1 1]) }}, ...
        'auxInit',    {{ struct('name','q0','length',1,'slot',[1 1]) }}), ...
    'bounds',     {{ struct('name','c0','lower','0','upper','2') }}, ...
    'statements', {{ struct('type','assign','target','q0','primed',0, ...
                            'expr', B('*', NM('beta'), NM('k'))) }}, ...
    'equations',  {{ struct('expr', B('-', NM('c0'), NM('q0')), 'primed', 0) }});
end
```

Create `tests/codegen/tInitView.m`:

```matlab
classdef tInitView < matlab.unittest.TestCase
    methods (Test)
        function mapsInitSectionsToMainShape(tc)
            iv = gdsge.codegen.initView(gdsgefix.minimalIRWithInit());
            tc.verifyEqual(iv.variables.policy{1}.name, 'c0');
            tc.verifyEqual(iv.variables.aux{1}.name, 'q0');
            tc.verifyEmpty(iv.variables.interp);
            tc.verifyEmpty(iv.interp);
            tc.verifyEqual(iv.bounds{1}.name, 'c0');
            tc.verifyEqual(numel(iv.model.equations), 1);
            tc.verifyEqual(iv.hooks.preModel, '');
            tc.verifyFalse(isfield(iv, 'modelInit'));
        end
    end
end
```

Create `tests/codegen/tGenerateCxxInit.m`:

```matlab
classdef tGenerateCxxInit < matlab.unittest.TestCase
    methods (Test)
        function emitsTaskInitAndDefines(tc)
            work = tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            files = gdsge.codegen.generateCxx(gdsgefix.minimalIRWithInit(), work);
            cpp = fileread(files.cppFile);
            tc.verifyTrue(contains(cpp, 'void task_init'));
            tc.verifyTrue(contains(cpp, 'void task_inf_horizon'));
            comp = fileread(files.compileFile);
            tc.verifyTrue(contains(comp, '-DHAS_INIT'));
            % MAXDIM = max(main policy total, init policy total) + 4 = max(1,1)+4
            tc.verifyTrue(contains(comp, '-DMAXDIM=5'));
        end
        function noInitStaysClean(tc)
            work = tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            files = gdsge.codegen.generateCxx(gdsgefix.minimalIR(), work);
            comp = fileread(files.compileFile);
            tc.verifyFalse(contains(comp, '-DHAS_INIT'));
        end
    end
end
```

- [ ] **Step 2: Run them — verify they fail** (`initView` undefined; generateCxx still
  refuses modelInit).

- [ ] **Step 3: Implement**

Create `src/+gdsge/+codegen/initView.m`:

```matlab
function irv = initView(ir)
% INITVIEW  A main-model-shaped view of ir.modelInit so every existing emitter
%   (emitModel/emitModelBody/emitArgUnpack/emitDeclare/emitEquations/emitAux/
%   emitPop/emitInterp/emitBounds) is reused verbatim for the init task.
%   The init problem has no interpolants and no hook code, and its GDSGE_DATA
%   layout is identical to the main task's (dataLayout ignores policy/aux).
irv = ir;
irv.variables.policy = ir.modelInit.variables.policyInit;
irv.variables.aux    = ir.modelInit.variables.auxInit;
irv.variables.interp = {};
irv.interp = {};
irv.bounds = ir.modelInit.bounds;
irv.model  = struct('statements', {ir.modelInit.statements}, ...
                    'equations',  {ir.modelInit.equations});
irv.hooks  = struct('preModel','','preIter','','postIter','', ...
                    'preJacCode','','postJacCode','','cxx','');
irv = rmfield(irv, 'modelInit');
end
```

`emitTask.m` — add the parameter:

```matlab
function txt = emitTask(ir, taskName)
% EMITTASK  One task function from the task template. taskName selects which
%   (task_inf_horizon for the main model, task_init for ir-as-initView).
if nargin < 2; taskName = 'task_inf_horizon'; end
```

and at the bottom of the fill list: `'TASK_NAME', taskName});`

`generateCxx.m` — delete the modelInit refusal (lines 15-17) and replace the template
fill with:

```matlab
taskInitCode = '';
if isfield(ir, 'modelInit')
    taskInitCode = gdsge.codegen.cxx.emitTask(gdsge.codegen.initView(ir), 'task_init');
end
cpp = fillTemplate(readTemplate('mex.tpl.cpp'), { ...
    'GDSGE_OTHER_INCLUDE',    ''; ...
    'TASK_INIT_CODE',         taskInitCode; ...
    'TASK_INF_HORIZON_CODE',  gdsge.codegen.cxx.emitTask(ir)});
```

`emitCompile.m` — replace the maxDim/extraDef computation:

```matlab
maxDim = sum(cellfun(@(p) p.length, pol)) + 4;
extraDef = '';
if isfield(ir, 'modelInit')
    initDim = sum(cellfun(@(p) p.length, ir.modelInit.variables.policyInit)) + 4;
    maxDim = max(maxDim, initDim);     % old: max(num_policy_total, num_policy_init_total)+4
    extraDef = '-DHAS_INIT';
end
```

and pass `'EXTRA_DEF', extraDef` in the fill list.

- [ ] **Step 4: Run the codegen suite — verify PASS** (new tests + the HL1996 C++
  snapshots, which must be byte-identical since HL1996 has no modelInit):

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','codegen')), runtests(fullfile('tests','HeatonLucas1996','codegen','tSnapshotHL1996.m')), runtests(fullfile('tests','HeatonLucas1996','codegen','tSnapshotCxxHL1996.m'))]; disp(table(results)); exit(any([results.Failed]))"
```

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen tests/+gdsgefix tests/codegen
git commit -m "feat(codegen-cxx): task_init from initView; -DHAS_INIT and max() MAXDIM"
```

---

### Task 10: Generated iter init segment (`emitIterInit`)

**Files:**
- Create: `src/+gdsge/+codegen/+mat/emitIterInit.m`
- Modify: `src/+gdsge/+codegen/+mat/emitIter.m`
- Test: `tests/codegen/tEmitIterInit.m` (new)

- [ ] **Step 1: Write the failing tests**

Create `tests/codegen/tEmitIterInit.m`:

```matlab
classdef tEmitIterInit < matlab.unittest.TestCase
    methods (Test)
        function emptyWithoutModelInit(tc)
            tc.verifyEmpty(gdsge.codegen.mat.emitIterInit(gdsgefix.minimalIR()));
        end
        function emitsGuardedInitSolve(tc)
            txt = gdsge.codegen.mat.emitIterInit(gdsgefix.minimalIRWithInit());
            tc.verifyTrue(contains(txt, 'if ~SkipModelInit'));
            tc.verifyTrue(contains(txt, 'GDSGE_CFG.taskName = MEX_TASK_INIT;'));
            tc.verifyTrue(contains(txt, 'GDSGE_CFG.splineVec = [];'));
            tc.verifyTrue(contains(txt, 'c0 = GDSGE_SOL(1:1,:);'));
            tc.verifyTrue(contains(txt, 'q0 = GDSGE_AUX(1:1,:);'));
        end
        function iterFileIsValidMatlab(tc)
            txt = gdsge.codegen.mat.emitIter(gdsgefix.minimalIRWithInit());
            tc.verifyTrue(contains(txt, 'if ~SkipModelInit'));
            t = mtree(txt);
            tc.verifyEqual(count(mtfind(t, 'Kind', 'ERR')), 0, ...
                'generated iter file does not parse');
        end
    end
end
```

- [ ] **Step 2: Run it — verify it fails** (`emitIterInit` undefined).

- [ ] **Step 3: Implement** — create `src/+gdsge/+codegen/+mat/emitIterInit.m`:

```matlab
function txt = emitIterInit(ir)
% EMITITERINIT  The model-init (last-period) solve segment of iter_<model>.m.
%   Returns '' when the model has no model_init. Old iter_init_template.m
%   parity: gated ONLY by SkipModelInit (runs even under WarmUp; the warm-up
%   interpolants supersede its output later), allocates the init problem over
%   the SAME state tensor, solves task_init, and leaves the init policy/aux
%   values in the workspace as rows (consumed by the interp initial exprs).
%   Runs BEFORE the main bounds/space section, which then re-allocates every
%   GDSGE_* array it touches.
if ~isfield(ir, 'modelInit'); txt = ''; return; end
iv     = gdsge.codegen.initView(ir);
bounds = gdsge.codegen.mat.emitBounds(iv);
pack   = gdsge.codegen.mat.emitDataPack(ir);   % layout identical to the main task
m      = ir.modelName;
nAux   = 0;
for i = 1:numel(iv.variables.aux); nAux = max(nAux, iv.variables.aux{i}.slot(2)); end
stateList     = strjoin(ir.states.names, ',');
stateNameCell = ['{''' strjoin(ir.states.names, ''',''') '''}'];
stateGridCell = ['{' stateList '}'];

w = gdsge.codegen.codeWriter();
w.add('%%%% Model-init solve (last-period problem; old semantics: SkipModelInit only)');
w.add('if ~SkipModelInit');
w.in();
w.addRaw(indent4(bounds.init));
w.add('GDSGE_EQVAL = 1e20*ones(%d,GDSGE_NPROB);', bounds.maxDim);
w.add('GDSGE_F = 1e20*ones(1,GDSGE_NPROB);');
w.add('GDSGE_SOL = zeros(%d,GDSGE_NPROB);', bounds.maxDim);
w.add('GDSGE_X0 = rand(size(GDSGE_SOL)) .* (GDSGE_UB-GDSGE_LB) + GDSGE_LB;');
w.add('GDSGE_SOL(:) = GDSGE_X0;');
w.add('GDSGE_AUX = zeros(%d,GDSGE_NPROB);', max(nAux, 1));
w.add('GDSGE_DATA = zeros(%d,GDSGE_NPROB);', pack.maxData);
w.add('%s', pack.iterPack);
w.add('GDSGE_CFG = struct();');
w.add('GDSGE_CFG.tolSol = TolSol; GDSGE_CFG.tolFun = TolFun;');
w.add('GDSGE_CFG.solMaxIter = SolMaxIter; GDSGE_CFG.numThreads = NumThreads;');
w.add('GDSGE_CFG.debugEvalOnly = GDSGE_DEBUG_EVAL_ONLY;');
w.add('GDSGE_CFG.useBroyden = UseBroyden; GDSGE_CFG.finiteDiffDelta = FiniteDiffDelta;');
w.add('GDSGE_CFG.useBroydenNow = 0;');
w.add('GDSGE_CFG.taskName = MEX_TASK_INIT;');
w.add('GDSGE_CFG.splineVec = [];');
w.add('GDSGE_CFG.ppNames = {}; GDSGE_CFG.ppCell = {};');
w.add('GDSGE_CFG.maxMinorIter = MaxMinorIter;');
w.add('GDSGE_CFG.probSize = GDSGE_SIZE;');
w.add('GDSGE_CFG.useNearestNeighbor = true;');
w.add('GDSGE_CFG.verboseRetry = ~NoPrint;');
w.add('GDSGE_CFG.adaptInSol = [];');
w.add('[GDSGE_SOL,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL,GDSGE_OPT_INFO,GDSGE_DIAG] = gdsge.runtime.solveProblems(@mex_%s, GDSGE_SOL, GDSGE_LB, GDSGE_UB, GDSGE_DATA, GDSGE_F, GDSGE_AUX, GDSGE_EQVAL, GDSGE_CFG); %%#ok<ASGLU>', m);
w.add('if any(GDSGE_DIAG.needResolved)');
w.add('    warning(''gdsge:runtime:unconverged'', ''model_init: %%s'', gdsge.runtime.reportUnconverged(GDSGE_DIAG.needResolved, GDSGE_F, GDSGE_SIZE, %s, %s, 5));', ...
    stateNameCell, stateGridCell);
w.add('end');
for i = 1:numel(iv.variables.policy)
    p = iv.variables.policy{i};
    w.add('%s = GDSGE_SOL(%d:%d,:);', p.name, p.slot(1), p.slot(2));
end
for i = 1:numel(iv.variables.aux)
    a = iv.variables.aux{i};
    w.add('%s = GDSGE_AUX(%d:%d,:);', a.name, a.slot(1), a.slot(2));
end
w.out();
w.add('end');
txt = w.str();
end

function s = indent4(s)
if isempty(s); return; end
lines = strsplit(s, newline, 'CollapseDelimiters', false);
for i = 1:numel(lines)
    if ~isempty(lines{i}); lines{i} = ['    ' lines{i}]; end
end
s = strjoin(lines, newline);
end
```

In `emitIter.m`, insert between the state-space-tensors section and the `%% Bounds`
section (i.e. after the `GDSGE_SIZE_STATE = ...` line and its `w.blank()`):

```matlab
initSeg = gdsge.codegen.mat.emitIterInit(ir);
if ~isempty(initSeg)
    w.addRaw(initSeg);
    w.blank();
end
```

- [ ] **Step 4: Run it — verify PASS**, plus the HL1996 snapshot suite (no modelInit →
  byte-identical):

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = [runtests(fullfile('tests','codegen','tEmitIterInit.m')), runtests(fullfile('tests','HeatonLucas1996','codegen','tSnapshotHL1996.m'))]; disp(table(results)); exit(any([results.Failed]))"
```

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/+mat tests/codegen/tEmitIterInit.m
git commit -m "feat(codegen-mat): iter model-init solve segment (SkipModelInit-gated)"
```

---

### Task 11: SIMU_INTERP simulate variant

**Files:**
- Create: `src/+gdsge/+codegen/+mat/emitSimulateInterp.m`
- Modify: `src/+gdsge/+codegen/generateMatlab.m`
- Test: `tests/codegen/tEmitSimulateInterp.m` (new)

- [ ] **Step 1: Write the failing tests**

Create `tests/codegen/tEmitSimulateInterp.m`:

```matlab
classdef tEmitSimulateInterp < matlab.unittest.TestCase
    methods (Test)
        function emitsInterpEvalNotSolver(tc)
            ir = gdsgefix.minimalIR();
            ir.options.simuInterp = 1; ir.options.simuResolve = 0;
            txt = gdsge.codegen.mat.emitSimulateInterp(ir);
            tc.verifyTrue(contains(txt, 'GDSGE_PP = IterRslt.output_interp;'));
            tc.verifyTrue(contains(txt, 'output_var_index'));
            tc.verifyFalse(contains(txt, 'solveProblems'));
            t = mtree(txt);
            tc.verifyEqual(count(mtfind(t, 'Kind', 'ERR')), 0);
        end
        function generateMatlabPicksVariant(tc)
            work = tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            ir = gdsgefix.minimalIR();
            ir.options.simuInterp = 1; ir.options.simuResolve = 0;
            files = gdsge.codegen.generateMatlab(ir, work);
            tc.verifyTrue(contains(fileread(files.simulateFile), 'output_interp'));
            ir2 = gdsgefix.minimalIR();
            files2 = gdsge.codegen.generateMatlab(ir2, work);
            tc.verifyTrue(contains(fileread(files2.simulateFile), 'solveProblems'));
        end
    end
end
```

- [ ] **Step 2: Run it — verify it fails.**

- [ ] **Step 3: Implement** — create `src/+gdsge/+codegen/+mat/emitSimulateInterp.m`:

```matlab
function txt = emitSimulateInterp(ir)
% EMITSIMULATEINTERP  simulate_<model>.m, SIMU_INTERP variant: no per-period
%   re-solve. var_output values come from IterRslt.output_interp evaluated at
%   the realized (shock, state) points each period; var_simu fields and state
%   transitions are then assigned exactly as in the resolve variant. Old
%   parity: code_template/simulate_interp_template.m + gdsge_parser.m:1684-1737.
simu = gdsge.codegen.mat.emitResultSimu(ir);
m = ir.modelName;
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
w.add('%% ---- pull the solved model from IterRslt');
for i = 1:numel(ir.states.names)
    w.add('%s = IterRslt.var_state.%s;', ir.states.names{i}, ir.states.names{i});
end
w.add('shock_trans = IterRslt.shock_trans;');
w.add('shock_num = IterRslt.shock_num;');
w.add('output_var_index = IterRslt.output_var_index;');
w.add('GDSGE_PP = IterRslt.output_interp;');
w.blank();
base = {'SimuPrintFreq','SimuSaveFreq','num_samples','num_periods', ...
    'NumThreads','EnforceSimuStateInbound','GEN_SHOCK_START_PERIOD','init'};
wl = [base, ir.states.names];
w.add('GDSGE_OPTIONS_VALID = {%s};', ['''' strjoin(wl, ''',''') '''']);
w.add('if nargin < 2; GDSGE_OPTIONS = struct(); end');
w.add('gdsge.runtime.unpackOptions(GDSGE_OPTIONS, GDSGE_OPTIONS_VALID);');
w.blank();
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
w.add('SimuRslt.shock(:,GEN_SHOCK_START_PERIOD:end) = gen_discrete_markov_rn(shock_trans,num_samples,length(GEN_SHOCK_START_PERIOD:num_periods+1),...');
w.add('    SimuRslt.shock(:,GEN_SHOCK_START_PERIOD));');
w.blank();
w.add('GDSGE_SHOCK_VAR_INDEX_BASE = ([0:num_samples-1]'')*shock_num;');
w.add('GDSGE_timer = tic;');
w.add('for GDSGE_t=1:num_periods');
w.in();
w.add('if EnforceSimuStateInbound==1');
for i = 1:numel(ir.states.names)
    s = ir.states.names{i};
    w.add('    SimuRslt.%s(:,GDSGE_t)=max(min(IterRslt.var_state.%s),SimuRslt.%s(:,GDSGE_t));', s, s, s);
    w.add('    SimuRslt.%s(:,GDSGE_t)=min(max(IterRslt.var_state.%s),SimuRslt.%s(:,GDSGE_t));', s, s, s);
end
w.add('end');
w.blank();
w.add('if shock_num>1');
w.add('    GDSGE_INTERP_RESULTS = myppual_mex(int32(NumThreads),GDSGE_PP.breaks,GDSGE_PP.coefs,...');
w.add('        int32(GDSGE_PP.pieces),int32(GDSGE_PP.order),int32(GDSGE_PP.dim),''not-a-knot'',[SimuRslt.shock(:,GDSGE_t)'';%s],[],[],[]);', simuStateRows);
w.add('else');
w.add('    GDSGE_INTERP_RESULTS = myppual_mex(int32(NumThreads),GDSGE_PP.breaks,GDSGE_PP.coefs,...');
w.add('        int32(GDSGE_PP.pieces),int32(GDSGE_PP.order),int32(GDSGE_PP.dim),''not-a-knot'',[%s],[],[],[]);', simuStateRows);
w.add('end');
for i = 1:numel(ir.variables.output)
    n = ir.variables.output{i};
    w.add('%s=GDSGE_INTERP_RESULTS(output_var_index.%s,:);', n, n);
end
w.blank();
w.add('GDSGE_SHOCK_VAR_LINEAR_INDEX = SimuRslt.shock(:,GDSGE_t+1) + GDSGE_SHOCK_VAR_INDEX_BASE;');
w.addRaw(indent4(simu.assign));
w.blank();
w.add('if mod(GDSGE_t,SimuPrintFreq)==0');
w.add('    fprintf(''Periods: %%d\\n'', GDSGE_t);');
w.add('    fprintf(''elapsed:%%.1fs\\n'', toc(GDSGE_timer));');
w.add('    GDSGE_timer = tic;');
w.add('end');
w.add('if mod(GDSGE_t,SimuSaveFreq)==0');
w.add('    save([''SimuRslt_%s_'' num2str(GDSGE_t) ''.mat''], ''SimuRslt'');', m);
w.add('end');
w.out();
w.add('end');
w.add('end');
txt = w.str();
end

function s = indent4(s)
if isempty(s); return; end
lines = strsplit(s, newline, 'CollapseDelimiters', false);
for i = 1:numel(lines)
    if ~isempty(lines{i}); lines{i} = ['    ' lines{i}]; end
end
s = strjoin(lines, newline);
end
```

In `generateMatlab.m`, replace the simulate line with:

```matlab
if isfield(ir.options, 'simuInterp') && ir.options.simuInterp == 1
    simTxt = gdsge.codegen.mat.emitSimulateInterp(ir);
else
    simTxt = gdsge.codegen.mat.emitSimulate(ir);
end
files.simulateFile = fullfile(outDir, ['simulate_' ir.modelName '.m']);
gdsge.codegen.writeText(files.simulateFile, simTxt);
```

Note for the implementer: `minimalIR`'s `var_simu` is `{'c'}` and `c` is its only output
— but `minimalIR` has no `output_var_index` consumer at unit-test level; the mtree check
plus the GLSW/Mendoza gates carry the correctness load.

- [ ] **Step 4: Run it — verify PASS** (plus `tSnapshotHL1996` unchanged: HL1996 has
  simuInterp=0).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen tests/codegen/tEmitSimulateInterp.m
git commit -m "feat(codegen-mat): SIMU_INTERP simulate variant from IterRslt.output_interp"
```

---

### Task 12: Barro end-to-end gate

**Files:**
- Create: `tests/Barro_et_al_2017/codegen/tEndToEndSafeAssets.m`

Known risk to clear first: the model body indexes a future variable with a constant
(`Re_n(2)` inside `x1 = shr_x1*(Rf/(Rf - Re_n(2)));`). HL1996 never did this. Check
`src/+gdsge/+codegen/+cxx/emitExpr.m` for the index/call-node-on-futureVar case; if it
mis-emits, FIRST add a failing case to the emitExpr unit tests in `tests/codegen/`
(pattern: expression AST with an index node whose base is a future var; expected C++
`Re_n_GDSGE_GRID[int(2)-1]`, consistent with the POPNARRAY accessor convention), then fix.

- [ ] **Step 1: Write the gate**

Create `tests/Barro_et_al_2017/codegen/tEndToEndSafeAssets.m`:

```matlab
classdef tEndToEndSafeAssets < matlab.unittest.TestCase
    % PHASE 7a GATE (Barro_et_al_2017): public API end-to-end vs goldens.
    % Slow: MEX compile + full policy iteration + 1000-period simulation.
    properties (Constant)
        RelTol = 1e-4;
        AbsTol = 1e-4;
    end
    methods (Test, TestTags = {'Slow'})
        function publicApiMatchesGolden(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'safe_assets.gmod'), work);

            gdsge_codegen('safe_assets');
            tc.assertTrue(exist(fullfile(work, ['mex_safe_assets.' mexext]), 'file') == 3, ...
                'MEX did not compile');

            opts = struct('SaveFreq', inf, 'NoSave', 1);
            IterRslt = iter_safe_assets(opts);

            golden = load(fullfile(modelDir, 'golden', 'IterRslt.mat'));
            G = golden.IterRslt;
            tc.verifyLessThan(IterRslt.Metric, 1e-6);
            tc.verifyLessThan(abs(IterRslt.Iter - G.Iter), 0.2*G.Iter + 20, ...
                sprintf('Iter=%d vs golden %d', IterRslt.Iter, G.Iter));
            r = gdsgetest.compareNumericClose(IterRslt.shock_trans, G.shock_trans, 1e-12, 1e-15);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_state, G.var_state, 1e-12, 1e-15);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_policy, G.var_policy, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_aux, G.var_aux, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_interp, G.var_interp, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));

            % ---- simulate (same reduced protocol as the capture script)
            rng(0823);
            SimuRslt = simulate_safe_assets(IterRslt, ...
                struct('num_samples', 6, 'num_periods', 1000));
            gs = load(fullfile(modelDir, 'golden', 'SimuRslt.mat'));
            GS = gs.SimuRslt;
            tc.verifyEqual(SimuRslt.shock, GS.shock, ...
                'shock paths differ — check shock_trans bit-identity first');
            flds = {'omega1','Rf','K1','b1','expectedRe'};
            T0 = 100;
            for i = 1:numel(flds)
                a = SimuRslt.(flds{i}); b = GS.(flds{i});
                tc.verifyEqual(size(a), size(b), flds{i});
                r = gdsgetest.compareNumericClose(a(:,1:T0), b(:,1:T0), 1e-2, 1e-2);
                tc.verifyTrue(r.pass, sprintf('%s early path: %s', flds{i}, ...
                    strjoin(r.failures, newline)));
                tc.verifyLessThan(abs(mean(a(:)) - mean(b(:))), 5e-3, flds{i});
                tc.verifyLessThan(abs(std(a(:)) - std(b(:))), 5e-3, flds{i});
            end
        end
    end
end
```

- [ ] **Step 2: Run the gate (slow)**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','Barro_et_al_2017','codegen','tEndToEndSafeAssets.m')); disp(table(results)); exit(any([results.Failed]))"
```

Expected: PASS. If the MEX fails to compile or numbers mismatch, debug with
superpowers:systematic-debugging; commit each real fix with its own failing-test-first
cycle. Do NOT loosen tolerances.

- [ ] **Step 3: Commit**

```bash
git add tests/Barro_et_al_2017/codegen
git commit -m "test(codegen): Barro end-to-end gate green - GDSGE_MIN + setup-replay bounds"
```

---

### Task 13: Mendoza end-to-end gate

**Files:**
- Create: `tests/Mendoza2010/codegen/tEndToEndMendoza2010.m`

This is the heavy gate: model_init (task_init + iter segment), 2-D states, adaptive
bounds, no-arg iter, fine-grid WarmUp re-solve, SIMU_INTERP simulate.

- [ ] **Step 1: Write the gate**

Create `tests/Mendoza2010/codegen/tEndToEndMendoza2010.m`:

```matlab
classdef tEndToEndMendoza2010 < matlab.unittest.TestCase
    % PHASE 7a GATE (Mendoza2010): public API end-to-end vs goldens, driven
    % exactly like the old test.m (no-arg iter, then 80x80 WarmUp re-solve),
    % then the reduced seeded SIMU_INTERP simulate. Slow.
    properties (Constant)
        RelTol = 1e-4;
        AbsTol = 1e-4;
    end
    methods (Test, TestTags = {'Slow'})
        function publicApiMatchesGolden(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'mendoza2010.gmod'), work);

            gdsge_codegen('mendoza2010');
            tc.assertTrue(exist(fullfile(work, ['mex_mendoza2010.' mexext]), 'file') == 3, ...
                'MEX did not compile');

            % ---- coarse solve, no-arg call (frozen public surface)
            IterRslt1 = iter_mendoza2010;
            g1 = load(fullfile(modelDir, 'golden', 'IterRslt1.mat'));
            tc.verifyLessThan(IterRslt1.Metric, 1e-6);
            tc.verifyLessThan(abs(IterRslt1.Iter - g1.IterRslt1.Iter), ...
                0.2*g1.IterRslt1.Iter + 20);

            % ---- fine-grid WarmUp re-solve (test.m parity)
            opts = struct;
            opts.cTilde = linspace(IterRslt1.var_state.cTilde(1), ...
                IterRslt1.var_state.cTilde(end), 80);
            opts.k = linspace(IterRslt1.var_state.k(1), IterRslt1.var_state.k(end), 80);
            opts.WarmUp = IterRslt1;
            IterRslt = iter_mendoza2010(opts);

            golden = load(fullfile(modelDir, 'golden', 'IterRslt.mat'));
            G = golden.IterRslt;
            tc.verifyLessThan(IterRslt.Metric, 1e-6);
            tc.verifyLessThan(abs(IterRslt.Iter - G.Iter), 0.2*G.Iter + 20, ...
                sprintf('Iter=%d vs golden %d', IterRslt.Iter, G.Iter));
            r = gdsgetest.compareNumericClose(IterRslt.shock_trans, G.shock_trans, 1e-12, 1e-15);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_state, G.var_state, 1e-12, 1e-15);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_policy, G.var_policy, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_aux, G.var_aux, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_interp, G.var_interp, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));

            % ---- SIMU_INTERP simulate (same protocol as capture)
            rng(0823);
            SimuRslt = simulate_mendoza2010(IterRslt, ...
                struct('num_samples', 6, 'num_periods', 1000));
            gs = load(fullfile(modelDir, 'golden', 'SimuRslt.mat'));
            GS = gs.SimuRslt;
            tc.verifyEqual(SimuRslt.shock, GS.shock, ...
                'shock paths differ — check shock_trans bit-identity first');
            flds = {'cTilde','k','c','Y','q','inv','b','mu','bNext','v','nx', ...
                'gdp','b_gdp','nx_gdp','lev','wkcptl'};
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
```

(No long-run moment check here: 1000 interp-simulated periods of a crisis model have
fat tails; the bit-exact shock path + early-path tolerance carry the load, mirroring
what the golden actually pins.)

- [ ] **Step 2: Run the gate (slow — two full solves + compile; budget per the capture
  timings recorded in Task 2)**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','Mendoza2010','codegen','tEndToEndMendoza2010.m')); disp(table(results)); exit(any([results.Failed]))"
```

Expected: PASS. Known-unknowns to debug here if red (failing-test-first for each):
- init-task data/slot mismatches (`task_init` POP vs MATLAB pack — both come from
  `dataLayout(ir)`, so a mismatch means the init segment packed the wrong arrays);
- interp `initial` exprs consuming init outputs (`flow/flow_scale` must be a defined
  1×NPROB row when `INTERP_INITIALIZE` runs — i.e. init segment emitted BEFORE the
  interp-initialization block);
- `applyWarmUp` interpolation onto the 80×80 tensor (2-D `evalPoints`);
- `GDSGE_INTERP_VEC'(cTildeNext',kNext)` — mixed primed/unprimed interp args in C++.

- [ ] **Step 3: Commit**

```bash
git add tests/Mendoza2010/codegen
git commit -m "test(codegen): Mendoza end-to-end gate green - model_init + SIMU_INTERP + 2-D states"
```

---

### Task 14: GLSW end-to-end gate

**Files:**
- Create: `tests/GLSW2020/codegen/tEndToEndGLSW.m`

- [ ] **Step 1: Write the gate**

Create `tests/GLSW2020/codegen/tEndToEndGLSW.m`:

```matlab
classdef tEndToEndGLSW < matlab.unittest.TestCase
    % PHASE 7a GATE (GLSW2020): public API end-to-end vs goldens: TolEq=1e-8,
    % INTERP_ORDER=4, model_init, SIMU_INTERP, simulate init/num overrides. Slow.
    properties (Constant)
        RelTol = 1e-4;
        AbsTol = 1e-4;
    end
    methods (Test, TestTags = {'Slow'})
        function publicApiMatchesGolden(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'GLSW_interp.gmod'), work);

            gdsge_codegen('GLSW_interp');
            tc.assertTrue(exist(fullfile(work, ['mex_GLSW_interp.' mexext]), 'file') == 3, ...
                'MEX did not compile');

            IterRslt = iter_GLSW_interp;

            golden = load(fullfile(modelDir, 'golden', 'IterRslt.mat'));
            G = golden.IterRslt;
            tc.verifyLessThan(IterRslt.Metric, 1e-8);
            tc.verifyLessThan(abs(IterRslt.Iter - G.Iter), 0.2*G.Iter + 20, ...
                sprintf('Iter=%d vs golden %d', IterRslt.Iter, G.Iter));
            r = gdsgetest.compareNumericClose(IterRslt.shock_trans, G.shock_trans, 1e-12, 1e-15);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_state, G.var_state, 1e-12, 1e-15);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_policy, G.var_policy, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_aux, G.var_aux, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_interp, G.var_interp, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));

            % ---- baseline SIMU_INTERP simulate
            rng(0823);
            SimuRslt = simulate_GLSW_interp(IterRslt, ...
                struct('num_samples', 6, 'num_periods', 1000));
            gs = load(fullfile(modelDir, 'golden', 'SimuRslt.mat'));
            GS = gs.SimuRslt;
            tc.verifyEqual(SimuRslt.shock, GS.shock);
            flds = {'a1','a2','P1','r','c1_shr'};
            T0 = 100;
            for i = 1:numel(flds)
                a = SimuRslt.(flds{i}); b = GS.(flds{i});
                tc.verifyEqual(size(a), size(b), flds{i});
                r = gdsgetest.compareNumericClose(a(:,1:T0), b(:,1:T0), 1e-2, 1e-2);
                tc.verifyTrue(r.pass, sprintf('%s early path: %s', flds{i}, ...
                    strjoin(r.failures, newline)));
            end

            % ---- init-override batch (pins options.init/num_samples/num_periods)
            irfOpts = struct;
            irfOpts.init = struct('shock', [ones(3,1); 2*ones(3,1)], 'a1', zeros(6,1));
            irfOpts.num_samples = 6;
            irfOpts.num_periods = 100;
            rng(0823);
            SimuRsltIrf = simulate_GLSW_interp(IterRslt, irfOpts);
            gi = load(fullfile(modelDir, 'golden', 'SimuRsltIrf.mat'));
            GI = gi.SimuRsltIrf;
            tc.verifyEqual(SimuRsltIrf.shock, GI.shock);
            tc.verifyEqual(SimuRsltIrf.shock(:,1), [1;1;1;2;2;2]);
            for i = 1:numel(flds)
                a = SimuRsltIrf.(flds{i}); b = GI.(flds{i});
                tc.verifyEqual(size(a), size(b), flds{i});
                r = gdsgetest.compareNumericClose(a, b, 1e-2, 1e-2);
                tc.verifyTrue(r.pass, sprintf('%s irf path: %s', flds{i}, ...
                    strjoin(r.failures, newline)));
            end
        end
    end
end
```

- [ ] **Step 2: Run the gate (slow)** — same command pattern as Task 13 with the GLSW
  path. Expected: PASS. Most machinery is shared with Mendoza; new exposure is
  TolEq=1e-8 convergence depth and the init-override path through
  `emitResultSimu.initOverwrite` (which writes `SimuRslt.shock(:,1:...)` BEFORE the
  shock-path generation that conditions on column 1 — order matters and matches the
  emitter layout from Task 11).

- [ ] **Step 3: Commit**

```bash
git add tests/GLSW2020/codegen
git commit -m "test(codegen): GLSW end-to-end gate green - TolEq 1e-8 + simulate overrides"
```

---

### Task 15: Full suite, PROGRESS.md, finish the branch

**Files:**
- Modify: `PROGRESS.md`

- [ ] **Step 1: Run the full suite**

```
matlab -batch "cd('tests'); run_tests"
```

Expected: exit 0; `tests/results/junit.xml` shows 0 failures. Budget: the previous
suite (~15-30 min with HL1996 Slow gates) plus the three new gates (use the Task 2
capture timings as the estimate).

- [ ] **Step 2: Update PROGRESS.md** — check off Phase 7a in the Phases list, replacing
  its line with the done-form (fill in the captured Iter values):

```markdown
- ☑ **Phase 7a — Spline-path widening** (done <date>)
  - ☑ Goldens captured (old toolbox): safe_assets, mendoza2010 (incl. 80×80 WarmUp
    re-solve), GLSW_interp (incl. init-override batch); reduced seeded simulates
  - ☑ IR: `setup` section (gmod declaration code replayed verbatim — old
    setParamsCode parity); transition `primed` flag; init-bound item
  - ☑ Parser: model_init + var_policy_init/var_aux_init/inbound_init; var_simu →
    var_output → var_aux promotion (old parity); SimuPrintFreq/SimuSaveFreq
  - ☑ Codegen: task_init via initView (emitter reuse), -DHAS_INIT + max() MAXDIM,
    iter init segment (SkipModelInit-gated, runs under WarmUp), SIMU_INTERP
    simulate variant
  - ☑ Gates green: tEndToEndSafeAssets, tEndToEndMendoza2010, tEndToEndGLSW
    (+ HL1996 gates stay green)
- ☐ **Phase 7b — ASG widening** (CaoKS2016, Bianchi2011_asg; ASG architecture,
  var_tensor, asg_*_struct result fields, SkipModelInit/shock-override options) — NEXT
```

and add a changelog entry at the top of the Changelog section summarizing the same
(mention the two latent-bug fixes: setup replay, unprimed-transition indexing).

- [ ] **Step 3: Commit**

```bash
git add PROGRESS.md
git commit -m "docs(progress): Phase 7a complete - three spline models green vs goldens"
```

- [ ] **Step 4: Finish the branch** — use the superpowers:finishing-a-development-branch
  skill (prior phases were merged into `main`).

---

## Self-review notes (already applied)

- **Spec coverage:** §6 capture protocol → Task 2 (full-convergence iter incl. Mendoza
  two-stage; ≤10×≤1000 seeded simulates; integrity tests; timings logged with the
  stop-and-escalate clause). §7 front-end → Tasks 3 (schema/setup/primed), 5 (options),
  6 (promotion — discovered requirement, old-parser parity), 7 (model_init), 8 (gates;
  "no hand-authored reference IRs" honored — gates assert structure + round-trip).
  §8.1 Barro → Task 12; §8.2 Mendoza → Tasks 9/10/11/13 (init task, iter segment with
  old SkipModelInit-only gating, SIMU_INTERP, 2-D states, no-arg iter — verified already
  supported, warm-up finer grids); §8.3 GLSW → Tasks 5/11/14 (INTERP_ORDER=4 verified
  already-default; overrides). §9 testing → unit tests in each task + gates Slow-tagged +
  HL1996 kept green at Tasks 4/9/10/11. §10 risks → Task 2 timing clause; init-task
  layout via shared dataLayout; per-path goldens (interp vs resolve) by construction.
- **Beyond-spec discoveries built in:** setup-text replay (Barro `Re_n`, GLSW `Ngrid`
  would be undefined otherwise); transition `primed` flag (Mendoza `k' = kNext`, GLSW
  `a1' = a1n` would index out of bounds); var_simu→output→aux promotion (Mendoza/Barro
  SimuRslt fields and old IterRslt.var_aux parity).
- **Type consistency:** `emitTask(ir, taskName)` matches both call sites (Task 9);
  `initView` field mapping matches schema `modelInit` (Task 3 initBoundItem name is
  fText, consumed by `emitBounds` which only reads name/lower/upper/adaptiveFactor);
  `emitIterInit` consumes `bounds.maxDim/pack.maxData` exactly as `emitIter` does;
  gates' simulate options literals match the capture scripts verbatim (6×1000;
  GLSW irf 6×100 with the same init struct).
- **Placeholder scan:** the only deliberately-unpinned numbers are golden Iter counts
  (unknowable before capture) — gates derive their windows from the loaded golden
  (`0.2*G.Iter + 20`) instead of hardcoding.

