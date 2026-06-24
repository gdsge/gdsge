# Cartesian SIMU_INTERP in-MEX simulation + myppual retirement — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make cartesian `SIMU_INTERP=1` simulation run its whole period loop in a generated C++ MEX (hans-style, in place), and delete `myppual` (`myppual.m`, `myppual_mex.mexw64`, `convert_to_interp_eval_array.m`, `get_scalar.m`) entirely.

**Architecture:** Rebuild `IterRslt.output_interp` and the SIMU_RESOLVE/warmup interpolants in the **stacked uniform-order** layout (shock = stacked vector index, uniform cubic spline over states only), built by the existing `interp_construct_mex` and evaluated by the existing `MatlabPp`/`InterpEval` headers. Add a generic MATLAB-callable evaluator `interp_eval_mex` for the resolve warm-start, `applyWarmUp`, and the SIMU_INTERP fallback; add a codegen'd whole-simulation MEX `simulate_<model>_mex` for the fast path.

**Tech Stack:** MATLAB R2025b, C++ MEX (`include/MatlabPp.h`, `include/InterpEval.h`, `include/interp_lite.h`, OpenMP), `matlab.unittest`.

**Spec:** `docs/superpowers/specs/2026-06-15-cartesian-simu-interp-mex-design.md`

---

## Orientation (read before starting)

- **Golden rule:** never persist a MATLAB path; each MATLAB run `addpath`s one source. Tests are driven by `matlab -batch "cd('tests'); run_tests"` (PowerShell 7 is **not** installed — do not use `tests/run.ps1`). Targeted run: `matlab -batch "cd('tests'); runtests('<relpath>')"`.
- **Run MATLAB sequentially** — never two MATLAB processes at once (OpenMP saturates all cores).
- **Pin `rng`** before any `iter_*` that feeds an end-to-end gate (safe_assets-style root-picking).
- Key existing files to mirror:
  - `src/+gdsge/+codegen/ensureSplineConstructMex.m` — kernel compile + `.cache` pattern.
  - `src/kernels/interp_construct_mex.cpp` — fused constructor; its `ppCell{i}` has fields `breaks/coefs/pieces/order/dim` and `dim = numArray` (the number of stacked vectors), exactly what `MatlabPp` reads.
  - `include/MatlabPp.h` — `MatlabPp<xDim,interpOrder>(pp)`, `search(double*,int*,double*)`, `nosearch_eval_2/4(xSiteToLeft,cellOfSite,vecIdx)`.
  - `src/+gdsge/+runtime/constructSplines.m` — calls `interp_construct_mex(breaks, valuesStruct, int32(orderVec), int32(extrapVec), int32(numThreads))`.
  - `src/+gdsge/+codegen/+mat/emitResultIter.m:40` — builds `output_interp` (to change).
  - `src/+gdsge/+codegen/+mat/emitSimulateInterp.m` — SIMU_INTERP simulate (to rewrite).
  - `src/+gdsge/+codegen/+mat/emitSimulate.m:54-83,140-147` — resolve warm-start construct + eval (to change).
  - `src/+gdsge/+runtime/applyWarmUp.m:33-35,51` — warmup construct + eval (to change).
  - `src/+gdsge/+codegen/generateCxx.m` / `codegen.m` — C++ emission + compile flow (to extend for the simulate MEX).

## File structure

**Create**
- `src/kernels/interp_eval_mex.cpp` — generic uniform-order evaluator MEX.
- `src/+gdsge/+codegen/ensureInterpEvalMex.m` — compile + cache for `interp_eval_mex`.
- `src/+gdsge/+codegen/+cxx/emitSimulateMex.m` — emits `simulate_<model>_mex.cpp` (codegen'd whole-sim MEX).
- `templates/simulate_mex.tpl.cpp` — the C++ template for the whole-sim MEX.
- Tests: `tests/kernels/tInterpEvalMex.m`, `tests/codegen/tEmitSimulateMex.m`, `tests/GLSW2020/codegen/tSimuInterpMexGLSW.m`, `tests/Mendoza2010/codegen/tSimuInterpMexMendoza.m`.

**Modify**
- `src/+gdsge/+codegen/+mat/emitResultIter.m` — stacked-uniform `output_interp` via `interp_construct_mex`; drop `get_scalar`.
- `src/+gdsge/+codegen/+mat/emitSimulateInterp.m` — rewrite (Phase 2 per-period; Phase 3 whole-sim).
- `src/+gdsge/+codegen/+mat/emitSimulate.m` — resolve warm-start construct/eval onto `interp_construct_mex`/`interp_eval_mex`.
- `src/+gdsge/+runtime/applyWarmUp.m` — construct/eval onto the new kernels.
- `src/+gdsge/+runtime/constructSplines.m` — drop legacy `myppual` branch.
- `src/+gdsge/+runtime/useFusedConstruct.m` — delete (toggle gone).
- `src/+gdsge/+codegen/generateCxx.m` — emit `simulate_<model>_mex.cpp` + add to `compile_<model>.m`; `ensureInterpEvalMex` for spline models.
- `templates/compile.tpl.m` — add a second `mex` call for `simulate_MODEL_NAME_mex.cpp` (spline + SIMU_INTERP only).
- `src/+gdsge/+runtime/Contents.m`, `src/+gdsge/+runtime/ensurePath.m` — drop `myppual` mentions.
- Existing codegen snapshot tests asserting the old `myppual(output_interp)` / `GDSGE_PP = IterRslt.output_interp` strings: `tests/codegen/tEmitResultIter.m:45`, `tests/codegen/tEmitSimulateInterp.m`.

**Delete**
- `src/kernels/myppual.m`, `src/kernels/myppual_mex.mexw64`, `src/kernels/convert_to_interp_eval_array.m`, `src/kernels/get_scalar.m`.

---

# Phase 1 — Generic uniform-order evaluator `interp_eval_mex`

Independent, no behavior change. Produces the evaluator that lets Phase 2 delete `myppual_mex`.

### Task 1.1: `interp_eval_mex` kernel + ensure/compile

**Files:**
- Create: `src/kernels/interp_eval_mex.cpp`
- Create: `src/+gdsge/+codegen/ensureInterpEvalMex.m`
- Test: `tests/kernels/tInterpEvalMex.m`

- [ ] **Step 1: Write the failing test**

`tests/kernels/tInterpEvalMex.m`:

```matlab
classdef tInterpEvalMex < matlab.unittest.TestCase
    % interp_eval_mex evaluates a stacked uniform-order pp (states-only spline,
    % shock = vector index) and must agree with a direct MatlabPp-equivalent
    % evaluation (here cross-checked against interp1 on a 1-D cubic build).
    methods (TestClassSetup)
        function addPath(tc)
            here = fileparts(mfilename('fullpath'));
            src  = fullfile(fileparts(fileparts(here)), 'src');
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(src));
            gdsge.runtime.ensurePath();
            gdsge.codegen.ensureSplineConstructMex();
            gdsge.codegen.ensureInterpEvalMex();
        end
    end
    methods (Test)
        function evalMatchesConstructValuesAtNodes(tc)
            % 1 state, 3 grid points (>=3: myppual 2-point bug), 2 stacked vectors.
            x = [0 1 2];
            v = [ 10 11 12 ;     % vector 1
                  20 19 18 ];    % vector 2   -> values [numArray=2, nx=3]
            [ppCell, ~] = interp_construct_mex({x}, struct('GDSGE_V1', v), ...
                int32(4), int32([]), int32(1));
            pp = ppCell{1};
            % evaluate both vectors at the grid nodes -> must reproduce v
            sites = [0 1 2];                      % 1 x 3 sites
            idx   = int32([0;1]);                 % both stacked vectors (0-based)
            vals  = interp_eval_mex(int32(1), pp, sites, idx);   % 2 x 3
            tc.verifyEqual(vals, v, 'AbsTol', 1e-9);
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "cd('tests'); runtests('kernels/tInterpEvalMex.m')"`
Expected: FAIL — `Undefined function 'interp_eval_mex'` (and `ensureInterpEvalMex`).

- [ ] **Step 3: Write `interp_eval_mex.cpp`**

`src/kernels/interp_eval_mex.cpp`:

```cpp
// interp_eval_mex.cpp — generic uniform-order evaluator for stacked pp structs.
//
//   vals = interp_eval_mex(int32 numThreads, pp, sites, int32 idx)
//
//   pp     : struct with breaks (1xD cell), coefs, order (int32 1xD, uniform),
//            dim (int32 numArray) — as produced by interp_construct_mex.
//   sites  : D x N matrix of evaluation points (column = one site).
//   idx    : K x 1 or K x N int32, 0-based stacked-vector indices to evaluate.
//            K x 1  -> same K vectors at every site;
//            K x N  -> per-site vector indices (one column per site).
//   vals   : K x N, vals(k,n) = pp vector idx(k[,n]) evaluated at sites(:,n).
//
// Uniform order only (every breaks dim shares order[0]); supports order 2 or 4.
// Evaluation core is include/MatlabPp.h (the same evaluator the model MEX uses).
#include "mex.h"
#include "MatlabPp.h"
#ifdef USE_OMP
#include <omp.h>
#endif

#define MAXDIM 10

template <int XDIM, int ORDER>
static void run(const mxArray* pp, const double* sites, int N,
                const int* idx, int K, bool perSite, double* out, int nThreads)
{
    MatlabPp<XDIM, ORDER> ev(pp);
#ifdef USE_OMP
    omp_set_num_threads(nThreads);
    #pragma omp parallel for
#endif
    for (int n = 0; n < N; ++n)
    {
        double xSite[MAXDIM], xLeft[MAXDIM];
        int cell[MAXDIM];
        for (int d = 0; d < XDIM; ++d) xSite[d] = sites[(size_t)n*XDIM + d];
        ev.search(xSite, cell, xLeft);
        const int* col = perSite ? idx + (size_t)n*K : idx;
        for (int k = 0; k < K; ++k)
        {
            double r = (ORDER == 4)
                ? ev.nosearch_eval_4(xLeft, cell, col[k])
                : ev.nosearch_eval_2(xLeft, cell, col[k]);
            out[(size_t)n*K + k] = r;
        }
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs != 4) mexErrMsgIdAndTxt("gdsge:interpEval:nargin",
        "interp_eval_mex(numThreads, pp, sites, idx)");
    int nThreads  = *((int*)mxGetData(prhs[0]));
    const mxArray* pp = prhs[1];
    const double* sites = mxGetPr(prhs[2]);
    int XDIM = (int)mxGetM(prhs[2]);
    int N    = (int)mxGetN(prhs[2]);
    const int* idx = (const int*)mxGetData(prhs[3]);
    int K = (int)mxGetM(prhs[3]);
    bool perSite = (mxGetN(prhs[3]) == (mwSize)N) && N > 1;

    int order0 = *((int*)mxGetData(mxGetField(pp, 0, "order")));

    plhs[0] = mxCreateDoubleMatrix(K, N, mxREAL);
    double* out = mxGetPr(plhs[0]);

    #define DISPATCH(D) \
        if (XDIM == D) { if (order0 == 4) run<D,4>(pp,sites,N,idx,K,perSite,out,nThreads); \
                         else              run<D,2>(pp,sites,N,idx,K,perSite,out,nThreads); return; }
    DISPATCH(1) DISPATCH(2) DISPATCH(3) DISPATCH(4) DISPATCH(5)
    #undef DISPATCH
    mexErrMsgIdAndTxt("gdsge:interpEval:dim", "unsupported state dim %d", XDIM);
}
```

- [ ] **Step 4: Write `ensureInterpEvalMex.m`** (mirror `ensureSplineConstructMex.m`)

`src/+gdsge/+codegen/ensureInterpEvalMex.m`:

```matlab
function ensureInterpEvalMex()
% ENSUREINTERPEVALMEX  Compile src/kernels/interp_eval_mex.cpp when needed
%   (hash cache interp_eval_mex.cache beside it, same mechanism as
%   ensureSplineConstructMex). Generic uniform-order evaluator over MatlabPp,
%   the myppual_mex replacement for MATLAB-side interpolation. The cache key
%   includes MatlabPp.h / InterpEval.h / interp_lite.h so header edits rebuild.
here      = fileparts(mfilename('fullpath'));   % src/+gdsge/+codegen
srcRoot   = fileparts(fileparts(here));         % src
repoRoot  = fileparts(srcRoot);
kernels   = fullfile(srcRoot, 'kernels');
includeDir = fullfile(repoRoot, 'include');
cppFile   = fullfile(kernels, 'interp_eval_mex.cpp');
cacheFile = fullfile(kernels, 'interp_eval_mex.cache');
mexFile   = fullfile(kernels, ['interp_eval_mex.' mexext]);
cppText   = fileread(cppFile);
for h = {'MatlabPp.h','InterpEval.h','interp_lite.h'}
    cppText = [cppText, fileread(fullfile(includeDir, h{1}))]; %#ok<AGROW>
end
if exist(mexFile, 'file') == 3 && ~gdsge.codegen.needsCompile(cppText, cacheFile)
    return;
end
fprintf('Compiling interp_eval_mex (cache-gated):\n');
oldCd = pwd; restore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(kernels);
mex(['-I' includeDir], '-DUSE_OMP', 'interp_eval_mex.cpp', ...
    'OPTIMFLAGS=/O2 /DNDEBUG', ...
    'COMPFLAGS=$COMPFLAGS /wd4267 /wd4068 /wd4091 /openmp');
gdsge.codegen.writeText(cacheFile, cppText);
end
```

- [ ] **Step 5: Run test to verify it passes**

Run: `matlab -batch "cd('tests'); runtests('kernels/tInterpEvalMex.m')"`
Expected: PASS (1 test).

- [ ] **Step 6: Commit**

```bash
git add src/kernels/interp_eval_mex.cpp src/+gdsge/+codegen/ensureInterpEvalMex.m tests/kernels/tInterpEvalMex.m
git commit -m "feat(kernels): interp_eval_mex generic uniform-order evaluator"
```

---

# Phase 2 — Stacked uniform-order interpolants; retire myppual (per-period eval)

This phase changes `output_interp`/resolve/warmup to the stacked uniform-order layout, routes all MATLAB-side construction/eval through `interp_construct_mex`/`interp_eval_mex`, and **deletes `myppual`**. SIMU_INTERP simulation still loops periods in MATLAB (now calling `interp_eval_mex`); Phase 3 makes it whole-sim. Ending all goldens green here proves the layout change is correct before adding the MEX loop.

### Task 2.1: Stacked-uniform `output_interp` in `emitResultIter`; drop `get_scalar`

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitResultIter.m:24-40,59-71`
- Modify: `tests/codegen/tEmitResultIter.m:45`

- [ ] **Step 1: Update the snapshot test first (it pins the old string)**

In `tests/codegen/tEmitResultIter.m`, replace the assertion at line 45:

```matlab
% OLD: tc.verifyTrue(contains(txt, 'IterRslt.output_interp=myppual(output_interp);'));
tc.verifyTrue(contains(txt, 'interp_construct_mex'));
tc.verifyTrue(contains(txt, 'IterRslt.output_interp ='));
tc.verifyFalse(contains(txt, 'myppual'));
tc.verifyFalse(contains(txt, 'get_scalar'));
```

- [ ] **Step 2: Run to verify it fails**

Run: `matlab -batch "cd('tests'); runtests('codegen/tEmitResultIter.m')"`
Expected: FAIL — generated text still contains `myppual(output_interp)` / `get_scalar`.

- [ ] **Step 3: Rewrite the output-interp + scalar emission**

In `src/+gdsge/+codegen/+mat/emitResultIter.m`, replace the block at lines 24-40 (the `outputVarStack`/`output_interp=struct...`/`myppual` lines) with the stacked-uniform build. `nComp = sum(output var lengths)`; values stacked `[shock_num*nComp, stateDims]` with array index `(shock-1)*nComp + comp` (comp fastest):

```matlab
w.add('outputVarStack = cat(1,%s);', strjoin(parts, ','));    % [nComp, shock_num*prod(state)]
w.add('GDSGE_NCOMP = size(outputVarStack,1);');
w.add('if shock_num>1');
w.add('    outputVarStack = reshape(outputVarStack, GDSGE_NCOMP, shock_num, GDSGE_SIZE_STATE{:});');
% collapse (comp,shock) -> stacked array dim, comp fastest == (shock-1)*nComp+comp
w.add('    GDSGE_OUTPUT_VALUES = reshape(outputVarStack, GDSGE_NCOMP*shock_num, GDSGE_SIZE_STATE{:});');
w.add('else');
w.add('    GDSGE_OUTPUT_VALUES = reshape(outputVarStack, GDSGE_NCOMP, GDSGE_SIZE_STATE{:});');
w.add('end');
w.add('[GDSGE_OUTPUT_PP_CELL, ~] = interp_construct_mex({%s}, ...', stateList);
w.add('    struct(''GDSGE_V1'', GDSGE_OUTPUT_VALUES), ...');
w.add('    int32(OutputInterpOrder*ones(1,length(%s))), int32([]), int32(NumThreads));', gridCell);
w.add('IterRslt.output_interp = GDSGE_OUTPUT_PP_CELL{1};');
w.add('IterRslt.output_interp.GDSGE_NCOMP = GDSGE_NCOMP;');   % needed by simulate
```

Then replace the two `get_scalar(...)` lines (currently 64 and 71) with plain element extraction (the only thing `get_scalar` did when the grid is `{1,0}`; for the real grid it returned the struct unchanged, so just assign the struct):

```matlab
% line ~64
w.add('IterRslt.var_policy = GDSGE_VAR_POLICY;');
% line ~71
w.add('IterRslt.var_aux = GDSGE_VAR_AUX;');
```

> Note: `get_scalar` only collapsed to element (1) for the degenerate `state_vector={1,0}` shape, which never occurs for a real model grid (it always returned the struct unchanged). Assigning the struct directly is equivalent for every supported model.

- [ ] **Step 4: Run to verify the snapshot test passes**

Run: `matlab -batch "cd('tests'); runtests('codegen/tEmitResultIter.m')"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/+mat/emitResultIter.m tests/codegen/tEmitResultIter.m
git commit -m "feat(codegen): stacked uniform-order output_interp; drop get_scalar"
```

### Task 2.2: Per-period `interp_eval_mex` SIMU_INTERP loop (stacked layout)

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitSimulateInterp.m:34,71-85`
- Modify: `tests/codegen/tEmitSimulateInterp.m` (strings referencing `myppual_mex` / `output_interp` access)

- [ ] **Step 1: Update the codegen-string test first**

In `tests/codegen/tEmitSimulateInterp.m`, replace assertions that pin `GDSGE_PP = IterRslt.output_interp;` / `myppual_mex` with the new shape (spline branch):

```matlab
tc.verifyTrue(contains(txt, 'interp_eval_mex'));
tc.verifyFalse(contains(txt, 'myppual'));
```

(Leave the ASG-branch assertions unchanged — ASG still uses `asg.eval_vec`.)

- [ ] **Step 2: Run to verify it fails**

Run: `matlab -batch "cd('tests'); runtests('codegen/tEmitSimulateInterp.m')"`
Expected: FAIL — still emits `myppual_mex`.

- [ ] **Step 3: Rewrite the spline branch of `emitSimulateInterp.m`**

Replace the spline (`~isAsg`) construct line (34) and the eval block (74-85) so the per-period eval uses `interp_eval_mex` on the stacked layout. The site is **states only**; the per-sample idx selects the realized shock's component block:

```matlab
% pull from IterRslt (spline branch)
w.add('GDSGE_PP = IterRslt.output_interp;');
w.add('GDSGE_NCOMP = IterRslt.output_interp.GDSGE_NCOMP;');
...
% inside the period loop, spline branch (replaces lines 74-85):
w.add('GDSGE_SITES = [%s];', simuStateRows);                          % nStates x num_samples
w.add('GDSGE_IDX = int32((SimuRslt.shock(:,GDSGE_t)''-1)*GDSGE_NCOMP + (0:GDSGE_NCOMP-1)'');'); % nComp x num_samples
w.add('GDSGE_INTERP_RESULTS = interp_eval_mex(int32(NumThreads), GDSGE_PP, GDSGE_SITES, GDSGE_IDX);');
```

The downstream `output_var_index` extraction (lines 82-85) and `simu.assign` (the primed/unprimed transitions) are unchanged — `GDSGE_INTERP_RESULTS` has the same `[nComp x num_samples]` shape as before, indexed by `output_var_index`. The `EnforceSimuStateInbound` clamp (63-69) stays.

- [ ] **Step 4: Run the codegen-string test**

Run: `matlab -batch "cd('tests'); runtests('codegen/tEmitSimulateInterp.m')"`
Expected: PASS.

- [ ] **Step 5: End-to-end gate — GLSW (unprimed) + Mendoza (primed)**

Run: `matlab -batch "cd('tests'); runtests('GLSW2020/codegen/tEndToEndGLSW.m')"`
Run: `matlab -batch "cd('tests'); runtests('Mendoza2010/codegen/tEndToEndMendoza2010.m')"`
Expected: PASS — `SimuRslt` matches goldens within tolerance (interpolants identical at integer shocks). If FAIL, debug the stacked idx ordering (`(shock-1)*nComp + comp`, comp fastest) before proceeding.

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+codegen/+mat/emitSimulateInterp.m tests/codegen/tEmitSimulateInterp.m
git commit -m "feat(codegen): SIMU_INTERP per-period eval via interp_eval_mex (stacked)"
```

### Task 2.3: Resolve warm-start + applyWarmUp onto the new kernels

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitSimulate.m:53-83,140-147`
- Modify: `src/+gdsge/+runtime/applyWarmUp.m:30-51`

- [ ] **Step 1: Rewrite the resolve warm-start construct + eval**

In `emitSimulate.m`, replace the `GDSGE_PP=struct(...);GDSGE_PP=myppual(...)` block (53-83) with a stacked-uniform build via `interp_construct_mex` of `IterRslt.GDSGE_PROB.GDSGE_SOL` (shape `[nSolComp, shock_num, stateDims]`), and replace the per-period `myppual_mex` eval (140-147) with `interp_eval_mex`:

```matlab
% construct (replaces 53-83): GDSGE_SOL stacked by shock as vector index, uniform order
w.add('GDSGE_NSOLCOMP = size(IterRslt.GDSGE_PROB.GDSGE_SOL,1);');
w.add('if shock_num>1');
w.add('    GDSGE_SOL_VALUES = reshape(IterRslt.GDSGE_PROB.GDSGE_SOL, GDSGE_NSOLCOMP*shock_num, %s);', ...
    'IterRslt.GDSGE_PROB.GDSGE_SIZE(2:end)');
w.add('else');
w.add('    GDSGE_SOL_VALUES = reshape(IterRslt.GDSGE_PROB.GDSGE_SOL, GDSGE_NSOLCOMP, %s);', ...
    'IterRslt.GDSGE_PROB.GDSGE_SIZE(2:end)');
w.add('end');
w.add('[GDSGE_PP_CELL, ~] = interp_construct_mex({%s}, struct(''GDSGE_V1'', GDSGE_SOL_VALUES), ...', stateList);
w.add('    int32(OutputInterpOrder*ones(1,length(%s))), int32([]), int32(NumThreads));', stateGridCell);
w.add('GDSGE_PP = GDSGE_PP_CELL{1};');

% eval (replaces 140-147):
w.add('GDSGE_IDX = int32((SimuRslt.shock(:,GDSGE_t)''-1)*GDSGE_NSOLCOMP + (0:GDSGE_NSOLCOMP-1)'');');
w.add('GDSGE_SOL = interp_eval_mex(int32(NumThreads), GDSGE_PP, [%s], GDSGE_IDX);', simuStateRows);
```

> The macOS-only `'pp'`/`myppual(myppual())` branch (68-83) is removed — the stacked-uniform path is platform-independent.

- [ ] **Step 2: Rewrite `applyWarmUp.m` construct + eval**

In `src/+gdsge/+runtime/applyWarmUp.m`, replace `myppual(...)` construction (line 51) and the three `myppual(solI/lbI/ubI, evalPoints)` evals (33-35) with `interp_construct_mex` + `interp_eval_mex`. The warmup interpolant is state-only uniform-order; `evalPoints` are the fine-grid state sites:

```matlab
% construct (replaces line 51): one pp over states, numArray = number of stacked rows
[ppCellWU, ~] = interp_construct_mex(breaks, struct('GDSGE_V1', values), ...
    int32(order), int32(extrap), int32(numThreads));
pp = ppCellWU{1};
% eval (replaces 33-35): all stacked vectors at each fine-grid point
nVec = pp.dim;
idx  = int32((0:double(nVec)-1)');
sol   = reshape(interp_eval_mex(int32(numThreads), ppSol, spec.evalPoints, idxSol), size(sol));
lbNew = reshape(interp_eval_mex(int32(numThreads), ppLb,  spec.evalPoints, idxLb ), size(lb));
ubNew = reshape(interp_eval_mex(int32(numThreads), ppUb,  spec.evalPoints, idxUb ), size(ub));
```

> Read the current `applyWarmUp.m` in full first and adapt variable names (`solI/lbI/ubI`, `values`, `breaks`, `order`, `extrap`) to its actual locals; the transformation is mechanical (construct→`interp_construct_mex`, eval→`interp_eval_mex`).

- [ ] **Step 3: Ensure kernels are built where these run**

`emitSimulate`/`applyWarmUp` run inside generated `simulate_*`/`iter_*`. Add `gdsge.codegen.ensureInterpEvalMex()` next to the existing `ensureSplineConstructMex()` call in `generateCxx.m` (the spline branch — see Task 3.4) so the evaluator exists before any generated file runs.

- [ ] **Step 4: End-to-end gates (resolve + warmup paths)**

Run: `matlab -batch "cd('tests'); runtests('HeatonLucas1996/codegen/tEndToEndHL1996.m')"`  (resolve + warmup)
Run: `matlab -batch "cd('tests'); runtests('Mendoza2010/codegen/tEndToEndMendoza2010.m')"`  (warmup re-solve)
Run: `matlab -batch "cd('tests'); runtests('Barro_et_al_2017/codegen/tEndToEndSafeAssets.m')"`  (pin rng)
Expected: PASS within tolerance. (Resolve warm-start only seeds the solver; small interp differences converge to the same solution.)

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/+mat/emitSimulate.m src/+gdsge/+runtime/applyWarmUp.m src/+gdsge/+codegen/generateCxx.m
git commit -m "feat: resolve warm-start + applyWarmUp onto interp_construct_mex/interp_eval_mex"
```

### Task 2.4: Delete myppual and prune the legacy construct path

**Files:**
- Modify: `src/+gdsge/+runtime/constructSplines.m` (drop legacy branch)
- Delete: `src/+gdsge/+runtime/useFusedConstruct.m`
- Delete: `src/kernels/myppual.m`, `src/kernels/myppual_mex.mexw64`, `src/kernels/convert_to_interp_eval_array.m`, `src/kernels/get_scalar.m`
- Modify: `src/+gdsge/+runtime/Contents.m`, `src/+gdsge/+runtime/ensurePath.m` (comment text)
- Test: `tests/runtime/tNoMyppualReferences.m` (new guard)

- [ ] **Step 1: Write the failing guard test**

`tests/runtime/tNoMyppualReferences.m`:

```matlab
classdef tNoMyppualReferences < matlab.unittest.TestCase
    % After retirement, no source under src/ may reference myppual /
    % convert_to_interp_eval_array / get_scalar / useFusedConstruct.
    methods (Test)
        function srcIsClean(tc)
            here = fileparts(mfilename('fullpath'));
            src  = fullfile(fileparts(fileparts(here)), 'src');
            files = dir(fullfile(src, '**', '*.m'));
            files = [files; dir(fullfile(src, '**', '*.cpp'))];
            bad = {};
            pats = {'myppual','convert_to_interp_eval_array','get_scalar','useFusedConstruct'};
            for i = 1:numel(files)
                t = fileread(fullfile(files(i).folder, files(i).name));
                for p = pats
                    if contains(t, p{1}); bad{end+1} = [files(i).name ':' p{1}]; end %#ok<AGROW>
                end
            end
            tc.verifyEmpty(bad, sprintf('stale references: %s', strjoin(bad, ', ')));
        end
    end
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `matlab -batch "cd('tests'); runtests('runtime/tNoMyppualReferences.m')"`
Expected: FAIL — many references remain (constructSplines legacy branch, kernels, comments).

- [ ] **Step 3: Prune `constructSplines.m`**

Replace the whole body of `src/+gdsge/+runtime/constructSplines.m` with the fused-only path (drop the `if gdsge.runtime.useFusedConstruct()` toggle and the `else` legacy `myppual`/`convert_to_interp_eval_array` branch):

```matlab
function [ppCell, splineVec] = constructSplines(valueCell, breaks, sizeState, interpOrder, extrapOrder, numThreads)
% CONSTRUCTSPLINES  Build interpolants + GDSGE_SPLINE_VEC via the fused
%   interp_construct_mex kernel (constructs all vars and writes the eval-order
%   array in one call). myppual is retired.
orderVec = interpOrder * ones(1, numel(sizeState));
if interpOrder == 4
    extrapVec = extrapOrder * ones(1, numel(sizeState));
else
    extrapVec = [];
end
nv = numel(valueCell);
if nv == 0
    ppCell = {}; splineVec = struct(); return;
end
values = struct();
for i = 1:nv
    values.(sprintf('GDSGE_V%d', i)) = reshape(valueCell{i}, [], sizeState{:});
end
[ppCell, splineVec] = interp_construct_mex(breaks, values, ...
    int32(orderVec), int32(extrapVec), int32(numThreads));
end
```

- [ ] **Step 4: Delete the files and clean comments**

```bash
git rm src/+gdsge/+runtime/useFusedConstruct.m \
       src/kernels/myppual.m src/kernels/myppual_mex.mexw64 \
       src/kernels/convert_to_interp_eval_array.m src/kernels/get_scalar.m
```

Edit `src/+gdsge/+runtime/Contents.m` and `src/+gdsge/+runtime/ensurePath.m` to remove the words `myppual` from their comment lines (replace with `interp_construct_mex` / `interp_eval_mex` as appropriate). Check `tests/runtime/tConstructSplinesEquivalence.m` — if it exercised the legacy branch via `useFusedConstruct(false)`, delete that case (the branch no longer exists).

- [ ] **Step 5: Run the guard + the construct equivalence + kernel tests**

Run: `matlab -batch "cd('tests'); runtests('runtime/tNoMyppualReferences.m')"`
Run: `matlab -batch "cd('tests'); runtests('runtime/tConstructSplines.m')"`
Run: `matlab -batch "cd('tests'); runtests('kernels/tKernels.m')"`
Expected: PASS (guard clean; constructs still work).

- [ ] **Step 6: Full suite green (myppual fully retired, still per-period eval)**

Run: `matlab -batch "cd('tests'); run_tests"` then check `tests/results/junit.xml` (exit 0 = all pass).
Expected: all pass. This is the Phase 2 milestone — myppual is gone, behavior preserved.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "refactor: delete myppual/get_scalar/convert_to_interp_eval_array; fused-only constructSplines"
```

---

# Phase 3 — Whole-simulation `simulate_<model>_mex` (the fast path)

Now move the SIMU_INTERP period loop into a codegen'd MEX. Per-period `interp_eval_mex` (Phase 2) stays as the documented fallback for compound transitions.

### Task 3.1: The C++ template `simulate_mex.tpl.cpp`

**Files:**
- Create: `templates/simulate_mex.tpl.cpp`

- [ ] **Step 1: Write the template**

`templates/simulate_mex.tpl.cpp` (placeholders filled by the emitter in Task 3.2: `GDSGE_NSTATES`, `GDSGE_ORDER`, `RECORD_AND_TRANSITION_CODE`, `STATE_FIELD_SETUP`, `VARSIMU_FIELD_SETUP`, `PRINT_FIELDS_CODE`):

```cpp
// Generated by gdsge.codegen.cxx.emitSimulateMex — do not edit.
// simulate_MODEL_NAME_mex(int32 numThreads, output_interp_pp, int32 nComp,
//                         shockMat, simuStruct, int32 [t0 t1], int32 num_periods,
//                         int32 EnforceSimuStateInbound, int32 SimuPrintFreq)
// Whole SIMU_INTERP period loop, in place. Stacked uniform-order interp:
// component c at realized shock s is vector index (s-1)*nComp + c (0-based).
#include "mex.h"
#include "MatlabPp.h"
#include <omp.h>
#define MAXDIM 10

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    int numThreads = *((int*)mxGetData(prhs[0]));
    omp_set_num_threads(numThreads);
    const mxArray* pp = prhs[1];
    int nComp = *((int*)mxGetData(prhs[2]));
    const double* shockMat = mxGetPr(prhs[3]);   // num_samples x (num_periods+1)
    mxArray* simu = (mxArray*)prhs[4];
    const int* trange = (const int*)mxGetData(prhs[5]);
    int t0 = trange[0], t1 = trange[1];
    int num_periods = *((int*)mxGetData(prhs[6]));
    int enforceInbound = *((int*)mxGetData(prhs[7]));
    int printFreq = *((int*)mxGetData(prhs[8]));

    int num_samples = (int)mxGetM(prhs[3]);
    const int XDIM = GDSGE_NSTATES;

    MatlabPp<GDSGE_NSTATES, GDSGE_ORDER> ev(pp);
    // state grid endpoints for inbound clamping
    double gmin[MAXDIM], gmax[MAXDIM];
    for (int d = 0; d < XDIM; ++d) { gmin[d] = ev.xGrid[d][0]; gmax[d] = ev.xGrid[d][ev.xPts[d]-1]; }

    // field pointers (in place)
    STATE_FIELD_SETUP        // double* state[XDIM]; each = mxGetPr(field) (column t = +t*num_samples)
    VARSIMU_FIELD_SETUP      // double* varsimu[...]; each = mxGetPr(field)

    for (int t = t0; t <= t1; ++t)
    {
        const double* shockT  = shockMat + (size_t)(t-1)*num_samples;
        const double* shockT1 = shockMat + (size_t)(t)*num_samples;   // next period
        #pragma omp parallel for
        for (int i = 0; i < num_samples; ++i)
        {
            double xSite[MAXDIM], xLeft[MAXDIM];
            int cell[MAXDIM];
            for (int d = 0; d < XDIM; ++d) {
                double s = state[d][(size_t)(t-1)*num_samples + i];
                if (enforceInbound) { if (s < gmin[d]) s = gmin[d]; if (s > gmax[d]) s = gmax[d]; }
                state[d][(size_t)(t-1)*num_samples + i] = s;     // matches MATLAB clamp-write
                xSite[d] = s;
            }
            ev.search(xSite, cell, xLeft);
            int base = ((int)shockT[i]-1) * nComp;
            int nextShock = (int)shockT1[i];                     // 1-based
            RECORD_AND_TRANSITION_CODE
        }
        if (printFreq > 0 && (t % printFreq == 0)) {
            mexPrintf("Periods: %d\n", t);
            PRINT_FIELDS_CODE
            mexEvalString("drawnow;");
        }
    }
}
```

- [ ] **Step 2: Commit the template (no test yet — exercised in 3.2/3.3)**

```bash
git add templates/simulate_mex.tpl.cpp
git commit -m "feat(templates): simulate_mex whole-loop C++ template"
```

### Task 3.2: The emitter `emitSimulateMex.m`

**Files:**
- Create: `src/+gdsge/+codegen/+cxx/emitSimulateMex.m`
- Test: `tests/codegen/tEmitSimulateMex.m`

Helper facts: a var's stacked row is `output_var_index.<name>` start (0-based). `nosearch_eval_<order>(xLeft, cell, base + row)`. For each `ir.simulate.varSimu{k}` (length-1 row `r_k`): `varsimu[k][(t-1)*num_samples+i] = ev.nosearch_eval_<order>(xLeft, cell, base + r_k);`. For each `ir.simulate.transitions{k}` with state slot `d`: unprimed → `state[d][t*num_samples+i] = ev.nosearch_eval(base + r);`; primed → `... base + r0 + (nextShock-1)` guarded by `if (t < num_periods)`.

- [ ] **Step 1: Write the failing test**

`tests/codegen/tEmitSimulateMex.m`:

```matlab
classdef tEmitSimulateMex < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addPath(tc)
            here = fileparts(mfilename('fullpath'));
            src  = fullfile(fileparts(fileparts(here)), 'src');
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(src));
        end
    end
    methods (Test)
        function emitsLoopAndTransitions(tc)
            ir = jsondecode(fileread(fullfile(tc.glswIrPath())));   % helper below
            txt = gdsge.codegen.cxx.emitSimulateMex(ir);
            tc.verifyClass(txt, 'char');
            tc.verifyTrue(contains(txt, 'omp parallel for'));
            tc.verifyTrue(contains(txt, 'nosearch_eval'));
            tc.verifyTrue(contains(txt, 'GDSGE_NSTATES'));   % template placeholder filled
            tc.verifyFalse(contains(txt, 'GDSGE_NSTATES;'),  'placeholder not substituted');
        end
    end
    methods
        function p = glswIrPath(tc) %#ok<MANU>
            here = fileparts(mfilename('fullpath'));
            p = fullfile(fileparts(here), 'GLSW2020','ir','GLSW_interp.gdsge.json');
        end
    end
end
```

> Replace `jsondecode(...)` with `gdsge.ir.decode(fileread(...))` if the IR JSON needs the typed decoder (check `tests/GLSW2020/ir/`). Use whatever the other codegen tests use to load an IR.

- [ ] **Step 2: Run to verify it fails**

Run: `matlab -batch "cd('tests'); runtests('codegen/tEmitSimulateMex.m')"`
Expected: FAIL — `emitSimulateMex` undefined.

- [ ] **Step 3: Write `emitSimulateMex.m`**

`src/+gdsge/+codegen/+cxx/emitSimulateMex.m`:

```matlab
function txt = emitSimulateMex(ir)
% EMITSIMULATEMEX  simulate_<model>_mex.cpp from templates/simulate_mex.tpl.cpp.
%   Whole SIMU_INTERP period loop in C++ (stacked uniform-order interp).
%   Only emitted for spline (non-ASG) SIMU_INTERP models whose simulate block
%   is MEX-expressible (var_simu = output vars; transitions = single output-var
%   reference). Caller (generateCxx) checks isMexExpressible first.
import gdsge.codegen.cxx.fillTemplate
import gdsge.codegen.cxx.readTemplate
nStates = numel(ir.states.names);
order   = ir.options.interpOrder;                 % 2 or 4

% 0-based stacked row of each output var (sum of preceding output lengths)
rowOf = containers.Map('KeyType','char','ValueType','double');
lenOf = containers.Map('KeyType','char','ValueType','double');
row = 0;
for i = 1:numel(ir.variables.output)
    n = ir.variables.output{i};
    L = outputLen(ir, n);
    rowOf(n) = row; lenOf(n) = L; row = row + L;
end

% state name -> 0-based slot index (column order of SimuRslt states)
stateSlot = containers.Map(ir.states.names, num2cell(0:nStates-1));

% RECORD_AND_TRANSITION_CODE
L = {};
for k = 1:numel(ir.simulate.varSimu)
    n = ir.simulate.varSimu{k};
    L{end+1} = sprintf('varsimu[%d][(size_t)(t-1)*num_samples+i] = ev.nosearch_eval_%d(xLeft, cell, base + %d);', ...
        k-1, order, rowOf(n)); %#ok<AGROW>
end
for k = 1:numel(ir.simulate.transitions)
    tr = ir.simulate.transitions{k};
    d  = stateSlot(tr.state);
    n  = tr.expr;                                  % single output-var name (validated)
    % Next-state column is 0-based index t (= MATLAB column t+1); written for
    % every t in [t0,t1] including t=num_periods (state fields are preallocated
    % with num_periods+1 columns to match MATLAB's auto-grow — see Task 3.4).
    if tr.primed
        L{end+1} = sprintf(['state[%d][(size_t)t*num_samples+i] = ' ...
            'ev.nosearch_eval_%d(xLeft, cell, base + %d + (nextShock-1));'], d, order, rowOf(n)); %#ok<AGROW>
    else
        L{end+1} = sprintf(['state[%d][(size_t)t*num_samples+i] = ' ...
            'ev.nosearch_eval_%d(xLeft, cell, base + %d);'], d, order, rowOf(n)); %#ok<AGROW>
    end
end
recordCode = strjoin(L, sprintf('\n            '));

% STATE_FIELD_SETUP / VARSIMU_FIELD_SETUP  (field order = states then varSimu; see emitResultSimu)
sf = {sprintf('double* state[%d];', nStates)};
for j = 1:nStates
    sf{end+1} = sprintf('state[%d] = mxGetPr(mxGetField(simu, 0, "%s"));', j-1, ir.states.names{j}); %#ok<AGROW>
end
stateSetup = strjoin(sf, sprintf('\n    '));
vf = {sprintf('double* varsimu[%d];', max(numel(ir.simulate.varSimu),1))};
for k = 1:numel(ir.simulate.varSimu)
    vf{end+1} = sprintf('varsimu[%d] = mxGetPr(mxGetField(simu, 0, "%s"));', k-1, ir.simulate.varSimu{k}); %#ok<AGROW>
end
varsimuSetup = strjoin(vf, sprintf('\n    '));

% PRINT_FIELDS_CODE — sample-1 values of each field (states then varSimu), parity-ish
pf = {};
allNames = [ir.states.names(:); ir.simulate.varSimu(:)];
for j = 1:numel(allNames)
    pf{end+1} = sprintf('mexPrintf("%%8.4g", mxGetPr(mxGetField(simu,0,"%s"))[(size_t)(t-1)*num_samples]);', allNames{j}); %#ok<AGROW>
end
printCode = strjoin([pf, {'mexPrintf("\\n");'}], sprintf('\n            '));

txt = readTemplate('simulate_mex.tpl.cpp');
txt = strrep(txt, 'MODEL_NAME', ir.modelName);
txt = fillTemplate(txt, { ...
    'GDSGE_NSTATES', num2str(nStates); ...
    'GDSGE_ORDER',   num2str(order); ...
    'STATE_FIELD_SETUP',   stateSetup; ...
    'VARSIMU_FIELD_SETUP', varsimuSetup; ...
    'RECORD_AND_TRANSITION_CODE', recordCode; ...
    'PRINT_FIELDS_CODE', printCode});
end

function L = outputLen(ir, name)
% output var length from the IR (1 for scalars; shock_num for primed/future vars).
L = 1;
for i = 1:numel(ir.variables.policy)
    if strcmp(ir.variables.policy{i}.name, name); L = ir.variables.policy{i}.length; return; end
end
for i = 1:numel(ir.variables.aux)
    if strcmp(ir.variables.aux{i}.name, name); L = ir.variables.aux{i}.length; return; end
end
end
```

> `fillTemplate` does word-boundary substitution; confirm `GDSGE_NSTATES`/`GDSGE_ORDER` appear in the template where boundaries hold (`<GDSGE_NSTATES, GDSGE_ORDER>` and `const int XDIM = GDSGE_NSTATES;`). If `fillTemplate` misses the `<...>` use, fall back to `strrep` for those two tokens. Verify `outputLen` against the real IR — if output lengths live elsewhere (e.g. `ir.variables.output` carries lengths), adapt; the Phase-2 `emitResultIter` already computes `lenOf` the same way (mirror it).

- [ ] **Step 4: Run to verify it passes**

Run: `matlab -batch "cd('tests'); runtests('codegen/tEmitSimulateMex.m')"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/+cxx/emitSimulateMex.m tests/codegen/tEmitSimulateMex.m
git commit -m "feat(codegen): emitSimulateMex — whole-loop simulate MEX emitter"
```

### Task 3.3: MEX-expressibility check + generateCxx wiring + compile template

**Files:**
- Create: `src/+gdsge/+codegen/+cxx/isSimuMexExpressible.m`
- Modify: `src/+gdsge/+codegen/generateCxx.m`
- Modify: `templates/compile.tpl.m`

- [ ] **Step 1: Write the expressibility predicate**

`src/+gdsge/+codegen/+cxx/isSimuMexExpressible.m`:

```matlab
function tf = isSimuMexExpressible(ir)
% ISSIMUMEXEXPRESSIBLE  true if a SIMU_INTERP spline model's simulate block can
%   run in simulate_<model>_mex: var_simu entries resolve to output vars, and
%   every transition RHS is a single output-var name (primed or unprimed).
%   Otherwise the model uses the per-period interp_eval_mex fallback.
tf = false;
if strcmp(ir.options.interpMethod, 'asg'); return; end
if ~(isfield(ir.options,'simuInterp') && ir.options.simuInterp==1); return; end
outNames = ir.variables.output;
isOut = @(n) any(strcmp(n, outNames));
for k = 1:numel(ir.simulate.varSimu)
    if ~isOut(ir.simulate.varSimu{k}); return; end
end
for k = 1:numel(ir.simulate.transitions)
    e = strtrim(ir.simulate.transitions{k}.expr);
    if isempty(regexp(e, '^[A-Za-z_]\w*$', 'once')) || ~isOut(e); return; end
end
tf = true;
end
```

- [ ] **Step 2: Wire `generateCxx.m`** — emit the simulate MEX + ensure the eval kernel

In `src/+gdsge/+codegen/generateCxx.m`, in the non-ASG (spline) branch where `ensureSplineConstructMex()` is called, add:

```matlab
gdsge.codegen.ensureSplineConstructMex();
gdsge.codegen.ensureInterpEvalMex();         % evaluator for resolve/warmup/fallback
```

After writing `files.cppFile`/`files.compileFile`, add:

```matlab
files.simuMexFile = '';
if gdsge.codegen.cxx.isSimuMexExpressible(ir)
    files.simuMexFile = fullfile(outDir, ['simulate_' ir.modelName '_mex.cpp']);
    gdsge.codegen.writeText(files.simuMexFile, gdsge.codegen.cxx.emitSimulateMex(ir));
end
```

And pass a flag to `emitCompile` so the compile script builds the simulate MEX too (Step 3). Change the `emitCompile` call to:

```matlab
gdsge.codegen.writeText(files.compileFile, ...
    gdsge.codegen.cxx.emitCompile(ir, includeDir, ~isempty(files.simuMexFile)));
```

- [ ] **Step 3: Add the second `mex` call to the compile template + emitCompile**

In `templates/compile.tpl.m`, before the final `end`, add a guarded second compile (placeholder `SIMU_MEX_COMPILE`):

```matlab
SIMU_MEX_COMPILE
end
```

In `src/+gdsge/+codegen/+cxx/emitCompile.m`, add a third argument `withSimuMex` (default false) and fill `SIMU_MEX_COMPILE`:

```matlab
function txt = emitCompile(ir, includeDir, withSimuMex)
...
if nargin < 3; withSimuMex = false; end
...
simuMexCompile = '';
if withSimuMex
    simuMexCompile = [ ...
        'simuCpp = fullfile(current_folder,''simulate_' ir.modelName '_mex.cpp'');' newline ...
        'simuFlags = '' -DMAXDIM=' num2str(maxDim) ' -DINTERP_ORDER=' num2str(ir.options.interpOrder) ''';' newline ...
        'eval([mexCommand simuFlags flag2 flag3 sprintf('' "%s"'',simuCpp) link_to_lib '' -outdir "'' current_folder ''"'' '' -I"'' include_folder ''"'']);'];
end
txt = fillTemplate(txt, { ...
    'INCLUDE_FOLDER',      strrep(includeDir, '\', '/'); ...
    'GDSGE_MAXDIM',        num2str(maxDim); ...
    'GDSGE_INTERP_ORDER',  num2str(ir.options.interpOrder); ...
    'EXTRA_DEF',           extraDef; ...
    'SIMU_MEX_COMPILE',    simuMexCompile});
```

> `simulate_<model>_mex.cpp` includes only `mex.h`/`MatlabPp.h`/`<omp.h>` — it does NOT need adept/Eigen/MKL defines, but compiling with the same `flag2/flag3` (OpenMP, `link_to_lib`) is harmless and reuses the established flags. `MatlabPp.h` pulls `interp_lite.h` (header-only); no extra link needed beyond `essential_blas.lib` already in `link_to_lib`.

- [ ] **Step 4: Update `codegen.m` cache to also gate the simulate MEX**

In `src/+gdsge/+codegen/codegen.m`, the cache currently hashes only `mex_<model>.cpp`. Extend the cache text to include the simulate MEX source so a change to it triggers recompile:

```matlab
cppText = fileread(files.cppFile);
if isfield(files,'simuMexFile') && ~isempty(files.simuMexFile)
    cppText = [cppText, fileread(files.simuMexFile)];
end
```

(Leave the rest of the compile-gate logic unchanged; `compile_<model>` builds both MEX files.)

- [ ] **Step 5: Smoke-compile via a model codegen**

Run: `matlab -batch "cd('tests'); runtests('GLSW2020/codegen/tEndToEndGLSW.m')"`
Expected: PASS — `gdsge_codegen('GLSW_interp')` now also writes+compiles `simulate_GLSW_interp_mex.cpp`. (The simulate `.m` doesn't call it yet — Task 3.4 — so this only proves emission+compile.)

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+codegen/+cxx/isSimuMexExpressible.m src/+gdsge/+codegen/generateCxx.m src/+gdsge/+codegen/codegen.m src/+gdsge/+codegen/+cxx/emitCompile.m templates/compile.tpl.m
git commit -m "feat(codegen): emit+compile simulate_<model>_mex for expressible SIMU_INTERP models"
```

### Task 3.4: Switch `emitSimulateInterp` to call the whole-sim MEX (with fallback)

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitSimulateInterp.m`
- Test: `tests/GLSW2020/codegen/tSimuInterpMexGLSW.m`, `tests/Mendoza2010/codegen/tSimuInterpMexMendoza.m`

- [ ] **Step 1: Write the failing parity tests**

`tests/GLSW2020/codegen/tSimuInterpMexGLSW.m` (mirror the e2e harness; assert SimuRslt vs golden AND that the generated simulate `.m` calls the MEX):

```matlab
classdef tSimuInterpMexGLSW < matlab.unittest.TestCase
    methods (Test, TestTags={'Slow'})
        function simuMexMatchesGolden(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir,'GLSW_interp.gmod'), work);
            tc.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(work));
            gdsge_codegen('GLSW_interp');
            simTxt = fileread(fullfile(work,'simulate_GLSW_interp.m'));
            tc.verifyTrue(contains(simTxt, 'simulate_GLSW_interp_mex'), 'simulate .m must call the MEX');
            rng(0); IterRslt = iter_GLSW_interp;
            rng(0); SimuRslt = simulate_GLSW_interp(IterRslt);
            g = load(fullfile(modelDir,'golden','SimuRslt.mat'));   % adapt to actual golden path
            tc.verifyEqual(SimuRslt.a1, g.SimuRslt.a1, 'RelTol', 1e-6, 'AbsTol', 1e-6);
            for f = {'a2','P1','r','c1_shr'}
                tc.verifyEqual(SimuRslt.(f{1}), g.SimuRslt.(f{1}), 'RelTol', 1e-6, 'AbsTol', 1e-6);
            end
        end
    end
end
```

> Check the real golden location/shape for GLSW SimuRslt (look at `tests/GLSW2020/codegen/tEndToEndGLSW.m` for how it compares SimuRslt today) and mirror it. Same structure for `tSimuInterpMexMendoza.m`, comparing `cTilde`,`k` (states, primed) + a couple `var_simu` fields.

- [ ] **Step 2: Run to verify they fail**

Run: `matlab -batch "cd('tests'); runtests('GLSW2020/codegen/tSimuInterpMexGLSW.m')"`
Expected: FAIL — generated simulate `.m` does not yet call `simulate_GLSW_interp_mex`.

- [ ] **Step 3: Rewrite `emitSimulateInterp.m` spline branch to call the MEX (with fallback)**

In `src/+gdsge/+codegen/+mat/emitSimulateInterp.m`, for the spline branch, branch on `gdsge.codegen.cxx.isSimuMexExpressible(ir)`:
- **Expressible** → emit the chunked whole-sim call (no per-period MATLAB loop). The shock path is pre-generated; loop chunks of `SimuSaveFreq`:

```matlab
if gdsge.codegen.cxx.isSimuMexExpressible(ir)
    w.add('GDSGE_NCOMP = IterRslt.output_interp.GDSGE_NCOMP;');
    % State fields need num_periods+1 columns (MATLAB auto-grows them via the
    % (:,t+1) transition write at t=num_periods; the MEX cannot grow, so size
    % up front to match the golden SimuRslt state shape). var_simu stay num_periods.
    for i = 1:numel(ir.states.names)
        w.add('SimuRslt.%s(:,num_periods+1) = 0;', ir.states.names{i});
    end
    w.add('GDSGE_T0 = 1;');
    w.add('while GDSGE_T0 <= num_periods');
    w.add('    GDSGE_T1 = min(GDSGE_T0 + SimuSaveFreq - 1, num_periods);');
    w.add('    simulate_%s_mex(int32(NumThreads), IterRslt.output_interp, int32(GDSGE_NCOMP), ...', m);
    w.add('        SimuRslt.shock, SimuRslt, int32([GDSGE_T0 GDSGE_T1]), int32(num_periods), ...');
    w.add('        int32(EnforceSimuStateInbound), int32(SimuPrintFreq));');
    w.add('    if mod(GDSGE_T1,SimuSaveFreq)==0');
    w.add('        save([''SimuRslt_%s_'' num2str(GDSGE_T1) ''.mat''], ''SimuRslt'');', m);
    w.add('    end');
    w.add('    GDSGE_T0 = GDSGE_T1 + 1;');
    w.add('end');
else
    % per-period interp_eval_mex loop (Phase 2 code path) — keep verbatim
    ... existing Task 2.2 loop ...
end
```

> `SimuRslt` is passed to the MEX, which writes its state/var_simu columns in place (hans pattern: `mxGetPr` of each field). Because `SimuRslt` is both the 5th input and the returned variable in the caller, the in-place write is safe under MATLAB copy-on-write. Keep the ASG branch and the Phase-2 per-period branch unchanged.

- [ ] **Step 4: Run the parity tests**

Run: `matlab -batch "cd('tests'); runtests('GLSW2020/codegen/tSimuInterpMexGLSW.m')"`
Run: `matlab -batch "cd('tests'); runtests('Mendoza2010/codegen/tSimuInterpMexMendoza.m')"`
Expected: PASS — SimuRslt matches goldens; simulate `.m` calls the MEX.

- [ ] **Step 5: Full regression suite**

Run: `matlab -batch "cd('tests'); run_tests"` → check `tests/results/junit.xml` (exit 0).
Expected: all pass — including CaoKS2016_simu_interp (ASG, unchanged), HL1996/safe_assets (resolve+warmup), and the new MEX paths.

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+codegen/+mat/emitSimulateInterp.m tests/GLSW2020/codegen/tSimuInterpMexGLSW.m tests/Mendoza2010/codegen/tSimuInterpMexMendoza.m
git commit -m "feat(codegen): SIMU_INTERP spline simulation runs whole loop in simulate_<model>_mex"
```

### Task 3.5: Perf check + docs

**Files:**
- Modify: `PROGRESS.md`, `docs/perf-report.md`

- [ ] **Step 1: Measure GLSW + Mendoza simulate wall-clock before/after**

In a scratch script (not committed), time `simulate_<model>` with the MEX path vs. the Phase-2 per-period path (temporarily force the fallback). Record the speedup.

- [ ] **Step 2: Record results**

Append a `docs/perf-report.md` entry (model, num_periods×num_samples, per-period-eval vs whole-sim-MEX seconds) and tick the phase in `PROGRESS.md`.

- [ ] **Step 3: Commit**

```bash
git add PROGRESS.md docs/perf-report.md
git commit -m "docs(perf): record SIMU_INTERP whole-sim MEX speedup"
```

---

## Self-review notes (filled during writing)

- **Spec coverage:** Component 1 (simulate MEX) → Tasks 3.1–3.4; Component 2 (`interp_eval_mex`) → Task 1.1; Component 3 (simulate `.m`) → Tasks 2.2, 3.4; Component 4 (construct/eval swaps) → Tasks 2.1, 2.3; Component 5 (build plumbing) → Task 3.3; Component 6 (fallback) → Tasks 3.3 (`isSimuMexExpressible`), 3.4 (else branch); deletions → Task 2.4; stacked-layout enabler → Tasks 2.1, 2.3; testing → tests in every task + guard 2.4/Step 1.
- **Ordering rationale:** Phase 2 retires myppual with a low-risk per-period eval and proves the layout change against goldens *before* Phase 3 adds the MEX loop. Each phase ends with a green suite.
- **Open verifications for the implementer (do these as you reach the task, not upfront):**
  - Confirm `OutputInterpOrder` is in scope in `emitResultIter`/`emitSimulate` setup (it is used today at lines 29/57); if named differently, use the actual option local.
  - Confirm the exact GLSW/Mendoza SimuRslt golden paths/fields from the existing e2e tests and mirror them in 3.4/Step 1.
  - Confirm `fillTemplate` substitutes the `<GDSGE_NSTATES, GDSGE_ORDER>` tokens; else `strrep` them.
  - Confirm `ir.simulate.transitions{k}.primed`/`.expr`/`.state` field names (from `parseSimulate.m`) — they are `state`, `expr`, `primed`.
  - Confirm the golden SimuRslt **state** fields have `num_periods+1` columns (MATLAB auto-grows them) while `var_simu` fields have `num_periods`; the Phase-3 prealloc and the MEX's `t=num_periods` next-state write depend on this. If a golden differs, adjust the prealloc and the chunk's final `t1`.
```
