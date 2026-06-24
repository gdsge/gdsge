# In-MEX Nearby-Warmup Resolve Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the nearby-neighbor warm-start sweep from MATLAB (`solveProblems.m` Phase 2) into the C++ MEX (`task.tpl.cpp`), so the whole sweep is one MEX call per outer iteration instead of `~2·num_dims·K` round-trips — with bit-exact parity to today's goldens.

**Architecture:** `task.tpl.cpp` reads an optional per-dimension `GDSGE_PROBLEM_STRIDES` vector from the caller workspace; when present, after the initial parallel solve it runs the neighbor sweep internally with deterministic snapshot semantics. `solveProblems.m` keeps its exact control flow but replaces the *body* of its MATLAB sweep with a single strides-enabled MEX call, guarded by the existing persistent `numNeedResolvedAfter` entry-gate. The switch is an undocumented `UseMexResolve` runtime option (default 1). ASG and the init task pass no strides → unchanged.

**Tech Stack:** MATLAB R2025b (`matlab.unittest`, headless via `matlab -batch`), C++ MEX (MSVC 2022, OpenMP, adept/CoDoSol), template-substitution codegen.

**Reference spec:** `docs/superpowers/specs/2026-06-14-in-mex-nearby-resolve-design.md`.

**Ground rules for every task:**
- Run MATLAB **one process at a time** (each saturates all cores via OpenMP; concurrent runs distort everything). Never launch two `matlab -batch` at once.
- The authoritative pass/fail signal is the `runtests` result object printed to stdout (or `results/junit.xml`); `results.tap` accumulates stale lines across runs — ignore it.
- Targeted test run (from repo root):
  ```bash
  matlab -batch "cd('tests'); addpath(pwd,'../src','../src/kernels'); r=runtests('REL/PATH/TO/tClass.m'); disp(r); exit(double(any([r.Failed])))"
  ```
- Full suite: `matlab -batch "cd('tests'); run_tests"` (exit 0 = all pass).
- The SymPy backend tests need the `uv` Python env (`uv sync` under `pyext/`); they `assumeTrue(sympyAvailable)` and skip if absent.

---

## File Structure

| File | Responsibility | Task |
|------|----------------|------|
| `templates/cxx/task.tpl.cpp` | per-point solve as a reusable lambda; Pass-0 loop; optional in-MEX snapshot resolve driver | 1 |
| `templates/cxx/call_fmin.tpl.cpp` | drop the `sol0→sol` memcpy (warm start pre-placed by caller) | 1 |
| `templates/cxx/call_fmin_sympy.tpl.cpp` | same | 1 |
| `tests/HeatonLucas1996/codegen/golden/mex_HL1996_golden_cpp.txt` | regenerated C++ snapshot | 1 |
| `src/+gdsge/+runtime/solveProblems.m` | `GDSGE_PROBLEM_STRIDES` contract local; `cfg.useMexResolve`; dedicated sweep call w/ persistent entry-gate | 2 |
| `tests/runtime/InMexResolveRecorder.m` | test helper: handle class that records MEX-call structure | 2 |
| `tests/runtime/tInMexResolveDriver.m` | fast unit test of the MATLAB driver's call structure | 2 |
| `src/+gdsge/+codegen/+mat/emitSetup.m` | default `UseMexResolve = 1;` | 3 |
| `src/+gdsge/+codegen/+mat/optionsWhitelist.m` | add `'UseMexResolve'` to the whitelist | 3 |
| `src/+gdsge/+codegen/+mat/emitIter.m` | `GDSGE_CFG.useMexResolve = UseMexResolve;` | 3 |
| `tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt` | regenerated iter `.m` snapshot | 3 |
| `tests/Barro_et_al_2017/codegen/tInMexResolveSafeAssets.m` | end-to-end A/B bit-identical gate | 4 |
| `tests/perf/inmex_resolve_bench.m` | informational round-trip/timing benchmark script | 5 |
| `PROGRESS.md` | changelog entry | 5 |

---

## Task 1: C++ in-MEX resolve driver (no-strides behavior unchanged)

This task changes the shared C++ task template and both `call_fmin` templates. Because nothing supplies `GDSGE_PROBLEM_STRIDES` yet (that lands in Tasks 2–3), the resolve driver is dormant (`num_dims == 0`), so the only observable effect must be **zero numeric change** — the lambda refactor + memcpy relocation is behavior-preserving. The C++ snapshot text changes and is regenerated.

**Files:**
- Modify: `templates/cxx/task.tpl.cpp` (full rewrite below)
- Modify: `templates/cxx/call_fmin.tpl.cpp`
- Modify: `templates/cxx/call_fmin_sympy.tpl.cpp`
- Regenerate: `tests/HeatonLucas1996/codegen/golden/mex_HL1996_golden_cpp.txt`

- [ ] **Step 1: Rewrite `templates/cxx/task.tpl.cpp`** to factor the per-point body into a reusable lambda, add the optional strides read, and append the snapshot resolve driver. Replace the entire file with:

```cpp
// This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE
void TASK_NAME(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  GET_DBL(TolFun);
  GET_DBL(TolSol);
  GET_INT(SolMaxIter);
  GET_INT(NumThreads);
  GET_INT(GDSGE_DEBUG_EVAL_ONLY);
  GET_INT(UseBroyden);
  GET_DBL(FiniteDiffDelta);
  GET_INT(GDSGE_USE_BROYDEN_NOW);

  #ifndef NO_OMP
  omp_set_num_threads(NumThreads);
  #endif

  //
  int GDSGE_NPROB = mxGetN(prhs[3]);
  int GDSGE_NUM_SOL0 = mxGetN(prhs[0]);

  // Input
  GET_DM_VIEW_FROM_MX(GDSGE_SOL0,prhs[0]);
  GET_DM_VIEW_FROM_MX(GDSGE_LB,prhs[1]);
  GET_DM_VIEW_FROM_MX(GDSGE_UB,prhs[2]);
  GET_DM_VIEW_FROM_MX(GDSGE_DATA,prhs[3]);
  GET_DV_VIEW_FROM_MX(GDSGE_SKIP,prhs[4]);
  GET_DV_VIEW_FROM_MX(GDSGE_F0,prhs[5]);
  GET_DM_VIEW_FROM_MX(GDSGE_AUX0,prhs[6]);
  GET_DM_VIEW_FROM_MX(GDSGE_EQVAL0,prhs[7]);

  // Output
  plhs[0] = mxCreateDoubleMatrix(NUM_EQUATIONS,GDSGE_NPROB,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,GDSGE_NPROB,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(NUM_AUX,GDSGE_NPROB,mxREAL);
  plhs[3] = mxCreateDoubleMatrix(NUM_EQUATIONS,GDSGE_NPROB,mxREAL);
  plhs[4] = mxCreateDoubleMatrix(1,GDSGE_NPROB,mxREAL);

  GET_DM_VIEW_FROM_MX(GDSGE_SOL,plhs[0]);
  GET_DV_VIEW_FROM_MX(GDSGE_F,plhs[1]);
  GET_DM_VIEW_FROM_MX(GDSGE_AUX,plhs[2]);
  GET_DM_VIEW_FROM_MX(GDSGE_EQVAL,plhs[3]);
  GET_DV_VIEW_FROM_MX(GDSGE_OPT_INFO,plhs[4]);

  // Optional analytic-Jacobian debug output (Phase 8 cross-check): per problem,
  // NUM_EQUATIONS*NUM_EQUATIONS column-major, filled only when nlhs>=6 and
  // GDSGE_DEBUG_EVAL_ONLY==2. Existing 5-output callers are unaffected.
  double* GDSGE_JAC_OUT = 0;
  if (nlhs >= 6) {
    plhs[5] = mxCreateDoubleMatrix(NUM_EQUATIONS*NUM_EQUATIONS, GDSGE_NPROB, mxREAL);
    GDSGE_JAC_OUT = mxGetPr(plhs[5]);
  }

  // Input and output
  INTERP_GET_CODE;

  // Optional in-MEX nearby-warmup resolve. Read the per-dimension problem
  // strides from the caller workspace; absent/empty => no resolve (Pass 0 only),
  // which keeps the init task, ASG, and the MATLAB-sweep fallback on today's path.
  int GDSGE_NUM_DIMS = 0;
  double* GDSGE_STRIDES = 0;
  {
    const mxArray* GDSGE_STRIDES_MX = mexGetVariablePtr("caller", "GDSGE_PROBLEM_STRIDES");
    if (GDSGE_STRIDES_MX != 0 && mxIsDouble(GDSGE_STRIDES_MX) && mxGetNumberOfElements(GDSGE_STRIDES_MX) > 0) {
      GDSGE_NUM_DIMS = (int) mxGetNumberOfElements(GDSGE_STRIDES_MX);
      GDSGE_STRIDES = mxGetPr(GDSGE_STRIDES_MX);
    }
  }

  // Per-grid-point solve, reusable across the initial pass and each resolve pass.
  // CONTRACT: the warm-start guess must already be in GDSGE_SOL(:,GDSGE_I) before
  // this runs. Pass 0 copies it from GDSGE_SOL0; a resolve pass copies it from a
  // converged neighbor. The solve/eval operate in place on GDSGE_SOL.
  auto GDSGE_solve_one = [&](int GDSGE_I)
  {
    START_LOOP_CODE;

    //Initiate adept
    Stack _stack;
    int GDSGE_EVAL_GRAD_FLAG = 0;

    INTERP_GET_THREAD_CODE;

    double* GDSGE_data = &GDSGE_DATA(1, GDSGE_I);
    double* GDSGE_aux = &GDSGE_AUX(1, GDSGE_I);
    double* GDSGE_eqval = &GDSGE_EQVAL(1, GDSGE_I);

    int GDSGE_POP_I = 0;
#define POPN(var) double var = GDSGE_data[GDSGE_POP_I++]
#define POPNARRAY(var, length) double* var = &GDSGE_data[GDSGE_POP_I]; GDSGE_POP_I+=length
#define POPNARRAY_BASE1(var, length) double* var = &GDSGE_data[GDSGE_POP_I-1]; GDSGE_POP_I+=length

    POP_CODE;

    PRE_MODEL_CODE;

    MODEL_CODE;

    CALL_FMIN_CODE;

    FINISH_LOOP_CODE;
#undef POPN
#undef POPNARRAY
#undef POPNARRAY_BASE1
  };

  // Pass 0: solve every non-skipped point from its supplied guess.
#pragma omp parallel for schedule(dynamic)
  for (int GDSGE_I = 1; GDSGE_I <= GDSGE_NPROB; GDSGE_I++)
  {
    GDSGE_F(GDSGE_I) = GDSGE_F0(GDSGE_I);
    if (GDSGE_SKIP(GDSGE_I) == 1)
    {
      memcpy(&GDSGE_SOL(1,GDSGE_I), &GDSGE_SOL0(1,GDSGE_I), sizeof(double)*NUM_EQUATIONS);
      memcpy(&GDSGE_EQVAL(1,GDSGE_I), &GDSGE_EQVAL0(1,GDSGE_I), sizeof(double)*NUM_EQUATIONS);
      memcpy(&GDSGE_AUX(1,GDSGE_I), &GDSGE_AUX0(1,GDSGE_I), sizeof(double)*NUM_AUX);
      continue;
    }
    memcpy(&GDSGE_SOL(1,GDSGE_I), &GDSGE_SOL0(1,GDSGE_I), sizeof(double)*NUM_EQUATIONS);
    GDSGE_solve_one(GDSGE_I);
  }

  // In-MEX nearby-warmup resolve (deterministic snapshot semantics; bit-exact
  // with solveProblems.m's MATLAB sweep). Runs only when strides were supplied.
  if (GDSGE_NUM_DIMS > 0)
  {
    vector<char> GDSGE_SOLVED0(GDSGE_NPROB + 1);
    auto GDSGE_count_unconv = [&]() -> int {
      int c = 0;
      for (int i = 1; i <= GDSGE_NPROB; i++) if (!(GDSGE_F(i) <= TolSol)) c++;
      return c;
    };
    int GDSGE_BEFORE = GDSGE_count_unconv();
    int GDSGE_AFTER = 0x7fffffff;
    while (GDSGE_AFTER != GDSGE_BEFORE)
    {
      GDSGE_BEFORE = GDSGE_count_unconv();
      for (int GDSGE_DIM = 0; GDSGE_DIM < GDSGE_NUM_DIMS; GDSGE_DIM++)
      {
        int GDSGE_STRIDE = (int) GDSGE_STRIDES[GDSGE_DIM];

        // Lower neighbor: snapshot the converged set, then fused copy + solve.
        for (int i = 1; i <= GDSGE_NPROB; i++)
          GDSGE_SOLVED0[i] = (GDSGE_F(i) <= TolSol) ? 1 : 0;
#pragma omp parallel for schedule(dynamic)
        for (int GDSGE_I = 1; GDSGE_I <= GDSGE_NPROB; GDSGE_I++)
        {
          int GDSGE_J = GDSGE_I - GDSGE_STRIDE;
          if (!GDSGE_SOLVED0[GDSGE_I] && GDSGE_J >= 1 && GDSGE_SOLVED0[GDSGE_J])
          {
            memcpy(&GDSGE_SOL(1,GDSGE_I), &GDSGE_SOL(1,GDSGE_J), sizeof(double)*NUM_EQUATIONS);
            GDSGE_solve_one(GDSGE_I);
          }
        }

        // Upper neighbor.
        for (int i = 1; i <= GDSGE_NPROB; i++)
          GDSGE_SOLVED0[i] = (GDSGE_F(i) <= TolSol) ? 1 : 0;
#pragma omp parallel for schedule(dynamic)
        for (int GDSGE_I = 1; GDSGE_I <= GDSGE_NPROB; GDSGE_I++)
        {
          int GDSGE_J = GDSGE_I + GDSGE_STRIDE;
          if (!GDSGE_SOLVED0[GDSGE_I] && GDSGE_J <= GDSGE_NPROB && GDSGE_SOLVED0[GDSGE_J])
          {
            memcpy(&GDSGE_SOL(1,GDSGE_I), &GDSGE_SOL(1,GDSGE_J), sizeof(double)*NUM_EQUATIONS);
            GDSGE_solve_one(GDSGE_I);
          }
        }
      }
      GDSGE_AFTER = GDSGE_count_unconv();
    }
  }
}
```

Notes for the implementer: `START_LOOP_CODE … FINISH_LOOP_CODE` are template placeholders substituted by the emitters (per-backend / per-region) — keep them verbatim and in this order; only their *enclosing* structure changed (now a lambda + an explicit Pass-0 loop). The `vector<char>` uses the `<vector>` already included by `mex.tpl.cpp`. `memcpy` is already used in the old template. The C++11 lambda + `auto` are within MSVC 2022's default `/std:c++14` (Eigen/adept already require C++11+), so `compile.tpl.m` needs no standard-flag change — the recompile in Step 7 confirms it builds.

- [ ] **Step 2: Edit `templates/cxx/call_fmin.tpl.cpp`** — remove both `sol0→sol` memcpys (the warm start is now pre-placed in `GDSGE_SOL` by the caller). Replace the whole file with:

```cpp
// This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE
if (MODEL_CONDITION) {
  if (GDSGE_DEBUG_EVAL_ONLY==0)
  {
    // Warm start already in GDSGE_SOL (placed by the caller); solve in place.
    GDSGE_F(GDSGE_I) = CoDoSol::solve(NUM_EQUATIONS, &GDSGE_SOL(1, GDSGE_I), &GDSGE_LB(1, GDSGE_I), &GDSGE_UB(1, GDSGE_I), UseBroyden, TolFun, 0, SolMaxIter, GDSGE_OBJ_MODEL_NUMBER, &GDSGE_OPT_INFO(GDSGE_I));
  }

  // Evaluate at the solution
  double* GDSGE_x = &GDSGE_SOL(1, GDSGE_I);
  #if MAXDIM>MAX_STACK_DIM
  vector<double> GDSGE_EQ(NUM_EQUATIONS);
  #else
  double GDSGE_EQ[NUM_EQUATIONS];
  #endif

  if (GDSGE_DEBUG_EVAL_ONLY==2 && GDSGE_JAC_OUT)
  {
    GDSGE_OBJ_MODEL_NUMBER(&GDSGE_x[0], &GDSGE_EQ[0], &GDSGE_JAC_OUT[(GDSGE_I-1)*NUM_EQUATIONS*NUM_EQUATIONS]);
  }
  else
  {
    GDSGE_FUNC_MODEL_NUMBER_double(&GDSGE_x[0], &GDSGE_EQ[0], 1);
  }
}
```

- [ ] **Step 3: Edit `templates/cxx/call_fmin_sympy.tpl.cpp`** — same memcpy removal (note the eval call has no trailing `, 1`). Replace the whole file with:

```cpp
// This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE
// SymPy backend: same CoDoSol drive as call_fmin.tpl.cpp, but the model
// function is all-double (x,f,jac) — no adept eval flag. Value-only eval omits
// the jac pointer; the 6th-output cross-check fills it via the OBJ functor.
if (MODEL_CONDITION) {
  if (GDSGE_DEBUG_EVAL_ONLY==0)
  {
    // Warm start already in GDSGE_SOL (placed by the caller); solve in place.
    GDSGE_F(GDSGE_I) = CoDoSol::solve(NUM_EQUATIONS, &GDSGE_SOL(1, GDSGE_I), &GDSGE_LB(1, GDSGE_I), &GDSGE_UB(1, GDSGE_I), UseBroyden, TolFun, 0, SolMaxIter, GDSGE_OBJ_MODEL_NUMBER, &GDSGE_OPT_INFO(GDSGE_I));
  }

  // Evaluate at the solution
  double* GDSGE_x = &GDSGE_SOL(1, GDSGE_I);
  #if MAXDIM>MAX_STACK_DIM
  vector<double> GDSGE_EQ(NUM_EQUATIONS);
  #else
  double GDSGE_EQ[NUM_EQUATIONS];
  #endif

  if (GDSGE_DEBUG_EVAL_ONLY==2 && GDSGE_JAC_OUT)
  {
    GDSGE_OBJ_MODEL_NUMBER(&GDSGE_x[0], &GDSGE_EQ[0], &GDSGE_JAC_OUT[(GDSGE_I-1)*NUM_EQUATIONS*NUM_EQUATIONS]);
  }
  else
  {
    GDSGE_FUNC_MODEL_NUMBER_double(&GDSGE_x[0], &GDSGE_EQ[0]);
  }
}
```

- [ ] **Step 4: Run the C++ snapshot test to confirm it now fails (drift)**

Run:
```bash
matlab -batch "cd('tests'); addpath(pwd,'../src','../src/kernels'); r=runtests('HeatonLucas1996/codegen/tSnapshotCxxHL1996.m'); disp(r); exit(double(any([r.Failed])))"
```
Expected: `cppMatchesSnapshot` FAILS with "cpp snapshot drifted" (the generated C++ now uses the lambda + resolve driver). `compileMatchesSnapshot` PASSES (compile flags unchanged).

- [ ] **Step 5: Regenerate the snapshot and review the diff**

Run:
```bash
matlab -batch "cd('tests/HeatonLucas1996/codegen'); addpath(fullfile(pwd,'..','..'),fullfile(pwd,'..','..','..','src'),fullfile(pwd,'..','..','..','src','kernels')); regen_snapshots"
```
Then inspect: `git diff tests/HeatonLucas1996/codegen/golden/mex_HL1996_golden_cpp.txt`. Confirm the diff shows ONLY: the lambda wrapper, the new Pass-0 loop, the strides read, the resolve driver, and the memcpy moved out of `call_fmin`. The model body, POP code, and interp code must be unchanged. Confirm `iter_HL1996_golden.txt`/`compile_HL1996_golden.txt`/`simulate_HL1996_golden.txt` are unchanged (`git status` shows only the cpp golden modified).

- [ ] **Step 6: Run the snapshot test to confirm it passes**

Run the Step-4 command again.
Expected: both `cppMatchesSnapshot` and `compileMatchesSnapshot` PASS.

- [ ] **Step 7: Verify no-strides equivalence end-to-end (both backends)**

This recompiles the HL1996 MEX from the new template and runs full iter + simulate vs the golden. Because no code supplies `GDSGE_PROBLEM_STRIDES` yet, the resolve driver stays dormant; the result must match the golden exactly.

Run:
```bash
matlab -batch "cd('tests'); addpath(pwd,'../src','../src/kernels'); r=runtests('HeatonLucas1996/codegen/tEndToEndHL1996.m'); disp(r); exit(double(any([r.Failed])))"
```
Expected: PASS (Iter=209, golden match). If the SymPy `uv` env is available, also run `tEndToEndHL1996Sympy.m` the same way and expect PASS; if Python is absent it will be filtered/skipped.

- [ ] **Step 8: Commit**

```bash
git add templates/cxx/task.tpl.cpp templates/cxx/call_fmin.tpl.cpp templates/cxx/call_fmin_sympy.tpl.cpp tests/HeatonLucas1996/codegen/golden/mex_HL1996_golden_cpp.txt
git commit -m "feat(cxx): in-MEX nearby-warmup resolve driver in task.tpl.cpp (dormant until strides supplied)

Per-point solve refactored into a reusable lambda; Pass-0 loop copies the
supplied guess into the working GDSGE_SOL then solves; an optional snapshot
resolve driver (gated on a non-empty GDSGE_PROBLEM_STRIDES from the caller
workspace) propagates converged neighbors and re-solves until no progress.
Both call_fmin templates drop the sol0->sol memcpy (warm start now pre-placed
by the caller). No code supplies strides yet, so behavior is unchanged:
HL1996 end-to-end matches the golden bit-for-bit on both backends.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: MATLAB driver — dedicated in-MEX sweep call

Add `cfg.useMexResolve` handling to `solveProblems.m`: keep the control flow byte-for-byte, but when on, replace the body of the MATLAB neighbor sweep with one strides-enabled MEX call guarded by the persistent `numNeedResolvedAfter` entry-gate. Fast TDD via a fake MEX that records the call structure (no compile).

**Files:**
- Create: `tests/runtime/InMexResolveRecorder.m`
- Create: `tests/runtime/tInMexResolveDriver.m`
- Modify: `src/+gdsge/+runtime/solveProblems.m`

- [ ] **Step 1: Write the test helper `tests/runtime/InMexResolveRecorder.m`**

A handle class whose `mex` method has the MEX signature and is passed as `@rec.mex` so `solveProblems` calls it directly (its `evalin('caller',…)` then reads `solveProblems`'s workspace, exactly as the real MEX reads the caller-workspace contract).

```matlab
classdef InMexResolveRecorder < handle
    % Test double for the model MEX. Records each call's strides + skip and
    % simulates convergence so solveProblems' resolve loops terminate.
    %   call 1            : Pass 0 — converge only problem 1 (rest stay unconverged)
    %   strides non-empty : in-MEX sweep — converge everything
    %   else              : MATLAB-sweep direction / random restart — converge
    %                       the non-skipped (skip==0) problems
    properties
        calls = {};   % each: struct('stridesEmpty',logical,'strides',[],'skip',[])
        n = 0;
    end
    methods
        function [sol, f, aux, eqVal, optInfo] = mex(obj, sol, lb, ub, data, skip, f, aux, eqVal) %#ok<INUSL>
            obj.n = obj.n + 1;
            strides = evalin('caller', 'GDSGE_PROBLEM_STRIDES');
            obj.calls{end+1} = struct('stridesEmpty', isempty(strides), ...
                'strides', strides(:)', 'skip', skip(:)');
            if obj.n == 1
                f(:) = 1e20; f(1) = 0;                 % Pass 0
            elseif ~isempty(strides)
                f(:) = 0;                              % in-MEX sweep converges all
            else
                f(skip == 0) = 0;                      % MATLAB sweep / restart
            end
            optInfo = zeros(1, numel(f));
        end
    end
end
```

- [ ] **Step 2: Write the failing test `tests/runtime/tInMexResolveDriver.m`**

```matlab
classdef tInMexResolveDriver < matlab.unittest.TestCase
    % Fast unit test of solveProblems' MATLAB driver under cfg.useMexResolve.
    methods (Static)
        function cfg = baseCfg(useMexResolve)
            cfg = struct('tolFun',1e-8,'tolSol',1e-8,'solMaxIter',200, ...
                'numThreads',1,'debugEvalOnly',0,'useBroyden',0, ...
                'finiteDiffDelta',1e-6,'useBroydenNow',0,'taskName',1, ...
                'splineVec',[],'ppNames',{{}},'ppCell',{{}},'maxMinorIter',50, ...
                'probSize',[2 3],'useNearestNeighbor',true,'verboseRetry',false, ...
                'adaptInSol',[],'useMexResolve',useMexResolve);
        end
        function [sol,lb,ub,data,f,aux,eqVal] = batch()
            n = 6;                       % prod([2 3])
            sol = zeros(1,n); lb = zeros(1,n); ub = ones(1,n);
            data = zeros(1,n); f = 1e20*ones(1,n); aux = zeros(1,n); eqVal = zeros(1,n);
        end
    end
    methods (Test)
        function inMexPassesStridesOnceAndConverges(tc)
            rec = InMexResolveRecorder();
            cfg = tc.baseCfg(true);
            [s,l,u,d,f,a,e] = tc.batch();
            [~, fOut] = gdsge.runtime.solveProblems(@rec.mex, s,l,u,d,f,a,e, cfg);
            tc.verifyEqual(fOut, zeros(1,6), 'all problems should end converged');
            withStrides = find(~cellfun(@(c) c.stridesEmpty, rec.calls));
            tc.verifyNumElements(withStrides, 1, ...
                'exactly one MEX call should carry strides (the dedicated sweep)');
            sweep = rec.calls{withStrides};
            tc.verifyEqual(sweep.strides, [1 2], ...
                'strides must be cumprod([1, probSize(1:end-1)])');
            tc.verifyEqual(sweep.skip, ones(1,6), ...
                'the sweep call must skip Pass-0 (skip all 1)');
        end
        function matlabPathNeverPassesStrides(tc)
            rec = InMexResolveRecorder();
            cfg = tc.baseCfg(false);
            [s,l,u,d,f,a,e] = tc.batch();
            gdsge.runtime.solveProblems(@rec.mex, s,l,u,d,f,a,e, cfg);
            anyStrides = any(~cellfun(@(c) c.stridesEmpty, rec.calls));
            tc.verifyFalse(anyStrides, ...
                'useMexResolve=false must never set GDSGE_PROBLEM_STRIDES');
            tc.verifyGreaterThan(numel(rec.calls), 2, ...
                'the MATLAB sweep should issue multiple MEX calls');
        end
    end
end
```

- [ ] **Step 3: Run the test to confirm it fails**

Run:
```bash
matlab -batch "cd('tests'); addpath(pwd,'../src','../src/kernels'); r=runtests('runtime/tInMexResolveDriver.m'); disp(r); exit(double(any([r.Failed])))"
```
Expected: `inMexPassesStridesOnceAndConverges` FAILS — `solveProblems` does not yet read `cfg.useMexResolve`, so it never sets `GDSGE_PROBLEM_STRIDES` and the `evalin('caller','GDSGE_PROBLEM_STRIDES')` errors (undefined variable) or no call carries strides. (`matlabPathNeverPassesStrides` may already pass.)

- [ ] **Step 4: Edit `src/+gdsge/+runtime/solveProblems.m`** — three edits.

(a) In the doc comment `cfg:` list (line ~13), append `useMexResolve` to the documented fields:
```matlab
%        maxMinorIter probSize useNearestNeighbor verboseRetry adaptInSol asgHandle
%        useMexResolve
```

(b) In the contract block, after the `GDSGE_SPLINE_VEC = cfg.splineVec;` line (line ~31), add the strides contract local (always defined; empty by default so Pass-0 and restart calls do no in-MEX resolve):
```matlab
GDSGE_SPLINE_VEC = cfg.splineVec;             %#ok<NASGU>
GDSGE_PROBLEM_STRIDES = [];                   %#ok<NASGU>  % MEX in-resolve gate (set per sweep call)
useMexResolve = isfield(cfg,'useMexResolve') && cfg.useMexResolve;
if useMexResolve
    GDSGE_RESOLVE_STRIDES = cumprod([1, cfg.probSize(1:end-1)]);  % prod(probSize(1:d-1)), d=1..D
end
```

(c) Replace the existing `if cfg.useNearestNeighbor … end` block (today's lines 51–85) with the in-MEX branch in front of it. The full new block:
```matlab
    if useMexResolve
        % In-MEX neighbor sweep. Replicate the OLD inner-while ENTRY gate exactly
        % (numNeedResolvedAfter persists across outer iterations, init = inf), then
        % let the MEX run the entire progress loop internally in ONE call.
        needResolved = (f > cfg.tolSol) | isnan(f);
        numNeedResolved = sum(needResolved);
        if numNeedResolvedAfter ~= numNeedResolved
            skip(:) = 1;                                 % Pass 0 no-ops; strides drive the sweep
            GDSGE_PROBLEM_STRIDES = GDSGE_RESOLVE_STRIDES; %#ok<NASGU>  % enable in-MEX resolve
            [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);
            GDSGE_PROBLEM_STRIDES = [];                  %#ok<NASGU>  % off for the restart call below
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
```
Leave everything else untouched: `numNeedResolvedAfter = inf;` before the outer `while`, the outer `while` condition, the random-restart block, the `adaptInSol` block, and `diag`.

- [ ] **Step 5: Run the test to confirm it passes**

Run the Step-3 command.
Expected: both tests PASS.

- [ ] **Step 6: Run the runtime suite to confirm no regression**

Run:
```bash
matlab -batch "cd('tests'); addpath(pwd,'../src','../src/kernels'); r=runtests('runtime','IncludingSubfolders',true); disp(r); exit(double(any([r.Failed])))"
```
Expected: PASS. (No generated file sets `cfg.useMexResolve` yet, so `solveProblems` callers without the field take the `elseif useNearestNeighbor` branch — current behavior.)

- [ ] **Step 7: Commit**

```bash
git add src/+gdsge/+runtime/solveProblems.m tests/runtime/InMexResolveRecorder.m tests/runtime/tInMexResolveDriver.m
git commit -m "feat(runtime): solveProblems cfg.useMexResolve — dedicated in-MEX sweep call

Keeps the control flow byte-for-byte; when useMexResolve is on, replaces the
body of the MATLAB neighbor sweep with one strides-enabled MEX call, guarded
by the persistent numNeedResolvedAfter entry-gate so it fires at exactly the
iterations the MATLAB inner-while would. GDSGE_PROBLEM_STRIDES is a new
caller-workspace contract local (empty except during the sweep call). Random
restart, RNG draws, and convergence checks are unchanged. Fast unit test via a
recording fake MEX proves the call structure for both modes.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Codegen — turn the in-MEX path on by default

Wire the `UseMexResolve` runtime option (default 1) through the generated `iter_<model>.m`. After this task the in-MEX resolve is the default path for every spline/pchip model; the full slow suite must still match goldens on both backends.

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitSetup.m`
- Modify: `src/+gdsge/+codegen/+mat/optionsWhitelist.m`
- Modify: `src/+gdsge/+codegen/+mat/emitIter.m`
- Regenerate: `tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt`

- [ ] **Step 1: Edit `src/+gdsge/+codegen/+mat/emitSetup.m`** — add the default next to `UseAdaptiveBoundInSol`. After the line `w.add('UseAdaptiveBoundInSol = 0;');` (line ~42) insert:
```matlab
w.add('UseMexResolve = 1;');
```

- [ ] **Step 2: Edit `src/+gdsge/+codegen/+mat/optionsWhitelist.m`** — add `'UseMexResolve'` to the `base` list. Change the last line of the `base` cell (line ~12) from:
```matlab
    'CONSTRUCT_OUTPUT','NumThreads','WarmUp'};
```
to:
```matlab
    'CONSTRUCT_OUTPUT','NumThreads','WarmUp','UseMexResolve'};
```

- [ ] **Step 3: Edit `src/+gdsge/+codegen/+mat/emitIter.m`** — emit the cfg field. After the line `w.add('GDSGE_CFG.useNearestNeighbor = true;');` (line ~173) insert:
```matlab
w.add('GDSGE_CFG.useMexResolve = UseMexResolve;');
```

- [ ] **Step 4: Run the iter `.m` snapshot test to confirm drift**

Run:
```bash
matlab -batch "cd('tests'); addpath(pwd,'../src','../src/kernels'); r=runtests('HeatonLucas1996/codegen/tSnapshotHL1996.m'); disp(r); exit(double(any([r.Failed])))"
```
Expected: FAILS — the generated `iter_HL1996.m` now contains `UseMexResolve = 1;` and `GDSGE_CFG.useMexResolve = UseMexResolve;`.

- [ ] **Step 5: Regenerate snapshots and review**

Run:
```bash
matlab -batch "cd('tests/HeatonLucas1996/codegen'); addpath(fullfile(pwd,'..','..'),fullfile(pwd,'..','..','..','src'),fullfile(pwd,'..','..','..','src','kernels')); regen_snapshots"
```
Inspect `git diff tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt`: the only additions must be the `UseMexResolve = 1;` default line and the `GDSGE_CFG.useMexResolve = UseMexResolve;` cfg line. Confirm the cpp golden is unchanged (`git status`).

- [ ] **Step 6: Run the iter snapshot test to confirm it passes**

Run the Step-4 command. Expected: PASS.

- [ ] **Step 7: Run the FULL suite — the integration gate**

Every end-to-end model now runs the in-MEX resolve by default; all must still match their goldens.

Run:
```bash
matlab -batch "cd('tests'); run_tests"
```
Expected: exit 0. Decisive cases: `tEndToEndSafeAssets` (multi-root resolve — the real exercise of the sweep) and its SymPy sibling, plus all HL1996/Mendoza/GLSW/Cao2011EZ/CaoNie2016 gates and the ASG gates (CaoKS2016/Bianchi2011 — ASG passes no strides, must be untouched). If `safe_assets` diverges, STOP and use systematic-debugging: compare against the A/B harness in Task 4 to localize whether it is the C++ sweep or the driver gate.

- [ ] **Step 8: Commit**

```bash
git add src/+gdsge/+codegen/+mat/emitSetup.m src/+gdsge/+codegen/+mat/optionsWhitelist.m src/+gdsge/+codegen/+mat/emitIter.m tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt
git commit -m "feat(codegen): default UseMexResolve=1 — in-MEX resolve is the default path

emitSetup defaults UseMexResolve=1, optionsWhitelist accepts it, emitIter wires
GDSGE_CFG.useMexResolve. Generated iter_<model>.m now drives the in-MEX neighbor
sweep by default. Full suite green on both backends — all goldens match,
including safe_assets (multi-root) and the ASG models (unaffected).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: End-to-end A/B differential gate

Prove the in-MEX resolve produces **bit-identical** `IterRslt` to the MATLAB sweep, via the undocumented `UseMexResolve` toggle, on `safe_assets` (the decisive multi-root model) and HL1996 (sanity).

**Files:**
- Create: `tests/Barro_et_al_2017/codegen/tInMexResolveSafeAssets.m`

- [ ] **Step 1: Write the failing test `tests/Barro_et_al_2017/codegen/tInMexResolveSafeAssets.m`**

```matlab
classdef tInMexResolveSafeAssets < matlab.unittest.TestCase
    % A/B differential: in-MEX resolve (UseMexResolve=1) must produce a
    % bit-identical IterRslt to the MATLAB sweep (UseMexResolve=0). safe_assets
    % is the decisive case — it actually exercises the neighbor sweep + random
    % restart (clean models converge in Pass 0 and never enter the sweep).
    methods (Test, TestTags = {'Slow'})
        function inMexEqualsMatlabSweep(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'safe_assets.gmod'), work);

            gdsge_codegen('safe_assets');   % one compile; the toggle changes only the MATLAB path

            optsOn  = struct('SaveFreq', inf, 'NoSave', 1, 'UseMexResolve', 1);
            optsOff = struct('SaveFreq', inf, 'NoSave', 1, 'UseMexResolve', 0);
            R1 = iter_safe_assets(optsOn);
            R0 = iter_safe_assets(optsOff);

            tc.verifyEqual(R1.Iter, R0.Iter, 'iteration count diverged');
            tc.verifyEqual(R1.Metric, R0.Metric, 'final metric diverged');
            flds = {'var_policy','var_aux','var_interp'};
            for i = 1:numel(flds)
                a = R1.(flds{i}); b = R0.(flds{i});
                fn = fieldnames(a);
                for k = 1:numel(fn)
                    tc.verifyTrue(isequal(a.(fn{k}), b.(fn{k})), ...
                        sprintf('%s.%s not bit-identical between in-MEX and MATLAB resolve', ...
                        flds{i}, fn{k}));
                end
            end
        end
    end
end
```

- [ ] **Step 2: Run it to confirm it passes**

(After Task 3 the toggle and the C++ resolve both exist, so this should pass immediately — it is a regression lock, not red-then-green. If it FAILS, that is a real bit-exactness defect; debug before proceeding.)

Run:
```bash
matlab -batch "cd('tests'); addpath(pwd,'../src','../src/kernels'); r=runtests('Barro_et_al_2017/codegen/tInMexResolveSafeAssets.m'); disp(r); exit(double(any([r.Failed])))"
```
Expected: PASS — `Iter`, `Metric`, and every `var_policy`/`var_aux`/`var_interp` field bit-identical.

- [ ] **Step 3: Write the HL1996 sanity A/B `tests/HeatonLucas1996/codegen/tInMexResolveHL1996.m`**

HL1996 converges in Pass 0 (never enters the sweep), so this is the "the toggle changes nothing on clean models" check.

```matlab
classdef tInMexResolveHL1996 < matlab.unittest.TestCase
    % A/B differential: in-MEX resolve (UseMexResolve=1) vs MATLAB sweep
    % (UseMexResolve=0) must produce a bit-identical IterRslt. HL1996 converges
    % in Pass 0 and never enters the sweep, so this is the clean-model sanity.
    methods (Test, TestTags = {'Slow'})
        function inMexEqualsMatlabSweep(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'HL1996.gmod'), work);

            gdsge_codegen('HL1996');

            optsOn  = struct('SaveFreq', inf, 'NoSave', 1, 'UseMexResolve', 1);
            optsOff = struct('SaveFreq', inf, 'NoSave', 1, 'UseMexResolve', 0);
            R1 = iter_HL1996(optsOn);
            R0 = iter_HL1996(optsOff);

            tc.verifyEqual(R1.Iter, R0.Iter, 'iteration count diverged');
            tc.verifyEqual(R1.Metric, R0.Metric, 'final metric diverged');
            flds = {'var_policy','var_aux','var_interp'};
            for i = 1:numel(flds)
                a = R1.(flds{i}); b = R0.(flds{i});
                fn = fieldnames(a);
                for k = 1:numel(fn)
                    tc.verifyTrue(isequal(a.(fn{k}), b.(fn{k})), ...
                        sprintf('%s.%s not bit-identical between in-MEX and MATLAB resolve', ...
                        flds{i}, fn{k}));
                end
            end
        end
    end
end
```

Verify the HL1996 gmod filename first: the model dir is `tests/HeatonLucas1996/`; confirm the gmod is `HL1996.gmod` (it is, per the corpus) before relying on the `copyfile`.

- [ ] **Step 4: Run the HL1996 A/B to confirm it passes**

```bash
matlab -batch "cd('tests'); addpath(pwd,'../src','../src/kernels'); r=runtests('HeatonLucas1996/codegen/tInMexResolveHL1996.m'); disp(r); exit(double(any([r.Failed])))"
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add tests/Barro_et_al_2017/codegen/tInMexResolveSafeAssets.m tests/HeatonLucas1996/codegen/tInMexResolveHL1996.m
git commit -m "test: A/B differential gate — in-MEX resolve == MATLAB sweep (bit-identical)

iter_<model>(UseMexResolve=1) vs (=0) produce isequal IterRslt for safe_assets
(multi-root, exercises sweep+restart) and HL1996 (clean, sanity). Locks the
bit-exact parity the design promises.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Round-trip characterization (informational) + PROGRESS

Demonstrate the perf win (fewer MEX round-trips) and record the work. Not a hard gate.

**Files:**
- Create: `tests/perf/inmex_resolve_bench.m`
- Modify: `PROGRESS.md`

- [ ] **Step 1: Write the benchmark script `tests/perf/inmex_resolve_bench.m`**

A plain script (not a unittest) that compiles `safe_assets` once, then times `iter_safe_assets` under `UseMexResolve` 1 vs 0 and prints iteration time. Run manually; record the numbers in PROGRESS.

```matlab
function inmex_resolve_bench()
% Manual benchmark: in-MEX resolve vs MATLAB sweep on safe_assets.
% Run from repo root: matlab -batch "cd('tests'); addpath(pwd,'../src','../src/kernels'); cd('perf'); inmex_resolve_bench"
here = fileparts(mfilename('fullpath'));
work = tempname; mkdir(work);
copyfile(fullfile(here, '..', 'Barro_et_al_2017', 'safe_assets.gmod'), work);
oldDir = cd(work); cleanup = onCleanup(@() cd(oldDir));
gdsge_codegen('safe_assets');
base = struct('SaveFreq', inf, 'NoSave', 1, 'NoPrint', 1);
for useMex = [1 0]
    opts = base; opts.UseMexResolve = useMex;
    t = tic; R = iter_safe_assets(opts); el = toc(t);
    fprintf('UseMexResolve=%d : Iter=%d  Metric=%.3e  time=%.2fs\n', ...
        useMex, R.Iter, R.Metric, el);
end
end
```

- [ ] **Step 2: Run the benchmark and capture output**

```bash
matlab -batch "cd('tests'); addpath(pwd,'../src','../src/kernels'); cd('perf'); inmex_resolve_bench"
```
Expected: both rows show the same `Iter`/`Metric`; the `UseMexResolve=1` row should have a lower `time` (fewer MEX round-trips). Record both lines.

- [ ] **Step 3: Update `PROGRESS.md`** — add a changelog entry at the top of the `## Changelog` section summarizing the feature (in-MEX nearby-warmup resolve, default on, bit-exact A/B verified, the safe_assets timing from Step 2), and note that the Phase-10 `safe_assets` perf follow-up is addressed. Mention the spec/plan paths and the branch `inmex-nearby-resolve`.

- [ ] **Step 4: Commit**

```bash
git add tests/perf/inmex_resolve_bench.m PROGRESS.md
git commit -m "docs(perf): in-MEX resolve benchmark + PROGRESS changelog

safe_assets benchmark script (UseMexResolve 1 vs 0, same Iter/Metric, fewer MEX
round-trips on the in-MEX path); PROGRESS changelog for the feature.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final verification

- [ ] Run the full suite once more and confirm exit 0:
```bash
matlab -batch "cd('tests'); run_tests"
```
- [ ] Confirm both backends are covered: `tEndToEndSafeAssets` + `tEndToEndSafeAssetsSympy` green (the SymPy one requires the `uv` env).
- [ ] Confirm the ASG gates (`tEndToEndCaoKS2016`, `tEndToEndBianchi2011`) are green — they prove the no-strides path is untouched.
- [ ] `git log --oneline` shows the five focused commits on branch `inmex-nearby-resolve`.
- [ ] Hand off via `superpowers:finishing-a-development-branch` (merge / PR decision).
