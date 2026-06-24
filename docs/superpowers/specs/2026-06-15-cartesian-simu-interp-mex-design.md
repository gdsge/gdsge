# Design — in-MEX cartesian SIMU_INTERP simulation + myppual retirement

**Date:** 2026-06-15
**Status:** approved (pending spec review)
**Anchors:** GLSW2020 `GLSW_interp` (unprimed transition), Mendoza2010 (primed transition)

## Problem

The cartesian (spline) `SIMU_INTERP=1` simulation path (`emitSimulateInterp.m`) runs the
period loop in **MATLAB**, calling `myppual_mex` once per period to evaluate
`IterRslt.output_interp`, then doing `var_simu` / state-transition assignments in MATLAB.
For long horizons (e.g. Mendoza `num_periods=50000`) this is 10k–50k MATLAB↔MEX
round-trips plus per-period array churn.

The hans toolbox's `simulate_stochastic_mex.cpp` shows the fast pattern: a MEX that takes
the constructed `pp` struct and the preallocated result struct, loops samples under OpenMP,
**searches the interpolant once per site** and `nosearch_eval`s every output, writing
results (and next-period states) **in place** into the struct.

This work ports that pattern to the refactored cartesian SIMU_INTERP path and, in the same
pass, **retires `myppual` entirely**: `myppual.m`, `myppual_mex.mexw64`,
`convert_to_interp_eval_array.m`, and `get_scalar.m` are deleted, with zero `myppual`
references remaining.

ASG SIMU_INTERP already evaluates in-MEX via `asg.eval_vec`; it is **out of scope** and
must stay green.

## Where `myppual` lives today (the five roles)

| Role | Site | Use | Replacement |
|------|------|-----|-------------|
| R1 | `emitResultIter.m:40` | construct `IterRslt.output_interp` | `interp_construct_mex` |
| R2 | `emitSimulate.m` (SIMU_RESOLVE) | construct + per-period `myppual_mex` eval of warm-start `GDSGE_PP` | `interp_construct_mex` + `interp_eval_mex` |
| R3 | `emitSimulateInterp.m` | per-period `myppual_mex` eval | `simulate_<model>_mex` (whole loop) |
| R4 | `applyWarmUp.m` | construct + eval warm-up interpolants | `interp_construct_mex` + `interp_eval_mex` |
| R5 | `emitResultIter.m` (`get_scalar`), `constructSplines.m` legacy branch (`convert_to_interp_eval_array`) | reshape helper / eval-array build | plain `reshape` / delete dead branch |

Construction already has a bit-exact drop-in (`interp_construct_mex`, produces `== myppual(pp)`).
Evaluation has **no** MATLAB-callable replacement except `myppual_mex`, so R2/R4 force the
creation of `interp_eval_mex`.

## Interpolant layout (the key enabler)

`IterRslt.output_interp` is rebuilt in the **stacked uniform-order** layout the model MEX and
hans already use, instead of today's mixed-order `[2, OutputInterpOrder…]` (linear in shock,
cubic in states):

- **Today:** shock is a spline *dimension* (order 2); states order `OutputInterpOrder`. Mixed
  order — only `myppual_mex` can evaluate it.
- **New:** shock is the **stacked vector index**, a single uniform cubic spline over **states
  only**. For a realized shock `s`, output component `c` is vector index
  `(s-1)*nComp + c`, where `nComp` = total output rows (sum of output-var lengths). Built by
  `interp_construct_mex` (one value field, `numArray = shock_num*nComp`), read by the existing
  `MatlabPp<nStates,OutputInterpOrder>` — no new evaluator, no `myppual_mex.cpp` source.

This is **numerically identical at integer shocks**: each shock slice is the same per-shock
cubic state-spline, and a not-a-knot spline reproduces node values exactly in the shock
direction regardless of its order there. `output_interp` is internal (only the simulate path
reads it; no golden compares its coefs), so the layout change is safe; `SimuRslt` is
unchanged. The SIMU_RESOLVE warm-start interpolant (`GDSGE_PP`) and `applyWarmUp` move to the
same stacked uniform-order layout.

## Architecture

Two new C++ MEX functions sharing one evaluation core (the existing
`include/InterpEval.h` / `include/MatlabPp.h` headers — the same uniform-order evaluator the
model MEX and hans use), plus construction routed onto `interp_construct_mex`.

```
gen_discrete_markov_rn (MATLAB, RNG parity)
        │  shock path
        ▼
simulate_<model>.m ── preallocate SimuRslt, init, chunk loop ──► simulate_<model>_mex (whole period loop, in place)
                                                                       │ uses
IterRslt.output_interp (built by interp_construct_mex) ────────────────┘ MatlabPp eval core

interp_eval_mex (generic) ◄── SIMU_RESOLVE warm-start, applyWarmUp, SIMU_INTERP fallback
```

### Component 1 — `simulate_<model>_mex` (codegen'd, per model)

The fast path: the entire SIMU_INTERP period loop in C++, written in place.

**Inputs**
- `output_interp` — the `pp` struct (`breaks`, `coefs`, `pieces`, `order`, `dim`).
- shock matrix — `num_samples × (num_periods+1)`, **pre-generated in MATLAB** (preserves
  `gen_discrete_markov_rn` RNG parity; see the safe_assets RNG-restart note).
- the preallocated `SimuRslt` field arrays (state columns + `var_simu` columns), written in
  place.
- period range `[t0, t1]` for chunking (see Component 3).
- `num_threads`, `EnforceSimuStateInbound`, `shock_num`, `SimuPrintFreq`.

Output-variable index ranges (`output_var_index`) and the `var_simu` / transition→output
mappings are **baked into the generated C++** from the IR slot layout at codegen time.

**Per period `t = t0..t1`, `#pragma omp parallel for` over samples `i`:**
1. If `EnforceSimuStateInbound==1`, clamp each current state `states[j][i]` to its grid
   endpoints (`breaks` min/max — equals `[min(var_state), max(var_state)]`, matching
   `emitSimulateInterp.m:63-69`).
2. Build the site from **states only**: `xSite = [states]` (shock is the vector index, not a
   spline dim).
3. `pp.search(xSite, cellOfSite, xSiteToLeft)` — **once** per sample.
4. Let `base = (shock[i,t]-1)*nComp`. For each needed output row `r`,
   `nosearch_eval_<order>(xSiteToLeft, cellOfSite, base + r)`.
5. Record `var_simu`: `SimuRslt.<name>(i,t) = nosearch_eval(base + <name's row>)` (var_simu
   rows are length 1).
6. Advance states (`t < num_periods`):
   - unprimed transition `s' = outvar` → `next_states[s][i] = nosearch_eval(base + <outvar's
     single row>)`.
   - primed transition `s' = outvar'` → `outvar` occupies `shock_num` consecutive output
     rows; pick the realized next-shock row:
     `next_states[s][i] = nosearch_eval(base + <outvar's row 0> + (shock[i,t+1]-1))`
     (the C++ equivalent of `outvar(GDSGE_SHOCK_VAR_LINEAR_INDEX)`).
7. Print every `SimuPrintFreq` periods via `mexPrintf` (period number + sample-1 field
   values + elapsed), matching the current MATLAB print block.

Row offsets (`nComp`, each output var's row range, primed/unprimed flags) are baked into the
generated C++ from the IR slot layout.

Writes directly into the caller's `SimuRslt` arrays (hans-style: `mxGetPr` of each field +
period offset).

### Component 2 — `interp_eval_mex` (generic, prebuilt once)

A thin MATLAB-callable evaluator over the `MatlabPp`/`InterpEval` core, for the
**stacked uniform-order** interpolants (states-only spline; shock = vector index). It accepts
a `pp` (or `ppCell{1}`) from `interp_construct_mex`, the per-sample state sites, and a
per-sample vector-index array (`(shock-1)*nComp + (0:nComp-1)`), and returns the evaluated
rows:

```
vals = interp_eval_mex(int32(numThreads), pp, sites, int32(idx))
```

`pp` carries `breaks`/`coefs`/`pieces`/`order`/`dim`; `idx` selects the stacked vectors per
site (empty = all). Compiled once by `ensureInterpEvalMex` (mirrors `ensureSplineConstructMex`:
`/openmp`, `-DUSE_OMP`, `-Iinclude`, MSVC default `/fp:precise`). Used by:
- SIMU_RESOLVE warm-start eval (`emitSimulate`),
- `applyWarmUp` eval,
- the SIMU_INTERP backward-compat fallback (Component 6).

### Component 3 — `simulate_<model>.m` (rewritten `emitSimulateInterp`)

Pure-MATLAB scaffolding around one or a few MEX calls:
1. setup / pull `output_interp`, `shock_trans`, `shock_num` from `IterRslt`;
2. preallocate `SimuRslt` (shock, states, var_simu);
3. initial conditions + `GDSGE_OPTIONS.init` overwrite; initial-shock-bound check;
4. `SimuRslt.shock(:) = gen_discrete_markov_rn(...)` — RNG stays in MATLAB;
5. **chunked** MEX calls: loop in chunks of `SimuSaveFreq` periods (one call when saving is
   off), calling `simulate_<model>_mex` with `[t0,t1]`; `save(...)` between chunks. The MEX
   prints every `SimuPrintFreq` internally.

No `myppual` anywhere.

### Component 4 — construction / evaluation swaps (the rest of myppual)

- `emitResultIter`: `IterRslt.output_interp` built via `interp_construct_mex` in the stacked
  uniform-order layout (state breaks only; one value field of `shock_num*nComp` stacked rows);
  `get_scalar(GDSGE_VAR_POLICY, ...)` / `get_scalar(GDSGE_VAR_AUX, ...)` → plain `reshape` to
  `{[1:shock_num], states}` shape (the only thing `get_scalar` did).
- `emitSimulate` (SIMU_RESOLVE): build `GDSGE_PP` via `interp_construct_mex` in the same
  stacked uniform-order layout; per-period `myppual_mex` eval → `interp_eval_mex` with
  `idx = (shock-1)*nSolComp + (0:nSolComp-1)`. (The macOS `'pp'`/`myppual(myppual())` branch
  collapses to the same construct path.)
- `applyWarmUp`: construct via `interp_construct_mex`, eval via `interp_eval_mex` (state-only
  uniform layout).
- `constructSplines`: delete the dead legacy `myppual` + `convert_to_interp_eval_array`
  branch and the `useFusedConstruct` toggle (fused is already the only live path).

### Component 5 — build & codegen plumbing

`gdsge.codegen.generateMatlab` emits `simulate_<model>_mex.cpp`, appends it to
`compile_<model>.m`, and gives it a `.cache` keyed like the model MEX (so a second
`gdsge_codegen` run skips recompilation). `ensureInterpEvalMex` builds the generic evaluator
once. `ensurePath` / `Contents.m` / doc comments updated to drop `myppual` mentions.

### Component 6 — backward-compat safety net

Supported simulate blocks: `var_simu` entries resolve to output vars (always true via
`resolveOutputs`), and each transition RHS is a **single output-var reference** (primed or
unprimed). This covers the entire corpus (GLSW `a1'=a1n`; Mendoza `cTilde'=cTildeNext'`,
`k'=kNext`; CaoKS `K'=Kp`, `X'=Xp`).

If codegen encounters a compound transition RHS (arithmetic, indexing), it **falls back** to
a MATLAB period loop that calls `interp_eval_mex` and does the transition in MATLAB — still
myppual-free, so no existing `.gmod` breaks (honoring the hard backward-compat constraint).
No current corpus model triggers the fallback.

## Deletions

`src/kernels/myppual.m`, `src/kernels/myppual_mex.mexw64`,
`src/kernels/convert_to_interp_eval_array.m`, `src/kernels/get_scalar.m`. Final gate asserts
no `myppual` / `convert_to_interp_eval_array` / `get_scalar` references remain in `src/`.

## Testing (TDD)

- **Parity goldens (write failing first):** new-MEX `SimuRslt` vs existing goldens for
  GLSW_interp (unprimed) and Mendoza2010 (primed), within the models' existing tolerances.
- **Differential check (dev-time):** new-MEX `SimuRslt` vs the old MATLAB+`myppual_mex`
  `SimuRslt`, within tight tolerance, before deleting `myppual`. The interpolants agree at
  integer shocks (same per-shock state spline), so `SimuRslt` should match to round-off; pin
  `rng` so the shock paths are identical.
- **Regression:** CaoKS2016_simu_interp (ASG) stays green; HL1996 / safe_assets /
  CaoKS2016 (resolve + warmup paths now on `interp_construct_mex` / `interp_eval_mex`) stay
  green.
- **Build:** second `gdsge_codegen` run does not recompile `simulate_<model>_mex` (cache).
- **Final:** zero `myppual` references in `src/`.

Numerics: pin `rng` before any iter that feeds these gates (safe_assets-style), and run
MATLAB processes **sequentially** (one-MATLAB-at-a-time).

## Non-goals

- ASG simulation (already in-MEX).
- Changing `IterRslt` / `SimuRslt` struct shapes (frozen).
- Generating the shock path inside the MEX (kept in MATLAB for RNG parity).
- Optimizing the SIMU_RESOLVE solve loop itself (only its interpolation calls move off
  myppual).
