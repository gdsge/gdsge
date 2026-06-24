# In-MEX randomized restart (`UseMexRandomize`) — design

**Date:** 2026-06-15
**Status:** approved (brainstorming → plan)
**Branch:** `inmex-randomize`

## 1. Overview & goal

Move the random-restart trial loop — and the outer-iteration-1 initial guess — out of
MATLAB and into the C++ MEX, driven by a deterministic, fixed-seedable C++ RNG. The MEX runs
up to `MexRandomizeBatch` (default **100**) full minor-iterations internally, returning to
MATLAB only at batch boundaries so MATLAB keeps ownership of **diagnostics**, the
**`MaxMinorIter` cap**, and the **converge / error / exit** decision.

Today (cartesian path) each "minor iteration" inside `gdsge.runtime.solveProblems` costs two
MATLAB↔MEX round-trips: an in-MEX neighbor-resolve fixpoint (already in C++ via
`UseMexResolve=1`) **then** one MATLAB random restart (`x0Rand = rand(...)` → one MEX call).
On restart-heavy models this dominates: **safe_assets** needs >1000 restarts on its first
outer iteration. Batching collapses *two round-trips per trial* into *one MEX call per ~100
trials*.

This deliberately **breaks byte-identity** with the current goldens (a different — C++,
deterministic — random sequence picks different restart paths). That is accepted: after a
tolerance-equivalence check against the old goldens, the new output becomes the new
bit-exact benchmark.

### Scope

- **In scope:** the **iter** cartesian (spline/pchip) solve path —
  `gdsge.runtime.solveProblems` + `templates/cxx/task.tpl.cpp`, driven by `emitIter`.
- **Out of scope (keep today's MATLAB randomize path):** `simulate` (SIMU_RESOLVE) and the
  `model_init` task. Their per-call resolves converge in ~0 restarts from the warm start, so
  they carry none of the cost; keeping them on the MATLAB path leaves their goldens
  untouched. They share `solveProblems`/`task.tpl.cpp`, so they can adopt the mechanism
  later for free.
- **Out of scope:** the ASG path (`solveProblemsAsg` / ASG task), mirroring how the in-MEX
  neighbor-resolve kept ASG out of scope. Tracked as a deferred follow-up.

## 2. Architecture & data flow

The mechanism lives in the two files the in-MEX neighbor-resolve already touches.

### 2.1 `templates/cxx/task.tpl.cpp`

1. Refactor the existing neighbor-resolve block (current lines 122–171) into a reusable
   lambda `GDSGE_neighbor_fixpoint()` so both the standalone in-MEX-resolve call and the new
   batched loop invoke the identical sweep. (`GDSGE_solve_one` is already a lambda.)
2. Add a batched randomize loop, active only when a new caller-workspace input
   `GDSGE_RANDOMIZE_BATCH > 0` is present:

```
Pass 0 (existing): solve each non-skip point from its supplied guess
if GDSGE_RANDOMIZE_BATCH > 0:
    GDSGE_neighbor_fixpoint()                           // propagate warm-started solutions
    for r in 0 .. BATCH-1:
        if all converged: break
        #pragma omp parallel for                        // unconverged points only (f > TolSol)
            GDSGE_SOL(:,i) = counter-RNG(seed, salt, globalPointIdx_i, trialOffset + r)
            GDSGE_solve_one(i)
        GDSGE_neighbor_fixpoint()
```

The batch calls pass `skip(:) = 1`, so the MEX's Pass 0 is a pass-through copy of the current
state (no wasted solves); the restart loop selects points by residual (`f > TolSol`), not by
`skip`. The separate first MEX call (before the loop) is the real Pass 0 from the supplied
guess.

The exact ordering (initial neighbor sweep, then per-trial *restart → solve → neighbor
sweep*) mirrors the spirit of the current cascade (neighbor before random; each random
restart followed by a neighbor re-sweep) so convergence behavior stays comparable and the
tolerance gate passes. Byte-parity is **not** a goal, so the ordering may be tuned during
TDD against the tolerance gate.

3. The MEX returns the **same 5 outputs** (`sol, f, aux, eqVal, optInfo`). **No extra output
   is needed for trial accounting:** the loop breaks only when *every* point is converged, so
   if MATLAB sees any unconverged point after a call, the full `BATCH` was consumed; if all
   converged, the exact count is moot (the run is done). MATLAB therefore advances its trial
   bookkeeping by `BATCH` whenever the returned state is still unconverged.

### 2.2 `src/+gdsge/+runtime/solveProblems.m`

When `cfg.useMexRandomize` is set (and `cfg.adaptInSol` is empty and cartesian strides are
present), the `while` body becomes one MEX call per batch:

```
f(:) = 1e20; skip(:) = 0;
[sol,f,...] = mexFn(...)                                % Pass 0 from supplied warm start
minorIter = 0; trialOffset = 0;
while (still unconverged) && minorIter < maxMinorIter
    batch = min(randomizeBatch, maxMinorIter - minorIter)   % cap-aware (inf => batch = N)
    set caller-workspace: GDSGE_RANDOMIZE_BATCH = batch,
        GDSGE_RANDOM_SEED, GDSGE_RANDOMIZE_SALT, GDSGE_RANDOMIZE_TRIAL_OFFSET = trialOffset,
        GDSGE_PROBLEM_STRIDES = strides; skip(:) = 1
    [sol,f,aux,eqVal,optInfo] = mexFn(...)              % runs up to `batch` minor-iters in C++
    if converged: break
    minorIter  += batch
    trialOffset += batch
    diagnostics heartbeat (now at batch granularity)
end
% final-diagnostics + error(iter)/warn(simulate) block: UNCHANGED
```

The MEX caller-workspace contract (variables read via `mexGetVariablePtr`/`mexGetVariable`)
gains the new `GDSGE_RANDOMIZE_*` names alongside the existing `GDSGE_PROBLEM_STRIDES`. As
with the existing contract, all `mexFn(...)` calls and these locals stay in the body of
`solveProblems` (never a subfunction).

The `useMexRandomize=false` branch keeps today's code path **verbatim** (in-MEX/MATLAB
neighbor sweep + one MATLAB `rand` restart per minor iter) so it stays the A/B reference.

## 3. C++ RNG design (the "fix seed" core)

**Counter-based, per-grid-point deterministic.** A single shared (or per-thread) generator
would make the random sequence depend on OpenMP scheduling, so a "fixed seed" would not
reproduce a run. Instead every draw is a pure function of integer coordinates:

> `guess(point i, component c, trial t)` = a uniform on `[lb_c, ub_c]` derived by hashing
> `(MexRandomSeed, salt, globalPointIndex_i, t, c)`.

Implementation: mix the integer key with splitmix64 (and/or seed a per-draw `std::mt19937_64`
from it), convert 64 random bits to a double in `[0,1)` via the standard `(bits >> 11) *
2^-53`, then `guess = lb_c + u * (ub_c - lb_c)` — matching the semantics of the old
`rand .* (ub-lb) + lb`. Exact primitive/constants are an implementation detail pinned in TDD;
the **semantics** below are the contract:

- **Thread-count independent.** No shared RNG state → identical results at any `NumThreads`,
  so goldens are portable across machines (important: `NumThreads` defaults to
  `feature('numcores')`).
- **Batch-resumable.** `trialOffset` is an integer MATLAB passes; batch *k+1* continues the
  stream seamlessly with no carried RNG state.
- **Fresh per outer iteration.** `salt = GDSGE_Iter`, so each outer iteration draws
  independent (still reproducible) streams — matching the old global `rand` never repeating
  draws across iterations.
- **Restart trials** are indexed `0,1,2,…`, continued across batches by `trialOffset`. Outer
  iterations ≥2 warm-start from the carried-over converged solution (no randomness), exactly
  as today.

Default seed is **fixed** (`MexRandomSeed = 0`) so the default run is reproducible out of the
box (more deterministic than the old unpinned MATLAB `rand`). Goldens are captured under this
fixed default.

### Initial-guess handling

To honor "move the first guess to C++" — i.e. *MATLAB `rng` no longer affects iter results;
reproducibility depends on the C++ seed alone* — the generated `iter_<model>.m` replaces the
unconditional `GDSGE_X0 = rand(...)` with a runtime branch on `UseMexRandomize` (see §4):
when on, `GDSGE_SOL` is seeded to the **deterministic bound midpoint** `(GDSGE_LB+GDSGE_UB)/2`
(no MATLAB `rand`); the first *random* guess for any point that doesn't converge from the
midpoint is restart trial 0, drawn in the MEX from `MexRandomSeed`. When off, today's MATLAB
random initial guess is used, byte-identical to current behavior. Outer iterations ≥2 reuse
the carried-over converged solution as the warm start, exactly as today.

## 4. Codegen & options

- **`gdsge.codegen.mat.optionsWhitelist`** — add `UseMexRandomize`, `MexRandomizeBatch`,
  `MexRandomSeed` to the frozen option surface.
- **`gdsge.codegen.mat.emitSetup`** — defaults: `UseMexRandomize = 1;`,
  `MexRandomizeBatch = 100;`, `MexRandomSeed = 0;`.
- **`gdsge.codegen.mat.emitIter`** — thread the new cfg fields into `GDSGE_CFG`
  (`useMexRandomize`, `randomizeBatch = MexRandomizeBatch`, `randomSeed = MexRandomSeed`,
  `randomizeSalt = GDSGE_Iter`). Replace the unconditional `GDSGE_X0 = rand(...)` init with a
  runtime branch:
  ```matlab
  if UseMexRandomize
      GDSGE_SOL(:) = (GDSGE_LB + GDSGE_UB) / 2;   % deterministic midpoint; restarts drawn in MEX
  else
      GDSGE_X0 = rand(size(GDSGE_SOL)) .* (GDSGE_UB-GDSGE_LB) + GDSGE_LB;
      GDSGE_SOL(:) = GDSGE_X0;
  end
  ```
- **`emitSimulate` / `emitIterInit`** — unchanged (leave `useMexRandomize` unset → MATLAB
  randomize path). Out of scope.
- **Incompatibility guard:** if `UseAdaptiveBoundInSol == 1`, `solveProblems` falls back to
  the MATLAB randomize path (in-MEX batching cannot adjust bounds mid-loop). No corpus model
  sets this by default.

## 5. Backward-compat, defaults, goldens

- **Default-on** (`UseMexRandomize = 1`). The `UseMexRandomize = 0` branch preserves today's
  exact MATLAB randomize behavior; set it to reproduce the *old* `.mat` goldens bit-exactly.
- **Re-baseline procedure:**
  1. Tolerance-equivalence gate: the new default output lands on the same equilibrium as the
     old goldens within tolerance (§6).
  2. Regenerate the **iter** goldens — the functional `IterRslt*.mat` (re-captured under the
     fixed default seed) **and** the changed `iter_*_golden.txt` text snapshots.
- **Simulate goldens are untouched:** simulate keeps the MATLAB path, and shock draws via
  `rng(0823)` are independent of the solve randomization regardless.
- IR schema unchanged (only new public options + generated-file text).

## 6. Validation / tests

- **Determinism gate** — drive the HL1996 MEX directly (via
  `tests/+gdsgetest/buildHL1996MexInputs.m`) with a fixed `BATCH` + `MexRandomSeed` at
  `NumThreads = 1` vs many → assert bit-identical `sol` (proves per-point / thread
  independence).
- **Driver unit test** (à la `tInMexResolveDriver`, nested fake MEX) — asserts the call
  structure: batch sizing (cap-aware), `trialOffset` advancing across batches,
  `salt = GDSGE_Iter`, `skip(:) = 1` on batch calls, and the early-exit-on-converge
  trial-count inference.
- **Tolerance-equivalence gates** — HL1996 + safe_assets (+ Mendoza / GLSW) iter under
  `UseMexRandomize = 1` vs the OLD goldens, within tolerance (see §7 for safe_assets).
- **A/B gate** — `UseMexRandomize = 0` still reproduces the OLD goldens bit-exactly (the
  MATLAB path is unchanged at runtime).
- Full suite stays green; restart-heavy gates keep their `rng`/`MaxMinorIter=inf` handling
  where they currently need it (now also pinning `MexRandomSeed`).

## 7. Risks

- **safe_assets multi-root points.** A different restart path can select a different (still
  valid) root, so pointwise policy may differ beyond tight tolerance at those points.
  Mitigation: the tolerance gate checks pointwise where the root is unique and falls back to
  economic / simulated-moment agreement at multi-root points; since the new output is a valid
  equilibrium, it becomes the benchmark. Owner (domain expert) signs off on the re-baseline.
- **Diagnostic granularity** coarsens from `DiagnoseMinorIter` (10) to the batch size (100)
  when the feature is on — the heartbeat fires at each batch boundary. Documented; with the
  default `MexRandomizeBatch = 100` being a multiple of `DiagnoseMinorIter = 10`, the existing
  `mod(minorIter, diagnoseAt) == 0` check fires every boundary unchanged.
- **`UseAdaptiveBoundInSol = 1`** is incompatible with in-MEX batching → falls back to the
  MATLAB randomize path (§4).
- **Cap semantics** are now enforced at batch granularity: `minorIter` advances by `BATCH`,
  so the effective cap rounds up to a multiple of the batch. Acceptable; restart-heavy gates
  use `MaxMinorIter = inf`.

## 8. Alternatives considered (rejected)

- **Global / per-thread RNG** — not reproducible under OpenMP; defeats "fix seed".
- **Batch only the restarts, keep the neighbor sweep MATLAB-driven** — still round-trips
  MATLAB per trial; forgoes most of the speedup.
- **Opt-in (default-off)** — forgoes the speed win by default; the owner wants speed by
  default and accepts re-baselining.
- **Nondeterministic default seed** — rejected in favor of a fixed default for out-of-the-box
  reproducibility.

## 9. Non-goals

ASG path; simulate / init in-MEX randomize; changing the CoDoSol solver or the
neighbor-resolve math; any IR schema change.
