# Bounded retries + tiered convergence diagnostics — design

**Date:** 2026-06-14
**Status:** approved (brainstorming)
**Scope:** cartesian (spline/pchip) solve path in `gdsge.runtime.solveProblems`

## Problem

When a model fails to converge, the inner resolve loop in
`gdsge.runtime.solveProblems` retries random initial guesses **without bound**
(`MaxMinorIter` defaults to `inf`). On a genuinely non-converging model this
loop spins forever, and the only output is a bare per-retry line. The user gets
no actionable signal about *why* it is failing.

There are two nested loops; this design touches only the inner one:

1. **Outer (major) loop** in the generated `iter_<model>.m`: `while(~stopFlag)`,
   stops at `GDSGE_Metric < TolEq` or `GDSGE_Iter >= MaxIter`. **Out of scope.**
2. **Inner (resolve) loop** in `solveProblems`: for grid points whose residual
   `GDSGE_F` exceeds `TolSol`, runs nearest-neighbor warm start → **random
   restart**, looping `while (max(f) > tolSol || any isnan(f)) && minorIter <
   MaxMinorIter`. **This is what we bound and instrument.**

## Goal

Bound the retry budget and emit **tiered, named diagnostics** so a
non-converging cartesian solve stops with an actionable report instead of
running forever. The diagnostics point at the three most common root causes:
bounds too tight, NaN in the model, and a mis-normalized/buggy equation.

Behavior is centralized in `solveProblems`, so it covers both `iter` and
`simulate`'s resolve path, and both Jacobian backends (autodiff + SymPy).

## Decisions (locked during brainstorming)

- **On exhaustion: hard error on the `iter` path; warn-and-continue on
  `simulate`.** When the retry budget is exhausted with points still
  unconverged, both paths print the full diagnostic block. `iter` then
  `error(...)`, halting the whole run; `simulate` warns for that period and
  continues to the next (now *bounded* by `MaxMinorIter` instead of looping
  forever). This is a `cfg.errorOnNonconvergence` flag: `iter` sets it `true`,
  `simulate` sets it `false`.
- **Schedule:** `MaxMinorIter` default `inf → 200`; a diagnostic **heartbeat**
  prints every `DiagnoseMinorIter` retries (default `10`) while still stuck; the
  full block + error fire at the cap. All tunable. *(Revised during
  implementation from the original `cap=20, one-shot at 10`: the corpus showed
  restart-heavy models need a much larger budget, and a recurring heartbeat is
  more useful than a single print across a 200-retry budget — see Risks.)*
- **Report detail:** named per-policy-component bounds breakdown + NaN + dominant
  equation + the existing worst-K grid-point table.
- **Architecture:** in `gdsge.runtime.solveProblems`, with a new
  `gdsge.runtime.diagnoseConvergence` formatter; codegen threads the labels and
  schedule into `GDSGE_CFG`.
- **Scope:** cartesian path first. ASG (`solveProblemsAsg`) is a follow-up.

## Retry schedule (per outer iteration's solve)

The inner loop does one random restart per `minorIter` (occasionally two when the
`UseAdaptiveBoundInSol` hook is active). New behavior:

| When | Action |
|------|--------|
| Every retry (existing) | Light line `resolve round k: N unconverged, worst residual g`, gated on `~NoPrint`. |
| Every `DiagnoseMinorIter` retries (default 10) while stuck — `mod(minorIter, DiagnoseMinorIter) == 0` and `minorIter < MaxMinorIter` | Print the **compact summary** heartbeat (`level = 'summary'`), then keep retrying. |
| Loop exits still unconverged at `MaxMinorIter` (default 200) | Print the **full** block (`level = 'full'`), then `error('gdsge:runtime:solveNotConverged', …)` if `cfg.errorOnNonconvergence` (iter), else return `diag` and let the caller warn (simulate). |

A converging solve exits the loop before retry 10, so it produces **no new
output and never errors** — behavior is byte-identical to today for every model
in the corpus.

The heartbeat uses `mod(minorIter, DiagnoseMinorIter) == 0` and is skipped at
`minorIter == MaxMinorIter` (the final block covers the cap iteration), so there
is no duplicate at exhaustion. Setting `DiagnoseMinorIter = inf` suppresses the
heartbeat entirely (`mod(k, inf) = k ≠ 0`), which the restart-heavy corpus gates
use together with `MaxMinorIter = inf` to restore the old quiet, unbounded
regime.

## Diagnostic content

Computed from `sol`, `lb`, `ub`, `f` (1×nProb worst-residual), `eqVal`
(nEq×nProb per-equation residuals), and the threaded labels. Two populations:

- **unconverged:** `f > TolSol | isnan(f)`
- **converged:** the complement

Signals:

1. **NaN** — count of problems with `isnan(f)`. Flagged prominently; usually a
   `log` of a negative or a divide-by-zero in the model body.
2. **Bounds-pinned (unconverged, per named component)** — for each solution row
   `r` (named `solNames{r}`), among unconverged points, count those with `sol`
   within `pinTol·span` of `lb` ("at LB") or `ub` ("at UB"), where
   `span = max(1, ub-lb)`, `pinTol = 1e-6`. Hints that bounds are too tight and
   are blocking convergence. Reported sorted by count, e.g.
   `w1n(3): 45/120 at UB`.
3. **Near-bound (converged, per named component)** — among converged points,
   count those within `nearBand·span` of a bound, `nearBand = 0.01` (1%). Hints
   an active/binding bound even where the solve succeeded. e.g.
   `w1n(3): 30 near UB`.
4. **Dominant equation** — over unconverged points, rank equation rows of
   `eqVal` by `max(abs(eqVal(e, unconverged)))`; report the top `topEq = 3` rows
   with the ratio of the top value to the **median row's** value. A ratio `>> 1`
   indicates one equation dominates the residual — a bug or a scaling/
   normalization issue. e.g. `eq#7: 1.4e+02 (142× median) <-- dominant`.
5. **Worst-K points** — reuse `gdsge.runtime.reportUnconverged`: top
   `topK = 5` unconverged grid points by `f`, with shock + state coordinates.

**Tier split:**
- `summary` (heartbeat, every `DiagnoseMinorIter` retries): NaN count + top
  pinned components + the dominant equation — the headlines.
- `full` (at the cap, before the error): all five signals in full tables.

## Naming / labels

- **Policy components** are named from the IR: for each `ir.variables.policy{i}`
  with `.name` and `.slot = [lo hi]`, row `r` in `[lo,hi]` is labeled
  `name` (scalar var, `lo==hi`) or `name(r-lo+1)` (array var, e.g. `w1n(3)`).
  Built once at codegen as a length-`maxDim` cell `solNames` and threaded via
  `cfg`.
- **Equations** have no user-given names in the IR (each `equations` item is a
  residual `expr` + a `primed` flag). They are labeled by **index**: `eq#k`
  for the k-th equation row, `eq#k[shock j]` for a primed expansion. Index
  labels are honest and still actionable ("equation #7 is 142× the median —
  check its scale"). AST→text rendering of the equation is an *optional* nicety,
  not required for this design; equations are reported by index in the runtime,
  so no equation-label array needs threading.

## Thresholds

`pinTol = 1e-6`, `nearBand = 0.01`, `topK = 5`, `topEq = 3` are **internal
constants** in `diagnoseConvergence`, deliberately kept off the public option
surface (YAGNI). Only the schedule numbers (`MaxMinorIter`, `DiagnoseMinorIter`)
are options.

## Architecture

### New: `gdsge.runtime.diagnoseConvergence`

```
function msg = diagnoseConvergence(sol, lb, ub, f, eqVal, cfg, level)
```

A **pure formatter** that returns the report string for `level ∈
{'summary','full'}`. Reads `cfg.solNames`, `cfg.probSize`, `cfg.tolSol`, and —
**only when present** — `cfg.stateNames` / `cfg.stateGrids`. The `full` level
appends the `reportUnconverged` point table **only when
`isfield(cfg,'stateGrids') && ~isempty(cfg.stateGrids)`** (iter supplies them;
simulate does not, because its "problems" are samples, not a state grid — see
note below). No side effects (the caller decides `fprintf` vs `error`).
Absent/empty `cfg.solNames` falls back to `x1, x2, …` component labels so the
formatter never errors on a minimally-populated `cfg`.

### Changed: `gdsge.runtime.solveProblems`

All new behavior is **gated on `isfield(cfg,'diagnoseAt')`**, so callers that do
not supply the diagnostic fields (the bare `tSolveProblems` unit tests using
`baseCfg`) keep their exact current behavior: cap → return `diag`, no printing,
no throw.

- Inside the inner `while`: after `f` is updated, a recurring heartbeat block:
  ```
  if isfield(cfg,'diagnoseAt') && cfg.diagnoseAt > 0 ...
          && mod(minorIter, cfg.diagnoseAt) == 0 && minorIter < cfg.maxMinorIter ...
          && (max(isnan(f)) || max(f(:)) > cfg.tolSol)
      fprintf('[gdsge] convergence diagnostics (retry %d, continuing):\n%s\n', ...
          minorIter, gdsge.runtime.diagnoseConvergence(sol,lb,ub,f,eqVal,cfg,'summary'));
  end
  ```
- After the loop, the exhaustion check (only when diagnostics are enabled):
  ```
  if isfield(cfg,'diagnoseAt') && (max(isnan(f)) || max(f(:)) > cfg.tolSol)
      fprintf('%s\n', gdsge.runtime.diagnoseConvergence(sol,lb,ub,f,eqVal,cfg,'full'));
      if isfield(cfg,'errorOnNonconvergence') && cfg.errorOnNonconvergence
          error('gdsge:runtime:solveNotConverged', ...
              'Solve did not converge: %d/%d points after %d retries. See diagnostics above.', ...
              nnz((f>cfg.tolSol)|isnan(f)), numel(f), minorIter);
      end
  end
  ```
- The light per-retry `verboseRetry` line is kept.
- The existing `diag` struct return is unchanged (still populated for the
  caller). On the `simulate` path (`errorOnNonconvergence=false`) the caller's
  existing per-period `reportUnconverged` warning still fires.

### Changed: codegen threading

`emitIter` builds `GDSGE_CFG` once before the loop. Add:
- `GDSGE_CFG.solNames = {…};` — the length-`maxDim` component-name cell,
  built from the policy slot layout (a new small helper, e.g.
  `gdsge.codegen.mat.solComponentNames(ir)`).
- `GDSGE_CFG.stateNames = {…};` and `GDSGE_CFG.stateGrids = {…};` — already
  available in `emitIter` as `stateNameCell` / `stateGridCell`.
- `GDSGE_CFG.diagnoseAt = DiagnoseMinorIter;`
- `GDSGE_CFG.errorOnNonconvergence = true;` (iter). `emitSimulate` sets it
  `false`.
- `maxMinorIter`, `probSize`, `tolSol` already on `cfg`.

The existing `reportUnconverged` call currently emitted **after**
`solveProblems` returns in the generated iter file moves into the formatter
(`full` level). The post-return `NeedResolved` warning block in `emitIter` is
removed (its job is now done inside `solveProblems`; a successful solve has no
leftover unconverged points because exhaustion now errors).

`emitSimulate`'s solve-cfg gets `solNames` + `diagnoseAt` +
`errorOnNonconvergence = false`, but **not** `stateNames`/`stateGrids` — so the
formatter prints the per-component bounds / NaN / dominant-equation signals
(well-defined over samples) and skips the grid-point table. Simulate keeps its
existing per-period `reportUnconverged` warning untouched as the period locator.

## Option surface

- `MaxMinorIter` — existing public option; default literal in `emitSetup`
  changes `inf → 200`. Set to `inf` to restore the old unbounded behavior (the
  restart-heavy corpus gates do this).
- `DiagnoseMinorIter` — **new** option, default `10`. The heartbeat interval (a
  summary prints every `DiagnoseMinorIter` retries while stuck); set to `inf` to
  silence the heartbeat. Added to `emitSetup` (default),
  `gdsge.codegen.mat.optionsWhitelist` (so it is accepted), and threaded into
  `cfg`.
- Both flow through the existing `gdsge.runtime.unpackOptions` whitelist; an
  unknown field still errors as today.

## Error semantics

On exhaustion, the full block is written with `fprintf` to stdout (so it is
visible even if the caller wraps the solve in `try/catch`). Then, **on the iter
path** (`cfg.errorOnNonconvergence = true`),
`error('gdsge:runtime:solveNotConverged', …)` is thrown; the outer iter loop
does **not** catch it, so it propagates to the user and halts the run. **On the
simulate path** (`cfg.errorOnNonconvergence = false`), no error is thrown — the
block is printed, `solveProblems` returns, and the generated simulate file's
existing per-period `reportUnconverged` warning fires; the simulation continues
to the next period. The error message names the unconverged count and retry
budget and points at the printed diagnostics. (Optional refinement: thread
`cfg.outerIter = GDSGE_Iter` for a sharper "at outer iteration I" message.)

## Backward compatibility

- **Converging models:** the inner loop exits before retry 10, so none of the
  new code runs and no error is raised. Functional `IterRslt`/`SimuRslt` goldens
  are **bit-exact**.
- **Default change `MaxMinorIter inf → 200`** only changes behavior for solves
  that fail to converge within 200 retries — previously these hung; now they
  error with diagnostics. This is the intended improvement. Restart-heavy models
  (safe_assets, Mendoza-sympy) set `MaxMinorIter = inf` to keep the old regime.
- **Generated-file text** changes (new `cfg` fields, new `DiagnoseMinorIter`
  default, removed post-return warning block), so all `iter_<model>.m` and
  `simulate_<model>.m` **text snapshots** regenerate via the existing
  `regen_snapshots.m`. No functional golden changes.

## Testing (TDD — failing test first)

1. **`tDiagnoseConvergence`** (unit, pure, fast) — synthetic `sol/lb/ub/f/eqVal`
   with: one component pinned at UB on several unconverged points, one NaN
   residual, and one equation row inflated ~100×. Assert the `full` report names
   the pinned component (`w1n(3)`-style), states the NaN count, and flags the
   inflated row as the dominant equation with a plausible ratio. Assert the
   `summary` report is the headline subset.
2. **`tSolveProblemsCap`** (driver, fast) — a nested fake `mexFn` that never
   drives one point below `tolSol`, with a full `cfg` including `diagnoseAt` and
   `errorOnNonconvergence=true`. Assert `solveProblems` throws
   `gdsge:runtime:solveNotConverged` after exactly `MaxMinorIter` retries, and
   that the heartbeat fired at the `DiagnoseMinorIter` interval (capture via
   `evalc`). A `recurringDiagnosticsEveryInterval` case asserts the heartbeat
   prints once per interval (retries 2,4,6 for `diagnoseAt=2`, cap 7). A second
   case with `errorOnNonconvergence=false` asserts it instead returns `diag`
   (all `needResolved`) after printing the full block (the simulate path).
   Precedent: `tInMexResolveDriver`.
3. **Existing `tSolveProblems` stays green** — its `baseCfg` supplies no
   `diagnoseAt`, so `maxMinorIterCapsAndReportsDiag` still observes the
   cap → return path unchanged (no printing, no throw). This is the regression
   guard for the `isfield(cfg,'diagnoseAt')` gate.
4. **Backward-compat** — full corpus gates stay green; regenerate the text
   snapshots and confirm functional goldens are unchanged. Run via
   `matlab -batch "cd('tests'); run_tests"`.

## Risks

- **Restart-heavy corpus models need an unbounded budget (resolved during
  implementation).** **safe_assets** (multi-root grid points) and the
  **Mendoza2010 sympy** gate genuinely need a heavy-tailed, RNG-dependent number
  of random restarts — the old default was `MaxMinorIter = inf`, and measurement
  showed the golden-reproducing seed (`rng('default')`) needs **>1000** restarts
  on the very first outer iteration of safe_assets. So *no* fixed cap with a hard
  error is safe for these models across their ~1271 outer iterations — it would
  fire probabilistically. **Resolution:** the default cap was raised to `200`
  (comfortable for normal models, which converge in 1–3 restarts), and the 5
  restart-heavy gates (`tEndToEndSafeAssets`, `tEndToEndSafeAssetsSympy`,
  `tFusedConstructSafeAssets`, `tInMexResolveSafeAssets`,
  `tEndToEndMendoza2010Sympy`) pass `MaxMinorIter = inf` + `DiagnoseMinorIter =
  inf` (and `rng('default')` where they did not already pin), restoring the exact
  quiet, unbounded regime their goldens were captured under. Real users of such
  models raise `MaxMinorIter` themselves — the thrown diagnostics tell them to.
  Full suite **482/482** green after this resolution.
- **Snapshot churn** — every `iter_*`/`simulate_*` text snapshot regenerates.
  Use `regen_snapshots.m`; verify the diff is confined to the expected `cfg` /
  defaults / removed-warning-block lines.

## Out of scope

- **ASG path** (`solveProblemsAsg`) — sparse-grid points have no fixed cartesian
  `LB`/`UB` layout; a diagnostics variant is a separate follow-up.
- **Outer-loop `MaxIter = inf` stalls** — when every inner solve succeeds but the
  metric never reaches `TolEq`. That loop already prints progress every
  `PrintFreq`, so it is not silent; a finite `MaxIter` default or stall detection
  is a separate enhancement.
- **Exposing `pinTol`/`nearBand`/`topK`/`topEq` as options** — internal constants
  for now.
