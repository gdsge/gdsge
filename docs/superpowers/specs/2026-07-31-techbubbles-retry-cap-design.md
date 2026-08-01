# Design: uncapped resolve retries, sample-path diagnostics, and an install self-check

**Date:** 2026-07-31
**Status:** Approved
**Motivating model:** `base_package/gdsge/tests/TechBubblesV4nblogKRAv3/TechBubblesV4nblogKRAv3.gmod`
(runs on the old toolbox; `simulate` produces garbage on the new one).

## 1. Problem and goal

A two-state, two-shock technology-bubble model that the old toolbox solves and simulates
correctly. On the new toolbox `iter` is already right — it converges at **Iter=355,
Metric=9.93e-07** against the old toolbox's **Iter=355, Metric=9.75e-07**, with
`var_policy` agreeing to 3.9e-08 and identical `GDSGE_PROB` bounds. Parsing, codegen, and
the MEX are all correct.

`simulate` is not. Every period leaves 6–77 of 100 samples unconverged and the paths walk
off: `q` reaches 987 against a true maximum of 2.5, `qd` reaches 2.8e+04 against 10.3.

Investigating that surfaced three more defects, plus one unrelated failure reported from a
macOS install. **Goal:** `TechBubblesV4nblogKRAv3` runs end-to-end on the new toolbox and
matches the old toolbox's results, with the failure modes uncovered along the way closed
and regression-covered.

**Non-goals:** no IR schema change; no new gmod syntax; no change to the solver itself.

## 2. Root causes

### 2.1 A — the resolve retry cap (the reason simulate fails)

`src/+gdsge/+codegen/+mat/emitSetup.m:23` emits `MaxMinorIter = 200;` with the comment
*"bounded retries; was inf (could loop forever)"*. The old toolbox emits
`MaxMinorIter = inf` and its simulate retry loop is unbounded — it does not even reference
the variable.

This model needs that. Instrumenting the old toolbox's retry counter over 20 periods:

```
OLDPERIOD 1: firstbad=100 retries=4920
OLDPERIOD 3: firstbad=94  retries=5253
OLDPERIOD 4: firstbad=89  retries=5120
OLDPERIOD 8: firstbad=16  retries=3787
```

Thousands of randomized restarts per period is normal operation here, not a pathology. The
cap stops at 200 and — because the simulate path sets `errorOnNonconvergence=false` — the
unconverged solutions are written into `SimuRslt` and fed forward as the next period's
state. The cap does not fail loudly; it silently truncates results.

Confirmed by overriding the option alone, no source change, `MaxMinorIter=inf`:

| field | max abs diff vs old |
|---|---|
| shock | 0 (bit-exact) |
| logK | 8.4e-11 |
| b | 1.3e-15 |
| q | 6.6e-11 |
| qn | 5.7e-11 |
| qd | 9.1e-10 |
| ER / r | 9.7e-11 |
| IK | 3.7e-12 |

Cost is at parity with the old toolbox: 44.7 s (old) vs 52.5 s (new) for 20 periods × 100
samples.

The same default feeds the iter path (where `errorOnNonconvergence=true` turns exhaustion
into a hard error) and both ASG paths via `solveProblemsAsg`.

### 2.2 B — `reportUnconverged` crashes on the simulate path

`src/+gdsge/+runtime/reportUnconverged.m:23` indexes `subs{s+1}` once per state name.
`emitSimulate.m:157` passes `probSize = [1 GDSGE_NPROB]` — two dimensions, where the
second is the sample index — together with the model's two state names. `ind2sub` yields
two subscripts, the loop asks for `subs{3}`, and the call dies with
`MATLAB:badsubscript`. The simulate "warn and continue" path crashes whenever it fires.

### 2.3 C — unbounded diagnostics heartbeat

`solveProblems.m:93` and `:178` print a full `diagnoseConvergence(...,'summary')` block
whenever `mod(minorIter, cfg.diagnoseAt) == 0`. With `DiagnoseMinorIter = 10` and an
uncapped loop that is ~530 blocks per period, ×10,000 periods.

The `mod` test is also wrong for the in-MEX randomize branch, which advances `minorIter` by
`actualTrials` (up to `MexRandomizeBatch`) per call: a step that jumps across a multiple of
`diagnoseAt` without landing on it skips the trigger entirely.

### 2.4 D — a corrupt source file gives an unreadable error

Reported from a macOS install:

```
Error using gdsge.runtime.unpackOptions
Execution of script unpackOptions as a function is not supported:
.../src/+gdsge/+runtime/unpackOptions.m
```

Reproduced exactly (`MATLAB:scriptNotAFunction`, message verbatim) with an **empty or
comment-only** `unpackOptions.m`. MATLAB finds the file, sees no `function` line,
classifies it as a script, and refuses to call it with arguments.

Ruled out: filename case (a case-only mismatch gives `MATLAB:undefinedVarOrClass`,
"Unable to resolve the name"), UTF-8 BOM (no error), and CR-only line endings (no error).
Our sources are clean — the file is intact in `D:\refactor_gdsge` and `D:\gdsge_publish`,
correctly cased in every commit of both histories, and no file under `src/` has a
filename/function-name mismatch. This was a corrupt working copy on one machine; the defect
is that the toolbox reports it unintelligibly.

## 3. Fix A — restore the uncapped default

`emitSetup.m:23` → `MaxMinorIter = inf;`, with the comment recording why the cap was
removed: a finite default silently truncates results on models that legitimately need
thousands of restarts, which is worse than the runaway loop it was guarding against.

`errorOnNonconvergence = true` on the iter path stays as-is. It now fires only when a user
sets a finite `MaxMinorIter`, which is the correct contract: opting into a bound means
opting into being told when it is hit.

**Golden impact: none.**

- No corpus model exceeds 200 retries today (the suite is green with the cap in place).
- The in-MEX randomize batch is `min(cfg.randomizeBatch, cfg.maxMinorIter - minorIter)`;
  with `maxMinorIter = inf` this is 100, exactly as it is under the cap for the first two
  batches. The RNG stream is unchanged.
- `MaxMinorIter` is a generated-file default, not an IR option — **no IR schema change, and
  no existing model's IR snapshot changes**. The `adding-an-IR-option` fan-out does not
  apply. (§7 adds one *new* IR snapshot, for the new corpus model only.)
- Only `tests/HeatonLucas1996/codegen/golden/{iter,simulate}_HL1996_golden.txt` change,
  regenerated with the existing `tests/HeatonLucas1996/codegen/regen_snapshots.m`.
- `optionsWhitelist.m` already lists `MaxMinorIter` and `DiagnoseMinorIter`; unchanged.

## 4. Fix B — a sibling reporter for the sample path

`reportUnconverged` is grid-shaped by construction: it maps a linear problem index through
`ind2sub(probSize, ...)` onto shock index plus per-dimension grid coordinates. In simulate
the problems are samples, not grid points, and `shock`, `logK`, `b` are already per-sample
column vectors in scope at the call site.

Rather than overload one function with two coordinate systems, add:

```matlab
gdsge.runtime.reportUnconvergedSamples(needResolved, f, shockVec, stateNames, stateVals, topK)
%   3 of 100 problems unconverged; worst 3:
%     [sample=17, shock=2, logK=1.234, b=-0.31] residual=1.02
```

indexing `shockVec` and each `stateVals{s}` directly by the problem index. `emitSimulate.m`
calls this instead. `reportUnconverged` is untouched, and both stay small and
independently testable.

`emitSimulateAsg.m` does not emit an unconverged warning today, so it needs no change.

## 5. Fix C — geometric backoff, and an aggregate for simulate

### 5.1 Backoff

In `solveProblems`, replace the `mod(minorIter, cfg.diagnoseAt) == 0` test in both branches
with a `nextDiagnose` threshold initialised to `cfg.diagnoseAt` and doubled after each
print, advancing past the current `minorIter` so a large batch step cannot leave it behind:

```matlab
if diagnoseAt > 0 && minorIter >= nextDiagnose && minorIter < cfg.maxMinorIter && stillStuck
    fprintf(...);
    while nextDiagnose <= minorIter; nextDiagnose = nextDiagnose * 2; end
end
```

Triggers at 10, 20, 40, 80, … — 9 blocks for a 5,300-retry solve instead of 530,
deterministic, and threshold-crossing also fixes the batch-skipping bug in §2.3.

### 5.2 Simulate goes quiet by default

`emitSetup(ir, 'simulate')` emits `DiagnoseMinorIter = 0;` (iter keeps `10`). That matches
old-toolbox silence for a model where heavy retrying is normal. The `diagnoseAt` field is
still emitted, so the `isfield`-gated FINAL diagnostics block still prints on cap
exhaustion, and a user can set `DiagnoseMinorIter` to opt back in.

### 5.3 A better replacement: per-checkpoint retry aggregate

The `retryStatsInit` / `retryStatsAdd` / `retryStatsEndIter` helpers are already generic
counters — "iter" is only a label — so the simulate period loop reuses them verbatim:
`retryStatsInit` before the loop, `retryStatsAdd(..., GDSGE_DIAG.minorIters)` after
`solveProblems`, `retryStatsEndIter(..., GDSGE_t)` at the end of each period.

A new `gdsge.runtime.printSimuRetryStats(stats)` prints one line at the existing
`SimuPrintFreq` checkpoint and returns a fresh accumulator, mirroring what
`printIterProgress` does for iter:

```
  resolve: 312/1000 periods retried, 1243905 rounds total, worst period 3 (5253 rounds)
```

Generated files stay thin over a unit-tested runtime library, per the Phase 4 principle.

## 6. Fix D — install self-check

New `gdsge.runtime.checkInstall()`, called first by `ensurePath()` — already the single
choke point ("called first by `gdsge_codegen` and by every generated
`iter_<model>`/`simulate_<model>`"). A `persistent` flag limits it to once per MATLAB
session; the scan is ~70 small files, a few milliseconds.

For every `.m` under `src/+gdsge/**`, plus `src/*.m` and `src/kernels/*.m`, it requires the
first non-comment, non-blank token to be `function` or `classdef`. Files named
`Contents.m` are excluded — they are comment-only by design. Failures are collected and
raised **once**:

```
gdsge:install:corrupt
GDSGE install looks incomplete or corrupted — these files define no function:
  src/+gdsge/+runtime/unpackOptions.m
Restore them (e.g. `git checkout -- src/`) or re-clone the toolbox.
```

Fatal rather than a warning: the toolbox cannot run in this state, and a clear message is
the entire point. It cannot protect against `checkInstall.m` or `ensurePath.m` themselves
being corrupt; that limitation is accepted.

## 7. Corpus gate

`tests/TechBubblesV4nblogKRAv3/`, following the Phase 9a/9b precedent:

- `TechBubblesV4nblogKRAv3.gmod` (copied from `base_package/gdsge/tests/`)
- `ir/` IR JSON snapshot, `codegen/` end-to-end test
- `golden/{IterRslt,SimuRslt}.mat` captured from the **old** toolbox
- `tests/golden/capture_TechBubblesV4nblogKRAv3.m`
- `tGoldenTechBubblesV4nblogKRAv3.m`

**Simulate size: 6 samples × 100 periods.** Measured with `MaxMinorIter=inf`: 6×20 = 17.7 s,
6×100 = 21.0 s, 20×20 = 28.7 s. Cost is dominated by the hard early periods, so 6×100 buys
much better coverage for ~3 s. With iter at ~8 s the gate is ~30 s, in line with the other
Slow gates.

**Tolerances, not bit-exactness.** This model has multi-root grid points, so the randomized
restart is order-sensitive — the `safe_assets` / `CaoNie2016` fragility. `rng(0823)` is
pinned before iter *and* before simulate. Old-vs-new `var_policy` agrees to 3.9e-08 and
simulate fields to ≤9.1e-10 given the same `IterRslt`, so the gate compares at ~1e-6
relative with the shock path required bit-exact.

This is the only test in the suite that exercises a >200-retry solve; without it, Fix A can
silently regress.

## 8. Verification

Test-first, per the repo's TDD constraint:

1. `reportUnconvergedSamples` — shape and indexing, including the 2-state simulate case
   that crashes today.
2. Backoff schedule — triggers at 10/20/40/80, and a batch step that jumps across several
   thresholds fires once and resumes correctly (the §2.3 bug).
3. `checkInstall` — fires on an empty file and on a comment-only non-`Contents` file, does
   *not* fire on `Contents.m` or a healthy tree, and reports every bad file in one error.
4. Emitted-default regression: generated `iter_*`/`simulate_*` carry `MaxMinorIter = inf`,
   and simulate carries `DiagnoseMinorIter = 0`.
5. The new end-to-end gate.
6. Full suite green (600 tests as of the last run); all existing goldens byte-identical
   except the two regenerated HL1996 codegen snapshots.

## 9. Out of scope

**In-MEX randomize for the simulate path.** `UseMexRandomize` is wired only for cartesian
iter (`emitIter.m:185-188`); simulate falls back to the MATLAB loop, one MEX call per
retry. Batching restarts inside the MEX would likely cut this model's ~7 h full run
(10,000 periods × 100 samples at ~2.6 s/period) substantially. It is performance work
beyond restoring parity and is cleanly separable — its own cycle if wanted.

**Parser warnings on this gmod.** Four warnings fire (`"lAI" is defined in the "var_shock"
block but is associated with "parameters"`, likewise `ql`, `qh`, and `shock_trans` under
`var_state`) because the model assigns those names after a later declaration keyword.
Values flow correctly — `IterRslt.params.lAI` is the intended 0.0569 — and the old toolbox
is silent. Left alone rather than folding a style-nag change into this fix.
