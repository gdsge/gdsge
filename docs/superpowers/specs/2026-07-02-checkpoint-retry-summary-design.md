# Checkpoint retry summary — design

**Date:** 2026-07-02
**Status:** approved

## Problem

When an iteration has unconverged grid points, `solveProblems` /
`solveProblemsAsg` print a resolve summary line (`[iter N] resolve converged
after M rounds`) at the end of every major iteration that needed retries. On
models where most iterations retry (e.g. Mendoza), this doubles the console
output and buries the `PrintFreq` checkpoint lines.

## Goal

Record retry activity per major iteration and print it **compactly at the
PrintFreq checkpoint** instead of per iteration. Non-convergence behavior
(diagnostics + error/warning) is unchanged.

## Decisions (from brainstorming)

- Keep the throttled in-solve round heartbeat (every `ResolvePrintFreq`
  rounds) as a liveness signal for a long-stuck solve.
- Summary format: **one extra indented line** under the checkpoint line.
- Scope: **both iter paths** (cartesian `solveProblems`, ASG
  `solveProblemsAsg`). Simulate already runs `verboseRetry=false` — untouched.
- **No new gmod option** (avoids the IR-option fan-out). Mid-solve visibility
  remains available by lowering the existing `ResolvePrintFreq`.
- Implementation: **runtime accumulator threaded through `printIterProgress`**
  (pure stats struct; no handle state, no logic in generated string templates).

## User-visible behavior

- The per-iteration `[iter N] resolve converged after M rounds` lines are gone
  from both iter paths.
- At each checkpoint (`mod(Iter, PrintFreq)==0` or the final iteration), if
  any covered iteration needed retries, one indented line prints:

  ```
  Iter:80, Metric:0.00123, maxF:8.2e-09, unconverged:0, elapsed:5.1s
    resolve: 7/10 iters retried, 134 rounds total, worst iter 74 (58 rounds)
  ```

  If no iteration since the last printed checkpoint retried, nothing extra
  prints.
- **Kept exactly as today:** the round heartbeat (every `ResolvePrintFreq`
  rounds inside one solve), the `DiagnoseMinorIter` diagnostics heartbeat, and
  the whole exhaustion path — the `stopped at MaxMinorIter` line, FINAL
  diagnostics, `error` on the iter path / warn-and-continue on simulate.
  `NoPrint=1` still silences everything.
- Model-init solves (cartesian init and ASG init refinement) have no
  checkpoint; their converged-retry chatter simply disappears. A stuck init
  still surfaces via the round heartbeat and the exhaustion path.

## Components

1. **`src/+gdsge/+runtime/retryStatsAdd.m`** —
   `stats = retryStatsAdd(stats, minorIters)`. Passing `[]` initializes
   (`itersSeen`, `itersRetried`, `totalRounds`, `worstIter`, `worstRounds`,
   `curRounds` all zero). Adds `minorIters` into the current iteration's
   running bucket (`curRounds`); ASG calls this once per refinement level.
2. **`src/+gdsge/+runtime/retryStatsEndIter.m`** —
   `stats = retryStatsEndIter(stats, iter)`. Called once per major iteration
   after all solves: increments `itersSeen`; if `curRounds > 0` increments
   `itersRetried`, adds into `totalRounds`, updates `worstIter`/`worstRounds`;
   clears `curRounds`.
3. **`printIterProgress.m`** — optional 9th argument and 2nd return:
   `[printed, stats] = printIterProgress(iter, metric, maxF, nUnconverged,
   elapsedSec, printFreq, noPrint, stopFlag, stats)`. When the checkpoint line
   prints and `stats.itersRetried > 0`, prints the `resolve:` line
   (`itersRetried/itersSeen iters retried, totalRounds rounds total, worst
   iter worstIter (worstRounds rounds)`). Whenever the checkpoint line prints
   (with or without retries), the returned stats are reset so `itersSeen`
   counts only the window since the last printed checkpoint. Existing 8-arg
   callers are unchanged.
4. **`solveProblems.m` / `solveProblemsAsg.m`** — the end-of-loop
   `printResolveProgress('summary', ...)` call fires only when the cap was
   hit; the converged path prints nothing. `printResolveProgress.m` drops the
   now-dead "converged after N rounds" wording and updates its header comment.
5. **`emitIter.m` / `emitIterAsg.m`** — three emitted lines each:
   - `GDSGE_RETRY_STATS = [];` before the main loop;
   - `GDSGE_RETRY_STATS = gdsge.runtime.retryStatsAdd(GDSGE_RETRY_STATS,
     GDSGE_DIAG.minorIters);` after each in-loop solve call (inside the ASG
     refinement loop);
   - at the checkpoint: `GDSGE_RETRY_STATS = gdsge.runtime.retryStatsEndIter(
     GDSGE_RETRY_STATS, GDSGE_Iter);` then the two-output `printIterProgress`
     call (timer reset keeps riding the `printed` return).
   Init segments (`emitIterInit`, ASG init) are untouched.
6. **Codegen goldens** regenerate (text change only; no IR/schema/option
   change, so no model-JSON snapshot churn).

## Edge cases

- **Denominator correctness:** `itersSeen` counts major iterations, not solve
  calls — the Add/EndIter split makes `7/10` right even when ASG solves
  several refinement levels per iteration.
- **Skipped checkpoints:** stats accumulate until a checkpoint actually
  prints; the reset rides the existing `printed` return, same as the elapsed
  timer.
- **Non-convergence:** unchanged and immediate — `solveProblems` still prints
  FINAL diagnostics and errors mid-iteration on `MaxMinorIter` exhaustion;
  pending stats are superseded by the error.
- **Frozen surfaces:** no gmod option, no IR change, no `IterRslt`/`SimuRslt`
  shape change; simulate path untouched; `printIterProgress` stays
  backward compatible for bare 8-arg callers.

## Testing (TDD)

- Unit tests (no MEX) for `retryStatsAdd`/`retryStatsEndIter`: init from `[]`,
  multiple adds within one iteration, worst-iteration tracking, zero-retry
  iterations leave `itersRetried` at 0.
- Unit tests for `printIterProgress` with stats: summary line prints only at a
  freq boundary with retries; stats reset after printing and accumulate across
  non-printing iterations; silent under `NoPrint`; 8-arg call still works.
- `evalc` test on `solveProblems` with a fake MEX (pre-seed
  `GDSGE_RANDOMIZE_TRIALS_USED`, per the known fake-MEX gotcha): a converged
  retry run prints no summary line; a capped run still prints the
  `stopped at MaxMinorIter` line.
- Regenerate codegen goldens; full suite
  `matlab -batch "cd('tests'); run_tests"` green; one end-to-end run with
  retries to eyeball the new checkpoint output.
