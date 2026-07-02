# Design — throttled resolve-round printing + end-of-loop summary

Date: 2026-06-16
Status: Approved — partially superseded by
`2026-07-02-checkpoint-retry-summary-design.md` (the converged
"resolve converged after N rounds" summary line described below no longer
prints; converged retries are now summarized at the PrintFreq checkpoint)

## Problem

During the solve/iter phase, the resolve loop prints a progress line on **every**
minor iteration whenever printing is enabled. This is noisy. Two changes:

1. Throttle the per-round line to print only every `ResolvePrintFreq` rounds
   (default **100**).
2. Print a single resolve **summary** when the loop ends — either because all
   points converged (report the total round count) or because it hit the
   `MaxMinorIter` cap (report how many points remain unconverged).

## Affected code sites

The per-round resolve print exists at three sites, with inconsistent behavior today:

- `src/+gdsge/+runtime/solveProblems.m:159-162` (cartesian, legacy resolve loop) —
  prints `resolve round %d` every minor iter, gated by `cfg.verboseRetry`.
- `src/+gdsge/+runtime/solveProblems.m` `useMexRandomize` branch (the default fast
  path) — prints **no** per-round line, and advances `minorIter` in **batches**
  (`+batch` per MEX call) rather than `+1`.
- `src/+gdsge/+runtime/solveProblemsAsg.m:103-106` — prints `asg resolve round %d`
  every minor iter, same gating.

`cfg.verboseRetry` is set to `~NoPrint` in the iter emitters and to `false` in the
simulate emitters.

## Design

### New option: `ResolvePrintFreq`

- Default **100**.
- Wired through the same chain as the existing `PrintFreq` option:
  - `src/+gdsge/+parser/defaultSetupCode.m` — add `ResolvePrintFreq = 100;`
  - `src/+gdsge/+parser/resolveOptions.m` — `o.resolvePrintFreq = getf(ws, 'ResolvePrintFreq', 100);`
  - `src/+gdsge/+ir/schema.m` — add `resolvePrintFreq` (`opt(fScalar())`).
  - `src/+gdsge/+codegen/+mat/emitSetup.m` — emit `ResolvePrintFreq = <value>;`
  - `src/+gdsge/+codegen/+mat/optionsWhitelist.m` — add `'ResolvePrintFreq'`.
  - `src/+gdsge/+codegen/+mat/emitIterInit.m`, `emitIter.m`, `emitIterAsg.m` —
    set `GDSGE_CFG.resolvePrintFreq = ResolvePrintFreq;` alongside the existing
    `GDSGE_CFG.verboseRetry = ~NoPrint;`.

### Gating

- `cfg.verboseRetry` remains the **master on/off** switch (`~NoPrint`). When false,
  nothing prints — per-round lines or summary.
- The simulate paths (`emitSimulate.m`, `emitSimulateAsg.m`) already set
  `verboseRetry = false`, so they stay silent and need no change. `ResolvePrintFreq`
  affects only the iter/solve paths.

### Behavior — applied uniformly to all three sites

1. **Throttled per-round line.** Print only when a `ResolvePrintFreq` window
   boundary is crossed this round. Rule:

   ```
   print if floor(minorIter / freq) ~= floor((minorIter - step) / freq)
   ```

   where `step` is this round's increment to `minorIter` (`1` for the `+1` paths,
   `batch` for the `useMexRandomize` path). This prints **at most once per window**
   even when a single batch advances `minorIter` by more than `freq`. Text is
   unchanged from today:

   ```
   resolve round %d: %d unconverged, worst residual %g
   ```

   ASG keeps its `asg ` prefix: `asg resolve round %d: ...`.

2. **End-of-loop summary.** Printed once after the loop exits, only when
   `verboseRetry && minorIter > 0` (no summary on the happy path where the main
   solve converged everything with zero resolve rounds):

   - All converged:
     ```
     resolve converged after %d rounds
     ```
   - Hit the cap (`minorIter >= maxMinorIter` with points still unconverged):
     ```
     resolve stopped at MaxMinorIter=%d: %d points still unconverged, worst residual %g
     ```

   ASG carries the same `asg ` prefix on the summary.

### Shared helper

Introduce `src/+gdsge/+runtime/printResolveProgress.m`, mirroring the existing
`printIterProgress.m`, to hold the window-crossing logic and the summary text so all
three call sites are one-liners and behave identically.

Proposed shape (final signature decided during implementation):

```matlab
function printResolveProgress(mode, label, minorIter, step, nUnconverged, worstF, ...
    freq, verboseRetry, hitCap)
% mode = 'round'   -> throttled per-round line
% mode = 'summary' -> end-of-loop summary (converged vs cap, chosen by hitCap)
% label = '' for cartesian, 'asg ' for ASG
```

Call sites:

- Per round: replace each `if cfg.verboseRetry; fprintf(...); end` block with a
  single `gdsge.runtime.printResolveProgress('round', label, minorIter, step, ...)`.
- After the loop in both files: one
  `gdsge.runtime.printResolveProgress('summary', label, minorIter, ...)` call,
  guarded by `verboseRetry && minorIter > 0`.
- Add the per-round call to the `useMexRandomize` branch, which is currently silent;
  pass `step = batch` there.

## Testing

TDD — write the failing unit test first. A focused test that drives
`printResolveProgress` directly (capturing output with `evalc`):

- **Silent between boundaries:** `mode='round'` with `minorIter`/`step` that do not
  cross a `freq` boundary produces no output.
- **One line per crossed window:** crossing a boundary prints exactly once; a single
  `step > freq` batch still prints exactly once (not once per window inside the
  batch).
- **Master switch:** `verboseRetry = false` produces no output in any mode.
- **Summary, converged:** `mode='summary'`, `hitCap=false` prints the
  `resolve converged after N rounds` line; `minorIter == 0` prints nothing.
- **Summary, cap hit:** `mode='summary'`, `hitCap=true` prints the
  `resolve stopped at MaxMinorIter=...` line with the unconverged count and worst
  residual.
- **ASG label:** `label='asg '` produces the `asg ` prefix.

Run via `matlab -batch "cd('tests'); run_tests"` (PowerShell 7 is unavailable, per
project notes); authoritative result is `tests/results/junit.xml` + exit code.

## Backward compatibility

- New option defaults to 100; existing `.gmod` models gain it automatically through
  `defaultSetupCode`. No public API or result-struct change.
- `verboseRetry`/`NoPrint` semantics are preserved (master mute unchanged).
- Visible change for users who left printing on:
  - ASG and cartesian-legacy paths: far fewer resolve lines (throttled to every
    `ResolvePrintFreq` rounds) plus one summary line at the end.
  - The cartesian `useMexRandomize` path (previously silent per-round): gains
    throttled per-round lines and the summary — a small, bounded increase.
