# In-MEX Nearby-Warmup Resolve — Design Spec

- **Date:** 2026-06-14
- **Status:** Approved (design); pending plan
- **Owner:** Wenlan Luo
- **Scope:** Post-Phase-10 performance enhancement (see `PROGRESS.md`; parent design
  `docs/superpowers/specs/2026-06-11-refactor-gdsge-design.md`). Not a numbered phase.

---

## 1. Context

Each GDSGE outer iteration solves every `(state × shock)` grid point as an independent
bounded nonlinear system. When a point fails to converge from its supplied initial guess,
the toolbox retries via a **resolve cascade**:

1. **Initial solve** — all grid points from the supplied guess (one MEX call).
2. **Nearby-neighbor warm-start sweep** — for each grid dimension, copy a *converged*
   neighbor's solution into the still-unconverged points and re-solve, looping until a full
   sweep produces no newly-converged point.
3. **Randomized restarts** — draw `rand·(ub−lb)+lb` for whatever remains and re-solve,
   bounded by `MaxMinorIter`.

In the current refactor this **entire cascade runs in MATLAB**
(`src/+gdsge/+runtime/solveProblems.m`). The C++ MEX (`templates/cxx/task.tpl.cpp`) solves
each point **exactly once per call**, gated by a per-point `GDSGE_SKIP` mask; the MATLAB
driver re-calls it for every sweep direction. A single outer iteration therefore fires roughly
`1 + 2·num_dims·K_progress + K_restart` MEX round-trips, each paying fixed marshalling cost
(`mexGetVariable` caller-workspace reads, output `mxArray` allocation, interp re-fetch). The
Phase-10 perf report attributed `safe_assets`' ~3× iter time to exactly this: ~12.5 MEX
solve-calls per iteration from the resolve cascade, versus ~1 on models that converge cleanly.

The **hans/dpopt** toolbox (the analytic-Jacobian reference, also authored by the owner) solves
this differently. With `solver_use_resolve=true`, its MATLAB driver
(`auto_dpopt/source/template/solve_vfi_template.m:127-133`) passes per-dimension
`ProblemStrides` into the solver and sets `max_minor_iter=0`, disabling the MATLAB sweep: the
**neighbor-warmup sweep runs inside the compiled solver** (`dpopt_lib/libdpopt.cpp:1545-1624`,
`:1684-1766`), and MATLAB re-engages only for the random trial. The in-MEX sweep is the same
progress-loop / per-dimension lower-then-upper structure that GDSGE's `solveProblems.m` Phase 2
is a parity port of.

This spec ports that idea: **move the neighbor-warmup sweep (cascade step 2) into GDSGE's own
C++ MEX**, collapsing all of step 2's MATLAB↔MEX round-trips into a single call, while leaving
the randomized restart (step 3) in MATLAB.

## 2. Goals

1. Run the nearby-neighbor warm-start sweep **inside `task.tpl.cpp`**, so an outer iteration's
   neighbor-propagation phase costs **one** MEX call instead of `~2·num_dims·K_progress`.
2. Cover the **regular cartesian-grid path** (spline / pchip), which serves **both** the adept
   autodiff and SymPy backends through the shared `task.tpl.cpp` — one implementation, both
   backends.
3. **Bit-exact parity** with every existing golden, on both backends. Reuse the captured
   goldens; capture none. Verification is the existing end-to-end gates plus a new in-process
   A/B differential test.
4. Keep the MATLAB-side sweep (`solveProblems.m` Phase 2) **reachable behind an internal
   `cfg` switch** (default = in-MEX), both as a safety fallback and to drive the A/B test.
5. Keep the **randomized restart in MATLAB** (parent §step 3), drawn from MATLAB's RNG, so the
   `safe_assets` multi-root RNG behavior is preserved exactly.

## 3. Non-goals

- **ASG interpolation.** ASG proposes an irregular adaptive sparse grid level-by-level; its
  warm-start is *interpolation*-based (`solInterpNew`/`solInterpOld`, evaluated in MATLAB), not
  cartesian stride-neighbors. The sweep concept does not map to it. `solveProblemsAsg.m` stays
  exactly as today. (A separate in-MEX scheme for ASG is out of scope.)
- **Randomized restart in C++.** Stays in MATLAB by design (Goal 5).
- **A documented user feature.** The switch is internal: an **undocumented runtime option**
  `UseMexResolve` (default `1`) added to the existing options whitelist, flowing into
  `cfg.useMexResolve`. It is a solver knob like `UseAdaptiveBoundInSol`/`MaxMinorIter`, not a
  gmod-language feature; it adds no IR field, changes no result shape, and is omitted from the
  user guide. Its only consumers are the in-MEX default and the A/B differential test.
  (Mirroring hans's `solver_use_resolve` as a *documented* public toggle was explicitly
  declined.)
- **Algorithmic change to results.** This is a *where-it-runs* refactor, not a *what-it-computes*
  one. Any speedup is reported informationally; parity is the bar.

## 4. The determinism decision (snapshot, not live-read)

libdpopt's in-MEX sweep reads the **live** `exitflag[]` array inside each parallel pass
(`libdpopt.cpp:1576` — `check_opt_solved(exitflag[i_warmup])` while sibling OpenMP threads write
`exitflag[]`). This is a benign race: a point that converges early in a pass *may or may not*
be visible as a warm-start source to a higher-indexed point in the **same** pass, depending on
thread scheduling — so libdpopt's resolve is **non-deterministic across runs**.

GDSGE's current MATLAB cascade instead **snapshots** the converged set, copies neighbors, then
solves — one propagation layer per direction-pass, fully deterministic. That snapshot semantics
produced every existing golden.

**Decision: port the algorithm but keep GDSGE's deterministic snapshot semantics**, not
libdpopt's live read. Per direction-pass:

1. Freeze the converged predicate into a local snapshot array `GDSGE_SOLVED0[i] =
   (GDSGE_F(i) <= TolSol && !isnan(GDSGE_F(i)))` for all `i`.
2. In one OpenMP `parallel for`, each point `i` that is unconverged in the snapshot and whose
   stride-neighbor is converged in the snapshot copies the neighbor's working solution and
   re-solves.

Eligibility is provably race-free: a warm-start **source** must be snapshot-converged, an
**eligible** point must be snapshot-unconverged, so the two sets are disjoint — no point both
supplies and receives a warm start in the same pass, and each eligible point writes only its own
solution column. The fused copy+solve is therefore both efficient (no separate serial copy pass)
**and** bit-identical to the MATLAB snapshot cascade. The divergence from libdpopt (snapshot vs
live read) is deliberate, documented, and the reason we get bit-exact verification.

## 5. Architecture

### 5.1 Gate: presence of strides

The resolve runs **iff** a non-empty per-dimension stride vector is supplied to the MEX. The
init task, ASG, and the MATLAB-fallback A/B mode simply do **not** supply it, so they hit the
exact current code path. No template fork, no `#ifdef`.

Strides follow the dpopt formula: for grid problem-shape `GDSGE_SIZE = [shock_num, n_state1,
n_state2, …]`, `stride(d) = prod(GDSGE_SIZE(1:d-1))` (so dim 1 = stride 1, i.e. the shock
dimension; dim 2 = stride `shock_num`; …). This matches `solveProblems.m`'s existing
`prod(cfg.probSize(1:iDim-1))` and `solve_vfi_template.m:128-130`.

### 5.2 `templates/cxx/task.tpl.cpp` (shared; the core change)

- **Read the optional strides** from the caller workspace via the existing `mexGetVariable`
  contract (a new contract variable, e.g. `GDSGE_PROBLEM_STRIDES`, a `double` row vector;
  `num_dims = numel`). Absent or empty → `num_dims = 0` → no resolve.
- **Refactor the per-point solve body into a reusable unit** — a `[&]`-capturing lambda
  `GDSGE_solve_one(int GDSGE_I)` wrapping the existing `START_LOOP_CODE … POP_CODE …
  PRE_MODEL_CODE … MODEL_CODE … CALL_FMIN_CODE … FINISH_LOOP_CODE` block. The adept `Stack` and
  `INTERP_GET_THREAD_CODE` stay **inside** the lambda so each invocation/thread gets its own
  (unchanged per-call semantics). **Contract:** the lambda assumes the warm-start guess is
  already in `GDSGE_SOL(:,GDSGE_I)` and solves it in place — the `sol0→sol` copy moves *out* of
  `call_fmin` to the call sites (§5.3).
- **`GDSGE_SOL` becomes the persistent working array**: it is both the warm-start source and the
  output. Pass 0 copies `GDSGE_SOL0(:,i)→GDSGE_SOL(:,i)` before solving each non-skipped point;
  the resolve passes copy a converged neighbor's `GDSGE_SOL(:,j)→GDSGE_SOL(:,i)` before solving.
- **Pass 0** (semantics unchanged): `#pragma omp parallel for` over all `GDSGE_I`, honoring
  `GDSGE_SKIP`; for a non-skipped point, `memcpy(sol0→sol)` then `GDSGE_solve_one(i)`. With
  `GDSGE_SKIP(:)=1` (the dedicated-sweep call, §5.4) Pass 0 no-ops and only the resolve driver
  runs.
- **Resolve driver** (only when `num_dims > 0`), mirroring `libdpopt.cpp:1561-1623` but with
  snapshot semantics:

  ```
  int before = count_unconverged();            // GDSGE_F(i) > TolSol || isnan
  int after  = BIG;
  while (after != before) {
      before = count_unconverged();
      for (int d = 0; d < num_dims; ++d) {
          int s = strides[d];
          snapshot SOLVED0[] from GDSGE_F;      // freeze, then lower pass
          #pragma omp parallel for schedule(dynamic)
          for i: if (!SOLVED0[i] && i-s>=0 && SOLVED0[i-s]) {
                     memcpy(&GDSGE_SOL(1,i), &GDSGE_SOL(1,i-s), NUM_EQ);
                     GDSGE_solve_one(i);         // updates GDSGE_F(i), aux, eqval
                 }
          snapshot SOLVED0[] from GDSGE_F;      // re-freeze, then upper pass
          #pragma omp parallel for schedule(dynamic)
          for i: if (!SOLVED0[i] && i+s<NPROB && SOLVED0[i+s]) {
                     memcpy(&GDSGE_SOL(1,i), &GDSGE_SOL(1,i+s), NUM_EQ);
                     GDSGE_solve_one(i);
                 }
      }
      after = count_unconverged();
  }
  ```

  Snapshot-before-each-direction matches `solveProblems.m`'s "recompute `needResolved` before
  the lower copy and again before the upper copy" (the documented parity quirk).

### 5.3 `templates/cxx/call_fmin.tpl.cpp` and `call_fmin_sympy.tpl.cpp`

Today both **begin** with `memcpy(GDSGE_SOL ← GDSGE_SOL0)` (in the `==0` solve branch and again
in the `==1||==2` eval-only branch). **Remove both memcpys.** The warm-start guess is now placed
into `GDSGE_SOL(:,i)` by the caller (Pass 0 copies from `GDSGE_SOL0`; a resolve pass copies from
the converged neighbor), so `CoDoSol::solve` / the eval simply operate in place on `GDSGE_SOL`.
This is what makes the per-point body a reusable unit whose warm start can vary by call.

The skip-path memcpys (`task.tpl.cpp:64-66`, copying prior `sol0/aux/eqval` for skipped points)
are preserved. The `GDSGE_DEBUG_EVAL_ONLY` 1/2 branches still eval in place on `GDSGE_SOL`
(which Pass 0 set from `GDSGE_SOL0`), so eval-only/analytic-Jacobian-debug behavior is unchanged;
those callers pass no strides, so the resolve driver never runs in debug-eval mode.

### 5.4 `src/+gdsge/+runtime/solveProblems.m`

**Design principle: preserve the current control flow exactly; replace only the *body* of the
MATLAB neighbor sweep with one strides-enabled MEX call.** Do **not** fold the sweep into Pass-0
or the restart calls — that would move the convergence check and break RNG-draw parity (§7).

Add a `GDSGE_PROBLEM_STRIDES = []` local in the MEX caller-workspace contract block (so Pass 0
and the restart calls see an empty strides vector → no internal resolve). Add `cfg.useMexResolve`
(default `true` via the generated cfg). Compute `strides = cumprod([1, cfg.probSize(1:end-1)])`
once when `useMexResolve` is on (this is `prod(probSize(1:d-1))` for `d=1..D`).

Replace the existing `if cfg.useNearestNeighbor … end` block (today's lines 51–85) with:

```matlab
if cfg.useMexResolve
    % In-MEX neighbor sweep. Replicate the OLD inner-while ENTRY gate exactly
    % (numNeedResolvedAfter persists across outer iterations, init = inf), then
    % let the MEX run the entire progress loop internally in ONE call.
    needResolved = (f > cfg.tolSol) | isnan(f);
    numNeedResolved = sum(needResolved);
    if numNeedResolvedAfter ~= numNeedResolved
        skip(:) = 1;                       % Pass 0 no-ops; strides drive the sweep
        GDSGE_PROBLEM_STRIDES = strides;    %#ok<NASGU>  % enable in-MEX resolve
        [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);
        GDSGE_PROBLEM_STRIDES = [];         %#ok<NASGU>  % off for the restart call below
        numNeedResolvedAfter = sum((f > cfg.tolSol) | isnan(f));
    end
elseif cfg.useNearestNeighbor
    ... today's MATLAB sweep, verbatim (the inner while + per-dim lower/upper) ...
end
```

Everything else is **unchanged**: `numNeedResolvedAfter = inf` before the outer while; the
outer `while ((max(isnan(f)) || max(f(:)) > tolSol) && minorIter < maxMinorIter)`; Phase 3
random restart with its single per-outer-iteration `rand(size(sol))` draw; the optional
`adaptInSol` tightening; `diag`. The dedicated sweep call carries `GDSGE_USE_BROYDEN_NOW = 0`
(set after Pass 0, as today). When `useMexResolve=false`, the `elseif` runs today's MATLAB sweep
and `GDSGE_PROBLEM_STRIDES` stays `[]` throughout — byte-for-byte the current behavior.

The persistent `if numNeedResolvedAfter ~= numNeedResolved` gate is what makes the MEX call fire
at *exactly* the outer iterations the MATLAB inner-while would have run (including the obscure
case where a stale count equals the fresh count and the sweep is skipped). Inside the MEX, the
`while (after != before)` progress loop reproduces the MATLAB inner-while's rounds (§7).

### 5.5 `src/+gdsge/+codegen/+mat/{emitSetup,optionsWhitelist,emitIter}.m`

- **`emitSetup.m`** (defaults block): add `UseMexResolve = 1;` next to `UseAdaptiveBoundInSol = 0;`.
- **`optionsWhitelist.m`** (`base` list): add `'UseMexResolve'` so the runtime override is
  accepted by `unpackOptions`.
- **`emitIter.m`** (cfg spine, ~line 173): add `GDSGE_CFG.useMexResolve = UseMexResolve;` and
  keep `GDSGE_CFG.useNearestNeighbor = true;` (now the fallback selector). `emitIterAsg.m` is
  untouched — ASG never sets `useMexResolve`, so its `solveProblemsAsg` cascade is unaffected
  (and `solveProblemsAsg` does not read `cfg.useMexResolve`).

## 6. Components touched (summary)

| File | Change |
|------|--------|
| `templates/cxx/task.tpl.cpp` | per-point body → lambda; `GDSGE_SOL` working array; read optional `GDSGE_PROBLEM_STRIDES`; snapshot resolve driver; Pass-0 `sol0→sol` copy moves here |
| `templates/cxx/call_fmin.tpl.cpp` | drop the `sol0→sol` memcpys (warm start now pre-placed in `GDSGE_SOL`) |
| `templates/cxx/call_fmin_sympy.tpl.cpp` | same |
| `src/+gdsge/+runtime/solveProblems.m` | `GDSGE_PROBLEM_STRIDES` contract local; `cfg.useMexResolve`; dedicated sweep call w/ persistent entry-gate; MATLAB sweep kept under `elseif useNearestNeighbor` |
| `src/+gdsge/+codegen/+mat/emitSetup.m` | default `UseMexResolve = 1;` |
| `src/+gdsge/+codegen/+mat/optionsWhitelist.m` | add `'UseMexResolve'` |
| `src/+gdsge/+codegen/+mat/emitIter.m` | `GDSGE_CFG.useMexResolve = UseMexResolve;` |
| C++ codegen snapshot goldens | regenerate (template text changed); EOL-normalized diff review |
| tests | new A/B differential gate; round-trip-count characterization |

`emitIterAsg.m`, `solveProblemsAsg.m`, the ASG templates, the IR schema, and all public API are
**unchanged**.

## 7. Why bit-exact parity holds

The design deliberately keeps `solveProblems.m`'s control flow byte-for-byte and swaps only the
*implementation* of the neighbor sweep from MATLAB-orchestrated calls to one in-MEX call. So:

- **Identical control flow ⇒ identical RNG stream.** The outer `while`, the single
  `rand(size(sol))` draw per outer iteration, the random-restart application, and the
  convergence checks are all in the same place as today. The dedicated sweep call does not draw
  randomness and does not move any check. Therefore the RNG stream is bit-identical — the
  `safe_assets` multi-root selection is preserved. *(This is precisely what the rejected
  "folded" variant got wrong: folding the sweep into the solve call moved the post-sweep
  convergence check ahead of the restart, so a sweep that fully converged a previously-incomplete
  state would skip a `rand` draw the old path makes, desyncing every later iteration.)*
- **Persistent entry-gate ⇒ the MEX fires at the same iterations.** The `if numNeedResolvedAfter
  ~= numNeedResolved` gate (with `numNeedResolvedAfter` persisting across outer iterations,
  `inf` initially) reproduces the MATLAB inner-while's entry condition exactly — including the
  obscure stale-count==fresh-count case where the sweep is skipped entirely.
- **Internal progress loop ⇒ same rounds.** Inside the MEX, `while (after != before)` over the
  per-dimension lower/upper passes reproduces the MATLAB inner-while's rounds in the same order,
  with the same `solved = (f <= TolSol)` predicate and the same `TolSol`.
- **Snapshot semantics ⇒ same per-round propagation.** Freezing `solved[]` before each direction
  pass (§4) reproduces the MATLAB sweep's "recompute `needResolved`, then copy from a converged
  neighbor, then solve" — one propagation layer per pass. The fused copy+solve is bit-identical
  to MATLAB's copy-all-then-solve-all because a warm-start source is always snapshot-converged
  (hence not modified in the pass) and each eligible point writes only its own column.
- **Same arithmetic.** The per-point solve is the identical `CoDoSol::solve` over identical
  inputs (same warm start, `lb/ub/data`, tolerances). Grid points are independent — no
  cross-point floating-point reduction — so OpenMP scheduling never affects a result.

No edge-case deviation remains; the A/B differential test (§8) is expected to be bit-identical
on every model, not merely within tolerance.

## 8. Verification

1. **Existing goldens, both backends.** All end-to-end gates — HL1996, safe_assets, Mendoza2010,
   GLSW, Cao2011EZ, CaoNie2016 — pass unchanged under `useMexResolve=true`, on adept **and**
   SymPy. (safe_assets, with its multi-root restart, is the key witness.)
2. **A/B differential test (new).** For at least HL1996 and safe_assets: run
   `iter_<model>(struct('UseMexResolve',1,…))` and `iter_<model>(struct('UseMexResolve',0,…))`
   on identical options and assert `IterRslt` `var_policy`/`var_aux`/`var_interp`/`Metric`/`Iter`
   are **bit-identical** (`isequal`, not tolerance). This directly proves the in-MEX resolve
   reproduces the MATLAB sweep. safe_assets is the decisive case (it actually exercises the
   sweep + random restart; clean models converge in Pass 0 and never enter the sweep).
   **Required:** pin `rng(SEED)` before *each* run. The iter loop draws `rand()` for the initial
   guess and every random restart from the global RNG; without an identical seed the two runs
   draw different sequences and diverge to ~1e-12 (single-root) or different multi-roots — a
   *test* artifact, not a resolve difference (verified: identically seeded, both models are
   bit-for-bit identical). This is the same RNG-fragility the end-to-end gates already pin for.
3. **Round-trip characterization (new, informational).** Count MEX invocations per outer
   iteration; assert the new path issues strictly fewer on `safe_assets` (expected ≈1–2 vs
   ≈12.5). Reported, not a hard gate.
4. **C++ snapshot goldens** regenerated and reviewed (template text changed); functional gates
   are the real check.

## 9. Risks & mitigations

- **Hot-path restructuring.** The lambda + working-array change is in the most performance- and
  correctness-sensitive C++. Mitigation: keep Pass 0 byte-for-byte equivalent to today; the
  snapshot bookkeeping is serial-but-O(NPROB)-cheap; the A/B test catches any semantic drift.
- **Caller-workspace contract growth.** `GDSGE_PROBLEM_STRIDES` is a new `mexGetVariable`
  contract variable; it must exist as a local in `solveProblems.m`'s body (the contract requires
  every read variable to be a caller local). Mitigation: follow the existing `TolSol`/
  `GDSGE_SPLINE_VEC` pattern exactly; absent/empty handled as "no resolve."
- **pchip / extrapolation-order variants** share `task.tpl.cpp` and the cartesian grid, so they
  inherit the resolve for free — but none is gated end-to-end today. Mitigation: the gate is
  spline-family-agnostic (it keys on strides, not interp order); the A/B test on the gated
  spline models covers the mechanism.
- **safe_assets fragility.** The known multi-root RNG sensitivity is *preserved*, not
  introduced, because restart RNG stays in MATLAB and the solve sequence is identical; the A/B
  test is the guard.

## 10. Open items / future

- An in-MEX scheme for **ASG** (interpolation-warm-start propagation inside the MEX) is a
  possible later enhancement; out of scope here.
- If profiling later shows the random-restart round-trips themselves dominate on some model,
  moving restart into the MEX (behind the same stride gate) could be revisited — but it would
  forfeit RNG parity, so it is deliberately excluded now.
