# Performance report (new vs old toolbox)

Point-in-time measurement on `PCWIN64`, 8 compute threads. Iter time is the median of 3 reps; SaveFreq=inf (no disk dumps). **Re-measured 2026-06-16** on the current toolbox (after the in-MEX `SIMU_INTERP`, `myppual` retirement, and the ASG-kernel optimization). Regenerate the headline table with `tests/perf/run_perf.ps1` or `matlab -batch "cd('tests/perf'); run_perf"`; the breakdown tables come from `scratch/profile_breakdown.m`.

| Model | Iter count (old / AD / SymPy) | old iter s | AD iter s | SymPy iter s | AD/old | codegen+compile s (AD) |
|---|---|---|---|---|---|---|
| HL1996 | 209 / 209 / 209 | 1.22 | 0.79 | 0.74 | 0.65 | 6.3 |
| GLSW_interp | 903 / 903 / 903 | 0.53 | 0.35 | 0.41 | 0.65 | 7.4 |
| CaoKS2016 | 281 / 281 / 281 | 8.88 | 4.60 | 4.58 | 0.52 | 6.1 |

> **`safe_assets` is omitted this run** — it did not produce a clean measurement. Its multi-root grid point makes the randomized restart fragile, and the perf worker does not pin `rng` before iter, so the run either fails to converge or errors at `MaxMinorIter`. This is the known `safe-assets-rng-restart-fragility` issue and the long-standing follow-up #1 below; the e2e gates pin `rng` and pass. Pin `rng` in the worker to re-measure it.

**Reading:** the autodiff backend emits the same adept/CoDoSol/kernel structure as the old toolbox, so AD/old ≈ 1 is the no-regression bar. The SymPy column quantifies the analytic-Jacobian backend vs autodiff. Iter counts must match across columns (same convergence).

## Findings

**The new toolbox is now faster than the old on every cleanly-measured model** (AD/old 0.65 / 0.65 / 0.52) — a step up from the Phase-10 measurement, where the spline models sat near parity (≈1.0) and only the ASG model led. The gains since then come from the in-MEX randomize/resolve path, the fused `interp_construct_mex` spline constructor (myppual retired), and the ASG-kernel optimization. The **SymPy** analytic-Jacobian backend is within noise of autodiff on HL1996 (0.74 vs 0.79s) and CaoKS2016 (4.58 vs 4.60s), and slightly behind on the tiny GLSW body (0.41 vs 0.35s) where recording an adept tape is already cheap. All iter counts and metrics match across columns.

**Where the spline-path iter time goes** (MATLAB profiler self-time, iter call only; compilation excluded):

| profiled | MEX solve | `iter_<model>` loop | spline rebuild\* | `solveProblems` | metric+print |
|---|---|---|---|---|---|
| HL1996 (AD) | 49% | 41% | 5% | 4% | <1% |
| HL1996 (SymPy) | 46% | 43% | 5% | 4% | <1% |
| GLSW (AD) | 48% | 14% | 14% | 19% | ~2% |
| GLSW (SymPy) | 40% | 17% | 17% | 21% | ~2% |
| Mendoza (SymPy) | 74% | 15% | 8% | 3% | <1% |

\* `constructSplines` + `interp_construct_mex` + `interp_eval_mex`

- The **generated `iter_<model>.m` loop** (per-iteration data marshalling for the MEX) is now a
  co-dominant MATLAB overhead on the cheap-solve models (HL1996 ~41%), and is **backend-agnostic** —
  autodiff and SymPy self-times match within noise. As the MEX solve got faster, this fixed overhead
  became a larger *fraction*.
- **Per-iteration spline reconstruction dropped sharply** versus Phase 10 (GLSW ~30% → ~14%): the fused
  `interp_construct_mex` constructor — which builds every cartesian interpolant in one call and writes
  the eval-order vector directly — replaced the per-variable `myppual` + `convert_to_interp_eval_array`
  rebuild that used to dominate.
- **`solveProblems` rises on models with a resolve cascade** (GLSW ~19%); on Mendoza the cascade now runs
  largely inside the MEX (in-MEX resolve), so its MATLAB-side `solveProblems` self-time is small (~3%)
  and the **MEX solve dominates (~74%)**. The number of real MEX solves per iteration is still the driver
  of the "slow" models.
- **`computeMetric` / `printIterProgress` are always negligible (<2%).**

**ASG (CaoKS2016) — diagnosis and optimization record (2026-06-15; top-level split re-confirmed 2026-06-16).**
The 2026-06-16 re-profile gives MEX solve **43% (AD) / 45% (SymPy)**, `iter_<model>` loop **~13%**, and the
ASG interpolation machinery **~40%** — i.e. the MEX solve is now the single largest bucket, exactly the
post-optimization state predicted below. The subsection that follows is the dated **before → after** record
of how that came to be; its first table is the *pre-optimization* diagnosis that motivated the work.

The ASG iter splits very differently from the spline path (profiler self-time, iter call only;
compilation is excluded). **Pre-optimization diagnosis:**

| profiled | MEX solve | `iter_<model>` loop | ASG construct/eval† | `solveProblemsAsg` | metric+print |
|---|---|---|---|---|---|
| CaoKS2016 (AD)    | 22% | 8% | 68% | 2% | <1% |
| CaoKS2016 (SymPy) | 24% | 7% | 67% | 2% | <1% |

\† the `asg` class (`construct` / `push_eval_results_at_grids` / `get_eval_grids` / `eval_vec`) plus
the `asg_mex` kernel — the adaptive-grid build/refine/evaluate, all MATLAB-side and **backend-agnostic**.

- At diagnosis the dominant ~67% was the **ASG interpolation machinery**, identical across backends; the
  **MEX solve** was only ~22–24%, and the generated **iter-file marshalling just ~7–8%** — much smaller
  than on the spline models, where the fixed grid makes the marshal loop + spline rebuild the bulk.
- This is why **SymPy ≈ autodiff on ASG**: the analytic-Jacobian backend only touches the MEX-solve slice,
  and there its solve is **no faster than adept** (the per-solve cost is dominated by the in-MEX ASG interp
  eval + chain rule, not the Jacobian assembly). Unlike the spline models (where SymPy's solve is fastest),
  SymPy gets no net edge on ASG.
- **The optimization lever for ASG was the `asg`/`asg_mex` construct-eval path** — not the Jacobian backend
  and not the iter-file marshalling.

**Inside the ASG bucket (asg_mex C++ time by command, backend-agnostic; pre-optimization).** Every `asg`
method funnels through `asg_mex(<command>)`, so attributing `asg_mex` time by calling method gives the
per-operation cost:

| asg_mex command | % of ASG bucket | what it does |
|---|---|---|
| `push_eval_results_at_grids` | ~38% | store the solved policies into the **sol-interp** (the next-iter warm-start source) |
| `push_eval_results_at_valid` | ~34% | store the `var_interp` values **and compute the hierarchical surplus** that drives refinement |
| `delete` | ~11% | destruct the two per-iteration `asg` C++ objects |
| `get_eval_grids` | ~9% | **propose** new grid points to refine |
| `eval_vec` | ~4% | evaluate the interpolant |
| `get_grids_info` | ~4% | feed `compute_inf_metric` |
| `get_current_level` | ~1% | inner-loop condition checks (~19k calls) |

`eval_vec` (~4%) splits by caller as `compute_inf_metric` ~56% / iter-loop warm-start ~43% /
`solveProblemsAsg` cascade ~1%.

- The dominant **~72% was the two `push` operations** — inserting solved values into the sparse grid and
  computing the surplus. That is the ASG construction itself (each iter runs ~11 refinement levels →
  ~3100 pushes), inherent to the method.
- **Warm-start eval is *not* a cost** (~2% of the ASG bucket); neither is grid proposal (~9%). The
  surprise was **object `delete` at ~11%** — per-iteration churn of fresh `asg` objects.
- Levers, in order: the **push/surplus kernel** in `asg_mex`, then the **per-iteration object churn** —
  not warm-start, eval, grid proposal, or the Jacobian backend.

**ASG kernel optimizations (applied 2026-06-15).** Acting on the diagnosis above, four changes to the
`asg`/`asg_mex` kernel (all **bit-identical** — CaoKS2016 + Bianchi2011_asg goldens on both backends + the
sympy↔adept↔FD Jacobian cross-check unchanged):

1. **Enable OpenMP** (`-DUSE_OMP`) — the parallel-for loops in `asg.h` (the eval inside `push_*` + the
   surplus loops) were guarded by `#ifdef USE_OMP`, which the compile never defined (`/openmp` alone is
   inert). *But instrumentation showed the eval is only ~0.2 s* — not the bottleneck.
2. **Incremental level-combination set** — `update_level_combinations` rescanned the whole grid map on
   every push (O(N²) over a run); now maintained incrementally. Removes an algorithmic wart; small on
   the corpus.
3. **Chunked surplus pool (the big one)** — the real cost was the **serial per-grid insertion**: `info`
   stored `GridInfo` *by value* with `surplus` sized to `ASG_MAX_NVEC=100` (an **840-byte** value), so
   every insert copied 840 B and every rehash re-copied all N. Moving `surplus` to a `double*` into a
   bulk-allocated chunked pool shrinks `GridInfo` to **48 B** (no per-point heap alloc, no `ASG_MAX_NVEC`
   ceiling change; readers unchanged since `surplus[i_vec]` works on a pointer).
4. **Ping-pong buffer reuse (the other big one)** — the loop allocated a fresh NEW interp every iteration
   and deleted the old OLD. Keep **two persistent buffers per interp** (value + sol), `reset()` them
   (cleared, **capacity retained**) and swap each iteration. The decisive saving is **avoiding hash-map
   rehashing**: re-inserting into a capacity-retained map skips the repeated rehash as a fresh map grows
   0→N, which was the bulk of the push cost (push_at_grids 1.27→0.35, push_at_valid 1.16→0.24, delete
   0.42→~0). The warmup/init OLD is used once then dropped, never recycled — recycling it corrupted
   two-stage WarmUp solves, whose warmup interp carries the prior stage's (different) `stateRange` that
   `reset()` does not touch.

**Cumulative:** the `asg_mex` bucket **6.62 s → 1.40 s** (push_at_grids 2.51→0.35, push_at_valid
2.28→0.24, delete 0.71→~0); CaoKS2016 iter **9.26 s → ~4.6 s (AD/SymPy)** — the new toolbox is now
**~1.9× faster than the old** on ASG (AD/old 0.52). The MEX solve is now the largest bucket (~43%, re-confirmed
2026-06-16) — the ASG bookkeeping is no longer the bottleneck. SymPy stays ≈ AD (the optimized machinery is
backend-agnostic).

**Follow-ups** (do not block the perf snapshot, which measures rather than optimizes): (1) **pin `rng`
before iter for safe_assets in the perf worker** — without it the model fails to produce a clean
measurement (confirmed this run); (2) if optimizing the spline path, the levers are the **generated
iter-loop marshalling** and `solveProblems` on cascade models — not the Jacobian backend (SymPy's is
already the fast part) and not the convergence metric.

---

## SIMU_INTERP whole-loop MEX (2026-06-16)

Cartesian `SIMU_INTERP=1` simulation now runs the entire period loop inside the
generated `simulate_<model>_mex` (search-once per site, OpenMP over samples,
in-place writes). A/B vs the per-period `interp_eval_mex` fallback, both at
`num_periods=100000`, `num_samples=100`, 8 threads (informational):

| Model | whole-loop MEX | per-period fallback | speedup |
|---|---|---|---|
| GLSW_interp (1 state, unprimed) | 2.59 s | 6.43 s | 2.5× |
| mendoza2010 (2 states, primed)  | 6.77 s | 14.89 s | 2.2× |

The win is avoiding `num_periods` MATLAB↔MEX round-trips and the per-period
site/index array churn; results are numerically identical (the interpolant is
evaluated at the same integer-shock sites). Reproduce: `matlab -batch
"addpath('src','src/kernels'); cd('scratch'); bench_simu_interp"`.
