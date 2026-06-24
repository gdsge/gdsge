# Fused In-MEX Spline Construction (Direct Eval-Order Output) — Design Spec

- **Date:** 2026-06-14
- **Status:** Approved (design); pending plan
- **Owner:** Wenlan Luo
- **Scope:** Post-Phase-10 performance enhancement (see `PROGRESS.md`; parent design
  `docs/superpowers/specs/2026-06-11-refactor-gdsge-design.md`). Not a numbered phase.
  Sibling to the in-MEX nearby-resolve enhancement
  (`docs/superpowers/specs/2026-06-14-in-mex-nearby-resolve-design.md`).

---

## 1. Context

Every GDSGE outer iteration rebuilds the interpolants of the policy/value functions from the
just-solved grid values, then hands them to the solver MEX for the next sweep of grid-point
solves. On the **cartesian-grid spline/linear path** this rebuild runs entirely in MATLAB
(`src/+gdsge/+runtime/constructSplines.m`) and costs, per interp variable:

1. `myppual(pp)` — dispatch through the 1176-line `myppual.m` gateway into `myppual_mex`
   (one MEX call **per interp variable**) to construct natural-order coefficients.
2. `convert_to_interp_eval_array(ppCell)` — **(a)** loop-copy each pp's `.coefs(:)` into one
   stacked array (a concat copy), then **(b)** `permute(coefs, newPermuteOrder)` into the
   evaluation order the solver MEX requires.

The construction math is already **ours**: `myppual_mex.cpp` (a ~70 KB MEX — *not* Intel MKL)
is a thin `dfdConstruct1D` shim (`include/mkl_dummy_interp.h`) over
`include/cubic_spline_notaknot.h` (tridiagonal not-a-knot cubic) and
`include/linear_interp_construct.h`. The MKL names are vestigial.

The Phase-10 perf report attributes up to **~30 %** of cheap/high-iteration model runtime
(GLSW) to this per-iteration rebuild. Three overheads are removable:

- the `myppual.m` gateway dispatch before each construct,
- **N** separate construct MEX calls (one per interp variable),
- the `convert_to_interp_eval_array` concat-copy **+** `permute`.

The `permute` is the load-bearing one: it exists only because construction emits coefficients
in **natural (curve-fit) order** `[i_array, order, pieces]` while the solver's vectorized
evaluator (`include/InterpEval.h` via `include/MatlabInterpEval.h`) consumes them in
**eval order** `[i_vec, order…, pieces…, i_array]`. A permute is a pure reordering of the same
doubles — so a constructor that writes directly into eval order produces **byte-identical**
output, with no permute and no copy.

## 2. Goals

1. Collapse the per-iteration cartesian spline/linear rebuild into **one** dedicated construct
   MEX call (`interp_construct_mex`) that:
   - takes **all** interp `Values` arrays via a **struct** (one field per variable) — MATLAB
     shares the buffers, so there is **no concat/copy**;
   - runs the **identical** `construct_*` math;
   - writes the `GDSGE_SPLINE_VEC.coefs` buffer **directly in eval order** (no `permute`);
   - in the **same call** also returns the per-variable natural-order pp structs.
2. **Bit-exact parity** with the current path on every existing golden. Reuse captured
   goldens; capture none. Verification = the existing end-to-end gates plus a new in-process
   A/B differential test and a synthetic driver test.
3. Leave the **C++ solver MEX untouched**: identical `GDSGE_SPLINE_VEC` field set and identical
   `GDSGE_PP_<name>` per-variable pp contract. Only the MATLAB-side constructor changes.
4. **Ship a minimal, self-contained construction kernel.** The new `interp_construct_mex.cpp`
   carries **construction only** — `mkl_start` / `mkl_start_with_extrap` + helpers + the
   eval-order packing + struct I/O. None of `myppual_mex.cpp`'s legacy machinery (the
   `eval2`/`eval4`/`evalN` variants, partial construct-eval, search-index in/out, reductions,
   `griddedInterpolant`, rational/`prp` splines) is carried over. Committed with its prebuilt
   MEX and a cache-gated auto-compile (the `asg_mex` precedent). The full legacy
   `myppual_mex.cpp` source is **not** vendored.

## 3. Non-goals

- **ASG interpolation.** ASG builds an irregular adaptive sparse grid and has its own
  construct/eval path (`asg`/`asg_mex`); it does not use `myppual`/`GDSGE_SPLINE_VEC`.
  `solveProblemsAsg.m` and the ASG codegen stay exactly as today.
- **pchip.** The `mkl_dummy_interp.h` constructor supports cubic not-a-knot (order 4) and
  linear (order 2) only. pchip is out of scope and keeps its current path.
- **The simulate path.** `emitSimulate`/`emitSimulateInterp` and `emitResultIter`'s
  `IterRslt.output_interp` keep using `myppual`/`myppual_mex` on natural-order pp structs.
  Those are not the iteration hot loop. The existing **prebuilt `myppual_mex.<mexext>`** stays
  committed and unchanged to serve this evaluation path; only its *source* is not vendored
  (it is legacy for the new toolbox, which constructs via `interp_construct_mex`).
- **Folding construction into the solver MEX.** Construction stays a dedicated MEX call between
  solves (the owner's "ask `myppual_mex` to directly output" framing). An in-solver-MEX
  construction is a possible future, not this spec.
- **Redesigning `GDSGE_SPLINE_VEC`.** The struct field set is frozen so the solver MEX is
  untouched and the A/B gate stays bit-exact. A leaner eval struct is a possible future.

## 4. Architecture & data flow

```
                       per outer iteration (cartesian spline/linear path)
  BEFORE:
    for each interp var i:  myppual(pp_i) ──► myppual_mex (construct, natural order)   [N MEX calls]
    convert_to_interp_eval_array(ppCell):  stack-copy + permute ──► GDSGE_SPLINE_VEC   [copy + permute]
    ppCell (natural-order pp per var) ─────────────────────────────► GDSGE_PP_<name>, IterRslt.pp

  AFTER:
    interp_construct_mex(breaks, GDSGE_INTERP_VALUES struct, order, extrap, threads)    [1 MEX call]
        ├─ construct_* (identical math)  ──► per-var natural-order coefs (intermediate)
        ├─ write directly in eval order  ──► GDSGE_SPLINE_VEC.coefs   (no permute)
        └─ wrap intermediates            ──► GDSGE_PP_CELL (natural-order pp per var)
    returns [GDSGE_PP_CELL, GDSGE_SPLINE_VEC]
```

The solver MEX (`templates/cxx/interp_spline_*.tpl.cpp` → `MatlabInterpEval`/`MatlabPp`) reads
the **same two contracts** it reads today and is **not modified**.

## 5. New component — `interp_construct_mex`

Flat kernel name matching `myppual_mex` / `asg_mex`. No old-toolbox name collides (the MATLAB
golden-path rule is satisfied).

**Signature (MATLAB call):**

```matlab
[ppCell, splineVec] = interp_construct_mex(breaks, values, order, extrapOrder, numThreads)
```

- `breaks` — `1×D` cell of state grid vectors (the shared grid, as in `constructSplines` today).
- `values` — a **struct** with one field per interp variable; each field is that variable's
  value array `reshape(valueCell{i}, [], sizeState{:})` (first dim = the array/shock dimension
  `numArray`, remaining dims = the state grid). Passed as a struct so MATLAB does not copy the
  underlying arrays (`mxGetField`/`mxGetPr` return the existing data pointer). Field **order**
  defines interp-variable order `i_vec = 1..numVec`.
- `order`, `extrapOrder`, `numThreads` — scalars/vectors exactly as `constructSplines` computes
  them today (`orderVec`, `extrapVec`, `NumThreads`).

**Outputs:**

- `ppCell` — `1×numVec` cell of natural-order pp structs, **identical in shape** to what
  `myppual(pp)` returns today (`form/breaks/Values/coefs/pieces/order/extrap_order/dim/Method/
  orient/thread` as consumed by `MatlabPp` and stored in `IterRslt.pp`). Feeds the solver's
  named-scalar interp calls (safe_assets/GLSW/Cao2011EZ) and `IterRslt.pp`.
- `splineVec` — the `GDSGE_SPLINE_VEC` struct, **identical field set** to
  `convert_to_interp_eval_array`'s output (`coefs, fullVecEvalCoefsLength,
  singleVecEvalCoefsLength, xDim, order, pieces, xPts, breaks, dim, arrayOffset`), with `coefs`
  laid out directly in eval order.

**Construction reuse (bit-exactness mechanism).** The exact multi-dimensional construction
driver that `myppual_mex.cpp` uses (`mkl_start` / `mkl_start_with_extrap` + the per-dimension
tensor sweep + the extrapolation break-extension, plus the `intmemcpy` / `doublememcpy` /
`vectorProd` helpers) is copied **verbatim** into `interp_construct_mex.cpp`, which `#include`s
`mkl_dummy_interp.h` (pulling `cubic_spline_notaknot.h` / `linear_interp_construct.h` and the
`mkl_domatcopy` transpose). Identical source → identical floating-point operations → bit-identical
coefficients. The driver test (§10) asserts these coefs equal the shipped prebuilt `myppual_mex`'s
coefs, catching any drift in the verbatim copy.

## 6. Eval-order layout contract

`GDSGE_SPLINE_VEC.coefs` must be **byte-identical** to today's
`convert_to_interp_eval_array(ppCell).coefs`. That function is the normative reference; the C++
writer reproduces its index map:

- natural per-var coefs have logical shape `[numArray, order_1, pieces_1, …, order_D, pieces_D]`
  (curve-fit order; `dim = numArray`);
- eval order is `[i_vec, order_D..order_1, …, pieces_1..pieces_D, i_array]`, i.e. interp
  variable `i_vec` is the fastest dimension, the per-piece order blocks come next (reversed per
  the MKL highest-to-lowest convention `InterpEval.h` decrements through), then the piece grid,
  with the array index `i_array` outermost at stride `arrayOffset`.

The implementation plan derives the exact strides; the **contract** is "equal to
`convert_to_interp_eval_array`," enforced bit-for-bit by the driver test. `fullVecEvalCoefsLength`,
`singleVecEvalCoefsLength`, and `arrayOffset` are computed identically (they are pure functions
of `order`, `pieces`, `numVec`, `numArray`).

## 7. MATLAB integration

- **`gdsge.runtime.constructSplines`** is rewritten to build the `values` struct from the input
  `valueCell` and make the single `interp_construct_mex` call, returning `[ppCell, splineVec]`.
  **Signature unchanged** — so `emitIter.m`'s three call-sites (main loop, WarmUp branch, init
  branch) and the generated `iter_<model>.m` are structurally unchanged.
- **`gdsge.codegen.ensureSplineConstructMex`** — a clone of `ensureAsgMex`: hash-cached
  (`interp_construct_mex.cache` beside the source via `gdsge.codegen.needsCompile`), compiles
  with `myppual`'s flags (§8), no-ops when current. Invoked wherever `ensureAsgMex` is invoked
  at codegen entry, so the kernel exists before the first solve. (Spline models currently need
  no kernel compile; this adds one, matching the ASG improvement over ship-prebuilt-only.)
- **`convert_to_interp_eval_array.m`, `myppual.m`, `myppual_mex`** remain vendored and on the
  path — still used by the simulate path and `IterRslt.output_interp` — just off the iteration
  hot path. No deletion (backward compatibility; frozen public surface).

## 8. Packaging & build

Committed into the package:

- `src/kernels/interp_construct_mex.cpp` — the new minimal, **self-contained** construction
  kernel (construction driver copied verbatim from `myppual_mex.cpp` + struct I/O + eval-order
  packing; no legacy eval/search/partial/reduction/`griddedInterpolant`/`prp` code).
- Prebuilt `src/kernels/interp_construct_mex.<mexext>` + `interp_construct_mex.cache`
  (the `myppual`/`asg` committed-prebuilt precedent).

Retained unchanged (already vendored): the prebuilt `src/kernels/myppual_mex.<mexext>` for the
simulate / `output_interp` eval path. The full legacy `myppual_mex.cpp` source is **not**
vendored. The three construction headers (`mkl_dummy_interp.h`, `cubic_spline_notaknot.h`,
`linear_interp_construct.h`) are already vendored under `include/` and verified
**byte-identical** to the owner's source tree — no reconciliation needed.

**Compile recipe** (cache-gated, from `compile_myppual.m`):
- Windows: `mex interp_construct_mex.cpp -DUSE_OMP -I<repo>/include COMPFLAGS="$COMPFLAGS /openmp"`
- Mac: `mex interp_construct_mex.cpp -I<repo>/include COMPFLAGS="$COMPFLAGS -std=c++14"`

**Bit-exactness & compile flags.** The A/B gate compares the new constructor against the shipped
prebuilt `myppual_mex` + `convert_to_interp_eval_array`. To keep coefs bit-identical, the
construct MEX compiles with MSVC default **`/fp:precise`** (never `/fp:fast`) and `/O2` — the
same effective floating-point model as the shipped `myppual_mex` (whose `mex` default is `/O2`,
`/fp:precise`). `/openmp` + `-DUSE_OMP` match `compile_myppual.m`; OpenMP only parallelizes
independent vector functions, not FP reduction order, so it does not perturb results.

## 9. Error handling & edge cases

- **Extrapolation** (order 4 with non-empty `extrapVec`): replicate `myppual_mex`'s
  break-extension (`pieces + 2`, `mkl_start_with_extrap`) byte-for-byte via the shared header.
- **Linear** (order 2): no extrapolation (matches today's `extrapVec = []` for `order ≠ 4`).
- **Empty interp list** (e.g. the init view with 0 interp variables): return an empty
  `splineVec` struct and empty `ppCell` — `constructSplines` / `emitInterp` already special-case
  `numel == 0`; the MEX must too.
- **Grid size:** the `myppual` 2-point-grid kernel bug is unaffected — construction math is
  unchanged, and tests use ≥3 points per dimension (existing convention).
- **Validation:** the MEX checks that struct-field count and each field's array shape are
  consistent with `breaks`/`order`; informative `gdsge:kernels:interpConstruct*` errors on
  mismatch (no silent misindex).

## 10. Testing (bit-exact, mirrors `tInMexResolve`)

1. **Driver unit test** (`tInterpConstructMex`, fast): for synthetic 1-D / 2-D / 3-D grids
   (≥3 pts/dim), cubic and linear, with and without extrapolation, multiple interp variables
   and `numArray > 1`, assert the fused MEX output is **bit-identical** to
   `myppual` + `convert_to_interp_eval_array`:
   - `splineVec.coefs` and every scalar field (`arrayOffset`, lengths, …) equal;
   - each `ppCell{i}.coefs` equal to `myppual(pp_i).coefs`.
   This also pins the §5 "shared header == shipped `myppual_mex`" invariant.
2. **A/B end-to-end differential gate** (`tFusedConstruct{HL1996,SafeAssets,...}`): run `iter`
   with the fused constructor vs the legacy `constructSplines` path and assert
   **bit-identical** results — `maxabsdiff = 0` on every `var_policy`/`var_aux`/`var_interp`
   field and `Metric` — with `rng` pinned before each run (the documented RNG-fragility of
   safe_assets' multi-root restart; same pin the end-to-end gates use). The legacy path stays
   reachable behind an internal switch in `constructSplines` (mirrors `UseMexResolve`) to drive
   this test and as a fallback.
3. **Existing corpus end-to-end gates** stay green unchanged (the constructor is internal):
   HL1996, safe_assets, Mendoza2010, GLSW, Cao2011EZ, CaoNie2016, both backends.
4. **Perf micro-bench** (reuse `tests/perf/`): construct-time and total iter-time before/after
   on GLSW and safe_assets (rng-pinned, identical Iter/Metric), logged to the perf report.

## 11. Risks & mitigations

- **Driver-copy drift** (the verbatim driver copy inside `interp_construct_mex.cpp` vs the
  shipped prebuilt `myppual_mex`): the §10.1 driver test asserts bit-identical coefs, failing CI
  on any divergence. The legacy kernel is frozen, so drift is unlikely.
- **Compile-flag FP divergence**: pinned to `/fp:precise /O2` (§8); the A/B gate would catch any
  ULP drift immediately.
- **Struct copy-on-pass**: confirm `mxGetField`/`mxGetPr` expose the existing buffer (no copy);
  the values struct is read-only in the MEX. If MATLAB ever copies-on-write, correctness is
  unaffected (only the copy-avoidance benefit is lost).
- **Eval-order packing bug on multi-D**: fully covered by the §10.1 2-D/3-D bit-exact cases
  against the normative `convert_to_interp_eval_array`.

## 12. File manifest

New:
- `src/kernels/interp_construct_mex.cpp` (minimal, self-contained), prebuilt `.<mexext>` + `.cache`
- `src/+gdsge/+codegen/ensureSplineConstructMex.m`
- `tests/kernels/tInterpConstructMex.m`
- `tests/<model>/.../tFusedConstruct*.m` (A/B gates; ≥ HL1996 + safe_assets)

Modified:
- `src/+gdsge/+runtime/constructSplines.m` (single fused call; internal legacy switch)
- `src/+gdsge/+codegen/generateCxx.m` (call `ensureSplineConstructMex` alongside `ensureAsgMex`)
- `docs/architecture.md` (kernel inventory + compile note), `PROGRESS.md` changelog

Unchanged (explicitly): the prebuilt `src/kernels/myppual_mex.<mexext>` (kept for the simulate /
`output_interp` eval path — source not vendored); the C++ solver MEX templates
(`interp_spline_*.tpl.cpp`, `MatlabInterpEval.h`, `MatlabPp.h`, `InterpEval.h`); `emitIter` call
shape; the simulate path; all ASG/pchip code; `IterRslt`/`SimuRslt` shapes.
