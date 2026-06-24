# GDSGE user guide

How to author and solve a model with the refactored toolbox. This is self-contained for the
package; it complements the published documentation at gdsge.com. The toolbox is **backward
compatible**: every existing `.gmod` runs unchanged and the `IterRslt`/`SimuRslt` shapes are frozen.

## Quickstart

From a directory containing `<model>.gmod`, with the package on the path:

```matlab
addpath('src');
gdsge_codegen('HL1996');          % parse -> IR -> generate iter_/simulate_/mex_ + compile
IterRslt = iter_HL1996;           % solve the model
SimuRslt = simulate_HL1996(IterRslt);
```

Or the one-call orchestrator, which codegens, solves (with caching), and simulates:

```matlab
eq = gdsge('HL1996');             % eq.model / eq.IterRslt / eq.SmltRslt
```

A worked example model is `tests/HeatonLucas1996/HL1996.gmod`.

## Authoring a model

A `.gmod` file is a **MATLAB superset**: grid and parameter lines (e.g.
`k = exp(linspace(log(kMin),log(kMax),kPts));`) are real MATLAB, `eval`'d during parsing to obtain
the grids and parameters. The model-specific syntax sits on top.

### Blocks

| Block | Purpose |
|---|---|
| `parameters` | scalar model parameters (real MATLAB assignments) |
| `var_state` | endogenous state variables (define grids) |
| `var_shock` | exogenous shocks; with `shock_num` and `shock_trans` |
| `var_policy[/_init]` | unknowns solved per grid point (`_init` = for the `model_init` warm-up) |
| `var_aux[/_init]` | auxiliary variables computed from the solution |
| `var_interp` | variables carried across iterations by interpolation |
| `var_tensor` | ndgrid tensors over states (MATLAB-side: feed bounds / initial interp) |
| `var_output` | variables written to `IterRslt` |
| `var_others` | extra names passed through to `IterRslt` |
| `model[/_init]` | the residual equations; `model(<cond>)` for conditional regions |
| `equations` | the equation block (supports `if/else/end` on the shock value) |
| `simulate` | the simulation step (policy updates, recorded paths) |
| `pre_model`, `pre_iter`, `post_iter` | hooks run around the iteration |
| `pre_jac_code`, `post_jac_code` | hooks around the Jacobian evaluation |

### Declarations

- `inbound[/_init]` — solver bounds per unknown; `adaptive(factor)` widens them across iterations.
- `initial` — initial guesses for the unknowns.
- `shock_num` — number of exogenous shock states; `shock_trans` — the transition matrix.

### Operators, markers, built-ins

- Reductions over shocks: `GDSGE_EXPECT{...}`, `GDSGE_MIN{...}`, `GDSGE_MAX{...}`, `GDSGE_PROD{...}`.
- Interpolation: `GDSGE_INTERP_VEC(...)` and the **primed** form `GDSGE_INTERP_VEC'(...)` (and
  named interp calls).
- Markers: a trailing `'` means **next-state / future**; `[N]` declares a **shock-indexed array**
  variable (e.g. `w1n[8]` — one value per shock), indexed `w1n(j)`.
- Built-ins available in model code: `GDSGE_Iter`, `TASK`, `OUTPUT_CONSTRUCT_CODE`.
- **Conditional regions** (Phase 9b): `model(X>0)` / `model(X==0)` declare independent square systems
  selected per grid point by the state guard; inside `equations`, `if/else/end` on the shock value
  fills an equation slot by branch.

### Macros

`#define`, `#for`, `#foreach`, `#if`, `#mat{}`, `#strcat_comma`, `include`, `cinclude` — expanded
before parsing. `cinclude` injects a C++ include into the generated MEX.

## Options reference

Set options in the gmod (e.g. `TolEq=1e-8;`) or pass them to `iter_<model>(struct(...))`.

| Option | Meaning | Default |
|---|---|---|
| `UseAutoDiff` | `1` = adept autodiff; `0` = SymPy analytic Jacobian; unset = auto-detect (see [Choosing the C++ Jacobian backend](#choosing-the-c-jacobian-backend)) | auto |
| `USE_SPLINE` / `USE_ASG` / `USE_PCHIP` | interpolation method | spline |
| `INTERP_ORDER` | spline interpolation order | per method |
| `ExtrapOrder` | extrapolation order outside the grid | per method |
| `SIMU_RESOLVE` / `SIMU_INTERP` | simulate by re-solving each period vs interpolating policies | resolve |
| `TolEq` | convergence tolerance on the iteration metric | `1e-6` |
| `MaxIter` | maximum iterations | model-set |
| `SaveFreq` | iterations between intermediate `IterRslt` saves (`inf` = never) | model-set |
| `PrintFreq` | iterations between progress prints | model-set |
| `NumThreads` | OpenMP threads (`0` = all cores) | `0` |
| `UseMexRandomize` | `1` = run the randomized-restart loop inside the C++ MEX (C++ RNG); `0` = legacy MATLAB restart path | `1` |
| `MexRandomizeBatch` | minor-iterations (restart trials) per MEX call before returning to MATLAB for diagnostics/cap/exit | `100` |
| `MexRandomSeed` | fixed seed for the in-MEX restart RNG (per-grid-point, thread-independent → reproducible) | `0` |
| `num_samples`, `num_periods` | simulate panel size (passed to `simulate_<model>`) | per call |
| `SimuPrintFreq`, `SimuSaveFreq` | simulate print/save frequency | per call |
| `AsgMinLevel`, `AsgMaxLevel`, `AsgThreshold` | adaptive sparse grid refinement controls | per model |
| `AsgOutputMaxLevel`, `AsgOutputThreshold` | ASG output-grid controls | per model |

### Randomized-restart reproducibility (`UseMexRandomize`)

By default the solver's randomized restarts run inside the C++ MEX (cartesian spline/pchip
iter path), using a per-grid-point counter RNG seeded by `MexRandomSeed`. This makes a run
**reproducible from the seed alone and independent of `NumThreads`** — MATLAB's `rng` state no
longer affects iter results (it still drives `simulate`'s shock draws). The MEX runs up to
`MexRandomizeBatch` restart trials per call before returning to MATLAB for diagnostics, the
`MaxMinorIter` cap, and the converge/exit decision. Set `UseMexRandomize=0` to fall back to the
legacy MATLAB restart path (which draws MATLAB `rand`); `UseAdaptiveBoundInSol=1` also forces
that fallback. The iter goldens were captured under the default (`MexRandomSeed=0`); the ASG and
`simulate` (SIMU_RESOLVE) paths keep the MATLAB restart path.

## Choosing the C++ Jacobian backend

The C++ MEX solver can build its per-grid-point Jacobian two ways: **adept autodiff** (records
a tape) or a **SymPy analytic Jacobian** (CSE-shared symbolic gradient). The backend is decided
at codegen time by this precedence (`gdsge.codegen.resolveBackend`):

1. **In-gmod `UseAutoDiff`** — authoritative when set. `UseAutoDiff=1;` ⇒ adept; `UseAutoDiff=0;`
   ⇒ SymPy.
2. **`GDSGE_BACKEND` environment variable** — `adept`/`autodiff` ⇒ adept, `sympy` ⇒ SymPy,
   `auto` (or unset) ⇒ auto-detect. An unrecognized value errors.
3. **Auto-detect** (nothing pinned): **Python present ⇒ SymPy**, with a codegen-time fallback to
   adept if SymPy generation fails; **Python absent ⇒ adept**.

Two consequences worth knowing:

- With nothing pinned, **installing the `uv` Python env flips the default to SymPy**. To keep a
  Python-free adept build regardless, pin it: `UseAutoDiff=1;` in the gmod, or
  `GDSGE_BACKEND=adept`. (The test harness pins `GDSGE_BACKEND=adept` for exactly this reason.)
- The auto-mode SymPy→adept fallback only fires in **auto** mode. An **explicit** `UseAutoDiff` or
  `GDSGE_BACKEND` pin does **not** fall back — it errors if that backend cannot generate.

The driver prints a one-line `Backend: …` message every run so the resolved choice is visible.

### Using the SymPy backend

Add `UseAutoDiff=0;` to the gmod (or supply it as the first line). This switches the C++ backend to
an analytic Jacobian generated with SymPy. Requirements and limits:

- Needs the `uv`-managed Python env: `uv sync --project pyext`.
- **Spline and ASG interpolation** are supported. **pchip** under `UseAutoDiff=0` still raises a
  clear error (`gdsge:codegen:sympyInterpUnsupported` / `gdsge:codegen:unsupported`) — see
  `docs/deferred-features.md` (it has no C++ backend at all).
- Results match the autodiff backend within numerical tolerance (proven end-to-end on six models).

The default autodiff path needs **no Python** — only MATLAB and a configured C++ MEX compiler.

## Reading results

- **`IterRslt`** — the converged solution: `Iter`, `Metric`, `shock_num`, `shock_trans`,
  `var_state`, `var_policy`, `var_aux`, `var_interp`, plus the `var_output` names. Field shapes are
  frozen.
- **`SimuRslt`** — the simulated panel: `shock` (the shock-index path) and one field per recorded
  variable, each `num_samples × (num_periods+1)`.

## Deferred features

A few constructs and option combinations are not yet generated and fail fast with a clear error
rather than producing broken code. See `docs/deferred-features.md` for the full map (C++-body
`var_tensor`, the SymPy backend for pchip, `GenCodeSegment`, …).
