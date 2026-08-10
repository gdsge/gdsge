# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.3] - 2026-08-10

### Added

- **Global var_aux via `GDSGE_DER(v, s)`** — declare a plain `var_aux` and define it in
  the var_aux declaration section (`MPK = GDSGE_DER(Y, K);`): the analytic partial
  derivative, w.r.t. one state, of the spline fitted to a policy/aux each iteration,
  evaluated at every grid point per shock. Constructed after each iteration (never
  assigned in the model block); feeds back through an ordinary var_interp update
  (`MPK_interp = MPK;`); chaining earlier globals gives higher derivatives. Tracked in
  the IR as `variables.auxGlobal` (minor bump to 2.1.0), merged into `IterRslt.var_aux`.
  New `interp_der_mex` kernel (analytic spline gradient). Cartesian splines only;
  ASG and sized (`[N]`) variables raise namespaced errors. Example + gates:
  `tests/rbc/` (spline MPK matches the closed form to ~2.5e-9 on the grid).
- An ordinary var_aux that is never assigned in the model block is now a clear parse
  error (`gdsge:parser:auxNotAssigned`) instead of an opaque C++ compile failure.

### Changed

- The simulate MEX (`simulate_<model>_mex`) is generated for every MEX-expressible
  spline model regardless of `SIMU_RESOLVE`/`SIMU_INTERP`, making the compile cache
  invariant to the toggle: flipping it only rewrites `simulate_<model>.m`, with no
  recompilation once cached.
- `SIMU_INTERP` models whose simulate block is not MEX-expressible now warn at codegen
  (`gdsge:codegen:simuInterpFallback`) instead of silently using the slower per-period
  MATLAB loop. Under `SIMU_RESOLVE`, global var_aux listed in `var_simu` warn that they
  are interpolated from `IterRslt.output_interp`.

## [0.2.2] - 2026-07-31

Driven by `TechBubblesV4nblogKRAv3`, an asset-pricing model with growth shocks whose
simulation needs thousands of randomized solver restarts per period.

### Fixed

- Resolve retries are unbounded again by default (`MaxMinorIter = inf`, as in the
  pre-refactor toolbox). The previous default of 200 did not fail loudly: `simulate`
  wrote the unconverged solution into `SimuRslt` and fed it forward as the next
  period's state, so a model needing more restarts produced plausible-looking but
  wholly wrong paths. Set `MaxMinorIter` yourself to restore a bound — doing so also
  arms the existing hard error on the iteration path.
- `simulate` no longer crashes with `MATLAB:badsubscript` when it reports unconverged
  samples. The report was grid-shaped and mis-indexed for any model with more than one
  state; it now names the sample and its own shock and state values.

### Changed

- `simulate` reports resolve activity once per `SimuPrintFreq` checkpoint
  (`resolve: N/M periods retried, ... worst period P`) instead of a full diagnostics
  block every few retries. Set `DiagnoseMinorIter` to re-enable the per-solve block.
- The diagnostics heartbeat backs off geometrically (at `DiagnoseMinorIter`, then each
  doubling) rather than firing at a fixed interval, so a long solve reports about
  `log2(N)` times instead of `N/interval`.

### Added

- An install self-check runs once per session and reports a truncated or empty toolbox
  source file as `gdsge:install:corrupt`, naming every damaged file and the remedy,
  instead of MATLAB's `MATLAB:scriptNotAFunction` pointing at an internal file.
- `simulate` batches randomized solver restarts inside the MEX (`UseMexRandomize`, on by
  default), instead of one MEX call per restart. Restart-heavy simulations get roughly
  1.1–1.8× faster depending on `num_samples` — the gain is largest when samples are few,
  because that is when per-call overhead dominates. Results are unchanged (agreement with
  the pre-refactor toolbox stays at ~1e-9); set `UseMexRandomize=0` for the MATLAB restart
  loop, which draws from MATLAB's RNG instead.

## [0.2.1] - 2026-07-03

Driven by a large CCP (discrete-choice) model: 13 states, adaptive sparse grids, and
~2,500 statements of nested-`max` logsumexp algebra — the largest model the toolbox
has compiled to date, and the first to run multithreaded on Windows.

### Added

- The gmod declaration region accepts arbitrary MATLAB again, as the pre-refactor
  toolbox did: multi-output assignments (`[a,b] = ndgrid(...)`), bare function
  calls, struct-field assignments, `for`/`if` blocks, and declaration blocks that
  re-open (`parameters ... var_shock ... parameters ...`). A new guard catches
  misspelled declaration keywords (`paramters`) at parse time with a suggestion
  instead of a confusing downstream error.
- Relational operators (`<`, `>`, `<=`, `>=`, `==`, `~=`) in model expressions compile on
  the C++ (adept) path with C++-correct precedence and parenthesization (`~=` emits as
  `!=`). Auto backend selection treats comparisons like `max`/`min`/`abs` — non-smooth,
  so SymPy is never auto-picked for models that use them.
- MATLAB row-vector literals in model expressions: `[a b c]` (space- or comma-separated)
  broadcasts elementwise through operators and functions and collapses to a scalar under
  `sum` / `prod` / `max` / `min`, matching the pre-refactor symbolic expansion.
- Common-subexpression hoisting in the parser: repeated subexpressions are computed once
  into generated locals, shrinking generated code, compile time, autodiff tape work
  (~10x faster iteration on CCP-style models), and stack usage. Values are unchanged.

### Fixed

- Auto backend selection now picks adept autodiff for models that call `max`/`min`/`abs`:
  SymPy analytic Jacobians of non-smooth models blow up symbolically.
- Windows: very large models compile (`/bigobj`) and run multithreaded — generated stack
  frames now fit the 1 MB OpenMP worker stacks, and the OpenMP runtime is pinned so
  `clear mex` after a parallel solve no longer crashes MATLAB.
- `0*x` terms fold to zero as the old symbolic pipeline did, so models that reference
  undeclared names inside zero-multiplied terms keep running.
- The parser no longer warns about `initial`/`var_interp` assignments placed after the
  `model_init` block — the canonical layout when initial interp values come from
  `model_init` results.
- The SymPy backend copes with a Python interpreter that other code (e.g. `startup.m`,
  another toolbox) loaded into MATLAB first: an out-of-process interpreter without SymPy
  is replaced with the GDSGE environment automatically; an in-process one cannot be
  replaced, so the toolbox now says exactly that (restart MATLAB) instead of
  "Python not detected", and `gdsge_setup_sympy` validates the session at setup time.

## [0.2.0] - 2026-07-02

Ground-up refactor of the toolbox. Drop-in backward compatible: every existing `.gmod`
runs unchanged, the public API is preserved, and `IterRslt` / `SimuRslt` shapes are frozen.
See [refactor_summary.md](refactor_summary.md) and [PROGRESS.md](PROGRESS.md) for detail.
The original toolbox is archived on the
[`legacy-0.1.x`](https://github.com/gdsge/gdsge/tree/legacy-0.1.x) branch.

### Added

- SymPy analytic-Jacobian backend, auto-selected when a Python env is present.
- Clearer diagnostics on solver errors and unsolved equations.

### Changed

- Rebuilt around an explicit, versioned JSON intermediate representation (IR): a modular
  MATLAB parser emits the IR; interchangeable backends generate the runtime MATLAB and the
  C++ MEX solver.
- Replaced Intel MKL with Eigen for linear algebra (MKL dependency and license removed).

## [0.1.5]

### Added

- Support for MacOS silicon processors
- A set of new macros to write symmetric models. (See the multi-country example.)

### Changed

- Replace Intel MKL with Eigen for linear algebra and own-implemented linear and cubic spline interpolation. This enables silicon processor support. Still keep MKL linear algebra for Windows version for performance.

## [0.1.2]

### Added

- Editor plugins
- Return more model information (now parameter values in .gmod files are returned)

### Changed

- Add instructions for compiling a gmod file in README
- Add instructions for setting up compilers on a macOS
- Fix include dependence bugs
- Fix bugs in using var_tensor
- Fix bugs in initiating var_interp with _INTERP affix

## [0.1.1] - 2022-09-30

### Changed

- Remove dependencies on v2struct in convert_to_interp_eval_array.m. This should improve performance

## [0.1.0] - 2022-09-09

### Added

- Changelog
- cinclude <> 
- Blocks with snippets added to mex_xxx.cpp: start_loop, finish_loop
- Blocks with snippets added to iter_xxx.cpp: pre_mex_call, pre_jac_code, post_jac_code

### Changed

- Ending of cxx; block is changed to endcxx;

### Removed

- Figure in Mendoza example that produces broken titles
