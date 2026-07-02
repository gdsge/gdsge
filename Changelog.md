# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
