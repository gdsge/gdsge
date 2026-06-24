# GDSGE refactor — summary

A ground-up refactor of the **GDSGE** MATLAB toolbox, which solves global nonlinear DSGE economic
models. The toolbox is rebuilt around an explicit, versioned **JSON intermediate representation
(IR)**: a modular MATLAB parser emits the IR, and interchangeable backends generate the runtime
MATLAB and the C++ MEX solver.

It is **drop-in backward compatible** — every existing `.gmod` runs unchanged, the public API is
preserved, and the `IterRslt` / `SimuRslt` result shapes are frozen.

## Two C++ Jacobian backends

- **adept autodiff** — needs only MATLAB and a configured C++ MEX compiler. No Python.
- **SymPy analytic Jacobian** — analytic Jacobians via an `uv`-managed Python env under `pyext/`.
  Spline interpolants only today (ASG/pchip are autodiff-only — see
  [docs/deferred-features.md](docs/deferred-features.md)).

When `UseAutoDiff` is not set in the gmod, the backend is auto-detected: a configured Python env
selects SymPy, otherwise adept autodiff. Force a backend with `UseAutoDiff=1;` (adept) or
`UseAutoDiff=0;` (SymPy) in the gmod.

## Quickstart

From a directory containing `<model>.gmod`, with the package on the path:

```matlab
addpath('src');
gdsge_codegen('HL1996');          % parse -> IR -> generate iter_/simulate_/mex_ + compile
IterRslt = iter_HL1996;           % solve
SimuRslt = simulate_HL1996(IterRslt);
```

Or the one-call orchestrator (codegen → cache-gated solve → simulate):

```matlab
eq = gdsge('HL1996');             % eq.model / eq.IterRslt / eq.SmltRslt
```

A worked example model is `tests/HeatonLucas1996/HL1996.gmod`. See
[docs/user-guide.md](docs/user-guide.md) for authoring models, the options reference, and reading
results.

## Running the tests

```
matlab -batch "cd('tests'); run_tests"
```

Exit code 0 means all pass; the report is written to `tests/results/junit.xml`. (The convenience
wrapper `tests/run.ps1` needs PowerShell 7 / `pwsh`, which is absent on the dev machine — invoke
`matlab -batch` directly as above.) The alternative SymPy-backend Python tests:
`uv run --project pyext pytest pyext/tests`.

## Repository layout

| Path | What |
|---|---|
| `src/+gdsge/+parser` | gmod → IR front-end (lexer, macro engine, block split, expression parser) |
| `src/+gdsge/+ir` | the IR contract: schema, validate, encode/decode, doc-gen |
| `src/+gdsge/+codegen` | IR → MATLAB + C++ (`+mat`, `+cxx` adept emitters, `+cxx/+sympymodel`, `+sympy` bridge) |
| `src/+gdsge/+runtime` | the unit-tested library the generated `.m` files call |
| `src/{gdsge.m,gdsge_codegen.m,asg.m}`, `src/kernels/` | flat public shims + vendored kernels |
| `templates/cxx/` | C++ MEX templates |
| `pyext/` | `uv`-managed SymPy backend (`gdsge_sympy`) |
| `tests/` | per-model gates (`tGolden*`, `tFrontEnd*`, `tEndToEnd*[Sympy]`) + unit suites |
| `docs/` | design specs, architecture, user guide, IR schema |

## Architecture in one line

`model.gmod → MATLAB parser → JSON IR → { MATLAB codegen, C++ codegen (adept autodiff | SymPy
analytic Jacobian) }`. The IR is the contract; backends never re-parse the gmod.

## Documentation

- [docs/superpowers/specs/2026-06-11-refactor-gdsge-design.md](docs/superpowers/specs/2026-06-11-refactor-gdsge-design.md) — the approved design
- [docs/architecture.md](docs/architecture.md) — developer guide (pipeline, module map, how to add a model)
- [docs/user-guide.md](docs/user-guide.md) — authoring models, options, reading results
- [docs/ir-schema.md](docs/ir-schema.md) — full IR schema (auto-generated)
- [docs/ir-versioning.md](docs/ir-versioning.md) — IR stability policy + changelog
- [docs/deferred-features.md](docs/deferred-features.md) — honest-error map + open sub-phases
- [docs/perf-report.md](docs/perf-report.md) — performance vs the old toolbox
- [PROGRESS.md](PROGRESS.md) — phase-by-phase status
