# Architecture / developer guide

How the refactored GDSGE toolbox is put together, for maintainers and for agents extending it. The
end-user view (authoring `.gmod` models) is in `docs/user-guide.md`; this document is about the
internals.

## The pipeline

```
model.gmod
   │  gdsge.parser            (MATLAB front-end)
   ▼
<model>.gdsge.json  ──────────  the IR: the contract between front-end and backends
   │
   ├─ gdsge.codegen.generateMatlab ──▶  iter_<model>.m, simulate_<model>.m   (thin; call gdsge.runtime)
   │
   └─ gdsge.codegen.generateCxx ─────▶  mex_<model>.cpp, compile_<model>.m
                                          ├─ adept autodiff   (default; no Python)
                                          └─ SymPy analytic Jacobian   (UseAutoDiff=0)
```

The IR is the single contract. **Backends never re-parse the gmod** — they consume the validated
IR struct (or its JSON round-trip). The driver `gdsge.codegen.codegen` orchestrates parse → encode
JSON → generateMatlab/generateCxx → cache-gated MEX compile. The flat public entry points
(`gdsge_codegen.m`, `gdsge.m`, `asg.m`) are thin shims over the `+gdsge` packages.

## Module map

**`gdsge.parser`** — gmod → IR, in stages:
`preprocess` → `expandMacros` (`#define/#for/#foreach/#if/#mat/#strcat_comma/include/cinclude`) →
`splitBlocks` (extracts `parameters`, `var_*`, `model[/_init]`, `equations`, `simulate`, … and the
`model(<cond>)` regions) → `splitStatements` → `parseVarDecls` / `parseDeclarations` (slot layout) →
`evalSetup` (runs gmod MATLAB to obtain grids/params) → `parseSimulate` / `resolveOptions` /
`resolveOutputs` → the model-body parser `tokenize` → `parseExpr` (full MATLAB precedence) →
`parseModel` (statement classification, inline-reduction hoisting) → `analyzeModel` (name
resolution, square-system check, interp arity). `assemblePartialIR` + `parseFrontEnd` assemble the
whole IR.

**`gdsge.ir`** — the contract machinery, all driven by one declarative descriptor:
`schema` (single source of truth) → `validate` (shape/refs/slots) → `canonicalize` (type-aware
normalization) → `encode`/`decode`/`roundtrip` (JSON, orientation/non-finite quirks handled) →
`isequalIR` (golden comparisons) → `regionView` (Phase 9b region adapter) → `gendoc` (emits
`docs/ir-schema.md`, no-drift tested). AST nodes live under `+node`.

**`gdsge.codegen`** — IR → generated files:
`codegen` (driver + cache via `needsCompile`/`ensureAsgMex`/`ensureSplineConstructMex`), `generateMatlab` (+ the `+mat`
emitters), `generateCxx` (+ the `+cxx` emitters), the shared `dataLayout` (MATLAB packer + C++
`POP*` macros from one descriptor), `assertSupportedIR` (defense-in-depth deferred-error invariant),
`codeWriter`/`writeText`, `initView` (model_init reuse). The C++ emitters (`+cxx`): `emitExpr`
(AST→C++ printer), `emitModel`/`emitModelBody`/`emitEquations`/`emitTask`/`emitArgUnpack`/`emitPop`/
`emitDeclare`/`emitAux`/`emitInterp`/`emitCompile`, `modelVars`/`modelBody` (region-agnostic body
helpers), `lowerCondition` (region/branch guard → C++ boolean), `fillTemplate`/`readTemplate`
(template engine), and the SymPy model emitters under `+cxx/+sympymodel`. The SymPy bridge
(`+sympy`): `ensurePyenv` (an **out-of-process** pyenv at the uv venv —
`ExecutionMode='OutOfProcess'`, to avoid an in-process native crash) + `callSympy` (one JSON
in/out to `pyext/gdsge_sympy`).

**`gdsge.runtime`** — the unit-tested library the thin generated `.m` files call:
`unpackOptions` (whitelist + informative errors), `solveProblems` / `solveProblemsAsg` (the MEX
caller-workspace contract, resolve cascade), `applyWarmUp`, `constructSplines`, `computeMetric`,
`ensurePath`, `printIterProgress`, `reportUnconverged`.

**`src/kernels/`** — vendored flat MEX kernels. `interp_construct_mex` is the fused
spline/linear **constructor**: one call builds every cartesian interpolant from all `var_interp`
`Values` (passed via a struct, no concat) and writes the eval-order `GDSGE_SPLINE_VEC`
**directly** (no MATLAB permute), plus the per-variable natural-order pp structs — replacing the
per-variable `myppual` + `convert_to_interp_eval_array` rebuild in `constructSplines`'s hot loop.
Construction-only and self-contained (the `mkl_start[_with_extrap]` driver over
`mkl_dummy_interp.h`'s `construct_cubic_spline_notaknot`/`construct_linear_interp`); bit-exact
with the legacy path (the `tFusedConstruct*` A/B gates assert it). `extrap_order==4` reuses plain
cubic construction (the boundary cubic extrapolates naturally — no extra knots). Compiled
cache-gated by `ensureSplineConstructMex` (flags from the owner's `compile_myppual.m` + `/O2`,
MSVC `/fp:precise`). The simulate / `output_interp` evaluation path uses a **stacked
uniform-order** layout built by `interp_construct_mex` and evaluated by the generic
`interp_eval_mex` (both on the double-only `include/interp_eval_double.h`); cartesian
`SIMU_INTERP` runs the whole period loop in a generated `simulate_<model>_mex`. `myppual.m`,
`myppual_mex`, and `convert_to_interp_eval_array.m` are retired.

## How each backend consumes the IR

- **Shared data layout** — `gdsge.codegen.dataLayout` emits, from one descriptor, both the MATLAB
  side that packs grid-point inputs into the flat `GDSGE_DATA` array and the C++ `POPN(var)` /
  `POPNARRAY(var,len)` macros that pop them back into stack locals. Stack allocation is preserved
  throughout (`double helper_i;` temporaries).
- **Autodiff (default)** — `emitExpr` is a direct AST→C++ printer (no Symbolic Toolbox); the model
  function is `adouble`, the Jacobian is adept-recorded, `CoDoSol::solve` drives the per-grid-point
  solve.
- **SymPy analytic Jacobian** (`UseAutoDiff=0`) — MATLAB owns the fused shock-loop structure and a
  gradient registry (name → sparse `{slot,expr,templated}` rows, forward-mode chain rule); SymPy
  (`pyext/gdsge_sympy`, out-of-process via `pyenv`) differentiates each body and returns value+partials
  as shared-CSE C++. The reduction-fusion chain rule closes through the spline kernel's
  `search_eval_with_grad_vec_at_array` (pure double). Spline interpolants only today; ASG/pchip are
  deferred (Phase 8b — see `docs/deferred-features.md`).
- **Region-agnostic emitters** (Phase 9b) — `model.regions` and tagged `plain|conditional`
  equations are looped only by the **drivers** (`emitModel`, the SymPy `emitTask`, `analyzeModel`);
  the per-body emitters (`modelVars`/`emitModelBody`/`emitEquations`) stay region-agnostic via the
  `modelBody` helper, and `lowerCondition` lowers each region/branch guard to a C++ boolean.

## How to add a model

The corpus discipline (Phases 9a/9b), in order:

1. **Probe first.** Write `scratch/probe_<model>.m` that enumerates every expression form, operator,
   marker, and option the gmod uses. Check each against what the parser and both backends already
   cover (HL1996 / safe_assets / Mendoza / GLSW / Cao2011EZ / CaoNie2016 between them cover most).
   The output is a precise gap list — what's genuinely new vs already supported.
2. **Capture a golden** from the old toolbox (`tests/golden/capture_<model>.m`): one source on the
   path, full-convergence iter + a reduced seeded simulate; save `-v7`.
3. **Spec → plan → implement** the gap list only, bringing **both** backends to parity.
4. **Gate** with `tFrontEnd<Model>` (IR equality), `tGolden<Model>` (integrity), and
   `tEndToEnd<Model>[Sympy]` (public API vs golden, tolerance-based).

## The MATLAB-path golden rule

Old and new share flat public names (`gdsge_codegen`, `iter_<model>`, `simulate_<model>`, generated
file names). **Never** put a toolbox source on a persistent MATLAB path. Each run `addpath`s exactly
**one** source — old `base_package/gdsge/source` or new `src/` — in its **own `matlab -batch`
process** with a controlled `cd`, restoring the path on cleanup (`onCleanup`). Sharing a path makes
old and new shadow each other. The golden-capture scripts (`tests/golden/capture_*.m`) and the perf
worker (`tests/perf/perf_worker.m`) are the canonical pattern.

## Testing model

Native `matlab.unittest`, headless via `matlab -batch "cd('tests'); run_tests"` (exit 0 = all pass;
report `tests/results/junit.xml`). Per-model gates live in plain folders under `tests/<Model>/`
(`tGolden*`, `parser/tFrontEnd*`, `codegen/tEndToEnd*[Sympy]`); utility functions live in the
`+gdsgetest` package. Comparison is tolerance-based (not bit-exact) via
`gdsgetest.compareNumericClose`. SymPy tests `assumeTrue(gdsgetest.sympyAvailable())`, so the default
Python-free loop stays green. Python-side SymPy unit tests: `uv run --project pyext pytest pyext/tests`.

## Pointers

- `docs/ir-schema.md` — full field-by-field IR schema (auto-generated, no-drift tested).
- `docs/ir-versioning.md` — the IR stability policy + changelog.
- `docs/deferred-features.md` — the honest-error map + open deferred sub-phases.
- `docs/notes-for-agents.md` — concrete facts about the old codebase.
- `docs/superpowers/specs/` — the per-phase design specs (the design history of every feature).
