# GDSGE Toolbox Refactor — Design Spec

- **Date:** 2026-06-11
- **Status:** Approved (design); pending implementation
- **Owner:** Wenlan Luo
- **Source task:** `task20260611_refactor_gdsge.md`

---

## 1. Context & motivation

GDSGE is a MATLAB toolbox that solves global DSGE economic models. The user writes a
`.gmod` file (a domain-specific language that is, in practice, a *superset of MATLAB*),
and the toolbox compiles it into MATLAB iteration code plus a C++ MEX solver, then runs
policy iteration and simulation.

The current implementation works and is computationally fast, but the *generation*
pipeline is hard to maintain and extend:

- **The parser** (`base_package/gdsge/source/gdsge_parser.m`) is a single 3027-line
  function that walks straight from `.gmod` text to emitted MATLAB/C++ strings using
  many fixed regular expressions and string replacements. Extending it is very hard
  because there is no intermediate structure — every change risks a fragile regex.
- **The generated MATLAB** relies on the `v2struct` package to pack/unpack workspace
  variables. `v2struct` copies variables (overhead) and can silently overwrite them.
  There is no error checking, no reporting when a model fails to solve, and the console
  output is unstructured and uninformative.
- **Differentiation** of the model equations is done in C++ by the `adept` automatic
  differentiation library (tape-based). It works, but the user wants an alternative
  that emits a **stack-allocated analytic Jacobian** (via SymPy) for maximum efficiency.

The refactor keeps the proven, fast runtime architecture (per-grid-point parallel solve,
stack-allocated C++ locals, the `CoDoSol` nonlinear solver) while replacing the brittle
generation layer with a clean, modular, testable pipeline built around an explicit
**Intermediate Representation (IR)**.

## 2. Goals

1. Keep the entire **running stack in MATLAB + C MEX**. The default path must require
   **no Python**. Python (managed by `uv`) is used only by an *alternative* symbolic
   backend and during refactoring/tooling.
2. Rewrite the parser as a **modular MATLAB package** with distinct, independently
   testable stages, replacing fixed-regex shortcuts with a real tokenizer + AST for the
   model body. Organize parser and templates better — without over-engineering.
3. Introduce a **language-neutral IR** (JSON) that the parser emits and that all code
   generators consume. The IR is the contract; backends never re-parse `.gmod`. This
   makes future backends in other languages (Julia, Python, C++) possible.
4. Rewrite the generated MATLAB to **drop `v2struct`**, use direct struct dot-access
   (a view, not a copy), add **error checking**, **report unsolved problems**, and print
   **structured, informative** progress.
5. Add an **alternative symbolic backend**: SymPy-generated **stack-allocated analytic
   Jacobian** replacing adept autodiff. Requires Python; `cxx;…end;` raw-C++ blocks are
   disabled in this mode. Solve the reduction-operator-to-loop fusion problem.
6. Preserve **maximum efficiency** and the **stack-allocation convention** for C++.
7. Develop test-first (**TDD**) with a command-line MATLAB test harness.
8. Remain **backward compatible**: every existing `.gmod` file must continue to run, and
   the public API + result struct shapes must not change.

## 3. Non-goals

- Rewriting the numerical kernels (`CoDoSol`, `adept`, `Eigen`, ASG, splines). These are
  vendored and moved, not rewritten.
- Changing the `.gmod` language surface. We add nothing the user must adopt; existing
  syntax keeps working. (New *options* to select the SymPy backend are additive.)
- Building the non-MATLAB backends now. The IR is designed to allow them later; we do not
  implement them in this project.
- Performance *improvements* beyond preserving current performance (the analytic-Jacobian
  backend may improve it, but parity is the bar).

## 4. Constraints

- **MATLAB:** R2025b is the working version on this machine (R2020b–R2024a also present).
- **Python:** managed by `uv` (0.7.18 present), pinned local version, SymPy dependency.
- **No-Python-default:** parsing + autodiff path must run with MATLAB + MEX only.
- **Backward compatibility:** all six current test models must reproduce results within
  numerical tolerance of the old toolbox.
- **Stack allocation:** generated C++ uses stack-allocated locals as much as possible.

## 5. Architecture — IR-centered pipeline

The single most important change is inserting a typed IR between parsing and codegen.

```
   model.gmod ──▶ PARSER (MATLAB) ──▶ IR (model.gdsge.json)
                       │                      │
                       │            ┌─────────┴──────────┐
                       ▼            ▼                    ▼
                (no Python)  MATLAB CODEGEN       C++ CODEGEN
                             iter_*.m / sim_*.m   mex_*.cpp
                                                  ├─ autodiff (adept)   ← default, no Python
                                                  └─ analytic Jacobian  ← alt, needs Python/SymPy
                                                     (SymPy + CSE)
```

**Why IR-centered** (over the alternatives considered):

- A *strangler/shim* (wrap the old parser, swap one stage at a time) keeps a runnable
  system but preserves the regex internals we want gone and yields no IR.
- A *Python/Lark parser for everyone* (the `hans` approach) is cleaner to write but forces
  Python into the default path, violating the no-Python-default constraint.
- *IR-centered with a MATLAB parser* is more upfront work but is the only option that
  satisfies both "default needs no Python" and "IR for other languages," and it makes the
  two Jacobian backends drop-in interchangeable.

**The parser lives in MATLAB** because `.gmod` is a MATLAB superset: grid/parameter lines
such as `k = exp(linspace(...))` and `shock_trans = shock_trans ./ repmat(...)` are real
MATLAB that must be `eval`'d to obtain grids and parameters. A non-MATLAB parser cannot do
this without re-implementing MATLAB.

## 6. Package layout

The new package lives at the repository root. `base_package/` and throwaway refactoring
code are git-ignored.

```
src/+gdsge/            MATLAB package (namespaced)
  +parser/             preprocess, lex, block-split, declaration parser, model-expr parser
  +ir/                 IR struct definition, JSON encode/decode, validators
  +codegen/            IR -> MATLAB and IR -> C++ generators
  +runtime/            runtime helpers used by generated code (error reporting, printing)
templates/
  matlab/              cleaned .m templates
  cxx/                 cleaned .cpp/.h templates
pyext/                 uv-managed Python: SymPy analytic-Jacobian codegen (alt path only)
include/               vendored C++ (adept, Eigen, codosol, asg, splines) — moved as-is
tests/
  run_tests.m          builds + runs the MATLAB unittest suite, writes JUnit/TAP, sets exit code
  run.ps1              shell wrapper: find R2025b, matlab -batch, propagate exit code
  <Model>/golden/      captured IterRslt/SimuRslt reference results
docs/
  superpowers/specs/   this spec
  notes-for-agents.md  orientation for future agents
  ir-schema.md         IR schema reference
CLAUDE.md              project orientation for Claude sessions
PROGRESS.md            living progress log, phase by phase
```

## 7. Parser design

A modular MATLAB pipeline, each stage independently testable:

1. **Preprocess** — macro expansion: `#define`, `#for`, `#foreach`, `#if`, `#mat{}`,
   `#strcat_comma`, `include(...)`, `cinclude(...)`/`cinclude<...>`. Each macro is its own
   function with its own tests. Deprecated-keyword rewriting (`GNDSGE`→`GDSGE`, etc.).
2. **Lex + block-split** — segment the file into declaration lines and named blocks
   (`pre_model`, `model`, `model_init`, `simulate`, `pre_iter`, `post_iter`,
   `pre_jac_code`, `post_jac_code`, etc.).
3. **Declaration parser** — `parameters`, `var_state/shock/policy/aux/interp/tensor/output/
   others`, `var_policy_init`, `var_aux_init`, `inbound`/`inbound_init` (incl.
   `adaptive(factor)`), `initial`. MATLAB setup code (grids, `shock_trans`, params) is
   `eval`'d to obtain numeric values. Produces structured variable tables with the flat
   solution/data/aux **slot layout** (each variable's index range in the packed arrays).
4. **Model-expression parser** — a real recursive-descent parser over a tokenizer for the
   `model;…equations;…end;` body. Produces an **AST** for: assignments, the `'`
   future/next-state marker, `GDSGE_INTERP_VEC` (and primed form), the reduction operators
   `GDSGE_EXPECT/GDSGE_MIN/GDSGE_MAX/GDSGE_PROD` (with optional custom transition matrix),
   future-indexed equations, and the equation list. This replaces today's most fragile
   regex/`str2sym` code.
5. **Semantic analysis** — name resolution, slot assignment, validation: square system
   (#equations == #unknowns after array/prime expansion), bounds present for every policy,
   reserved-word checks, interp-argument consistency, interpolation-method invariants
   (exactly one of `USE_SPLINE`/`USE_ASG`/`USE_PCHIP`; ASG excludes `var_tensor`; PCHIP
   single-state; `INTERP_ORDER ∈ {2,4}`).
6. **Emit IR** — serialize the populated IR struct to JSON.

## 8. IR design

The IR is a **versioned**, language-neutral JSON document describing *what to compute*,
not *how*. Top-level sections (illustrative; schema finalized in Phase 1):

- `irVersion`, `modelName`, `options` (all toolbox flags resolved to values).
- `shocks`: names, count, `shock_trans` matrix, realized values per shock.
- `states`: names, grids, interpolation method/order.
- `variables`: policy / aux / interp / tensor / output / others — each with name, length
  (scalar or array), and **flat slot index range** in the packed solution/data/aux arrays.
- `bounds`: per policy variable lower/upper expressions, adaptive factors.
- `interp`: interpolation objects, their state arguments, initial expressions, update rules.
- `model`: ordered list of typed statements —
  - `assign{ target, exprAst }`
  - `interpCall{ targets[], primed, argAsts[] }`
  - `reduction{ kind: EXPECT|MIN|MAX|PROD, target, bodyAst, transRef }`
  - `equation{ exprAst, primed }`
- `modelInit`, `simulate`, and the verbatim hook blocks (`pre_iter`, `post_iter`,
  `pre_model`, `cxx`, …) carried as opaque text where appropriate.

Expression ASTs are a small node set (number, name, primed-name, binary/unary op, call,
index). The MATLAB struct form and the JSON form round-trip via `jsonencode`/`jsondecode`,
validated by `+ir` validators.

## 9. Backends

### 9.1 MATLAB runtime (`iter_*.m`, `simulate_*.m`)

- **No `v2struct`.** Codegen knows every variable name from the IR, so it emits explicit
  named locals and direct struct dot-assignment (`IterRslt.var_policy.c1 = ...`), which is
  a view in MATLAB — no copy, no silent overwrite.
- **Error checking + reporting.** After each MEX solve, inspect `CoDoSol` exit flags and
  per-point residual norms; if any grid point failed, report how many, which states/shocks,
  and the worst residual. A non-converging solve ends with a clear diagnostic.
- **Structured printing.** A clean iteration table (iter, metric, max-residual,
  #unconverged, elapsed), gated by existing `PrintFreq`/`NoPrint` options.
- **Unchanged public surface.** `iter_<model>(options)`,
  `simulate_<model>(IterRslt, options)`, and the `IterRslt`/`SimuRslt` field shapes stay
  identical.

### 9.2 C++ MEX — default (adept autodiff)

Generated from IR, preserving today's structure: per-grid-point problems solved in parallel
under OpenMP; data popped into **stack-allocated locals** via `POPN`/`POPNARRAY`;
`CoDoSol::solve` as the solver; `adouble` recording the Jacobian. This is the
**backward-compatibility anchor** and needs no Python. The rewrite here is mostly mechanical
(clean templates fed by IR), so its output can be diffed against the old generator.

### 9.3 C++ MEX — alternative (SymPy analytic Jacobian)

Selected by an additive option flag. Same outer structure (POPN locals, codosol, OpenMP,
stack allocation), but the Jacobian comes from **SymPy + CSE** instead of adept:

- IR expression ASTs → SymPy; gradients of each residual w.r.t. each unknown computed
  symbolically, run through `cse()`, emitted as `double helper_i;` stack temporaries via
  `ccode()` (the `hans` recipe, generalized — see `base_package/hans/hans_parser.py` and
  `template/vfi_template_sympy.h`).
- `cxx;…end;` raw-C++ blocks are **disabled** in this mode (cannot differentiate opaque
  C++); a clear parser error is raised if present.
- Interp evaluators expose, per call, both the value **and its derivative w.r.t. each
  interpolation argument**, so the chain rule closes.

## 10. Reduction-to-loop fusion (the "subtle" problem)

A reduction such as
`es1 = GDSGE_EXPECT{ g'^(1-γ)*(c1n'/c1)^(-γ)*(psn'+d')/ps }`
is mathematically `es1 = Σ_j trans(i,j)·body_j`, where `body_j` mixes **current unknowns**
(`c1, ps`) with **per-shock primed inputs**: `g',d'` (shock realizations indexed by `j`) and
`psn',c1n'` (interpolants evaluated at the future state, which is itself a policy unknown
`w1n`).

The bridge — used by both backends, and what makes the SymPy backend possible:

1. Parse the reduction body into an AST in which primed inputs are **per-shock symbolic
   placeholders** (`g_j`, `psn_j`, …). Differentiate the body **once**, symbolically, w.r.t.
   each unknown `x_m`. Where a primed input is itself an interpolant of an unknown
   (`psn_j = interp(w1n_j)`), the chain rule yields
   `∂psn_j/∂x_m = interp'(w1n_j)·∂w1n_j/∂x_m`, with `interp'` supplied by the interp
   evaluator's derivative output.
2. Emit a **single loop over `j = 1..shock_num`** that accumulates both the value and the
   gradient row, stack-allocated and CSE'd:

   ```cpp
   double es1 = 0.0, des1[NX] = {0.0};
   for (int j = 1; j <= shock_num; ++j) {
     /* load g_j,d_j; eval psn_j,c1n_j and their arg-derivatives; CSE'd body & dbody */
     es1 += trans_ij * body_j;
     for (m) des1[m] += trans_ij * dbody_j[m];
   }
   ```

   `EXPECT` accumulates with `trans`; `PROD` via the product rule; `MIN`/`MAX` track the
   extremal `j` and copy its sub-gradient. The accumulated value + gradient feed the residual
   and the analytic Jacobian assembly.

The autodiff backend uses the same loop shape but lets adept record the gradient instead of
the symbolic `dbody`, so the two backends share the lowering and only differ in how the
gradient is produced.

## 11. Error handling & reporting

- **Parser:** precise, located errors (which block/line) for: malformed declarations,
  non-square systems, missing bounds, reserved-word collisions, conflicting interpolation
  options, `cxx` blocks under the SymPy backend.
- **Runtime:** per-iteration convergence diagnostics; explicit unsolved-problem reports;
  no silent fallthrough on failure.
- **Codegen:** generated C++ must compile; snapshot tests guard structural regressions.

## 12. Testing & TDD

- **Harness:** native `matlab.unittest`. `tests/run_tests.m` builds a `TestSuite`, runs it
  through a `TestRunner` with `XMLPlugin.producingJUnitFormat` (+ TAP), writes
  `tests/results/junit.xml`, and calls `exit(any([results.Failed]))`. Shell entry:
  `matlab -batch "run_tests"` (exit 0 = pass). `tests/run.ps1` wraps this for CI/agents.
- **Path & isolation policy (explicit, never persistent).** The toolbox source — old or
  new — is **never added to the MATLAB path persistently**. Each test/run explicitly
  `addpath`s the source it wants (old `base_package/gdsge/source` for golden capture; new
  `src/` for the new pipeline) or calls with an explicit path, and `rmpath`s on teardown.
  This is required because old and new share **generated-file names** (`iter_<model>.m`,
  `simulate_<model>.m`, `mex_<model>`) and the **public flat entry points**
  (`gdsge.m`, `gdsge_codegen.m`) have the same names in both. The public flat names are
  preserved deliberately for backward compatibility (existing `test.m` files call
  `gdsge_codegen('HL1996')`); they are *thin shims* that forward to the namespaced
  implementation. Isolation between old and new therefore comes from **process + path
  control, not renaming**: (1) old and new internal modules are separated — the new
  *internal* modules live in a **`+gdsge` package namespace** (`gdsge.parser.*`,
  `gdsge.codegen.*`), so only the deliberate flat shims sit at top level; (2) golden
  capture (old) and new-pipeline runs execute in **separate `matlab -batch` processes**
  with controlled working directories and exactly one source `addpath`ed, so neither the
  shims nor the generated files ever clobber each other. The harness owns all
  `addpath`/`rmpath`; nothing depends on the user's saved path.
- **Unit tests first** for each stage: macro preprocessing, declaration parsing,
  expression-AST round-trips, IR schema validation, reduction lowering, codegen snapshots.
- **Golden integration tests:** capture `IterRslt`/`SimuRslt` from the *old* toolbox
  (Phase 0), then assert the new pipeline reproduces each model within tolerance
  (relative/absolute bounds on `var_policy`/`var_aux`/`var_interp` arrays, the convergence
  metric, and iteration count within a small window). Valid across both backends because
  both solve the same equations to the same `TolEq`.

## 13. Backward-compatibility surface

Must remain stable:

- **Public functions:** `gdsge(model)`, `gdsge_codegen(model[,options])`,
  `iter_<model>(options)`, `simulate_<model>(IterRslt, options)`.
- **Result structs:** `IterRslt` fields (`var_state`, `var_policy`, `var_aux`, `var_interp`,
  `var_others`, `params`, `shock_num`, `shock_trans`, `pp`, `GDSGE_PROB`, `Iter`, `Metric`,
  …) and `SimuRslt` fields.
- **gmod feature inventory** (full surface to support): blocks `parameters`, `var_state`,
  `var_shock`, `var_policy[/_init]`, `var_aux[/_init]`, `var_interp`, `var_tensor`,
  `var_output`, `var_others`, `model[/_init]`, `equations`, `simulate`, `pre_model`,
  `pre_iter`, `post_iter`, `pre_jac_code`, `post_jac_code`, hook blocks; declarations
  `inbound[/_init]` with `adaptive(factor)`, `initial`, `shock_num`, `shock_trans`;
  operators `GDSGE_EXPECT`, `GDSGE_MIN`, `GDSGE_MAX`, `GDSGE_PROD`, `GDSGE_INTERP_VEC`
  (and primed `'` form); the `'` future marker; `[N]` shock-indexed array variables;
  built-ins `GDSGE_Iter`, `TASK`, `OUTPUT_CONSTRUCT_CODE`; macros `#define`, `#for`,
  `#foreach`, `#if`, `#mat{}`, `#strcat_comma`, `include`, `cinclude`; options
  `USE_SPLINE/USE_ASG/USE_PCHIP`, `INTERP_ORDER`, `ExtrapOrder`, `SIMU_RESOLVE/SIMU_INTERP`,
  `TolEq`, `SaveFreq`, `PrintFreq`, `NumThreads`, `AsgMaxLevel`, `AsgThreshold`, etc.

## 14. Risks

- **Old toolbox must build here** to capture goldens (needs a MATLAB-configured C++
  compiler). Phase 0 verifies this *first*; if it fails, pause and decide (user supplies
  goldens, or fix the toolchain) before building on top.
- **MATLAB parser effort:** writing a real recursive-descent parser in MATLAB is more work
  than regex. Mitigated by the vertical-slice approach (prove it on one model first) and a
  deliberately *minimum* parser that grows only as models require.
- **SymPy reduction fusion** is the most novel piece; it is deferred to Phase 8 and
  cross-checked against autodiff goldens.
- **Numerical reproducibility:** floating-point and iteration-count differences require
  tolerance-based (not bit-exact) golden comparison with sensibly chosen bounds.

## 15. Phased plan

Vertical slice first (HeatonLucas1996, autodiff backend), then widen. Each phase is its own
spec → plan → agent cycle, tracked in `PROGRESS.md`.

- **Phase 0 — Infrastructure (this session).** git init; `.gitignore`; package skeleton;
  uv Python project; `CLAUDE.md`; `PROGRESS.md`; `docs/notes-for-agents.md` + IR schema
  stub; unittest runner + `run.ps1`; capture goldens (HL1996 first, then all six if the old
  toolbox builds).
- **Phase 1 — IR schema + MATLAB scaffolding.** Versioned IR (JSON schema + MATLAB struct),
  encode/decode round-trip, validators, tests.
- **Phase 2 — Parser front-end.** Macros + lexer + block-split + declaration parser →
  partial IR; tests against HL1996 declarations.
- **Phase 3 — Model-expression parser.** Tokenizer + recursive-descent AST (primes, interp
  calls, reductions, equations) + semantic analysis + slot layout → complete IR for HL1996.
- **Phase 4 — MATLAB codegen.** IR → `iter_HL1996.m` / `simulate_HL1996.m` (no v2struct,
  error reporting).
- **Phase 5 — C++ codegen (autodiff).** IR → `mex_HL1996.cpp` / `compile_*.m`; clean
  templates; stack allocation preserved.
- **Phase 6 — Vertical slice green.** HL1996 end-to-end; compile; run; match golden.
  Architecture proven.
- **Phase 7 — Widen for backward-compat.** Bring the other five models on (ASG, `var_tensor`,
  `model_init`, pre/post_iter, macros, adaptive bounds, custom `shock_trans` in reductions);
  each model a golden test until all six are green.
- **Phase 8 — SymPy analytic-Jacobian backend.** Add the alternative backend + reduction
  fusion; cross-check against autodiff goldens; disable `cxx` in this mode.
- **Phase 9 — Polish.** Docs, IR version freeze, performance check vs old, cleanup.

## 16. This session's deliverables (Phase 0)

1. `git init` and a `.gitignore` that excludes `base_package/`, scratch refactoring code,
   generated artifacts (`mex_*`, `iter_*.m`, `simulate_*.m`, `*.mexw64`, `*.cache`,
   `IterRslt_*.mat`, `SmltRslt_*`), uv venv, and editor cruft.
2. Package skeleton directories (Section 6).
3. uv Python project (pinned Python + SymPy) under `pyext/`.
4. `CLAUDE.md`, `PROGRESS.md`, `docs/notes-for-agents.md`, IR schema stub.
5. MATLAB unittest runner (`tests/run_tests.m`) + `tests/run.ps1`.
6. Verify the old toolbox builds/runs; capture HL1996 golden (and the rest if possible).
7. Commit.

The detailed Phase 1 plan is produced via the writing-plans skill after this spec is
approved.
