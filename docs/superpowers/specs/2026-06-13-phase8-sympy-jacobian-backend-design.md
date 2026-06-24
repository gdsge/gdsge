# Phase 8 — SymPy Analytic-Jacobian Backend — Design Spec

- **Date:** 2026-06-13
- **Status:** Approved (design); pending plan
- **Owner:** Wenlan Luo
- **Phase:** 8 (see `PROGRESS.md`; parent design
  `docs/superpowers/specs/2026-06-11-refactor-gdsge-design.md` §9.3, §10)

---

## 1. Context

The default GDSGE C++ backend (Phases 5/7) generates a MEX whose model function records its
Jacobian with the **adept** autodiff library: the model body lowers reductions to explicit
shock loops over `adouble` accumulators (`emitModelBody.m`), and `model_eval.tpl.cpp` runs an
adept recording pass to fill `GDSGE_jac`. It works and is the backward-compat anchor; it needs
no Python.

This phase adds the **alternative backend** promised by the parent design §9.3: the same MEX
outer structure (POPN stack locals, OpenMP per-grid-point, CoDoSol), but the Jacobian comes
from **SymPy + CSE** as stack-allocated `double helper_i;` temporaries instead of adept. The
genuinely novel piece is **reduction-to-loop fusion with the chain rule closed through
interpolant derivatives** (parent §10) — deferred to this phase precisely because it is the
hardest part.

Two facts discovered during design make the spline path tractable:

1. **The spline kernel already exposes a pure-`double` gradient evaluator**
   (`MatlabInterpEval.h:41`, `search_eval_with_grad_vec_at_array(arrayIdx, xSite, evalResult,
   grad, cellOfSite)`). So `interp'` — the derivative of each interpolant w.r.t. its evaluation
   point, the input the chain rule in §10 needs — is available with **zero kernel work**. (The
   adept path already uses the adouble sibling; the autodiff `prepare_space` template fills
   `GDSGE_INTERP_GRAD` from it.)
2. **A finite-difference Jacobian path already exists** (`#ifdef USE_FINITE_DIFF` in
   `model_eval.tpl.cpp`, `cfg.finiteDiffDelta`). It gives a free *third* witness for
   cross-checking adept vs SymPy.

## 2. Goals

1. Generate a correct analytic Jacobian for the **spline-interpolation** corpus models
   (HL1996, safe_assets, Mendoza2010, GLSW), selected by an additive codegen-time option,
   with the default (autodiff) path producing **byte-identical** output as today (no golden
   regeneration).
2. Solve the **reduction-to-loop fusion** problem (parent §10): differentiate each reduction
   body once, symbolically, and emit a single fused shock loop accumulating both value and
   gradient row, with the chain rule closed through interpolant derivatives.
3. Keep the symbolic work confined to Python/SymPy; MATLAB owns the loop structure and all
   IR→C++ lowering (Approach C, §5).
4. Cross-check the analytic Jacobian against adept and finite-difference, per grid point, to
   tolerance — plus full end-to-end golden match for every gated model.
5. Disable `cxx;…end;` raw-C++ blocks in this mode (cannot differentiate opaque C++); raise a
   clear error.

## 3. Non-goals

- **ASG and pchip interpolation** under the SymPy backend (→ Phase 8b; §8). ASG's
  interp-derivative availability is unproven, unlike the spline `with_grad` path.
- **`var_tensor`** (already a deferred error, Phase 7e).
- **Performance improvement.** Parity is the bar (parent non-goal #4); any speedup is reported
  informationally, not gated.
- **Changing the MATLAB backend** (`iter_*.m`/`simulate_*.m`). Data packing and the CoDoSol
  cascade are identical; only the generated `.cpp` model function differs.

## 4. Backend selection & guards

Selection is a **codegen-time IR option** (consistent with `interpMethod`, and because the
choice bakes a different `model_eval` into the MEX — there is no runtime switch):

- New IR field `options.jacobianBackend ∈ {'autodiff','sympy'}`, default `'autodiff'`. Added to
  `gdsge.ir.schema` + validator. Surfaced as an additive gmod/codegen option (resolved in
  `resolveOptions`). **The default path is unchanged**, so every existing model still generates
  identical autodiff code and no goldens regenerate.
- `generateCxx` branches on it; the MATLAB backend is untouched.

Guards (clear errors, no wrong code — mirroring the Phase 7c/7e deferral style):

| condition | error id |
|---|---|
| `sympy` + any `cxx;…end;` block | `gdsge:codegen:cxxUnderSympy` |
| `sympy` + `interpMethod ∈ {asg,pchip}` | `gdsge:codegen:sympyInterpUnsupported` |
| `sympy` + `var_tensor` | (already blocked, Phase 7e invariant) |
| `sympy` requested but Python/SymPy unavailable | `gdsge:codegen:sympyPythonUnavailable` (hint: `uv sync --project pyext`) |

## 5. Architecture — Approach C (MATLAB loops, Python bodies)

The Jacobian is **symbolic forward-mode AD with CSE and loop fusion**. MATLAB owns the loop
*structure* and an IR→C++ lowering it already has; Python's only job is symbolic
differentiation of one RHS at a time. The two backends literally share the loop shape and
differ only in how the body-gradient is produced (parent §10's stated goal).

### 5.1 The pyenv bridge

The MATLAB↔Python boundary is **one JSON string in, one JSON string out** (deeply-nested AST
cells convert awkwardly to `py.dict`; JSON reuses the IR's existing `jsonencode`/`jsondecode`
machinery and keeps the Python side unit-testable without MATLAB).

MATLAB side, `gdsge.codegen.sympy.*`:
- `ensurePyenv()` — idempotent; points `pyenv` at the uv venv interpreter
  (`pyext/.venv/Scripts/python.exe` on Windows), inserts `pyext/` into `py sys.path`, verifies
  `import gdsge_sympy` + sympy. Returns false (→ guard error / test skip) if the env is unsynced.
  Callable from `matlab -batch`.
- `callSympy(fnName, requestStruct)` — `jsonencode` → `py.gdsge_sympy.<fn>(jsonChar)` → `char`
  → `jsondecode`. The single chokepoint for every MATLAB→SymPy call.

Python side, new package `pyext/gdsge_sympy/` (path-imported; no install step):
- `__init__.py` — entry functions taking/returning JSON strings.
- `ast_to_sympy.py` — IR node (`num`/`name`/`primed`/`binop`/`unop`/`call`) → SymPy, with a
  fixed map for gmod built-in calls (`exp`,`log`,`pow`,…). Inverse of `emitExpr`.
- `diff_body.py` — the §5.3 core.
- `ccode_helpers.py` — the `hans` `my_ccode` recipe (custom `user_functions`, `pow`).

Tracked: package source + `uv.lock`; git-ignored: `.venv`.

### 5.2 The gradient registry

Every computed name — unknown, interp result, reduction result, reused aux — carries a
**gradient row**: a sparse map `slot → C++ expression` of its partials w.r.t. the unknowns
`GDSGE_x[0..NX-1]`. MATLAB maintains the registry while walking model statements in order.

Seed rows, built by MATLAB from the IR slot layout:

- **Unknown** (policy var `c1`, array element `w1n[j]`): row `{own slot: 1}`.
- **Shock realization** (`g'`,`d'`), **param**: empty row (constant w.r.t. unknowns).
- **Interp result** `psn_j = interp(w1n_j)`: chain row `{slot(w1n_j): GDSGE_INTERP_GRAD[psn]}`
  (one term per interp arg if multi-dimensional). This is where `interp'` enters.

### 5.3 The per-statement contract (the one thing Python does)

Given an RHS AST + the list of its free symbols, Python returns `∂RHS/∂s` for **each** free
symbol `s`, all CSE'd together with the value, as C++ (`helper_i` temporaries + one expr per
partial). Python never sees slot indices — it differentiates pure expressions.

MATLAB forms the new name's gradient row by the chain rule, assembling the sparse combination
in C++ (it knows which slots are nonzero):

```
grad_row[new][m]  =  Σ_s  (∂RHS/∂s) · grad_row[s][m]      // only nonzero s,m
```

The **two chain-rule levels** — interp-of-unknown inside a reduction loop, and reduction-result
`es1` outside — are the *same* registry mechanism at different scopes; no special-casing.

### 5.4 Worked example — HL1996 `es1 = GDSGE_EXPECT{ g'^(1-γ)·(c1n'/c1)^(-γ)·(psn'+d')/ps }`

- Inside the shock loop, the body's free symbols are `g_j, c1n_j, c1, psn_j, d_j, ps`. Python
  returns `∂body/∂{each}` (CSE'd with the value).
- MATLAB combines via the registry: `c1,ps` use seed rows `{slot:1}`; `c1n_j,psn_j` use their
  interp chain rows `{slot(w1n_j): dc1n_j}` / `{slot(w1n_j): dpsn_j}`; `g_j,d_j` drop out.
  Result `dbody[m]` is nonzero only at `m ∈ {slot(c1), slot(ps), slot(w1n_j)}`.
- The loop accumulates `es1 += trans_ij·body` and `des1[m] += trans_ij·dbody[m]`. `es1`
  registers with row `des1[·]`, which the downstream equation differentiation consumes the same
  way.

## 6. The emitted C++

A new emitter `gdsge.codegen.cxx.emitModelSympy` + template `model_sympy.tpl.cpp` (the autodiff
`emitModel`/`model.tpl.cpp`/adept block of `model_eval.tpl.cpp` are untouched). Computes
residuals + Jacobian in **one all-`double` pass — no adouble, no `_stack`**:

```cpp
GDSGE_FUNC_MODEL_1_double(double* GDSGE_x, double* GDSGE_f, double* GDSGE_jac) {
  // unpack unknowns + data (POPN, unchanged)
  // interp -> search_eval_with_grad_vec_at_array (value + GDSGE_INTERP_GRAD, double)
  // <Python helper_i temporaries; fused loops>     // value always; gradient if(GDSGE_jac)
  // GDSGE_f[i] = residual_i;   JAC(i,m) = drow_i[m];
}
```

- **Gradient storage = stack arrays, sparse-written.** `NX` denotes the unknown count, equal
  to `NUM_EQUATIONS` for the square system (`emitEquations` already asserts this). Each
  materialized intermediate gets `double d<name>[NX] = {0}`; only nonzero slots assigned.
  `JAC(i,m)=GDSGE_jac[i + NUM_EQUATIONS*m]` — the same column-major layout the adept/finite-diff
  paths use, so the solver side is identical. Reuses the existing `MAXDIM>MAX_STACK_DIM` → `vector<double>`
  fallback for large-NX models. Gradient accumulation guarded by `if (GDSGE_jac)` so value-only
  calls stay cheap.
- **Interp evaluator:** the SymPy path's double evaluator calls `search_eval_with_grad_vec_at_array`
  (value + `GDSGE_INTERP_GRAD`), replacing the value-only `search_eval_vec_at_array` the
  autodiff double-evaluator uses.

**Per-reduction loop shapes** (all share the scaffolding; only the accumulate differs — §10):

| kind | value | gradient row |
|---|---|---|
| EXPECT | `acc += trans_ij·body` | `dacc[m] += trans_ij·dbody[m]` |
| PROD | `acc *= body_j` | `dacc[m] = dacc[m]·body_j + acc·dbody_j[m]` *(before updating acc)* |
| MIN/MAX | `if(body_j ≷ acc) acc=body_j` | in that branch, `dacc[m]=dbody_j[m]` (extremal-`j` subgradient) |

MIN/MAX copy the argmin/argmax shock's sub-gradient — the same branch adept's tape would
record, so cross-check parity holds at smooth points (the kink is measure-zero; both backends
pick a side). The corpus exercises **EXPECT** (all 4 models) + **MIN** (safe_assets); PROD/MAX
share the lowering and get a synthetic unit test but are not corpus-gated.

## 7. Testing (TDD, layered)

- **Python unit tests** (`pyext/tests/`, `uv run --project pyext pytest`, no MATLAB):
  `ast_to_sympy` round-trips; `diff_body` vs hand-computed derivatives; CSE/`ccode` output
  shape; the chain-rule combination. The symbolic core pinned cheaply.
- **MATLAB unit tests:** `ensurePyenv` gating, `callSympy` JSON round-trip, `emitModelSympy`
  snapshot, and every §4 guard error.
- **Cross-check gate (per spline model) — the decisive one:** a debug evaluation mode returns
  the analytic Jacobian at a supplied point; the gate asserts **SymPy Jacobian == adept Jacobian
  == finite-difference reconstruction**, per grid point, to tolerance. Three independent
  witnesses directly verify the §5/§6 fusion. (Exposing the Jacobian likely needs a small
  task-template debug hook analogous to the existing `debugEvalOnly`.)
- **End-to-end gate (per spline model):** the SymPy-backend MEX converges to the same golden
  (Iter/Metric within window; `var_policy`/`aux`/`interp` within tolerance; simulate shock
  path) — `jacobianBackend='sympy'` variants of the four existing `tEndToEnd*` gates.
- **Python-availability gating:** all SymPy gates `assumeTrue(ensurePyenv())`, so the default
  correctness loop stays green on a Python-free machine (CLAUDE.md constraint) and runs fully
  where the uv env is synced.

## 8. Scope boundaries & deferrals

- **ASG → Phase 8b.** ASG interp-derivative availability is unproven (no spline-style
  `with_grad` path confirmed); CaoKS2016 / Bianchi2011 wait for it.
- **pchip → unsupported** (parity with the autodiff backend, which also defers pchip in
  `generateCxx`).
- **`cxx` blocks → hard error** (§4). **`var_tensor` → already blocked** (Phase 7e).
- **PROD/MAX → implemented** (share the lowering) but synthetic-tested only; EXPECT + MIN are
  corpus-gated.

## 9. Risks

- The sparse gradient-registry + per-body CSE is the novel core → mitigated by the three-way
  Jacobian cross-check (adept + FD witnesses) per model and the Python unit tests.
- `pyenv`-at-uv-venv under `matlab -batch` can be finicky on Windows → idempotent `ensurePyenv`
  + clear unavailability guard + skippable gates.
- Stack pressure from many materialized `d<name>[NX]` rows on larger-NX models (Mendoza) →
  sparse writes + the existing `vector<double>` fallback; worth measuring.
- MIN/MAX non-differentiability at the kink → parity argument with adept; synthetic test at
  smooth points.

## 10. Deliverables

1. IR option `jacobianBackend` (schema + validator + `resolveOptions`); default unchanged.
2. `gdsge.codegen.sympy.{ensurePyenv,callSympy}` + the `pyext/gdsge_sympy/` package
   (`ast_to_sympy`, `diff_body`, `ccode_helpers`).
3. `gdsge.codegen.cxx.emitModelSympy` + `templates/cxx/model_sympy.tpl.cpp`; `generateCxx`
   branch + guards; double with-grad interp evaluator.
4. Debug Jacobian-return hook (task template) for the cross-check.
5. Tests: Python unit suite; MATLAB unit + guard tests; per-model cross-check + end-to-end
   `sympy`-variant gates (HL1996, safe_assets, Mendoza2010, GLSW), Python-availability-gated.
6. `PROGRESS.md` + changelog update; note Phase 8b (ASG) as the deferred sub-phase.
```
