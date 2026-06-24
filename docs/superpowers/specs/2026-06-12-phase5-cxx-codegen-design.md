# Phase 5 — C++ Codegen (adept autodiff) Design

- **Date:** 2026-06-12
- **Status:** Approved (design); pending implementation
- **Parent spec:** `2026-06-11-refactor-gdsge-design.md` (§9.2, §10, §15 Phase 5)
- **Scope:** IR → `mex_<model>.cpp` / `compile_<model>.m` for HL1996 (adept autodiff
  backend), cleaned C++ templates, vendored C++ includes, and the shared
  `GDSGE_DATA` layout contract.

---

## 1. Goal

Consume the validated HL1996 IR and generate the C++ MEX source plus its compile
script, preserving the old runtime architecture exactly: per-grid-point problems
solved in parallel under OpenMP, inputs popped into **stack-allocated locals** via
`POPN`/`POPNARRAY`, `CoDoSol::solve` as the solver, adept `adouble` recording the
Jacobian. The new MEX is a **drop-in replacement** for the old binary: same
positional interface, same caller-workspace contract.

Decisions made during brainstorming (with the user):

1. **Full functional gate.** Compile the new `mex_HL1996`, then run the Phase-4
   generated MATLAB files against it to convergence and match the committed
   goldens — the same shape as Phase 4's gate with the roles reversed. A direct
   old-MEX-vs-new-MEX comparison on identical packed inputs is added as a sharper
   diagnostic layer. Phase 6 shrinks to wiring the `gdsge_codegen`/`gdsge` shims.
2. **Templates + emitters.** Cleaned template files in `templates/cxx/` keep the
   static boilerplate; `gdsge.codegen.cxx.emit*` functions generate the hole
   contents from the IR. A direct AST→C++ printer replaces the old
   `str2sym`/`ccode` round-trip (no Symbolic Toolbox). Matches parent spec §9.2
   ("clean templates fed by IR").

The generated C++ stays **structurally close** to the old generator's output
(same macro and naming conventions: `POPN`, `*_GRID` accessors,
`GDSGE_REDUCTION_VALUE<n>`, `GDSGE_INTERP_RSLT`) so it can be debugged by
side-by-side comparison, but is **not byte-diffable**; equivalence is proven
numerically by the MEX-equivalence test and the functional gate.

## 2. Architecture

```
IR (validated struct)
  └─ gdsge.codegen.generateCxx(ir, outDir)
       ├─ fillTemplate(mex.tpl.cpp, holes from emit*)   ──▶ mex_<model>.cpp
       └─ emitCompile(ir)  (compile.tpl.m)              ──▶ compile_<model>.m

mex_<model>.cpp at runtime
  ├─ #includes vendored headers from include/ (adept, codosol, Eigen, splines)
  └─ honors the Phase-4 caller-workspace contract (gdsge.runtime.solveProblems
     is the direct caller; no change on the MATLAB side)
```

### 2.1 Components

```
templates/cxx/                      cleaned static shells (from old code_template/)
  mex.tpl.cpp                       includes, mexFunction dispatch, task selection
  task.tpl.cpp                      OpenMP per-point driver, POPN macro defs, I/O views
  model.tpl.cpp                     adouble/double model lambda skeleton
  call_fmin.tpl.cpp                 CoDoSol invocation + eval-at-solution
  interp_spline_construct.tpl.cpp   spline pp ingestion        (ASG variants → Phase 7)
  interp_spline_get.tpl.cpp         spline evaluator lambdas
  interp_spline_prepare_space.tpl.cpp
  compile.tpl.m                     compile_<model>.m shell

include/                            vendored C++ moved AS-IS from
                                    base_package/gdsge/source/include/:
                                    adept*.h, Eigen/, codosol.h, InterpEval.h,
                                    MatlabInterpEval.h, MatlabPp.h, interp_lite.h,
                                    cubic_spline_notaknot.h, wl_math.h, mm_lite.h,
                                    flat_hash_map.hpp, … (+ the essential_blas
                                    link-time closure found during implementation)

src/+gdsge/+codegen/
  generateCxx.m                     driver: IR → two file texts, writes to outDir
  dataLayout.m                      THE single source of truth for the GDSGE_DATA
                                    layout (ordered descriptor: name, kind
                                    scalar/array, length, grid-accessor flag)
  +cxx/
    emitExpr.m                      AST → C++ expression printer
    emitPop.m                       POP_CODE: POPN/POPNARRAY sequence + *_GRID
                                    #defines, generated from dataLayout(ir)
    emitArgUnpack.m                 ARG_CODE: adouble c1 = GDSGE_x[0]; … (slot layout)
    emitDeclare.m                   DECLARE_CODE: adouble temps (assign targets,
                                    reduction targets, aux, per-shock interp locals)
    emitModelBody.m                 MODEL_CODE: assigns, primed interp-vec loops,
                                    reduction loops
    emitEquations.m                 EQUATION_CODE: GDSGE_EQ[i] = …;
    emitAux.m                       HEADER_AUX_ASSIGN_CODE
    emitInterp.m                    INTERP_GET_CODE / INTERP_GET_THREAD_CODE /
                                    prepare-space from the IR interp section (spline)
    emitCompile.m                   compile_<model>.m: EXTRA_DEF, GDSGE_MAXDIM,
                                    GDSGE_INTERP_ORDER, include path, model name
    fillTemplate.m                  named-placeholder substitution; errors on any
                                    unfilled placeholder or unknown hole name
```

### 2.2 The shared `GDSGE_DATA` layout contract

The layout is currently encoded once on the MATLAB side (Phase 4's
`gdsge.codegen.mat.emitDataPack`). Phase 5 introduces
`gdsge.codegen.dataLayout(ir)` — an ordered descriptor of the packed data
vector (shock count, params, `shock_trans`, shock grids, shock index, state
values) — and refactors `emitDataPack` to consume it while `emitPop` consumes
the same descriptor for the C++ side. The MATLAB packer and the C++ POP
sequence therefore **cannot drift**. The refactor is behavior-preserving and
pinned by the existing Phase-4 snapshot tests plus a new cross-test.

### 2.3 Caller contract (unchanged)

The new MEX reads, via `mexGetVariable("caller", …)`, exactly the names the
Phase-4 runtime already defines in `gdsge.runtime.solveProblems`'s workspace:
the 11 scalars (`MEX_TASK_NAME`, `TolFun`, `TolSol`, `SolMaxIter`,
`NumThreads`, `GDSGE_DEBUG_EVAL_ONLY`, `UseBroyden`, `FiniteDiffDelta`,
`GDSGE_USE_BROYDEN_NOW`, `MEX_TASK_INIT`, `MEX_TASK_INF_HORIZON`) plus
`GDSGE_SPLINE_VEC` and the per-interp `GDSGE_PP_*` structs. Positional
interface unchanged:
`[SOL,F,AUX,EQVAL,OPT_INFO] = mex_<model>(SOL,LB,UB,DATA,SKIP,F,AUX,EQVAL)`.

## 3. Lowering rules

### 3.1 Expression printer (`emitExpr`)

- Direct AST→C++: `^` → `pow(…)`; precedence-aware minimal parenthesization;
  unary minus; function calls passed through; numeric literals printed `%.17g`.
- `name` nodes → C++ locals of the same name; `primed` nodes are only legal
  inside per-shock loop contexts and print as the per-shock local/accessor
  (`name_GDSGE_GRID[GDSGE_iter-1]` locals or `name_GRID(GDSGE_iter)` macro,
  following the old conventions per context).
- No Symbolic Toolbox anywhere (removes the old `str2sym` R2024a workaround
  class of bugs).

### 3.2 Statement lowering (`emitModelBody`)

- **assign** → `target = <expr>;` (`adouble`, declared via `emitDeclare`).
- **interpCall, primed** (`[psn',…] = GDSGE_INTERP_VEC'(w1n')`) → set
  `INTERP_VEC_FLAG` for the participating interps, loop
  `for (GDSGE_iter = 1..shock_num)` evaluating the vector interp at the future
  state and assigning each target's per-shock local from `GDSGE_INTERP_RSLT` —
  the old mechanism, spline path.
- **reduction** → explicit loop with stack accumulator:

  ```cpp
  adouble GDSGE_REDUCTION_VALUE<n> = <init>;   // EXPECT 0, PROD 1, MIN 1e20, MAX -1e20
  for (int GDSGE_iter = 1; GDSGE_iter <= shock_num; GDSGE_iter++) {
    // EXPECT: += trans(shock, GDSGE_iter) * body;  (transRef honored)
    // PROD:   *= body;   MIN/MAX: fmin/fmax against body
  }
  target = GDSGE_REDUCTION_VALUE<n>;
  ```

- **equation** → `GDSGE_EQ[i] = <expr>;` in declared order; future-indexed
  equations expand across the shock loop exactly as the old generator did
  (HL1996's `w1_consis'` → 8 equations; the 19×19 system).

### 3.3 HL1996-shaped holes (real `if`s for Phase 7)

`MODEL_CONDITION = 1` (single model block); no `HAS_INIT` (no `model_init`);
empty `GDSGE_OTHER_INCLUDE` (no `cinclude`); no `cxx` blocks; spline interp
only. Each is a branch on IR presence in the emitters, not a template fork.

### 3.4 `compile_<model>.m` (`emitCompile`)

Multi-platform parity with the old `compile_template.m` (MSVC / MinGW64 /
Intel / g++ / clang branches, OpenMP flags, `-DUSE_MKL` where applicable,
DEBUG branch), modernized in form only: `v2struct(GDSGE_OPTIONS)` replaced with
explicit field reads. Filled holes: `MODEL_NAME`, include path (repo
`include/`), `EXTRA_DEF` from IR options, `GDSGE_MAXDIM`
(= `max(num_policy_total, num_policy_init_total) + 4`, old formula),
`GDSGE_INTERP_ORDER`. Links `essential_blas.lib` on
MSVC as before (vendored with the include closure).

## 4. Error handling

- **`fillTemplate`**: any placeholder left unfilled, or any hole name not
  present in the template → hard error naming the placeholder and template.
  No silent partial substitution (the old generator's failure mode).
- **`generateCxx`**: validates the IR first (existing validator); refuses
  unsupported sections this phase (ASG/PCHIP options, `model_init`, `cxx`
  blocks) with a clear "Phase 7/8" error rather than emitting wrong code.
- **Generated C++ must compile**: the compile gate is part of the test suite,
  not a manual step.

## 5. Testing (TDD order)

1. **`emitExpr` unit tests**: precedence and associativity (incl. MATLAB's
   left-assoc `^` → nested `pow`), unary minus, parenthesization, literals
   round-trip, primed-name printing in loop context, error on primed outside
   loop context.
2. **Per-emitter unit tests**: each `emit*` vs expected text for HL1996 IR
   sections (POP sequence order from `dataLayout`, ARG/DECLARE from slot
   layout, one reduction loop, the primed interp-vec loop, equations block).
3. **`dataLayout` cross-test**: MATLAB-side pack expression (refactored
   `emitDataPack`) and C++-side `emitPop` walk the same descriptor; Phase-4
   snapshots unchanged (pins the refactor as behavior-preserving).
4. **`fillTemplate` tests**: completeness errors both directions.
5. **Snapshot test**: full generated `mex_HL1996.cpp` + `compile_HL1996.m`
   equal committed golden text.
6. **Compile gate**: the generated MEX builds cleanly (MSVC 2022) in an
   isolated work dir with `include/`.
7. **MEX-equivalence test**: `capture_HL1996_mex.m` extended to also stash the
   old generated `mex_HL1996.cpp` / `compile_HL1996.m` as
   `*_old_reference.txt` (debug reference, like the old `.m` files). The test
   drives the OLD and NEW MEX binaries on **identical packed inputs**
   (same SOL/LB/UB/DATA/SKIP from a controlled state): first
   `GDSGE_DEBUG_EVAL_ONLY = 1` (pure residual evaluation — catches any
   expression-printer divergence immediately), then one full solve round;
   outputs (F/SOL/AUX/EQVAL) compared within tight tolerance. Fails loudly
   with "run capture_HL1996_mex first" if the old binary is absent.
8. **Functional phase gate**: isolated process with `src/` + `src/kernels/`;
   `generateMatlab` + `generateCxx` → `compile_HL1996` → `iter_HL1996` to
   convergence → `rng(0823); simulate_HL1996`; compared to the committed
   goldens via `compareNumericClose` (same tolerances and Iter window as
   `tGoldenHL1996` / Phase 4's gate).

## 6. Out of scope (explicit)

- ASG / PCHIP interp templates, `var_tensor`, `model_init` task body, macros,
  `cxx;…end;` user blocks — Phase 7 (emitters already branch on IR presence).
- SymPy analytic-Jacobian backend and reduction-gradient fusion — Phase 8.
- The flat `gdsge_codegen.m` / `gdsge.m` shims and one-command wiring — Phase 6.
- Any change to the vendored numerical C++ — `adept`, `Eigen`, `codosol.h`,
  spline/ASG headers are moved verbatim into `include/`.

## 7. Risks

- **Layout drift** between the MATLAB packer and the C++ POP sequence — killed
  by the shared `dataLayout` descriptor and the cross-test.
- **Expression-printer correctness** (precedence/associativity subtleties) —
  unit tests plus the eval-only MEX-equivalence comparison catch numeric
  divergence at the first residual.
- **Vendor closure**: the exact link-time set for `essential_blas` (`.lib` /
  `.dll`) and any header the templates pull in beyond the listed set is
  discovered at the compile gate; the design accepts adding files to
  `include/`, as Phase 4 did for `src/kernels/`.
- **Old-binary coupling** in the MEX-equivalence test: depends on the local
  `oldmex` artifact; mitigated the same way as Phase 4 (capture script
  regenerates it; the test fails loudly when missing, never skips silently).
