# Phase 8b (ASG slice) — SymPy analytic-Jacobian backend for ASG

**Date:** 2026-06-15
**Status:** approved (design)
**Relationship:** the ASG half of the deferred **Phase 8b**. pchip (the other half) stays
deferred — see §10.

> **Update (2026-06-15, during implementation):** the **CaoNie2016_asg** capstone (§6) was
> **dropped**. The model does not converge cleanly under the reference (old) toolbox, so it is not
> a reliable golden/gate. The SymPy-ASG slice ships validated by **CaoKS2016 + Bianchi2011_asg**
> end-to-end gates plus the **sympy↔adept↔FD Jacobian cross-check** (all green). Consequence: the
> **ASG × Phase-9b conditional-region** combination is left unexercised by the corpus — no other
> model uses it, so backward compatibility is unaffected, but that intersection is not yet proven
> on either backend (tracked as a future item should a converging conditional-region ASG model appear).

---

## 1. Context & motivation

The refactor ships two C++ Jacobian backends: the default **adept autodiff** and the
alternative **SymPy analytic Jacobian** (`UseAutoDiff=0` → `options.jacobianBackend='sympy'`,
Phase 8). The SymPy backend was proven on the four **spline** corpus models
(HL1996 / safe_assets / Mendoza2010 / GLSW). The **ASG** (adaptive sparse grid) interpolant —
widened to autodiff parity in Phase 7b — has only ever run on autodiff. Today
`gdsge.codegen.generateCxx` raises `gdsge:codegen:sympyInterpUnsupported` when
`jacobianBackend=='sympy'` and `interpMethod=='asg'`.

This sub-phase removes that gap: ASG models run on the SymPy backend with the **same support**
the spline path already has, so the cartesian/spline and ASG paths reach backend parity.

## 2. Goal & non-goals

**Goal.** Lift `sympyInterpUnsupported` for `interpMethod=='asg'`. ASG models generate, compile,
and converge on the SymPy analytic-Jacobian backend, matching their autodiff goldens.

**In scope.**
- SymPy backend support for the ASG interpolant (value + analytic Jacobian via the ASG
  interpolant's derivative).
- Three ASG end-to-end SymPy gates: **CaoKS2016**, **Bianchi2011_asg**, **CaoNie2016_asg**
  (a new corpus model — §6).
- A SymPy↔adept↔finite-difference Jacobian cross-check on an ASG model.

**Non-goals (unchanged behavior).**
- **pchip under any C++ backend** stays deferred (separate gap; §10).
- No IR schema change (`interpMethod` already carries `asg`; `jacobianBackend` already optional).
- No change to the MATLAB ASG drivers `emitIterAsg` / `emitSimulateAsg` or to
  `gdsge.runtime.solveProblemsAsg` — they call the MEX and are backend-agnostic. The autodiff
  ASG output stays byte-identical (the SymPy path is only emitted when `jacobianBackend='sympy'`).
- No in-MEX randomize / resolve / convergence-diagnostics work for ASG — those are the other
  three ASG-parity sub-phases (B/C/D), tracked separately.

## 3. Key enabling fact — the ASG double gradient already exists

The SymPy chain rule closes through the interpolant's derivative. For splines that is the kernel's
`search_eval_with_grad_vec_at_array` (value + gradient, pure `double`); the emitted lambda
`GDSGE_INTERP_VEC_double_grad` calls it and `addInterpCall.m` threads the gradient into the
forward-mode registry.

The ASG interpolant **already computes the same value+gradient in pure double**.
`include/asg_adouble.h::eval_vec_adept(...)` accumulates the value (`p_totalSurplus`) **and** the
gradient (`p_gradient[i_dim] += surplus * gradientProd[i_dim]`) as `double`. The adouble evaluator
`AsgInterpArrayAdoubleEvaluator::eval_vec_adouble(...)` is just a thin wrapper: it calls
`eval_vec_adept`, scales `gradient[i] /= stateRange[i]`, then injects the gradient into the adept
tape via `add_derivative_dependence`. **No new ASG math is required** — only a pure-double entry
point that stops short of the adept injection, plus the codegen wiring.

## 4. Architecture & change surface

Four edits. The downstream SymPy model emitters (`+cxx/+sympymodel/*`, the registry chain rule in
`addInterpCall.m`) are **untouched** because the new ASG path reuses the exact same lambda name and
signature the spline path emits.

### 4.1 `include/asg_adouble.h` — one pure-double method (vendored, additive)

Add to `AsgInterpArrayAdoubleEvaluator` a method mirroring `eval_vec_adouble`'s
`evalGradFlag==1` branch but staying in `double` and skipping `add_derivative_dependence`:

```cpp
void eval_vec_with_grad(int i_array, double* site, double* cell, double* ratio, double* slope,
    double* rslt, double* gradient)            // rslt[numVec]; gradient native [i_vec*numDim+i_dim]
{
    double siteScaled[ASG_MAX_DIM];
    for (int d = 0; d < numDim; d++)
        siteScaled[d] = (site[d] - stateMin[d]) / stateRange[d];
    AsgInterp& interp = interps[i_array];
    search_adept_with_slope(siteScaled, cell, ratio, slope, currentLevel, numDim);   // slope needed for grad
    eval_vec_adept(numVec, cell, ratio, slope, currentLevel, numDim,
                   interp.info, interp.setOfLevelCombinations, rslt, gradient);
    for (int v = 0; v < numVec; v++)
        for (int d = 0; d < numDim; d++)
            gradient[v*numDim + d] /= stateRange[d];   // -> derivative w.r.t. the real site
}
```

Purely additive; mirrors a proven method; no existing path changes. (The VEC evaluator covers
named scalar interp calls too — the spline SymPy path also routes named calls through the VEC
evaluator and picks an index, so no per-interp ASG double-grad method is needed.)

### 4.2 New emitter — SymPy ASG interp section

Make `gdsge.codegen.cxx.emitInterpSympy` ASG-aware (branch on `interpMethod=='asg'`, or add a
sibling `emitInterpSympyAsg` it delegates to). It returns `frag.getCode` (task scope) and
`frag.threadCode` (per-thread), like the spline version:

- **getCode:** the ASG construct block — reuse `templates/cxx/interp_asg_construct.tpl.cpp`
  (converts `GDSGE_ASG_HANDLE` → the evaluator). Identical to the autodiff ASG getCode.
- **threadCode:** a **double-only** scratch block (no adouble scratch) — `GDSGE_ASG_CELL/RATIO/SLOPE`,
  a native-layout gradient scratch, plus the `GDSGE_INTERP_VEC_double_grad` lambda. Author a small
  SymPy-ASG prepare-space rather than reusing `interp_asg_prepare_space.tpl.cpp` verbatim (that one
  declares adouble scratch the SymPy path does not use).

The lambda has the **same signature** the spline SymPy path emits:

```cpp
auto GDSGE_INTERP_VEC_double_grad =
  [&GDSGE_CPP_ASG, &GDSGE_ASG_CELL, &GDSGE_ASG_RATIO, &GDSGE_ASG_SLOPE, &GDSGE_ASG_GRAD_NATIVE]
  (int shockIdx, /* double site0, site1, ... */ double* GDSGE_out, double* GDSGE_grad)
{
    double xSite[] = { /* site0, site1, ... */ };
    GDSGE_CPP_ASG.eval_vec_with_grad(shockIdx-1, xSite,
        GDSGE_ASG_CELL, GDSGE_ASG_RATIO, GDSGE_ASG_SLOPE, GDSGE_out, GDSGE_ASG_GRAD_NATIVE);
    // transpose native [v*DIM+d] -> spline layout [v + NVEC*d] (see §5)
    for (int v = 0; v < NVEC; v++)
        for (int d = 0; d < DIM; d++)
            GDSGE_grad[v + NVEC*d] = GDSGE_ASG_GRAD_NATIVE[v*DIM + d];
};
```

### 4.3 `generateCxx.m` — lift the guard

Change the SymPy interp guard from "spline only" to "spline **or** asg". Keep the pchip rejection
and the `cxx`-hook-under-SymPy rejection. `ensureAsgMex()` is already called for ASG regardless of
backend, and `emitCompile` already appends `-DUSE_ASG -DASG_MAX_DIM/NVEC/LEVEL` purely on
`interpMethod=='asg'` (backend-agnostic) — so the SymPy ASG compile gets the right defines and
kernel for free.

### 4.4 Tests — §7.

**Why this is small:** `emitCompile` (ASG defines), `mex.tpl.cpp` (shared shell), the SymPy
`emitTask` (already assembles its interp section via `emitInterpSympy` and already loops
`ir.model.regions`), and the whole registry/chain-rule layer are all already in place. The new code
is one C++ method + one emitter branch + one guard edit.

## 5. The gradient contract (the single correctness-critical point)

The native ASG gradient and the layout `addInterpCall.m` consumes differ on **two** axes — both must
be reconciled or the Jacobian is silently wrong:

1. **Scaling.** `eval_vec_adept` returns the gradient w.r.t. the **scaled** `[0,1]` site. The real-site
   gradient requires `/= stateRange[d]`. Handled inside `eval_vec_with_grad` (§4.1), mirroring
   `eval_vec_adouble`.
2. **Layout.** `addInterpCall.m` reads `GDSGE_grad[ix + numInterp*dim]` (**vec-fastest, dim-major**).
   `eval_vec_adept` writes `gradient[vec*numDim + dim]` (**dim-fastest, vec-major**). The emitted
   lambda transposes native → spline layout (§4.2).

Both are pure, local transforms with no downstream coupling. The Jacobian cross-check (§7) pins them
exactly — this is the same class of bug (a mis-assembled interp gradient row) that the Phase 8
cross-check caught on the spline side.

## 6. CaoNie2016_asg — the capstone gate (new corpus model)

`base_package/.../CaoNie2016_asg/CaoNie2016.nmod` is the **same model body** as the cartesian
**CaoNie2016** (Phase 9b, green on **both** backends) with only the interpolation axis flipped:
`USE_SPLINE=0; USE_ASG=1; AsgThreshold=1e-4` (no `InterpOrder`). It is the strongest possible
SymPy-ASG gate and the first model to combine ASG with Phase 9b's control structures.

**Why it de-risks rather than expands parser/SymPy risk:**
- The model body — conditional regions `model(X>0)` / `model(X==0)`, `if/else` conditional equations
  on the shock value, `Xp[3]` array policy, `qp(1..3)`/`Xp(1..3)` const-indexing, multi-output
  `GDSGE_INTERP_VEC'`, the `consis/Xp(1)-1` mid-expression case — is already parsed and proven on
  both backends. Zero parser risk; the SymPy machinery already handles every construct.
- An **old-toolbox golden already exists**: `base_package/.../CaoNie2016_asg/IterRslt_CaoNie2016_1077.mat`
  (Iter=1077). The new `.gmod` is the existing `tests/CaoNie2016/CaoNie2016.gmod` plus the ASG option
  diff (flip `USE_SPLINE`/`USE_ASG`, `AsgThreshold`, drop `InterpOrder`; keep the empty `var_tensor;`,
  which already parses).

**The genuinely new axis it exercises — on both backends:**
- **ASG × conditional regions/equations.** ASG (Phase 7b) was only tested on non-conditional models;
  conditional regions (Phase 9b) only on cartesian. The intersection is *structurally* supported —
  autodiff `emitTask` builds its body from the region-aware `emitModel` + ASG `emitInterp`; the SymPy
  `emitTask` already loops regions and will use the new ASG SymPy interp — but **untested**, so
  bringing CaoNie2016_asg up may surface **autodiff-side** ASG-region wiring bugs, independent of any
  SymPy work.

**Therefore CaoNie2016_asg is brought up autodiff-first:** add it as a corpus model
(`tests/CaoNie2016_asg/` with `.gmod` + IR JSON snapshot), get it green on the **autodiff** ASG path
against the golden (fixing any ASG×region wiring), then add its SymPy gate.

## 7. Verification strategy

Mirrors Phase 8's "prove the chain rule, then prove convergence on the corpus."

1. **Jacobian cross-check (strongest).** Clone `tSympyJacobianHL1996.analyticJacMatchesAdeptAndFD`
   for an ASG model (CaoKS2016 is the natural fixture): build both adept + SymPy MEX, evaluate the
   Jacobian at the converged solution via the 6th output (`GDSGE_DEBUG_EVAL_ONLY==2`, already in the
   shared `task.tpl.cpp`), assert **sympy == adept** (rtol 1e-6) and **sympy == finite-difference**
   (rtol 1e-4) per grid point. This directly validates the §5 contract.
2. **Three ASG end-to-end SymPy gates**, each under `UseAutoDiff=0`, matching the existing autodiff
   golden (Iter / Metric + `var_policy`/`var_aux`/`var_interp` within the project's iter tolerance):
   - **CaoKS2016** and **Bianchi2011_asg** — already autodiff-green; isolate pure SymPy-ASG risk.
   - **CaoNie2016_asg** — the capstone (after its autodiff bring-up), green on **both** backends.
3. **Python gating.** All SymPy tests stay `assumeTrue(gdsgetest.sympyAvailable())`, so the default
   Python-free correctness loop is unaffected.
4. **Structural snapshot.** A fast generate-only test asserting the SymPy ASG `.cpp` emits the ASG
   construct + `GDSGE_INTERP_VEC_double_grad` + `JAC(` and no adouble model lambda (mirrors the fast
   `generatesAnalyticModel` test).

The autodiff goldens stay byte-identical (no autodiff codegen path changes); only new tests and the
new corpus model's snapshots are added.

## 8. Risks & mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| §5 gradient contract (scaling/layout) wrong → subtly wrong Jacobian | medium | The Jacobian cross-check (§7.1) pins it exactly, at machine precision vs adept. |
| ASG × conditional-region intersection needs autodiff-side wiring | medium | CaoNie2016_asg brought up **autodiff-first**; if it balloons, split into a follow-on sub-phase (8b-2) without blocking the core (CaoKS2016 + Bianchi2011 are non-conditional). |
| Vendored-header edit destabilizes the autodiff ASG path | low | Addition is purely additive (new method, no edit to existing methods); the autodiff goldens are a regression tripwire. |
| A SymPy construct unexercised by the spline corpus appears in an ASG model | low | All three model bodies already pass SymPy on the cartesian side (CaoNie2016) or are simple reductions+interp (CaoKS2016/Bianchi2011). Implementation step 1 is an IR-coverage probe (reduction kinds, interp/region usage) per the project's "probe before the gap list" habit. |

## 9. Sequencing / phasing

1. **IR-coverage probe** of the three ASG models (reduction kinds, interp arity, regions/conditional
   eqs) — confirm no SymPy construct gap before coding.
2. **SymPy-ASG core:** `eval_vec_with_grad` (§4.1) + the emitter (§4.2) + the guard lift (§4.3).
3. **Gate the core** on CaoKS2016 + Bianchi2011_asg (already autodiff-green) + the Jacobian
   cross-check.
4. **CaoNie2016_asg bring-up:** add the corpus model; green on **autodiff** ASG vs the golden (fix any
   ASG×region wiring); then add the SymPy gate.
5. Full suite green; update `PROGRESS.md` / `docs/deferred-features.md` (the
   `sympyInterpUnsupported` row now covers pchip only).

Each step keeps commits small; the failing test precedes the code.

## 10. Deferred (Phase 8b remainder)

**pchip under the C++ backend** stays deferred. pchip is a *cartesian* interpolant with **no** C++
backend at all (`generateCxx` raises `gdsge:codegen:unsupported` for pchip on *either* Jacobian
backend) — a distinct gap from ASG, with no corpus model exercising it. After this sub-phase, the
`gdsge:codegen:sympyInterpUnsupported` deferral covers **pchip only**; `docs/deferred-features.md`
is updated accordingly.
