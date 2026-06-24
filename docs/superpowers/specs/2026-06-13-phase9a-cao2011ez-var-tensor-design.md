# Phase 9a Design — Cao2011EZ (`var_tensor` + per-shock `_GRID` indexing)

- **Date:** 2026-06-13
- **Status:** Approved (design); pending plan
- **Parent spec:** `2026-06-11-refactor-gdsge-design.md` (§13 backward-compat surface, §15 Phase 7
  widening); continues the one-model-per-sub-phase precedent of 7a/7b
- **Related:** `2026-06-13-phase7e-var-tensor-error-design.md` (the deferred-error this phase
  partially lifts; §7 "full `var_tensor`" sketch)
- **Owner:** Wenlan Luo

---

## 1. Context

Phases 7a/7b brought the six original corpus models green; Phase 8 added the SymPy
analytic-Jacobian backend for the four spline models. The owner is adding two further example
models from the GDSGE distribution — **Cao2011EZ** and **CaoNie2016** — to the backward-compat
corpus before moving to polish. Both already converge in the old toolbox (their
`IterRslt_*.mat` artifacts exist under `base_package/gdsge/tests/`), so goldens are capturable.

Per the owner's scoping decision, the two models are tackled as **two sequenced sub-phases**,
**Cao2011EZ first** (lower-risk feature set), each with its own spec → plan → golden-capture →
implement cycle. Both C++ Jacobian backends — **adept autodiff** (the default / backward-compat
anchor) and **SymPy** (Phase 8) — are brought to parity on these models.

**This spec covers Cao2011EZ only.** CaoNie2016 (conditional model regions `model(X>0)` /
`model(X==0)`, `if/else/end` conditional equations, `var_others`) gets its own brainstorm/spec
as Phase 9b.

### 1.1 What Cao2011EZ exercises that the refactor lacks

`base_package/gdsge/tests/Cao2011EZ/Cao2011EZ.gmod` — a two-agent Epstein–Zin asset-pricing
model. 1 state `W` (wealth share, 100 pts on `[0,1]`), 2 shocks (`d e1 e2`), 13 unknowns
(`p q c1 c2 lambda1 lambda2 mu1 mu2 kp1 a1 a2` + vector `Wp[2]`), 13 equations (square). It
reuses constructs already supported — a *single* piped `GNDSGE_EXPECT{…|shock_trans1}` against a
square `shock_num×shock_num` parameter matrix (7b precedent), `GNDSGE_MIN{}`, vector policy
`Wp[2]` with `Wp(n)` indexing, the whole-RHS named interp `c1p' = c1Future'(Wp')`, primed shock
`d'`, non-integer powers.

**Four constructs are new** (the three model-body items were confirmed empirically by a parser
probe — see §1.2):

1. **Multiple inline reductions per statement** — the Epstein–Zin Euler equations carry **two**
   `GNDSGE_EXPECT{}` terms each (`eq1 eq2 eq5 eq6`), e.g.
   `eq1 = -1 + … + beta1*(GNDSGE_EXPECT{(v1p')^(1-sigma1)|shock_trans1})^θ * GNDSGE_EXPECT{qp'*…|shock_trans1}/…`.
   `parseModel` today hoists **one** inline reduction per statement and hard-errors
   `gdsge:parser:multipleReductions` on the second.

2. **Embedded named-interp call** — `qp' = qFuture'(Wp')+d';` uses a named interp call as a
   *sub-expression* of a larger primed RHS. `parseModel`/`parseExpr` accept a named interp call
   only as the **entire** RHS and hard-error `gdsge:parser:badInterpCall` otherwise.

3. **`qp_GRID(n)` per-shock indexing** — `qp'` (a primed per-shock assignment) is collapsed by
   `qCC = GNDSGE_MIN{qp'};`, and the consistency residuals index the **uncollapsed** per-shock
   value:
   ```matlab
   consis1 = (qp_GRID(1)*kp1+bp1) / qp_GRID(1) - Wp(1);
   consis2 = (qp_GRID(2)*kp1+bp1) / qp_GRID(2) - Wp(2);
   ```
   `qp_GRID(idx)` is the old toolbox's internal macro for the per-shock value of a future variable
   (`#define qp_GRID(idx) qp_GDSGE_GRID[int(idx)-1]`). The gmod author wrote it explicitly; other
   models (e.g. CaoNie2016) write the bare `qp(n)` form. Both denote the same thing. This **parses**
   today (as a `call` node `fn='qp_GRID'`) but generates a call to a non-existent C++ function —
   a codegen/semantic gap, not a parse gap.

4. **`var_tensor`** — `budget1 budget2 vbudget1 vbudget2 etotal`, each assigned a broadcast
   expression over states/shocks:
   ```matlab
   budget1 = e1+d.*W;      % shock (d, 1×2) × state (W, 1×100) → 2×100 tensor
   budget2 = e2+d.*(1-W);
   vbudget1 = budget1.^(1-beta1);
   vbudget2 = budget2.^(1-beta2);
   etotal = e1+e2+d;       % shock-only → 1×2
   ```
   These are **not used in the `model;` body**. They feed only:
   - `inbound c1 1e-12 etotal;` and `inbound c2 1e-12 etotal;` — a tensor-valued **upper bound**;
   - `initial c1Future budget1;` `initial c2Future budget2;` `initial v1Future vbudget1;`
     `initial v2Future vbudget2;` — tensor-valued **initial interpolation values**.

   Today a non-commented `var_tensor` raises `gdsge:parser:varTensorUnsupported` (Phase 7e). This
   phase lifts that error for the **MATLAB-side** usage Cao2011EZ needs.

### 1.2 Probe findings (2026-06-13, empirical)

Before planning, a one-source-per-process probe (`scratch/probe_cao2011ez.m`) ran the new
parser's model-body stages on Cao2011EZ, testing every body statement and equation in isolation.
Confirmed results (the source of the inventory above):

- `qp' = qFuture'(Wp')+d'` → `gdsge:parser:badInterpCall` (gap 2).
- `eq1 eq2 eq5 eq6` → `gdsge:parser:multipleReductions` (gap 1). `v1 v2` (one `EXPECT` each, wrapped
  in powers) parse cleanly, so the gap is specifically *>1 reduction per statement*.
- `consis1 consis2` (with `qp_GRID(1)`) **parse OK** → gap 3 is codegen-level, not parse-level.
- All other body statements (`qCC`, `c1p'…v2p'`, `bp1 bp2`, `eq3 eq4 eq7 eq8`, `mc`) and all 13
  equations parse cleanly.

This is the authoritative gap list; the plan is built against it.

## 2. Key architectural insight — the new work is parser-level and backend-neutral

Two facts keep "both backends at parity" cheap:

1. **`var_tensor` is MATLAB-side only here.** No tensor variable appears in the `model;` residual,
   so the C++ solver never sees a tensor. The old toolbox pops every tensor into C++ as a `double`
   (`POPN(<tensor>)`, `nData += num_tensor`), but for Cao2011EZ those popped values would be dead.
   Tensors are constructed in MATLAB and consumed in MATLAB (bounds + initial-interp arrays) — **no
   C++ codegen work** for either backend, and no `IterRslt.var_tensor` unless the golden demands it.

2. **Gaps 1–2 are parser-level hoisting.** Multiple inline reductions and embedded interp calls are
   both resolved by **hoisting** to generated statements (the existing single-reduction hoist
   already does exactly this for one reduction; we generalize it and add an interp-call analogue).
   After hoisting, the model body is ordinary reduction/primed-assign statements that **both
   backends already handle** (Phase 5/7a autodiff; Phase 8 SymPy). No backend-specific reduction or
   interp logic is added.

The **only** backend-specific (autodiff + SymPy) work is gap 3 — `qp_GRID(n)`/`qp(n)` per-shock
indexing of a primed-assign target in scalar scope — which mirrors Phase 8's `Re_n(2)` lift and
the autodiff future-var indexed-call branch (`emitExpr.m`), and is verified by the cross-check.

## 3. Scope decisions

**In scope:**

1. **Multiple inline reductions per statement** (parser) — generalize `parseModel`'s single-reduction
   hoist to hoist *each* reduction in a statement to its own `GDSGE_<KIND>_<n>` generated statement,
   dropping the `multipleReductions` error. The nested-reduction guard stays.
2. **Embedded named-interp call** (parser) — hoist a named interp call used as a sub-expression to a
   generated primed target (`GDSGE_INTERP_<n>'`), then reference it as a primed name in the host
   statement. Mirrors the reduction-hoist pattern; drops the whole-RHS-only restriction.
3. **`qp_GRID(n)` / `qp(n)` per-shock primed indexing** — parser normalizes the `_GRID` suffix to
   indexed access of the primed variable; both backends materialize the per-shock primed array and
   index it (autodiff future-var call branch; SymPy `liftIndexed`).
4. **`var_tensor`, MATLAB-side only** — parser captures each tensor's assignment RHS as opaque
   MATLAB text; IR enriches `variables.tensor` from a name list to `{name, expr}` items; the
   MATLAB iter codegen constructs `GDSGE_TENSOR_<t>` broadcast arrays (after the existing
   state/shock `ndgrid` block) and makes tensor names available to bound and initial-interp
   rewrites as `GDSGE_TENSOR_<t>(:)'`.
5. **Tensor-valued bounds and initial-interp values** — `emitBounds` and the initial-interp
   construction rewrite tensor names like they already rewrite state/shock names.
6. **`SIMU_RESOLVE` per-period tensor recompute** — bounds in the per-period re-solve need the
   tensors; recompute them per simulated `(state, shock)` (the 7e sketch's strip-`GDSGE_TENSOR_`
   transform). Trivial for Cao2011EZ (`num_periods=1`) but implemented correctly for generality.
7. **Goldens + gates** for Cao2011EZ on both backends; `PROGRESS.md` bookkeeping.

**Deferred (recorded in `PROGRESS.md`, unchanged from the 7e §7 sketch where applicable):**

- **C++-body `var_tensor` (the `POPN` path)** — tensors that appear inside the `model;` residual,
  requiring `dataLayout`/`emitDataPack`/`emitPop` rows and a popped `double` in C++ (and, under
  SymPy, registration as a zero-gradient constant input). No corpus model needs it; the IR
  enrichment in this phase is forward-compatible with it.
- **`IterRslt.var_tensor` result field** — the old toolbox returns tensors in `IterRslt`; Cao2011EZ
  is verified against the old golden, so whether this field is needed is settled by the golden
  diff (see §6). If the golden carries it, we emit it; if not, deferred. (To be confirmed at
  capture time — see Risk 8.1.)
- **`USE_ASG + var_tensor` guard** — remains a blanket error (ASG is orthogonal; Cao2011EZ is
  spline).

**Out of scope:** CaoNie2016 and its features (Phase 9b); performance work beyond parity
(Phase polish); any change to the `var_tensor`/`_GRID` `.gmod` surface syntax.

## 4. Goal & success criteria

Phase 9a is done when:

1. **Old-toolbox golden captured first** (honoring "verify the original toolbox runs"):
   Cao2011EZ converges in `base_package`, captured in an isolated `matlab -batch` process under
   the one-source-per-process path policy; `tests/Cao2011EZ/golden/` holds `IterRslt.mat` +
   `SimuRslt.mat` + the capture script.
2. **Autodiff end-to-end gate** green: the public API (`gdsge_codegen` → `iter_Cao2011EZ` →
   `simulate_Cao2011EZ`, or the `gdsge`/`gndsge` orchestrator as the old `test.m` uses) reproduces
   the golden — iteration count and metric bit-match; `IterRslt`/`SimuRslt` fields within
   tolerance.
3. **SymPy end-to-end gate** green: same model under `UseAutoDiff=0` reproduces the same golden;
   guarded `assumeTrue(sympyAvailable)` so the default Python-free loop stays green.
4. **Front-end gate** green: `parseFrontEnd(Cao2011EZ.gmod)` yields a schema-valid IR with the
   enriched `variables.tensor`; IR JSON snapshot committed.
5. **The 7e gate stays meaningful**: a non-commented `var_tensor` that this phase does *not*
   support (i.e. one used in the model body) still errors clearly; the lifted path is exactly the
   MATLAB-side usage. (Implementation: the parser gate narrows from "any tensor" to "a tensor
   referenced in the model body", or the codegen invariant takes over for the body case — decided
   in the plan; see §5.5.)
6. All existing gates (HL1996, safe_assets, Mendoza2010, GLSW, CaoKS2016, Bianchi2011, macro /
   orchestration / sympy suites) stay green.

## 5. The change (detail)

### 5.1 Parser — capture tensor assignments

`parseVarDecls.m` already collects `var_tensor` names into `decls.tensorNames` (pure collector).
Today the *assignments* (`budget1 = e1+d.*W;`) fall through Pass B into `decls.setupStmts` because
`budget1` is neither an interp name nor a state name — and they would **fail `evalSetup`** (plain
`eval` of `e1+d.*W` is a dimension-mismatch: `d` is `1×2`, `W` is `1×100`). The fix: in Pass B,
when `lhs ∈ decls.tensorNames`, route the statement to a new `decls.tensorAssign{}` list as
`{name, expr}` (opaque RHS text), **not** to `setupStmts`. Tensor assignments are thus never
plain-eval'd; they are emitted as broadcast `ndgrid` expressions by codegen.

`parseDeclarations.m` assembles `decls.tensorAssign` into the enriched IR `variables.tensor`
(§5.3). The Phase 7e parser gate is re-scoped (§5.5).

### 5.2 Parser — `_GRID` suffix normalization

In the model-body tokenizer/expression parser, an indexed reference `name_GRID(idx)` is
normalized to the same AST node as `name(idx)` — an indexed access of the (primed) variable
`name`. Concretely: recognize the `_GRID` suffix on an indexed call and strip it to the base name,
so `qp_GRID(1)` and `qp(1)` produce identical AST. This keeps the AST/IR free of the old toolbox's
internal macro name. (Validated by parsing both forms to the same node in a unit test.)

### 5.3 IR schema — enrich `variables.tensor`

`variables.tensor` changes from `fList(fText())` (name list) to a list of items, each
`{ name: fText, expr: fText }` — `expr` is the opaque-MATLAB-text assignment RHS (the same
opaque-text world as grids/bounds/setup, **not** the model AST). This is an additive schema
change: empty tensor lists (every existing corpus model) encode to `[]` exactly as before, so
**no existing IR golden regenerates**. The schema descriptor, validator, round-trip, and doc-gen
update together (single-source-of-truth in `gdsge.ir.schema`); the no-drift doc test regenerates
`docs/ir-schema.md`.

### 5.4 MATLAB codegen — tensor construction, bounds, initial, simulate

- **Construction (`emitIter`):** after the existing `GDSGE_TENSOR_<state>` / `GDSGE_TENSOR_shockIdx`
  / `GDSGE_TENSOR_<shock>` `ndgrid` block, emit one line per tensor:
  `GDSGE_TENSOR_<t> = <expr with state/shock/earlier-tensor names rewritten to GDSGE_TENSOR_*>;`
  (tensors may reference earlier tensors — `vbudget1 = budget1.^(1-beta1)` — so emit in
  declaration order and include prior tensor names in the rewrite set). Result: `[shock_num ×
  state-grid]` (or broadcastable) arrays.
- **Bounds (`emitBounds`):** add tensor names to the `names`/`reps` rewrite lists so a bound like
  `etotal` becomes `GDSGE_TENSOR_etotal(:)'`. `etotal` (shock-only, `1×2`) and `budget1`
  (`2×100`) both flatten correctly via `(:)'` against the `[shock × state]` `NPROB` ordering.
- **Initial interp (`initial <interp> <tensor>`):** the initial-interp construction rewrites the
  tensor name to its `GDSGE_TENSOR_<t>` array, sized `[shock_num × state-grid]` — exactly the
  shape of the interp data for one named interpolant.
- **Simulate (`SIMU_RESOLVE`):** emit a per-period tensor recompute (strip the `GDSGE_TENSOR_`
  prefix — scalar/vector ops on the current simulated `state`/`shock`) ahead of the per-period
  bound construction. Cao2011EZ's `num_periods=1` exercises a single step; the code path is
  written for the general case.

### 5.5 Re-scoping the Phase 7e error

Phase 7e raises `gdsge:parser:varTensorUnsupported` for **any** non-commented `var_tensor`. This
phase supports the MATLAB-side usage but **not** the C++-body usage (deferred, §3). The gate is
therefore narrowed rather than removed:

- The **parser** accepts `var_tensor` declarations + assignments unconditionally (they are now a
  supported MATLAB-side construct).
- A new check rejects a tensor **referenced inside the `model;` body** with a clear identified
  error (e.g. `gdsge:parser:varTensorInBodyUnsupported`), pointing to the deferred C++-pop path.
  This is the analysis-time analogue: after `analyzeModel` resolves body names, any name in
  `variables.tensor` that appears as a body reference triggers it.
- The **codegen invariant** `assertSupportedIR` is updated in lockstep: it no longer rejects a
  non-empty `variables.tensor` outright, but still rejects (a) a hand-authored IR whose tensor is
  consumed by the body, and (b) `USE_ASG + var_tensor`.

The plan decides the precise placement (analyze-time vs codegen-time) following the per-stage-owner
precedent; the user-facing guarantee is: **MATLAB-side tensors work; body tensors fail fast with a
clear message.**

### 5.6 Backends — `qp_GRID(n)` / `qp(n)`

`qp'` is a primed assignment already materialized as a per-shock array in both backends (Phase 7a
primed-assign path; Phase 8 SymPy per-shock value+grad arrays). Indexing it by a constant
(`qp(1)`) is the same shape as Phase 8's `Re_n(2)` *array-policy* indexing, but `qp` is a
**primed-interp-derived** array rather than an unknown slot — so the existing path may or may not
generalize. The work:

- **Autodiff (`emitExpr`):** ensure a constant-indexed reference to a primed assignment emits the
  array element (`qp_GDSGE_GRID[idx-1]` analogue) — extend the future-var indexed-call path
  (Phase 7a) if it only covers interp/equation expansion, not standalone scalar use in another
  residual.
- **SymPy (`+sympymodel`):** ensure a constant-indexed primed value is liftable to a scalar with
  the correct per-shock gradient row (analogous to `liftIndexed` for `Re_n(2)`), so the chain rule
  through `consis1`/`consis2` is exact.

The probe (§1.2) confirmed `qp` is a primed assign and `qp_GRID(1)` reaches codegen; the autodiff
`emitExpr` already has a future-var indexed-`call` branch (`emitExpr.m:68-92`) and SymPy has
`liftIndexed` (`Re_n(2)`), so after §5.2 normalization (`qp_GRID(1)`→`qp(1)`), both are expected to
work once `qp` is registered as a per-shock base in each backend's context. The plan verifies this
with a focused codegen-diff before the full build.

### 5.7 Parser — multiple inline reductions per statement

`parseModel`'s `parseStatement` currently finds all reduction openers (`findReductions`) and
hard-errors when `numel(redIdx) > 1`. Generalize: hoist **each** reduction (left-to-right) to its
own `GDSGE_<KIND>_<counter>` generated reduction statement, substitute its single-token name back
into the host token stream, then parse the rewritten host. The existing single-reduction branch is
the `n==1` case of this loop. Keep:

- the **nested-reduction** guard (a reduction inside another reduction's braces stays an error);
- the **whole-RHS** fast path (if the entire RHS is exactly one reduction, it becomes the statement
  directly — no hoist, preserving existing IR snapshots for HL1996/safe_assets/etc.).

Each hoisted reduction carries its own `transRef` (so `eq1`'s two `|shock_trans1` EXPECTs become two
independent reduction statements). Order: hoist in source order so generated names are stable and
the front-end IR snapshot is deterministic. After hoisting, the statements are ordinary reductions —
both backends already emit them (Phase 5/7a/8).

### 5.8 Parser — embedded named-interp call hoisting

`parsePostfix` (`parseExpr.m:75-78`) errors when a primed name is followed by `(` mid-expression,
and `parseModel`'s named-interp detection only fires when the call is the whole RHS. Add a
pre-pass in `parseStatement` (before reduction hoisting) that detects a named-interp call
`name'(args)` embedded in a larger RHS and hoists it to a generated primed target
`GDSGE_INTERP_<counter>' = name'(args)` (an `interpCall` statement, `interpRef=name`), substituting
the primed token `GDSGE_INTERP_<counter>'` back into the host. For `qp' = qFuture'(Wp')+d'` this
yields:

```
GDSGE_INTERP_1' = qFuture'(Wp')        % interpCall stmt, interpRef='qFuture'
qp'             = GDSGE_INTERP_1' + d'  % assign stmt, primed; refs the hoisted primed name
```

`analyzeModel` already treats interp-call targets as `primeable` (`analyzeModel.m:70-71`), so the
host's `GDSGE_INTERP_1'` resolves. Both backends already materialize primed assigns/interp targets
as per-shock arrays, so no backend change. Detection must distinguish the *whole-RHS* named call
(kept as today, so existing snapshots are stable) from the *embedded* case (hoisted). Multi-output
`GDSGE_INTERP_VEC'` is always whole-RHS and is unaffected.

## 6. Golden capture (first plan step)

Honoring "first verify the original toolbox runs": before any new-pipeline code, capture
Cao2011EZ from `base_package` in an isolated process.

- `tests/golden/capture_Cao2011EZ.m` — `addpath` old `base_package/gdsge/source` only, `cd` a
  scratch dir, run the model to convergence, save `IterRslt` and the (1-period) `SimuRslt`, restore
  path. **Pin `rng` before iter.** Cao2011EZ has KKT complementarity (`lambda*kp1`, `mu*a`) — almost
  certainly multi-root grid points, so the iter's randomized restart is order-sensitive (the
  documented safe_assets fragility — both backends affected). Capture and both end-to-end gates must
  pin the same `rng` seed before `iter_Cao2011EZ` so the converged root is reproducible; the old
  `test.m` does not seed (its simulate is a trivial 1 period), so we choose and record a seed.
- `tests/Cao2011EZ/golden/{IterRslt.mat, SimuRslt.mat}` committed (expected small: 100-pt grid, 2
  shocks). `tGoldenCao2011EZ` integrity test pins iteration count + metric + field checksums.
- The capture **also resolves the deferred questions**: whether the golden `IterRslt` carries a
  `var_tensor` field (→ §3 result-field decision) and the exact converged iteration/metric the
  gates assert.

## 7. Testing

Harness unchanged: headless `matlab.unittest` via `matlab -batch "cd('tests'); run_tests"`;
`junit.xml` + exit code authoritative; one-source-per-process path policy. New suites mirror the
existing per-model layout (`tests/Cao2011EZ/...`):

- **`tGoldenCao2011EZ`** (`tests/Cao2011EZ/`) — golden integrity (iteration, metric, field shapes).
- **`tFrontEndCao2011EZ`** (`tests/Cao2011EZ/parser/`) — `parseFrontEnd` → schema-valid IR with the
  enriched `variables.tensor`; IR JSON snapshot.
- **`tEndToEndCao2011EZ`** (`tests/Cao2011EZ/codegen/`) — autodiff: public API reproduces the golden
  (iteration + metric bit-match, fields within tolerance); second-run compile-skip; artifact
  presence; IR round-trip.
- **`tEndToEndCao2011EZSympy`** (`tests/Cao2011EZ/codegen/`) — `UseAutoDiff=0` reproduces the same
  golden; `assumeTrue(sympyAvailable)`. Where feasible, a three-way Jacobian cross-check
  (sympy == adept to ~1e-6 == finite-diff to ~1e-4) at the golden solution, reusing the Phase 8
  6th-output analytic-Jacobian harness.
- **Unit tests** for the new pieces: multiple inline reductions hoist to N `GDSGE_<KIND>_<n>`
  statements (and the whole-RHS single-reduction fast path is unchanged); an embedded named-interp
  call hoists to a `GDSGE_INTERP_<n>'` interpCall + a host referencing it (and the whole-RHS named
  call is unchanged); `_GRID`/bare-index parse to the same AST; parser routes tensor assignments to
  `tensorAssign` and not `setupStmts`; `emitBounds`/initial rewrite tensor names; the narrowed 7e
  gate (MATLAB-side tensor OK; body-tensor errors).
- **Regression:** full suite stays green (no existing IR golden regenerates — §5.3; the hoist fast
  paths keep HL1996/safe_assets/Mendoza/GLSW/CaoKS/Bianchi IR snapshots byte-identical).

## 8. Risks

1. **`IterRslt.var_tensor` result field.** The old golden may include a `var_tensor` struct.
   Mitigation: the §6 capture settles it; if present and the gate's field comparison needs it, emit
   it (small additive MATLAB-side assembly); if absent, stays deferred. Low risk.
2. **`qp_GRID(n)` not covered by the `Re_n(2)` path.** A primed-interp-derived array indexed in a
   sibling residual may need a small emitter extension in both backends (§5.6). Mitigation: probe
   immediately post-capture with a focused codegen-diff before broad work. Medium risk, well
   isolated.
3. **Tensor broadcast shape mismatch.** `etotal` (`1×2`, shock-only) vs `budget1` (`2×100`,
   shock×state) must both flatten correctly via `(:)'` under the `[shock × state]` `NPROB` ordering.
   Mitigation: the autodiff golden diff catches any transpose/broadcast error immediately;
   shock-only tensors broadcast against the state dim in `ndgrid`. Medium risk.
4. **SymPy parity.** Since tensors are MATLAB-side, the only SymPy-specific work is `qp(n)`
   (Risk 8.2). If that generalizes, SymPy parity is nearly free. Low–medium risk.
5. **IR schema change drift.** `variables.tensor` shape change touches schema/validator/round-trip/
   doc-gen. Mitigation: single-source descriptor + no-drift doc test; additive (empty encodes
   identically), so no golden regen. Low risk.
6. **Hoist fast-path regressions.** Generalizing the reduction hoist and adding the interp-call
   hoist must leave the *whole-RHS* single-reduction and whole-RHS named-interp cases byte-identical,
   or every existing IR snapshot regenerates. Mitigation: explicit fast-path branches + the full
   existing front-end suite (6 models) as the regression guard; assert no snapshot churn. Medium risk.
7. **Generated-name collisions.** Reduction hoist (`GDSGE_<KIND>_<n>`) and interp hoist
   (`GDSGE_INTERP_<n>`) share the statement `counter`. Mitigation: one monotone counter across both
   hoist kinds per model body; unit-test a statement that triggers both. Low risk.

## 9. Result-struct & API contract

Public entry points unchanged. `IterRslt`/`SimuRslt` shapes for all existing models untouched. A
`var_tensor` result field is introduced for Cao2011EZ **only if** the captured old golden carries
it and the gate requires it (§3, §8.1); otherwise unchanged. The `.gmod` surface syntax for
`var_tensor`, `inbound`, `initial`, and `_GRID` indexing is unchanged (we are matching the old
toolbox, not redesigning it).

## 10. Sequencing (for the plan)

0. **Done (2026-06-13):** model-body gap probe (§1.2) — confirmed the gap inventory.
1. Capture old-toolbox Cao2011EZ golden; pin `rng`+iteration/metric; resolve §8.1/§3 questions.
2. Parser hoisting: multiple inline reductions (§5.7) + embedded named-interp call (§5.8) — these
   unblock `parseModel(Cao2011EZ)` entirely (gaps 1–2). Guard the whole-RHS fast paths.
3. Parser: `_GRID` normalization (§5.2) + tensor-assignment capture (§5.1) + narrowed 7e gate (§5.5).
4. IR: enrich `variables.tensor` (schema/validator/round-trip/doc-gen) (§5.3) + front-end gate.
5. MATLAB codegen: tensor construction + bounds + initial + simulate recompute (§5.4).
6. Backend `qp(n)` indexing verify/extend (autodiff `emitExpr`, then SymPy `liftIndexed`) (§5.6) —
   focused codegen-diff first.
7. Autodiff end-to-end gate green → SymPy end-to-end gate green.
8. `PROGRESS.md` + branch merge.

CaoNie2016 follows as Phase 9b (separate brainstorm/spec).
