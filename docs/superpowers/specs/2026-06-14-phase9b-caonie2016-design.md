# Phase 9b Design — CaoNie2016 (conditional model regions + if/else equations)

- **Date:** 2026-06-14
- **Status:** Approved (design); pending plan
- **Parent spec:** `2026-06-11-refactor-gdsge-design.md` (§13 backward-compat surface, §15 Phase 7
  widening); continues the one-model-per-sub-phase precedent of 7a/7b/9a
- **Related:** `2026-06-13-phase9a-cao2011ez-var-tensor-design.md` (the prior corpus sub-phase;
  shares the model-body-probe-first discipline and the both-backends-at-parity goal)
- **Owner:** Wenlan Luo

---

## 1. Context

Phase 9a added **Cao2011EZ** to the backward-compat corpus, green end-to-end on both C++ Jacobian
backends. **CaoNie2016** is the second of the two example models the owner asked to add before
polish (Phase 10). Its source lives at `base_package/gdsge/tests/CaoNie2016/CaoNie2016.gmod`, with
a committed `IterRslt_CaoNie2016_1050.mat` (an iteration-1050 snapshot at `SaveFreq=10`).

CaoNie2016 is a two-agent (entrepreneur / household) production economy with a collateral
constraint and aggregate-productivity shocks. 1 state `X` (entrepreneur wealth share, 701 pts on
`[0,1]`), 3 shocks `A` (`[0.97 1 1.03]` with a 3×3 transition), 11 unknowns
(`p q c h a cp L mu` + vector `Xp[3]`), square 11×11 system per `(state × shock)` grid point.

Per the owner's two scoping decisions for this phase:

- **Both backends in one phase** — adept autodiff (default / backward-compat anchor) and SymPy
  (Phase 8) are brought to parity together, matching 9a. The new work is piecewise *control
  structure*, not new symbolic math, so the SymPy cost is structural emission rather than novel
  differentiation (§2).
- **Clean schema evolution** — the IR `model` section is restructured to carry regions and
  conditional equations as first-class shapes, bumping `irVersion` to **1.1.0** and regenerating
  every model's IR JSON snapshot **once** (a reviewed, mechanical migration), rather than bolting on
  additive-optional fields to keep existing goldens byte-identical (§5.1).

### 1.1 What CaoNie2016 exercises that the refactor lacks

`CaoNie2016.gmod` reuses constructs already supported — array policy `Xp[3]` with `Xp(n)` indexing
(HL1996 `w1n[8]`), per-shock future-var indexing `qp(n)` (9a), multi-output primed interp
`[cNext',cpNext',qp'] = GNDSGE_INTERP_VEC'(Xp')` and the named form `qp' = qFuture'(Xp')` (Phase
7a), `GNDSGE_EXPECT{}` reductions with reuse of a hoisted `Eqp` across equations (9a's
multi-reduction generalization), `model_init` (GLSW/9a), non-integer powers.

**Two constructs are new** (both confirmed empirically by a parser probe — see §1.2), plus one
wiring task and a set of verify items:

1. **Conditional model regions** — the gmod has **two** model bodies, guarded on the *state*:
   ```matlab
   model(X>0);  ... equations; ... end; end;
   model(X==0); ... equations; ... end; end;
   ```
   The new pipeline structurally assumes a **single** model body: `splitBlocks` matches a bare
   `model` opener only (it strips the trailing `;`, so `model(X>0);` is not recognized — the inner
   `equations;`/`end;` desync and it errors `gdsge:parser:unterminatedBlock`), the IR `model`
   section is one `{statements, equations}` struct, and `parseFrontEnd` calls `parseModel` exactly
   once. Each region is a separate square system selected per grid point by its state condition
   (`X==0` is the single boundary grid point `XMin=0`; `X>0` is every other).

2. **`if/else/end` conditional equations** — inside `equations;`, the consistency residual for a
   slot is selected at runtime on the *shock* value:
   ```matlab
   if A==AGood
     Xp(1);
   else
     consis1;
   end
   consis2;
   if A==ABad
     Xp(3);
   else
     consis3;
   end
   ```
   `parseModel` chokes on `if A==AGood` (`gdsge:parser:badExpr` — two names, no operator). The
   if/else/end lines carry no `;`, so `splitStatements` cannot group them; they need recognition as
   control lines that wrap one-or-more equation slots per branch.

3. **Wiring — `var_others shock_trans`** — `var_others` exports a setup-eval'd variable (here the
   3×3 transition matrix) into `IterRslt`. `parseVarDecls` already collects it into
   `decls.otherNames` → `ir.variables.others` (probe-confirmed), and the schema already has
   `variables.others`. The remaining work is the **MATLAB codegen** packing those names into
   `IterRslt` (the old toolbox's `RSLT_VAR_OTHERS`). No C++ work.

4. **Verify (probe-confirmed parse-clean; confirm they flow through analyze → both backends):**
   array policy `Xp[3]`; const-index `qp(1)`/`Xp(1)` in the body; `var_output Xp` (an array policy
   in the output set); `X' = Xp'` transition (array next-state selected by realized shock — HL1996
   `w1' = w1n'` precedent); empty `var_tensor;` (probe: `tensorNames = {}`, a no-op, no Phase 7e
   error).

### 1.2 Probe findings (2026-06-14, empirical)

Before this spec's gap inventory — honoring 9a's lesson that static reading mis-scoped the gap
list — a one-source-per-process probe (`scratch/probe_caonie2016.m`) ran the new parser's stages on
CaoNie2016, carving the two regions manually (since `splitBlocks` cannot) and testing every body
statement and equation in isolation. Confirmed results (the source of the inventory above):

- `splitBlocks` → `gdsge:parser:unterminatedBlock` at line 170 (gap 1, structural).
- `parseModel` on each region body → `gdsge:parser:badExpr` on the `if A==AGood` equation line
  (gap 2).
- **All 20 body statements of region 1 parse OK** — including `[cNext',cpNext',qp'] =
  GDSGE_INTERP_VEC'(Xp')`, `Eqp = GDSGE_EXPECT{qp'}`, the four equations carrying an inline
  `GDSGE_EXPECT` plus a reference to `Eqp`, and `consis1 = (qp(1)*h+b)/(qp(1)*H)/Xp(1) - 1`. No new
  expression-level work.
- **All plain equations parse OK**; only the `if A==AGood` / `if A==ABad` lines fail (gap 2).
- `parseVarDecls` on the declaration region: `otherNames = {shock_trans}` (gap 3 captured),
  `tensorNames = {}` (empty `var_tensor` is a no-op — no Phase 7e error), `policy` includes the
  array `Xp`, `output` includes `Xp`.
- Full `parseFrontEnd` → `gdsge:parser:unterminatedBlock` (the gap-1 structural failure fires
  first, as expected).

This is the authoritative gap list; the plan is built against it.

## 2. Key architectural insight — the new work is control structure, not new math

Both new constructs are **piecewise control flow** over machinery that already exists:

1. **A region is just a model body with a guard.** `parseModel` already turns a body into
   `{statements, equations}`; we call it per region and tag each with its (opaque-text) state
   condition. Each backend already emits one body; it now emits each region's body inside a guard
   `if (<state cond>) { ... }`. adept records through the taken branch; SymPy differentiates each
   region independently with its existing `diff_body` machinery (called per region, emitted
   guarded). **No new symbolic capability** — the per-region Jacobian is the same analytic Jacobian
   Phase 8 already produces, emitted twice under a runtime branch.

2. **An if/else equation is a per-slot conditional residual.** Each branch contributes the same
   number of equation slots (validated equal across branches); the residual filling slot *k* is
   chosen at runtime on the shock condition. Autodiff emits `if (<shock cond>) { GDSGE_EQ[k]=…; }
   else { GDSGE_EQ[k]=…; }`. SymPy differentiates *each branch's* expression separately and emits
   the value + gradient rows under the same guard. The condition is opaque text lowered to C++ by a
   small translator (state/shock C++ locals; params are `#define`d constants), mirroring the old
   `gdsge_parser`'s `replace_grid_variable` + `MODEL_CONDITION`.

Consequently the SymPy-specific work is **emitting per-region / per-branch value+gradient under a
runtime guard**, not deriving anything new — unlike Phase 8's reduction-fusion chain rule or Phase
8b's deferred ASG-interp derivatives. This is what makes "both backends in one phase" tractable
here, and it is cross-checked against adept via the Phase-8 6th-output analytic Jacobian.

## 3. Scope decisions

**In scope:**

1. **Conditional model regions** (parser + IR + both backends) — `splitBlocks` recognizes
   `model(<cond>);`; the IR `model` section becomes a `regions` list; `parseFrontEnd` iterates
   regions; `analyzeModel` checks each region as an independent square system; both backends emit
   guarded region bodies. State-condition lowering to C++.
2. **`if/else/end` conditional equations** (parser + IR + both backends) — the equations parser
   recognizes `if/elseif/else/end` control lines and groups them into a `conditional` equation node
   (cases, each with one-or-more plain equations, equal count across cases); both backends emit the
   guarded per-slot residual. Shock-condition lowering to C++.
3. **`var_others` → `IterRslt`** (MATLAB codegen) — pack `ir.variables.others` names into the
   result struct (old `RSLT_VAR_OTHERS`).
4. **IR schema evolution to 1.1.0** — `model.regions`, tagged equation list (`plain | conditional`),
   `modelInit.equations` likewise tagged; one reviewed regeneration of every model's IR JSON golden
   (§5.1).
5. **Verify items (§1.1.4)** flow through analyze + both backends (array policy in output,
   const-index, array transition, empty tensor no-op).
6. **Goldens + gates** for CaoNie2016 on both backends; `PROGRESS.md` bookkeeping.

**Deferred (recorded in `PROGRESS.md`):**

- **Multi-equation-per-branch conditionals beyond what the corpus needs** — the representation
  supports N plain equations per branch (equal across branches), but **nested** `if` inside a branch
  is rejected with a clear error (`gdsge:parser:nestedConditionalEquation`); no corpus model nests.
- **`if/else` in the model *body*** (not the equations block) — the old `gdsge_parser` supports it
  (lines 2242–2281); CaoNie2016 uses if/else only in `equations`. A body-level `if` raises a clear
  deferred error. (Forward-compatible: the same condition-lowering translator would serve it.)
- **`var_others` for non-setup expressions** — only setup-eval'd names (the old behavior) are
  exported; arbitrary derived expressions are out of scope (none needed).

**Out of scope:** performance work beyond parity (Phase 10); any change to the `.gmod` surface
syntax for `model(cond)`, `if/else`, or `var_others` (we match the old toolbox, not redesign it);
ASG for CaoNie2016 (it is a spline model — `USE_SPLINE=1, USE_ASG=0`).

## 4. Goal & success criteria

Phase 9b is done when:

1. **Old-toolbox golden captured first** (honoring "verify the original toolbox runs"): CaoNie2016
   regenerated and converged in `base_package` in an isolated `matlab -batch` process under the
   one-source-per-process path policy; `tests/CaoNie2016/golden/` holds `IterRslt.mat` + a reduced
   seeded `SimuRslt.mat` + the capture script. (Primary path; fallback in §6.)
2. **Autodiff end-to-end gate** green: the public API (`gdsge_codegen` → `iter_CaoNie2016` →
   `simulate_CaoNie2016`, or the `gdsge`/`gndsge` orchestrator) reproduces the golden — iteration
   count and metric bit-match; `IterRslt`/`SimuRslt` fields within tolerance; `IterRslt.var_others`
   (`shock_trans`) present and equal.
3. **SymPy end-to-end gate** green: same model under `UseAutoDiff=0` reproduces the same golden;
   guarded `assumeTrue(sympyAvailable)` so the default Python-free loop stays green. Where feasible,
   a per-region Jacobian cross-check (sympy == adept to ~1e-6 == finite-diff to ~1e-4) at the golden
   solution, reusing the Phase 8 6th-output analytic-Jacobian harness.
4. **Front-end gate** green: `parseFrontEnd(CaoNie2016.gmod)` yields a schema-valid IR with two
   `model.regions` and a `conditional` equation node in each; IR JSON snapshot committed.
5. **The IR migration is clean**: every existing model's regenerated IR golden parses `isequalIR`
   its pre-migration IR after the mechanical wrap (regions[0], plain-tagged equations); round-trip
   and no-drift doc tests updated; `docs/ir-schema.md` regenerated.
6. All existing gates (HL1996, safe_assets, Mendoza2010, GLSW, CaoKS2016, Bianchi2011, Cao2011EZ,
   macro / orchestration / sympy suites) stay green after the migration.

## 5. The change (detail)

### 5.1 IR schema — regions + tagged equations (irVersion 1.1.0)

`gdsge.ir.schema` is the single source of truth (drives validator, round-trip, doc-gen). Changes:

- **`model`** becomes `{ regions: fList(regionItem) }`, where
  `regionItem = { condition: fText, statements: fStmtList, equations: fList(eqEntry) }`.
  `condition` is opaque MATLAB-text on the states; `''` means unconditional (always). Existing
  single-body models migrate to a 1-element `regions` list with `condition = ''`.
- **Equation entries become a tagged list** over a new `eqkinds` registry (mirroring the `stmts`
  registry mechanism):
  - `plain`   → `{ kind:'plain', expr: fNode, primed: fScalar }` (today's `eqItem` + a `kind` tag).
  - `conditional` → `{ kind:'conditional', cases: fList(caseItem) }`, where
    `caseItem = { cond: fText, equations: fList(plainEqItem) }`; the final case with `cond=''` is the
    `else`. `equations` inside a case are `plain` only (nested conditionals rejected at parse — §3).
- **`modelInit.equations`** uses the same tagged list (always one unconditional region; `modelInit`
  keeps `statements`/`equations` directly — model_init is never region-split in the corpus).
- **`irVersion` → 1.1.0.** The validator, type-aware round-trip (`canonicalize`/`encode`/`decode`/
  `isequalIR`), and doc generator update together. The no-drift doc test regenerates
  `docs/ir-schema.md`.

**Migration of existing goldens (one reviewed pass):** every committed model IR JSON snapshot is
regenerated by re-running `parseFrontEnd` (which now emits the regions/tagged shape) and the change
is reviewed as a mechanical wrap — `model.statements/equations` → `model.regions[0]` with
`condition=''`; each equation gains `kind:'plain'`. A differential test asserts each model's
regenerated IR `isequalIR` a programmatic wrap of its old IR (no semantic drift, only shape).

### 5.2 Parser — `splitBlocks` multi-region

`splitBlocks` learns the `model(<cond>);` opener (regex `^model(\(.*\))?;$`, capturing the
condition), and emits an ordered `out.modelRegions = [{condition, body}, …]` (depth-tracked exactly
as today so nested `equations;…end;` stays inside the region; a bare `model;` → one region,
`condition=''`). The singleton blocks (`simulate`, `model_init`, hooks) are unchanged. `model_init`
remains a single block (no conditional model_init in the corpus; a `model_init(<cond>)` raises a
clear deferred error).

### 5.3 Parser — `parseModel` per region + if/else equations

- `parseFrontEnd` iterates `out.modelRegions`, calling `parseModel` on each body and assembling
  `model.regions = [{condition, statements, equations}, …]`.
- `parseModel`'s **body** handling is unchanged (probe-confirmed: all body statements parse).
- `parseModel`'s **equations** handling learns `if/elseif/else/end` control lines (no `;`):
  scanning the equations region line-wise, an `if <cond>` opens a conditional block, `elseif
  <cond>` / `else` start new cases, `end` closes it; `expr;` lines are plain equations belonging to
  the current case (or the top level). The result is an ordered list of `plain | conditional`
  entries (mirrors old `gdsge_parser` lines 2635–2730). Each `conditional` block's cases must hold
  the same number of plain equations (validated). Conditions are kept as opaque text; nested `if`
  inside a case errors (`gdsge:parser:nestedConditionalEquation`).

### 5.4 Parser — `analyzeModel` per region

- **Name resolution** runs per region (each region's `statements` define names consumed by its own
  `equations`); cross-region names are independent.
- **Square check** runs per region: count equation slots, where a `conditional` entry counts as its
  per-case plain-equation count; the total must equal the number of unknowns (11). A region whose
  count differs errors clearly, naming the region condition.
- **Condition validation**: a region `condition` may reference only state names (+ params); an
  equation-`cond` may reference only shock names (+ params). Unknown identifiers error.
- The verify items (§1.1.4) are exercised here: an array policy in `var_output`, const-index of a
  primed/array var, and the `X' = Xp'` array transition resolve as today (HL1996 precedent); the
  empty tensor list is a no-op.

### 5.5 Both backends — guarded regions + guarded equation slots

A shared **condition-lowering** helper translates an opaque condition to a C++ boolean expression,
mapping state/shock identifiers to their per-grid-point C++ locals and leaving params (which are
`#define`d constants, e.g. `AGood`) as-is — the new-pipeline analogue of the old
`replace_grid_variable`. Used by both backends for both region and equation conditions.

- **Autodiff (C++):** the model function emits each region's body inside `if (<state cond>) { … }
  else if (…) { … }` (an unconditional region is the trailing `else`/`{1}` arm). Within a region's
  equation assembly, a `conditional` entry emits `if (<shock cond>) { GDSGE_EQ[k]=…; } else if (…)
  { … }`. adept records the taken branch; the recorded Jacobian is exact for the active region/case.
- **SymPy (C++):** the registry/emitters in `gdsge.codegen.cxx.+sympymodel` run **per region** —
  each region's value+gradient C++ is generated by the existing `diff_body` machinery and wrapped in
  the region guard. A `conditional` equation differentiates each branch's residual expression
  separately and emits the value + gradient rows under the shock guard, so the analytic Jacobian
  matches whichever branch is active. Shared-CSE (`helper_i`) stays per emitted block.

This is the only non-mechanical backend work; it is verified by the §7 cross-check (sympy == adept).

### 5.6 MATLAB codegen — `var_others`

`generateMatlab` packs each name in `ir.variables.others` into the assembled `IterRslt` (old
`RSLT_VAR_OTHERS`: the setup-eval'd value is available in the iter workspace and dot-assigned to the
result struct). Absent/empty `others` (every other corpus model) emits nothing — no change to their
generated code or result shape.

## 6. Golden capture (first plan step)

Honoring "first verify the original toolbox runs": before any new-pipeline code, capture CaoNie2016
from `base_package` in an isolated process.

- **Primary — regenerate via the vendored old toolbox.** Its `gdsge_parser` supports
  `model(<cond>)` (`extract_model_segment`), `if/else` equations (lines 2635–2730), and `var_others`
  (lines 306, 1560); and CaoNie2016 has **≤1 literal `GNDSGE_EXPECT` per line** (`Eqp` is a hoisted
  variable, not an inline reduction), so the 9a two-`EXPECT`-per-line blocker does **not** apply —
  the old parser is expected to regenerate, compile (MSVC, proven Phase 0), and solve the model
  cleanly, exactly like 7a/7b. `tests/golden/capture_CaoNie2016.m` `addpath`s old
  `base_package/gdsge/source` only, `cd`s a scratch dir, runs to convergence, saves `IterRslt` and a
  **reduced seeded** `SimuRslt` (the gmod's `100000 × 100` is impractical for a gate — capture
  `6 × 1000` per the 7a precedent), restores the path.
- **Pin `rng`.** CaoNie2016 has a collateral-constraint KKT structure (`mu*a`, the `X==0` corner)
  — almost certainly multi-root grid points, so the iter's randomized restart is order-sensitive
  (the documented safe_assets fragility; both backends affected). Capture and both end-to-end gates
  pin the same seed before `iter_CaoNie2016` and before `simulate_CaoNie2016`; the chosen seed is
  recorded.
- **Cross-check / fallback — the committed `IterRslt_CaoNie2016_1050.mat`.** Verify it is converged
  (metric `< TolEq`) and use it to cross-check the regenerated `IterRslt`; if old-toolbox
  regeneration proves infeasible or too slow (Risk 8.2), fall back to comparing the new pipeline's
  `IterRslt` against this committed artifact (iteration count unknown in that case → the gate
  asserts field-level tolerance, not a bit-matched iteration).
- `tests/CaoNie2016/golden/{IterRslt.mat, SimuRslt.mat}` committed; `tGoldenCaoNie2016` integrity
  test pins iteration count + metric + field checksums. The capture also settles the exact converged
  iteration/metric the gates assert.

## 7. Testing

Harness unchanged: headless `matlab.unittest` via `matlab -batch "cd('tests'); run_tests"`;
`junit.xml` + exit code authoritative; one-source-per-process path policy. New suites mirror the
existing per-model layout (`tests/CaoNie2016/...`):

- **`tGoldenCaoNie2016`** (`tests/CaoNie2016/`) — golden integrity (iteration, metric, field shapes;
  `var_others.shock_trans` present).
- **`tFrontEndCaoNie2016`** (`tests/CaoNie2016/parser/`) — `parseFrontEnd` → schema-valid IR with two
  `model.regions` (conditions `X>0`, `X==0`) and a `conditional` equation entry in each; IR JSON
  snapshot.
- **`tEndToEndCaoNie2016`** (`tests/CaoNie2016/codegen/`) — autodiff: public API reproduces the
  golden (iteration + metric bit-match where regenerated, fields within tolerance); second-run
  compile-skip; artifact presence; IR round-trip.
- **`tEndToEndCaoNie2016Sympy`** (`tests/CaoNie2016/codegen/`) — `UseAutoDiff=0` reproduces the same
  golden; `assumeTrue(sympyAvailable)`. Where feasible, a per-region three-way Jacobian cross-check
  (sympy == adept to ~1e-6 == finite-diff to ~1e-4) at the golden solution, reusing the Phase 8
  6th-output analytic-Jacobian harness, exercising both the `X>0` and `X==0` regions and an active
  if/else branch.
- **Unit / differential tests** for the new pieces:
  - `splitBlocks` returns two `modelRegions` with the right conditions; a bare `model;` → one region
    with `condition=''`.
  - `parseModel` equations parser groups `if/else/end` into a `conditional` entry with the right
    cases; unequal per-case counts and nested `if` error clearly.
  - `analyzeModel` square check counts a `conditional` block as one slot-group; a wrong count errors
    naming the region.
  - condition lowering: `X>0`/`X==0` and `A==AGood`/`A==ABad` translate to the expected C++ booleans.
  - `var_others` appears in a generated `IterRslt` (a small `generateMatlab` snapshot).
  - **IR migration differential**: each existing model's regenerated IR `isequalIR` a programmatic
    wrap of its old IR (regions[0] + plain-tagged equations) — proves the schema bump is shape-only.
- **Regression:** full suite stays green after the one-time IR golden regeneration (every model
  snapshot regenerated and reviewed; round-trip + no-drift doc tests updated).

## 8. Risks

1. **IR migration drift.** The schema bump regenerates every model's IR JSON golden and updates
   round-trip/validator/doc expectations. Mitigation: single-source descriptor; the §7 migration
   differential proves the change is shape-only (regions[0] + plain tag), so any semantic drift fails
   a gate; review the regenerated snapshots once. Medium risk (broad but mechanical).
2. **Old-toolbox regeneration feasibility / runtime.** The 701×3 grid over ~1050 iterations may be
   slow to capture, and the old parser must actually emit/compile CaoNie2016. Mitigation: the old
   parser provably supports every construct (§6) and ≤1 EXPECT/line clears the 9a blocker; capture
   once and commit; the committed `IterRslt_..._1050.mat` is the fallback reference (§6). Reduce the
   simulate to `6×1000` for the gate. Medium risk, mitigated by the fallback.
3. **SymPy per-region / per-branch gradient emission.** The only non-mechanical backend work:
   emitting value+gradient under runtime guards for two regions and the if/else branches.
   Mitigation: the per-body `diff_body` machinery is reused unchanged (it differentiates one body);
   the guards are structural; the §7 adept cross-check (per region + an active branch) catches any
   mis-emitted row. Medium risk, well isolated.
4. **Condition lowering correctness.** `X==0` and `A==AGood` rely on exact double equality; the grid
   places `X=0` exactly at `XMin` and the shock values are exactly `ABad/AGood`, so equality holds
   (matching the old toolbox, which does the same). Mitigation: the autodiff golden diff catches any
   mis-selected region/branch immediately (a wrong region at the boundary point diverges visibly).
   Low–medium risk.
5. **Square-system count with conditional blocks.** `analyzeModel` must count a conditional block as
   one slot-group (per-case count), not as its flattened lines, or the square check mis-fires.
   Mitigation: explicit per-region count using the `conditional` entry's case arity; unit test on
   both regions (11 == 11). Low risk.
6. **`var_others` result placement.** The setup-eval'd `shock_trans` must be available in the iter
   workspace at result-assembly time and dot-assigned. Mitigation: it is a setup name (the
   `setupNames` / setup-replay machinery from 7a/7b already evaluates it); small additive MATLAB-side
   assembly verified by a `generateMatlab` snapshot + the golden's `var_others` field. Low risk.
7. **`model_init` / array-output interactions.** `var_output Xp` (array policy) and `model_init`'s
   own equations (now tagged `plain`) must round-trip and codegen as before. Mitigation: covered by
   the front-end + autodiff gates and the migration differential. Low risk.

## 9. Result-struct & API contract

Public entry points unchanged. `IterRslt`/`SimuRslt` shapes for all existing models untouched.
CaoNie2016's `IterRslt` gains a `var_others` member (`shock_trans`) to match the old golden — this
is the documented `var_others` behavior, not a new shape for existing models. The `.gmod` surface
syntax for `model(cond)`, `if/else` equations, and `var_others` is unchanged (we match the old
toolbox). The IR JSON **format** changes (irVersion 1.1.0, regions + tagged equations) — an internal
contract, regenerated once and reviewed; no public API consumes the IR JSON directly.

## 10. Sequencing (for the plan)

0. **Done (2026-06-14):** model-body gap probe (§1.2) — confirmed the gap inventory.
1. Capture old-toolbox CaoNie2016 golden (regenerate → compile → solve → reduced seeded simulate);
   pin `rng`+iteration/metric; verify the committed `..._1050.mat` as cross-check/fallback (§6).
2. IR schema evolution to 1.1.0 (§5.1): `model.regions`, tagged `plain|conditional` equations,
   validator + round-trip + doc-gen; the migration differential + one reviewed golden regen.
3. Parser: `splitBlocks` multi-region (§5.2) + `parseModel` per-region & if/else equations (§5.3) +
   `analyzeModel` per-region square check & condition validation (§5.4) → front-end gate green.
4. Both backends: shared condition lowering + guarded region bodies + guarded equation slots
   (§5.5) — autodiff first, then SymPy (per-region/per-branch emission), with the per-region
   cross-check as the SymPy guard.
5. MATLAB codegen: `var_others` → `IterRslt` (§5.6).
6. Autodiff end-to-end gate green → SymPy end-to-end gate green.
7. `PROGRESS.md` + branch merge.

This is the last corpus model before Phase 10 (polish).
