# Phase 7c Design — Legacy Orchestration & ASG-interp Simulate

- **Date:** 2026-06-13
- **Status:** Approved (design); pending implementation
- **Parent spec:** `2026-06-11-refactor-gdsge-design.md` (§13, §15 Phase 7); the original
  Phase 7c grab-bag is split here (decision below)
- **Owner:** Wenlan Luo

---

## 1. Context

Phases 7a (spline widening) and 7b (ASG widening) brought all six corpus models green
end-to-end through the public API. What remains of the backward-compat surface is a set of
**legacy items no corpus model exercises** — `PROGRESS.md` calls this the "golden-less
legacy surface" (Phase 7c): the macro engine, `gdsge.m`, `GenCodeSegment`, the legacy
5-output `gdsge_codegen` signature, `var_tensor`, and the ASG + `SIMU_INTERP` simulate
variant.

These six items are genuinely independent and vary widely in size. Following the 7a/7b
precedent of focused sub-phases (own spec → plan → branch → merge), the owner approved
(2026-06-13) **splitting Phase 7c into three sub-phases**:

- **7c (this spec) — legacy orchestration & ASG-interp simulate:** `gdsge.m`, the 5-output
  `gdsge_codegen` signature, `GenCodeSegment`, and ASG + `SIMU_INTERP`. The smaller items;
  the first three reuse existing corpus goldens, the fourth needs one small new capture.
- **7d — macro engine** (`#define/#for/#foreach/#if/#mat{}/#strcat_comma/include/cinclude`).
- **7e — `var_tensor`** (new variable kind; spline-path only).

`PROGRESS.md` is updated to reflect this split.

**Old-toolbox reference points** (all under `base_package/gdsge/source/`):
- `gdsge.m` — the orchestrator being mirrored (codegen → cache-gated solve → simulate →
  `eq = v2struct(model,IterRslt,SmltRslt)`).
- `gdsge_codegen.m` — the 5-output signature `[model,iterCode,cppCache,cppCode,codeSegment]`
  (line 25); `GenCodeSegment` writes five intermediate `.m` files (lines 68–90).
- `gdsge_parser.m` — `codeSegment = v2struct(setParamsCode,iterInitCode,prepareSpaceCode,
  constructSplineCode,solveAndAssignCode)` (line 1868); simulate-template selection
  (lines 1684–1742), where ASG + `SIMU_INTERP` reuses the **generic**
  `simulate_interp_template.m` with construction
  `GDSGE_ASG_INTERP = asg.construct_from_struct(IterRslt.asg_output_struct)` and per-period
  `GDSGE_ASG_INTERP.eval_vec(shock(:,t)', [states])` (lines 1686–1688).

**Groundwork already in the new package:**
- `src/gdsge_codegen.m` — current single-output shim over `gdsge.codegen.codegen`; its
  docstring already names the deferred outputs.
- `src/+gdsge/+codegen/codegen.m` — raises a deliberate "deferred to Phase 7c" error for
  `options.GenCodeSegment` (≈lines 42–46); owns the `mex_<model>.cache` gating (Phase 6).
- `src/+gdsge/+codegen/+mat/emitSimulateInterp.m` — the spline `SIMU_INTERP` emitter (7a),
  to be extended with an ASG branch.
- `src/+gdsge/+parser/resolveOptions.m:33–35` — raises the 7b "SIMU_INTERP with USE_ASG is
  not supported yet (deferred to Phase 7c)" error, to be removed.
- The vendored public `asg` class and ASG runtime helpers (7b).

## 2. Scope decisions

**In scope (four items):**

1. `gdsge_codegen` legacy **5-output signature** `[model,iterCode,cppCache,cppCode,
   codeSegment]`.
2. **`GenCodeSegment`** option → informative error (the new architecture has no equivalent
   segment decomposition).
3. **`gdsge.m`** top-level orchestrator (flat shim, no `v2struct`).
4. **ASG + `SIMU_INTERP`** simulate variant.

**Deferred to later sub-phases:** macro engine (7d); `var_tensor` (7e); SymPy backend
(Phase 8); performance work beyond parity (Phase 9).

## 3. Goal & success criteria

Phase 7c is done when:

1. `gdsge_codegen('<model>')` returns the **5-output** form, and existing single-output
   callers (`model = gdsge_codegen(...)`) are unaffected. `iterCode`/`cppCache`/`cppCode`/
   `codeSegment` carry the real generated artifacts (§5).
2. `options.GenCodeSegment==1` raises a clear, identified error (§5).
3. `gdsge('<model>')` drives codegen → cache-gated solve → simulate and returns
   `eq` with fields `model`, `IterRslt`, `SmltRslt`, reproducing the existing model's
   golden `IterRslt`/`SmltRslt`; a second call with unchanged code **skips the solve**
   (IterRslt-cache hit). No `v2struct`.
4. A gmod combining `USE_ASG` with `SIMU_INTERP` parses (no error), and its generated
   `simulate_<model>.m` evaluates `asg_output_struct` per period (§6); the result matches a
   freshly captured old-toolbox golden (shock path bit-exact, fields within tolerance).
5. All HL1996 / 7a / 7b gates stay green throughout.

## 4. Approach

The established TDD + differential discipline, sequenced low-risk first:

1. **5-output signature → GenCodeSegment error** — mechanical surface on `gdsge_codegen` /
   `gdsge.codegen.codegen`; structural tests only.
2. **`gdsge.m`** — reuse an existing fast model's golden (HL1996); no new capture.
3. **ASG + `SIMU_INTERP`** — the one behavior-bearing item: capture a single small
   old-toolbox golden (CaoKS2016 with `SIMU_INTERP` flipped on; it already declares
   `var_output Kp Xp kp1 kp2`), then extend `emitSimulateInterp` and differential-test.

Own branch `phase7c-legacy-orchestration`, merged when green.

## 5. Orchestration components (items 1–3)

### 5.1 `gdsge_codegen` 5-output signature

`src/gdsge_codegen.m` widens to `[model,iterCode,cppCache,cppCode,codeSegment]`. This is a
**backward-compatible** change: MATLAB lets callers request fewer outputs, so
`model = gdsge_codegen('HL1996')` keeps working unchanged.

- `model` — IR struct (unchanged).
- `iterCode` — text written to `iter_<model>.m`.
- `cppCache` — prior `mex_<model>.cache` content, read **before** this run overwrites it
  (`''` if none). This lets `gdsge.m` detect whether the C++ changed.
- `cppCode` — current `mex_<model>.cpp` content.
- `codeSegment` — struct carrying the **real generated strings** (e.g. `iterCode`,
  `simulateCode`, `cppCode`, `compileCode`). Deliberately **not** the old five fragments
  (`setParamsCode`/`iterInitCode`/`prepareSpaceCode`/`constructSplineCode`/
  `solveAndAssignCode`), which have no counterpart in the thin-file architecture.

`gdsge.codegen.codegen` already owns the cache and writes the files; it surfaces these
strings (returns them and the pre-overwrite cache content) so the flat shim stays thin and
does no file I/O of its own beyond what `codegen` already does.

### 5.2 `GenCodeSegment` → informative error

When `options.GenCodeSegment==1`, raise an error
(identifier `gdsge:codegen:genCodeSegmentUnsupported`, or sibling under
`gdsge:codegen:*`) explaining that the refactored thin-file architecture (explicit named
locals over `gdsge.runtime` helpers, no `v2struct`) has **no equivalent segment
decomposition**, and pointing the user at the generated `iter_<model>.m` and
`mex_<model>.cpp`. This replaces the current "deferred to Phase 7c" stub in `codegen.m`.

### 5.3 `gdsge.m` orchestrator

New flat shim `src/gdsge.m`, mirroring the old driver's control flow:

1. `[model,iterCode,cppCache,cppCode] = gdsge_codegen(modelName);`
2. **IterRslt cache gate** (old `gdsge.m` semantics): read `iter_<model>.cache`; if it
   equals `iterCode` **and** `cppCache` equals `cppCode` **and** `IterRslt_<model>.mat`
   exists, `load` it; otherwise solve (`iter_<model>`), `save` the mat, and write the iter
   cache.
3. `SmltRslt = simulate_<model>(IterRslt);`
4. Return `eq` assembled **explicitly** — `eq.model`, `eq.IterRslt`, `eq.SmltRslt` — the
   same field shape the old `v2struct(model,IterRslt,SmltRslt)` produced, with no
   `v2struct`.

Console output stays in the spirit of the old driver (parse / solve-or-load / simulate
banners), consistent with the structured-printing direction of the refactor.

## 6. ASG + `SIMU_INTERP` (item 4)

Extend `src/+gdsge/+codegen/+mat/emitSimulateInterp.m` (currently spline-only) with an ASG
branch, mirroring the old generic `simulate_interp_template.m` specialization
(`gdsge_parser.m:1684–1736`):

- **Construct:** `GDSGE_ASG_INTERP = asg.construct_from_struct(IterRslt.asg_output_struct);`
- **Per-period eval:** `GDSGE_INTERP_RESULTS = GDSGE_ASG_INTERP.eval_vec(shock(:,GDSGE_t)',
  [<states>]);`

in place of the spline path's `GDSGE_PP = IterRslt.output_interp;` + `myppual_mex(...)`. The
output-index assignment (`output_var_index`) and the simulate state/post-assign machinery
are shared with the existing spline interp path.

- **Parser:** remove the 7b error at `resolveOptions.m:33–35`; `USE_ASG` + `SIMU_INTERP`
  now resolves cleanly (with `USE_SPLINE`/`USE_PCHIP` exclusivity and `SIMU_RESOLVE` ⊕
  `SIMU_INTERP` still enforced).
- **Runtime:** reuse the ASG reconstruction helpers vendored in 7b.
- **No C++/MEX change:** interp-simulate does **not** re-solve the model, so there is no
  per-period model evaluation and no new MEX contract — distinct from the ASG
  `SIMU_RESOLVE` path delivered in 7b.

## 7. Testing

Harness unchanged: headless `matlab.unittest` via `matlab -batch "cd('tests');
run_tests"`; `junit.xml` + exit code authoritative; one-source-per-process path policy.

- **5-output signature** — structural unit test: all five outputs present and correct
  (`iterCode` equals `iter_<model>.m` content; `cppCache`/`cppCode` semantics; `codeSegment`
  carries the real strings); the single-output call `model = gdsge_codegen(...)` still
  returns the IR.
- **GenCodeSegment** — asserts the informative error (identifier + message) when the option
  is set.
- **`gdsge.m`** — Slow end-to-end on HL1996: `eq.IterRslt`/`eq.SmltRslt` match the existing
  golden; a second call hits the IterRslt cache and **skips the solve** (assert via
  iteration/timing or a solve sentinel); `eq` field shape verified. **Reuses existing
  goldens — no new capture.**
- **ASG + `SIMU_INTERP`** —
  - *Front-end gate:* the fixture gmod parses (no error), IR carries `interpMethod='asg'`
    with the interp-simulate selection; JSON round-trips; IR snapshot committed once the
    end-to-end gate is green.
  - *End-to-end gate (Slow):* `gdsge_codegen` + `iter` + the new interp `simulate` match the
    freshly captured golden — seeded shock path bit-exact, `SimuRslt` fields within the
    comparison utility's tolerance.
- All HL1996 / 7a / 7b gates stay green throughout.

## 8. Golden capture protocol (one new capture)

One old-toolbox capture, in its **own `matlab -batch` process** with controlled cwd and
only `base_package/gdsge/source` on the path (per the standing policy):

- **CaoKS2016 + `SIMU_INTERP`:** the CaoKS2016 gmod copied into the new test fixture with
  `SIMU_INTERP=1; SIMU_RESOLVE=0;` (the model already declares `var_output`). Iterate to the
  same settings the 7b CaoKS gate uses (so the converged `IterRslt`/`asg_output_struct` is
  comparable), then a **seeded reduced simulate** (6 samples × 1,000 periods, matching the
  7a/7b sizes) through the interp path. The committed capture script records exact settings
  and wall-clock timing.
- Golden lands under `tests/CaoKS2016_simu_interp/golden/` (a sibling fixture, to keep the
  existing `tests/CaoKS2016/` golden untouched) with a golden-integrity test. Size checked
  before commit (expected small, like the other ASG goldens — far under Mendoza's 99.8 MB
  outlier). If the capture proves prohibitively long, pause and revisit before building the
  gate on it.

## 9. Result-struct & API contract (frozen)

- `gdsge(modelName)` returns `eq` with fields `model`, `IterRslt`, `SmltRslt` — identical
  shape to the old `v2struct`-built struct.
- `gdsge_codegen` outputs 1..5 match the old positional contract; `IterRslt`/`SimuRslt`
  field shapes are unchanged from the spline/ASG paths.
- ASG `IterRslt` continues to carry `asg_output_struct`/`output_var_index` (7b); the interp
  simulate consumes them.

## 10. Risks

- **IterRslt-cache fidelity in `gdsge.m`.** The new `gdsge.codegen.codegen` already does
  MEX-cache gating; `gdsge.m` adds a second (IterRslt) layer keyed on iter+cpp code. Risk:
  the two layers interacting unexpectedly (e.g. a recompile that doesn't change `cppCode`).
  Mitigation: mirror the old key logic exactly and assert the skip-solve path in the gate.
- **`codeSegment` shape expectations.** No corpus caller inspects it, but the field names
  shouldn't masquerade as the old five fragments. Mitigation: name them after the real
  generated artifacts and document the deviation in `PROGRESS.md`.
- **ASG interp-simulate fixture validity.** Flipping CaoKS2016 to `SIMU_INTERP` is a
  synthetic configuration the model's authors didn't ship. Mitigation: it's a pure
  evaluate-the-output-interpolant path (no re-solve), differential-tested against the old
  toolbox running the identical flipped gmod, so both sides exercise the same code path.

## 11. Out of scope

Macro engine (7d); `var_tensor` (7e); SymPy backend and reduction-loop fusion (Phase 8);
`GDSGE_ASG_FIX_GRID` / unexercised ASG option branches; PCHIP; hook/`cxx` raw blocks;
performance work beyond parity (Phase 9).
