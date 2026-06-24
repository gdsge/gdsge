# Phase 7a Design — Spline-Path Widening (Barro_et_al_2017, Mendoza2010, GLSW2020)

- **Date:** 2026-06-12
- **Status:** Approved (design); pending implementation
- **Parent spec:** `2026-06-11-refactor-gdsge-design.md` (§15, Phase 7)
- **Owner:** Wenlan Luo

---

## 1. Context

Phases 0–6 proved the IR-centered architecture end-to-end on HeatonLucas1996: the public
`gdsge_codegen('HL1996')` → `iter_HL1996` → `simulate_HL1996` chain reproduces the old
toolbox's goldens through the new parser, IR, MATLAB codegen, and C++ (adept) codegen.

The parent spec's Phase 7 ("widen for backward-compat") covers the five remaining test
models. Exploration of those models' `.gmod` files and `test.m` drivers showed the phase
as specced bundles two very different kinds of work, plus items with no test coverage:

- **ASG** (`USE_ASG=1`: Bianchi2011_asg, CaoKS2016) is a structurally different iteration
  architecture — adaptive grid refinement inside the iteration loop, a different MEX data
  layout, new result fields (`asg_interp_struct` / `asg_output_struct`), and the public
  `asg` MATLAB class.
- The **spline-path models** (Barro_et_al_2017, Mendoza2010, GLSW2020) need only
  incremental parser/codegen/runtime widening.
- **No `.gmod` file in the entire corpus** (all six test models, everything under
  `base_package/`) uses a single macro (`#define`, `#for`, `#foreach`, `#if`, `#mat{}`,
  `#strcat_comma`, `include`, `cinclude`), nor `pre_iter`/`post_iter`/`pre_model` hooks,
  nor `cxx;…end;` blocks. There is no golden model to validate a macro engine against.
- The items Phase 6 deferred (`gdsge.m` orchestrator, `GenCodeSegment`, the legacy
  5-output `gdsge_codegen` signature) are used by **none** of the five remaining
  `test.m` drivers.

## 2. Decision — Phase 7 splits into 7a / 7b / 7c

| Sub-phase | Contents | Why |
|---|---|---|
| **7a (this spec)** | Barro_et_al_2017, Mendoza2010, GLSW2020 green end-to-end vs fresh goldens | Incremental widening; lands value soonest; de-risks 7b |
| **7b** | CaoKS2016, Bianchi2011_asg: ASG architecture, `var_tensor`, `asg_*_struct` result fields, `SkipModelInit`/shock-override options exercised by their drivers | ASG deserves its own focused spec → plan cycle |
| **7c** | Golden-less legacy surface: macro engine, `gdsge.m` orchestrator, `GenCodeSegment`, legacy 5-output signature; validated by synthetic fixtures (and a real macro-using example model if the owner supplies one) | All are frozen public surface (parent spec §13) but have no golden model; keeping them out of 7a/7b keeps those phases golden-testable |

`pre_iter`/`post_iter`/`pre_model`/`pre_jac_code`/`post_jac_code` hooks and `cxx` blocks
stay deferred alongside 7c unless a model that uses them enters the corpus.

## 3. Goal & success criteria

Phase 7a is done when, for each of the three models, a Slow-tagged gate drives the public
API exactly as the model's `test.m` does and matches a freshly captured golden:

1. `gdsge_codegen('<model>')` produces `iter_*.m`, `simulate_*.m`, `mex_*.cpp`,
   `<model>.gdsge.json` (schema-valid, round-trips), and a compiled MEX with cache-gated
   recompilation — same artifact contract as `tEndToEndHL1996`.
2. `iter_<model>` converges with iteration count and metric matching the golden within
   the established tolerance discipline, including each driver's structurally distinctive
   steps (Mendoza: no-arg call, then fine-grid warm-up re-solve; GLSW: `TolEq=1e-8`,
   `INTERP_ORDER=4`).
3. `simulate_<model>` reproduces the golden's reduced seeded simulation (shock path
   bit-exact; fields within tolerance), including GLSW's `init.*`-override batch.

## 4. Feature delta (what 7a must add)

From the gmod/driver inventory:

| Model | New surface exercised |
|---|---|
| Barro_et_al_2017 (`safe_assets.gmod`) | `GDSGE_MIN` reduction (parser + C++ lowering already exist; gate confirms) |
| Mendoza2010 | `model_init` block; `var_policy_init`/`var_aux_init`/`inbound_init`; 2-D state space; `SIMU_RESOLVE=0; SIMU_INTERP=1`; no-arg `iter_*()`; warm-up onto user-supplied finer grids; `SaveFreq=inf` |
| GLSW2020 | `model_init` + `*_init`; `INTERP_ORDER=4`; `SIMU_INTERP`; simulate option overrides `init.<var>`, `num_samples`, `num_periods`, `SimuPrintFreq`, `SimuSaveFreq` |

Already supported (verified in `src/`): all four reduction kinds lower in
`gdsge.codegen.cxx.emitModelBody`; `var_output` (HL1996 uses it); `adaptive()` bounds
(HL1996 uses `adaptive(1.5)`); the IR schema already defines the optional `modelInit`
section — the parser just never populates it.

## 5. Approach — batch goldens, front-end once, then vertical per model

Chosen over (a) pure model-by-model verticals (would design the shared `model_init`
front-end work mid-stream, twice) and (b) layer-by-layer across all three models
(end-to-end feedback arrives last). Hybrid order:

1. **Capture all three goldens** in one old-toolbox session (§6).
2. **Settle the parser/IR front-end for all three models** (§7) — the front-end gaps
   (`model_init`, `*_init` declarations) are shared by Mendoza and GLSW.
3. **Vertical codegen/runtime per model**, easiest first: Barro → Mendoza → GLSW (§8),
   each closed by its end-to-end golden gate.

## 6. Golden capture protocol

One capture session; each model in its **own `matlab -batch` process** with controlled
cwd and only `base_package/gdsge/source` on the path (parent-spec path policy). Committed
capture scripts (like HL1996's) record the exact settings and log wall-clock timings.

Per model: **full-convergence iter, reduced seeded simulate.**

- **Iter** runs exactly as `test.m` does, to full convergence — convergence is the
  correctness anchor. Mendoza's two-stage solve (default grids → 80×80 fine-grid
  `WarmUp` re-solve) is captured as-is because warm-up-onto-finer-grids is itself a
  feature under test.
- **Simulate** is shrunk to a small fixed size — at most 10 samples × 1,000 periods
  (smaller if the model's gmod defaults are smaller) — with `rng` seeded (HL1996
  precedent: `rng(0823)`); the exact values live in the committed capture script and are
  reused verbatim by the gate. The `test.m` plotting and moment-table code is skipped. GLSW additionally captures a tiny `init.*`-override
  batch (small `num_samples`, `num_periods=100`) to pin the override feature.
- Barro's `test.m` never calls simulate, but `simulate_safe_assets.m` is generated public
  surface, so a small seeded simulate is captured for it too.
- Captured `IterRslt`/`SimuRslt` land under `tests/<Model>/golden/` with a
  golden-integrity test each, mirroring HL1996.

If a full-convergence run proves prohibitively long (unknown until tried; timings are
logged), pause and revisit the protocol for that model before building gates on it.

## 7. Front-end work (parser → IR), settled once

- **`model_init` block:** `splitBlocks` learns the `model_init;…end;` block; its body is
  parsed by the **existing** model-statement parser (`parseModel`/`analyzeModel` reused
  with the init-variable tables); results populate the IR's existing `modelInit` section
  (`variables.policyInit/auxInit`, `bounds`, `statements`, `equations`).
- **Init declarations:** `parseVarDecls`/`parseDeclarations` learn `var_policy_init`,
  `var_aux_init`, `inbound_init` (same syntax as their non-init forms, incl. `adaptive()`
  on `inbound_init`), with an **independent slot layout** for the init problem (the init
  system is its own square system: #init equations == #init unknowns; validated like the
  main system).
- **Options:** `resolveOptions` and the IR options section carry `SIMU_INTERP` /
  `SIMU_RESOLVE` (mutually exclusive; exactly one active), `SimuPrintFreq`,
  `SimuSaveFreq`, and `INTERP_ORDER=4` (already constrained to {2,4} by the parent spec).
  Whatever of these the options whitelist already passes through stays as-is; gaps are
  filled with the same informative-error discipline.
- **`GDSGE_MIN`:** no parser work expected (Phase 3 implemented all four reductions);
  Barro's front-end gate is the confirmation.
- **Front-end gates (one per model):** `parseFrontEnd(<model>.gmod)` produces an IR that
  (a) passes the schema validator, (b) round-trips through JSON
  (`encode`/`decode`/`isequalIR`), and (c) satisfies targeted structural assertions —
  state/shock/policy/equation counts, `modelInit` presence and its square system, slot
  layout invariants. **No hand-authored reference IRs** — that was a bootstrap device for
  HL1996; here the validator + structural assertions carry the load, and the emitted IR
  JSON is committed as a snapshot regression artifact once the end-to-end gate is green.

## 8. Codegen + runtime work, vertical per model

### 8.1 Barro_et_al_2017 (first — validates the widening loop cheaply)

Expected near-free: `GDSGE_MIN` already lowers in C++ (`MIN(acc, body)` accumulation) and
the model uses no other new surface. Run the full pipeline, fix whatever actually breaks,
close the golden gate.

### 8.2 Mendoza2010 (the init task + simulate-interp + 2-D states)

- **C++ init task:** emit a second model body from `IR.modelInit` through the existing
  emitters; add `task_init` dispatch beside `task_inf_horizon` (old convention:
  `MEX_TASK_INIT` task id selected via the `MEX_TASK_NAME` scalar in the caller-workspace
  contract); compile with `-DHAS_INIT`; the data-buffer sizing follows the old formula
  `GDSGE_MAXDIM = max(num_policy_total, num_policy_init_total) + 4`. The shared
  `gdsge.codegen.dataLayout` descriptor gains the init-problem layout so the MATLAB
  packer and C++ `POP*` macros stay generated from one source.
- **Generated iter init segment:** before the main loop, gated **only** by
  `if ~SkipModelInit`, solve the init problem over the state tensor and use its
  policies/aux to seed `var_interp`. Old-toolbox semantics are replicated exactly: the
  init solve runs even when `WarmUp` is supplied (Mendoza's fine-grid re-solve does
  this), and the warm-up interpolants then supersede its output in the main loop —
  rebuilt on
  `gdsge.runtime` helpers (explicit locals, no v2struct), reporting through the same
  structured progress/unconverged machinery as the main loop. `SkipModelInit` joins the
  iter options whitelist now (its heavy use arrives with Bianchi in 7b).
- **`SIMU_INTERP` simulate variant:** generated `simulate_*` interpolates policies from
  the converged splines instead of re-solving each period (old
  `simulate_interp_template.m` vs `simulate_resolve_template.m` branch). Codegen selects
  the variant from the IR options; the resolve path stays the default.
- **2-D state space:** the spline stack (`myppual`, `constructSplines`, the C++ interp
  evaluator) is general-dimension by design; this model is where that claim is proven.
- **No-arg entry:** `iter_mendoza2010` (and all generated iter/simulate files) callable
  with zero arguments.
- **Warm-up onto finer grids:** `options.<state> = linspace(...)` + `options.WarmUp`
  re-solves on new grids seeded from the old solution's interpolants (existing
  `applyWarmUp`/`constructSplines` widened as needed).

### 8.3 GLSW2020 (interp order 4 + simulate overrides)

- **`INTERP_ORDER=4`:** order-4 not-a-knot spline construction in MATLAB and the
  `GDSGE_INTERP_ORDER` compile-time define in the MEX build (old `compile_template`
  convention), end-to-end.
- **Simulate runtime overrides:** `options.init.<var>` (initial states/shocks per
  sample), `options.num_samples`, `options.num_periods`, `SimuPrintFreq`, `SimuSaveFreq`
  resolve through the simulate options path with whitelist errors for unknowns;
  `SimuRslt` array shapes follow the overridden sizes (frozen field names).
- GLSW also consumes `model_init` and `SIMU_INTERP` built in §8.2.

## 9. Testing

- Harness unchanged: headless `matlab.unittest` via `matlab -batch "cd('tests'); run_tests"`,
  `junit.xml` + exit code authoritative; path policy one-source-per-process.
- New tests per model: golden-integrity, front-end gate (§7), end-to-end gate (§3) tagged
  Slow. Unit tests accompany each new mechanism at its own layer (init-declaration
  parsing, init slot layout, task dispatch emission, simulate-interp emission, order-4
  spline construction, simulate option resolution).
- Existing HL1996 gates must stay green throughout — widening must not disturb the
  vertical slice.
- Tolerances follow the established discipline: shock paths bit-exact (seeded), result
  fields within the comparison utility's relative/absolute bounds, iteration count and
  metric matched to the golden capture.

## 10. Risks

- **Unknown wall-clock for Mendoza/GLSW full-convergence solves** (and the 80×80
  re-solve). Mitigation: capture scripts log timings; if prohibitive, pause and revisit
  the protocol (§6) before gates are built. Gates are Slow-tagged either way.
- **Init-task data layout** is the deepest old-convention dependency (buffer sizing,
  task dispatch, caller-workspace contract). Mitigation: it flows through the existing
  shared `dataLayout` descriptor, and the Mendoza golden gate catches contract drift.
- **`SIMU_INTERP` fidelity:** the interp-simulate path produces different numerics than
  resolve-simulate by construction; goldens are captured per-path from the old toolbox,
  so comparisons stay like-for-like.

## 11. Out of scope

ASG and `var_tensor` (7b); macro engine, `gdsge.m`, `GenCodeSegment`, legacy 5-output
signature (7c); hook blocks and `cxx` blocks (deferred with 7c); SymPy backend (Phase 8);
performance work beyond parity (Phase 9).
