# Phase 2 — Parser Front-End → Partial IR — Design Spec

- **Date:** 2026-06-12
- **Status:** Approved (design); pending implementation plan
- **Owner:** Wenlan Luo
- **Parent design:** `docs/superpowers/specs/2026-06-11-refactor-gdsge-design.md` (§7 Parser design)
- **Builds on:** `docs/superpowers/specs/2026-06-12-phase1-ir-schema-design.md` (the IR contract)
- **Phase tracker:** `PROGRESS.md` (Phase 2)

---

## 1. Purpose & scope

Phase 2 builds the **front-end** of the new parser: the stages that turn `.gmod` text into a
**partial IR** — everything *except* the `model;…equations;…end;` body, which becomes an AST
in Phase 3. It is the first phase that actually *reads `.gmod`* (Phase 1's reference IR was
hand-authored), so it establishes the parser package, its eval sandbox, and the
section-by-section golden test that Phase 3 will complete.

The vertical-slice target is **HeatonLucas1996** (`HL1996.gmod`) — a clean "core" model: no
macros, no ASG, no `var_tensor`, no `model_init`, no hook blocks, no `cxx`. Phase 2 produces
the partial IR for this model and proves it against the (corrected) Phase-1 reference IR.

### In scope

- A `gdsge.parser` package with independently testable stages: `preprocess`, `splitBlocks`,
  `parseDeclarations`, `parseSimulate`, `resolveOptions`, `assemblePartialIR`, plus a
  top-level `parseFrontEnd` that chains them.
- The **eval sandbox**: running pre-model setup code (params, shock values, `shock_trans`,
  grids) in a controlled MATLAB workspace seeded with toolbox defaults, mirroring the old
  `gdsge_precode.m` mechanism.
- **Slot layout** assignment for policy and aux variables (fat front-end).
- The **`simulate` block** parsed into opaque-text IR fields.
- **Options resolution** from the real toolbox defaults (`default_mod.nmod` +
  `params_template.m`) with gmod overrides, mapped to the curated IR `options` subset.
- One **additive schema change**: a top-level `params` section (name → value), giving
  parameter values a home in the IR.
- **Tests:** a unit test per stage on small hand-built inputs, plus a golden partial-IR test
  that parses the real `HL1996.gmod` and matches the corrected reference IR section by
  section.

### Out of scope (Phase 2 non-goals)

- **No model-body parsing.** `model.statements` and `model.equations` stay empty lists;
  the tokenizer + recursive-descent AST is Phase 3.
- **No model-semantic validation** (square system, name resolution against the model body,
  reserved-word collisions in expressions). Phase 2 validates only what the partial IR can
  support (schema shape, slot layout, option invariants, bounds present for every policy).
- **No full `#…` macro engine.** Phase 2 ships the preprocess *stage* with a documented
  pass-through stub; `#define/#for/#foreach/#if/#mat/#strcat_comma/include/cinclude` are
  Phase 7 (where the models that use them come on board).
- **No codegen.** Phases 4+.
- **No new gmod features** beyond what HL1996 exercises: `model_init`, `var_tensor`, ASG, and
  hook-block *bodies* are recognized by `splitBlocks` only enough to not break; their full
  handling is Phase 7.

---

## 2. Decisions taken in brainstorming

1. **Macros — minimal preprocess only.** Phase 2 builds the preprocess stage (comments,
   line-continuation, deprecated-keyword rewrite, default-prepend) but defers the `#…` macro
   implementations to Phase 7. Rationale: HL1996 uses no macros, the phase plan lists
   "macros" under Phase 7's widen step, and this matches the spec's "minimum parser that
   grows only as models require."
2. **Fat front-end — slots + simulate in Phase 2.** The flat slot layout is mechanical
   (variable order + lengths) and needs no model body, so it is assigned here. The `simulate`
   block's fields are opaque text (no AST), so it is parsed here. Phase 3 is left with only
   the model AST and the final square-system check. This resolves the conflict between
   parent-spec §7 stage 3 (declaration parser produces slots) and the phase-plan line (slot
   layout under Phase 3) in favor of the declaration parser.
3. **Options — resolve correctly, fix the fixture.** Phase 2 resolves options from the actual
   toolbox defaults and corrects the committed `buildHL1996IR.m` + `HL1996.gdsge.json`, so
   the IR drives a numerically faithful Phase 6.
4. **New `params` IR section.** Parameter values currently have nowhere to live (the schema
   has no `params`; the model AST references `beta`/`gamma`/`Kb` as *name* nodes, not inlined
   constants). Phase 2 adds an additive `params` section.

---

## 3. The schema change: a `params` section

### 3.1 Why it is needed

In the old toolbox, `get_parameters` evals the pre-model setup and captures
`params = v2struct(beta, gamma, Kb)` — a name→value struct later `POPN`'d into the C++
solver and emitted into the generated MATLAB. The Phase-1 IR omitted this: the model AST in
`buildHL1996IR` references `beta`/`gamma`/`Kb` by *name*, but no IR section carries their
*values*. Without a home for values, no backend can emit `beta = 0.95;`.

### 3.2 Shape

Add to the schema root:

```
params : fList( fStruct( name: fText(), value: fMatrix() ) )
```

- `value` is `fMatrix` because parameters may be arrays (the old toolbox tracks
  `parameters_size = numel(param)`); scalars are the 1×1 case and round-trip via the
  existing matrix-wrapping rule.
- The list preserves declaration order.
- Empty list is valid (a model with no `parameters`), exactly like `variables.others`.

For HL1996: `params = { {name:'beta', value:0.95}, {name:'gamma', value:1.5},
{name:'Kb', value:-0.05} }`.

### 3.3 Ripple

- `src/+gdsge/+ir/schema.m` — add the `params` field-spec to `s.root`.
- `docs/ir-schema.md` — regenerated from the schema (the no-drift test enforces this).
- `tests/+gdsgefix/minimalIR.m` — add a `params` field (may be a single param or empty).
- `tests/HeatonLucas1996/ir/buildHL1996IR.m` — add the three params.
- `tests/HeatonLucas1996/ir/HL1996.gdsge.json` — regenerated golden.

`irVersion` stays `1.0.0` (the freeze is Phase 9; nothing external consumes the IR yet).

### 3.4 Name resolution note (forward to Phase 3)

Once `params` exists, parameter names become a resolvable pool for the model AST. Phase 2
does **not** wire that check (it has no model AST). Phase 3's semantic analysis adds params
to the name-resolution pools alongside states, shocks, and policy/aux vars.

---

## 4. Options reconciliation

### 4.1 Two defaults sources

The old toolbox resolves options from two places:

- **`code_template/default_mod.nmod`** — prepended to the gmod and eval'd at parse time.
  Supplies parse-time flags: `INTERP_ORDER=4`, `EXTRAP_ORDER=2`, `USE_SPLINE=1`,
  `USE_ASG=0`, `USE_PCHIP=0`, `SIMU_INTERP=0`, `SIMU_RESOLVE=1`, `AsgMaxLevel=10`,
  `AsgThreshold=1e-2`, `shock_num=1`, `shock_trans=1`, `NumThreads=feature('numcores')`, …
- **`code_template/params_template.m`** — emitted into the generated `iter_*.m`. Supplies
  runtime defaults: `TolEq=1e-6`, `PrintFreq=10`, `NoPrint=0`, `SaveFreq=10`,
  `SimuPrintFreq=1000`, `SimuSaveFreq=inf`, `MaxIter=inf`, …

A gmod may override any of these by assigning the flag; the runtime `options` struct passed to
`iter_*` overrides the `params_template` values at run time (not at parse time).

### 4.2 Resolution in Phase 2

`resolveOptions` reads the eval'd flag workspace (which already contains `default_mod`
defaults plus any gmod overrides) and the `params_template` runtime defaults, then maps to
the curated IR `options` subset:

| IR option       | Source                                                            |
|-----------------|-------------------------------------------------------------------|
| `interpMethod`  | `USE_SPLINE/USE_ASG/USE_PCHIP` → `'spline'/'asg'/'pchip'`          |
| `interpOrder`   | `INTERP_ORDER`                                                     |
| `extrapOrder`   | `EXTRAP_ORDER`                                                     |
| `asgMaxLevel`   | `AsgMaxLevel` (only when `interpMethod='asg'`)                     |
| `asgThreshold`  | `AsgThreshold` (only when `interpMethod='asg'`)                    |
| `tolEq`         | `TolEq` (`params_template`)                                        |
| `numThreads`    | `NumThreads` (machine-derived `feature('numcores')`)              |
| `simuResolve`   | `SIMU_RESOLVE`                                                     |
| `simuInterp`    | `SIMU_INTERP`                                                      |
| `printFreq`     | `PrintFreq` (`params_template`)                                    |
| `saveFreq`      | `SaveFreq` (`params_template`)                                     |

### 4.3 Fixture corrections

The Phase-1 reference IR's `options` were hand-authored approximations. Phase 2 corrects them
to the resolved defaults:

- `simuResolve`: `0 → 1`
- `simuInterp`: `1 → 0`
- `printFreq`: `100 → 10`
- `saveFreq`: `Inf → 10`
- `numThreads`: **not** changed in the fixture. The fixture keeps a fixed value (`8`) so its
  own deterministic JSON golden (`HL1996.gdsge.json`, compared byte-for-byte by
  `tIrHL1996`) stays reproducible. Only the *parse* of `HL1996.gmod` yields a machine-derived
  `feature('numcores')`, so the **Phase-2 front-end test** (§7.2) compares `numThreads`
  loosely (positive integer) — it is the one field where the parsed value legitimately
  differs from the fixture across machines.

The runtime `options.SaveFreq = inf` in `test.m` is a *run-time* override, not a parse-time
value, and does not change the parsed IR.

---

## 5. Pipeline architecture

All modules live under `src/+gdsge/+parser/`. Each is a single function file with focused
unit tests. A top-level `gdsge.parser.parseFrontEnd(gmodText)` chains them and returns a
validated partial IR.

```
gmod text
   │  preprocess(text)            comments, continuations, deprecated rewrite, prepend defaults
   ▼
preprocessed text + default-flag workspace seed
   │  splitBlocks(lines)          declaration/setup lines  ┊  named blocks (model, simulate, hooks)
   ▼
{ declLines, blocks }
   │  parseDeclarations(declLines)  two-pass: learn names, then classify & route assignments
   │     ├─ eval sandbox  → params, shock values, shock_trans, option flags, grid sizes
   │     └─ text capture  → grids, interp initial/update exprs
   ▼
{ params, shocks, states, variables(+slots), bounds, interp, flags }
   │  parseSimulate(blocks.simulate)   num_periods, num_samples, initial, var_simu, transitions
   │  resolveOptions(flags)            curated IR options
   ▼
assemblePartialIR(...)            build struct, model.statements/equations = {}, gdsge.ir.validate
   ▼
partial IR
```

### 5.1 `preprocess(text) → {text, ...}`

- Strip `%` comments (respecting that `%` inside strings/`#mat{}` is not a comment — for
  HL1996 the simple line-comment rule suffices; the stage is structured so Phase 7 can refine).
- Join MATLAB line-continuations (`...`).
- Deprecated-keyword rewrite (`GNDSGE → GDSGE`, etc.) via word-boundary replacement.
- Prepend `default_mod.nmod` content so its defaults are part of the setup code that the eval
  sandbox runs.
- **Macro hook:** a pass-through that scans for any `#…` directive and raises a clear
  "macro support arrives in Phase 7" error if one is present. HL1996 has none.

### 5.2 `splitBlocks(lines) → {declLines, blocks}`

- Recognize block openers (`model;`, `simulate;`, and the hook names `pre_model`, `pre_iter`,
  `post_iter`, `pre_jac_code`, `post_jac_code`, `model_init`) and their matching `end;`.
- Everything outside a block is a declaration/setup line.
- Each block and line keeps its **source line number** for diagnostics.
- For HL1996 only `model` and `simulate` blocks appear; hook blocks are captured as raw text
  into `hooks` (all empty for HL1996) but their bodies are not interpreted in Phase 2.

### 5.3 `parseDeclarations(declLines) → declaration tables`

Two passes:

**Pass 1 — learn names.** Scan declaration keywords to collect the name sets:
`parameters`, `var_shock`, `var_state`, `var_policy` (with `[N]` array suffix →
length `N`), `var_aux`, `var_interp`, `var_output`, `var_others`, and the `inbound` /
`initial` declarations. (`var_tensor`, `var_policy_init`/`var_aux_init`, and `inbound_init`
are recognized only as known keywords so an HL1996-shaped file does not trip
`unknownDeclaration`; their full handling — together with `model_init` — is Phase 7. HL1996
uses none.)

**Pass 2 — classify & route assignment lines** by LHS name:

| LHS class                      | Action                                                      |
|--------------------------------|-------------------------------------------------------------|
| param / shock var / `shock_num` / `shock_trans` | add to **pre-model setup script** (eval'd)  |
| state name                     | add to setup script **and** capture RHS **text** as grid    |
| `var_interp` name              | capture RHS **text** as `updateExpr`; **exclude from eval**  |
| (declaration keyword)          | structural; handled in pass 1                                |

`initial <interp> <expr>` lines supply each interp object's `initialExpr` (text).
Interp `args` default to `states.names`.

Outputs the structured tables: `params`, `shocks` (names, count, values, transitions),
`states` (names, grid text), `variables` (policy/aux/interp/tensor/output/others with
lengths), `bounds` (lower/upper/adaptiveFactor), `interp` objects, and the eval'd
option-flag workspace.

### 5.4 The eval sandbox (the subtle part)

The pre-model setup script (collected in pass 2) is `eval`'d in a **dedicated function
workspace** seeded by `default_mod.nmod`, reproducing the old `gdsge_precode.m`:

- **Captured numerically:** parameter values (`beta=0.95`, arrays allowed), shock realized
  values (`g`, `d`, `eta1`), and `shock_trans` (after the gmod's own normalization
  `./ repmat(sum(...))`). The normalized matrix is the IR's `shocks.transitions.shock_trans`.
- **Captured as text, not value:** state grids. The RHS source of `w1 = linspace(...)` is
  stored verbatim as `states.grids.w1`; the eval'd vector is used only to know the grid is
  valid / its size, not stored.
- **Never eval'd:** interp-update assignments (`ps_future = ps`) — they reference policy
  unknowns absent at parse time. Their RHS is captured as text in pass 2.
- **Error handling:** eval failures are reported with the offending line echoed (mirroring the
  old `get_parameters` catch that prints numbered pre-code lines), wrapped in a
  `gdsge:parser:setupEvalFailed` error.

Isolation: the sandbox runs in its own function scope (no `assignin('base',…)`), so it never
touches the caller workspace or the MATLAB path (consistent with the project's path policy).

### 5.5 Slot layout

- Policy vars: walk in declaration order, accumulating `slot = [start stop]` from each
  variable's length (scalars length 1; `w1n[8]` length 8). HL1996 →
  `c1[1 1] … pb[11 11] w1n[12 19]`.
- Aux vars: independent slot space starting at 1 → `equity_premium[1 1]`.
- This is exactly the layout `gdsge.ir.validate`'s slot checks already enforce (contiguous,
  starts at 1, length matches span).

### 5.6 `parseSimulate(block) → simulate section`

Parse the `simulate;…end;` block lines:

- `num_periods = 10000` → `numPeriods` (scalar)
- `num_samples = 24` → `numSamples` (scalar)
- `initial w1 0.5`, `initial shock 1` → `initial` list of `{var, value(text)}`
- `var_simu c1 c2 ps pb equity_premium` → `varSimu` (refs into policy/aux)
- `w1' = w1n'` → `transitions` list of `{state:'w1', expr:'w1n'}` (RHS text, prime stripped)

### 5.7 `assemblePartialIR(...) → partial IR`

Assemble the full struct with `model.statements = {}` and `model.equations = {}`, set
`irVersion`/`modelName`, then run `gdsge.ir.validate`. The partial IR must validate (the
schema permits empty statement/equation lists; the square-system check is parser-side and
deferred to Phase 3, so an empty model body is not a validation error here).

---

## 6. The Phase 2 ↔ Phase 3 boundary

| Concern                              | Phase 2 | Phase 3 |
|--------------------------------------|:-------:|:-------:|
| Preprocess / block-split             |   ✅    |         |
| Declarations, params, shocks, states |   ✅    |         |
| Slot layout                          |   ✅    |         |
| Bounds, interp objects               |   ✅    |         |
| `simulate` section                   |   ✅    |         |
| Options resolution                   |   ✅    |         |
| `model;…end;` body → AST statements  |         |   ✅    |
| `equations;…end;` → AST + prime expand |       |   ✅    |
| Square-system / name-resolution check |        |   ✅    |
| Full `#…` macros                     |         | Phase 7 |

---

## 7. Testing & TDD

Following the project's TDD discipline (failing test first, small commits), with the Phase-0
headless harness (`tests/run_tests.m`, `tests/run.ps1`).

### 7.1 Stage unit tests (`tests/parser/`)

- `tPreprocess` — comment strip, line-continuation join, deprecated rewrite, default prepend,
  `#`-directive → clear Phase-7 error.
- `tSplitBlocks` — separates decl lines from `model`/`simulate` blocks; line numbers correct;
  unterminated block → error.
- `tParseDeclarations` — small fragments: a `var_policy` with a `[N]` array → right length;
  `inbound … adaptive(f)` → `adaptiveFactor`; LHS classification routes a param vs an interp
  update correctly; missing bound for a policy var → error.
- `tEvalSandbox` — `shock_trans` normalization reproduced; a param array captured with right
  size; a deliberately broken setup line → `gdsge:parser:setupEvalFailed` with the line echoed.
- `tParseSimulate` — the five simulate constructs parsed to the right fields.
- `tResolveOptions` — flag-workspace → curated options mapping; ASG path pulls
  `asgMaxLevel`/`asgThreshold`; spline path omits them.

### 7.2 Golden partial-IR test (`tests/HeatonLucas1996/parser/`)

- `tFrontEndHL1996` — `parseFrontEnd(fileread('HL1996.gmod'))`, then for each populated
  section assert `gdsge.ir.isequalIR` against the corresponding section of the corrected
  `buildHL1996IR()` with `model.statements`/`equations` blanked. `numThreads` is checked
  loosely (positive integer). The whole partial IR must `gdsge.ir.validate`.

Single source of truth stays the reference IR: the golden test reads sections off
`buildHL1996IR`, so correcting the fixture (§3, §4) automatically updates the target.

### 7.3 Regression

The existing Phase-1 IR tests (`tIrHL1996`, round-trip, validate, gendoc no-drift) must stay
green after the `params` schema change and the fixture corrections.

---

## 8. Error handling

Parser errors are precise and located (block/line where possible), using a
`gdsge:parser:*` identifier namespace:

- `gdsge:parser:macroUnsupported` — a `#…` directive (Phase 7).
- `gdsge:parser:unterminatedBlock` — a block opener with no matching `end;`.
- `gdsge:parser:setupEvalFailed` — pre-model eval error, with the numbered setup line echoed.
- `gdsge:parser:missingBound` — a policy variable with no `inbound`.
- `gdsge:parser:unknownDeclaration` — an unrecognized declaration keyword.

Model-body and square-system diagnostics belong to Phase 3.

---

## 9. File plan

Created in `src/+gdsge/+parser/`:

- `parseFrontEnd.m`, `preprocess.m`, `splitBlocks.m`, `parseDeclarations.m`,
  `evalSetup.m` (the sandbox), `parseSimulate.m`, `resolveOptions.m`, `assemblePartialIR.m`,
  `Contents.m`.

Modified:

- `src/+gdsge/+ir/schema.m` — add `params`.
- `docs/ir-schema.md` — regenerated.
- `tests/+gdsgefix/minimalIR.m`, `tests/HeatonLucas1996/ir/buildHL1996IR.m`,
  `tests/HeatonLucas1996/ir/HL1996.gdsge.json` — `params` + options corrections.

Tests:

- `tests/parser/tPreprocess.m`, `tSplitBlocks.m`, `tParseDeclarations.m`, `tEvalSandbox.m`,
  `tParseSimulate.m`, `tResolveOptions.m`.
- `tests/HeatonLucas1996/parser/tFrontEndHL1996.m`.

Reference inputs (read-only, already present): `tests/HeatonLucas1996/HL1996.gmod`,
`base_package/gdsge/source/code_template/default_mod.nmod`, `…/params_template.m`.

---

## 10. Risks

- **Comment/string edge cases.** The minimal comment stripper assumes `%` always starts a
  comment. HL1996 holds, but a `%` inside a string literal would break it. Mitigated by
  keeping `preprocess` a seam Phase 7 can harden; flagged here honestly.
- **Grid-text vs eval.** Storing grids as text but eval'ing them for validation could diverge
  if a grid references a parameter defined later. HL1996's grid is self-contained; the
  single-pass setup script (params before grids in source order) preserves this. General
  ordering is a Phase 7 concern.
- **Machine-dependent `numThreads`.** Pinning it would make the golden machine-specific;
  the loose check (positive integer) avoids that while still exercising the path.
- **Fixture correction churn.** Changing the committed golden is deliberate (call #3) and
  gated by the regenerated JSON + the Phase-1 regression tests staying green.

---

## 11. Acceptance criteria

1. `gdsge.parser.parseFrontEnd(fileread('HL1996.gmod'))` returns a struct that
   `gdsge.ir.validate` passes.
2. Every populated section of that struct `gdsge.ir.isequalIR`-matches the corrected
   `buildHL1996IR()` reference (model body blanked; `numThreads` loose).
3. All stage unit tests green; all Phase-0/Phase-1 tests still green.
4. `docs/ir-schema.md` regenerated; no-drift test green.
5. `pwsh -File tests/run.ps1` exits 0.
