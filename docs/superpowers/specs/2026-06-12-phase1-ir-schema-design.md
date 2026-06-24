# Phase 1 — IR Schema + MATLAB Scaffolding — Design Spec

- **Date:** 2026-06-12
- **Status:** Approved (design); pending implementation plan
- **Owner:** Wenlan Luo
- **Parent design:** `docs/superpowers/specs/2026-06-11-refactor-gdsge-design.md` (§8 IR design)
- **Phase tracker:** `PROGRESS.md` (Phase 1)

---

## 1. Purpose & scope

Phase 1 builds the **IR contract** — the versioned, language-neutral intermediate
representation that the parser (Phases 2–3) emits and every backend (Phases 4–5, 8)
consumes. The IR is the seam of the whole refactor: backends never re-parse `.gmod`, so the
IR's shape, its MATLAB ↔ JSON round-trip, and its validator are what every later phase is
built on.

This phase delivers **only the contract**: the `gdsge.ir` package, its schema descriptor,
serialization, validators, doc generator, AST-node API, and a hand-authored HL1996 reference
IR that doubles as Phase 3's acceptance target.

### In scope

- A single declarative **schema descriptor** (`gdsge.ir.schema`) that is the source of truth.
- **Validation** (`gdsge.ir.validate`) driven by the descriptor.
- **Serialization** — `encode`/`decode`/`canonicalize` — descriptor-driven, type-aware,
  defeating MATLAB's JSON round-trip quirks.
- **Doc generation** (`gdsge.ir.gendoc`) → `docs/ir-schema.md`, so the doc cannot drift.
- **AST node API** (`gdsge.ir.node.*` constructors + accessors).
- **Full-feature-inventory** schema: representations for every gmod feature in the
  backward-compat surface (not just what HL1996 uses).
- **Tests:** synthetic fragment tests for every rule/quirk, plus a complete hand-authored
  **HL1996 reference IR** that validates, round-trips, and is committed as a JSON golden.

### Out of scope (Phase 1 non-goals)

- **No `.gmod` parsing.** The HL1996 reference IR is hand-authored, not generated.
- **No codegen** consuming the IR.
- **No model-semantic validation** (square system, bounds-completeness, reserved-word
  collisions, name resolution). Those need gmod/parse context and belong to the parser's
  semantic-analysis stage (Phase 3). See §6 for the precise boundary.

---

## 2. Decisions taken in brainstorming

1. **Schema ambition — full feature inventory now.** The schema defines representations for
   the entire backward-compat surface (var_tensor, ASG, model_init, hook blocks, custom
   transition matrices, multi-state, all options), even though only HL1996 is exercised
   end-to-end in this phase. Features no test model uses are *schema-defined but not yet
   round-trip-proven against a real model* — flagged honestly in §5.7.
2. **Source of truth — MATLAB validator + auto-generated doc.** A declarative descriptor in
   MATLAB is canonical; `docs/ir-schema.md` is generated from it. No separate hand-maintained
   JSON-Schema file to drift. The default path stays pure MATLAB (no Python).
3. **Round-trip — descriptor-driven, type-aware (Approach C).** The same descriptor that
   drives validation and doc-gen also drives (de)serialization, so the round-trip quirks are
   handled once, type-aware, in a single small walker rather than scattered special-cases.
4. **Test target — both.** Fine-grained synthetic fragments *and* a hand-authored complete
   HL1996 reference IR (the parser's eventual golden target).

---

## 3. Module surface & package layout

```
src/+gdsge/+ir/
  schema.m          % THE source of truth. Returns the descriptor: every section, its
                    %   fields, each field's KIND, required/optional, enum domains, and
                    %   the current irVersion.
  validate.m        % validate(ir) -> report{pass, errors{}}. Interprets schema.m. Located,
                    %   human-readable errors ("model.statements{3}.target: missing").
  encode.m          % ir struct -> JSON text (pre-normalize, then jsonencode).
  decode.m          % JSON text -> canonical ir struct (jsondecode, then canonicalize).
  canonicalize.m    % schema-guided fixups for the round-trip quirks (§5).
  isequalIR.m       % structural IR equality, orientation/array-kind tolerant (for tests).
  gendoc.m          % schema.m -> docs/ir-schema.md (generated; never hand-edited).
  roundtrip.m       % convenience: decode(encode(ir)); used heavily in tests.
  +node/            % AST node API:
      num.m  name.m  primed.m  unop.m  binop.m  call.m  index.m   % constructors
      kindOf.m  children.m                                        % accessors
  Contents.m
```

`+node` is a sub-package of small constructors because it is the API the Phase 2–3 parser
will call thousands of times (`gdsge.ir.node.binop('+', a, b)`). Fixing the node shape now,
behind constructors, means the round-trip and validator have one node definition to target,
and the parser writes intention-revealing code instead of raw struct literals.

---

## 4. IR document shape

Two principles cut the schema's complexity roughly in half.

### Principle A — two expression worlds

The IR carries expressions in exactly two forms and never conflates them:

- **Model AST** — expressions *inside* `model;…equations;…end;` and `model_init`. These
  lower to C++ (and MATLAB) and are differentiated. Parsed to the 7-node AST (§4.2).
- **Opaque MATLAB text** — grids, params, `shock_trans`, `inbound` bounds, `var_interp`
  `initial`/update rules, `simulate` transitions. These are `eval`'d in MATLAB at
  setup/runtime, never differentiated. Carried verbatim as strings (`text` kind).

This is why grids/bounds/initial-expressions need no expression parser: they ride along as
text, exactly as the old toolbox `eval`s them. Only the model body — the part that becomes
C++ and gets differentiated — needs an AST.

### Principle B — reductions are flat statements; nested ones are hoisted

HL1996 has `equity_premium = GDSGE_EXPECT{…} - 1/pb;` — a reduction used as a *sub-term*.
Rather than make reductions nestable AST nodes, the parser (Phase 3) **hoists** each
reduction to its own statement with a generated target (`__r1 = GDSGE_EXPECT{…}`) and
references the temp in the enclosing expression. So in the IR, every reduction is a
top-level `reduction` statement with an explicit `target`. This keeps `model.statements` a
flat list and matches the spec §10 "one reduction → one fused loop" lowering. (This refines
parent-spec §8, which implied reductions are only ever a whole RHS. Hoisting is a parser
normalization; it does not change the schema, which simply requires `reduction.target`.)

### 4.1 Top-level sections

| Section | Kind | Shape (HL1996 unless noted) |
|---|---|---|
| `irVersion` | text | `"1.0.0"` (semver; see §7) |
| `modelName` | text | `"HL1996"` |
| `options` | struct | resolved flags, concrete values after defaults: `interpMethod:"spline"`, `interpOrder`, `extrapOrder`, `asgMaxLevel`, `asgThreshold`, `tolEq`, `numThreads`, `simuResolve`, `simuInterp`, `printFreq`, `saveFreq`, `simuPrintFreq`, `simuSaveFreq`, … |
| `shocks` | struct | `names:["g","d","eta1"]`, `count:8`, `values:{g:[…],d:[…],eta1:[…]}` (each 1×count), **`transitions`** registry `{shock_trans: 8×8, …}` |
| `states` | struct | `names:["w1"]`, `grids:{w1:"linspace(-0.05,1.05,201)"}` (text); multi-state = N named grids (Mendoza2010 `cTilde k`, CaoKS2016 `K X`) |
| `variables` | struct | category lists `policy / aux / interp / tensor / output / others`; each policy/aux entry `{name, length, slot:[start,stop]}`; `output` = name refs |
| `bounds` | list | `[{name:"c1",lower:"0",upper:"1",adaptiveFactor:null},{name:"ps",lower:"0",upper:"3",adaptiveFactor:1.5},…]` (lower/upper = MATLAB text) |
| `interp` | list | var_interp objects `[{name:"ps_future",args:["w1"],initialExpr:"0.0",updateExpr:"ps"},{name:"c1_future",…,initialExpr:"w1.*d+eta1",updateExpr:"c1"},…]` (text) |
| `model` | struct | `statements:[…]` (assign/interpCall/reduction) + `equations:[…]` (AST) |
| `modelInit` | optional struct | `{variables:{policyInit,auxInit}, bounds, statements, equations}` (GLSW2020-shaped); absent for HL1996 |
| `simulate` | struct | `numPeriods:10000`, `numSamples:24`, `initial:[{var:"w1",value:"0.5"},{var:"shock",value:"1"}]`, `varSimu:["c1","c2","ps","pb","equity_premium"]`, `transitions:[{state:"w1",expr:"w1n"}]` (text) |
| `hooks` | struct | opaque strings `preModel, preIter, postIter, preJacCode, postJacCode, cxx` (all `""` for HL1996) |

Notes:

- **`variables` vs. the `interp` section.** `variables.policy`/`aux` carry the slot layout
  (`{name, length, slot}`). `variables.interp`/`output`/`others`/`tensor` are **name
  references only** — their definitions live elsewhere: the `interp` *section* holds the
  var_interp object bodies (args/initialExpr/updateExpr); `output`/`others` names point back
  into `policy`/`aux`. The validator (§6) enforces that these references resolve.
- **Transition-matrix registry.** Reductions may name a custom matrix
  (`GDSGE_EXPECT{ expr | shock_trans2 }`, CaoKS2016). `shocks.transitions` is therefore a
  *registry* of named N×N matrices (default key `shock_trans` plus any named extras);
  `reduction.transRef` names one.
- **ASG.** `USE_ASG=1` coexists with the `var_state` grids (the grid is the base level the
  adaptive grid refines from). ASG parameters (`asgMaxLevel`, `asgThreshold`) live in
  `options`; nothing special is needed in `states` beyond the base grid text.
- **Slots.** `slot:[start,stop]` is a 1-based inclusive index range within the variable's
  category-packed array; `length` (1 = scalar, N = array such as `w1n[8]`) must match the
  span. The exact numbers are computed during Phase 3 semantic analysis; the hand-authored
  HL1996 reference IR carries the real layout (`c1→[1,1] … w1n→[12,19]`, a 19-wide policy
  vector).

### 4.2 Model statement types (the AST world)

- `assign{ target, primed:bool, expr:AST }` — `b1p = nb1p + Kb`; `primed:true` for
  future-indexed targets like `w1_consis'`.
- `interpCall{ targets[], primed:bool, args[], interpRef }` —
  `[psn',pbn',c1n',c2n'] = GDSGE_INTERP_VEC'(w1n')`.
- `reduction{ kind:EXPECT|MIN|MAX|PROD, target, body:AST, transRef }` — `transRef` names a
  key in `shocks.transitions` (default `"shock_trans"`).
- `equation{ expr:AST, primed:bool }` — `-1+beta*es1+ms1`; `primed:true` (`w1_consis'`)
  expands at codegen to `count` equations (one per shock).

`model.statements` holds the body (assign/interpCall/reduction, ordered); `model.equations`
holds the residual list — mirroring the gmod `model; <body> equations; <residuals> end;`.

### 4.3 Expression AST node set (7 nodes, tagged by `kind`)

`num{value}` · `name{id}` · `primed{id}` (trailing `'` future/next-state) · `unop{op,arg}` ·
`binop{op,lhs,rhs}` · `call{fn,args[]}` · `index{base,args[]}`.

Reduction bodies are ordinary ASTs over this set — e.g.
`g'^(1-gamma)*(c1n'/c1)^(-gamma)*(psn'+d')/ps` mixes `primed` nodes (`g`,`c1n`,`psn`,`d`)
with `name` nodes (`c1`,`ps`,`gamma`). No dedicated `reduce` node exists (Principle B).

---

## 5. Serialization & the round-trip quirks

### 5.1 Field kinds

`schema.m` tags every field with one kind; the generic walker in `encode`/`decode`/
`canonicalize`/`validate` understands this small vocabulary:

`scalar` · `matrix` · `enum` · `text` · `ref` · `list` (array of homogeneous records) ·
`nodeList` (array of heterogeneous AST nodes/statements) · `struct` (nested section) ·
`optional` (wrapper marking a field that may be absent).

### 5.2 Quirk 1 — orientation collapse

`jsondecode(jsonencode([1 2 3]))` returns a **3×1 column**, silently transposing grids and
`shock_trans`. Fix: `matrix` fields **always** encode as an array-of-rows
(`[1 2 3]` → `[[1,2,3]]`; an 8×8 stays nested; an N×1 → `[[…],[…],…]`). `jsondecode`
reconstructs the shape faithfully, and the JSON stays human-readable.

### 5.3 Quirk 2 — array-kind ambiguity

A JSON array of objects decodes to a **struct array** when objects share fields, but a
**cell array** when they differ — so an AST node-list is sometimes one, sometimes the other.
Fix: every object-array (`list`, `nodeList`) canonicalizes to **one** form — a **cell array
of scalar structs**. Codegen and validators always iterate `for k = 1:numel(cellArr)`.

### 5.4 Quirk 3 — scalar/empty flattening

The descriptor knows each field's expected kind, so `canonicalize` coerces: `scalar` → 1×1;
empty `text` → `""`; missing `optional` → `[]`; a singleton `list` stays a 1-element cell.

### 5.5 Quirk 4 — non-finite numerics

`jsonencode(Inf)` silently becomes `null` (and `SaveFreq = inf` is real — Mendoza2010,
GLSW2020). Fix: the walker encodes non-finite numbers as tagged strings
(`"Inf"`, `"-Inf"`, `"NaN"`) and `decode` restores them to numeric.

### 5.6 Equality

`isequalIR(a, b)` compares two IRs through the canonical lens (orientation- and
array-kind-tolerant, non-finite-aware), so round-trip tests assert genuine structural
equality rather than incidental encoding differences.

### 5.7 Honesty about coverage

The schema *defines* the full inventory, but features no test model exercises — `var_tensor`,
`var_others`, `GDSGE_PROD`/`GDSGE_MAX`, `USE_PCHIP`, `pre_iter`/`post_iter`,
`pre_jac_code`/`post_jac_code`, `cxx`, macros (which never reach the IR; they expand in the
Phase 2 preprocess step) — are schema-defined and unit-testable with synthetic fragments,
but are **not** round-trip-proven against a real model until the relevant model arrives
(Phase 7+). This is stated in `docs/ir-schema.md` so no one mistakes schema coverage for
validated coverage.

---

## 6. Validation boundary

The rule: **if it is checkable from the IR data alone, `gdsge.ir.validate` checks it; if it
needs gmod/parse context, the parser checks it (Phase 3).**

`gdsge.ir.validate(ir)` → `{pass, errors[]}` with located messages, checks:

- **Shape** — every required section/field present and matching its declared kind; AST nodes
  recursively well-formed (`binop` has `lhs`+`rhs`, `call` has `fn`+`args`, etc.).
- **Internal references** — `reduction.transRef` ∈ `shocks.transitions`;
  `variables.output`, `simulate.varSimu`, and `simulate.transitions[].state` names resolve;
  slot ranges well-formed (start ≤ stop, contiguous and non-overlapping within a category,
  `length` matches span).
- **Enum / option domains** — `reduction.kind` ∈ {EXPECT,MIN,MAX,PROD}; `node.kind` valid;
  `interpOrder` ∈ {2,4}; exactly one of spline/asg/pchip selected; ASG excludes `var_tensor`.

**Deferred to the parser** (needs gmod context, not pure IR): reserved-word collisions,
name resolution against declarations, and the *model-soundness* checks. In particular,
**square-system** (#equations == #unknowns after prime/array expansion) and
**bounds-completeness** (every policy has bounds) are kept **parser-side primary** — the
parser must check them to emit a correct IR anyway — and are deliberately *not* duplicated
in the Phase 1 `+ir` validator, to keep this phase tight. (They may be added later as a
defense-in-depth net if a need appears; that is a future decision, not Phase 1 scope.)

---

## 7. Versioning

- `irVersion` is **semver** (`MAJOR.MINOR.PATCH`), starting `"1.0.0"`.
- `schema.m` owns the current version. `decode`/`validate` require a matching **MAJOR**
  (breaking changes bump MAJOR and are rejected by older readers with a clear, located
  error); additive **MINOR** changes (new optional sections) are accepted.
- No backend consumes the IR until Phase 4, so versions churn freely through Phases 2–6; the
  real freeze is Phase 9. The *mechanism* matters now, the *number* does not.

---

## 8. Testing strategy

Native `matlab.unittest`, picked up by the existing `tests/run_tests.m`
(`fromFolder`, `IncludingSubfolders`). TDD throughout: a failing test per stage first
(`schema`/`validate`/`encode`/`decode`/`node`/`gendoc`), then implement; small, frequent
commits.

### 8.1 Synthetic fragment tests (`tests/ir/`)

- `tSchemaRoundtrip.m` — each quirk in isolation: a `matrix` field stays 1×8 (not 8×1);
  8×8 `shock_trans` survives; a mixed `nodeList` → cell-of-structs; `SaveFreq=inf` survives;
  `NaN`, empty list, and empty text round-trip.
- `tValidate.m` — pass/fail per check: missing required field → located error; bad `enum`;
  dangling `transRef`; overlapping slots; `interpOrder=3`; two interp methods selected; a
  clean fragment passes.
- `tNodeConstructors.m` — `gdsge.ir.node.binop('+',a,b)` builds the right tagged struct;
  `kindOf`/`children` accessors behave; a small AST round-trips.
- `tGendoc.m` — `gendoc()` covers every descriptor section, and the committed
  `docs/ir-schema.md` **equals** fresh `gendoc()` output (fails if stale → the doc cannot
  drift).

### 8.2 HL1996 reference IR (`tests/HeatonLucas1996/ir/`)

- `buildHL1996IR.m` is the **source** (a MATLAB fixture builder producing the full IR struct:
  all sections, the 19×19 system, the `GDSGE_EXPECT` reductions, the primed vector interp
  `GDSGE_INTERP_VEC'`, the future-indexed equation `w1_consis'`, `adaptive(1.5)` bounds, the
  KKT complementarity equations).
- `HL1996.gdsge.json` is its committed, generated golden.
- `tIrHL1996.m` asserts: `validate` passes clean; `isequalIR(ir, decode(encode(ir)))`;
  and `encode(ir)` matches the committed JSON golden (so any schema/encoding change surfaces
  in the diff).
- This reference IR is the **forward target**: Phase 3's parser test will assert its output
  `isequalIR` this file. We are authoring Phase 3's acceptance criterion now.

### 8.3 Harness note

`tests/run_tests.m` gains an `addpath(<repoRoot>/src)` so `gdsge.ir.*` resolves. This is
consistent with the path policy: the *one* source a new-pipeline test process adds is `src/`;
old-toolbox golden capture remains in its own separate `matlab -batch` process. No
persistent path; the harness owns add/restore.

---

## 9. Done-when

- `pwsh -File tests/run.ps1` exits 0 with all `ir/` fragment tests and the HL1996
  reference-IR tests green (alongside the existing Phase 0 tests).
- `gdsge.ir.validate` passes on the HL1996 reference IR; round-trip `isequalIR` holds; the
  JSON golden matches.
- `docs/ir-schema.md` is regenerated from `schema.m`, committed, and the gendoc test confirms
  no drift.
- `PROGRESS.md` marks Phase 1 done and Phase 2 next; the four brainstorming decisions (§2)
  are logged.

---

## 10. Risks

- **Round-trip quirks beyond the four catalogued.** MATLAB's `jsonencode`/`jsondecode` has
  more corners (e.g. struct field-order, integer vs double, very large arrays). Mitigated by
  the descriptor-driven walker (one place to add a rule) and the HL1996 reference IR exercising
  a real model end-to-end, not just fragments.
- **Schema over-fitting to HL1996.** Designing the full inventory while exercising one model
  risks shapes that don't survive the other five. Mitigated by §5.7 honesty, the
  versioning mechanism (cheap to evolve pre-Phase-4), and grounding the full-inventory shapes
  in the actual gmod survey (transition registry, multi-state grids, ASG coexistence,
  model_init, custom-trans reductions are all taken from real models).
- **gendoc/validator descriptor drift from intent.** The descriptor is both the validator's
  rulebook and the doc's source; a wrong descriptor entry is wrong in two places at once.
  Mitigated because that single-source property is also the goal — fixing the descriptor fixes
  both — and the `tGendoc` no-drift test plus `tValidate` coverage pin the behaviour.
```
