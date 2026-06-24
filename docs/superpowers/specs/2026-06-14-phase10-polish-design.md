# Phase 10 Design — Polish (docs, IR version freeze, perf check, cleanup)

- **Date:** 2026-06-14
- **Status:** Approved (design); pending plan
- **Parent spec:** `2026-06-11-refactor-gdsge-design.md` (§15 "Phase 9 — Polish": *docs, IR
  version freeze, performance check vs old, cleanup* — renumbered Phase 10 in `PROGRESS.md`
  after the corpus phases 7a/7b/9a/9b were split out)
- **Related:** all prior phase specs (this phase documents and freezes the system they built)
- **Owner:** Wenlan Luo

---

## 1. Context

Phase 9b completed the backward-compat corpus: all nine models (HeatonLucas1996, safe_assets,
Mendoza2010, GLSW2020, CaoKS2016, Bianchi2011_asg, Cao2011EZ, CaoNie2016, plus the
CaoKS2016+SIMU_INTERP variant) run green end-to-end through the public API on the autodiff
backend, and six of them on the SymPy backend too (the ASG-path gates — CaoKS2016, its
SIMU_INTERP variant, and Bianchi2011_asg — wait on the deferred Phase 8b). The full suite is
**447/447**. The architecture — `gmod → parser → JSON IR → {MATLAB
codegen, C++ codegen (adept autodiff | SymPy analytic Jacobian)}` — is proven.

Phase 10 is the finishing phase. No new model, no new language feature, no new backend. It
**documents, freezes, measures, and tidies** what exists, so the package is something a
maintainer (or a future agent) can pick up without re-reading nine phase specs, and so the IR
contract is explicitly stable.

The owner approved all four buckets, the fuller scope within each, and the structure below.

### 1.1 What this phase does *not* do

These remain deferred with their own future spec→gate cycles; Phase 10 only records them in the
deferred-error/version docs, it does not implement them:

- **Phase 8b** — SymPy analytic-Jacobian backend for ASG + pchip interpolation (CaoKS2016,
  Bianchi2011_asg under `UseAutoDiff=0`).
- **Phase 7e-full (remainder)** — C++-body `var_tensor` (the `POPN` path) and the
  `IterRslt.var_tensor` result field. The MATLAB-side `var_tensor` subset landed in 9a.

No corpus model exercises either, so neither blocks "the refactor is done."

---

## 2. Goals / non-goals

**Goals**

1. A maintainer-facing **README** + **architecture guide** + **user guide** that make the package
   self-explanatory.
2. The **IR contract is explicitly frozen** for major version 1, with a documented versioning
   policy, an IR changelog, and a test that enforces the major-version gate.
3. A **reusable performance harness** and a **committed parity report** confirming the autodiff
   path does not regress vs the old toolbox, and quantifying SymPy vs autodiff.
4. A **clean tree**: the 97 MB Mendoza golden trimmed, a TODO/dead-code pass done, the deferred-
   error taxonomy consolidated.

**Non-goals**

- No performance *improvement* target (parity is the bar, per parent spec §3); the report
  measures, it does not optimize.
- No pass/fail **perf gate** in the unittest suite (timing on a shared machine is noisy and would
  flake). The harness is a deliberately-run script; the committed report is the artifact.
- No `.gmod` language change, no IR shape change (the freeze is a *declaration of stability*, not
  a new schema version).
- `scratch/` is git-ignored (working-dir only) and left untouched.

---

## 3. Structure & ordering

One spec, one plan, one branch `phase10-polish`, four workstreams committed independently and
executed in this order:

1. **Cleanup** — first, so the tree being documented is the final one.
2. **IR version freeze** — small, self-contained.
3. **Perf harness + report** — produces numbers the docs cite.
4. **Docs** — last, describing the final state and linking the perf report.

Rationale (vs splitting into 10a/10b/10c/10d sub-phases): the buckets are independent enough to
commit separately but small enough — freeze and cleanup especially — that four full
spec→plan→implement cycles would be pure overhead. A single ordered spec matches how the
pre-corpus phases (0–6) ran.

---

## 4. Workstream 1 — Cleanup

### 4.1 Trim the Mendoza2010 golden (97 MB → small)

`tests/Mendoza2010/golden/IterRslt.mat` is 97 MB. `tEndToEndMendoza2010.m` compares only these
fields of the loaded golden `G`:

- `G.Iter`, `G.Metric` (scalars)
- `G.shock_trans`, `G.var_state` (compared at 1e-12 / 1e-15)
- `G.var_policy`, `G.var_aux`, `G.var_interp` (compared at the suite's RelTol/AbsTol)

**Plan:** load the golden, inspect per-field byte sizes to find what dominates, then rewrite the
stored file as a struct holding **only the compared set**, saved `-v7` (zlib-compressed). Whatever
the old toolbox additionally stored in `IterRslt` (full output arrays, `var_others`, duplicated
state, etc.) that the test never reads is dropped. If `var_interp` (the 80×80 spline pp
coefficients) is itself the bulk and *is* compared, it stays — document the residual size and why.

- Update the golden-integrity test (`tests/Mendoza2010/...`) to expect the trimmed field set.
- Re-run `tEndToEndMendoza2010` + the integrity test: green.
- Target: the repo is `git clone`-able without LFS; ideally the golden is < a few MB.

The sibling `IterRslt1.mat` (13 MB, the coarse first-stage WarmUp source — `g1.IterRslt1.Iter`
and `IterRslt1.var_state` are read) gets the same treatment if its compared surface is small;
otherwise left as-is (it is already 7× smaller).

### 4.2 TODO / dead-code pass

Grep `src/`, `templates/`, `pyext/gdsge_sympy/` for `TODO|FIXME|XXX|HACK` and obviously-unreferenced
helpers. Triage each:

- **Fix now** if trivial and in-scope.
- **Record** in `PROGRESS.md` (or the deferred-error doc, §5) if it is a real future item.
- **Delete** if stale/obsolete.

Output: zero unexplained TODO markers; a short list in the plan of what was fixed/recorded/deleted.

### 4.3 Deferred-error taxonomy consolidation

Collect every architecture-honest deferral the codebase raises into one reference table
(`docs/deferred-features.md`, linked from the architecture guide):

| Error ID | Raised by | Meaning / workaround | Tracked in |
|---|---|---|---|
| `gdsge:codegen:genCodeSegmentUnsupported` | codegen | `GenCodeSegment` option | Phase 7c |
| `gdsge:parser:varTensorUnsupported` | parser | non-commented `var_tensor` (body path) | 7e / 7e-full |
| `gdsge:parser:varTensorInBodyUnsupported` | analyzeModel | tensor enters residual | 7e-full |
| `gdsge:codegen:varTensorUnsupported` | codegen invariant | hand-authored tensor IR | 7e |
| (sympy + asg/pchip) | sympy backend | `UseAutoDiff=0` with ASG/pchip | 8b |

Audit for consistent IDs/wording; fix any drift. (Exact rows confirmed during implementation by
grepping the `error('gdsge:...')` sites.)

**Testing for WS1:** full suite stays 447/447; Mendoza gate + integrity green after the trim.

---

## 5. Workstream 2 — IR version freeze

The IR is the contract (`CLAUDE.md`: "The IR is the contract; backends never re-parse gmod").
`gdsge.ir.schema.irVersion` is `1.1.0`. `decode.m` already rejects an incompatible **major**
version (`gdsge:ir:incompatibleVersion`). This workstream makes the stability guarantee explicit
and tested; it does **not** change the schema.

### 5.1 Versioning policy doc — `docs/ir-versioning.md`

- The IR JSON is the stable interface between the MATLAB front-end and the code generators.
- **Major** version (`N.x.x`): incremented only on a breaking shape change; `decode` rejects a
  mismatched major (already implemented). Generated artifacts and goldens are tied to a major.
- **Minor** version (`1.N.x`): additive / backward-compatible enrichment (e.g. 1.0.0 → 1.1.0 added
  `model.regions` + tagged plain/conditional equations in Phase 9b; older 1.0.0 consumers would
  ignore-or-fail-soft, newer consumers read both).
- **Patch**: editorial/no-semantic.
- "Freeze" = the major-1 schema is stable; further corpus work adds minors, never silently breaks
  major 1.

### 5.2 IR changelog

A `## IR changelog` section (in `docs/ir-versioning.md` or appended to `docs/ir-schema.md`):

- **1.0.0** — initial IR (Phases 1–8): sections params/setup/setupNames/variables/model/
  equations/simulate/options/hooks; single unconditional model body.
- **1.1.0** — Phase 9b: `model` becomes a `regions` list (`{condition, ...}`); `equations` become a
  tagged `plain | conditional` list (`eqkinds` registry). Existing single-region models migrate as
  `regions[0]` / condition `''` / `plain` (proven shape-only by
  `tRoundtrip/existingModelsAreSingleUnconditionalRegion`).

### 5.3 Freeze test

`tests/ir/tIrVersionFreeze.m`:

- `decode` of a synthetic IR with `irVersion = '2.0.0'` raises `gdsge:ir:incompatibleVersion`.
- `decode` of a `'1.x.y'` IR (e.g. the current `1.1.0` and a hypothetical `'1.9.0'`) succeeds.
- Asserts `gdsge.ir.schema().irVersion` begins `1.` (a tripwire: bumping to major 2 forces a
  deliberate edit of this test, i.e. a conscious break of the freeze).

**Testing for WS2:** the new freeze test green; no golden regeneration (schema unchanged).

---

## 6. Workstream 3 — Performance harness + committed report

### 6.1 The golden-rule constraint

Old and new share flat public names (`gdsge_codegen`, `iter_<model>`, …). Per `CLAUDE.md`, each
run must `addpath` **exactly one** source in its **own `matlab -batch` process**. The harness is
therefore a *driver* that spawns separate processes — it never loads both toolboxes in one MATLAB.

### 6.2 Harness — `tests/perf/run_perf.ps1` (+ MATLAB workers)

For each model in the subset, the driver runs, each in its own process with a controlled `cd` and
a single source on the path:

- **NEW (autodiff):** `addpath src/` → `gdsge_codegen` (cold: time codegen + MEX compile once) →
  `iter_<model>` (time iter, median of *k* reps) → `simulate_<model>` (time simulate). Records
  `Iter`, `Metric`.
- **NEW (SymPy):** same, with `UseAutoDiff=0` (spline models only; ASG models skip this column).
- **OLD:** `addpath base_package/gdsge/source` → old `gdsge_codegen`/`iter`/`simulate` (or the old
  `gdsge` orchestrator) → time the same stages. Records `Iter`, `Metric`.

The old toolbox is known to build here (Phase 0 risk gate: MSVC 2022). A model worker is a small
`.m` that takes a model name + backend flag, returns a struct of timings as JSON on stdout so the
PowerShell driver can aggregate without sharing a MATLAB session.

**Methodology notes recorded in the report:**

- The **iter wall-time** is the headline metric (the dominant, steady-state cost). Codegen+compile
  is one-time and cache-gated; reported but secondary.
- Median of a few reps; note machine + thread count (`NumThreads`).
- Convergence equality (`Iter`, `Metric`) is reported alongside time, so a "faster" number that
  converged differently is visible.

### 6.3 Subset (representative)

| Model | Why | Backends timed |
|---|---|---|
| HeatonLucas1996 | core spline, the vertical slice | old, new-autodiff, new-sympy |
| safe_assets | MIN reduction + named interp + `Re_n(2)` | old, new-autodiff, new-sympy |
| GLSW2020 | model_init + SIMU_INTERP | old, new-autodiff, new-sympy |
| CaoKS2016 | ASG path | old, new-autodiff |

Skips the slow Mendoza fine-grid re-solve and the two-stage Bianchi driver to keep the harness
tractable; the four span spline/ASG and single/multi-feature. (Subset is a config list at the top
of the harness — trivially extended.)

### 6.4 Report — `docs/perf-report.md`

A committed table: `model × {old iter s, new-autodiff iter s, new-sympy iter s, Iter (old/new),
ratio new/old}`, plus codegen+compile times and a short prose reading. Expected story, to be
confirmed by the numbers:

- **autodiff ≈ old** — the new autodiff backend emits the same adept/CoDoSol/kernel structure, so
  iter-time parity is expected (this is the no-regression bar).
- **SymPy vs autodiff** — analytic Jacobian may differ; the report quantifies it.

The report header notes it is a point-in-time measurement on the owner's machine, regenerable via
`tests/perf/run_perf.ps1`.

**Testing for WS3:** the harness runs to completion and writes the report; correctness is already
guaranteed by the existing end-to-end gates (the harness reuses the same public API). No new
unittest assertion on timing.

---

## 7. Workstream 4 — Docs

Written last, describing the final tree and citing the perf report.

### 7.1 Root `README.md`

- **What this is** — a ground-up refactor of the GDSGE toolbox around an explicit JSON IR;
  drop-in backward compatible (every existing `.gmod` runs; `IterRslt`/`SimuRslt` shapes frozen).
- **Quickstart** — `gdsge('HL1996')` (or the `gdsge_codegen` + `iter_*`/`simulate_*` path);
  pointer to a model under `tests/`.
- **The two backends** — autodiff (default, MATLAB + C MEX, no Python) vs SymPy analytic Jacobian
  (`UseAutoDiff=0`, needs the `uv` Python env); when to use which.
- **Running the tests** — `matlab -batch "cd('tests'); run_tests"` (note: `tests/run.ps1` needs
  pwsh, which is absent here — see the memory note); SymPy tests via `uv run`.
- **Repo layout** — `src/+gdsge/{parser,ir,codegen,runtime}`, `templates/`, `tests/`, `pyext/`,
  `docs/`.
- **Links** — design spec, `docs/architecture.md`, `docs/ir-schema.md`, `docs/perf-report.md`.

### 7.2 Developer / architecture guide — `docs/architecture.md`

Consolidates knowledge currently scattered across the phase specs + `notes-for-agents.md`:

- The pipeline `gmod → parser → IR → {MATLAB, C++ autodiff|SymPy}`, with the module map
  (`gdsge.parser.*`, `gdsge.ir.*`, `gdsge.codegen.*`, `gdsge.runtime.*`).
- **The IR is the contract** — backends never re-parse gmod; link the versioning policy.
- How each backend consumes the IR (the `dataLayout`/`POP*` contract; the SymPy registry +
  pyenv bridge; the region-agnostic per-body emitters from 9b).
- **How to add a model** — the model-body probe-first discipline (9a/9b): probe the `.gmod`,
  enumerate which constructs are already covered, write the gap list, then spec.
- **The MATLAB-path golden rule** — one source per process; why.
- Pointer to `docs/deferred-features.md` (§4.3).

### 7.3 User guide — `docs/user-guide.md`

Self-contained for the package (complements, does not replace, gdsge.com):

- **Authoring a model** — the gmod blocks (`parameters`, `var_state`, `var_shock`,
  `var_policy[/_init]`, `var_aux`, `var_interp`, `var_tensor`, `var_output`, `var_others`,
  `model[/_init]`, `equations`, `simulate`, `pre_*`/`post_*`), declarations (`inbound[/_init]` +
  `adaptive`, `initial`, `shock_num`, `shock_trans`), operators (`GDSGE_EXPECT/MIN/MAX/PROD`,
  `GDSGE_INTERP_VEC` + primed `'`), the `'` future marker and `[N]` shock-indexed arrays, the 9b
  conditional `model(<cond>)` regions + `if/else` equations, and macros
  (`#define/#for/#foreach/#if/#mat/#strcat_comma/include/cinclude`).
- **Options reference** — interpolation (`USE_SPLINE/USE_ASG/USE_PCHIP`, `INTERP_ORDER`,
  `ExtrapOrder`), simulate (`SIMU_RESOLVE/SIMU_INTERP`), solver/iteration (`TolEq`, `MaxIter`,
  `SaveFreq`, `PrintFreq`, `NumThreads`), ASG (`AsgMinLevel`/`AsgMaxLevel`/`AsgThreshold`/…), and
  `UseAutoDiff` to select the backend.
- **Reading results** — `IterRslt` and `SimuRslt` field shapes (frozen, backward compatible).
- A short list of currently-deferred features (link `docs/deferred-features.md`).

Much of the raw inventory already exists in `docs/notes-for-agents.md` §"Full gmod feature
inventory" and the parent spec §13 — the user guide turns it into user-facing prose.

**Testing for WS4:** the `ir-schema.md` no-drift test stays green; docs are prose, so verification
is a build/links sanity pass (no broken internal links, code snippets match the real API), not a
unittest.

---

## 8. Risks

- **Mendoza trim changes the golden bytes.** Mitigated: the trim keeps the exact compared field
  *values*; only un-compared fields are dropped and the file is recompressed. The end-to-end gate
  re-runs to prove the comparison still passes. If `var_interp` alone is large, the file stays
  larger and the report documents it (no LFS, no remote yet).
- **Old-toolbox timing variance.** Mitigated: median of reps, fixed thread count, machine noted;
  the report is explicitly point-in-time, not a gate.
- **Doc drift over time.** Accepted: prose docs are a snapshot; the no-drift-tested `ir-schema.md`
  and the regenerable perf report are the load-bearing pieces.

## 9. Deliverables

1. Trimmed `tests/Mendoza2010/golden/IterRslt.mat` (+ updated integrity test); TODO pass;
   `docs/deferred-features.md`.
2. `docs/ir-versioning.md` (policy + changelog); `tests/ir/tIrVersionFreeze.m`.
3. `tests/perf/run_perf.ps1` + MATLAB worker(s); `docs/perf-report.md`.
4. `README.md`; `docs/architecture.md`; `docs/user-guide.md`.
5. `PROGRESS.md` updated: Phase 10 ☑, changelog entry; the two deferred sub-phases (8b, 7e-full
   remainder) restated as the only open items.

Each deliverable is its own commit on `phase10-polish`; the detailed task breakdown comes from the
writing-plans skill after this spec is approved.
