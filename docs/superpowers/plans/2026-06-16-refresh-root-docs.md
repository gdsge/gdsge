# Refresh root `docs/*.md` to current status — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring the three stale root docs (`architecture.md`, `user-guide.md`, `perf-report.md`) back in line with the current `src/` and PROGRESS changelog, including freshly re-measured performance numbers.

**Architecture:** Docs-only change plus one run of the existing perf harness and a small profiling script. No `src/` behavior changes. For prose edits the TDD loop is grep-driven: a grep proving the stale text is present is the "failing test"; a grep proving it is gone and the correction is in place is the "passing test". `ir-schema.md`, `ir-versioning.md`, `deferred-features.md`, `notes-for-agents.md` are verified-current and untouched.

**Tech Stack:** Markdown; MATLAB R2025b (`matlab -batch`) for the perf harness and profiling; `tests/perf/run_perf`.

**Spec:** `docs/superpowers/specs/2026-06-16-refresh-root-docs-design.md`

---

## File structure

- Modify: `docs/architecture.md` — three fact corrections (Task 1)
- Modify: `docs/user-guide.md` — ASG-under-SymPy fix + new backend-selection subsection (Task 2)
- Modify: `docs/perf-report.md` — re-measured headline table + re-profiled breakdown + corrected footnote (Tasks 3–5)
- Create (throwaway): `scratch/profile_breakdown.m` — re-profiling helper (Task 4; `scratch/` is git-ignored)

---

## Task 1: `architecture.md` — three fact corrections

**Files:**
- Modify: `docs/architecture.md`

- [ ] **Step 1: Confirm the stale text is present (failing test)**

Run:
```
rg -n "useFusedConstruct|myppual_mex is retained|in-process via `pyenv`|in-process" docs/architecture.md
```
Expected: matches on the runtime-module line (`useFusedConstruct`), the kernels paragraph (`myppual_mex` retained), and the SymPy bridge "in-process" mentions.

- [ ] **Step 2: Remove `useFusedConstruct` from the runtime module list**

In the `gdsge.runtime` bullet (currently ends `... reportUnconverged, useFusedConstruct (internal A/B switch).`), change to end at `reportUnconverged`:

Old:
```
`ensurePath`, `printIterProgress`, `reportUnconverged`, `useFusedConstruct` (internal A/B switch).
```
New:
```
`ensurePath`, `printIterProgress`, `reportUnconverged`.
```

- [ ] **Step 3: Rewrite the kernels paragraph closing claim**

The `src/kernels/` paragraph currently ends:
```
Compiled
cache-gated by `ensureSplineConstructMex` (flags from the owner's `compile_myppual.m` + `/O2`,
MSVC `/fp:precise`). The prebuilt `myppual_mex` is retained for the simulate / `output_interp`
evaluation path only.
```
Replace the final sentence so it reflects the retirement of `myppual`:
```
Compiled
cache-gated by `ensureSplineConstructMex` (flags from the owner's `compile_myppual.m` + `/O2`,
MSVC `/fp:precise`). The simulate / `output_interp` evaluation path uses a **stacked
uniform-order** layout built by `interp_construct_mex` and evaluated by the generic
`interp_eval_mex` (both on the double-only `include/interp_eval_double.h`); cartesian
`SIMU_INTERP` runs the whole period loop in a generated `simulate_<model>_mex`. `myppual.m`,
`myppual_mex`, and `convert_to_interp_eval_array.m` are retired.
```

- [ ] **Step 4: Correct the SymPy "in-process" claims to out-of-process**

There are two mentions. First, the module-map `+sympy` line:
```
The SymPy bridge
(`+sympy`): `ensurePyenv` (pyenv at the uv venv) + `callSympy` (one JSON in/out to
`pyext/gdsge_sympy`).
```
Change to:
```
The SymPy bridge
(`+sympy`): `ensurePyenv` (an **out-of-process** pyenv at the uv venv —
`ExecutionMode='OutOfProcess'`, to avoid an in-process native crash) + `callSympy` (one JSON
in/out to `pyext/gdsge_sympy`).
```

Second, the backend-consumption bullet:
```
SymPy
  (`pyext/gdsge_sympy`, in-process via `pyenv`) differentiates each body and returns value+partials
  as shared-CSE C++.
```
Change `in-process via `pyenv`` to:
```
SymPy
  (`pyext/gdsge_sympy`, out-of-process via `pyenv`) differentiates each body and returns value+partials
  as shared-CSE C++.
```

- [ ] **Step 5: Confirm the corrections (passing test)**

Run:
```
rg -n "useFusedConstruct|myppual_mex is retained|in-process" docs/architecture.md
```
Expected: **no matches** (the only `myppual` mentions now say "are retired"; the only `pyenv` mentions say "out-of-process").

Run:
```
rg -n "interp_construct_mex|interp_eval_mex|OutOfProcess|out-of-process" docs/architecture.md
```
Expected: matches present (corrections landed).

- [ ] **Step 6: Commit**

```
git add docs/architecture.md
git commit -m "docs(architecture): retire myppual + useFusedConstruct, mark pyenv out-of-process"
```

---

## Task 2: `user-guide.md` — ASG-under-SymPy fix + backend-selection subsection

**Files:**
- Modify: `docs/user-guide.md`

- [ ] **Step 1: Confirm the stale text is present (failing test)**

Run:
```
rg -n "Spline interpolation only today|SymPy backend for ASG/pchip" docs/user-guide.md
```
Expected: matches at the "Selecting the SymPy backend" limits bullet and the "Deferred features" sentence.

- [ ] **Step 2: Correct the ASG-under-SymPy limit bullet**

In "Selecting the SymPy backend", replace the limits bullet:

Old:
```
- **Spline interpolation only** today. ASG or pchip under `UseAutoDiff=0` raises a clear error
  (`gdsge:codegen:sympyInterpUnsupported` / `gdsge:codegen:unsupported`) — see
  `docs/deferred-features.md` (Phase 8b).
```
New:
```
- **Spline and ASG interpolation** are supported. **pchip** under `UseAutoDiff=0` still raises a
  clear error (`gdsge:codegen:sympyInterpUnsupported` / `gdsge:codegen:unsupported`) — see
  `docs/deferred-features.md` (it has no C++ backend at all).
```

- [ ] **Step 3: Trim the deferred-features sentence to pchip**

Old:
```
rather than producing broken code. See `docs/deferred-features.md` for the full map (C++-body
`var_tensor`, the SymPy backend for ASG/pchip, `GenCodeSegment`, …).
```
New:
```
rather than producing broken code. See `docs/deferred-features.md` for the full map (C++-body
`var_tensor`, the SymPy backend for pchip, `GenCodeSegment`, …).
```

- [ ] **Step 4: Rename the section and add the backend-selection subsection**

Rename the heading `## Selecting the SymPy backend` to `## Choosing the C++ Jacobian backend`, and insert the following **before** the existing "Add `UseAutoDiff=0;`…" paragraph (keep that paragraph and the limits bullets that follow it):

```markdown
## Choosing the C++ Jacobian backend

The C++ MEX solver can build its per-grid-point Jacobian two ways: **adept autodiff** (records
a tape) or a **SymPy analytic Jacobian** (CSE-shared symbolic gradient). The backend is decided
at codegen time by this precedence (`gdsge.codegen.resolveBackend`):

1. **In-gmod `UseAutoDiff`** — authoritative when set. `UseAutoDiff=1;` ⇒ adept; `UseAutoDiff=0;`
   ⇒ SymPy.
2. **`GDSGE_BACKEND` environment variable** — `adept`/`autodiff` ⇒ adept, `sympy` ⇒ SymPy,
   `auto` (or unset) ⇒ auto-detect. An unrecognized value errors.
3. **Auto-detect** (nothing pinned): **Python present ⇒ SymPy**, with a codegen-time fallback to
   adept if SymPy generation fails; **Python absent ⇒ adept**.

Two consequences worth knowing:

- With nothing pinned, **installing the `uv` Python env flips the default to SymPy**. To keep a
  Python-free adept build regardless, pin it: `UseAutoDiff=1;` in the gmod, or
  `GDSGE_BACKEND=adept`. (The test harness pins `GDSGE_BACKEND=adept` for exactly this reason.)
- The auto-mode SymPy→adept fallback only fires in **auto** mode. An **explicit** `UseAutoDiff` or
  `GDSGE_BACKEND` pin does **not** fall back — it errors if that backend cannot generate.

The driver prints a one-line `Backend: …` message every run so the resolved choice is visible.
```

- [ ] **Step 5: Point the options-table row at the new subsection**

Old (options-reference table row):
```
| `UseAutoDiff` | `1` = adept autodiff backend; `0` = SymPy analytic Jacobian | `1` |
```
New:
```
| `UseAutoDiff` | `1` = adept autodiff; `0` = SymPy analytic Jacobian; unset = auto-detect (see [Choosing the C++ Jacobian backend](#choosing-the-c-jacobian-backend)) | auto |
```

- [ ] **Step 6: Confirm the corrections (passing test)**

Run:
```
rg -n "Spline interpolation only today|SymPy backend for ASG/pchip" docs/user-guide.md
```
Expected: **no matches**.

Run:
```
rg -n "Choosing the C\+\+ Jacobian backend|GDSGE_BACKEND|Spline and ASG interpolation" docs/user-guide.md
```
Expected: matches present.

- [ ] **Step 7: Commit**

```
git add docs/user-guide.md
git commit -m "docs(user-guide): ASG-under-SymPy supported; document backend auto-selection"
```

---

## Task 3: Re-measure headline perf numbers

**Files:**
- Modify (regenerated): `docs/perf-report.md`

> Note: `run_perf` **overwrites** `docs/perf-report.md` with only the headline table + a one-line "Reading". The full narrative (Findings, breakdown tables, ASG-optimization record, SIMU_INTERP section) is preserved in the spec/plan context and is re-authored in Task 5. Run this task first so the fresh table exists; do not commit the half-file it produces.

- [ ] **Step 1: Run the perf harness (background, sequential — minutes)**

Run:
```
matlab -batch "cd('tests/perf'); run_perf"
```
Expected: console lines `[run_perf] <model> / <mode> ...` for HL1996/safe_assets/GLSW_interp/CaoKS2016 across old/autodiff/sympy, then `[run_perf] wrote docs/perf-report.md (N results)`. A child may print `(child exited … after writing result — likely MEX teardown; result kept)` — that is the known Windows MEX/OpenMP teardown artifact, **not** a failure.

- [ ] **Step 2: Capture the regenerated headline table**

Run:
```
rg -n "^\| " docs/perf-report.md
```
Expected: the freshly written headline table (Model / Iter count / old·AD·SymPy iter s / AD/old / codegen+compile s). Record these numbers — they are the authoritative headline figures used in Task 5. Do not commit yet.

---

## Task 4: Re-profile the breakdown percentages

**Files:**
- Create: `scratch/profile_breakdown.m` (git-ignored throwaway)

- [ ] **Step 1: Write the profiling helper**

Create `scratch/profile_breakdown.m`. It runs one model's `iter` under the MATLAB profiler (one source on the path, controlled `cd`, restored on cleanup — the golden-rule pattern), then prints self-time for the buckets used in the report: MEX solve (`mex_<model>`), the generated `iter_<model>` loop, spline rebuild (`constructSplines` + `interp_construct_mex` + `interp_eval_mex`), `solveProblems`/`solveProblemsAsg`, and `computeMetric`+`printIterProgress`.

```matlab
function profile_breakdown(model, gmodPath, backend)
% profile_breakdown('HL1996','../tests/HeatonLucas1996/HL1996.gmod','adept')
% Run from scratch/ with one source per process. Prints self-time per bucket.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(fullfile(repoRoot, 'src'));
setenv('GDSGE_BACKEND', backend);
tmp = tempname; mkdir(tmp); cd(tmp);
copyfile(gmodPath, pwd);
[~, base, ext] = fileparts(gmodPath);
gdsge_codegen([base ext]);                      % parse -> codegen -> compile
profile clear; profile on;
IterRslt = feval(['iter_' model]); %#ok<NASGU>
profile off;
S = profile('info');
names = {S.FunctionTable.FunctionName};
self  = [S.FunctionTable.TotalTime] - arrayfun(@(f) sum([f.Children.TotalTime]), S.FunctionTable);
get = @(n) sum(self(contains(names, n)));
total = sum(self);
fprintf('--- %s / %s : total self %.2fs ---\n', model, backend, total);
fprintf('MEX solve (mex_%s)      : %5.1f%%\n', model, 100*get(['mex_' model])/total);
fprintf('iter_%s loop            : %5.1f%%\n', model, 100*get(['iter_' model])/total);
fprintf('spline rebuild          : %5.1f%%\n', 100*(get('constructSplines')+get('interp_construct_mex')+get('interp_eval_mex'))/total);
fprintf('solveProblems[Asg]      : %5.1f%%\n', 100*(get('solveProblems'))/total);
fprintf('metric+print            : %5.1f%%\n', 100*(get('computeMetric')+get('printIterProgress'))/total);
end
```

- [ ] **Step 2: Profile the spline models (AD + SymPy) and the ASG model**

Run each in its own process (sequential — one MATLAB at a time):
```
matlab -batch "cd('scratch'); profile_breakdown('HL1996','../tests/HeatonLucas1996/HL1996.gmod','adept')"
matlab -batch "cd('scratch'); profile_breakdown('HL1996','../tests/HeatonLucas1996/HL1996.gmod','sympy')"
matlab -batch "cd('scratch'); profile_breakdown('GLSW_interp','../tests/GLSW2020/GLSW_interp.gmod','adept')"
matlab -batch "cd('scratch'); profile_breakdown('GLSW_interp','../tests/GLSW2020/GLSW_interp.gmod','sympy')"
matlab -batch "cd('scratch'); profile_breakdown('mendoza2010','../tests/Mendoza2010/mendoza2010.gmod','sympy')"
matlab -batch "cd('scratch'); profile_breakdown('CaoKS2016','../tests/CaoKS2016/CaoKS2016.gmod','adept')"
matlab -batch "cd('scratch'); profile_breakdown('CaoKS2016','../tests/CaoKS2016/CaoKS2016.gmod','sympy')"
```
Expected: each prints a per-bucket percentage block. Record the percentages — they replace the spline-breakdown table (lines ~36–44) and confirm or update the ASG-breakdown table (lines ~62–65). If the ASG percentages match the existing table within a few points, keep the existing text; otherwise update.

> If a model's gmod filename differs (e.g. exact case), confirm with `ls tests/<Model>/*.gmod` first. Mendoza is the largest golden but profiling runs codegen+iter fresh, independent of the golden `.mat`.

---

## Task 5: Re-author `perf-report.md` with fresh numbers + corrected footnote

**Files:**
- Modify: `docs/perf-report.md`

- [ ] **Step 1: Rewrite the file**

Reconstruct `docs/perf-report.md` from the fresh measurements, preserving the existing narrative structure. Concretely:

1. **Headline table + "Reading" line** — use the Task 3 numbers (the file already holds the regenerated table from `run_perf`; restore the richer "Reading" line that notes the `safe_assets` RNG sensitivity).
2. **Findings narrative** — keep the structure, but update any figure that changed. In particular fix the known internal inconsistency: the CaoKS2016 AD/old ratio in the narrative must match the headline table (the post-optimization value, ~0.43, not the stale 0.94).
3. **Spline profiling-breakdown table (~lines 36–44)** — replace the percentages with the Task 4 numbers and **correct the footnote** from
   `\* `constructSplines` + `myppual` + `myppual_mex` + `convert_to_interp_eval_array``
   to
   `\* `constructSplines` + `interp_construct_mex` + `interp_eval_mex``.
4. **ASG breakdown tables and the ASG-optimization record (~lines 58–133)** — keep as-is unless Task 4 shows a material drift; this section is a dated record of the 2026-06-15 optimization (it legitimately reports historical before/after numbers).
5. **`## SIMU_INTERP whole-loop MEX (2026-06-16)` section (~lines 143–158)** — keep verbatim; it is current.

- [ ] **Step 2: Confirm the corrections (passing test)**

Run:
```
rg -n "myppual" docs/perf-report.md
```
Expected: **no matches** (the spline-rebuild footnote no longer names retired kernels).

Run:
```
rg -n "interp_construct_mex|interp_eval_mex" docs/perf-report.md
```
Expected: present in the corrected footnote.

- [ ] **Step 3: Commit**

```
git add docs/perf-report.md
git commit -m "docs(perf): re-measure on current code; correct retired-kernel footnote"
```

---

## Task 6: Final verification

**Files:** none (verification only)

- [ ] **Step 1: Cross-file stale-term sweep**

Run:
```
rg -n "myppual|useFusedConstruct|in-process" docs/architecture.md docs/user-guide.md docs/perf-report.md
```
Expected: **no matches**, except deliberate "are retired" phrasing in `architecture.md`.

- [ ] **Step 2: Run the doc/IR tests to confirm nothing regressed**

Run:
```
matlab -batch "cd('tests'); runtests('ir')"
```
Expected: all green (notably `tGendoc` — confirms `ir-schema.md` is still drift-free and was correctly left untouched). No `src/` code changed, so the broader suite is unaffected.

- [ ] **Step 3: Confirm the untouched docs are unchanged**

Run:
```
git status --porcelain docs/
```
Expected: only `architecture.md`, `user-guide.md`, `perf-report.md` modified (already committed); `ir-schema.md`, `ir-versioning.md`, `deferred-features.md`, `notes-for-agents.md` untouched.

---

## Self-review notes

- **Spec coverage:** Task 1 ⇄ architecture.md (3 corrections); Task 2 ⇄ user-guide.md (ASG fix + backend subsection); Tasks 3–5 ⇄ perf-report.md (re-measure + re-profile + footnote); Task 6 ⇄ verification incl. the untouched-files check. All spec sections covered.
- **Measured-value note:** Tasks 3–5 reference numbers produced by Steps in Task 3/4 rather than hard-coding them — these are measurements, not placeholders; the commands that produce them are exact.
- **Naming consistency:** `interp_construct_mex` / `interp_eval_mex` / `resolveBackend` / `ensurePyenv` used consistently across tasks and match `src/`.
