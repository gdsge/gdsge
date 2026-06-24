# Phase 7e — `var_tensor` Deferred Error — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make a non-commented `var_tensor` declaration fail fast with a clear, identified error instead of today's silent broken-codegen path, and defer the full feature to a future phase.

**Architecture:** Two guards. (1) A parser gate in `parseDeclarations.m` raises `gdsge:parser:varTensorUnsupported` the moment `parseVarDecls` reports a non-empty tensor list — before `evalSetup`, so the diagnostic is purely structural and takes precedence over generic setup errors. (2) A defense-in-depth invariant `gdsge.codegen.assertSupportedIR(ir)`, called at the top of both `generateMatlab` and `generateCxx`, raises `gdsge:codegen:varTensorUnsupported` if a hand-authored IR carries a tensor (the parser normally prevents this). The IR schema is unchanged — `variables.tensor` stays the always-empty forward-compat placeholder.

**Tech Stack:** MATLAB R2025b; `matlab.unittest`; the `+gdsge.parser` / `+gdsge.codegen` packages. No Python, no C++ compile, no golden capture.

**Spec:** `docs/superpowers/specs/2026-06-13-phase7e-var-tensor-error-design.md`

**Branch:** `phase7e-var-tensor-error` (already created).

---

## Conventions used in this plan

- **Run one suite (fast TDD loop)** — prints pass/fail counts:
  ```
  matlab -batch "addpath('src'); addpath(fullfile('src','kernels')); addpath('tests'); r = runtests('tests/parser/tVarTensorUnsupported.m'); fprintf('Passed=%d Failed=%d\n', nnz(~[r.Failed]), nnz([r.Failed]))"
  ```
- **Run the full gate (authoritative):**
  ```
  matlab -batch "cd('tests'); run_tests"
  ```
  Exit code 0 = all pass (per memory: `junit.xml` + exit code are authoritative; `results.tap` accumulates stale lines across runs — ignore it). After the run, check `echo $?` / `$LASTEXITCODE`.
- New test classes are plain `.m` files under `tests/` subfolders; `run_tests` discovers them recursively (`fromFolder(..., 'IncludingSubfolders', true)`).

---

## File Structure

- **Create** `tests/parser/tVarTensorUnsupported.m` — the cohesive suite for this feature: parser gate (positive), no-over-fire (commented), and codegen invariant (both generators).
- **Modify** `src/+gdsge/+parser/parseDeclarations.m` — add the pre-eval gate.
- **Create** `src/+gdsge/+codegen/assertSupportedIR.m` — the shared IR invariant.
- **Modify** `src/+gdsge/+codegen/generateMatlab.m` — call `assertSupportedIR` at the top.
- **Modify** `src/+gdsge/+codegen/generateCxx.m` — call `assertSupportedIR` at the top.
- **Modify** `PROGRESS.md` — mark 7e done (re-scoped), add a deferred entry + changelog.

---

## Task 1: Parser gate — `var_tensor` → `gdsge:parser:varTensorUnsupported`

**Files:**
- Create: `tests/parser/tVarTensorUnsupported.m`
- Modify: `src/+gdsge/+parser/parseDeclarations.m` (insert after line 7)

- [ ] **Step 1: Write the failing test (new suite, parser-gate + no-over-fire cases)**

Create `tests/parser/tVarTensorUnsupported.m`:

```matlab
classdef tVarTensorUnsupported < matlab.unittest.TestCase
    % Phase 7e: var_tensor is deferred — declaring one must raise a clear,
    % identified error (parser gate) without over-firing on commented examples.
    methods (Test)
        function varTensorDeclarationErrors(tc)
            % A non-commented var_tensor must be rejected structurally,
            % before any setup eval (no other declarations needed).
            declText = 'var_tensor r w;';
            tc.verifyError(@() gdsge.parser.parseDeclarations(declText), ...
                'gdsge:parser:varTensorUnsupported');
        end
        function errorMessageNamesTheVariables(tc)
            declText = 'var_tensor budget1 budget2;';
            try
                gdsge.parser.parseDeclarations(declText);
                tc.verifyFail('expected parseDeclarations to error on var_tensor');
            catch err
                tc.verifyEqual(err.identifier, 'gdsge:parser:varTensorUnsupported');
                tc.verifyTrue(contains(err.message, 'budget1'));
                tc.verifyTrue(contains(err.message, 'budget2'));
            end
        end
        function commentedVarTensorDoesNotFire(tc)
            % preprocess strips the % comment, so the gate never sees it and a
            % model with only a commented var_tensor parses with an empty tensor
            % list. Guards against the gate over-firing (e.g. live CaoKS2016).
            raw = sprintf(['%% var_tensor r w budget1 budget2;\n' ...
                'var_state k;\nk = linspace(0,1,3);']);
            out = gdsge.parser.parseDeclarations(gdsge.parser.preprocess(raw));
            tc.verifyEmpty(out.variables.tensor);
        end
    end
end
```

- [ ] **Step 2: Run the suite to verify the gate tests fail**

Run:
```
matlab -batch "addpath('src'); addpath(fullfile('src','kernels')); addpath('tests'); r = runtests('tests/parser/tVarTensorUnsupported.m'); fprintf('Passed=%d Failed=%d\n', nnz(~[r.Failed]), nnz([r.Failed]))"
```
Expected: `varTensorDeclarationErrors` and `errorMessageNamesTheVariables` FAIL (no error is raised today — `parseDeclarations` evals defaults and assembles `variables.tensor = {'r','w'}`). `commentedVarTensorDoesNotFire` PASSES (the comment is already stripped). So: `Passed=1 Failed=2`.

- [ ] **Step 3: Add the gate to `parseDeclarations.m`**

In `src/+gdsge/+parser/parseDeclarations.m`, insert immediately after line 7
(`decls = gdsge.parser.parseVarDecls(declText);`), before the `setupBody` line:

```matlab
% Phase 7e: var_tensor is a deferred variable kind. Reject it here —
% structurally, before evalSetup — so the diagnostic is specific and fires
% ahead of any generic setup error. Full support is a future phase
% (see docs/superpowers/specs/2026-06-13-phase7e-var-tensor-error-design.md).
if ~isempty(decls.tensorNames)
    error('gdsge:parser:varTensorUnsupported', ...
        ['var_tensor is not yet supported by the refactored GDSGE pipeline ' ...
         '(deferred to a future phase). Found tensor variable(s): %s. ' ...
         'Remove the var_tensor declaration and inline the precomputed ' ...
         'state/shock quantity into the model equations.'], ...
        strjoin(decls.tensorNames, ', '));
end
```

- [ ] **Step 4: Run the suite to verify all three pass**

Run:
```
matlab -batch "addpath('src'); addpath(fullfile('src','kernels')); addpath('tests'); r = runtests('tests/parser/tVarTensorUnsupported.m'); fprintf('Passed=%d Failed=%d\n', nnz(~[r.Failed]), nnz([r.Failed]))"
```
Expected: `Passed=3 Failed=0`.

- [ ] **Step 5: Commit**

```bash
git add tests/parser/tVarTensorUnsupported.m src/+gdsge/+parser/parseDeclarations.m
git commit -m "feat(parser): reject var_tensor with gdsge:parser:varTensorUnsupported"
```

---

## Task 2: Codegen invariant — `gdsge.codegen.assertSupportedIR`

**Files:**
- Create: `src/+gdsge/+codegen/assertSupportedIR.m`
- Modify: `src/+gdsge/+codegen/generateMatlab.m` (top of function)
- Modify: `src/+gdsge/+codegen/generateCxx.m` (top of function)
- Modify (add tests): `tests/parser/tVarTensorUnsupported.m`

- [ ] **Step 1: Write the failing tests (codegen invariant, both generators)**

Add a `TestClassSetup` (to put `tests/HeatonLucas1996/ir` on the path so `buildHL1996IR` resolves) and two test methods to `tests/parser/tVarTensorUnsupported.m`. The class becomes:

```matlab
classdef tVarTensorUnsupported < matlab.unittest.TestCase
    % Phase 7e: var_tensor is deferred — declaring one must raise a clear,
    % identified error (parser gate) without over-firing on commented examples;
    % a hand-authored IR carrying a tensor is caught by the codegen invariant.
    methods (TestClassSetup)
        function addIrHelper(tc)
            here = fileparts(mfilename('fullpath'));            % tests/parser
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(here), 'HeatonLucas1996', 'ir')));
        end
    end
    methods (Test)
        function varTensorDeclarationErrors(tc)
            declText = 'var_tensor r w;';
            tc.verifyError(@() gdsge.parser.parseDeclarations(declText), ...
                'gdsge:parser:varTensorUnsupported');
        end
        function errorMessageNamesTheVariables(tc)
            declText = 'var_tensor budget1 budget2;';
            try
                gdsge.parser.parseDeclarations(declText);
                tc.verifyFail('expected parseDeclarations to error on var_tensor');
            catch err
                tc.verifyEqual(err.identifier, 'gdsge:parser:varTensorUnsupported');
                tc.verifyTrue(contains(err.message, 'budget1'));
                tc.verifyTrue(contains(err.message, 'budget2'));
            end
        end
        function commentedVarTensorDoesNotFire(tc)
            raw = sprintf(['%% var_tensor r w budget1 budget2;\n' ...
                'var_state k;\nk = linspace(0,1,3);']);
            out = gdsge.parser.parseDeclarations(gdsge.parser.preprocess(raw));
            tc.verifyEmpty(out.variables.tensor);
        end
        function generateMatlabRejectsTensorIR(tc)
            ir = buildHL1996IR();
            ir.variables.tensor = {'r'};
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            tc.verifyError(@() gdsge.codegen.generateMatlab(ir, work), ...
                'gdsge:codegen:varTensorUnsupported');
        end
        function generateCxxRejectsTensorIR(tc)
            ir = buildHL1996IR();
            ir.variables.tensor = {'r'};
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            tc.verifyError(@() gdsge.codegen.generateCxx(ir, work), ...
                'gdsge:codegen:varTensorUnsupported');
        end
    end
end
```

- [ ] **Step 2: Run the suite to verify the two new tests fail**

Run:
```
matlab -batch "addpath('src'); addpath(fullfile('src','kernels')); addpath('tests'); r = runtests('tests/parser/tVarTensorUnsupported.m'); fprintf('Passed=%d Failed=%d\n', nnz(~[r.Failed]), nnz([r.Failed]))"
```
Expected: the three Task-1 tests PASS; `generateMatlabRejectsTensorIR` and `generateCxxRejectsTensorIR` FAIL (today both generators ignore `variables.tensor` and write files without error). So: `Passed=3 Failed=2`.

- [ ] **Step 3: Create the shared invariant helper**

Create `src/+gdsge/+codegen/assertSupportedIR.m`:

```matlab
function assertSupportedIR(ir)
% ASSERTSUPPORTEDIR  Defense-in-depth invariant: refuse IR features the
%   refactored pipeline does not yet generate code for. Called at the top of
%   the IR-consuming generators (generateMatlab/generateCxx) so a hand-authored
%   or hand-edited IR cannot silently slip an unsupported feature past codegen.
%   Phase 7e: var_tensor is deferred (the parser already rejects it; this
%   guards the direct-generator-call path).
if ~isempty(ir.variables.tensor)
    error('gdsge:codegen:varTensorUnsupported', ...
        ['IR carries var_tensor variable(s) (%s), which the refactored ' ...
         'pipeline does not yet generate code for. var_tensor is deferred ' ...
         'to a future phase.'], strjoin(ir.variables.tensor, ', '));
end
end
```

- [ ] **Step 4: Wire it into `generateMatlab.m`**

In `src/+gdsge/+codegen/generateMatlab.m`, add the call right after the
`if nargin < 2; outDir = pwd; end` line and before `files = struct();`:

```matlab
gdsge.codegen.assertSupportedIR(ir);
```

- [ ] **Step 5: Wire it into `generateCxx.m`**

In `src/+gdsge/+codegen/generateCxx.m`, add the call right after the
`if ~exist(outDir, 'dir'); mkdir(outDir); end` line and before the existing
`% --- refuse what this phase does not support` block:

```matlab
gdsge.codegen.assertSupportedIR(ir);
```

- [ ] **Step 6: Run the suite to verify all five pass**

Run:
```
matlab -batch "addpath('src'); addpath(fullfile('src','kernels')); addpath('tests'); r = runtests('tests/parser/tVarTensorUnsupported.m'); fprintf('Passed=%d Failed=%d\n', nnz(~[r.Failed]), nnz([r.Failed]))"
```
Expected: `Passed=5 Failed=0`.

- [ ] **Step 7: Commit**

```bash
git add src/+gdsge/+codegen/assertSupportedIR.m src/+gdsge/+codegen/generateMatlab.m src/+gdsge/+codegen/generateCxx.m tests/parser/tVarTensorUnsupported.m
git commit -m "feat(codegen): assertSupportedIR invariant rejects var_tensor IR"
```

---

## Task 3: PROGRESS bookkeeping + full-suite gate

**Files:**
- Modify: `PROGRESS.md`

- [ ] **Step 1: Run the full suite (regression gate) BEFORE editing docs**

Run:
```
matlab -batch "cd('tests'); run_tests"
```
Then check the exit code (`echo $?` in bash / `$LASTEXITCODE` in PowerShell).
Expected: exit code `0` — every existing suite (HL1996, safe_assets, Mendoza2010, GLSW, CaoKS2016, CaoKS2016_simu_interp, Bianchi2011, macros, codegen, ir, runtime, …) stays green, plus the 5 new `tVarTensorUnsupported` tests. No corpus model declares a non-commented `var_tensor`, so nothing regresses.

If anything fails, stop and diagnose before proceeding (use superpowers:systematic-debugging).

- [ ] **Step 2: Update the Phase 7e line in `PROGRESS.md`**

In `PROGRESS.md`, replace the current Phase 7e entry:

```markdown
- ☐ **Phase 7e — `var_tensor`** — **NEXT** — (new variable kind across decl → IR → MATLAB+C++ codegen →
  runtime; spline-path only, ⊥ ASG; synthetic differential golden)
```

with:

```markdown
- ☑ **Phase 7e — `var_tensor` deferred error** (done 2026-06-13) — re-scoped (owner-approved):
  instead of building the full variable kind, a non-commented `var_tensor` now raises a clear
  parse-time error (`gdsge:parser:varTensorUnsupported`, gate in `parseDeclarations` before
  `evalSetup`), with a defense-in-depth codegen invariant
  (`gdsge.codegen.assertSupportedIR` → `gdsge:codegen:varTensorUnsupported`, called by both
  generators). Mirrors Phase 7c's `GenCodeSegment` deferral; replaces today's silent
  broken-codegen path. No corpus gmod declares a non-commented tensor (CaoKS2016's is a
  stripped comment), so nothing regresses; zero golden capture. IR schema unchanged
  (`variables.tensor` stays the always-empty placeholder). Verified by `tVarTensorUnsupported`
  (parser gate + no-over-fire + both-generator invariant). Design:
  `docs/superpowers/specs/2026-06-13-phase7e-var-tensor-error-design.md`.
- ☐ **Phase 7e-full (deferred) — `var_tensor` full support** — the real variable kind
  (decl → IR enrichment → MATLAB+C++ codegen → `SIMU_RESOLVE` recompute → `IterRslt.var_tensor`;
  spline-path only, ⊥ ASG; synthetic differential golden). Future-phase sketch in the 7e spec §7.
```

- [ ] **Step 3: Add a changelog entry at the top of the `## Changelog` list in `PROGRESS.md`**

Insert immediately under the `## Changelog` heading, above the existing
`2026-06-13: **Phase 7d complete.**` entry:

```markdown
- 2026-06-13: **Phase 7e complete (re-scoped).** `var_tensor` deferred to a future phase;
  declaring one now fails fast instead of silently generating broken code. Parser gate in
  `parseDeclarations` (after `parseVarDecls`, before `evalSetup`) raises
  `gdsge:parser:varTensorUnsupported`, naming the offending variables and pointing to the
  inline-into-equations workaround — purely structural, so it fires ahead of any setup-eval
  error. Defense-in-depth: a shared `gdsge.codegen.assertSupportedIR(ir)` invariant, called at
  the top of `generateMatlab` and `generateCxx`, raises `gdsge:codegen:varTensorUnsupported`
  if a hand-authored IR carries a tensor (`dataLayout`/`emitDataPack`/`emitPop` only handle
  states). Mirrors the Phase 7c `GenCodeSegment` architecture-honest deferral. No shipped gmod
  declares a non-commented `var_tensor` (CaoKS2016's `% var_tensor ...` is stripped by
  `preprocess` before parsing), so the corpus stays green; IR schema unchanged
  (`variables.tensor` remains the always-empty forward-compat placeholder — no golden
  regeneration). Zero old-toolbox capture. New suite `tVarTensorUnsupported` (5 tests: parser
  gate, message-names-vars, commented-no-over-fire, both-generator invariant). Full suite green.
  All on branch `phase7e-var-tensor-error`. Next: Phase 8 (SymPy analytic-Jacobian backend),
  with full `var_tensor` support tracked as a deferred sub-phase (7e spec §7).
```

- [ ] **Step 4: Commit**

```bash
git add PROGRESS.md
git commit -m "docs(progress): Phase 7e complete (var_tensor deferred error)"
```

---

## Done criteria

- `tVarTensorUnsupported` (5 tests) green: parser gate fires + names variables, commented `var_tensor` does not over-fire, both generators reject a tensor-carrying IR.
- Full suite (`matlab -batch "cd('tests'); run_tests"`) exits `0`.
- `PROGRESS.md` reflects the re-scope + the deferred full-feature entry.
- Branch `phase7e-var-tensor-error` ready to merge (handle via superpowers:finishing-a-development-branch).
