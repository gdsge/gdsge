# Phase 0 — Infrastructure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the repository scaffolding, the command-line MATLAB test harness, the uv-managed Python environment, a tested golden-comparison utility, and captured backward-compatibility golden results for HeatonLucas1996 — so that Phases 1–9 can be implemented test-first by future agents.

**Architecture:** New package at the repo root (`src/+gdsge/...`), `base_package/` git-ignored. Tests use native `matlab.unittest` run headless via `matlab -batch`. The toolbox source (old or new) is never persistently on the MATLAB path — each run `addpath`s exactly one source in its own process. Golden results are captured from the *old* toolbox and compared with a tolerance-based utility.

**Tech Stack:** MATLAB R2025b (`matlab.unittest`, `matlab -batch`), PowerShell wrapper, `uv` + Python 3.12 + SymPy (alt backend only), git.

---

## File structure (created/modified in this phase)

- `src/+gdsge/` — package skeleton: `+parser/`, `+ir/`, `+codegen/`, `+runtime/` (placeholders only this phase)
- `templates/matlab/`, `templates/cxx/` — empty template homes (placeholders)
- `pyext/` — `pyproject.toml`, `tests/test_smoke.py` (uv project; `.venv`/`.python-version` git-ignored)
- `tests/run_tests.m` — builds + runs the unittest suite, writes JUnit/TAP, sets exit code
- `tests/run.ps1` — finds R2025b, runs `run_tests` headless, propagates exit code
- `tests/smoke/tSmoke.m` — trivial harness smoke test
- `tests/+gdsgetest/compareNumericClose.m` — recursive tolerance comparison of structs (package function)
- `tests/utils/tCompareNumericClose.m` — unit tests for the comparison utility (plain folder)
- `tests/golden/capture_HL1996.m` — captures golden IterRslt/SimuRslt from the OLD toolbox
- `tests/HeatonLucas1996/golden/IterRslt.mat`, `SimuRslt.mat` — captured goldens (tracked)
- `tests/HeatonLucas1996/tGoldenHL1996.m` — golden-integrity test now; extended to compare new pipeline in Phase 6
- `CLAUDE.md` — project orientation for Claude sessions
- `PROGRESS.md` — living progress log / broad plan tracker
- `docs/notes-for-agents.md` — orientation + architecture facts discovered during exploration
- `docs/ir-schema.md` — IR schema stub

Already done before this plan (committed): `git init` on `main`, `.gitignore`, the design spec under `docs/superpowers/specs/`.

---

## Task 1: Package skeleton

**Files:**
- Create: `src/+gdsge/Contents.m`
- Create: `src/+gdsge/+parser/.gitkeep`, `src/+gdsge/+ir/.gitkeep`, `src/+gdsge/+codegen/.gitkeep`, `src/+gdsge/+runtime/.gitkeep`
- Create: `templates/matlab/.gitkeep`, `templates/cxx/.gitkeep`

- [ ] **Step 1: Create the package directories and placeholders**

`src/+gdsge/Contents.m`:
```matlab
% GDSGE  Refactored Global DSGE toolbox (package namespace for internal modules).
%
% Subpackages:
%   +parser   - gmod -> IR front end (preprocess, lex, declarations, model expressions)
%   +ir       - IR struct definition, JSON encode/decode, validators
%   +codegen  - IR -> MATLAB and IR -> C++ generators
%   +runtime  - helpers used by generated MATLAB (error reporting, printing)
%
% Public flat entry points (gdsge.m, gdsge_codegen.m) are thin shims added later;
% they are kept at top level for backward compatibility.
```

Create empty `.gitkeep` files in each subpackage and template folder so git tracks the empty dirs.

- [ ] **Step 2: Verify the tree exists**

Run: `Get-ChildItem -Recurse -Directory src, templates | Select-Object FullName`
Expected: lists `src\+gdsge`, the four subpackages, `templates\matlab`, `templates\cxx`.

- [ ] **Step 3: Commit**

```bash
git add src templates
git commit -m "chore: scaffold package skeleton (+gdsge subpackages, templates)"
```

---

## Task 2: CLAUDE.md

**Files:**
- Create: `CLAUDE.md`

- [ ] **Step 1: Write CLAUDE.md**

Content must include these sections (write real prose, not headers only):
- **What this is:** refactor of the GDSGE MATLAB toolbox; new package at repo root; `base_package/` is the old reference (git-ignored).
- **Read first:** `docs/superpowers/specs/2026-06-11-refactor-gdsge-design.md` (the design), `PROGRESS.md` (current phase), `docs/notes-for-agents.md` (codebase facts).
- **Golden rule — MATLAB path:** never add toolbox source to a persistent path; each run `addpath`s exactly one source (old `base_package/gdsge/source` or new `src/`) in its own `matlab -batch` process. Public flat names (`gdsge_codegen`, `iter_<model>`) exist in both old and new — isolate by process + cwd.
- **How to run tests:** `pwsh tests/run.ps1` (or `matlab -batch "cd('tests'); run_tests"`). Exit 0 = pass; report at `tests/results/junit.xml`.
- **How to run Python (alt backend only):** `uv run --project pyext pytest pyext/tests`.
- **Constraints:** default path needs no Python; keep C++ stack-allocated; backward compatible (all gmod run); TDD.
- **MATLAB:** R2025b at `C:\Program Files\MATLAB\R2025b\bin\matlab.exe`.

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: add CLAUDE.md project orientation"
```

---

## Task 3: PROGRESS.md

**Files:**
- Create: `PROGRESS.md`

- [ ] **Step 1: Write PROGRESS.md as the living phase tracker**

Include: a one-paragraph project summary; a status legend; and a checklist mirroring spec Section 15 (Phases 0–9) with sub-bullets for Phase 0 tasks. Mark Phase 0 items as they complete. Include a "Decisions log" section capturing the four brainstorming decisions (native MATLAB unittest; vertical slice = HeatonLucas1996; autodiff backend first; spec + plan + scaffold this session) and a "Risks/Open" section (old toolbox must build to capture goldens; SymPy reduction fusion deferred to Phase 8).

- [ ] **Step 2: Commit**

```bash
git add PROGRESS.md
git commit -m "docs: add PROGRESS.md living phase tracker"
```

---

## Task 4: notes-for-agents.md and ir-schema.md stub

**Files:**
- Create: `docs/notes-for-agents.md`
- Create: `docs/ir-schema.md`

- [ ] **Step 1: Write docs/notes-for-agents.md**

Capture the concrete facts discovered during exploration so future agents need not re-derive them:
- **Old toolbox map:** `base_package/gdsge/source/gdsge_parser.m` (3027-line monolith; macro→block→declaration→model→template stages), `gdsge_codegen.m` (file I/O + compile trigger), `gdsge.m` (orchestrates codegen→iter→simulate). Templates in `source/code_template/`. C++ uses `adept` autodiff; solver is `CoDoSol::solve`; per-grid-point OpenMP parallel; data popped to stack locals via `POPN`/`POPNARRAY`.
- **Reference for analytic Jacobian:** `base_package/hans/hans_parser.py` (+ `template/vfi_template_sympy.h`): SymPy + `cse()` → stack-allocated `double helper_i;` via `ccode()`. Uses Lark grammar (`hans/lark/*.lark`). EXPECT treated as identity there — GDSGE needs true loop fusion (spec Section 10).
- **gmod = MATLAB superset:** grid/param lines (`k = exp(linspace(...))`, `shock_trans = ...`) are `eval`'d by MATLAB; this is why the parser must be MATLAB.
- **Test models:** `base_package/gdsge/tests/{HeatonLucas1996,Bianchi2011_asg,Mendoza2010,CaoKS2016,GLSW2020,Barro_et_al_2017}`; each has `<model>.gmod` + `test.m`. `tests/runtests.m` (old) drives them. HeatonLucas1996 = first vertical slice (1 state w1, 8 shocks, array policy `w1n[8]`, `GDSGE_INTERP_VEC'`, `GDSGE_EXPECT`, future-indexed equation `w1_consis'`, adaptive bounds, KKT complementarity; converges ~iter 209; `rng(0823)` before simulate).
- **gmod feature inventory:** copy the list from spec Section 13 so agents know the full backward-compat surface.
- **How to drive MATLAB headless:** `matlab -batch "<expr>"`; exit code reflects errors / explicit `exit(n)`.

- [ ] **Step 2: Write docs/ir-schema.md stub**

A short document stating the IR is versioned JSON, listing the planned top-level sections from spec Section 8 (`irVersion`, `modelName`, `options`, `shocks`, `states`, `variables`, `bounds`, `interp`, `model`, `modelInit`, `simulate`, hook blocks) and the expression-AST node set (number, name, primed-name, unary/binary op, call, index). Mark it "Stub — finalized in Phase 1."

- [ ] **Step 3: Commit**

```bash
git add docs/notes-for-agents.md docs/ir-schema.md
git commit -m "docs: add notes-for-agents and IR schema stub"
```

---

## Task 5: MATLAB test harness (smoke)

**Files:**
- Create: `tests/run_tests.m`
- Create: `tests/run.ps1`
- Create: `tests/smoke/tSmoke.m`

- [ ] **Step 1: Write the smoke test (the failing test stands in for "harness works")**

`tests/smoke/tSmoke.m`:
```matlab
classdef tSmoke < matlab.unittest.TestCase
    methods (Test)
        function harnessCollectsResults(testCase)
            testCase.verifyTrue(true, 'Smoke test should pass.');
        end
    end
end
```

- [ ] **Step 2: Write the runner**

`tests/run_tests.m`:
```matlab
function run_tests()
% Build and run the MATLAB unittest suite headless; write JUnit + TAP; set exit code.
% Usage: matlab -batch "cd('tests'); run_tests"
import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.XMLPlugin
import matlab.unittest.plugins.TAPPlugin
import matlab.unittest.plugins.ToFile

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);   % so the +gdsgetest package (utility fns) resolves inside tests
resultsDir = fullfile(thisDir, 'results');
if ~exist(resultsDir, 'dir'); mkdir(resultsDir); end

% NOTE: fromFolder does NOT descend into namespace (+pkg) folders, so test
% CLASSES live in plain folders; only utility FUNCTIONS live in +gdsgetest.
suite = TestSuite.fromFolder(thisDir, 'IncludingSubfolders', true);

runner = TestRunner.withTextOutput();
runner.addPlugin(XMLPlugin.producingJUnitFormat(fullfile(resultsDir, 'junit.xml')));
runner.addPlugin(TAPPlugin.producingVersion13(ToFile(fullfile(resultsDir, 'results.tap'))));

results = runner.run(suite);

failed = ~isempty(results) && any([results.Failed]);
exit(double(failed));
end
```

`tests/run.ps1`:
```powershell
# Run the MATLAB unittest suite headless and propagate the exit code.
$ErrorActionPreference = 'Stop'
$matlab = 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe'
if (-not (Test-Path $matlab)) { throw "MATLAB not found at $matlab" }
$testsDir = $PSScriptRoot
& $matlab -batch "cd('$testsDir'); run_tests"
exit $LASTEXITCODE
```

- [ ] **Step 3: Run the harness and verify it passes**

Run: `pwsh -File tests/run.ps1`
Expected: MATLAB launches, runs `tSmoke`, prints "1 Passed", `tests/results/junit.xml` is written, and `$LASTEXITCODE` is 0.

- [ ] **Step 4: Verify failure path sets nonzero exit (temporary check)**

Temporarily change `verifyTrue(true,...)` to `verifyTrue(false,...)`, run `pwsh -File tests/run.ps1`, confirm exit code is 1 and junit.xml records a failure. Then revert to `true`.

- [ ] **Step 5: Commit**

```bash
git add tests/run_tests.m tests/run.ps1 tests/smoke/tSmoke.m
git commit -m "test: add headless MATLAB unittest harness with smoke test"
```

---

## Task 6: Golden comparison utility (real TDD)

**Files:**
- Create: `tests/+gdsgetest/compareNumericClose.m` (utility function — lives in the package)
- Test: `tests/utils/tCompareNumericClose.m` (test class — plain folder so `fromFolder` finds it)

- [ ] **Step 1: Write the failing tests**

`tests/utils/tCompareNumericClose.m`:
```matlab
classdef tCompareNumericClose < matlab.unittest.TestCase
    methods (Test)
        function identicalStructsPass(testCase)
            a = struct('x', [1 2 3], 'y', 4.0);
            r = gdsgetest.compareNumericClose(a, a, 1e-9, 1e-12);
            testCase.verifyTrue(r.pass);
            testCase.verifyEmpty(r.failures);
        end

        function withinTolerancePasses(testCase)
            a = struct('x', 1.0);
            b = struct('x', 1.0 + 1e-10);
            r = gdsgetest.compareNumericClose(a, b, 1e-8, 1e-12);
            testCase.verifyTrue(r.pass);
        end

        function outsideToleranceFails(testCase)
            a = struct('x', 1.0);
            b = struct('x', 1.1);
            r = gdsgetest.compareNumericClose(a, b, 1e-8, 1e-12);
            testCase.verifyFalse(r.pass);
            testCase.verifyNotEmpty(r.failures);
        end

        function nestedStructsCompared(testCase)
            a = struct('inner', struct('v', [1 2; 3 4]));
            b = struct('inner', struct('v', [1 2; 3 4] + 1e-11));
            r = gdsgetest.compareNumericClose(a, b, 1e-8, 1e-12);
            testCase.verifyTrue(r.pass);
        end

        function missingFieldFails(testCase)
            a = struct('x', 1.0, 'y', 2.0);
            b = struct('x', 1.0);
            r = gdsgetest.compareNumericClose(a, b, 1e-8, 1e-12);
            testCase.verifyFalse(r.pass);
        end

        function sizeMismatchFails(testCase)
            a = struct('x', [1 2 3]);
            b = struct('x', [1 2]);
            r = gdsgetest.compareNumericClose(a, b, 1e-8, 1e-12);
            testCase.verifyFalse(r.pass);
        end
    end
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `matlab -batch "cd('tests'); addpath(pwd); results = runtests(fullfile('utils','tCompareNumericClose.m')); disp(table(results))"`
Expected: errors/fails — `gdsgetest.compareNumericClose` not found.

- [ ] **Step 3: Implement the utility**

`tests/+gdsgetest/compareNumericClose.m`:
```matlab
function report = compareNumericClose(actual, expected, relTol, absTol)
% Recursively compare numeric leaves of two structs within tolerance.
% report.pass (logical) and report.failures (cellstr of messages).
if nargin < 3 || isempty(relTol); relTol = 1e-6; end
if nargin < 4 || isempty(absTol); absTol = 1e-9; end
failures = {};
failures = compareNode(actual, expected, '', relTol, absTol, failures);
report = struct('pass', isempty(failures), 'failures', {failures});
end

function failures = compareNode(a, b, path, relTol, absTol, failures)
if isstruct(a) && isstruct(b)
    fa = sort(fieldnames(a)); fb = sort(fieldnames(b));
    if ~isequal(fa, fb)
        failures{end+1} = sprintf('%s: field sets differ {%s} vs {%s}', ...
            iff(isempty(path),'<root>',path), strjoin(fa',','), strjoin(fb',','));
        return;
    end
    for i = 1:numel(fa)
        f = fa{i};
        failures = compareNode(a.(f), b.(f), [path '.' f], relTol, absTol, failures);
    end
elseif isnumeric(a) && isnumeric(b)
    if ~isequal(size(a), size(b))
        failures{end+1} = sprintf('%s: size %s vs %s', path, mat2str(size(a)), mat2str(size(b)));
        return;
    end
    d = abs(double(a(:)) - double(b(:)));
    tol = absTol + relTol .* abs(double(b(:)));
    bad = d > tol;
    if any(bad)
        failures{end+1} = sprintf('%s: max abs diff %g exceeds tol (relTol=%g, absTol=%g)', ...
            path, max(d), relTol, absTol);
    end
else
    if ~isequaln(a, b)
        failures{end+1} = sprintf('%s: non-numeric leaves differ', path);
    end
end
end

function out = iff(cond, a, b)
if cond; out = a; else; out = b; end
end
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `matlab -batch "cd('tests'); addpath(pwd); results = runtests(fullfile('utils','tCompareNumericClose.m')); assert(all(~[results.Failed]))"`
Expected: all pass, exit 0.

- [ ] **Step 5: Commit**

```bash
git add tests/+gdsgetest/compareNumericClose.m tests/utils/tCompareNumericClose.m
git commit -m "test: add tolerance-based struct comparison utility (TDD)"
```

---

## Task 7: uv Python environment (alt SymPy backend)

**Files:**
- Create: `pyext/pyproject.toml`
- Create: `pyext/tests/test_smoke.py`
- Create: `pyext/README.md`

- [ ] **Step 1: Create the uv project files**

`pyext/pyproject.toml`:
```toml
[project]
name = "gdsge-sympy"
version = "0.0.0"
description = "SymPy analytic-Jacobian code generation for the GDSGE alt backend"
requires-python = ">=3.12"
dependencies = ["sympy>=1.12"]

[dependency-groups]
dev = ["pytest>=8"]

[tool.uv]
package = false
```

`pyext/tests/test_smoke.py`:
```python
def test_sympy_diff_and_ccode():
    import sympy
    x = sympy.symbols("x", real=True)
    assert sympy.diff(x**2, x) == 2 * x
    # ccode is the codegen primitive the alt backend depends on
    assert "x" in sympy.ccode(x * x)
```

`pyext/README.md`: one paragraph — this env is only for the alternative SymPy analytic-Jacobian backend (Phase 8); the default toolbox path needs no Python. Run tests with `uv run --project pyext pytest pyext/tests`.

- [ ] **Step 2: Pin Python and sync**

Run:
```
uv python pin 3.12
uv sync --project pyext
```
Expected: creates `pyext/.venv` (git-ignored) and a `uv.lock`.

- [ ] **Step 3: Run the smoke test to verify it passes**

Run: `uv run --project pyext pytest pyext/tests -q`
Expected: 1 passed.

- [ ] **Step 4: Commit**

```bash
git add pyext/pyproject.toml pyext/tests/test_smoke.py pyext/README.md pyext/uv.lock
git commit -m "build: add uv-managed Python env with SymPy for alt backend"
```

---

## Task 8: Old-toolbox build verification + HL1996 golden capture (RISK GATE)

**Files:**
- Create: `tests/golden/capture_HL1996.m`
- Create (output, tracked): `tests/HeatonLucas1996/golden/IterRslt.mat`, `tests/HeatonLucas1996/golden/SimuRslt.mat`
- Create: `tests/HeatonLucas1996/HL1996.gmod` (copied from base_package so the model travels with the new tests)

- [ ] **Step 1: Copy the HL1996 model into the new tests tree**

Run: `Copy-Item base_package/gdsge/tests/HeatonLucas1996/HL1996.gmod tests/HeatonLucas1996/HL1996.gmod`
(The model dir already holds `tGoldenHL1996.m` after Task 9; create the dir now if needed.)

- [ ] **Step 2: Write the golden-capture script (adds OLD source explicitly, runs in a temp dir)**

`tests/golden/capture_HL1996.m`:
```matlab
function capture_HL1996()
% Capture golden IterRslt/SimuRslt for HL1996 from the OLD toolbox.
% Explicitly adds ONLY the old source to the path; runs in a temp working dir.
here     = fileparts(mfilename('fullpath'));          % tests/golden
repoRoot = fileparts(fileparts(here));                % repo root
oldSrc   = fullfile(repoRoot, 'base_package', 'gdsge', 'source');
modelDir = fullfile(repoRoot, 'tests', 'HeatonLucas1996');
goldenDir = fullfile(modelDir, 'golden');
if ~exist(goldenDir, 'dir'); mkdir(goldenDir); end

work = tempname; mkdir(work);
copyfile(fullfile(modelDir, 'HL1996.gmod'), work);

oldPath = path; restore = onCleanup(@() path(oldPath));
addpath(oldSrc);
oldCd = pwd; cdRestore = onCleanup(@() cd(oldCd));
cd(work);

gdsge_codegen('HL1996');                 % compiles the OLD mex — this is the build gate
options = struct; options.SaveFreq = inf;
IterRslt = iter_HL1996(options);
rng(0823);
SimuRslt = simulate_HL1996(IterRslt);

save(fullfile(goldenDir, 'IterRslt.mat'), 'IterRslt', '-v7');
save(fullfile(goldenDir, 'SimuRslt.mat'), 'SimuRslt', '-v7');
fprintf('Golden captured: Iter=%d Metric=%g\n', IterRslt.Iter, IterRslt.Metric);
end
```

- [ ] **Step 3: Run capture and verify the old toolbox builds (RISK GATE)**

Run: `matlab -batch "addpath('tests/golden'); capture_HL1996"`
Expected: the old MEX compiles, iteration converges (≈209 iters), and `tests/HeatonLucas1996/golden/IterRslt.mat` + `SimuRslt.mat` are written.
**If the MEX fails to compile:** STOP. Report the compiler error; do not proceed. Resolution is a user decision (configure the MATLAB C++ compiler, or the user supplies golden `.mat` files).

- [ ] **Step 4: Commit goldens + script**

```bash
git add tests/golden/capture_HL1996.m tests/HeatonLucas1996/HL1996.gmod tests/HeatonLucas1996/golden/IterRslt.mat tests/HeatonLucas1996/golden/SimuRslt.mat
git commit -m "test: capture HL1996 golden results from old toolbox"
```

---

## Task 9: HL1996 golden-integrity test (extended in Phase 6)

**Files:**
- Create: `tests/HeatonLucas1996/tGoldenHL1996.m`

- [ ] **Step 1: Write the golden-integrity test**

`tests/HeatonLucas1996/tGoldenHL1996.m`:
```matlab
classdef tGoldenHL1996 < matlab.unittest.TestCase
    % Phase 0: verify golden files exist, load, and have the expected shape.
    % Phase 6: extend to run the NEW pipeline and compare with
    %          gdsgetest.compareNumericClose against these goldens.
    methods (Test)
        function goldenLoadsAndHasShape(testCase)
            here = fileparts(mfilename('fullpath'));
            g = load(fullfile(here, 'golden', 'IterRslt.mat'));
            testCase.verifyTrue(isfield(g, 'IterRslt'));
            R = g.IterRslt;
            testCase.verifyEqual(R.shock_num, 8);
            testCase.verifyTrue(isfield(R, 'var_policy'));
            testCase.verifyTrue(isfield(R.var_policy, 'c1'));
            testCase.verifyTrue(isfield(R, 'var_state'));
            testCase.verifyEqual(numel(R.var_state.w1), 201);
        end

        function selfComparisonPasses(testCase)
            here = fileparts(mfilename('fullpath'));
            g = load(fullfile(here, 'golden', 'IterRslt.mat'));
            r = gdsgetest.compareNumericClose(g.IterRslt.var_policy, ...
                                              g.IterRslt.var_policy, 1e-9, 1e-12);
            testCase.verifyTrue(r.pass);
        end
    end
end
```

- [ ] **Step 2: Run the full suite headless to confirm everything is green**

Run: `pwsh -File tests/run.ps1`
Expected: `tSmoke`, `tCompareNumericClose` (6 tests), and `tGoldenHL1996` (2 tests) all pass; exit 0; `tests/results/junit.xml` written.

- [ ] **Step 3: Commit**

```bash
git add tests/HeatonLucas1996/tGoldenHL1996.m
git commit -m "test: add HL1996 golden-integrity test (extended in Phase 6)"
```

---

## Phase 0 done-when

- `pwsh -File tests/run.ps1` exits 0 with all tests passing and writes `tests/results/junit.xml`.
- `uv run --project pyext pytest pyext/tests` passes.
- `tests/HeatonLucas1996/golden/IterRslt.mat` + `SimuRslt.mat` exist and load.
- `CLAUDE.md`, `PROGRESS.md`, `docs/notes-for-agents.md`, `docs/ir-schema.md` exist and are committed.
- `base_package/` and refactoring scratch remain git-ignored.
- PROGRESS.md marks Phase 0 complete and Phase 1 as next.

## Self-review notes

- **Spec coverage (Section 16 deliverables):** git/gitignore (done pre-plan) ✓; skeleton (Task 1) ✓; uv+SymPy (Task 7) ✓; CLAUDE/PROGRESS/notes/IR-stub (Tasks 2–4) ✓; harness+run.ps1 (Task 5) ✓; old-build verify + golden capture (Task 8) ✓; commit (every task) ✓.
- **Path policy (spec Section 12):** capture script `addpath`s only the old source in its own process and restores on cleanup ✓.
- **No placeholders:** all code files have complete content; doc tasks list required sections explicitly.
- **Type consistency:** `gdsgetest.compareNumericClose(actual, expected, relTol, absTol)` returning `report.pass`/`report.failures` is used identically in Tasks 6 and 9.
