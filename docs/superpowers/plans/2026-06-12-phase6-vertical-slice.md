# Phase 6 — Vertical Slice Green: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** HL1996 runs end-to-end through the public API — `gdsge_codegen('HL1996')` → `iter_HL1996` → `simulate_HL1996` — and matches the committed goldens.

**Architecture:** A flat backward-compat shim `src/gdsge_codegen.m` forwards to a namespaced driver `gdsge.codegen.codegen` that parses the gmod from cwd, writes the IR contract artifact `<model>.gdsge.json`, calls the existing Phase 4/5 generators into cwd, and compiles the MEX gated by an old-style text cache (`mex_<model>.cache`). A new slow gate test drives everything through the shim and replaces the Phase 5 gate.

**Tech Stack:** MATLAB R2025b (`matlab -batch`, native `matlab.unittest`), existing `gdsge.parser` / `gdsge.ir` / `gdsge.codegen` / `gdsge.runtime` packages, MSVC for MEX.

**Spec:** `docs/superpowers/specs/2026-06-12-phase6-vertical-slice-design.md`

---

## Conventions for every task

- Work from the repo root `D:\refactor_gdsge`.
- `matlab` means `C:\Program Files\MATLAB\R2025b\bin\matlab.exe`. In PowerShell invoke it as
  `& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "<expr>"`.
  (PowerShell 7 / `pwsh` is **not** installed — never use `tests/run.ps1`.)
- Single-test run pattern (replicates `run_tests`'s path setup; exit code 0 = pass):

  ```
  matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','<folder>','<TestFile>.m')); disp(table(results)); exit(any([results.Failed]))"
  ```

- Full suite: `matlab -batch "cd('tests'); run_tests"` — exit 0 = all pass;
  `tests/results/junit.xml` is the authoritative report (ignore `results.tap`, it appends).
- The MATLAB path policy is sacred: never `addpath` persistently; the commands above add
  paths only inside their own `matlab -batch` process.

## File structure

| File | Action | Responsibility |
|---|---|---|
| `src/+gdsge/+codegen/needsCompile.m` | Create | Pure cache decision: does the cpp text differ from the cache file? |
| `src/+gdsge/+codegen/codegen.m` | Create | The public driver: gmod → IR → JSON artifact → generators → cache-gated compile |
| `src/gdsge_codegen.m` | Create | Flat backward-compat shim forwarding to the driver |
| `tests/codegen/tNeedsCompile.m` | Create | Fast unit tests for the cache decision |
| `tests/codegen/tCodegenDriver.m` | Create | Fast unit tests for driver error paths, skip path, shim forwarding |
| `tests/HeatonLucas1996/codegen/tEndToEndHL1996.m` | Create | The Phase 6 gate (Slow): public API end-to-end vs goldens |
| `tests/HeatonLucas1996/codegen/tFunctionalCxxHL1996.m` | Delete | Subsumed by the new gate (per spec §2.3) |
| `PROGRESS.md` | Modify | Check off Phase 6, add changelog entry |

---

### Task 1: Branch setup

**Files:** none (git only)

- [ ] **Step 1: Create the phase branch**

```bash
git checkout -b phase6-vertical-slice
```

Expected: `Switched to a new branch 'phase6-vertical-slice'`. (`git status` must be clean first; if not, stop and report.)

---

### Task 2: `needsCompile` — the pure cache decision

**Files:**
- Create: `src/+gdsge/+codegen/needsCompile.m`
- Test: `tests/codegen/tNeedsCompile.m`

- [ ] **Step 1: Write the failing tests**

Create `tests/codegen/tNeedsCompile.m`:

```matlab
classdef tNeedsCompile < matlab.unittest.TestCase
    methods (Test)
        function trueWhenNoCache(tc)
            tc.verifyTrue(gdsge.codegen.needsCompile('int main(){}', ...
                fullfile(tempname, 'no-such.cache')));
        end
        function falseWhenCacheMatches(tc)
            f = [tempname '.cache'];
            fid = fopen(f, 'w'); fwrite(fid, 'int main(){}'); fclose(fid);
            cleaner = onCleanup(@() delete(f)); %#ok<NASGU>
            tc.verifyFalse(gdsge.codegen.needsCompile('int main(){}', f));
        end
        function trueWhenCacheDiffers(tc)
            f = [tempname '.cache'];
            fid = fopen(f, 'w'); fwrite(fid, 'old text'); fclose(fid);
            cleaner = onCleanup(@() delete(f)); %#ok<NASGU>
            tc.verifyTrue(gdsge.codegen.needsCompile('new text', f));
        end
    end
end
```

- [ ] **Step 2: Run the tests — verify they fail**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','codegen','tNeedsCompile.m')); disp(table(results)); exit(any([results.Failed]))"
```

Expected: FAIL, all three tests, with `Unable to resolve the name 'gdsge.codegen.needsCompile'`.

- [ ] **Step 3: Implement `needsCompile`**

Create `src/+gdsge/+codegen/needsCompile.m`:

```matlab
function tf = needsCompile(cppText, cacheFile)
% NEEDSCOMPILE  True when cppText differs from the cached copy (or no cache).
%   Replicates the old gdsge_codegen cache semantics: a text-identical
%   mex_<model>.cache means the MEX binary is already up to date.
tf = true;
if exist(cacheFile, 'file')
    tf = ~strcmp(fileread(cacheFile), cppText);
end
end
```

- [ ] **Step 4: Run the tests — verify they pass**

Same command as Step 2. Expected: PASS (3/3), exit 0.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/needsCompile.m tests/codegen/tNeedsCompile.m
git commit -m "feat(codegen): needsCompile - pure mex-cache decision"
```

---

### Task 3: Driver error paths (options validation + missing gmod)

**Files:**
- Create: `src/+gdsge/+codegen/codegen.m`
- Test: `tests/codegen/tCodegenDriver.m`

- [ ] **Step 1: Write the failing tests**

Create `tests/codegen/tCodegenDriver.m`:

```matlab
classdef tCodegenDriver < matlab.unittest.TestCase
    % Fast (no-compile) tests of the public codegen driver.
    methods (Test)
        function unknownOptionErrors(tc)
            tc.verifyError(@() gdsge.codegen.codegen('AnyModel', struct('Bogus', 1)), ...
                'gdsge:codegen:unknownOption');
        end
        function genCodeSegmentDeferred(tc)
            tc.verifyError(@() gdsge.codegen.codegen('AnyModel', struct('GenCodeSegment', 1)), ...
                'gdsge:codegen:deferred');
        end
        function missingGmodErrors(tc)
            tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            tc.verifyError(@() gdsge.codegen.codegen('NoSuchModel'), ...
                'gdsge:codegen:gmodNotFound');
        end
    end
end
```

Notes for the implementer:
- Options are validated **before** the gmod lookup, so the two option tests need no
  working-folder fixture.
- `WorkingFolderFixture` creates a fresh temp folder, `cd`s into it, and restores
  the original cwd on teardown — this keeps the repo root free of stray files.

- [ ] **Step 2: Run the tests — verify they fail**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','codegen','tCodegenDriver.m')); disp(table(results)); exit(any([results.Failed]))"
```

Expected: FAIL, all three, with `Unable to resolve the name 'gdsge.codegen.codegen'`.

- [ ] **Step 3: Implement validation + gmod resolution (happy path lands in Task 4)**

Create `src/+gdsge/+codegen/codegen.m`:

```matlab
function ir = codegen(modelName, options)
% CODEGEN  Public codegen driver (error paths; happy path lands in Task 4).
if nargin < 2; options = []; end
validateOptions(options);

gmodFile = fullfile(pwd, [modelName '.gmod']);
if ~exist(gmodFile, 'file')
    error('gdsge:codegen:gmodNotFound', ...
        '%s.gmod not found in %s', modelName, pwd);
end
ir = [];
end

function validateOptions(options)
if isempty(options); return; end
if ~isstruct(options)
    error('gdsge:codegen:optionsNotAStruct', 'options must be a struct.');
end
fn = fieldnames(options);
if any(strcmp(fn, 'GenCodeSegment'))
    error('gdsge:codegen:deferred', ...
        'options.GenCodeSegment is not supported yet (deferred to Phase 7).');
end
if ~isempty(fn)
    error('gdsge:codegen:unknownOption', ...
        'Unknown option field(s): %s. gdsge_codegen accepts no options yet.', ...
        strjoin(fn', ', '));
end
end
```

- [ ] **Step 4: Run the tests — verify they pass**

Same command as Step 2. Expected: PASS (3/3), exit 0.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/codegen.m tests/codegen/tCodegenDriver.m
git commit -m "feat(codegen): driver error paths - options whitelist, gmod lookup"
```

---

### Task 4: Driver happy path (parse → JSON → generate → cache-gated compile)

**Files:**
- Modify: `src/+gdsge/+codegen/codegen.m` (replace entirely with the version below)
- Test: `tests/codegen/tCodegenDriver.m` (add one test method)

- [ ] **Step 1: Write the failing test**

Add this method inside the `methods (Test)` block of `tests/codegen/tCodegenDriver.m`:

```matlab
        function skipsCompileAndWritesArtifactsWhenCacheMatches(tc)
            % Pre-seed mex_HL1996.cache with the exact cpp the driver will
            % generate, so the driver takes the skip path and this test never
            % invokes mex (stays fast: two parses, no compile).
            here = fileparts(mfilename('fullpath'));          % tests/codegen
            modelDir = fullfile(fileparts(here), 'HeatonLucas1996');
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'HL1996.gmod'), work);

            gmodText = fileread(fullfile(modelDir, 'HL1996.gmod'));
            ir0 = gdsge.parser.parseFrontEnd(gmodText, 'HL1996');
            files0 = gdsge.codegen.generateCxx(ir0, work);
            movefile(files0.cppFile, fullfile(work, 'mex_HL1996.cache'));
            delete(fullfile(work, 'compile_HL1996.m'));

            out = evalc('ir = gdsge.codegen.codegen(''HL1996'');');
            tc.verifyTrue(contains(out, 'skip compiling'), ...
                sprintf('expected the skip path, got output:\n%s', out));

            artifacts = {'iter_HL1996.m', 'simulate_HL1996.m', ...
                'mex_HL1996.cpp', 'compile_HL1996.m', 'HL1996.gdsge.json'};
            for i = 1:numel(artifacts)
                tc.verifyTrue(exist(fullfile(work, artifacts{i}), 'file') == 2, ...
                    artifacts{i});
            end
            tc.verifyFalse(exist(fullfile(work, ['mex_HL1996.' mexext]), 'file') == 3, ...
                'driver compiled despite a matching cache');

            % the written IR is the contract: it must round-trip
            decoded = gdsge.ir.decode(fileread(fullfile(work, 'HL1996.gdsge.json')));
            tc.verifyTrue(gdsge.ir.isequalIR(decoded, ir), ...
                'IR JSON does not round-trip');
        end
```

- [ ] **Step 2: Run the test — verify it fails**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','codegen','tCodegenDriver.m')); disp(table(results)); exit(any([results.Failed]))"
```

Expected: the new test FAILS (driver returns `[]` and writes nothing — the `skip compiling` assertion fails); the three Task 3 tests still pass.

- [ ] **Step 3: Implement the full driver**

Replace `src/+gdsge/+codegen/codegen.m` entirely with:

```matlab
function ir = codegen(modelName, options)
% CODEGEN  Public codegen driver: <model>.gmod in cwd -> generated files + MEX.
%   Writes iter_<model>.m, simulate_<model>.m, mex_<model>.cpp,
%   compile_<model>.m and <model>.gdsge.json into the current directory, then
%   compiles the MEX — skipped when the C++ matches mex_<model>.cache. The
%   cache is written only after a successful compile, so failed builds retry.
%   Returns the IR struct. Called by the flat shim gdsge_codegen.
if nargin < 2; options = []; end
validateOptions(options);

gmodFile = fullfile(pwd, [modelName '.gmod']);
if ~exist(gmodFile, 'file')
    error('gdsge:codegen:gmodNotFound', ...
        '%s.gmod not found in %s', modelName, pwd);
end

fprintf('Parsing gmod file: ');
ir = gdsge.parser.parseFrontEnd(fileread(gmodFile), modelName);
writeText(fullfile(pwd, [modelName '.gdsge.json']), gdsge.ir.encode(ir));
gdsge.codegen.generateMatlab(ir, pwd);
files = gdsge.codegen.generateCxx(ir, pwd);
fprintf('ok\n');

cppText   = fileread(files.cppFile);
cacheFile = fullfile(pwd, ['mex_' modelName '.cache']);
if gdsge.codegen.needsCompile(cppText, cacheFile)
    fprintf('Compile mex file:\n');
    gdsge.runtime.ensurePath();
    feval(['compile_' modelName]);
    writeText(cacheFile, cppText);
else
    fprintf('C++ source unchanged, skip compiling.\n');
end
end

function validateOptions(options)
if isempty(options); return; end
if ~isstruct(options)
    error('gdsge:codegen:optionsNotAStruct', 'options must be a struct.');
end
fn = fieldnames(options);
if any(strcmp(fn, 'GenCodeSegment'))
    error('gdsge:codegen:deferred', ...
        'options.GenCodeSegment is not supported yet (deferred to Phase 7).');
end
if ~isempty(fn)
    error('gdsge:codegen:unknownOption', ...
        'Unknown option field(s): %s. gdsge_codegen accepts no options yet.', ...
        strjoin(fn', ', '));
end
end

function writeText(p, txt)
fid = fopen(p, 'w');
if fid < 0
    error('gdsge:codegen:cannotWrite', 'Cannot write %s', p);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, txt);
end
```

Implementation notes:
- `feval(['compile_' modelName])` resolves the freshly written `compile_<model>.m`
  because cwd always has top resolution priority, even from inside a package function.
  This matches how the old driver (`eval(['compile_' modelName])`) and the Phase 5
  gate test invoked it.
- `gdsge.runtime.ensurePath()` before compiling mirrors the proven Phase 5 gate
  sequence (it puts the `essential_blas.dll` folder on the system PATH; idempotent).

- [ ] **Step 4: Run the tests — verify they pass**

Same command as Step 2. Expected: PASS (4/4), exit 0.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/codegen.m tests/codegen/tCodegenDriver.m
git commit -m "feat(codegen): driver happy path - IR json artifact + cache-gated compile"
```

---

### Task 5: The flat backward-compat shim

**Files:**
- Create: `src/gdsge_codegen.m`
- Test: `tests/codegen/tCodegenDriver.m` (add one test method)

- [ ] **Step 1: Write the failing test**

Add this method inside the `methods (Test)` block of `tests/codegen/tCodegenDriver.m`:

```matlab
        function shimForwardsToDriver(tc)
            tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            tc.verifyError(@() gdsge_codegen('NoSuchModel'), ...
                'gdsge:codegen:gmodNotFound');
        end
```

- [ ] **Step 2: Run the test — verify it fails**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','codegen','tCodegenDriver.m')); disp(table(results)); exit(any([results.Failed]))"
```

Expected: the new test FAILS with `Unrecognized function or variable 'gdsge_codegen'`; the other four still pass.

- [ ] **Step 3: Implement the shim**

Create `src/gdsge_codegen.m`:

```matlab
function model = gdsge_codegen(modelName, options)
% GDSGE_CODEGEN  Backward-compatible flat entry point (frozen public API).
%   Thin shim over gdsge.codegen.codegen: parses <modelName>.gmod from the
%   current directory, writes the generated files there, and compiles the MEX
%   when its C++ changed. Returns the model IR struct.
if nargin < 2; options = []; end
model = gdsge.codegen.codegen(modelName, options);
end
```

- [ ] **Step 4: Run the tests — verify they pass**

Same command as Step 2. Expected: PASS (5/5), exit 0.

- [ ] **Step 5: Commit**

```bash
git add src/gdsge_codegen.m tests/codegen/tCodegenDriver.m
git commit -m "feat(codegen): flat gdsge_codegen shim - frozen public entry point"
```

---

### Task 6: The Phase 6 gate — end-to-end through the public API

**Files:**
- Create: `tests/HeatonLucas1996/codegen/tEndToEndHL1996.m`
- Delete: `tests/HeatonLucas1996/codegen/tFunctionalCxxHL1996.m`

- [ ] **Step 1: Write the gate test**

Create `tests/HeatonLucas1996/codegen/tEndToEndHL1996.m`. The golden assertions are
carried over verbatim from `tFunctionalCxxHL1996.m`; new are the shim entry, the
artifact checks, the JSON round-trip, and the second-run skip check.

```matlab
classdef tEndToEndHL1996 < matlab.unittest.TestCase
    % PHASE 6 GATE: the public surface (gdsge_codegen -> iter_HL1996 ->
    % simulate_HL1996), driven exactly like the old toolbox's test.m,
    % reproduces the committed goldens. Replaces the Phase 5 gate
    % tFunctionalCxxHL1996 (subsumed). Slow (~minutes): MEX compile + full
    % policy iteration + 10k-period simulation.
    properties (Constant)
        RelTol = 1e-4;
        AbsTol = 1e-4;
    end
    methods (Test, TestTags = {'Slow'})
        function publicApiMatchesGolden(tc)
            here = fileparts(mfilename('fullpath'));     % tests/HeatonLucas1996/codegen
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'HL1996.gmod'), work);

            % ---- codegen through the flat shim ------------------------------
            ir = gdsge_codegen('HL1996');

            mexFile = fullfile(work, ['mex_HL1996.' mexext]);
            artifacts = {'iter_HL1996.m', 'simulate_HL1996.m', 'mex_HL1996.cpp', ...
                'compile_HL1996.m', 'HL1996.gdsge.json', 'mex_HL1996.cache'};
            for i = 1:numel(artifacts)
                tc.assertTrue(exist(fullfile(work, artifacts{i}), 'file') == 2, ...
                    artifacts{i});
            end
            tc.assertTrue(exist(mexFile, 'file') == 3, 'new MEX did not compile');

            % the written IR is the contract: it must round-trip
            decoded = gdsge.ir.decode(fileread(fullfile(work, 'HL1996.gdsge.json')));
            tc.verifyTrue(gdsge.ir.isequalIR(decoded, ir), ...
                'IR JSON does not round-trip');

            % ---- second run takes the skip path ------------------------------
            before = dir(mexFile);
            gdsge_codegen('HL1996');
            after = dir(mexFile);
            tc.verifyEqual(after.datenum, before.datenum, ...
                'second gdsge_codegen run recompiled despite unchanged C++');

            % ---- iterate ------------------------------------------------------
            opts = struct('SaveFreq', inf, 'NoSave', 1);
            IterRslt = iter_HL1996(opts);

            golden = load(fullfile(modelDir, 'golden', 'IterRslt.mat'));
            G = golden.IterRslt;
            tc.verifyLessThan(IterRslt.Metric, 1e-6);
            tc.verifyGreaterThan(IterRslt.Iter, 100);
            tc.verifyLessThan(IterRslt.Iter, 400);     % golden converged at 209
            r = gdsgetest.compareNumericClose(IterRslt.shock_trans, G.shock_trans, 1e-12, 1e-15);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_state, G.var_state, 1e-12, 1e-15);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_policy, G.var_policy, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_aux, G.var_aux, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_interp, G.var_interp, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));

            % ---- simulate -----------------------------------------------------
            rng(0823);                                  % same seed as the golden capture
            SimuRslt = simulate_HL1996(IterRslt);
            gs = load(fullfile(modelDir, 'golden', 'SimuRslt.mat'));
            GS = gs.SimuRslt;
            % identical seed + transition matrix => identical shock path
            tc.verifyEqual(SimuRslt.shock, GS.shock, ...
                'shock paths differ — check shock_trans bit-identity first');
            flds = {'w1','c1','c2','ps','pb','equity_premium'};
            T0 = 100;   % early path: policy diffs have not accumulated yet
            for i = 1:numel(flds)
                a = SimuRslt.(flds{i}); b = GS.(flds{i});
                tc.verifyEqual(size(a), size(b), flds{i});
                r = gdsgetest.compareNumericClose(a(:,1:T0), b(:,1:T0), 1e-2, 1e-2);
                tc.verifyTrue(r.pass, sprintf('%s early path: %s', flds{i}, ...
                    strjoin(r.failures, newline)));
                % long-run moments (paths decorrelate; the distribution must not)
                tc.verifyLessThan(abs(mean(a(:)) - mean(b(:))), 5e-3, flds{i});
                tc.verifyLessThan(abs(std(a(:)) - std(b(:))), 5e-3, flds{i});
            end
        end
    end
end
```

- [ ] **Step 2: Delete the superseded Phase 5 gate**

```bash
git rm tests/HeatonLucas1996/codegen/tFunctionalCxxHL1996.m
```

- [ ] **Step 3: Run the gate (slow, ~5–10 minutes: compile + 209 iterations + simulation)**

```
matlab -batch "addpath('tests'); addpath('src'); addpath(fullfile('src','kernels')); results = runtests(fullfile('tests','HeatonLucas1996','codegen','tEndToEndHL1996.m')); disp(table(results)); exit(any([results.Failed]))"
```

Expected: PASS, exit 0. Iteration output should show convergence around Iter=209,
Metric≈9.58e-07. If the compile fails or goldens mismatch, STOP and debug with the
superpowers:systematic-debugging skill — do not loosen tolerances.

- [ ] **Step 4: Commit**

```bash
git add tests/HeatonLucas1996/codegen/tEndToEndHL1996.m
git commit -m "test(codegen): Phase 6 gate - public API end-to-end vs goldens, replaces Phase 5 gate"
```

(The `git rm` from Step 2 is already staged and lands in this commit.)

---

### Task 7: Full suite green + PROGRESS.md

**Files:**
- Modify: `PROGRESS.md`

- [ ] **Step 1: Run the full suite**

```
matlab -batch "cd('tests'); run_tests"
```

Expected: exit 0; `tests/results/junit.xml` shows 0 failures. This includes the new
gate and the remaining slow tests (`tFunctionalHL1996`, `tMexEquivalenceHL1996`) —
budget ~15–30 minutes.

- [ ] **Step 2: Update PROGRESS.md**

In the Phases list, replace the Phase 6 line:

```markdown
- ☐ **Phase 6 — Vertical slice green** (HL1996 end-to-end vs golden — architecture proven) — NEXT
```

with:

```markdown
- ☑ **Phase 6 — Vertical slice green** (done 2026-06-12)
  - ☑ `src/gdsge_codegen.m` flat shim + `gdsge.codegen.codegen` driver (parse →
    `<model>.gdsge.json` → generateMatlab/Cxx into cwd → cache-gated MEX compile)
  - ☑ Old cache semantics replicated (`mex_<model>.cache`; cache written only after a
    successful compile); options whitelist errors; `GenCodeSegment` deferred error
  - ☑ Gate green: `tEndToEndHL1996` drives the public API like the old `test.m`,
    asserts artifacts + IR JSON round-trip + second-run compile skip + goldens;
    replaces `tFunctionalCxxHL1996`
  - ☑ Deferred to Phase 7: `gdsge.m` orchestrator, `GenCodeSegment`, legacy 5-output
    signature
```

and update the next-phase marker by changing the Phase 7 line's suffix from nothing to ` — NEXT`:

```markdown
- ☐ **Phase 7 — Widen for backward-compat** (Bianchi2011_asg, Mendoza2010, CaoKS2016, GLSW2020, Barro_et_al_2017; ASG/var_tensor/model_init/pre-post_iter/macros) — NEXT
```

Add at the top of the Changelog section:

```markdown
- 2026-06-12: **Phase 6 complete.** Public wiring: flat `gdsge_codegen` shim +
  `gdsge.codegen.codegen` driver (gmod from cwd → IR → `<model>.gdsge.json` via
  `gdsge.ir.encode` → existing generators → `mex_<model>.cache`-gated compile,
  cache written only after success). Options whitelist with informative errors;
  `GenCodeSegment` raises a deferred-to-Phase-7 error. New Slow gate
  `tEndToEndHL1996` proves the architecture through the public surface (artifacts,
  IR round-trip, second-run skip, golden match) and replaces `tFunctionalCxxHL1996`.
  All on branch `phase6-vertical-slice`. Next: Phase 7 (widen for backward-compat).
```

- [ ] **Step 3: Commit**

```bash
git add PROGRESS.md
git commit -m "docs(progress): Phase 6 complete - vertical slice green via public API"
```

- [ ] **Step 4: Finish the branch**

Use the superpowers:finishing-a-development-branch skill to decide merge vs PR
(prior phases were merged into `main`).

---

## Self-review notes (already applied)

- **Spec coverage:** §3.1 shim → Task 5; §3.2 driver steps 1–6 → Tasks 3–4 (+ Task 2
  for the cache decision); §4 error handling → Task 3 (validation, gmod) and Task 4
  (cache written only after successful compile — exercised by the gate's skip check);
  §5.1 gate items 1–5 → Task 6; §5.2 fast tests → Tasks 2–5; §5.3 done-means →
  Task 7; §6 deferred → nothing to implement.
- **Type consistency:** `needsCompile(cppText, cacheFile)` (Task 2) is what Task 4
  calls; `generateCxx` returns `files.cppFile`/`files.compileFile` (verified against
  `src/+gdsge/+codegen/generateCxx.m`); `gdsge.ir.encode/decode/isequalIR` signatures
  verified against `src/+gdsge/+ir/`.
- **Failed-compile-no-cache** (§4) is implemented by ordering (`feval` before
  `writeText(cacheFile, …)`) — an error in `compile_<model>` propagates and skips the
  cache write. No dedicated test: forcing a MEX failure deterministically would need
  a broken toolchain fixture; the ordering is structurally evident and the retry path
  is the same code the gate exercises.
