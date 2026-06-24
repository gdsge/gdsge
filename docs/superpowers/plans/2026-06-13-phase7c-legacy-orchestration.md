# Phase 7c — Legacy Orchestration & ASG-interp Simulate — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the four small deferred legacy-surface items — the `gdsge.m` orchestrator, the legacy 5-output `gdsge_codegen` signature, the `GenCodeSegment` informative error, and the ASG + `SIMU_INTERP` simulate variant — so the public API matches the old toolbox's surface.

**Architecture:** `gdsge.codegen.codegen` gains an optional second output carrying the generated code strings; the flat `gdsge_codegen` shim re-exposes the old 5 positional outputs from it; a new flat `gdsge.m` mirrors the old orchestrator (codegen → cache-gated solve → simulate → `eq`) without `v2struct`. `GenCodeSegment` raises a clear "no segment decomposition in the thin-file architecture" error. ASG + `SIMU_INTERP` is delivered by extending the existing spline `emitSimulateInterp` with an ASG branch (construct from `asg_output_struct`, `eval_vec` per period) and removing the 7b parser guard; validated differentially against a freshly captured old-toolbox golden (CaoKS2016 with `SIMU_INTERP` flipped on).

**Tech Stack:** MATLAB R2025b (`matlab.unittest`), the `+gdsge` package (`gdsge.parser.*`, `gdsge.codegen.*`, `gdsge.runtime.*`), the vendored `asg` class, C MEX (only for the ASG iter, not for interp-simulate).

**Spec:** `docs/superpowers/specs/2026-06-13-phase7c-legacy-orchestration-design.md`
**Branch:** `phase7c-legacy-orchestration` (already created; spec + PROGRESS split committed).

---

## Conventions for every task

- **MATLAB binary:** `matlab` resolves to `C:\Program Files\MATLAB\R2025b\bin\matlab.exe`. If `matlab` is not on PATH, use the full path.
- **Run one test class** (fast loop), from the repo root `D:\refactor_gdsge`:
  ```
  matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/<dir>/<Class>.m'); disp(r); assert(~any([r.Failed]))"
  ```
  Exit code 0 = pass; a non-zero exit (assertion) = fail.
- **Run the full suite** (final gate): `matlab -batch "cd('tests'); run_tests"` — exit 0 = all pass; authoritative report `tests/results/junit.xml`.
- **Commit trailer:** every commit message ends with:
  ```
  Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
  ```
- Work stays on branch `phase7c-legacy-orchestration`.

---

## Task 1: `gdsge.codegen.codegen` exposes generated artifacts (2nd output)

Add an optional second output `gen` carrying the generated code strings and the pre-overwrite cache content. Single-output callers are unaffected.

**Files:**
- Modify: `src/+gdsge/+codegen/codegen.m`
- Test: `tests/codegen/tCodegenDriver.m` (add one fast, no-compile method)

- [ ] **Step 1: Write the failing test** — add this method inside the `methods (Test)` block of `tests/codegen/tCodegenDriver.m` (after `skipsCompileAndWritesArtifactsWhenCacheMatches`):

```matlab
        function exposesGeneratedArtifacts(tc)
            % [ir,gen] returns the real generated strings. Pre-seed the cache so
            % the driver skips compiling and the test stays fast (no mex).
            here = fileparts(mfilename('fullpath'));          % tests/codegen
            modelDir = fullfile(fileparts(here), 'HeatonLucas1996');
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'HL1996.gmod'), work);

            ir0 = gdsge.parser.parseFrontEnd( ...
                fileread(fullfile(modelDir, 'HL1996.gmod')), 'HL1996');
            files0 = gdsge.codegen.generateCxx(ir0, work);
            copyfile(files0.cppFile, fullfile(work, 'mex_HL1996.cache'));

            [ir, gen] = gdsge.codegen.codegen('HL1996'); %#ok<ASGLU>
            tc.verifyEqual(gen.iterCode,  fileread(fullfile(work, 'iter_HL1996.m')));
            tc.verifyEqual(gen.cppCode,   fileread(fullfile(work, 'mex_HL1996.cpp')));
            tc.verifyEqual(gen.compileCode, fileread(fullfile(work, 'compile_HL1996.m')));
            tc.verifyEqual(gen.cppCache,  gen.cppCode);   % skip path: prior cache matched
            tc.verifyTrue(isstruct(gen.codeSegment));
            tc.verifyTrue(isfield(gen.codeSegment, 'iterCode'));
            tc.verifyTrue(isfield(gen.codeSegment, 'cppCode'));
        end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tCodegenDriver.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: FAIL — `codegen` currently returns only `ir`, so `[ir,gen] = ...` errors ("Too many output arguments") or `gen` is undefined.

- [ ] **Step 3: Implement the 2nd output** — replace the body of `src/+gdsge/+codegen/codegen.m` (keep `validateOptions` below unchanged for now; Task 3 edits it):

```matlab
function [ir, gen] = codegen(modelName, options)
% CODEGEN  Public codegen driver: <model>.gmod in cwd -> generated files + MEX.
%   Writes iter_<model>.m, simulate_<model>.m, mex_<model>.cpp,
%   compile_<model>.m and <model>.gdsge.json into the current directory, then
%   compiles the MEX — skipped when the C++ matches mex_<model>.cache. The
%   cache is written only after a successful compile, so failed builds retry.
%   Returns the IR struct. With a second output, also returns gen, a struct of
%   the generated code strings (iterCode/simulateCode/cppCode/compileCode), the
%   pre-overwrite cache content (cppCache), and codeSegment (the real generated
%   strings — not the old template fragments). Called by the flat shim
%   gdsge_codegen, which re-exposes these as its legacy positional outputs.
if nargin < 2; options = []; end
validateOptions(options);

gmodFile = fullfile(pwd, [modelName '.gmod']);
if ~exist(gmodFile, 'file')
    error('gdsge:codegen:gmodNotFound', ...
        '%s.gmod not found in %s', modelName, pwd);
end

fprintf('Parsing gmod file: ');
ir = gdsge.parser.parseFrontEnd(fileread(gmodFile), modelName);
gdsge.codegen.writeText(fullfile(pwd, [modelName '.gdsge.json']), gdsge.ir.encode(ir));
matFiles = gdsge.codegen.generateMatlab(ir, pwd);
files    = gdsge.codegen.generateCxx(ir, pwd);
fprintf('ok\n');

cppText   = fileread(files.cppFile);
cacheFile = fullfile(pwd, ['mex_' modelName '.cache']);
cppCache  = '';
if exist(cacheFile, 'file'); cppCache = fileread(cacheFile); end
if gdsge.codegen.needsCompile(cppText, cacheFile)
    fprintf('Compile mex file:\n');
    gdsge.runtime.ensurePath();
    feval(['compile_' modelName]);
    gdsge.codegen.writeText(cacheFile, cppText);
else
    fprintf('C++ source unchanged, skip compiling.\n');
end

if nargout > 1
    iterCode     = fileread(matFiles.iterFile);
    simulateCode = fileread(matFiles.simulateFile);
    compileCode  = fileread(fullfile(pwd, ['compile_' modelName '.m']));
    gen = struct('iterCode', iterCode, 'simulateCode', simulateCode, ...
        'cppCode', cppText, 'compileCode', compileCode, 'cppCache', cppCache);
    gen.codeSegment = struct('iterCode', iterCode, 'simulateCode', simulateCode, ...
        'cppCode', cppText, 'compileCode', compileCode);
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tCodegenDriver.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: PASS (all methods in `tCodegenDriver`, including the existing ones, green).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/codegen.m tests/codegen/tCodegenDriver.m
git commit -m "feat(codegen): expose generated code strings as codegen 2nd output

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: `gdsge_codegen` legacy 5-output signature

Re-expose the old positional outputs `[model,iterCode,cppCache,cppCode,codeSegment]`. Single-output callers keep working.

**Files:**
- Modify: `src/gdsge_codegen.m`
- Test: `tests/codegen/tCodegenDriver.m` (add one fast method)

- [ ] **Step 1: Write the failing test** — add this method to `tests/codegen/tCodegenDriver.m`:

```matlab
        function shimReturnsLegacyFiveOutputs(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fullfile(fileparts(here), 'HeatonLucas1996');
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'HL1996.gmod'), work);

            ir0 = gdsge.parser.parseFrontEnd( ...
                fileread(fullfile(modelDir, 'HL1996.gmod')), 'HL1996');
            files0 = gdsge.codegen.generateCxx(ir0, work);
            copyfile(files0.cppFile, fullfile(work, 'mex_HL1996.cache'));

            % single-output back-compat
            m1 = gdsge_codegen('HL1996');
            tc.verifyTrue(isstruct(m1) && isfield(m1, 'modelName'));

            % legacy 5-output
            [model, iterCode, cppCache, cppCode, codeSegment] = gdsge_codegen('HL1996');
            tc.verifyTrue(isstruct(model) && isfield(model, 'modelName'));
            tc.verifyEqual(iterCode, fileread(fullfile(work, 'iter_HL1996.m')));
            tc.verifyEqual(cppCode,  fileread(fullfile(work, 'mex_HL1996.cpp')));
            tc.verifyEqual(cppCache, cppCode);   % skip path
            tc.verifyTrue(isstruct(codeSegment) && isfield(codeSegment, 'iterCode'));
        end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tCodegenDriver.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: FAIL — the current shim has only one output, so the 5-output destructure errors ("Too many output arguments").

- [ ] **Step 3: Implement the 5-output shim** — replace the whole `src/gdsge_codegen.m`:

```matlab
function [model, iterCode, cppCache, cppCode, codeSegment] = gdsge_codegen(modelName, options)
% GDSGE_CODEGEN  Backward-compatible flat entry point (frozen public API).
%   Thin shim over gdsge.codegen.codegen: parses <modelName>.gmod from the
%   current directory, writes the generated files there, and compiles the MEX
%   when its C++ changed. Returns the model IR struct. The legacy positional
%   outputs are re-exposed for the gdsge.m orchestrator:
%     iterCode    - text of the generated iter_<model>.m
%     cppCache    - mex_<model>.cache content before this run (''if none)
%     cppCode     - text of the generated mex_<model>.cpp
%     codeSegment - struct of the real generated strings (iterCode/simulateCode/
%                   cppCode/compileCode); NOT the old template fragments, which
%                   have no counterpart in the thin-file architecture.
if nargin < 2; options = []; end
if nargout <= 1
    model = gdsge.codegen.codegen(modelName, options);
else
    [model, gen] = gdsge.codegen.codegen(modelName, options);
    iterCode    = gen.iterCode;
    cppCache    = gen.cppCache;
    cppCode     = gen.cppCode;
    codeSegment = gen.codeSegment;
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tCodegenDriver.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/gdsge_codegen.m tests/codegen/tCodegenDriver.m
git commit -m "feat(codegen): restore legacy 5-output gdsge_codegen signature

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: `GenCodeSegment` informative error

Replace the placeholder "deferred to Phase 7c" error with the final, architecture-honest message.

**Files:**
- Modify: `src/+gdsge/+codegen/codegen.m` (the `validateOptions` local function)
- Test: `tests/codegen/tCodegenDriver.m` (flip the existing `genCodeSegmentDeferred` method)

- [ ] **Step 1: Update the test (it will fail against the old identifier)** — in `tests/codegen/tCodegenDriver.m`, replace the `genCodeSegmentDeferred` method with:

```matlab
        function genCodeSegmentUnsupported(tc)
            tc.verifyError(@() gdsge.codegen.codegen('AnyModel', struct('GenCodeSegment', 1)), ...
                'gdsge:codegen:genCodeSegmentUnsupported');
        end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tCodegenDriver.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: FAIL — the code still raises `gdsge:codegen:deferred`, not `gdsge:codegen:genCodeSegmentUnsupported`.

- [ ] **Step 3: Implement the error** — in `src/+gdsge/+codegen/codegen.m`, inside `validateOptions`, replace the `GenCodeSegment` branch:

```matlab
% GenCodeSegment is a recognized legacy option with no counterpart in the
% refactored architecture; give it a specific diagnostic before the generic
% unknown-field sweep.
if any(strcmp(fn, 'GenCodeSegment'))
    error('gdsge:codegen:genCodeSegmentUnsupported', ...
        ['options.GenCodeSegment is not supported: the refactored thin-file ' ...
         'architecture (explicit named locals over gdsge.runtime helpers, no ' ...
         'v2struct) has no equivalent segment decomposition. Inspect the ' ...
         'generated iter_<model>.m and mex_<model>.cpp directly.']);
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tCodegenDriver.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/codegen.m tests/codegen/tCodegenDriver.m
git commit -m "feat(codegen): GenCodeSegment raises an architecture-honest error

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: `gdsge.m` orchestrator

A new flat shim mirroring the old driver: codegen → cache-gated solve (persist `IterRslt_<model>.mat` + `iter_<model>.cache`) → simulate → `eq` (no `v2struct`).

**Files:**
- Create: `src/gdsge.m`
- Test: `tests/HeatonLucas1996/codegen/tGdsgeDriverHL1996.m`

- [ ] **Step 1: Write the failing test** — create `tests/HeatonLucas1996/codegen/tGdsgeDriverHL1996.m`:

```matlab
classdef tGdsgeDriverHL1996 < matlab.unittest.TestCase
    % PHASE 7c GATE: the gdsge() orchestrator (codegen -> cache-gated solve ->
    % simulate -> eq) reproduces HL1996's golden IterRslt and skips the solve on
    % a second run. Reuses the existing golden (no new capture). gdsge() does not
    % seed rng, so SmltRslt is checked structurally, not bit-for-bit. Slow.
    methods (Test, TestTags = {'Slow'})
        function orchestratesAndCaches(tc)
            here = fileparts(mfilename('fullpath'));          % tests/HeatonLucas1996/codegen
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'HL1996.gmod'), work);

            eq1 = gdsge('HL1996');
            tc.verifyTrue(all(isfield(eq1, {'model', 'IterRslt', 'SmltRslt'})));

            golden = load(fullfile(modelDir, 'golden', 'IterRslt.mat'));
            G = golden.IterRslt;
            R = eq1.IterRslt;
            tc.verifyLessThan(R.Metric, 1e-6);
            tc.verifyGreaterThan(R.Iter, 100);
            tc.verifyLessThan(R.Iter, 400);                  % golden converged at 209
            r = gdsgetest.compareNumericClose(R.var_policy, G.var_policy, 1e-4, 1e-4);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(R.var_interp, G.var_interp, 1e-4, 1e-4);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));

            % SmltRslt: structural only (no rng seeding inside gdsge())
            tc.verifyTrue(isstruct(eq1.SmltRslt));
            tc.verifyTrue(isfield(eq1.SmltRslt, 'shock'));
            tc.verifyTrue(all(isfinite(eq1.SmltRslt.w1(:))));

            % cache gate: tamper the saved IterRslt, re-run, prove the solve was skipped
            iterRsltName = fullfile(work, 'IterRslt_HL1996.mat');
            tc.assertTrue(exist(iterRsltName, 'file') == 2, 'IterRslt mat not saved');
            S = load(iterRsltName, 'IterRslt');
            IterRslt = S.IterRslt;
            IterRslt.GDSGE_CACHE_SENTINEL = 12345; %#ok<STRNU>
            save(iterRsltName, 'IterRslt');

            eq2 = gdsge('HL1996');
            tc.verifyTrue(isfield(eq2.IterRslt, 'GDSGE_CACHE_SENTINEL'), ...
                'second gdsge() re-solved instead of loading the cached IterRslt');
            tc.verifyEqual(eq2.IterRslt.GDSGE_CACHE_SENTINEL, 12345);
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/HeatonLucas1996/codegen/tGdsgeDriverHL1996.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: FAIL — `gdsge` (the function) does not exist in `src/`, so the call errors.

- [ ] **Step 3: Implement the orchestrator** — create `src/gdsge.m`:

```matlab
function eq = gdsge(modelName)
% GDSGE  Backward-compatible flat orchestrator (frozen public API).
%   Generates code (gdsge_codegen), solves policies (iter_<model>) with an
%   IterRslt cache, simulates (simulate_<model>), and returns eq with fields
%   model, IterRslt, SmltRslt. The solve is skipped when the generated iter and
%   C++ code are unchanged and IterRslt_<model>.mat already exists. No v2struct.
fprintf('GDSGE: A Toolbox for Solving Global DSGE Models:\n');

[model, iterCode, cppCache, cppCode] = gdsge_codegen(modelName);

iterCacheName = ['iter_' modelName '.cache'];
iterCache = '';
if exist(iterCacheName, 'file') ~= 0
    iterCache = fileread(iterCacheName);
end

iterRsltName = ['IterRslt_' modelName '.mat'];
if strcmp(iterCache, iterCode) && strcmp(cppCache, cppCode) && exist(iterRsltName, 'file') ~= 0
    fprintf('Iter file not changed. IterRslt found. Load results:\n');
    loaded = load(iterRsltName, 'IterRslt');
    IterRslt = loaded.IterRslt;
    fprintf('ok\n');
else
    fprintf('Solve policies:\n');
    IterRslt = feval(['iter_' modelName]);
    save(iterRsltName, 'IterRslt');
    fid = fopen(iterCacheName, 'w');
    fprintf(fid, '%s', iterCode);
    fclose(fid);
    fprintf('ok\n');
end

fprintf('Simulate:\n');
SmltRslt = feval(['simulate_' modelName], IterRslt);
fprintf('ok\n');

eq.model    = model;
eq.IterRslt = IterRslt;
eq.SmltRslt = SmltRslt;
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/HeatonLucas1996/codegen/tGdsgeDriverHL1996.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: PASS (Slow — full HL1996 codegen + solve + two simulates; minutes).

- [ ] **Step 5: Commit**

```bash
git add src/gdsge.m tests/HeatonLucas1996/codegen/tGdsgeDriverHL1996.m
git commit -m "feat(runtime): gdsge.m orchestrator (codegen->cache-gated solve->simulate)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Allow ASG + `SIMU_INTERP` in the parser

Remove the 7b guard so `USE_ASG` with `SIMU_INTERP` resolves cleanly (the `SIMU_RESOLVE ⊕ SIMU_INTERP` and method-exclusivity checks stay).

**Files:**
- Modify: `src/+gdsge/+parser/resolveOptions.m`
- Test: `tests/parser/tResolveOptions.m` (flip `asgRejectsSimuInterp`)

- [ ] **Step 1: Update the test** — in `tests/parser/tResolveOptions.m`, replace the `asgRejectsSimuInterp` method with:

```matlab
        function asgAcceptsSimuInterp(tc)
            ws = struct('USE_ASG', 1, 'USE_SPLINE', 0, ...
                'SIMU_INTERP', 1, 'SIMU_RESOLVE', 0);
            o = gdsge.parser.resolveOptions(ws);
            tc.verifyEqual(o.interpMethod, 'asg');
            tc.verifyEqual(o.simuInterp, 1);
            tc.verifyEqual(o.simuResolve, 0);
        end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/parser/tResolveOptions.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: FAIL — `resolveOptions` still raises `gdsge:parser:asgSimuInterp` for this workspace.

- [ ] **Step 3: Remove the guard** — in `src/+gdsge/+parser/resolveOptions.m`, delete these lines (currently 33–36):

```matlab
if strcmp(method, 'asg') && o.simuInterp ~= 0
    error('gdsge:parser:asgSimuInterp', ...
        'SIMU_INTERP with USE_ASG is not supported yet (deferred to Phase 7c).');
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/parser/tResolveOptions.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/resolveOptions.m tests/parser/tResolveOptions.m
git commit -m "feat(parser): allow USE_ASG with SIMU_INTERP (remove 7b guard)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: ASG branch in `emitSimulateInterp` + `generateMatlab` selection

Extend the spline interp-simulate emitter with an ASG branch (construct from `asg_output_struct`, `eval_vec` per period), and route ASG + `simuInterp` to it.

**Files:**
- Modify: `src/+gdsge/+codegen/+mat/emitSimulateInterp.m`
- Modify: `src/+gdsge/+codegen/generateMatlab.m`
- Test: `tests/codegen/tEmitSimulateInterp.m` (add ASG-branch + selection methods)

- [ ] **Step 1: Write the failing tests** — add these two methods to the `methods (Test)` block of `tests/codegen/tEmitSimulateInterp.m`:

```matlab
        function emitsAsgInterpEval(tc)
            ir = gdsgefix.minimalIRAsg();
            ir.options.simuInterp = 1; ir.options.simuResolve = 0;
            txt = gdsge.codegen.mat.emitSimulateInterp(ir);
            tc.verifyTrue(contains(txt, ...
                'asg.construct_from_struct(IterRslt.asg_output_struct)'));
            tc.verifyTrue(contains(txt, ...
                'GDSGE_ASG_OUTPUT_INTERP.eval_vec(SimuRslt.shock(:,GDSGE_t)'''));
            tc.verifyFalse(contains(txt, 'GDSGE_PP = IterRslt.output_interp;'));
            tc.verifyFalse(contains(txt, 'myppual_mex'));
            tc.verifyFalse(contains(txt, 'solveProblems'));
            t = mtree(txt);
            tc.verifyEqual(count(mtfind(t, 'Kind', 'ERR')), 0);
        end
        function generateMatlabPicksAsgInterpVariant(tc)
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            ir = gdsgefix.minimalIRAsg();
            ir.options.simuInterp = 1; ir.options.simuResolve = 0;
            files = gdsge.codegen.generateMatlab(ir, work);
            sim = fileread(files.simulateFile);
            tc.verifyTrue(contains(sim, 'asg_output_struct'));
            tc.verifyFalse(contains(sim, 'solveProblems'));   % interp = no re-solve
            % the default ASG path (simuResolve) still re-solves
            ir2 = gdsgefix.minimalIRAsg();
            files2 = gdsge.codegen.generateMatlab(ir2, work);
            tc.verifyTrue(contains(fileread(files2.simulateFile), 'solveProblems'));
        end
```

- [ ] **Step 2: Run the tests to verify they fail**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tEmitSimulateInterp.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: FAIL — `emitSimulateInterp` emits the spline `GDSGE_PP`/`myppual_mex` form for the ASG IR; `generateMatlab` routes ASG to `emitSimulateAsg` regardless of `simuInterp`.

- [ ] **Step 3a: Add the ASG branch in `emitSimulateInterp.m`.** In `src/+gdsge/+codegen/+mat/emitSimulateInterp.m`:

First, update the header comment (lines 1–6) to:

```matlab
function txt = emitSimulateInterp(ir)
% EMITSIMULATEINTERP  simulate_<model>.m, SIMU_INTERP variant: no per-period
%   re-solve. var_output values come from the output interpolant evaluated at
%   the realized (shock, state) points each period; var_simu fields and state
%   transitions are then assigned exactly as in the resolve variant. Spline:
%   IterRslt.output_interp via myppual_mex. ASG (interpMethod='asg'):
%   asg.construct_from_struct(IterRslt.asg_output_struct) via eval_vec. Old
%   parity: code_template/simulate_interp_template.m + gdsge_parser.m:1684-1737.
```

Then, right after `m = ir.modelName;` (line 8), add:

```matlab
isAsg = strcmp(ir.options.interpMethod, 'asg');
```

Then replace the line `w.add('GDSGE_PP = IterRslt.output_interp;');` (line 28) with:

```matlab
if isAsg
    w.add('GDSGE_ASG_OUTPUT_INTERP = asg.construct_from_struct(IterRslt.asg_output_struct);');
else
    w.add('GDSGE_PP = IterRslt.output_interp;');
end
```

Then replace the spline eval block (currently the `if shock_num>1 ... end` block at lines 64–70):

```matlab
w.add('if shock_num>1');
w.add('    GDSGE_INTERP_RESULTS = myppual_mex(int32(NumThreads),GDSGE_PP.breaks,GDSGE_PP.coefs,...');
w.add('        int32(GDSGE_PP.pieces),int32(GDSGE_PP.order),int32(GDSGE_PP.dim),''not-a-knot'',[SimuRslt.shock(:,GDSGE_t)'';%s],[],[],[]);', simuStateRows);
w.add('else');
w.add('    GDSGE_INTERP_RESULTS = myppual_mex(int32(NumThreads),GDSGE_PP.breaks,GDSGE_PP.coefs,...');
w.add('        int32(GDSGE_PP.pieces),int32(GDSGE_PP.order),int32(GDSGE_PP.dim),''not-a-knot'',[%s],[],[],[]);', simuStateRows);
w.add('end');
```

with:

```matlab
if isAsg
    w.add('GDSGE_INTERP_RESULTS = GDSGE_ASG_OUTPUT_INTERP.eval_vec(SimuRslt.shock(:,GDSGE_t)'',[%s]);', simuStateRows);
else
    w.add('if shock_num>1');
    w.add('    GDSGE_INTERP_RESULTS = myppual_mex(int32(NumThreads),GDSGE_PP.breaks,GDSGE_PP.coefs,...');
    w.add('        int32(GDSGE_PP.pieces),int32(GDSGE_PP.order),int32(GDSGE_PP.dim),''not-a-knot'',[SimuRslt.shock(:,GDSGE_t)'';%s],[],[],[]);', simuStateRows);
    w.add('else');
    w.add('    GDSGE_INTERP_RESULTS = myppual_mex(int32(NumThreads),GDSGE_PP.breaks,GDSGE_PP.coefs,...');
    w.add('        int32(GDSGE_PP.pieces),int32(GDSGE_PP.order),int32(GDSGE_PP.dim),''not-a-knot'',[%s],[],[],[]);', simuStateRows);
    w.add('end');
end
```

- [ ] **Step 3b: Route ASG + `simuInterp` to the interp emitter.** Replace the whole body of `src/+gdsge/+codegen/generateMatlab.m`:

```matlab
function files = generateMatlab(ir, outDir)
% GENERATEMATLAB  IR -> iter_<model>.m / simulate_<model>.m written to outDir.
%   interpMethod 'asg' -> emitIterAsg + (emitSimulateInterp when
%   options.simuInterp==1, else emitSimulateAsg). Otherwise the spline path:
%   emitIter + (emitSimulateInterp when options.simuInterp==1, else emitSimulate).
if nargin < 2; outDir = pwd; end
files = struct();
files.iterFile = fullfile(outDir, ['iter_' ir.modelName '.m']);
simuInterp = isfield(ir.options, 'simuInterp') && ir.options.simuInterp == 1;
if strcmp(ir.options.interpMethod, 'asg')
    gdsge.codegen.writeText(files.iterFile, gdsge.codegen.mat.emitIterAsg(ir));
    if simuInterp
        simTxt = gdsge.codegen.mat.emitSimulateInterp(ir);
    else
        simTxt = gdsge.codegen.mat.emitSimulateAsg(ir);
    end
else
    gdsge.codegen.writeText(files.iterFile, gdsge.codegen.mat.emitIter(ir));
    if simuInterp
        simTxt = gdsge.codegen.mat.emitSimulateInterp(ir);
    else
        simTxt = gdsge.codegen.mat.emitSimulate(ir);
    end
end
files.simulateFile = fullfile(outDir, ['simulate_' ir.modelName '.m']);
gdsge.codegen.writeText(files.simulateFile, simTxt);
end
```

- [ ] **Step 4: Run the tests to verify they pass**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tEmitSimulateInterp.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: PASS (all five methods, including the existing spline ones).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/+mat/emitSimulateInterp.m src/+gdsge/+codegen/generateMatlab.m tests/codegen/tEmitSimulateInterp.m
git commit -m "feat(codegen): ASG branch in emitSimulateInterp; route ASG+SIMU_INTERP

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: CaoKS2016 + `SIMU_INTERP` fixture, golden capture, integrity test

Build the differential fixture (CaoKS2016 with the simulate mode flipped), capture its golden from the **old** toolbox, and add a golden-integrity test.

**Files:**
- Create: `tests/CaoKS2016_simu_interp/CaoKS2016_simu_interp.gmod`
- Create: `tests/golden/capture_CaoKS2016_simu_interp.m`
- Create: `tests/CaoKS2016_simu_interp/golden/IterRslt.mat`, `tests/CaoKS2016_simu_interp/golden/SimuRslt.mat` (produced by the capture run)
- Create: `tests/CaoKS2016_simu_interp/tGoldenCaoKS2016SimuInterp.m`

- [ ] **Step 1: Create the fixture gmod** — from the repo root, in the Bash tool:

```bash
mkdir -p tests/CaoKS2016_simu_interp/golden tests/CaoKS2016_simu_interp/parser tests/CaoKS2016_simu_interp/codegen tests/CaoKS2016_simu_interp/ir
cp tests/CaoKS2016/CaoKS2016.gmod tests/CaoKS2016_simu_interp/CaoKS2016_simu_interp.gmod
sed -i '/^USE_SPLINE=0;/a SIMU_RESOLVE=0;\nSIMU_INTERP=1;' tests/CaoKS2016_simu_interp/CaoKS2016_simu_interp.gmod
grep -nE "USE_ASG|USE_SPLINE|SIMU_RESOLVE|SIMU_INTERP" tests/CaoKS2016_simu_interp/CaoKS2016_simu_interp.gmod
```
Expected `grep` output: `USE_ASG=1;`, `USE_SPLINE=0;`, then `SIMU_RESOLVE=0;`, `SIMU_INTERP=1;` (each appearing exactly once).

- [ ] **Step 2: Write the golden capture script** — create `tests/golden/capture_CaoKS2016_simu_interp.m`:

```matlab
function capture_CaoKS2016_simu_interp()
% Golden capture for CaoKS2016 with SIMU_INTERP from the OLD toolbox. The iter
% is identical to CaoKS2016 (SIMU_INTERP only changes simulate); simulate
% evaluates asg_output_struct per period (no re-solve), so it is fast. Reduced
% seeded simulate 6x1000 (matches the 7a/7b golden sizes).
here     = fileparts(mfilename('fullpath'));          % tests/golden
repoRoot = fileparts(fileparts(here));
oldSrc   = fullfile(repoRoot, 'base_package', 'gdsge', 'source');
modelDir = fullfile(repoRoot, 'tests', 'CaoKS2016_simu_interp');
goldenDir = fullfile(modelDir, 'golden');

work = tempname; mkdir(work);
copyfile(fullfile(modelDir, 'CaoKS2016_simu_interp.gmod'), work);

oldPath = path; restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(oldSrc);
oldCd = pwd; cdRestore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(work);

t0 = tic;
gdsge_codegen('CaoKS2016_simu_interp');
IterOptions.PrintFreq = 10;
IterOptions.SaveFreq = 100;
IterRslt = iter_CaoKS2016_simu_interp(IterOptions);
fprintf('iter wall-clock: %.1fs (Iter=%d)\n', toc(t0), IterRslt.Iter);

t0 = tic;
rng(0823);
SimuRslt = simulate_CaoKS2016_simu_interp(IterRslt, struct('num_samples', 6, 'num_periods', 1000));
fprintf('simulate wall-clock: %.1fs\n', toc(t0));

save(fullfile(goldenDir, 'IterRslt.mat'), 'IterRslt', '-v7');
save(fullfile(goldenDir, 'SimuRslt.mat'), 'SimuRslt', '-v7');
fprintf('Golden captured: Iter=%d Metric=%g\n', IterRslt.Iter, IterRslt.Metric);
end
```

- [ ] **Step 3: Run the capture (own process, old source only)** — from the repo root:

```
matlab -batch "addpath('tests/golden'); capture_CaoKS2016_simu_interp"
```
Expected: prints `iter wall-clock`, `simulate wall-clock`, and `Golden captured: Iter=... Metric=...` (Iter near 281, Metric < 1e-4). Two `.mat` files now exist under `tests/CaoKS2016_simu_interp/golden/`. Verify their size is small (well under 1 MB):
```bash
ls -la tests/CaoKS2016_simu_interp/golden/
```
If the run errors (old toolbox rejects the combination, etc.), STOP and report — do not fabricate a golden.

- [ ] **Step 4: Write the golden-integrity test** — create `tests/CaoKS2016_simu_interp/tGoldenCaoKS2016SimuInterp.m`:

```matlab
classdef tGoldenCaoKS2016SimuInterp < matlab.unittest.TestCase
    % Phase 7c: the CaoKS2016+SIMU_INTERP golden exists, loads, and has the
    % expected ASG / simulate shape.
    methods (Test)
        function goldenLoadsAndHasShape(tc)
            here = fileparts(mfilename('fullpath'));
            g = load(fullfile(here, 'golden', 'IterRslt.mat'));
            R = g.IterRslt;
            tc.verifyEqual(R.shock_num, 4);
            tc.verifyLessThan(R.Metric, 1e-4);
            tc.verifyTrue(isfield(R, 'asg_output_struct'));
            tc.verifyEqual(R.asg_output_struct.numVec, 4);   % Kp,Xp,kp1,kp2
            tc.verifyEqual(R.output_var_index.Kp, 1);
            s = load(fullfile(here, 'golden', 'SimuRslt.mat'));
            tc.verifyEqual(size(s.SimuRslt.shock), [6, 1001]);
            tc.verifyTrue(isfield(s.SimuRslt, 'kp1'));
            tc.verifyTrue(isfield(s.SimuRslt, 'K'));
        end
    end
end
```

- [ ] **Step 5: Run the integrity test**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/CaoKS2016_simu_interp/tGoldenCaoKS2016SimuInterp.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: PASS. If a field assertion fails, adjust it to the golden's actual shape (the golden is the source of truth) and re-run.

- [ ] **Step 6: Commit**

```bash
git add tests/CaoKS2016_simu_interp/CaoKS2016_simu_interp.gmod tests/golden/capture_CaoKS2016_simu_interp.m tests/CaoKS2016_simu_interp/golden tests/CaoKS2016_simu_interp/tGoldenCaoKS2016SimuInterp.m
git commit -m "test(golden): capture CaoKS2016+SIMU_INTERP golden from old toolbox

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: ASG + `SIMU_INTERP` front-end gate + IR snapshot

Prove the fixture parses to a valid, round-tripping IR with the expected options, and commit the IR JSON snapshot.

**Files:**
- Create: `tests/CaoKS2016_simu_interp/ir/CaoKS2016_simu_interp.gdsge.json` (generated)
- Create: `tests/CaoKS2016_simu_interp/parser/tFrontEndCaoKS2016SimuInterp.m`

- [ ] **Step 1: Generate the IR snapshot** — from the repo root:

```
matlab -batch "addpath('src','src/kernels'); ir=gdsge.parser.parseFrontEnd(fileread('tests/CaoKS2016_simu_interp/CaoKS2016_simu_interp.gmod'),'CaoKS2016_simu_interp'); gdsge.codegen.writeText('tests/CaoKS2016_simu_interp/ir/CaoKS2016_simu_interp.gdsge.json', gdsge.ir.encode(ir)); disp('IR written');"
```
Expected: prints `IR written`; the JSON exists. (This depends on Task 5 — the parser must accept ASG+SIMU_INTERP.)

- [ ] **Step 2: Write the front-end gate** — create `tests/CaoKS2016_simu_interp/parser/tFrontEndCaoKS2016SimuInterp.m`:

```matlab
classdef tFrontEndCaoKS2016SimuInterp < matlab.unittest.TestCase
    % Phase 7c front-end gate: ASG model with SIMU_INTERP parses, validates,
    % round-trips, and matches the committed IR snapshot.
    methods (Test)
        function frontEndProducesValidIR(tc)
            ir = parseIt();
            tc.verifyTrue(gdsge.ir.validate(ir).pass);
            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));

            tc.verifyEqual(ir.options.interpMethod, 'asg');
            tc.verifyEqual(ir.options.simuInterp, 1);
            tc.verifyEqual(ir.options.simuResolve, 0);
            tc.verifyEqual(ir.states.names, {'K','X'});
            tc.verifyEqual(ir.shocks.count, 4);
            tc.verifyEqual(ir.variables.output, {'Kp','Xp','kp1','kp2'});
        end
        function encodingMatchesSnapshot(tc)
            ir = parseIt();
            here = fileparts(mfilename('fullpath'));
            onDisk = fileread(fullfile(fileparts(here), 'ir', 'CaoKS2016_simu_interp.gdsge.json'));
            tc.verifyEqual(norm_(gdsge.ir.encode(ir)), norm_(onDisk), ...
                'CaoKS2016_simu_interp.gdsge.json is stale — regenerate it.');
        end
    end
end

function ir = parseIt()
here = fileparts(mfilename('fullpath'));
gmodPath = fullfile(fileparts(here), 'CaoKS2016_simu_interp.gmod');
ir = gdsge.parser.parseFrontEnd(fileread(gmodPath), 'CaoKS2016_simu_interp');
end

function s = norm_(s)
s = strtrim(strrep(s, sprintf('\r\n'), sprintf('\n')));
end
```

- [ ] **Step 3: Run the front-end gate**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/CaoKS2016_simu_interp/parser/tFrontEndCaoKS2016SimuInterp.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add tests/CaoKS2016_simu_interp/ir/CaoKS2016_simu_interp.gdsge.json tests/CaoKS2016_simu_interp/parser/tFrontEndCaoKS2016SimuInterp.m
git commit -m "test(parser): CaoKS2016+SIMU_INTERP front-end gate + IR snapshot

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 9: ASG + `SIMU_INTERP` end-to-end gate (Slow)

Drive the public API end-to-end and match the captured golden through the interp simulate path.

**Files:**
- Create: `tests/CaoKS2016_simu_interp/codegen/tEndToEndCaoKS2016SimuInterp.m`

- [ ] **Step 1: Write the end-to-end gate** — create `tests/CaoKS2016_simu_interp/codegen/tEndToEndCaoKS2016SimuInterp.m`:

```matlab
classdef tEndToEndCaoKS2016SimuInterp < matlab.unittest.TestCase
    % PHASE 7c GATE: ASG model with SIMU_INTERP, end-to-end through the public
    % API vs the freshly captured golden. iter is the standard ASG solve; the
    % simulate evaluates asg_output_struct per period (no re-solve). Slow.
    properties (Constant)
        RelTol = 1e-3;
        AbsTol = 1e-3;
    end
    methods (Test, TestTags = {'Slow'})
        function publicApiMatchesGolden(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'CaoKS2016_simu_interp.gmod'), work);

            ir = gdsge_codegen('CaoKS2016_simu_interp');

            mexFile = fullfile(work, ['mex_CaoKS2016_simu_interp.' mexext]);
            artifacts = {'iter_CaoKS2016_simu_interp.m', ...
                'simulate_CaoKS2016_simu_interp.m', 'mex_CaoKS2016_simu_interp.cpp', ...
                'compile_CaoKS2016_simu_interp.m', 'CaoKS2016_simu_interp.gdsge.json', ...
                'mex_CaoKS2016_simu_interp.cache'};
            for i = 1:numel(artifacts)
                tc.assertTrue(exist(fullfile(work, artifacts{i}), 'file') == 2, artifacts{i});
            end
            tc.assertTrue(exist(mexFile, 'file') == 3, 'MEX did not compile');
            decoded = gdsge.ir.decode(fileread(fullfile(work, 'CaoKS2016_simu_interp.gdsge.json')));
            tc.verifyTrue(gdsge.ir.isequalIR(decoded, ir), 'IR JSON does not round-trip');

            % ---- the generated simulate uses the interp path, not a re-solve
            simTxt = fileread(fullfile(work, 'simulate_CaoKS2016_simu_interp.m'));
            tc.verifyTrue(contains(simTxt, 'asg_output_struct'));
            tc.verifyFalse(contains(simTxt, 'solveProblems'));

            % ---- iter (test.m parity)
            IterOptions.PrintFreq = 10;
            IterOptions.SaveFreq = 100;
            IterRslt = iter_CaoKS2016_simu_interp(IterOptions);

            golden = load(fullfile(modelDir, 'golden', 'IterRslt.mat'));
            G = golden.IterRslt;
            tc.verifyLessThan(IterRslt.Metric, 1e-4);
            tc.verifyLessThan(abs(IterRslt.Iter - G.Iter), 0.2*G.Iter + 20, ...
                sprintf('Iter=%d vs golden %d', IterRslt.Iter, G.Iter));
            tc.verifyTrue(isfield(IterRslt, 'asg_output_struct'));
            tc.verifyEqual(IterRslt.output_var_index, G.output_var_index);

            % ---- seeded reduced simulate via the interp path (capture protocol)
            rng(0823);
            SimuRslt = simulate_CaoKS2016_simu_interp(IterRslt, ...
                struct('num_samples', 6, 'num_periods', 1000));
            gs = load(fullfile(modelDir, 'golden', 'SimuRslt.mat'));
            GS = gs.SimuRslt;
            tc.verifyEqual(SimuRslt.shock, GS.shock, ...
                'shock paths differ — check shock_trans bit-identity first');
            flds = {'K','X','kp1','kp2'};
            T0 = 100;
            for i = 1:numel(flds)
                a = SimuRslt.(flds{i}); b = GS.(flds{i});
                tc.verifyEqual(size(a), size(b), flds{i});
                r = gdsgetest.compareNumericClose(a(:,1:T0), b(:,1:T0), 1e-2, 1e-2);
                tc.verifyTrue(r.pass, sprintf('%s early path: %s', flds{i}, ...
                    strjoin(r.failures, newline)));
            end
        end
    end
end
```

- [ ] **Step 2: Run the gate to verify it passes**

Run:
```
matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/CaoKS2016_simu_interp/codegen/tEndToEndCaoKS2016SimuInterp.m'); disp(r); assert(~any([r.Failed]))"
```
Expected: PASS (Slow — compiles the ASG model MEX, runs the ~281-iter solve, then a fast interp simulate). If `flds` includes a field the golden lacks, replace `flds` with the var_output/state fields actually present in the golden's `SimuRslt` (inspect `fieldnames(GS)`), then re-run.

- [ ] **Step 3: Commit**

```bash
git add tests/CaoKS2016_simu_interp/codegen/tEndToEndCaoKS2016SimuInterp.m
git commit -m "test(gate): CaoKS2016+SIMU_INTERP green end-to-end vs golden

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 10: Full suite, PROGRESS update, finalize

- [ ] **Step 1: Run the full suite**

Run:
```
matlab -batch "cd('tests'); run_tests"
```
Expected: exit code 0. Confirm in `tests/results/junit.xml` that there are no failures and the new classes ran (`tGdsgeDriverHL1996`, `tEmitSimulateInterp` ASG methods, `tResolveOptions.asgAcceptsSimuInterp`, `tFrontEndCaoKS2016SimuInterp`, `tEndToEndCaoKS2016SimuInterp`, `tGoldenCaoKS2016SimuInterp`, `tCodegenDriver` new methods).

- [ ] **Step 2: Update `PROGRESS.md`** — mark Phase 7c done (change `☐` to `☑` and remove `**NEXT**` on the 7c line; mark 7d as the next phase), and add a changelog entry at the top of the Changelog section:

```markdown
- 2026-06-13: **Phase 7c complete.** Legacy orchestration surface + ASG-interp
  simulate. `gdsge.codegen.codegen` gained a 2nd output (generated code strings +
  pre-overwrite cache); the flat `gdsge_codegen` re-exposes the legacy 5 outputs
  `[model,iterCode,cppCache,cppCode,codeSegment]` (single-output callers
  unaffected); `codeSegment` carries the real generated strings (no old template
  fragments). New flat `gdsge.m` orchestrator: codegen -> cache-gated solve
  (IterRslt_<model>.mat + iter_<model>.cache) -> simulate -> eq{model,IterRslt,
  SmltRslt}, no v2struct; HL1996 gate proves the golden match and the second-run
  skip-solve. `options.GenCodeSegment` now raises an architecture-honest error
  (gdsge:codegen:genCodeSegmentUnsupported) — the thin-file architecture has no
  segment decomposition. ASG+SIMU_INTERP: removed the 7b parser guard; extended
  emitSimulateInterp with an ASG branch (asg.construct_from_struct(
  asg_output_struct) + per-period eval_vec, no re-solve) and routed it via
  generateMatlab; validated differentially against a fresh old-toolbox golden
  (CaoKS2016 with SIMU_INTERP flipped: own fixture/capture/integrity/front-end/
  end-to-end gates). Full suite green. All on branch phase7c-legacy-orchestration.
  Next: Phase 7d (macro engine).
```

- [ ] **Step 3: Commit the PROGRESS update**

```bash
git add PROGRESS.md
git commit -m "docs(progress): Phase 7c complete — legacy orchestration + ASG-interp

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 4: Finish the branch** — invoke the `superpowers:finishing-a-development-branch` skill to decide merge/PR for `phase7c-legacy-orchestration` (the 7a/7b precedent is a merge into `main`).

---

## Notes & risks (carried from the spec)

- **gdsge.m rng:** `gdsge()` does not seed before simulating, and the iter solve's randomized restarts can consume rng — so the `tGdsgeDriverHL1996` gate checks `eq.IterRslt` against the golden (deterministic) and `eq.SmltRslt` structurally, not bit-for-bit. Bit-exact simulate is already covered by `tEndToEndHL1996`.
- **`codeSegment` shape:** named after the real generated artifacts (`iterCode`/`simulateCode`/`cppCode`/`compileCode`), deliberately not the old five fragments. No corpus caller inspects it.
- **ASG-interp fixture:** flipping CaoKS2016 to `SIMU_INTERP` is a synthetic configuration; it is a pure evaluate-the-output-interpolant path (no re-solve), differential-tested against the old toolbox running the identical flipped gmod.
- **If the golden capture is prohibitively long or the old toolbox rejects the combination** (Task 7 Step 3): stop and revisit the protocol — do not fabricate a golden.
```
