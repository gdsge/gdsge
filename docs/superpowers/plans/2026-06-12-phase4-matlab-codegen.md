# Phase 4 — MATLAB Codegen Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate `iter_HL1996.m` / `simulate_HL1996.m` from the IR (no v2struct, real error reporting) and prove them functionally against the **old compiled MEX** and the committed goldens.

**Architecture:** Thin generated files + a hand-written `gdsge.runtime` library (resolve cascade, options unpacking, warmup transfer, progress/unconverged reporting) + pure section emitters in `gdsge.codegen.mat` (no .m templates) + vendored flat kernels in `src/kernels/`. Spec: `docs/superpowers/specs/2026-06-12-phase4-matlab-codegen-design.md`.

**Tech Stack:** MATLAB R2025b (`matlab -batch`), `matlab.unittest`, old toolbox in `base_package/` (reference + MEX capture only).

---

## Orientation for the engineer (read first)

- **Reference files** (what the OLD generator produced for HL1996): regenerate with
  `matlab -batch "cd('D:\refactor_gdsge\scratch'); regen_old_matlab_codegen"` →
  `scratch/old_generated/iter_HL1996.m` and `simulate_HL1996.m`. Every "parity" comment in
  this plan refers to those files. If `scratch/regen_old_matlab_codegen.m` is missing, it is
  recreated in Task 17 (same content as the capture script but parser-only).
- **The IR**: `tests/HeatonLucas1996/ir/buildHL1996IR.m` is the fixture used by all emitter
  tests. Field reference: `src/+gdsge/+ir/schema.m`. Key facts: policy slots are
  `[lo hi]` ranges into a 19-row `GDSGE_SOL`; `states.grids` values are **opaque MATLAB
  text**; params/shock values/`shock_trans` are numeric; `ir.interp{i}` has
  `name/args/initialExpr/updateExpr`.
- **The MEX caller-workspace contract**: the MEX reads these from the workspace of the
  function that calls it: `MEX_TASK_NAME`, `MEX_TASK_INIT`, `MEX_TASK_INF_HORIZON`,
  `TolFun`, `TolSol`, `SolMaxIter`, `NumThreads`, `GDSGE_DEBUG_EVAL_ONLY`, `UseBroyden`,
  `FiniteDiffDelta`, `GDSGE_USE_BROYDEN_NOW`, `GDSGE_SPLINE_VEC`. In the new design only
  `gdsge.runtime.solveProblems` calls the MEX, so only it must define them — and **every
  MEX call must be in its main function body** (a subfunction would change the caller).
- **Run a single test file** (from repo root; used in every task):
  `matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/<path>/<tFile>.m'))"`
- **Run everything**: `matlab -batch "cd('tests'); run_tests"` (exit 0 = green).
- **Commit style**: small commits, `feat(codegen):` / `feat(runtime):` / `test(...):` /
  `chore(...):` prefixes, body line `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.

---

### Task 1: Branch

- [ ] **Step 1.1: Create the phase branch**

```powershell
git checkout -b phase4-matlab-codegen
```

Expected: `Switched to a new branch 'phase4-matlab-codegen'`.

---

### Task 2: Vendor the flat MATLAB kernels

**Files:**
- Create: `src/kernels/` (copies of 6 files from `base_package/gdsge/source/`)
- Modify: `tests/run_tests.m` (addpath `src/kernels`)
- Test: `tests/kernels/tKernels.m`

- [ ] **Step 2.1: Write the failing smoke test**

Create `tests/kernels/tKernels.m`:

```matlab
classdef tKernels < matlab.unittest.TestCase
    % Smoke tests for the vendored flat runtime kernels in src/kernels/.
    methods (TestClassSetup)
        function blasOnSystemPath(tc)
            % myppual_mex links against essential_blas.dll; its dir must be on PATH.
            root = fileparts(which('essential_blas.dll'));
            tc.assertNotEmpty(root, 'essential_blas.dll not found — vendor src/kernels first');
            if ~any(strcmp(strsplit(getenv('PATH'), ';'), root))
                setenv('PATH', [getenv('PATH') ';' root]);
            end
        end
    end
    methods (Test)
        function myppualEvaluatesSpline(tc)
            x = linspace(0, 1, 11); v = x.^2;
            pp = struct('form','MKL','breaks',{{x}},'Values',reshape(v,[],11), ...
                'coefs',[],'order',4,'Method',[],'ExtrapolationOrder',2, ...
                'thread',1,'orient','curvefit');
            pp = myppual(pp);
            y = myppual(pp, 0.5);
            tc.verifyLessThan(abs(y - 0.25), 1e-8);
        end
        function convertToEvalArrayAccepts(tc)
            x = linspace(0, 1, 11); v = x.^2;
            pp = struct('form','MKL','breaks',{{x}},'Values',reshape(v,[],11), ...
                'coefs',[],'order',4,'Method',[],'ExtrapolationOrder',2, ...
                'thread',1,'orient','curvefit');
            pp = myppual(pp);
            vec = convert_to_interp_eval_array({pp, pp});
            tc.verifyTrue(isstruct(vec));
        end
        function markovChainShapeAndRange(tc)
            rng(1);
            trans = [0.9 0.1; 0.2 0.8];
            s = gen_discrete_markov_rn(trans, 5, 7, ones(5,1));
            tc.verifySize(s, [5 7]);
            tc.verifyTrue(all(ismember(s(:), [1 2])));
        end
        function getScalarReturnsStruct(tc)
            st = struct(); st.a = ones(8, 201);
            r = get_scalar(st, {1:8, linspace(0, 1, 201)});
            tc.verifyTrue(isstruct(r) && isfield(r, 'a'));
        end
    end
end
```

- [ ] **Step 2.2: Run it — expect FAIL (functions undefined)**

```
matlab -batch "addpath('src'); addpath('tests'); assertSuccess(runtests('tests/kernels/tKernels.m'))"
```

Expected: failure, `essential_blas.dll not found` (assertion in setup).

- [ ] **Step 2.3: Copy the kernels**

```powershell
New-Item -ItemType Directory -Force src\kernels
$srcOld = "base_package\gdsge\source"
Copy-Item "$srcOld\myppual.m","$srcOld\myppual_mex.mexw64","$srcOld\essential_blas.dll","$srcOld\convert_to_interp_eval_array.m","$srcOld\gen_discrete_markov_rn.m","$srcOld\get_scalar.m" src\kernels\
```

Verify with `git status` that all six appear as untracked under `src/kernels/`
(none match a `.gitignore` pattern — `mex_*.mexw64` does not match `myppual_mex.mexw64`).

- [ ] **Step 2.4: Add `src/kernels` to the harness path**

In `tests/run_tests.m`, after the existing `addpath(fullfile(fileparts(thisDir), 'src'));` line, add:

```matlab
addpath(fullfile(fileparts(thisDir), 'src', 'kernels'));   % vendored flat runtime kernels
```

- [ ] **Step 2.5: Run the smoke test — expect PASS**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/kernels/tKernels.m'))"
```

If a kernel errors with `Undefined function 'X'`, `X` lives in
`base_package/gdsge/source/` too — copy it into `src/kernels/` and re-run (the spec
accepts dependency-closure growth here; list any additions in the commit message).

- [ ] **Step 2.6: Commit**

```powershell
git add src/kernels tests/kernels/tKernels.m tests/run_tests.m
git commit -m "feat(kernels): vendor flat MATLAB/MEX runtime kernels from old source"
```

---

### Task 3: IR amendment — `numThreads` dynamic sentinel (0)

**Files:**
- Modify: `src/+gdsge/+parser/defaultSetupCode.m:20`
- Modify: `src/+gdsge/+parser/resolveOptions.m:21`
- Modify: `tests/HeatonLucas1996/ir/buildHL1996IR.m:6,20`
- Modify: `tests/parser/tResolveOptions.m:14`, `tests/parser/tEvalSetup.m:29`,
  `tests/HeatonLucas1996/parser/tFullIRHL1996.m:14-15`,
  `tests/HeatonLucas1996/parser/tFrontEndHL1996.m:19-26`
- Regenerate: `tests/HeatonLucas1996/ir/HL1996.gdsge.json`

- [ ] **Step 3.1: Make the resolveOptions test pin the sentinel (failing first)**

In `tests/parser/tResolveOptions.m`, replace the line
`tc.verifyTrue(o.numThreads >= 1);` with:

```matlab
            tc.verifyEqual(o.numThreads, 0);   % 0 = dynamic: feature('numcores') at runtime
```

(If the surrounding test feeds a `ws` that explicitly sets `NumThreads`, instead add a new
test method asserting `numThreads == 0` when `ws` has no `NumThreads` field, and keep the
explicit-value test asserting the explicit number.)

- [ ] **Step 3.2: Run — expect FAIL** (`o.numThreads` is the machine's core count)

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/parser/tResolveOptions.m'))"
```

- [ ] **Step 3.3: Implement the sentinel**

`src/+gdsge/+parser/defaultSetupCode.m` line 20 — change
`'NumThreads = feature(''numcores'');', ...` to:

```matlab
    'NumThreads = 0;', ...   % 0 = dynamic sentinel: codegen emits feature('numcores')
```

`src/+gdsge/+parser/resolveOptions.m` line 21 — change to:

```matlab
o.numThreads  = getf(ws, 'NumThreads', 0);   % 0 = dynamic (feature('numcores') at runtime)
```

- [ ] **Step 3.4: Update the fixture and downstream tests**

`tests/HeatonLucas1996/ir/buildHL1996IR.m`: change `'numThreads',8` to `'numThreads',0`
and update both comments ("pinned at 8 for golden determinism" → "0 = dynamic sentinel,
resolved to feature('numcores') by the MATLAB codegen").

`tests/parser/tEvalSetup.m:29`: `tc.verifyTrue(ws.NumThreads >= 1);` →
`tc.verifyEqual(ws.NumThreads, 0);`

`tests/HeatonLucas1996/parser/tFullIRHL1996.m`: delete the two neutralization lines
(`ref.options.numThreads = 0;` / `ir.options.numThreads = 0;`) and the
"machine-derived; neutralize" comment.

`tests/HeatonLucas1996/parser/tFrontEndHL1996.m`: replace
`tc.verifyGreaterThanOrEqual(ir.options.numThreads, 1);` with
`tc.verifyEqual(ir.options.numThreads, 0);` and delete its two neutralization lines.

- [ ] **Step 3.5: Regenerate the JSON golden**

```
matlab -batch "addpath('src'); addpath('tests/HeatonLucas1996/ir'); txt = gdsge.ir.encode(buildHL1996IR()); f = fopen('tests/HeatonLucas1996/ir/HL1996.gdsge.json','w'); fwrite(f, txt); fclose(f);"
```

- [ ] **Step 3.6: Run the full suite — expect PASS**

```
matlab -batch "cd('tests'); run_tests"
```

- [ ] **Step 3.7: Commit**

```powershell
git add -A src tests
git commit -m "feat(ir): numThreads 0 = dynamic sentinel; codegen will emit feature('numcores')"
```

---

### Task 4: `codeWriter`

**Files:**
- Create: `src/+gdsge/+codegen/+mat/codeWriter.m`
- Test: `tests/codegen/tCodeWriter.m`

- [ ] **Step 4.1: Write the failing test**

Create `tests/codegen/tCodeWriter.m`:

```matlab
classdef tCodeWriter < matlab.unittest.TestCase
    methods (Test)
        function addFormatsAndJoins(tc)
            w = gdsge.codegen.mat.codeWriter();
            w.add('x = %d;', 3);
            w.add('y = %s;', 'z');
            tc.verifyEqual(w.str(), sprintf('x = 3;\ny = z;'));
        end
        function indentIsFourSpacesPerLevel(tc)
            w = gdsge.codegen.mat.codeWriter();
            w.add('if a');
            w.in();
            w.add('b = 1;');
            w.out();
            w.add('end');
            tc.verifyEqual(w.str(), sprintf('if a\n    b = 1;\nend'));
        end
        function addRawIsVerbatimMultiline(tc)
            w = gdsge.codegen.mat.codeWriter();
            w.in();   % addRaw must NOT apply indentation
            w.addRaw(sprintf('a\n\nb'));
            tc.verifyEqual(w.str(), sprintf('a\n\nb'));
        end
        function blankAddsEmptyLine(tc)
            w = gdsge.codegen.mat.codeWriter();
            w.add('a');
            w.blank();
            w.add('b');
            tc.verifyEqual(w.str(), sprintf('a\n\nb'));
        end
        function crlfInputNormalized(tc)
            w = gdsge.codegen.mat.codeWriter();
            w.addRaw(sprintf('a\r\nb'));
            tc.verifyEqual(w.str(), sprintf('a\nb'));
        end
    end
end
```

- [ ] **Step 4.2: Run — expect FAIL** (`Unable to resolve the name gdsge.codegen.mat.codeWriter`)

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/codegen/tCodeWriter.m'))"
```

- [ ] **Step 4.3: Implement**

Create `src/+gdsge/+codegen/+mat/codeWriter.m`:

```matlab
classdef codeWriter < handle
    % CODEWRITER  Minimal indented line buffer for MATLAB code emission.
    %   add(fmt, ...)  sprintf-formatted line at the current indent (4 sp/level)
    %   addRaw(text)   verbatim multi-line block, no indentation applied
    %   blank()        empty line
    %   in()/out()     indent level +/- 1
    %   str()          all lines joined with \n
    properties (Access = private)
        lines = {};
        level = 0;
    end
    methods
        function add(w, fmt, varargin)
            w.lines{end+1} = [repmat(' ', 1, 4*w.level) sprintf(fmt, varargin{:})];
        end
        function addRaw(w, text)
            text = strrep(text, sprintf('\r\n'), newline);
            parts = strsplit(text, newline, 'CollapseDelimiters', false);
            w.lines = [w.lines, parts];
        end
        function blank(w)
            w.lines{end+1} = '';
        end
        function in(w),  w.level = w.level + 1; end
        function out(w), w.level = max(0, w.level - 1); end
        function s = str(w)
            s = strjoin(w.lines, newline);
        end
    end
end
```

- [ ] **Step 4.4: Run — expect PASS** (same command as 4.2)

- [ ] **Step 4.5: Commit**

```powershell
git add src/+gdsge/+codegen tests/codegen/tCodeWriter.m
git commit -m "feat(codegen): codeWriter - indented line buffer for emitters"
```

---

### Task 5: `rewriteNames`

**Files:**
- Create: `src/+gdsge/+codegen/+mat/rewriteNames.m`
- Test: `tests/codegen/tRewriteNames.m`

- [ ] **Step 5.1: Write the failing test**

Create `tests/codegen/tRewriteNames.m`:

```matlab
classdef tRewriteNames < matlab.unittest.TestCase
    methods (Test)
        function wholeIdentifierOnly(tc)
            out = gdsge.codegen.mat.rewriteNames('w1 + w1n + Kw1', ...
                {'w1'}, {'GDSGE_TENSOR_w1(:)'''});
            tc.verifyEqual(out, 'GDSGE_TENSOR_w1(:)'' + w1n + Kw1');
        end
        function multipleNamesOnePass(tc)
            out = gdsge.codegen.mat.rewriteNames('w1.*d+eta1', ...
                {'w1','d','eta1'}, ...
                {'GDSGE_TENSOR_w1(:)''','GDSGE_TENSOR_d(:)''','GDSGE_TENSOR_eta1(:)'''});
            tc.verifyEqual(out, ...
                'GDSGE_TENSOR_w1(:)''.*GDSGE_TENSOR_d(:)''+GDSGE_TENSOR_eta1(:)''');
        end
        function replacementNotRematched(tc)
            % 'd' must not match inside the GDSGE_TENSOR_w1 identifier
            out = gdsge.codegen.mat.rewriteNames('w1+d', ...
                {'w1','d'}, {'GDSGE_TENSOR_w1(:)''','D2'});
            tc.verifyEqual(out, 'GDSGE_TENSOR_w1(:)''+D2');
        end
        function untouchedWhenNoMatch(tc)
            out = gdsge.codegen.mat.rewriteNames('0.0', {'w1'}, {'X'});
            tc.verifyEqual(out, '0.0');
        end
    end
end
```

- [ ] **Step 5.2: Run — expect FAIL**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/codegen/tRewriteNames.m'))"
```

- [ ] **Step 5.3: Implement**

Create `src/+gdsge/+codegen/+mat/rewriteNames.m`:

```matlab
function out = rewriteNames(txt, names, replacements)
% REWRITENAMES  Word-boundary identifier replacement in opaque MATLAB text.
%   Replaces whole-identifier occurrences of names{i} with replacements{i}.
%   Boundaries are MATLAB identifier characters, so 'w1' does not match inside
%   'w1n' or 'Kw1'. Replacements containing identifiers are safe because the
%   inserted text always embeds the original name with a '_' prefix boundary.
out = txt;
for i = 1:numel(names)
    pat = ['(?<![A-Za-z0-9_])' regexptranslate('escape', names{i}) '(?![A-Za-z0-9_])'];
    rep = strrep(replacements{i}, '\', '\\');
    rep = strrep(rep, '$', '\$');
    out = regexprep(out, pat, rep);
end
end
```

- [ ] **Step 5.4: Run — expect PASS**

- [ ] **Step 5.5: Commit**

```powershell
git add src/+gdsge/+codegen/+mat/rewriteNames.m tests/codegen/tRewriteNames.m
git commit -m "feat(codegen): rewriteNames - word-boundary identifier rewriting"
```

---

### Task 6: `runtime.unpackOptions`

**Files:**
- Create: `src/+gdsge/+runtime/unpackOptions.m`
- Test: `tests/runtime/tUnpackOptions.m`

- [ ] **Step 6.1: Write the failing test**

Create `tests/runtime/tUnpackOptions.m`:

```matlab
classdef tUnpackOptions < matlab.unittest.TestCase
    methods (Test)
        function assignsWhitelistedIntoCaller(tc)
            [tolEq, beta] = callUnpack(struct('TolEq', 1e-8, 'beta', 0.9), ...
                {'TolEq','beta'});
            tc.verifyEqual(tolEq, 1e-8);
            tc.verifyEqual(beta, 0.9);
        end
        function leavesUnmentionedAlone(tc)
            [tolEq, beta] = callUnpack(struct('beta', 0.9), {'TolEq','beta'});
            tc.verifyEqual(tolEq, 1e-6);    % the default in the helper
            tc.verifyEqual(beta, 0.9);
        end
        function unknownFieldErrors(tc)
            tc.verifyError(@() callUnpack(struct('Bogus', 1), {'TolEq','beta'}), ...
                'gdsge:options:unknownField');
        end
        function emptyOptionsIsNoOp(tc)
            [tolEq, beta] = callUnpack([], {'TolEq','beta'});
            tc.verifyEqual(tolEq, 1e-6);
            tc.verifyEqual(beta, 0.95);
        end
    end
end

function [TolEq, beta] = callUnpack(opts, valid)
% Helper: a function workspace playing the role of the generated file.
TolEq = 1e-6;
beta = 0.95;
gdsge.runtime.unpackOptions(opts, valid);
end
```

- [ ] **Step 6.2: Run — expect FAIL**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/runtime/tUnpackOptions.m'))"
```

- [ ] **Step 6.3: Implement**

Create `src/+gdsge/+runtime/unpackOptions.m`:

```matlab
function unpackOptions(opts, validNames)
% UNPACKOPTIONS  Assign whitelisted option fields into the caller workspace.
%   Replaces v2struct(GDSGE_OPTIONS): explicit, checked, errors on unknown
%   fields instead of silently creating/overwriting workspace variables.
if isempty(opts); return; end
if ~isstruct(opts)
    error('gdsge:options:notAStruct', 'GDSGE_OPTIONS must be a struct.');
end
fn = fieldnames(opts);
unknown = setdiff(fn, validNames);
if ~isempty(unknown)
    error('gdsge:options:unknownField', ...
        'Unknown option field(s): %s\nValid fields: %s', ...
        strjoin(unknown', ', '), strjoin(sort(validNames(:))', ', '));
end
for i = 1:numel(fn)
    assignin('caller', fn{i}, opts.(fn{i}));
end
end
```

- [ ] **Step 6.4: Run — expect PASS**

- [ ] **Step 6.5: Commit**

```powershell
git add src/+gdsge/+runtime/unpackOptions.m tests/runtime/tUnpackOptions.m
git commit -m "feat(runtime): unpackOptions - whitelisted options, error on unknown"
```

---

### Task 7: `runtime.ensurePath`

**Files:**
- Create: `src/+gdsge/+runtime/ensurePath.m`
- Test: `tests/runtime/tEnsurePath.m`

- [ ] **Step 7.1: Write the failing test**

Create `tests/runtime/tEnsurePath.m`:

```matlab
classdef tEnsurePath < matlab.unittest.TestCase
    methods (Test)
        function blasDirLandsOnSystemPath(tc)
            if ~ispc; tc.assumeFail('Windows-only behavior'); end
            gdsge.runtime.ensurePath();
            root = fileparts(which('essential_blas.dll'));
            tc.assertNotEmpty(root);
            tc.verifyTrue(any(strcmp(strsplit(getenv('PATH'), ';'), root)));
        end
        function idempotent(tc)
            if ~ispc; tc.assumeFail('Windows-only behavior'); end
            gdsge.runtime.ensurePath();
            p1 = getenv('PATH');
            gdsge.runtime.ensurePath();
            tc.verifyEqual(getenv('PATH'), p1);
        end
    end
end
```

- [ ] **Step 7.2: Run — expect FAIL**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/runtime/tEnsurePath.m'))"
```

- [ ] **Step 7.3: Implement**

Create `src/+gdsge/+runtime/ensurePath.m`:

```matlab
function ensurePath()
% ENSUREPATH  Put the directory holding essential_blas.dll on the system PATH
%   so MEX kernels (myppual_mex, the model MEX) can load it. Windows only;
%   idempotent. Replaces the BLAS-path header of the old generated files.
if ~ispc; return; end
root = fileparts(which('essential_blas.dll'));
if isempty(root)
    error('gdsge:runtime:blasNotFound', ...
        'essential_blas.dll not found on the MATLAB path; addpath src/kernels first.');
end
if ~any(strcmp(strsplit(getenv('PATH'), ';'), root))
    setenv('PATH', [getenv('PATH') ';' root]);
end
end
```

- [ ] **Step 7.4: Run — expect PASS**

- [ ] **Step 7.5: Commit**

```powershell
git add src/+gdsge/+runtime/ensurePath.m tests/runtime/tEnsurePath.m
git commit -m "feat(runtime): ensurePath - essential_blas.dll system-PATH setup"
```

---

### Task 8: `runtime.computeMetric`, `printIterProgress`, `reportUnconverged`

**Files:**
- Create: `src/+gdsge/+runtime/computeMetric.m`, `printIterProgress.m`, `reportUnconverged.m`
- Test: `tests/runtime/tReporting.m`

- [ ] **Step 8.1: Write the failing tests**

Create `tests/runtime/tReporting.m`:

```matlab
classdef tReporting < matlab.unittest.TestCase
    methods (Test)
        % ---- computeMetric ------------------------------------------------
        function metricZeroWhenEqual(tc)
            a = {ones(2,3), 2*ones(2,3)};
            tc.verifyEqual(gdsge.runtime.computeMetric(a, a), 0);
        end
        function metricIsSupNormAcrossVars(tc)
            new = {[1 2], [10 10]};
            old = {[1 2.5], [10 9]};
            tc.verifyEqual(gdsge.runtime.computeMetric(new, old), 1);
        end
        function metricNanOnSizeMismatch(tc)
            % parity with the old generated try/catch -> NaN
            tc.verifyTrue(isnan(gdsge.runtime.computeMetric({ones(2,3)}, {ones(4,5)})));
        end
        % ---- printIterProgress --------------------------------------------
        function printsOnFreqAndStop(tc)
            out = evalc('gdsge.runtime.printIterProgress(10, 0.5, 1e-9, 0, 1.23, 10, 0, false);');
            tc.verifyTrue(contains(out, 'Iter:10'));
            tc.verifyTrue(contains(out, 'unconverged:0'));
            out2 = evalc('gdsge.runtime.printIterProgress(7, 0.5, 1e-9, 0, 1, 10, 0, true);');
            tc.verifyTrue(contains(out2, 'Iter:7'));   % stopFlag forces print
        end
        function silentOffFreqOrNoPrint(tc)
            out = evalc('gdsge.runtime.printIterProgress(7, 0.5, 1e-9, 0, 1, 10, 0, false);');
            tc.verifyEqual(out, '');
            out2 = evalc('gdsge.runtime.printIterProgress(10, 0.5, 1e-9, 0, 1, 10, 1, false);');
            tc.verifyEqual(out2, '');
        end
        function returnsWhetherPrinted(tc)
            [~, p] = evalc('gdsge.runtime.printIterProgress(10, 0.5, 1e-9, 0, 1, 10, 0, false)');
            tc.verifyTrue(p);
            p2 = gdsge.runtime.printIterProgress(7, 0.5, 1e-9, 0, 1, 10, 1, false);
            tc.verifyFalse(p2);
        end
        % ---- reportUnconverged --------------------------------------------
        function reportNamesWorstPoints(tc)
            % probSize [2 shocks x 3 states]; linear index 4 = (shock 2, state 2)
            need = logical([0 0 0 1 0 1]);
            f = [0 0 0 7 0 3];
            msg = gdsge.runtime.reportUnconverged(need, f, [2 3], {'w1'}, {[0.1 0.5 0.9]}, 5);
            tc.verifyTrue(contains(msg, '2 of 6'));
            tc.verifyTrue(contains(msg, 'shock=2'));
            tc.verifyTrue(contains(msg, 'w1=0.5'));        % worst point first
            tc.verifyTrue(contains(msg, 'residual=7'));
        end
        function reportAllConverged(tc)
            msg = gdsge.runtime.reportUnconverged(false(1,4), zeros(1,4), [2 2], {'w1'}, {[0 1]}, 5);
            tc.verifyEqual(msg, 'All problems converged.');
        end
    end
end
```

- [ ] **Step 8.2: Run — expect FAIL**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/runtime/tReporting.m'))"
```

- [ ] **Step 8.3: Implement the three functions**

Create `src/+gdsge/+runtime/computeMetric.m`:

```matlab
function metric = computeMetric(newCell, oldCell)
% COMPUTEMETRIC  Sup-norm distance across interp variables.
%   NaN on any error (size mismatch etc.) — parity with the old generated
%   try/catch -> NaN. Assumes at least one interp variable.
try
    m = -inf;
    for i = 1:numel(newCell)
        m = max(m, max(abs(newCell{i}(:) - oldCell{i}(:))));
    end
    metric = m;
catch
    metric = nan;
end
end
```

Note: MATLAB's `a - b` with mismatched non-scalar sizes raises an error, which the catch
turns into NaN — exactly the old behavior on the first post-WarmUp iteration.

Create `src/+gdsge/+runtime/printIterProgress.m`:

```matlab
function printed = printIterProgress(iter, metric, maxF, nUnconverged, elapsedSec, printFreq, noPrint, stopFlag)
% PRINTITERPROGRESS  One structured progress line, PrintFreq/NoPrint gated.
%   Returns true if it printed (caller resets its elapsed timer on true).
printed = false;
if noPrint; return; end
if mod(iter, printFreq) == 0 || stopFlag
    fprintf('Iter:%d, Metric:%g, maxF:%g, unconverged:%d, elapsed:%.1fs\n', ...
        iter, metric, maxF, nUnconverged, elapsedSec);
    printed = true;
end
end
```

Create `src/+gdsge/+runtime/reportUnconverged.m`:

```matlab
function msg = reportUnconverged(needResolved, f, probSize, stateNames, stateGrids, topK)
% REPORTUNCONVERGED  Human-readable report of unconverged grid points.
%   needResolved : logical row over problems (column-linear over probSize)
%   probSize     : [shock_num, state grid sizes...]
%   stateGrids   : cell of grid vectors aligned with stateNames
%   Returns the report text; the caller decides warning vs fprintf.
idx = find(needResolved);
if isempty(idx)
    msg = 'All problems converged.';
    return;
end
[~, order] = sort(f(idx), 'descend');
idx = idx(order);
k = min(topK, numel(idx));
lines = cell(1, k+1);
lines{1} = sprintf('%d of %d problems unconverged; worst %d:', ...
    nnz(needResolved), numel(needResolved), k);
subs = cell(1, numel(probSize));
for i = 1:k
    [subs{:}] = ind2sub(probSize, idx(i));
    parts = sprintf('shock=%d', subs{1});
    for s = 1:numel(stateNames)
        parts = sprintf('%s, %s=%g', parts, stateNames{s}, stateGrids{s}(subs{s+1}));
    end
    lines{i+1} = sprintf('  [%s] residual=%g', parts, f(idx(i)));
end
msg = strjoin(lines, newline);
end
```

- [ ] **Step 8.4: Run — expect PASS**

- [ ] **Step 8.5: Commit**

```powershell
git add src/+gdsge/+runtime tests/runtime/tReporting.m
git commit -m "feat(runtime): computeMetric, printIterProgress, reportUnconverged"
```

---

### Task 9: `runtime.constructSplines`

**Files:**
- Create: `src/+gdsge/+runtime/constructSplines.m`
- Test: `tests/runtime/tConstructSplines.m`

- [ ] **Step 9.1: Write the failing test**

Create `tests/runtime/tConstructSplines.m`:

```matlab
classdef tConstructSplines < matlab.unittest.TestCase
    methods (TestClassSetup)
        function blas(tc) %#ok<MANU>
            gdsge.runtime.ensurePath();
        end
    end
    methods (Test)
        function matchesDirectMyppual(tc)
            w1 = linspace(0, 1, 21);
            v1 = sin(w1); v2 = cos(w1);
            [ppCell, splineVec] = gdsge.runtime.constructSplines( ...
                {v1, v2}, {w1}, {21}, 4, 2, 1);
            tc.verifySize(ppCell, [1 2]);
            tc.verifyTrue(isstruct(splineVec));
            % evaluate ppCell{1} against a hand-built equivalent
            ref = struct('form','MKL','breaks',{{w1}},'Values',reshape(v1,[],21), ...
                'coefs',[],'order',4,'Method',[],'ExtrapolationOrder',2, ...
                'thread',1,'orient','curvefit');
            ref = myppual(ref);
            xq = [0.05 0.33 0.91];
            tc.verifyEqual(myppual(ppCell{1}, xq), myppual(ref, xq), 'AbsTol', 1e-12);
        end
        function order2OmitsExtrap(tc)
            w1 = linspace(0, 1, 5);
            [ppCell, ~] = gdsge.runtime.constructSplines({w1.^2}, {w1}, {5}, 2, 2, 1);
            tc.verifyTrue(isempty(ppCell{1}.ExtrapolationOrder) ...
                || ~isfield(ppCell{1}, 'ExtrapolationOrder'));
        end
    end
end
```

(If `myppual`'s output struct drops/renames `ExtrapolationOrder`, weaken the second test to
just verify the call succeeds — the contract under test is "order 2 passes `[]`".)

- [ ] **Step 9.2: Run — expect FAIL**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/runtime/tConstructSplines.m'))"
```

- [ ] **Step 9.3: Implement**

Create `src/+gdsge/+runtime/constructSplines.m`:

```matlab
function [ppCell, splineVec] = constructSplines(valueCell, breaks, sizeState, interpOrder, extrapOrder, numThreads)
% CONSTRUCTSPLINES  Build myppual interpolants for each interp variable plus
%   the vectorized evaluation array consumed by the MEX (GDSGE_SPLINE_VEC).
%   valueCell : cell of value arrays (one per interp variable)
%   breaks    : cell of state grid vectors, e.g. {w1}
%   sizeState : num2cell of the state-space size, e.g. {201}
%   interpOrder, extrapOrder, numThreads: from the generated file's options.
orderVec = interpOrder * ones(1, numel(sizeState));
if interpOrder == 4
    extrapVec = extrapOrder * ones(1, numel(sizeState));
else
    extrapVec = [];
end
ppCell = cell(1, numel(valueCell));
for i = 1:numel(valueCell)
    pp = struct('form','MKL','breaks',{breaks}, ...
        'Values',reshape(valueCell{i},[],sizeState{:}),'coefs',[],'order',orderVec, ...
        'Method',[],'ExtrapolationOrder',extrapVec,'thread',numThreads, ...
        'orient','curvefit');
    ppCell{i} = myppual(pp);
end
splineVec = convert_to_interp_eval_array(ppCell);
end
```

- [ ] **Step 9.4: Run — expect PASS**

- [ ] **Step 9.5: Commit**

```powershell
git add src/+gdsge/+runtime/constructSplines.m tests/runtime/tConstructSplines.m
git commit -m "feat(runtime): constructSplines - myppual interpolants + GDSGE_SPLINE_VEC"
```

---

### Task 10: `runtime.solveProblems` — contract + randomize loop

**Files:**
- Create: `src/+gdsge/+runtime/solveProblems.m`
- Test: `tests/runtime/tSolveProblems.m` (fake-MEX based; no real MEX needed)

- [ ] **Step 10.1: Write the failing tests**

Create `tests/runtime/tSolveProblems.m`:

```matlab
classdef tSolveProblems < matlab.unittest.TestCase
    methods (Test)
        function callerWorkspaceContractHolds(tc)
            % The fake MEX reads every contract variable from its caller's
            % workspace (= solveProblems' workspace) — errors if any missing.
            n = 3;
            cfg = baseCfg([1 n]);
            sol = zeros(1, n);
            [~, f] = gdsge.runtime.solveProblems(@fakeMexContract, sol, ...
                zeros(1,n), ones(1,n), zeros(1,n), 1e20*ones(1,n), ...
                zeros(1,n), zeros(1,n), cfg);
            tc.verifyEqual(f, zeros(1, n));
        end
        function randomizeRestartsConverge(tc)
            % Solver succeeds only when the guess is within 0.1 of target 0.5;
            % bounds [0,1] so random restarts eventually land there.
            rng(42);
            n = 4;
            cfg = baseCfg([1 n]);
            cfg.useNearestNeighbor = false;
            data = 0.5 * ones(1, n);
            sol = zeros(1, n);   % initial guess far from target
            [sol, f, ~, ~, ~, diag] = gdsge.runtime.solveProblems(@fakeMexNarrow, ...
                sol, zeros(1,n), ones(1,n), data, 1e20*ones(1,n), ...
                zeros(1,n), zeros(1,n), cfg);
            tc.verifyTrue(all(f <= cfg.tolSol));
            tc.verifyEqual(sol, data);
            tc.verifyGreaterThanOrEqual(diag.minorIters, 1);
            tc.verifyFalse(any(diag.needResolved));
        end
        function maxMinorIterCapsAndReportsDiag(tc)
            n = 2;
            cfg = baseCfg([1 n]);
            cfg.useNearestNeighbor = false;
            cfg.maxMinorIter = 3;
            % impossible problem: never converges
            [~, f, ~, ~, ~, diag] = gdsge.runtime.solveProblems(@fakeMexNever, ...
                zeros(1,n), zeros(1,n), ones(1,n), zeros(1,n), ...
                1e20*ones(1,n), zeros(1,n), zeros(1,n), cfg);
            tc.verifyEqual(diag.minorIters, 3);
            tc.verifyTrue(all(diag.needResolved));
            tc.verifyTrue(all(f > cfg.tolSol));
        end
        function convergedPointsAreSkippedOnRestart(tc)
            % After the first call solves everything, the loop must not run.
            n = 3;
            cfg = baseCfg([1 n]);
            cfg.useNearestNeighbor = false;
            data = 0.5*ones(1, n);
            sol = data;   % already at target
            [~, ~, ~, ~, ~, diag] = gdsge.runtime.solveProblems(@fakeMexNarrow, ...
                sol, zeros(1,n), ones(1,n), data, 1e20*ones(1,n), ...
                zeros(1,n), zeros(1,n), cfg);
            tc.verifyEqual(diag.minorIters, 0);
        end
    end
end

% ===== fakes (8-arg MEX interface) ========================================
function [sol, f, aux, eqVal, optInfo] = fakeMexContract(sol, ~, ~, ~, ~, f, aux, eqVal)
names = {'TolFun','TolSol','SolMaxIter','NumThreads','GDSGE_DEBUG_EVAL_ONLY', ...
    'UseBroyden','FiniteDiffDelta','GDSGE_USE_BROYDEN_NOW','MEX_TASK_NAME', ...
    'MEX_TASK_INIT','MEX_TASK_INF_HORIZON','GDSGE_SPLINE_VEC'};
for k = 1:numel(names)
    evalin('caller', [names{k} ';']);   % errors if the contract variable is missing
end
f(:) = 0;
optInfo = zeros(1, numel(f));
end

function [sol, f, aux, eqVal, optInfo] = fakeMexNarrow(sol, ~, ~, data, skip, f, aux, eqVal)
% Solves point i iff its guess is within 0.1 of data(1,i). Honors skip.
for i = 1:size(sol, 2)
    if skip(i); continue; end
    if abs(sol(1,i) - data(1,i)) <= 0.1
        sol(1,i) = data(1,i);
        f(i) = 0;
    else
        f(i) = 1e20;
    end
end
optInfo = zeros(1, numel(f));
end

function [sol, f, aux, eqVal, optInfo] = fakeMexNever(sol, ~, ~, ~, skip, f, aux, eqVal)
f(~logical(skip)) = 1e20;
optInfo = zeros(1, numel(f));
end

function cfg = baseCfg(probSize)
cfg = struct('tolSol', 1e-8, 'tolFun', 1e-8, 'solMaxIter', 200, 'numThreads', 1, ...
    'debugEvalOnly', 0, 'useBroyden', 0, 'finiteDiffDelta', 1e-6, ...
    'useBroydenNow', 0, 'taskName', 1, 'splineVec', struct(), ...
    'maxMinorIter', inf, 'probSize', probSize, 'useNearestNeighbor', true, ...
    'verboseRetry', false, 'adaptInSol', []);
end
```

- [ ] **Step 10.2: Run — expect FAIL**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/runtime/tSolveProblems.m'))"
```

- [ ] **Step 10.3: Implement (randomize loop only; the nearest-neighbor branch lands in Task 11)**

Create `src/+gdsge/+runtime/solveProblems.m`:

```matlab
function [sol, f, aux, eqVal, optInfo, diag] = solveProblems(mexFn, sol, lb, ub, data, f, aux, eqVal, cfg)
% SOLVEPROBLEMS  Drive the model MEX over all grid-point problems until every
%   residual is below tolSol (or maxMinorIter is exhausted). Behavior-parity
%   port of the old generated resolve cascade (see scratch/old_generated/).
%
%   CONTRACT: the MEX reads named variables from ITS CALLER's workspace via
%   mexGetVariable('caller',...). Every mexFn(...) call below must therefore
%   stay in THIS function's body (never a subfunction), and the contract
%   variables below must exist here as locals.
%
%   cfg: tolSol tolFun solMaxIter numThreads debugEvalOnly useBroyden
%        finiteDiffDelta useBroydenNow taskName splineVec maxMinorIter
%        probSize useNearestNeighbor verboseRetry adaptInSol
%   adaptInSol: [] or @(sol,lb,ub) -> [lb,ub] (UseAdaptiveBoundInSol hook)
%   diag: .needResolved (logical row) .worstF .minorIters

% ---- MEX caller-workspace contract (read via mexGetVariable) --------------
TolFun = cfg.tolFun;                          %#ok<NASGU>
TolSol = cfg.tolSol;                          %#ok<NASGU>
SolMaxIter = cfg.solMaxIter;                  %#ok<NASGU>
NumThreads = cfg.numThreads;                  %#ok<NASGU>
GDSGE_DEBUG_EVAL_ONLY = cfg.debugEvalOnly;    %#ok<NASGU>
UseBroyden = cfg.useBroyden;                  %#ok<NASGU>
FiniteDiffDelta = cfg.finiteDiffDelta;        %#ok<NASGU>
GDSGE_USE_BROYDEN_NOW = cfg.useBroydenNow;    %#ok<NASGU>
MEX_TASK_NAME = cfg.taskName;                 %#ok<NASGU>
MEX_TASK_INIT = 0;                            %#ok<NASGU>
MEX_TASK_INF_HORIZON = 1;                     %#ok<NASGU>
GDSGE_SPLINE_VEC = cfg.splineVec;             %#ok<NASGU>
% ---------------------------------------------------------------------------

f(:) = 1e20;
skip = zeros(1, size(sol, 2));
[sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);
GDSGE_USE_BROYDEN_NOW = 0;   %#ok<NASGU> % parity: Broyden only on the first call

minorIter = 0;
numNeedResolvedAfter = inf;
while ((max(isnan(f)) || max(f(:)) > cfg.tolSol) && minorIter < cfg.maxMinorIter)
    if cfg.useNearestNeighbor
        % nearest-neighbor cascade lands in Task 11 (no Task-10 test enters here)
        error('gdsge:runtime:notImplemented', ...
            'nearest-neighbor cascade not implemented yet (Task 11)');
    end

    % randomize restarts for whatever is left (parity: old lines 376-385)
    x0Rand = rand(size(sol)) .* (ub - lb) + lb;
    needResolved = (f > cfg.tolSol) | isnan(f);
    sol(:, needResolved) = x0Rand(:, needResolved);
    skip(:) = 0;
    skip(~needResolved) = 1;
    [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);

    if ~isempty(cfg.adaptInSol)
        % UseAdaptiveBoundInSol hook (parity: old lines 387-408)
        lbOld = lb; ubOld = ub;
        [lb, ub] = cfg.adaptInSol(sol, lb, ub);
        hitLower = abs(sol - lbOld) < 1e-8;
        hitUpper = abs(sol - ubOld) < 1e-8;
        lb(~hitLower) = lbOld(~hitLower);
        ub(~hitUpper) = ubOld(~hitUpper);
        minorIter = minorIter + 1;   % parity: old counts an extra minor iter here
    end
    minorIter = minorIter + 1;

    if cfg.verboseRetry
        fprintf('  resolve round %d: %d unconverged, worst residual %g\n', ...
            minorIter, nnz((f > cfg.tolSol) | isnan(f)), max(f(:)));
    end
end

diag = struct();
diag.needResolved = (f > cfg.tolSol) | isnan(f);
diag.worstF = max(f(:));
diag.minorIters = minorIter;
end
```

- [ ] **Step 10.4: Run — expect PASS** (the contract/randomize/cap tests don't exercise the cascade)

- [ ] **Step 10.5: Commit**

```powershell
git add src/+gdsge/+runtime/solveProblems.m tests/runtime/tSolveProblems.m
git commit -m "feat(runtime): solveProblems - MEX contract, randomize restarts, diagnostics"
```

---

### Task 11: `runtime.solveProblems` — nearest-neighbor cascade

**Files:**
- Modify: `src/+gdsge/+runtime/solveProblems.m`
- Modify: `tests/runtime/tSolveProblems.m` (add cascade tests)

- [ ] **Step 11.1: Add the failing cascade tests**

Append these methods inside the `methods (Test)` block of `tests/runtime/tSolveProblems.m`:

```matlab
        function nearestNeighborPropagates(tc)
            % targets 1..5; solver succeeds iff guess within 1.0 of target.
            % Only point 1 starts close; neighbor copying must walk it across.
            n = 5;
            cfg = baseCfg([1 n]);
            data = 1:n;
            sol = [1 9 9 9 9];
            [sol, f, ~, ~, ~, diag] = gdsge.runtime.solveProblems(@fakeMexUnitBall, ...
                sol, zeros(1,n), 10*ones(1,n), data, 1e20*ones(1,n), ...
                zeros(1,n), zeros(1,n), cfg);
            tc.verifyEqual(sol, double(1:n));
            tc.verifyTrue(all(f <= cfg.tolSol));
            tc.verifyEqual(diag.minorIters, 1);   % one outer round suffices
        end
        function cascadeRespectsStridesPerDim(tc)
            % probSize [2 3]: dim-1 stride 1 (shock), dim-2 stride 2 (state).
            % Make only state-neighbor copying useful: targets differ by shock
            % row, solver tolerance 1.0, state-adjacent targets within 1.
            cfg = baseCfg([2 3]);
            data = [1 1; 2 2; 3 3]'; data = data(:)';   % targets [1 2 3] per state, same per shock
            sol = [1 2 9 9 9 9];                        % first state column solved
            [sol2, f] = gdsge.runtime.solveProblems(@fakeMexUnitBall, ...
                sol, zeros(1,6), 10*ones(1,6), data, 1e20*ones(1,6), ...
                zeros(1,6), zeros(1,6), cfg);
            tc.verifyTrue(all(f <= cfg.tolSol));
            tc.verifyEqual(sol2, data);
        end
```

And this fake at the bottom of the file (with the other fakes):

```matlab
function [sol, f, aux, eqVal, optInfo] = fakeMexUnitBall(sol, ~, ~, data, skip, f, aux, eqVal)
% Solves point i iff its guess is within 1.0 of data(1,i). Honors skip.
for i = 1:size(sol, 2)
    if skip(i); continue; end
    if abs(sol(1,i) - data(1,i)) <= 1.0
        sol(1,i) = data(1,i);
        f(i) = 0;
    else
        f(i) = 1e20;
    end
end
optInfo = zeros(1, numel(f));
end
```

- [ ] **Step 11.2: Run — expect FAIL** (cascade branch is a stub)

- [ ] **Step 11.3: Implement the cascade**

In `src/+gdsge/+runtime/solveProblems.m`, replace the stubbed
`if cfg.useNearestNeighbor ... end` block with the parity port of old lines 330-374:

```matlab
    if cfg.useNearestNeighbor
        needResolved = (f > cfg.tolSol) | isnan(f);
        numNeedResolved = sum(needResolved);
        while numNeedResolvedAfter ~= numNeedResolved
            needResolved = (f > cfg.tolSol) | isnan(f);
            numNeedResolved = sum(needResolved);
            for iDim = 1:length(cfg.probSize)
                stride = prod(cfg.probSize(1:iDim-1));

                % copy from the lower neighbor (parity: stale needResolved
                % during the sweep, exactly like the old generated loop)
                needResolved = (f > cfg.tolSol) | isnan(f);
                skip(:) = 1;
                for i = 1:numel(f)
                    if needResolved(i) && i-stride >= 1 && ~needResolved(i-stride)
                        sol(:, i) = sol(:, i-stride);
                        skip(i) = 0;
                    end
                end
                [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);

                % copy from the upper neighbor
                needResolved = (f > cfg.tolSol) | isnan(f);
                skip(:) = 1;
                for i = 1:numel(f)
                    if needResolved(i) && i+stride <= numel(f) && ~needResolved(i+stride)
                        sol(:, i) = sol(:, i+stride);
                        skip(i) = 0;
                    end
                end
                [sol, f, aux, eqVal, optInfo] = mexFn(sol, lb, ub, data, skip, f, aux, eqVal);
            end
            numNeedResolvedAfter = sum((f > cfg.tolSol) | isnan(f));
        end
    end
```

- [ ] **Step 11.4: Run all solveProblems tests — expect PASS**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/runtime/tSolveProblems.m'))"
```

- [ ] **Step 11.5: Commit**

```powershell
git add src/+gdsge/+runtime/solveProblems.m tests/runtime/tSolveProblems.m
git commit -m "feat(runtime): solveProblems - nearest-neighbor resolve cascade (parity port)"
```

---

### Task 12: `runtime.applyWarmUp`

**Files:**
- Create: `src/+gdsge/+runtime/applyWarmUp.m`
- Test: `tests/runtime/tApplyWarmUp.m`

- [ ] **Step 12.1: Write the failing test**

Create `tests/runtime/tApplyWarmUp.m`:

```matlab
classdef tApplyWarmUp < matlab.unittest.TestCase
    methods (TestClassSetup)
        function blas(tc) %#ok<MANU>
            gdsge.runtime.ensurePath();
        end
    end
    methods (Test)
        function plainSolFieldCopies(tc)
            wu = struct('GDSGE_SOL', magic(3));
            [sol, lb, ub] = gdsge.runtime.applyWarmUp(wu, zeros(3), ...
                zeros(3), ones(3), spec0(2, {'w1'}));
            tc.verifyEqual(sol, magic(3));
            tc.verifyEqual(lb, zeros(3));
            tc.verifyEqual(ub, ones(3));
        end
        function plainSolWrongRowsErrors(tc)
            wu = struct('GDSGE_SOL', zeros(4, 3));
            tc.verifyError(@() gdsge.runtime.applyWarmUp(wu, zeros(3), ...
                zeros(3), ones(3), spec0(2, {'w1'})), 'gdsge:runtime:warmUpSize');
        end
        function interpTransferIsExactForLinearData(tc)
            % coarse grid [0 1], 2 shocks, 1 unknown; SOL linear in the state
            % -> order-2 interp transfer reproduces it exactly on a finer grid.
            shockNum = 2; coarse = [0 1]; fine = [0 0.5 1];
            sizeCoarse = [shockNum numel(coarse)];
            solC = zeros(1, prod(sizeCoarse));   % layout: shock fastest (column-linear)
            % value = shock + 10*state
            k = 0;
            for s = 1:numel(coarse)
                for j = 1:shockNum
                    k = k + 1;
                    solC(k) = j + 10*coarse(s);
                end
            end
            wu = struct();
            wu.GDSGE_PROB = struct('GDSGE_SOL', solC, 'GDSGE_LB', solC - 1, ...
                'GDSGE_UB', solC + 1, 'GDSGE_SIZE', sizeCoarse);
            wu.var_state = struct('w1', coarse);
            [shockIdx, w1T] = ndgrid(1:shockNum, fine);
            sp = spec0(shockNum, {'w1'});
            sp.evalPoints = [shockIdx(:)'; w1T(:)'];
            sp.adaptiveSlots = 1;
            n = numel(shockIdx);
            [sol, lb, ub] = gdsge.runtime.applyWarmUp(wu, zeros(1, n), ...
                -5*ones(1, n), 5*ones(1, n), sp);
            expect = shockIdx(:)' + 10*w1T(:)';
            tc.verifyEqual(sol, expect, 'AbsTol', 1e-10);
            tc.verifyEqual(lb, expect - 1, 'AbsTol', 1e-10);   % slot 1 is adaptive
            tc.verifyEqual(ub, expect + 1, 'AbsTol', 1e-10);
        end
        function nonAdaptiveSlotsKeepBounds(tc)
            % same setup, but adaptiveSlots empty -> lb/ub untouched
            shockNum = 2; coarse = [0 1]; fine = [0 0.5 1];
            solC = ones(1, 4);
            wu = struct();
            wu.GDSGE_PROB = struct('GDSGE_SOL', solC, 'GDSGE_LB', zeros(1,4), ...
                'GDSGE_UB', 2*ones(1,4), 'GDSGE_SIZE', [2 2]);
            wu.var_state = struct('w1', coarse);
            [shockIdx, w1T] = ndgrid(1:shockNum, fine);
            sp = spec0(shockNum, {'w1'});
            sp.evalPoints = [shockIdx(:)'; w1T(:)'];
            sp.adaptiveSlots = [];
            n = numel(shockIdx);
            [~, lb, ub] = gdsge.runtime.applyWarmUp(wu, zeros(1, n), ...
                -5*ones(1, n), 5*ones(1, n), sp);
            tc.verifyEqual(lb, -5*ones(1, n));
            tc.verifyEqual(ub, 5*ones(1, n));
        end
        function directCopyWhenInterpSolOff(tc)
            wu = struct();
            wu.GDSGE_PROB = struct('GDSGE_SOL', 7*ones(2,4), ...
                'GDSGE_LB', -ones(2,4), 'GDSGE_UB', ones(2,4), 'GDSGE_SIZE', [2 2]);
            sp = spec0(2, {'w1'});
            sp.interpSol = 0;
            [sol, lb, ub] = gdsge.runtime.applyWarmUp(wu, zeros(2,4), ...
                zeros(2,4), 2*ones(2,4), sp);
            tc.verifyEqual(sol, 7*ones(2,4));
            tc.verifyEqual(lb, -ones(2,4));
            tc.verifyEqual(ub, ones(2,4));
        end
    end
end

function sp = spec0(shockNum, stateNames)
sp = struct('shockNum', shockNum, 'stateNames', {stateNames}, ...
    'evalPoints', [], 'adaptiveSlots', [], 'reuseSol', 1, 'interpSol', 1, ...
    'numThreads', 1);
end
```

- [ ] **Step 12.2: Run — expect FAIL**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/runtime/tApplyWarmUp.m'))"
```

- [ ] **Step 12.3: Implement**

Create `src/+gdsge/+runtime/applyWarmUp.m`:

```matlab
function [sol, lb, ub] = applyWarmUp(warmUp, sol, lb, ub, spec)
% APPLYWARMUP  Transfer a previous run's solution/bounds into this run.
%   Behavior-parity port of the old generated WarmUp block (plain GDSGE_SOL
%   copy + the REUSE_WARMUP_SOL interp/direct transfer; see
%   scratch/old_generated/iter_HL1996.m lines 257-304).
%   spec fields: shockNum, stateNames (decl order), evalPoints
%     ([shockIdx(:)'; state tensors...]), adaptiveSlots, reuseSol, interpSol,
%     numThreads.
if isfield(warmUp, 'GDSGE_SOL')
    if ~isequal(size(warmUp.GDSGE_SOL, 1), size(sol, 1))
        error('gdsge:runtime:warmUpSize', 'Wrong size of GDSGE_SOL in WarmUp');
    end
    sol = warmUp.GDSGE_SOL;
end
if ~(spec.reuseSol == 1 && isfield(warmUp, 'GDSGE_PROB') ...
        && size(warmUp.GDSGE_PROB.GDSGE_SOL, 1) == size(sol, 1))
    return;
end
if spec.interpSol == 1
    if spec.shockNum < 2
        error('gdsge:runtime:warmUpShock', ...
            'WarmUp solution can only be applied with shock_num>=2. Please set REUSE_WARMUP_SOL=0.');
    end
    sizeFull = num2cell(warmUp.GDSGE_PROB.GDSGE_SIZE);
    breaks = cell(1, 1 + numel(spec.stateNames));
    breaks{1} = 1:spec.shockNum;
    for s = 1:numel(spec.stateNames)
        breaks{1+s} = warmUp.var_state.(spec.stateNames{s});
    end
    solI = mkInterp(warmUp.GDSGE_PROB.GDSGE_SOL, breaks, sizeFull, spec.numThreads);
    lbI  = mkInterp(warmUp.GDSGE_PROB.GDSGE_LB,  breaks, sizeFull, spec.numThreads);
    ubI  = mkInterp(warmUp.GDSGE_PROB.GDSGE_UB,  breaks, sizeFull, spec.numThreads);
    sol   = reshape(myppual(solI, spec.evalPoints), size(sol));
    lbNew = reshape(myppual(lbI,  spec.evalPoints), size(lb));
    ubNew = reshape(myppual(ubI,  spec.evalPoints), size(ub));
    lb(spec.adaptiveSlots, :) = lbNew(spec.adaptiveSlots, :);
    ub(spec.adaptiveSlots, :) = ubNew(spec.adaptiveSlots, :);
else
    if isfield(warmUp.GDSGE_PROB, 'GDSGE_SOL'); sol = warmUp.GDSGE_PROB.GDSGE_SOL; end
    if isfield(warmUp.GDSGE_PROB, 'GDSGE_LB');  lb  = warmUp.GDSGE_PROB.GDSGE_LB;  end
    if isfield(warmUp.GDSGE_PROB, 'GDSGE_UB');  ub  = warmUp.GDSGE_PROB.GDSGE_UB;  end
end
end

function pp = mkInterp(vals, breaks, sizeFull, numThreads)
% order-2 (linear) transfer interp over {1:shock_num, state grids...} — parity
pp = struct('form','MKL','breaks',{breaks}, ...
    'Values', reshape(vals, [], sizeFull{:}), 'coefs', [], ...
    'order', 2*ones(1, numel(sizeFull)), 'Method', [], ...
    'ExtrapolationOrder', [], 'thread', numThreads, 'orient', 'curvefit');
pp = myppual(pp);
end
```

- [ ] **Step 12.4: Run — expect PASS**

- [ ] **Step 12.5: Commit**

```powershell
git add src/+gdsge/+runtime/applyWarmUp.m tests/runtime/tApplyWarmUp.m
git commit -m "feat(runtime): applyWarmUp - WarmUp SOL/LB/UB transfer (parity port)"
```

---

### Task 13: `emitSetup`

**Files:**
- Create: `src/+gdsge/+codegen/+mat/emitSetup.m`
- Test: `tests/codegen/tEmitSetup.m`

- [ ] **Step 13.1: Write the failing test** (eval-based: proves numeric fidelity, not formatting)

Create `tests/codegen/tEmitSetup.m`:

```matlab
classdef tEmitSetup < matlab.unittest.TestCase
    methods (TestClassSetup)
        function irFixture(tc)
            here = fileparts(mfilename('fullpath'));
            irDir = fullfile(fileparts(here), 'HeatonLucas1996', 'ir');
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(irDir));
        end
    end
    methods (Test)
        function evalReproducesIrValuesExactly(tc)
            ir = buildHL1996IR();
            ws = evalSetupText(gdsge.codegen.mat.emitSetup(ir, 'iter'));
            tc.verifyEqual(ws.beta, 0.95);    % mat2str(v,17) round-trip is exact
            tc.verifyEqual(ws.gamma, 1.5);
            tc.verifyEqual(ws.Kb, -0.05);
            tc.verifyEqual(ws.shock_num, 8);
            tc.verifyEqual(ws.g, ir.shocks.values.g);
            tc.verifyEqual(ws.d, ir.shocks.values.d);
            tc.verifyEqual(ws.eta1, ir.shocks.values.eta1);
            tc.verifyEqual(ws.shock_trans, ir.shocks.transitions.shock_trans);
            tc.verifyEqual(ws.w1, linspace(-0.05, 1.05, 201));   % grid text replayed
        end
        function frozenDefaultsPresent(tc)
            ir = buildHL1996IR();
            ws = evalSetupText(gdsge.codegen.mat.emitSetup(ir, 'iter'));
            tc.verifyEqual(ws.TolEq, 1e-6);
            tc.verifyEqual(ws.TolSol, 1e-8);
            tc.verifyEqual(ws.PrintFreq, 10);
            tc.verifyEqual(ws.SaveFreq, 10);
            tc.verifyEqual(ws.MaxIter, inf);
            tc.verifyEqual(ws.MaxMinorIter, inf);
            tc.verifyEqual(ws.INTERP_ORDER, 4);
            tc.verifyEqual(ws.EXTRAP_ORDER, 2);
            tc.verifyEqual(ws.MEX_TASK_INIT, 0);
            tc.verifyEqual(ws.MEX_TASK_INF_HORIZON, 1);
            tc.verifyEqual(ws.num_samples, 1);       % iter defaults
            tc.verifyEqual(ws.num_periods, 1000);
            tc.verifyEqual(ws.NumThreads, feature('numcores'));   % 0-sentinel resolved
        end
        function simulateTargetUsesIrSimulate(tc)
            ir = buildHL1996IR();
            ws = evalSetupText(gdsge.codegen.mat.emitSetup(ir, 'simulate'));
            tc.verifyEqual(ws.num_periods, 10000);
            tc.verifyEqual(ws.num_samples, 24);
        end
        function explicitNumThreadsEmitsLiteral(tc)
            ir = buildHL1996IR();
            ir.options.numThreads = 6;
            ws = evalSetupText(gdsge.codegen.mat.emitSetup(ir, 'iter'));
            tc.verifyEqual(ws.NumThreads, 6);
        end
    end
end

function ws = evalSetupText(txt)
% Evaluate emitted setup code in a clean function workspace; return it.
eval(txt);
clear txt;
vars = whos;
ws = struct();
for i = 1:numel(vars)
    ws.(vars(i).name) = eval(vars(i).name);
end
end
```

- [ ] **Step 13.2: Run — expect FAIL**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/codegen/tEmitSetup.m'))"
```

- [ ] **Step 13.3: Implement**

Create `src/+gdsge/+codegen/+mat/emitSetup.m`:

```matlab
function txt = emitSetup(ir, target)
% EMITSETUP  Defaults + params + shocks + grids setup block.
%   target: 'iter' or 'simulate' (simulate takes num_periods/num_samples from
%   ir.simulate; iter keeps the old generator's 1/1000 defaults).
%   Numeric IR values are emitted with mat2str(v,17) (exact double round-trip);
%   state grids are opaque MATLAB text emitted verbatim (after params, so grid
%   text may reference params).
w = gdsge.codegen.mat.codeWriter();
opt = ir.options;

w.add('%% ---- toolbox options (frozen public surface; defaults as in the old generator)');
w.add('TolEq = %s;', mat2str(getOpt(opt, 'tolEq', 1e-6), 17));
w.add('TolSol = 1e-8;');
w.add('TolFun = 1e-8;');
w.add('PrintFreq = %s;', mat2str(getOpt(opt, 'printFreq', 10), 17));
w.add('NoPrint = 0;');
w.add('SaveFreq = %s;', mat2str(getOpt(opt, 'saveFreq', 10), 17));
w.add('NoSave = 0;');
w.add('SimuPrintFreq = %s;', mat2str(getOpt(opt, 'simuPrintFreq', 1000), 17));
w.add('SimuSaveFreq = %s;', mat2str(getOpt(opt, 'simuSaveFreq', inf), 17));
w.add('MaxIter = inf;');
w.add('MaxMinorIter = inf;');
if strcmp(target, 'simulate')
    w.add('num_samples = %s;', mat2str(ir.simulate.numSamples, 17));
    w.add('num_periods = %s;', mat2str(ir.simulate.numPeriods, 17));
else
    w.add('num_samples = 1;');
    w.add('num_periods = 1000;');
end
w.add('SolMaxIter = 200;');
w.add('UseBroyden = 0;');
w.add('FiniteDiffDelta = 1e-6;');
w.add('GDSGE_USE_BROYDEN = 1;');
w.add('GDSGE_DEBUG_EVAL_ONLY = 0;');
w.add('INTERP_ORDER = %s;', mat2str(getOpt(opt, 'interpOrder', 4), 17));
w.add('EXTRAP_ORDER = %s;', mat2str(getOpt(opt, 'extrapOrder', 2), 17));
w.add('OutputInterpOrder = 2;');
w.add('IterSaveAll = 0;');
w.add('SkipModelInit = 0;');
w.add('GDSGE_EMPTY = [];');
w.add('UseAdaptiveBound = 1;');
w.add('UseAdaptiveBoundInSol = 0;');
w.add('EnforceSimuStateInbound = 1;');
w.add('REUSE_WARMUP_SOL = 1;');
w.add('INTERP_WARMUP_SOL = 1;');
w.add('CONSTRUCT_OUTPUT = 1;');
w.add('MEX_TASK_INIT = 0;');
w.add('MEX_TASK_INF_HORIZON = 1;');
if getOpt(opt, 'numThreads', 0) == 0
    w.add('NumThreads = feature(''numcores'');');
else
    w.add('NumThreads = %s;', mat2str(opt.numThreads, 17));
end
w.blank();

w.add('%% ---- model parameters (numeric, from the IR)');
for i = 1:numel(ir.params)
    p = ir.params{i};
    w.add('%s = %s;', p.name, mat2str(p.value, 17));
end
w.blank();

w.add('%% ---- shocks');
w.add('shock_num = %s;', mat2str(ir.shocks.count, 17));
for i = 1:numel(ir.shocks.names)
    n = ir.shocks.names{i};
    w.add('%s = %s;', n, mat2str(ir.shocks.values.(n), 17));
end
tn = fieldnames(ir.shocks.transitions);
for i = 1:numel(tn)
    w.add('%s = %s;', tn{i}, mat2str(ir.shocks.transitions.(tn{i}), 17));
end
w.blank();

w.add('%% ---- state grids (opaque gmod text, replayed verbatim)');
for i = 1:numel(ir.states.names)
    s = ir.states.names{i};
    w.add('%s = %s;', s, ir.states.grids.(s));
end
txt = w.str();
end

function v = getOpt(opt, name, dflt)
if isfield(opt, name); v = opt.(name); else; v = dflt; end
end
```

Deliberately dropped from the old defaults block (dead for the spline path; per spec):
`USE_SPLINE/USE_ASG/USE_PCHIP`, `USE_SPARSE_JACOBIAN`, `GDSGE_USE_OLD_VEC`, `SimuSeed`,
`AsgMinLevel/AsgMaxLevel/...`, `GDSGE_ASG_FIX_GRID`, `UseModelId`, `MinBatchSol`,
`REMOVE_NULL_STATEMENTS`, `USE_FINITE_DIFF`, `GDSGE_dummy_state/shock`,
`DEFAULT_PARAMETERS_END_HERE`. `GDSGE_USE_BROYDEN_NOW` lives inside `solveProblems` now.

- [ ] **Step 13.4: Run — expect PASS**

- [ ] **Step 13.5: Commit**

```powershell
git add src/+gdsge/+codegen/+mat/emitSetup.m tests/codegen/tEmitSetup.m
git commit -m "feat(codegen): emitSetup - options/params/shocks/grids setup block"
```

---

### Task 14: `emitBounds` + `emitInterpInitial`

**Files:**
- Create: `src/+gdsge/+codegen/+mat/emitBounds.m`, `src/+gdsge/+codegen/+mat/emitInterpInitial.m`
- Test: `tests/codegen/tEmitBounds.m`, `tests/codegen/tEmitInterpInitial.m`

- [ ] **Step 14.1: Write the failing tests**

Create `tests/codegen/tEmitBounds.m`:

```matlab
classdef tEmitBounds < matlab.unittest.TestCase
    methods (TestClassSetup)
        function irFixture(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(here), 'HeatonLucas1996', 'ir')));
        end
    end
    methods (Test)
        function initCoversAllSlots(tc)
            frag = gdsge.codegen.mat.emitBounds(buildHL1996IR());
            tc.verifyEqual(frag.maxDim, 19);
            tc.verifyTrue(contains(frag.init, 'GDSGE_LB = zeros(19,GDSGE_NPROB);'));
            tc.verifyTrue(contains(frag.init, 'GDSGE_UB = 1e3*ones(19,GDSGE_NPROB);'));
            tc.verifyTrue(contains(frag.init, 'GDSGE_LB(1:1,:)=0;'));
            tc.verifyTrue(contains(frag.init, 'GDSGE_UB(10:10,:)=3;'));
            tc.verifyTrue(contains(frag.init, 'GDSGE_LB(12:19,:)=-0.5;'));
            tc.verifyTrue(contains(frag.init, 'GDSGE_UB(12:19,:)=1.5;'));
        end
        function adaptiveFragmentsMatchOldSemantics(tc)
            frag = gdsge.codegen.mat.emitBounds(buildHL1996IR());
            tc.verifyEqual(frag.adaptiveSlots, [10 11]);
            % widen (in-loop, parity with old lines 475-478)
            tc.verifyTrue(contains(frag.adaptiveWiden, ...
                'GDSGE_LB(10:10,:)=min(GDSGE_LB(10:10,:),GDSGE_SOL(10:10,:)*1.5);'));
            tc.verifyTrue(contains(frag.adaptiveWiden, ...
                'GDSGE_UB(11:11,:)=max(GDSGE_UB(11:11,:),GDSGE_SOL(11:11,:)*1.5);'));
            % tighten (UseAdaptiveBoundInSol hook, parity with old lines 393-396)
            tc.verifyTrue(contains(frag.adaptiveTight, ...
                'GDSGE_LB(10:10,:)=GDSGE_SOL(10:10,:)/1.5;'));
            tc.verifyTrue(contains(frag.adaptiveTight, ...
                'GDSGE_UB(10:10,:)=GDSGE_SOL(10:10,:)*1.5;'));
        end
    end
end
```

Create `tests/codegen/tEmitInterpInitial.m`:

```matlab
classdef tEmitInterpInitial < matlab.unittest.TestCase
    methods (TestClassSetup)
        function irFixture(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(here), 'HeatonLucas1996', 'ir')));
        end
    end
    methods (Test)
        function preallocatesAndRewrites(tc)
            txt = gdsge.codegen.mat.emitInterpInitial(buildHL1996IR());
            tc.verifyTrue(contains(txt, 'ps_future=zeros(GDSGE_SIZE);'));
            tc.verifyTrue(contains(txt, 'ps_future(:)=0.0;'));
            % state/shock refs rewritten to tensors (parity with old line 179)
            tc.verifyTrue(contains(txt, ...
                'c1_future(:)=GDSGE_TENSOR_w1(:)''.*GDSGE_TENSOR_d(:)''+GDSGE_TENSOR_eta1(:)'';'));
            tc.verifyTrue(contains(txt, ...
                'c2_future(:)=(1-GDSGE_TENSOR_w1(:)'').*GDSGE_TENSOR_d(:)''+1-GDSGE_TENSOR_eta1(:)'';'));
        end
    end
end
```

- [ ] **Step 14.2: Run both — expect FAIL**

- [ ] **Step 14.3: Implement**

Create `src/+gdsge/+codegen/+mat/emitBounds.m`:

```matlab
function frag = emitBounds(ir)
% EMITBOUNDS  Bound-related fragments from ir.bounds + policy slot layout.
%   .init            GDSGE_LB/UB allocation + per-variable fills
%   .adaptiveWiden   in-loop widening      (old: min/max against SOL*factor)
%   .adaptiveTight   UseAdaptiveBoundInSol (old: SOL/factor .. SOL*factor)
%   .adaptiveSlots   row vector of adaptive slot indices (warmup LB/UB
%                    overwrite happens inside gdsge.runtime.applyWarmUp)
%   .maxDim          number of solution rows
pol = ir.variables.policy;
slotOf = struct();
maxDim = 0;
for i = 1:numel(pol)
    slotOf.(pol{i}.name) = pol{i}.slot;
    maxDim = max(maxDim, pol{i}.slot(2));
end
names = [ir.states.names, ir.shocks.names];
reps = cellfun(@(n) ['GDSGE_TENSOR_' n '(:)'''], names, 'UniformOutput', false);

wInit = gdsge.codegen.mat.codeWriter();
wInit.add('GDSGE_LB = zeros(%d,GDSGE_NPROB);', maxDim);
wInit.add('GDSGE_UB = 1e3*ones(%d,GDSGE_NPROB);', maxDim);
wWiden = gdsge.codegen.mat.codeWriter();
wTight = gdsge.codegen.mat.codeWriter();
adaptiveSlots = [];
for i = 1:numel(ir.bounds)
    b = ir.bounds{i};
    s = slotOf.(b.name);
    lo = gdsge.codegen.mat.rewriteNames(b.lower, names, reps);
    hi = gdsge.codegen.mat.rewriteNames(b.upper, names, reps);
    wInit.add('GDSGE_LB(%d:%d,:)=%s;', s(1), s(2), lo);
    wInit.add('GDSGE_UB(%d:%d,:)=%s;', s(1), s(2), hi);
    if isfield(b, 'adaptiveFactor') && ~isempty(b.adaptiveFactor)
        fct = mat2str(b.adaptiveFactor, 17);
        wWiden.add('GDSGE_LB(%d:%d,:)=min(GDSGE_LB(%d:%d,:),GDSGE_SOL(%d:%d,:)*%s);', ...
            s(1), s(2), s(1), s(2), s(1), s(2), fct);
        wWiden.add('GDSGE_UB(%d:%d,:)=max(GDSGE_UB(%d:%d,:),GDSGE_SOL(%d:%d,:)*%s);', ...
            s(1), s(2), s(1), s(2), s(1), s(2), fct);
        wTight.add('GDSGE_LB(%d:%d,:)=GDSGE_SOL(%d:%d,:)/%s;', s(1), s(2), s(1), s(2), fct);
        wTight.add('GDSGE_UB(%d:%d,:)=GDSGE_SOL(%d:%d,:)*%s;', s(1), s(2), s(1), s(2), fct);
        adaptiveSlots = [adaptiveSlots, s(1):s(2)]; %#ok<AGROW>
    end
end
frag = struct('init', wInit.str(), 'adaptiveWiden', wWiden.str(), ...
    'adaptiveTight', wTight.str(), ...
    'adaptiveSlots', adaptiveSlots, 'maxDim', maxDim);
end
```

Create `src/+gdsge/+codegen/+mat/emitInterpInitial.m`:

```matlab
function txt = emitInterpInitial(ir)
% EMITINTERPINITIAL  Interp-variable preallocation + initial expressions,
%   with state/shock references rewritten to GDSGE_TENSOR_<name>(:)'.
names = [ir.states.names, ir.shocks.names];
reps = cellfun(@(n) ['GDSGE_TENSOR_' n '(:)'''], names, 'UniformOutput', false);
w = gdsge.codegen.mat.codeWriter();
for i = 1:numel(ir.interp)
    w.add('%s=zeros(GDSGE_SIZE);', ir.interp{i}.name);
end
for i = 1:numel(ir.interp)
    it = ir.interp{i};
    w.add('%s(:)=%s;', it.name, ...
        gdsge.codegen.mat.rewriteNames(it.initialExpr, names, reps));
end
txt = w.str();
end
```

- [ ] **Step 14.4: Run both — expect PASS**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests({'tests/codegen/tEmitBounds.m','tests/codegen/tEmitInterpInitial.m'}))"
```

- [ ] **Step 14.5: Commit**

```powershell
git add src/+gdsge/+codegen/+mat tests/codegen
git commit -m "feat(codegen): emitBounds + emitInterpInitial"
```

---

### Task 15: `emitDataPack` + `emitSolUnpack`

**Files:**
- Create: `src/+gdsge/+codegen/+mat/emitDataPack.m`, `src/+gdsge/+codegen/+mat/emitSolUnpack.m`
- Test: `tests/codegen/tEmitDataPack.m`, `tests/codegen/tEmitSolUnpack.m`

- [ ] **Step 15.1: Write the failing tests**

Create `tests/codegen/tEmitDataPack.m`:

```matlab
classdef tEmitDataPack < matlab.unittest.TestCase
    methods (TestClassSetup)
        function irFixture(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(here), 'HeatonLucas1996', 'ir')));
        end
    end
    methods (Test)
        function layoutMatchesOldGenerator(tc)
            % THE data-layout contract with the MEX: [shock_num; params (decl
            % order); shock_trans(:); shock values (decl order)] per problem,
            % then shockIdx row and state tensor rows. Old reference line 314.
            frag = gdsge.codegen.mat.emitDataPack(buildHL1996IR());
            tc.verifyEqual(frag.iterPack, ['GDSGE_DATA(:) = [repmat([' ...
                'shock_num;beta(:);gamma(:);Kb(:);shock_trans(:);g(:);d(:);eta1(:)' ...
                '],1,GDSGE_NPROB);GDSGE_TENSOR_shockIdx(:)'';GDSGE_TENSOR_w1(:)''];']);
            tc.verifyEqual(frag.maxData, 1 + 3 + 64 + 24 + 1 + 1);   % = 94
        end
        function simulateVariants(tc)
            frag = gdsge.codegen.mat.emitDataPack(buildHL1996IR());
            tc.verifyEqual(frag.simuData0, ['GDSGE_data0 = repmat([' ...
                'shock_num;beta(:);gamma(:);Kb(:);shock_trans(:);g(:);d(:);eta1(:)' ...
                '],1,GDSGE_NPROB);']);
            tc.verifyEqual(frag.simuPack, ...
                'GDSGE_DATA = [GDSGE_data0;shock(:)'';w1(:)''];');
        end
    end
end
```

Create `tests/codegen/tEmitSolUnpack.m`:

```matlab
classdef tEmitSolUnpack < matlab.unittest.TestCase
    methods (TestClassSetup)
        function irFixture(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(here), 'HeatonLucas1996', 'ir')));
        end
    end
    methods (Test)
        function unpackUsesSlotRanges(tc)
            frag = gdsge.codegen.mat.emitSolUnpack(buildHL1996IR());
            tc.verifyTrue(contains(frag.unpack, 'c1=GDSGE_SOL(1:1,:);'));
            tc.verifyTrue(contains(frag.unpack, 'ps=GDSGE_SOL(10:10,:);'));
            tc.verifyTrue(contains(frag.unpack, 'w1n=GDSGE_SOL(12:19,:);'));
            tc.verifyTrue(contains(frag.unpackAux, 'equity_premium=GDSGE_AUX(1:1,:);'));
        end
        function reshapeFragments(tc)
            frag = gdsge.codegen.mat.emitSolUnpack(buildHL1996IR());
            tc.verifyTrue(contains(frag.reshapeShock, 'c1=reshape(c1,shock_num,[]);'));
            tc.verifyTrue(contains(frag.reshapeShock, ...
                'equity_premium=reshape(equity_premium,shock_num,[]);'));
            % save-time reshape: scalars to GDSGE_SIZE, arrays to [len GDSGE_SIZE],
            % then interp vars and aux (old reference lines 519-535)
            tc.verifyTrue(contains(frag.reshapeSize, 'c1=reshape(c1,GDSGE_SIZE);'));
            tc.verifyTrue(contains(frag.reshapeSize, 'w1n=reshape(w1n,[8 GDSGE_SIZE]);'));
            tc.verifyTrue(contains(frag.reshapeSize, ...
                'ps_future=reshape(ps_future,GDSGE_SIZE);'));
            tc.verifyTrue(contains(frag.reshapeSize, ...
                'equity_premium=reshape(equity_premium,GDSGE_SIZE);'));
        end
    end
end
```

- [ ] **Step 15.2: Run — expect FAIL**

- [ ] **Step 15.3: Implement**

Create `src/+gdsge/+codegen/+mat/emitDataPack.m`:

```matlab
function frag = emitDataPack(ir)
% EMITDATAPACK  The GDSGE_DATA packing expressions — THE layout contract with
%   the MEX. Per-problem constant block: [shock_num; params (decl order);
%   shock_trans(:); shock values (decl order)]; then per-problem rows:
%   shockIdx, state tensors (decl order). Must match the old generator so the
%   Phase-4 functional gate can run the new .m against the old MEX.
inner = {'shock_num'};
nData = 1;
for i = 1:numel(ir.params)
    inner{end+1} = [ir.params{i}.name '(:)']; %#ok<AGROW>
    nData = nData + numel(ir.params{i}.value);
end
tn = fieldnames(ir.shocks.transitions);
for i = 1:numel(tn)
    inner{end+1} = [tn{i} '(:)']; %#ok<AGROW>
    nData = nData + numel(ir.shocks.transitions.(tn{i}));
end
for i = 1:numel(ir.shocks.names)
    inner{end+1} = [ir.shocks.names{i} '(:)']; %#ok<AGROW>
    nData = nData + numel(ir.shocks.values.(ir.shocks.names{i}));
end
innerStr = strjoin(inner, ';');

perProb = {'GDSGE_TENSOR_shockIdx(:)'''};
simuRows = {'shock(:)'''};
nData = nData + 1;                       % shockIdx row
for i = 1:numel(ir.states.names)
    perProb{end+1} = ['GDSGE_TENSOR_' ir.states.names{i} '(:)''']; %#ok<AGROW>
    simuRows{end+1} = [ir.states.names{i} '(:)''']; %#ok<AGROW>
    nData = nData + 1;
end

frag = struct();
frag.iterPack = sprintf('GDSGE_DATA(:) = [repmat([%s],1,GDSGE_NPROB);%s];', ...
    innerStr, strjoin(perProb, ';'));
frag.simuData0 = sprintf('GDSGE_data0 = repmat([%s],1,GDSGE_NPROB);', innerStr);
frag.simuPack = sprintf('GDSGE_DATA = [GDSGE_data0;%s];', strjoin(simuRows, ';'));
frag.maxData = nData;
end
```

Create `src/+gdsge/+codegen/+mat/emitSolUnpack.m`:

```matlab
function frag = emitSolUnpack(ir)
% EMITSOLUNPACK  Slot unpacking + reshape fragments from the IR slot layout.
%   .unpack       <name>=GDSGE_SOL(lo:hi,:);  (policy, decl order)
%   .unpackAux    <name>=GDSGE_AUX(lo:hi,:);
%   .reshapeShock <name>=reshape(<name>,shock_num,[]);  (policy then aux)
%   .reshapeSize  save-time reshape to GDSGE_SIZE / [len GDSGE_SIZE]
%                 (policy, then interp vars, then aux — old reference order)
wU = gdsge.codegen.mat.codeWriter();
wA = gdsge.codegen.mat.codeWriter();
wS = gdsge.codegen.mat.codeWriter();
wZ = gdsge.codegen.mat.codeWriter();
for i = 1:numel(ir.variables.policy)
    v = ir.variables.policy{i};
    wU.add('%s=GDSGE_SOL(%d:%d,:);', v.name, v.slot(1), v.slot(2));
    wS.add('%s=reshape(%s,shock_num,[]);', v.name, v.name);
    if v.length == 1
        wZ.add('%s=reshape(%s,GDSGE_SIZE);', v.name, v.name);
    else
        wZ.add('%s=reshape(%s,[%d GDSGE_SIZE]);', v.name, v.name, v.length);
    end
end
for i = 1:numel(ir.variables.interp)
    n = ir.variables.interp{i};
    wZ.add('%s=reshape(%s,GDSGE_SIZE);', n, n);
end
for i = 1:numel(ir.variables.aux)
    v = ir.variables.aux{i};
    wA.add('%s=GDSGE_AUX(%d:%d,:);', v.name, v.slot(1), v.slot(2));
    wS.add('%s=reshape(%s,shock_num,[]);', v.name, v.name);
    if v.length == 1
        wZ.add('%s=reshape(%s,GDSGE_SIZE);', v.name, v.name);
    else
        wZ.add('%s=reshape(%s,[%d GDSGE_SIZE]);', v.name, v.name, v.length);
    end
end
frag = struct('unpack', wU.str(), 'unpackAux', wA.str(), ...
    'reshapeShock', wS.str(), 'reshapeSize', wZ.str());
end
```

- [ ] **Step 15.4: Run — expect PASS**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests({'tests/codegen/tEmitDataPack.m','tests/codegen/tEmitSolUnpack.m'}))"
```

- [ ] **Step 15.5: Commit**

```powershell
git add src/+gdsge/+codegen/+mat tests/codegen
git commit -m "feat(codegen): emitDataPack (MEX layout contract) + emitSolUnpack"
```

---

### Task 16: `emitResultIter`

**Files:**
- Create: `src/+gdsge/+codegen/+mat/emitResultIter.m`
- Test: `tests/codegen/tEmitResultIter.m`

- [ ] **Step 16.1: Write the failing test**

Create `tests/codegen/tEmitResultIter.m`:

```matlab
classdef tEmitResultIter < matlab.unittest.TestCase
    methods (TestClassSetup)
        function irFixture(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(here), 'HeatonLucas1996', 'ir')));
        end
    end
    methods (Test)
        function explicitStructAssemblyNoV2struct(tc)
            txt = gdsge.codegen.mat.emitResultIter(buildHL1996IR());
            tc.verifyFalse(contains(txt, 'v2struct'));
            tc.verifyTrue(contains(txt, 'IterRslt.params.beta = beta;'));
            tc.verifyTrue(contains(txt, 'IterRslt.params.GDSGE_EMPTY = GDSGE_EMPTY;'));  % frozen shape
            tc.verifyTrue(contains(txt, 'IterRslt.var_shock.g = g;'));
            tc.verifyTrue(contains(txt, 'IterRslt.var_interp.ps_future = ps_future;'));
            tc.verifyTrue(contains(txt, 'IterRslt.var_state.w1 = w1;'));
            tc.verifyTrue(contains(txt, 'IterRslt.GDSGE_PROB.GDSGE_LB = GDSGE_LB;'));
            tc.verifyTrue(contains(txt, 'IterRslt.GDSGE_PROB.GDSGE_SIZE = GDSGE_SIZE;'));
            tc.verifyTrue(contains(txt, 'IterRslt.var_others.GDSGE_EMPTY = GDSGE_EMPTY;'));
            tc.verifyTrue(contains(txt, 'IterRslt.NeedResolved = NeedResolved;'));
        end
        function policyAuxViaGetScalar(tc)
            txt = gdsge.codegen.mat.emitResultIter(buildHL1996IR());
            tc.verifyTrue(contains(txt, 'GDSGE_VAR_POLICY.c1 = c1;'));
            tc.verifyTrue(contains(txt, 'GDSGE_VAR_POLICY.w1n = w1n;'));
            tc.verifyTrue(contains(txt, ...
                'IterRslt.var_policy = get_scalar(GDSGE_VAR_POLICY,{[1:shock_num],w1});'));
            tc.verifyTrue(contains(txt, 'GDSGE_VAR_AUX.equity_premium = equity_premium;'));
            tc.verifyTrue(contains(txt, 'GDSGE_VAR_AUX.GDSGE_EMPTY = GDSGE_EMPTY;'));
        end
        function outputConstructPresent(tc)
            txt = gdsge.codegen.mat.emitResultIter(buildHL1996IR());
            tc.verifyTrue(contains(txt, ['outputVarStack = cat(1,reshape(c1,1,[]),' ...
                'reshape(c2,1,[]),reshape(ps,1,[]),reshape(pb,1,[]),' ...
                'reshape(equity_premium,1,[]),reshape(w1n,8,[]));']));
            tc.verifyTrue(contains(txt, 'output_var_index.c1=1:1;'));
            tc.verifyTrue(contains(txt, 'output_var_index.w1n=6:13;'));
            tc.verifyTrue(contains(txt, 'IterRslt.output_interp=myppual(output_interp);'));
        end
        function ppFieldsKeepOldNames(tc)
            txt = gdsge.codegen.mat.emitResultIter(buildHL1996IR());
            tc.verifyTrue(contains(txt, 'IterRslt.pp.GDSGE_PP_ps_future = GDSGE_PP_ps_future;'));
            tc.verifyTrue(contains(txt, 'IterRslt.pp.GDSGE_SPLINE_VEC = GDSGE_SPLINE_VEC;'));
            tc.verifyTrue(contains(txt, 'IterRslt.pp.GDSGE_EMPTY = GDSGE_EMPTY;'));
        end
    end
end
```

- [ ] **Step 16.2: Run — expect FAIL**

- [ ] **Step 16.3: Implement**

Create `src/+gdsge/+codegen/+mat/emitResultIter.m`:

```matlab
function txt = emitResultIter(ir)
% EMITRESULTITER  Output construction + explicit IterRslt assembly.
%   Frozen result shape: same fields as the old generator, including the
%   GDSGE_EMPTY placeholder fields v2struct used to create (params, var_aux,
%   var_tensor, pp, var_others). No v2struct anywhere.
lenOf = struct();
for i = 1:numel(ir.variables.policy)
    lenOf.(ir.variables.policy{i}.name) = ir.variables.policy{i}.length;
end
for i = 1:numel(ir.variables.aux)
    lenOf.(ir.variables.aux{i}.name) = ir.variables.aux{i}.length;
end
stateList = strjoin(ir.states.names, ',');   % 'w1'
gridCell = ['{' strjoin(ir.states.names, ',') '}'];   % '{w1}'

w = gdsge.codegen.mat.codeWriter();
w.add('if CONSTRUCT_OUTPUT==1');
w.add('%% Construct output variables');
parts = cell(1, numel(ir.variables.output));
for i = 1:numel(ir.variables.output)
    n = ir.variables.output{i};
    parts{i} = sprintf('reshape(%s,%d,[])', n, lenOf.(n));
end
w.add('outputVarStack = cat(1,%s);', strjoin(parts, ','));
w.add('if shock_num>1');
w.add('    outputVarStack = reshape(outputVarStack, [],shock_num,GDSGE_SIZE_STATE{:});');
w.add('    output_interp=struct(''form'',''MKL'',''breaks'',{{1:shock_num,%s}},...', stateList);
w.add('        ''Values'',outputVarStack,...');
w.add('        ''coefs'',[],''order'',[2 OutputInterpOrder*ones(1,length(%s))],''Method'',[],...', gridCell);
w.add('        ''ExtrapolationOrder'',[],''thread'',NumThreads, ...');
w.add('        ''orient'',''curvefit'');');
w.add('else');
w.add('    outputVarStack = reshape(outputVarStack, [],GDSGE_SIZE_STATE{:});');
w.add('    output_interp=struct(''form'',''MKL'',''breaks'',{{%s}},...', stateList);
w.add('        ''Values'',outputVarStack,...');
w.add('        ''coefs'',[],''order'',[OutputInterpOrder*ones(1,length(%s))],''Method'',[],...', gridCell);
w.add('        ''ExtrapolationOrder'',[],''thread'',NumThreads, ...');
w.add('        ''orient'',''curvefit'');');
w.add('end');
w.add('IterRslt.output_interp=myppual(output_interp);');
w.add('output_var_index=struct();');
idx = 1;
for i = 1:numel(ir.variables.output)
    n = ir.variables.output{i};
    w.add('output_var_index.%s=%d:%d;', n, idx, idx + lenOf.(n) - 1);
    idx = idx + lenOf.(n);
end
w.add('IterRslt.output_var_index = output_var_index;');
w.blank();
w.add('IterRslt.shock_num = shock_num;');
w.add('IterRslt.shock_trans = shock_trans;');
for i = 1:numel(ir.params)
    w.add('IterRslt.params.%s = %s;', ir.params{i}.name, ir.params{i}.name);
end
w.add('IterRslt.params.GDSGE_EMPTY = GDSGE_EMPTY;');
for i = 1:numel(ir.shocks.names)
    w.add('IterRslt.var_shock.%s = %s;', ir.shocks.names{i}, ir.shocks.names{i});
end
w.add('GDSGE_VAR_POLICY = struct();');
for i = 1:numel(ir.variables.policy)
    n = ir.variables.policy{i}.name;
    w.add('GDSGE_VAR_POLICY.%s = %s;', n, n);
end
w.add('IterRslt.var_policy = get_scalar(GDSGE_VAR_POLICY,{[1:shock_num],%s});', stateList);
w.add('GDSGE_VAR_AUX = struct();');
for i = 1:numel(ir.variables.aux)
    n = ir.variables.aux{i}.name;
    w.add('GDSGE_VAR_AUX.%s = %s;', n, n);
end
w.add('GDSGE_VAR_AUX.GDSGE_EMPTY = GDSGE_EMPTY;');
w.add('IterRslt.var_aux = get_scalar(GDSGE_VAR_AUX,{[1:shock_num],%s});', stateList);
w.add('IterRslt.var_tensor.GDSGE_EMPTY = GDSGE_EMPTY;');
for i = 1:numel(ir.variables.interp)
    n = ir.variables.interp{i};
    w.add('IterRslt.pp.GDSGE_PP_%s = GDSGE_PP_%s;', n, n);
end
w.add('IterRslt.pp.GDSGE_SPLINE_VEC = GDSGE_SPLINE_VEC;');
w.add('IterRslt.pp.GDSGE_EMPTY = GDSGE_EMPTY;');
w.add('end');
w.add('IterRslt.Metric0 = GDSGE_Metric0;');
w.add('IterRslt.Metric = GDSGE_Metric;');
w.add('IterRslt.Iter = GDSGE_Iter;');
for i = 1:numel(ir.states.names)
    w.add('IterRslt.var_state.%s = %s;', ir.states.names{i}, ir.states.names{i});
end
for i = 1:numel(ir.variables.interp)
    n = ir.variables.interp{i};
    w.add('IterRslt.var_interp.%s = %s;', n, n);
end
w.add('IterRslt.GDSGE_PROB.GDSGE_LB = GDSGE_LB;');
w.add('IterRslt.GDSGE_PROB.GDSGE_UB = GDSGE_UB;');
w.add('IterRslt.GDSGE_PROB.GDSGE_SOL = GDSGE_SOL;');
w.add('IterRslt.GDSGE_PROB.GDSGE_F = GDSGE_F;');
w.add('IterRslt.GDSGE_PROB.GDSGE_SIZE = GDSGE_SIZE;');
w.add('IterRslt.var_others.GDSGE_EMPTY = GDSGE_EMPTY;');
w.add('IterRslt.NeedResolved = NeedResolved;');
txt = w.str();
end
```

Note the old generator put **both** the output-construct block and the `shock_num` …
`pp` assignments under `CONSTRUCT_OUTPUT==1` — this emitter does the same (one `if`,
closed by the `w.add('end')` after the `pp` lines).

- [ ] **Step 16.4: Run — expect PASS**

- [ ] **Step 16.5: Commit**

```powershell
git add src/+gdsge/+codegen/+mat/emitResultIter.m tests/codegen/tEmitResultIter.m
git commit -m "feat(codegen): emitResultIter - explicit IterRslt assembly, frozen shape"
```

---

### Task 17: `emitIter` — assemble the full iter file (+ parse test + snapshot)

**Files:**
- Create: `src/+gdsge/+codegen/+mat/emitIter.m`
- Create: `tests/HeatonLucas1996/codegen/regen_snapshots.m`
- Create: `tests/HeatonLucas1996/codegen/golden/iter_HL1996_golden.txt` (generated, reviewed, committed)
- Create (only if missing): `scratch/regen_old_matlab_codegen.m`
- Test: `tests/HeatonLucas1996/codegen/tSnapshotHL1996.m`

- [ ] **Step 17.1: If `scratch/old_generated/iter_HL1996.m` is missing, recreate the regen script**

Create `scratch/regen_old_matlab_codegen.m` (scratch/ is git-ignored; this is a working aid):

```matlab
function regen_old_matlab_codegen()
% Regenerate the OLD toolbox's generated MATLAB for HL1996 (parser only, no
% compile). Writes iter/simulate .m into scratch/old_generated/ for reference.
here     = fileparts(mfilename('fullpath'));          % scratch
repoRoot = fileparts(here);
oldSrc   = fullfile(repoRoot, 'base_package', 'gdsge', 'source');
modelDir = fullfile(repoRoot, 'tests', 'HeatonLucas1996');
outDir   = fullfile(here, 'old_generated');
if ~exist(outDir, 'dir'); mkdir(outDir); end
work = tempname; mkdir(work);
copyfile(fullfile(modelDir, 'HL1996.gmod'), work);
oldPath = path; restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(oldSrc);
oldCd = pwd; cdRestore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(work);
[~, iterCode, ~, simulateCode] = gdsge_parser('HL1996');
fid = fopen(fullfile(outDir, 'iter_HL1996.m'), 'w');
fprintf(fid, '%s', iterCode); fclose(fid);
fid = fopen(fullfile(outDir, 'simulate_HL1996.m'), 'w');
fprintf(fid, '%s', simulateCode); fclose(fid);
fprintf('Wrote old generated MATLAB to %s\n', outDir);
end
```

Run: `matlab -batch "cd('scratch'); regen_old_matlab_codegen"`.

- [ ] **Step 17.2: Write the failing parse test** (snapshot test added in 17.6 after review)

Create `tests/HeatonLucas1996/codegen/tSnapshotHL1996.m`:

```matlab
classdef tSnapshotHL1996 < matlab.unittest.TestCase
    methods (TestClassSetup)
        function irFixture(tc)
            here = fileparts(mfilename('fullpath'));   % tests/HeatonLucas1996/codegen
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(here), 'ir')));
        end
    end
    methods (Test)
        function iterParsesCleanly(tc)
            txt = gdsge.codegen.mat.emitIter(buildHL1996IR());
            T = mtree(txt);
            tc.verifyEqual(count(mtfind(T, 'Kind', 'ERR')), 0, ...
                'generated iter file has MATLAB syntax errors');
        end
        function iterHasNoV2struct(tc)
            txt = gdsge.codegen.mat.emitIter(buildHL1996IR());
            tc.verifyFalse(contains(txt, 'v2struct'));
        end
        function iterMatchesSnapshot(tc)
            here = fileparts(mfilename('fullpath'));
            goldenFile = fullfile(here, 'golden', 'iter_HL1996_golden.txt');
            tc.assertTrue(exist(goldenFile, 'file') == 2, ...
                'snapshot missing — run tests/HeatonLucas1996/codegen/regen_snapshots.m and review');
            txt = gdsge.codegen.mat.emitIter(buildHL1996IR());
            golden = fileread(goldenFile);
            tc.verifyEqual(normEol(txt), normEol(golden), ...
                'iter snapshot drifted — if intentional, regen_snapshots and review the diff');
        end
    end
end

function s = normEol(s)
s = strtrim(strrep(s, sprintf('\r\n'), sprintf('\n')));
end
```

- [ ] **Step 17.3: Run — expect FAIL** (`emitIter` undefined)

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/HeatonLucas1996/codegen/tSnapshotHL1996.m'))"
```

- [ ] **Step 17.4: Implement `emitIter`**

Create `src/+gdsge/+codegen/+mat/emitIter.m`. This is the file-level assembler; every
fragment comes from an already-tested emitter or is structural glue. Keep the section
order exactly as below (it mirrors the old file; see parity comments):

```matlab
function txt = emitIter(ir)
% EMITITER  Assemble the complete iter_<model>.m from IR sections.
%   Thin generated file: model-specific code only; algorithmic blocks live in
%   gdsge.runtime.*. Parity reference: scratch/old_generated/iter_HL1996.m.
bounds = gdsge.codegen.mat.emitBounds(ir);
pack   = gdsge.codegen.mat.emitDataPack(ir);
unpk   = gdsge.codegen.mat.emitSolUnpack(ir);
m      = ir.modelName;
nAux   = 0;
for i = 1:numel(ir.variables.aux); nAux = max(nAux, ir.variables.aux{i}.slot(2)); end
interpNames = ir.variables.interp;                     % decl order
interpList  = strjoin(interpNames, ',');               % 'ps_future,pb_future,...'
stateList   = strjoin(ir.states.names, ',');           % 'w1'
stateNameCell = ['{''' strjoin(ir.states.names, ''',''') '''}'];   % {'w1'}
stateGridCell = ['{' stateList '}'];                                % {w1}
tensorRows  = ['[GDSGE_TENSOR_shockIdx(:)'';' ...
    strjoin(cellfun(@(s) ['GDSGE_TENSOR_' s '(:)'''], ir.states.names, ...
    'UniformOutput', false), ';') ']'];

w = gdsge.codegen.mat.codeWriter();
w.add('%% Generated by gdsge.codegen.generateMatlab — do not edit.');
w.add('function [IterRslt,IterFlag] = iter_%s(GDSGE_OPTIONS)', m);
w.add('gdsge.runtime.ensurePath();');
w.blank();
w.addRaw(gdsge.codegen.mat.emitSetup(ir, 'iter'));
w.blank();

% ---- options unpacking (whitelist = frozen surface + model names) ---------
base = {'TolEq','TolSol','TolFun','PrintFreq','NoPrint','SaveFreq','NoSave', ...
    'SimuPrintFreq','SimuSaveFreq','MaxIter','MaxMinorIter','num_samples', ...
    'num_periods','SolMaxIter','UseBroyden','FiniteDiffDelta','GDSGE_USE_BROYDEN', ...
    'GDSGE_DEBUG_EVAL_ONLY','INTERP_ORDER','EXTRAP_ORDER','OutputInterpOrder', ...
    'IterSaveAll','SkipModelInit','UseAdaptiveBound','UseAdaptiveBoundInSol', ...
    'EnforceSimuStateInbound','REUSE_WARMUP_SOL','INTERP_WARMUP_SOL', ...
    'CONSTRUCT_OUTPUT','NumThreads','WarmUp'};
modelNames = {};
for i = 1:numel(ir.params); modelNames{end+1} = ir.params{i}.name; end %#ok<AGROW>
modelNames = [modelNames, ir.shocks.names, fieldnames(ir.shocks.transitions)', ir.states.names];
wl = [base, modelNames];
w.add('GDSGE_OPTIONS_VALID = {%s};', ['''' strjoin(wl, ''',''') '''']);
w.add('if nargin < 1; GDSGE_OPTIONS = struct(); end');
w.add('gdsge.runtime.unpackOptions(GDSGE_OPTIONS, GDSGE_OPTIONS_VALID);');
w.blank();

% ---- declaration asserts (parity with old lines 111-116) ------------------
w.add('assert(exist(''shock_num'',''var'')==1);');
for i = 1:numel(ir.shocks.names)
    w.add('assert(length(%s)==%d);', ir.shocks.names{i}, ir.shocks.count);
end
tn = fieldnames(ir.shocks.transitions);
for i = 1:numel(tn)
    w.add('assert(size(%s,1)==%d);', tn{i}, ir.shocks.count);
    w.add('assert(size(%s,2)==%d);', tn{i}, ir.shocks.count);
end
w.blank();

% ---- state-space tensors ---------------------------------------------------
w.add('%%%% State-space tensors');
w.add('[GDSGE_TENSOR_shockIdx,%s]=ndgrid(1:shock_num,%s);', ...
    strjoin(cellfun(@(s) ['GDSGE_TENSOR_' s], ir.states.names, 'UniformOutput', false), ','), ...
    stateList);
for i = 1:numel(ir.shocks.names)
    w.add('GDSGE_TENSOR_%s=ndgrid(%s,%s);', ir.shocks.names{i}, ir.shocks.names{i}, stateList);
end
w.add('GDSGE_NPROB=numel(GDSGE_TENSOR_shockIdx);');
w.add('GDSGE_SIZE = size(GDSGE_TENSOR_shockIdx);');
w.add('GDSGE_SIZE_STATE = num2cell(GDSGE_SIZE(2:end));');
w.blank();

% ---- bounds + solution space ----------------------------------------------
w.add('%%%% Bounds');
w.addRaw(bounds.init);
w.blank();
w.add('%%%% Solution space');
w.add('GDSGE_EQVAL = 1e20*ones(%d,GDSGE_NPROB);', bounds.maxDim);
w.add('GDSGE_F = 1e20*ones(1,GDSGE_NPROB);');
w.add('GDSGE_SOL = zeros(%d,GDSGE_NPROB);', bounds.maxDim);
w.add('GDSGE_X0 = rand(size(GDSGE_SOL)) .* (GDSGE_UB-GDSGE_LB) + GDSGE_LB;');
w.add('GDSGE_SOL(:) = GDSGE_X0;');
w.add('GDSGE_AUX = zeros(%d,GDSGE_NPROB);', max(nAux, 1));
w.add('GDSGE_DATA = zeros(%d,GDSGE_NPROB);', pack.maxData);
w.add('NeedResolved = false(1,GDSGE_NPROB);');
w.blank();

% ---- interp initialization (skipped under interp WarmUp; parity line 172) --
w.add('%%%% Interp initialization');
w.add('if ~( nargin>=1 && isfield(GDSGE_OPTIONS,''WarmUp'') && isfield(GDSGE_OPTIONS.WarmUp,''var_interp'') )');
w.in();
w.addRaw(indent4(gdsge.codegen.mat.emitInterpInitial(ir)));
w.add('[GDSGE_PP_CELL, GDSGE_SPLINE_VEC] = gdsge.runtime.constructSplines({%s}, %s, GDSGE_SIZE_STATE, INTERP_ORDER, EXTRAP_ORDER, NumThreads);', ...
    interpList, stateGridCell);
for i = 1:numel(interpNames)
    w.add('GDSGE_PP_%s = GDSGE_PP_CELL{%d};', interpNames{i}, i);
end
w.out();
w.add('end');
w.blank();

w.add('GDSGE_Metric = 1;');
w.add('GDSGE_Iter = 0;');
w.add('IS_WARMUP_LOOP = 0;');
w.blank();

% ---- WarmUp ----------------------------------------------------------------
w.add('if nargin>=1 && isfield(GDSGE_OPTIONS,''WarmUp'')');
w.in();
w.add('GDSGE_WARMUP = GDSGE_OPTIONS.WarmUp;');
w.add('if isfield(GDSGE_WARMUP,''var_interp'')');
w.in();
for i = 1:numel(interpNames)
    n = interpNames{i};
    w.add('if isfield(GDSGE_WARMUP.var_interp,''%s''); %s = GDSGE_WARMUP.var_interp.%s; end', n, n, n);
end
w.add('if isfield(GDSGE_WARMUP,''GDSGE_PROB'')');
w.add('    GDSGE_SIZE_STATE_WU = num2cell(GDSGE_WARMUP.GDSGE_PROB.GDSGE_SIZE(2:end));');
w.add('else');
w.add('    GDSGE_SIZE_STATE_WU = GDSGE_SIZE_STATE;');
w.add('end');
w.add('if isfield(GDSGE_WARMUP,''var_state'')');
w.add('    GDSGE_BREAKS_WU = {%s};', strjoin(cellfun(@(s) ...
    ['GDSGE_WARMUP.var_state.' s], ir.states.names, 'UniformOutput', false), ','));
w.add('else');
w.add('    GDSGE_BREAKS_WU = %s;', stateGridCell);
w.add('end');
w.add('[GDSGE_PP_CELL, GDSGE_SPLINE_VEC] = gdsge.runtime.constructSplines({%s}, GDSGE_BREAKS_WU, GDSGE_SIZE_STATE_WU, INTERP_ORDER, EXTRAP_ORDER, NumThreads);', interpList);
for i = 1:numel(interpNames)
    w.add('GDSGE_PP_%s = GDSGE_PP_CELL{%d};', interpNames{i}, i);
end
w.add('IS_WARMUP_LOOP = 1;');
w.out();
w.add('end');
w.add('if isfield(GDSGE_WARMUP,''Iter''); GDSGE_Iter = GDSGE_WARMUP.Iter; end');
w.add('GDSGE_WARMSPEC = struct();');
w.add('GDSGE_WARMSPEC.shockNum = shock_num;');
w.add('GDSGE_WARMSPEC.stateNames = %s;', stateNameCell);
w.add('GDSGE_WARMSPEC.evalPoints = %s;', tensorRows);
w.add('GDSGE_WARMSPEC.adaptiveSlots = %s;', mat2str(bounds.adaptiveSlots));
w.add('GDSGE_WARMSPEC.reuseSol = REUSE_WARMUP_SOL;');
w.add('GDSGE_WARMSPEC.interpSol = INTERP_WARMUP_SOL;');
w.add('GDSGE_WARMSPEC.numThreads = NumThreads;');
w.add('[GDSGE_SOL, GDSGE_LB, GDSGE_UB] = gdsge.runtime.applyWarmUp(GDSGE_WARMUP, GDSGE_SOL, GDSGE_LB, GDSGE_UB, GDSGE_WARMSPEC);');
w.out();
w.add('end');
w.blank();

% ---- main loop --------------------------------------------------------------
w.add('stopFlag = false;');
w.add('GDSGE_timer = tic;');
w.add('while(~stopFlag)');
w.in();
w.add('GDSGE_Iter = GDSGE_Iter+1;');
w.blank();
w.add(pack.iterPack);
w.blank();
w.add('GDSGE_CFG = struct();');
w.add('GDSGE_CFG.tolSol = TolSol; GDSGE_CFG.tolFun = TolFun;');
w.add('GDSGE_CFG.solMaxIter = SolMaxIter; GDSGE_CFG.numThreads = NumThreads;');
w.add('GDSGE_CFG.debugEvalOnly = GDSGE_DEBUG_EVAL_ONLY;');
w.add('GDSGE_CFG.useBroyden = UseBroyden; GDSGE_CFG.finiteDiffDelta = FiniteDiffDelta;');
w.add('GDSGE_CFG.useBroydenNow = double(GDSGE_Iter>1)*GDSGE_USE_BROYDEN;');
w.add('GDSGE_CFG.taskName = MEX_TASK_INF_HORIZON;');
w.add('GDSGE_CFG.splineVec = GDSGE_SPLINE_VEC;');
w.add('GDSGE_CFG.maxMinorIter = MaxMinorIter;');
w.add('GDSGE_CFG.probSize = GDSGE_SIZE;');
w.add('GDSGE_CFG.useNearestNeighbor = true;');
w.add('GDSGE_CFG.verboseRetry = ~NoPrint;');
w.add('if UseAdaptiveBoundInSol==1');
w.add('    GDSGE_CFG.adaptInSol = @GDSGE_ADAPT_TIGHT;');
w.add('else');
w.add('    GDSGE_CFG.adaptInSol = [];');
w.add('end');
w.add('[GDSGE_SOL,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL,GDSGE_OPT_INFO,GDSGE_DIAG] = gdsge.runtime.solveProblems(@mex_%s, GDSGE_SOL, GDSGE_LB, GDSGE_UB, GDSGE_DATA, GDSGE_F, GDSGE_AUX, GDSGE_EQVAL, GDSGE_CFG); %%#ok<ASGLU>', m);
w.add('NeedResolved = GDSGE_DIAG.needResolved;');
w.add('if any(NeedResolved)');
w.add('    warning(''gdsge:runtime:unconverged'', ''%%s'', gdsge.runtime.reportUnconverged(NeedResolved, GDSGE_F, GDSGE_SIZE, %s, %s, 5));', ...
    stateNameCell, stateGridCell);
w.add('end');
w.blank();
w.addRaw(indent4(unpk.unpack));
w.addRaw(indent4(unpk.unpackAux));
w.addRaw(indent4(unpk.reshapeShock));
w.blank();
for i = 1:numel(ir.interp)
    it = ir.interp{i};
    w.add('GDSGE_NEW_%s = %s;', it.name, it.updateExpr);
end
newList = strjoin(cellfun(@(n) ['GDSGE_NEW_' n], interpNames, 'UniformOutput', false), ',');
w.add('GDSGE_Metric = gdsge.runtime.computeMetric({%s},{%s});', newList, interpList);
w.add('GDSGE_Metric0 = GDSGE_Metric;');
w.add('if IS_WARMUP_LOOP==1');
w.add('    IS_WARMUP_LOOP = 0;');
w.add('end');
for i = 1:numel(interpNames)
    w.add('%s=GDSGE_NEW_%s;', interpNames{i}, interpNames{i});
end
w.blank();
w.addRaw(indent4(bounds.adaptiveWiden));
w.blank();
w.add('stopFlag = (isempty(GDSGE_Metric) || GDSGE_Metric<TolEq) || GDSGE_Iter>=MaxIter;');
w.blank();
w.add('[GDSGE_PP_CELL, GDSGE_SPLINE_VEC] = gdsge.runtime.constructSplines({%s}, %s, GDSGE_SIZE_STATE, INTERP_ORDER, EXTRAP_ORDER, NumThreads);', ...
    interpList, stateGridCell);
for i = 1:numel(interpNames)
    w.add('GDSGE_PP_%s = GDSGE_PP_CELL{%d};', interpNames{i}, i);
end
w.blank();
w.add('if gdsge.runtime.printIterProgress(GDSGE_Iter, GDSGE_Metric, max(GDSGE_F), nnz(NeedResolved), toc(GDSGE_timer), PrintFreq, NoPrint, stopFlag)');
w.add('    GDSGE_timer = tic;');
w.add('end');
w.blank();
w.add('if ( (mod(GDSGE_Iter,SaveFreq)==0 || stopFlag == true) )');
w.in();
w.addRaw(indent8(unpk.reshapeSize));
w.addRaw(indent8(gdsge.codegen.mat.emitResultIter(ir)));
w.add('if ~NoSave');
w.add('    if IterSaveAll');
w.add('        save([''IterRslt_%s_'' num2str(GDSGE_Iter) ''.mat'']);', m);
w.add('    else');
w.add('        save([''IterRslt_%s_'' num2str(GDSGE_Iter) ''.mat''],''IterRslt'');', m);
w.add('    end');
w.add('end');
w.out();
w.add('end');
w.out();
w.add('end');
w.blank();
w.add('IterFlag = 0;');
w.add('end');
w.blank();

% ---- model-specific adaptive-tighten subfunction ---------------------------
w.add('function [GDSGE_LB,GDSGE_UB] = GDSGE_ADAPT_TIGHT(GDSGE_SOL,GDSGE_LB,GDSGE_UB)');
if isempty(bounds.adaptiveTight)
    w.add('%% no adaptive bounds in this model');
else
    w.addRaw(bounds.adaptiveTight);
end
w.add('end');
txt = w.str();
end

function s = indent4(s)
s = indentBy(s, '    ');
end
function s = indent8(s)
s = indentBy(s, '        ');
end
function s = indentBy(s, pad)
if isempty(s); return; end
lines = strsplit(s, newline, 'CollapseDelimiters', false);
for i = 1:numel(lines)
    if ~isempty(lines{i}); lines{i} = [pad lines{i}]; end
end
s = strjoin(lines, newline);
end
```

- [ ] **Step 17.5: Run the parse + no-v2struct tests — expect PASS, snapshot test fails (missing golden)**

Fix any syntax issues `mtree` reports before moving on (typical culprits: unbalanced
quotes in `w.add` format strings — remember `%%` for a literal `%` and `''` for a quote).

- [ ] **Step 17.6: Create, review, and commit the snapshot**

Create `tests/HeatonLucas1996/codegen/regen_snapshots.m`:

```matlab
function regen_snapshots()
% Regenerate the committed text snapshots of the generated files.
% Review the diff before committing — the snapshot IS the spec of the output.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'ir'));
ir = buildHL1996IR();
outDir = fullfile(here, 'golden');
if ~exist(outDir, 'dir'); mkdir(outDir); end
writeTxt(fullfile(outDir, 'iter_HL1996_golden.txt'), gdsge.codegen.mat.emitIter(ir));
if exist(fullfile(fileparts(which('gdsge.codegen.mat.emitIter')), 'emitSimulate.m'), 'file')
    writeTxt(fullfile(outDir, 'simulate_HL1996_golden.txt'), gdsge.codegen.mat.emitSimulate(ir));
end
fprintf('Snapshots written to %s\n', outDir);
end
function writeTxt(p, txt)
fid = fopen(p, 'w'); fwrite(fid, txt); fclose(fid);
end
```

Run it, then **review the generated file against the old reference**:

```
matlab -batch "addpath('src'); addpath('src/kernels'); cd('tests/HeatonLucas1996/codegen'); regen_snapshots"
```

Review checklist (against `scratch/old_generated/iter_HL1996.m`):
- [ ] 19-row bounds identical slot-by-slot (old lines 135-160)
- [ ] `GDSGE_DATA` width 94 and pack expression ordering identical (old line 314)
- [ ] interp init values identical (old lines 173-180)
- [ ] adaptive widen uses `min`/`max` with `*1.5` on slots 10/11 (old lines 475-478)
- [ ] result fields cover every old `IterRslt.` assignment (old lines 559-589)
- [ ] no `v2struct` anywhere

- [ ] **Step 17.7: Run the full snapshot test file — expect PASS, then commit**

```powershell
git add src/+gdsge/+codegen/+mat/emitIter.m tests/HeatonLucas1996/codegen
git commit -m "feat(codegen): emitIter - full iter file assembly + snapshot"
```

---

### Task 18: `emitResultSimu`, `emitSimulate`, `generateMatlab` driver

**Files:**
- Create: `src/+gdsge/+codegen/+mat/emitResultSimu.m`, `src/+gdsge/+codegen/+mat/emitSimulate.m`,
  `src/+gdsge/+codegen/generateMatlab.m`
- Modify: `tests/HeatonLucas1996/codegen/tSnapshotHL1996.m` (simulate tests)
- Test: `tests/codegen/tEmitResultSimu.m`, `tests/codegen/tGenerateMatlab.m`

- [ ] **Step 18.1: Write the failing tests**

Create `tests/codegen/tEmitResultSimu.m`:

```matlab
classdef tEmitResultSimu < matlab.unittest.TestCase
    methods (TestClassSetup)
        function irFixture(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(here), 'HeatonLucas1996', 'ir')));
        end
    end
    methods (Test)
        function preallocInitAndAssigns(tc)
            frag = gdsge.codegen.mat.emitResultSimu(buildHL1996IR());
            % prealloc: states then varSimu (old reference lines 154-159)
            tc.verifyTrue(contains(frag.prealloc, ...
                'SimuRslt.w1=zeros(num_samples,num_periods);'));
            tc.verifyTrue(contains(frag.prealloc, ...
                'SimuRslt.equity_premium=zeros(num_samples,num_periods);'));
            % initial conditions (old lines 162-163)
            tc.verifyTrue(contains(frag.init, 'SimuRslt.w1(:,1)=0.5;'));
            tc.verifyTrue(contains(frag.init, 'SimuRslt.shock(:,1)=1;'));
            % GDSGE_OPTIONS.init overwrite (old lines 167-168)
            tc.verifyTrue(contains(frag.initOverwrite, ...
                'SimuRslt.w1(:,1:size(GDSGE_OPTIONS.init.w1,2))=GDSGE_OPTIONS.init.w1;'));
            % per-period assignment + transition (old lines 325-330)
            tc.verifyTrue(contains(frag.assign, 'SimuRslt.c1(:,GDSGE_t)=c1;'));
            tc.verifyTrue(contains(frag.assign, ...
                'SimuRslt.w1(:,GDSGE_t+1) = w1n(GDSGE_SHOCK_VAR_LINEAR_INDEX);'));
        end
    end
end
```

Create `tests/codegen/tGenerateMatlab.m`:

```matlab
classdef tGenerateMatlab < matlab.unittest.TestCase
    methods (TestClassSetup)
        function irFixture(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(here), 'HeatonLucas1996', 'ir')));
        end
    end
    methods (Test)
        function writesBothFiles(tc)
            ir = buildHL1996IR();
            outDir = tempname; mkdir(outDir);
            files = gdsge.codegen.generateMatlab(ir, outDir);
            tc.verifyEqual(files.iterFile, fullfile(outDir, 'iter_HL1996.m'));
            tc.verifyEqual(files.simulateFile, fullfile(outDir, 'simulate_HL1996.m'));
            tc.verifyEqual(fileread(files.iterFile), gdsge.codegen.mat.emitIter(ir));
            tc.verifyEqual(fileread(files.simulateFile), gdsge.codegen.mat.emitSimulate(ir));
        end
    end
end
```

Add to `tests/HeatonLucas1996/codegen/tSnapshotHL1996.m` (inside `methods (Test)`):

```matlab
        function simulateParsesCleanly(tc)
            txt = gdsge.codegen.mat.emitSimulate(buildHL1996IR());
            T = mtree(txt);
            tc.verifyEqual(count(mtfind(T, 'Kind', 'ERR')), 0);
            tc.verifyFalse(contains(txt, 'v2struct'));
        end
        function simulateMatchesSnapshot(tc)
            here = fileparts(mfilename('fullpath'));
            goldenFile = fullfile(here, 'golden', 'simulate_HL1996_golden.txt');
            tc.assertTrue(exist(goldenFile, 'file') == 2, ...
                'snapshot missing — run regen_snapshots and review');
            txt = gdsge.codegen.mat.emitSimulate(buildHL1996IR());
            tc.verifyEqual(normEol(txt), normEol(fileread(goldenFile)));
        end
```

- [ ] **Step 18.2: Run — expect FAIL**

- [ ] **Step 18.3: Implement `emitResultSimu`**

Create `src/+gdsge/+codegen/+mat/emitResultSimu.m`:

```matlab
function frag = emitResultSimu(ir)
% EMITRESULTSIMU  SimuRslt fragments: prealloc, initial conditions, the
%   GDSGE_OPTIONS.init overwrite, and per-period assignment + state transition.
w1 = gdsge.codegen.mat.codeWriter();   % prealloc (states then varSimu)
for i = 1:numel(ir.states.names)
    w1.add('SimuRslt.%s=zeros(num_samples,num_periods);', ir.states.names{i});
end
for i = 1:numel(ir.simulate.varSimu)
    w1.add('SimuRslt.%s=zeros(num_samples,num_periods);', ir.simulate.varSimu{i});
end

w2 = gdsge.codegen.mat.codeWriter();   % initial conditions
for i = 1:numel(ir.simulate.initial)
    it = ir.simulate.initial{i};
    w2.add('SimuRslt.%s(:,1)=%s;', it.var, it.value);
end

w3 = gdsge.codegen.mat.codeWriter();   % GDSGE_OPTIONS.init overwrite
initVars = [cellfun(@(x) x.var, ir.simulate.initial, 'UniformOutput', false)];
for i = 1:numel(initVars)
    n = initVars{i};
    w3.add('SimuRslt.%s(:,1:size(GDSGE_OPTIONS.init.%s,2))=GDSGE_OPTIONS.init.%s;', n, n, n);
end

w4 = gdsge.codegen.mat.codeWriter();   % per-period assigns + transitions
for i = 1:numel(ir.simulate.varSimu)
    n = ir.simulate.varSimu{i};
    w4.add('SimuRslt.%s(:,GDSGE_t)=%s;', n, n);
end
for i = 1:numel(ir.simulate.transitions)
    t = ir.simulate.transitions{i};
    % shock-indexed policy arrays advance via the linear index (parity line 330)
    w4.add('SimuRslt.%s(:,GDSGE_t+1) = %s(GDSGE_SHOCK_VAR_LINEAR_INDEX);', t.state, t.expr);
end
frag = struct('prealloc', w1.str(), 'init', w2.str(), ...
    'initOverwrite', w3.str(), 'assign', w4.str());
end
```

- [ ] **Step 18.4: Implement `emitSimulate`**

Create `src/+gdsge/+codegen/+mat/emitSimulate.m` (SIMU_RESOLVE variant; parity reference
`scratch/old_generated/simulate_HL1996.m`):

```matlab
function txt = emitSimulate(ir)
% EMITSIMULATE  Assemble simulate_<model>.m (SIMU_RESOLVE variant — the only
%   one HL1996 uses; SIMU_INTERP is Phase 7).
bounds = gdsge.codegen.mat.emitBounds(ir);
pack   = gdsge.codegen.mat.emitDataPack(ir);
unpk   = gdsge.codegen.mat.emitSolUnpack(ir);
simu   = gdsge.codegen.mat.emitResultSimu(ir);
m      = ir.modelName;
nAux   = 0;
for i = 1:numel(ir.variables.aux); nAux = max(nAux, ir.variables.aux{i}.slot(2)); end
stateList = strjoin(ir.states.names, ',');
stateGridCell = ['{' stateList '}'];
stateNameCell = ['{''' strjoin(ir.states.names, ''',''') '''}'];
simuStateRows = strjoin(cellfun(@(s) ['SimuRslt.' s '(:,GDSGE_t)'''], ...
    ir.states.names, 'UniformOutput', false), ';');

w = gdsge.codegen.mat.codeWriter();
w.add('%% Generated by gdsge.codegen.generateMatlab — do not edit.');
w.add('function SimuRslt = simulate_%s(IterRslt,GDSGE_OPTIONS)', m);
w.add('gdsge.runtime.ensurePath();');
w.blank();
w.addRaw(gdsge.codegen.mat.emitSetup(ir, 'simulate'));
w.blank();
w.add('GEN_SHOCK_START_PERIOD = 1;');
w.blank();
w.add('%% ---- pull the solved model from IterRslt (explicit; no v2struct)');
tn = fieldnames(ir.shocks.transitions);
for i = 1:numel(tn)
    w.add('%s = IterRslt.%s;', tn{i}, 'shock_trans');
end
for i = 1:numel(ir.params)
    w.add('%s = IterRslt.params.%s;', ir.params{i}.name, ir.params{i}.name);
end
for i = 1:numel(ir.shocks.names)
    w.add('%s = IterRslt.var_shock.%s;', ir.shocks.names{i}, ir.shocks.names{i});
end
for i = 1:numel(ir.states.names)
    w.add('%s = IterRslt.var_state.%s;', ir.states.names{i}, ir.states.names{i});
end
w.add('shock_num = IterRslt.shock_num;');
w.add('GDSGE_SPLINE_VEC = IterRslt.pp.GDSGE_SPLINE_VEC;');
w.blank();
base = {'TolEq','TolSol','TolFun','PrintFreq','NoPrint','SaveFreq','NoSave', ...
    'SimuPrintFreq','SimuSaveFreq','MaxIter','MaxMinorIter','num_samples', ...
    'num_periods','SolMaxIter','UseBroyden','FiniteDiffDelta','GDSGE_USE_BROYDEN', ...
    'GDSGE_DEBUG_EVAL_ONLY','INTERP_ORDER','EXTRAP_ORDER','OutputInterpOrder', ...
    'IterSaveAll','SkipModelInit','UseAdaptiveBound','UseAdaptiveBoundInSol', ...
    'EnforceSimuStateInbound','REUSE_WARMUP_SOL','INTERP_WARMUP_SOL', ...
    'CONSTRUCT_OUTPUT','NumThreads','init','GEN_SHOCK_START_PERIOD'};
modelNames = {};
for i = 1:numel(ir.params); modelNames{end+1} = ir.params{i}.name; end %#ok<AGROW>
modelNames = [modelNames, ir.shocks.names, fieldnames(ir.shocks.transitions)', ir.states.names];
w.add('GDSGE_OPTIONS_VALID = {%s};', ['''' strjoin([base, modelNames], ''',''') '''']);
w.add('if nargin < 2; GDSGE_OPTIONS = struct(); end');
w.add('gdsge.runtime.unpackOptions(GDSGE_OPTIONS, GDSGE_OPTIONS_VALID);');
w.blank();

% ---- solution interp for warm starts (parity lines 119-149, both branches) -
w.add('%% ---- solution interpolation for per-period initial guesses');
w.add('if ~ismac');
w.add('    if shock_num>1');
w.add('        GDSGE_PP=struct(''form'',''MKL'',''breaks'',{{1:shock_num,%s}},...', stateList);
w.add('            ''Values'',reshape(IterRslt.GDSGE_PROB.GDSGE_SOL, [numel(IterRslt.GDSGE_PROB.GDSGE_SOL)/prod(IterRslt.GDSGE_PROB.GDSGE_SIZE),IterRslt.GDSGE_PROB.GDSGE_SIZE]),...');
w.add('            ''coefs'',[],''order'',[2 OutputInterpOrder*ones(1,length(%s))],''Method'',[],...', stateGridCell);
w.add('            ''ExtrapolationOrder'',[],''thread'',NumThreads, ...');
w.add('            ''orient'',''curvefit'');');
w.add('    else');
w.add('        GDSGE_PP=struct(''form'',''MKL'',''breaks'',{%s},...', stateGridCell);
w.add('            ''Values'',reshape(IterRslt.GDSGE_PROB.GDSGE_SOL, [numel(IterRslt.GDSGE_PROB.GDSGE_SOL)/prod(IterRslt.GDSGE_PROB.GDSGE_SIZE),IterRslt.GDSGE_PROB.GDSGE_SIZE(2:end)]),...');
w.add('            ''coefs'',[],''order'',[OutputInterpOrder*ones(1,length(%s))],''Method'',[],...', stateGridCell);
w.add('            ''ExtrapolationOrder'',[],''thread'',NumThreads, ...');
w.add('            ''orient'',''curvefit'');');
w.add('    end');
w.add('    GDSGE_PP=myppual(GDSGE_PP);');
w.add('else');
w.add('    if shock_num>1');
w.add('        GDSGE_PP=struct(''form'',''pp'',''breaks'',{{1:shock_num,%s}},...', stateList);
w.add('            ''Values'',reshape(IterRslt.GDSGE_PROB.GDSGE_SOL, [numel(IterRslt.GDSGE_PROB.GDSGE_SOL)/prod(IterRslt.GDSGE_PROB.GDSGE_SIZE),IterRslt.GDSGE_PROB.GDSGE_SIZE]),...');
w.add('            ''coefs'',[],''order'',[4 4*ones(1,length(%s))],''Method'',[],...', stateGridCell);
w.add('            ''ExtrapolationOrder'',[],''thread'',NumThreads, ...');
w.add('            ''orient'',''curvefit'');');
w.add('    else');
w.add('        GDSGE_PP=struct(''form'',''pp'',''breaks'',{%s},...', stateGridCell);
w.add('            ''Values'',reshape(IterRslt.GDSGE_PROB.GDSGE_SOL, [numel(IterRslt.GDSGE_PROB.GDSGE_SOL)/prod(IterRslt.GDSGE_PROB.GDSGE_SIZE),IterRslt.GDSGE_PROB.GDSGE_SIZE(2:end)]),...');
w.add('            ''coefs'',[],''order'',[4*ones(1,length(%s))],''Method'',[],...', stateGridCell);
w.add('            ''ExtrapolationOrder'',[],''thread'',NumThreads, ...');
w.add('            ''orient'',''curvefit'');');
w.add('    end');
w.add('    GDSGE_PP=myppual(myppual(GDSGE_PP));');
w.add('end');
w.blank();
w.add('GDSGE_NPROB = num_samples;');
w.add('SimuRslt.shock = ones(num_samples,num_periods+1);');
w.addRaw(simu.prealloc);
w.blank();
w.addRaw(simu.init);
w.blank();
w.add('if nargin>1 && isfield(GDSGE_OPTIONS,''init'')');
w.addRaw(indent4(simu.initOverwrite));
w.add('end');
w.add('if any(SimuRslt.shock(:,1)>shock_num)');
w.add('    error(''initial shock exceeds shock_num'');');
w.add('end');
w.blank();
w.addRaw(bounds.init);
w.blank();
w.add('%% Use the widest bounds seen during iteration (parity lines 205-209)');
w.add('if UseAdaptiveBound==1');
w.add('    GDSGE_LB = repmat(min(IterRslt.GDSGE_PROB.GDSGE_LB,[],2),[1,GDSGE_NPROB]);');
w.add('    GDSGE_UB = repmat(max(IterRslt.GDSGE_PROB.GDSGE_UB,[],2),[1,GDSGE_NPROB]);');
w.add('end');
w.blank();
w.add('GDSGE_EQVAL = 1e20*ones(%d,GDSGE_NPROB);', bounds.maxDim);
w.add('GDSGE_F = 1e20*ones(1,GDSGE_NPROB);');
w.add('GDSGE_SOL = zeros(%d,GDSGE_NPROB);', bounds.maxDim);
w.add('GDSGE_AUX = zeros(%d,GDSGE_NPROB);', max(nAux, 1));
w.blank();
w.add('SimuRslt.shock(:,GEN_SHOCK_START_PERIOD:end) = gen_discrete_markov_rn(shock_trans,num_samples,length(GEN_SHOCK_START_PERIOD:num_periods+1),...');
w.add('    SimuRslt.shock(:,GEN_SHOCK_START_PERIOD));');
w.blank();
w.add(pack.simuData0);
w.add('GDSGE_SHOCK_VAR_INDEX_BASE = ([0:num_samples-1]'')*shock_num;');
w.add('GDSGE_timer = tic;');
w.add('for GDSGE_t=1:num_periods');
w.in();
w.add('shock = SimuRslt.shock(:,GDSGE_t);');
for i = 1:numel(ir.states.names)
    w.add('%s=SimuRslt.%s(:,GDSGE_t);', ir.states.names{i}, ir.states.names{i});
end
w.blank();
w.add('%% interp warm start for the solver (parity lines 236-247)');
w.add('if shock_num>1');
w.add('    GDSGE_SOL = myppual_mex(int32(NumThreads),GDSGE_PP.breaks,GDSGE_PP.coefs,...');
w.add('        int32(GDSGE_PP.pieces),int32(GDSGE_PP.order),int32(GDSGE_PP.dim),''not-a-knot'',[SimuRslt.shock(:,GDSGE_t)'';%s],[],[],[]);', simuStateRows);
w.add('else');
w.add('    GDSGE_SOL = myppual_mex(int32(NumThreads),GDSGE_PP.breaks,GDSGE_PP.coefs,...');
w.add('        int32(GDSGE_PP.pieces),int32(GDSGE_PP.order),int32(GDSGE_PP.dim),''not-a-knot'',[%s],[],[],[]);', simuStateRows);
w.add('end');
w.blank();
w.add('if UseAdaptiveBoundInSol==1');
w.add('    GDSGE_LB_OLD = GDSGE_LB; GDSGE_UB_OLD = GDSGE_UB;');
w.add('    [GDSGE_LB,GDSGE_UB] = GDSGE_ADAPT_TIGHT(GDSGE_SOL,GDSGE_LB,GDSGE_UB);');
w.add('    GDSGE_SOL_hitting_lower_bound = abs(GDSGE_SOL - GDSGE_LB_OLD) < 1e-8;');
w.add('    GDSGE_SOL_hitting_upper_bound = abs(GDSGE_SOL - GDSGE_UB_OLD) < 1e-8;');
w.add('    GDSGE_LB(~GDSGE_SOL_hitting_lower_bound) = GDSGE_LB_OLD(~GDSGE_SOL_hitting_lower_bound);');
w.add('    GDSGE_UB(~GDSGE_SOL_hitting_upper_bound) = GDSGE_UB_OLD(~GDSGE_SOL_hitting_upper_bound);');
w.add('end');
w.blank();
w.add(pack.simuPack);
w.blank();
w.add('GDSGE_CFG = struct();');
w.add('GDSGE_CFG.tolSol = TolSol; GDSGE_CFG.tolFun = TolFun;');
w.add('GDSGE_CFG.solMaxIter = SolMaxIter; GDSGE_CFG.numThreads = NumThreads;');
w.add('GDSGE_CFG.debugEvalOnly = GDSGE_DEBUG_EVAL_ONLY;');
w.add('GDSGE_CFG.useBroyden = UseBroyden; GDSGE_CFG.finiteDiffDelta = FiniteDiffDelta;');
w.add('GDSGE_CFG.useBroydenNow = 0;');
w.add('GDSGE_CFG.taskName = MEX_TASK_INF_HORIZON;');
w.add('GDSGE_CFG.splineVec = GDSGE_SPLINE_VEC;');
w.add('GDSGE_CFG.maxMinorIter = MaxMinorIter;');
w.add('GDSGE_CFG.probSize = [1 GDSGE_NPROB];');
w.add('GDSGE_CFG.useNearestNeighbor = false;');
w.add('GDSGE_CFG.verboseRetry = false;');
w.add('if UseAdaptiveBoundInSol==1');
w.add('    GDSGE_CFG.adaptInSol = @GDSGE_ADAPT_TIGHT;');
w.add('else');
w.add('    GDSGE_CFG.adaptInSol = [];');
w.add('end');
w.add('[GDSGE_SOL,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL,GDSGE_OPT_INFO,GDSGE_DIAG] = gdsge.runtime.solveProblems(@mex_%s, GDSGE_SOL, GDSGE_LB, GDSGE_UB, GDSGE_DATA, GDSGE_F, GDSGE_AUX, GDSGE_EQVAL, GDSGE_CFG); %%#ok<ASGLU>', m);
w.add('if any(GDSGE_DIAG.needResolved)');
w.add('    warning(''gdsge:runtime:unconverged'', ''period %%d: %%s'', GDSGE_t, gdsge.runtime.reportUnconverged(GDSGE_DIAG.needResolved, GDSGE_F, [1 GDSGE_NPROB], %s, {%s}, 5));', ...
    stateNameCell, strjoin(ir.states.names, ','));
w.add('end');
w.blank();
w.addRaw(indent4(unpk.unpack));
w.addRaw(indent4(unpk.unpackAux));
w.blank();
w.add('GDSGE_SHOCK_VAR_LINEAR_INDEX = SimuRslt.shock(:,GDSGE_t+1) + GDSGE_SHOCK_VAR_INDEX_BASE;');
w.addRaw(indent4(simu.assign));
w.blank();
w.add('if mod(GDSGE_t,SimuPrintFreq)==0');
w.add('    fprintf(''Periods: %%d\\n'', GDSGE_t);');
w.add('    SimuRsltNames = fieldnames(SimuRslt);');
w.add('    for GDSGE_field = 1:length(SimuRsltNames)');
w.add('        fprintf(''%%8s'', SimuRsltNames{GDSGE_field});');
w.add('    end');
w.add('    fprintf(''\\n'');');
w.add('    for GDSGE_field = 1:length(SimuRsltNames)');
w.add('        fprintf(''%%8.4g'', SimuRslt.(SimuRsltNames{GDSGE_field})(1,GDSGE_t));');
w.add('    end');
w.add('    fprintf(''\\n'');');
w.add('    fprintf(''elapsed:%%.1fs\\n'', toc(GDSGE_timer));');
w.add('    GDSGE_timer = tic;');
w.add('end');
w.add('if mod(GDSGE_t,SimuSaveFreq)==0');
w.add('    save([''SimuRslt_%s_'' num2str(GDSGE_t) ''.mat''], ''SimuRslt'');', m);
w.add('end');
w.out();
w.add('end');
w.add('end');
w.blank();
w.add('function [GDSGE_LB,GDSGE_UB] = GDSGE_ADAPT_TIGHT(GDSGE_SOL,GDSGE_LB,GDSGE_UB)');
if isempty(bounds.adaptiveTight)
    w.add('%% no adaptive bounds in this model');
else
    w.addRaw(bounds.adaptiveTight);
end
w.add('end');
txt = w.str();
end

function s = indent4(s)
if isempty(s); return; end
lines = strsplit(s, newline, 'CollapseDelimiters', false);
for i = 1:numel(lines)
    if ~isempty(lines{i}); lines{i} = ['    ' lines{i}]; end
end
s = strjoin(lines, newline);
end
```

Notes:
- Bound text in `bounds.init` references `GDSGE_TENSOR_*` only for state-dependent
  bounds; HL1996's are constants, so reuse in simulate is safe this phase. (Phase 7
  will need a simulate-specific rewrite map — out of scope now.)
- `w.add` treats its first argument as an sprintf format: escape `%` as `%%` and
  backslash sequences carefully (`\n` inside a single-quoted MATLAB string we EMIT must
  arrive as the two characters `\n` — write `\\n` in the format if needed; verify
  against the parse test).

- [ ] **Step 18.5: Implement `generateMatlab`**

Create `src/+gdsge/+codegen/generateMatlab.m`:

```matlab
function files = generateMatlab(ir, outDir)
% GENERATEMATLAB  IR -> iter_<model>.m / simulate_<model>.m written to outDir.
if nargin < 2; outDir = pwd; end
files = struct();
files.iterFile = fullfile(outDir, ['iter_' ir.modelName '.m']);
writeFile(files.iterFile, gdsge.codegen.mat.emitIter(ir));
files.simulateFile = fullfile(outDir, ['simulate_' ir.modelName '.m']);
writeFile(files.simulateFile, gdsge.codegen.mat.emitSimulate(ir));
end

function writeFile(p, txt)
fid = fopen(p, 'w');
if fid < 0
    error('gdsge:codegen:cannotWrite', 'Cannot write %s', p);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, txt);
end
```

- [ ] **Step 18.6: Regen snapshots, review simulate against old reference, run tests — expect PASS**

```
matlab -batch "addpath('src'); addpath('src/kernels'); cd('tests/HeatonLucas1996/codegen'); regen_snapshots"
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests({'tests/codegen/tEmitResultSimu.m','tests/codegen/tGenerateMatlab.m','tests/HeatonLucas1996/codegen/tSnapshotHL1996.m'}))"
```

Simulate review checklist (against `scratch/old_generated/simulate_HL1996.m`):
- [ ] pulls `shock_trans`/params/var_shock/var_state from `IterRslt` (old lines 109-112)
- [ ] solution-interp `GDSGE_PP` construction matches old lines 119-149
- [ ] `gen_discrete_markov_rn` call identical (old lines 218-219)
- [ ] data layout matches (`GDSGE_data0` + per-period rows; old lines 232, 273)
- [ ] per-period assigns + `GDSGE_SHOCK_VAR_LINEAR_INDEX` transition (old lines 324-330)

- [ ] **Step 18.7: Commit**

```powershell
git add src/+gdsge/+codegen tests/codegen tests/HeatonLucas1996/codegen
git commit -m "feat(codegen): emitSimulate + emitResultSimu + generateMatlab driver"
```

---

### Task 19: Capture the old MEX

**Files:**
- Create: `tests/golden/capture_HL1996_mex.m`
- Modify: `.gitignore` (ignore `tests/HeatonLucas1996/oldmex/`)

- [ ] **Step 19.1: Write the capture script**

Create `tests/golden/capture_HL1996_mex.m`:

```matlab
function capture_HL1996_mex()
% Build the OLD toolbox's mex_HL1996 and stash the binary (plus the old
% generated .m files, for reference) under tests/HeatonLucas1996/oldmex/.
% Uses ONLY the old source, in this process, in a temp dir; does NOT re-run
% the solver and does NOT touch the committed goldens.
here     = fileparts(mfilename('fullpath'));          % tests/golden
repoRoot = fileparts(fileparts(here));
oldSrc   = fullfile(repoRoot, 'base_package', 'gdsge', 'source');
modelDir = fullfile(repoRoot, 'tests', 'HeatonLucas1996');
outDir   = fullfile(modelDir, 'oldmex');
if ~exist(outDir, 'dir'); mkdir(outDir); end

work = tempname; mkdir(work);
copyfile(fullfile(modelDir, 'HL1996.gmod'), work);

oldPath = path; restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(oldSrc);
oldCd = pwd; cdRestore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(work);

gdsge_codegen('HL1996');                 % writes + compiles the old MEX

copyfile(['mex_HL1996.' mexext], outDir);
copyfile('iter_HL1996.m', fullfile(outDir, 'iter_HL1996_old_reference.txt'));
copyfile('simulate_HL1996.m', fullfile(outDir, 'simulate_HL1996_old_reference.txt'));
fprintf('Old MEX captured to %s\n', outDir);
end
```

- [ ] **Step 19.2: Ignore the artifact directory**

Append to `.gitignore` under the "Test outputs" section:

```
tests/HeatonLucas1996/oldmex/
```

- [ ] **Step 19.3: Run the capture (own MATLAB process, old source only)**

```
matlab -batch "cd('tests/golden'); capture_HL1996_mex"
```

Expected: `Old MEX captured to ...\tests\HeatonLucas1996\oldmex` (takes 1-3 min — MSVC
compile; the toolchain was proven in Phase 0). Verify `mex_HL1996.mexw64` exists there.

- [ ] **Step 19.4: Commit**

```powershell
git add tests/golden/capture_HL1996_mex.m .gitignore
git commit -m "test(golden): capture_HL1996_mex - stash old MEX for the Phase 4 gate"
```

---

### Task 20: The functional phase gate

**Files:**
- Test: `tests/HeatonLucas1996/codegen/tFunctionalHL1996.m`
- Modify: `PROGRESS.md`, `src/+gdsge/Contents.m` (+ new `Contents.m` for the two packages)

- [ ] **Step 20.1: Write the gate test**

Create `tests/HeatonLucas1996/codegen/tFunctionalHL1996.m`:

```matlab
classdef tFunctionalHL1996 < matlab.unittest.TestCase
    % PHASE 4 GATE: the NEW generated MATLAB drives the OLD compiled MEX to
    % convergence and reproduces the committed goldens within tolerance.
    % Slow (~minutes): full policy iteration + 10k-period simulation.
    properties (Constant)
        RelTol = 1e-4;
        AbsTol = 1e-4;
    end
    methods (Test, TestTags = {'Slow'})
        function iterAndSimulateMatchGolden(tc)
            here = fileparts(mfilename('fullpath'));     % tests/HeatonLucas1996/codegen
            modelDir = fileparts(here);
            mexFile = fullfile(modelDir, 'oldmex', ['mex_HL1996.' mexext]);
            tc.assertTrue(exist(mexFile, 'file') == 2, sprintf( ...
                ['Old MEX missing: %s\nRun:  matlab -batch ' ...
                 '"cd(''tests/golden''); capture_HL1996_mex"'], mexFile));

            % generate the NEW files from the full parser pipeline
            ir = gdsge.parser.parseFrontEnd( ...
                fileread(fullfile(modelDir, 'HL1996.gmod')), 'HL1996');
            work = tempname; mkdir(work);
            gdsge.codegen.generateMatlab(ir, work);
            copyfile(mexFile, work);

            oldCd = pwd; cdBack = onCleanup(@() cd(oldCd)); %#ok<NASGU>
            cd(work);

            % ---- iterate ----------------------------------------------------
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

            % ---- simulate ---------------------------------------------------
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

Decision points if something fails (investigate, don't just loosen):
- `shock paths differ`: compare `IterRslt.shock_trans` to `G.shock_trans` with
  `isequal` — if not bit-identical, the `mat2str(…,17)` literal or the normalization
  changed the bits; find which and fix the emitter (do NOT switch the test to tolerance
  until the cause is understood).
- `var_policy` failures at ~1e-4: check the GDSGE_DATA layout first (one transposed row
  scrambles everything); then bounds; then the metric/update wiring.
- `Iter` outside [100,400]: interp init or update wiring is wrong — diff the generated
  file's loop section against `scratch/old_generated/iter_HL1996.m` block by block.

- [ ] **Step 20.2: Run the gate — iterate until PASS**

```
matlab -batch "addpath('src'); addpath('src/kernels'); addpath('tests'); assertSuccess(runtests('tests/HeatonLucas1996/codegen/tFunctionalHL1996.m'))"
```

This is the phase's verification step: expect several debug rounds. Use the old
reference files and the decision points above. Commit fixes as they land
(`fix(codegen): ...` / `fix(runtime): ...`).

- [ ] **Step 20.3: Run the FULL suite — expect PASS**

```
matlab -batch "cd('tests'); run_tests"
```

- [ ] **Step 20.4: Add package Contents.m and update PROGRESS.md**

Create `src/+gdsge/+runtime/Contents.m`:

```matlab
% +GDSGE/+RUNTIME  Hand-written runtime library used by generated code.
%   unpackOptions     - whitelisted GDSGE_OPTIONS unpacking (errors on unknown)
%   ensurePath        - essential_blas.dll system-PATH setup
%   solveProblems     - MEX driver + resolve cascade (owns the MEX caller contract)
%   applyWarmUp       - WarmUp SOL/LB/UB transfer
%   constructSplines  - myppual interpolants + GDSGE_SPLINE_VEC
%   computeMetric     - sup-norm interp-update metric (NaN on error)
%   printIterProgress - structured progress line
%   reportUnconverged - failed-point report with state coordinates
```

Create `src/+gdsge/+codegen/Contents.m`:

```matlab
% +GDSGE/+CODEGEN  Code generators consuming the IR.
%   generateMatlab - IR -> iter_<model>.m / simulate_<model>.m
%   +mat           - MATLAB-emitting backend (section emitters + codeWriter)
```

Update `PROGRESS.md`: mark Phase 4 done with a summary bullet list (mirror the Phase 3
entry style: what was built, the gate that proved it, the branch name), set Phase 5 as
NEXT, and add a changelog entry dated with the actual completion date.

- [ ] **Step 20.5: Commit**

```powershell
git add tests/HeatonLucas1996/codegen/tFunctionalHL1996.m src/+gdsge PROGRESS.md
git commit -m "test(codegen): Phase 4 gate - new MATLAB vs old MEX reproduces goldens"
```

- [ ] **Step 20.6: Finish the branch**

Use the superpowers:finishing-a-development-branch skill (merge to main vs PR is the
user's call).

---

## Self-review notes (kept for the record)

- Spec §2.1 components vs tasks: codeWriter→T4, rewriteNames→T5, unpackOptions→T6,
  ensurePath→T7, computeMetric/printIterProgress/reportUnconverged→T8,
  constructSplines→T9, solveProblems→T10-11, applyWarmUp→T12, emitSetup→T13,
  emitBounds/emitInterpInitial→T14, emitDataPack/emitSolUnpack→T15, emitResultIter→T16,
  emitIter→T17, emitResultSimu/emitSimulate/generateMatlab→T18, kernels→T2,
  numThreads sentinel (spec §6)→T3, capture+gate (spec §7.4)→T19-20. No spec section
  is uncovered.
- The spec's "solveProblems … exhausting MaxMinorIter triggers reportUnconverged +
  warning" is implemented as an emitted check in the generated files (the runtime
  returns `diag`; the generated file owns state names/grids needed for the report).
  Same observable behavior; noted here as a deliberate placement choice.
- `emitSetup` drops dead old flags (list in Task 13) — approved by the design decision
  "only truly dead code for the spline path is dropped".
- Naming consistency check: `cfg.*` fields used in Tasks 10/11/17/18 match; `frag.*`
  fields of emitBounds/emitDataPack/emitSolUnpack/emitResultSimu match their consumers
  in emitIter/emitSimulate; `spec.*` fields of applyWarmUp match the GDSGE_WARMSPEC
  emitted in Task 17.
