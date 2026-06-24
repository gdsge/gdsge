# Phase 2 — Parser Front-End → Partial IR — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a `gdsge.parser` pipeline that turns `HL1996.gmod` text into a schema-valid **partial IR** (everything except the `model;…end;` body), proven against the corrected Phase-1 reference IR.

**Architecture:** Six small, independently-tested stages under `src/+gdsge/+parser/` — `preprocess` → `splitBlocks` → `parseDeclarations` (which uses `parseVarDecls`, `splitStatements`, `defaultSetupCode`, `evalSetup`) + `parseSimulate` + `resolveOptions` → `assemblePartialIR`, chained by `parseFrontEnd`. The IR gains one additive `params` section so parameter values have a home. The declaration region is eval'd in an isolated workspace (mirroring the old `gdsge_precode.m`); interp-update assignments and state grids are captured as text, never eval'd.

**Tech Stack:** MATLAB R2025b (`eval`, `regexp`, `matlab.unittest`), the Phase-0 headless harness (`tests/run_tests.m`, `tests/run.ps1`), the Phase-1 `gdsge.ir` package (schema/validate/canonicalize/encode/isequalIR).

**Spec:** `docs/superpowers/specs/2026-06-12-phase2-parser-frontend-design.md`.

---

## File structure

Created in `src/+gdsge/+parser/`:

- `preprocess.m` — clean gmod text (comments, continuations, deprecated rewrite); error on `#` macros.
- `splitBlocks.m` — separate declaration lines from named blocks (`model`, `simulate`, hooks).
- `splitStatements.m` — bracket-aware split of text into logical statements (multi-line matrices stay whole).
- `defaultSetupCode.m` — parser-owned default flag values (the IR-relevant subset of the old `default_mod.nmod` + `params_template.m`).
- `evalSetup.m` — eval a setup script in an isolated workspace; capture variables into a struct; located errors.
- `parseVarDecls.m` — pure structural parse of the declaration region (names, lengths, bounds, interp exprs, grid text, residual setup statements).
- `parseDeclarations.m` — orchestrates `parseVarDecls` + `evalSetup`, builds the `params`/`shocks`/`states`/`variables`(+slots)/`bounds`/`interp` IR sections.
- `parseSimulate.m` — the `simulate;…end;` block → IR `simulate` section (opaque text).
- `resolveOptions.m` — eval'd flag workspace → curated IR `options`.
- `assemblePartialIR.m` — combine sections into a partial IR, run semantic gate (bounds completeness) + `gdsge.ir.validate`.
- `parseFrontEnd.m` — top-level chain: gmod text → validated partial IR.
- `Contents.m` — package help.

Modified:

- `src/+gdsge/+ir/schema.m` — add the `params` root section.
- `tests/ir/tSchema.m` — add `'params'` to the expected root-field list.
- `tests/+gdsgefix/minimalIR.m` — add a `params` field.
- `tests/HeatonLucas1996/ir/buildHL1996IR.m` — add `params`; correct `options`.
- `docs/ir-schema.md` — regenerated from the schema.
- `tests/HeatonLucas1996/ir/HL1996.gdsge.json` — regenerated golden.

Tests (plain folders so `TestSuite.fromFolder` finds them; `src/` is already on the test path):

- `tests/ir/tParams.m`
- `tests/parser/tPreprocess.m`, `tSplitBlocks.m`, `tSplitStatements.m`, `tEvalSetup.m`,
  `tParseVarDecls.m`, `tParseDeclarations.m`, `tParseSimulate.m`, `tResolveOptions.m`
- `tests/HeatonLucas1996/parser/tFrontEndHL1996.m`

**No harness change needed:** `tests/run_tests.m` already adds `src/` to the path and runs
`TestSuite.fromFolder(..., 'IncludingSubfolders', true)`, so the new `tests/parser/` and
`tests/HeatonLucas1996/parser/` folders are auto-discovered.

**Commands** (run from repo root `D:\refactor_gdsge`):
- Full suite: `pwsh -File tests/run.ps1` (exit 0 = pass).
- Single test file: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tPreprocess.m')); disp(table(r)); assert(all(~[r.Failed]))"`
- MATLAB binary if needed explicitly: `C:\Program Files\MATLAB\R2025b\bin\matlab.exe`.

---

## Task 1: Add the `params` IR section + correct the HL1996 fixture/goldens

**Files:**
- Modify: `src/+gdsge/+ir/schema.m`, `tests/ir/tSchema.m`, `tests/+gdsgefix/minimalIR.m`, `tests/HeatonLucas1996/ir/buildHL1996IR.m`
- Regenerate: `docs/ir-schema.md`, `tests/HeatonLucas1996/ir/HL1996.gdsge.json`
- Test: `tests/ir/tParams.m`

- [ ] **Step 1: Write the failing test**

`tests/ir/tParams.m`:
```matlab
classdef tParams < matlab.unittest.TestCase
    methods (Test)
        function paramsPresentValidates(tc)
            ir = gdsgefix.minimalIR();
            r = gdsge.ir.validate(ir);
            tc.verifyTrue(r.pass, strjoin(r.errors, ' | '));
        end
        function missingParamsFails(tc)
            ir = gdsgefix.minimalIR();
            ir = rmfield(ir, 'params');
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
            tc.verifyTrue(any(contains(r.errors, 'params')));
        end
        function paramScalarRoundtrips(tc)
            ir = gdsgefix.minimalIR();
            tc.verifyEqual(ir.params{1}.name, 'beta');
            ir2 = gdsge.ir.roundtrip(ir);
            tc.verifyEqual(ir2.params{1}.name, 'beta');
            tc.verifyEqual(ir2.params{1}.value, 0.95, 'AbsTol', 1e-12);
        end
        function paramArrayRoundtrips(tc)
            ir = gdsgefix.minimalIR();
            ir.params = { struct('name','rho','value',[0.9 0.8 0.7]) };
            ir2 = gdsge.ir.roundtrip(ir);
            tc.verifyEqual(ir2.params{1}.value, [0.9 0.8 0.7], 'AbsTol', 1e-12);
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('ir','tParams.m')); disp(table(r))"`
Expected: FAIL — `minimalIR` has no `params` field; schema has no `params` section.

- [ ] **Step 3: Add `params` to the schema**

In `src/+gdsge/+ir/schema.m`, after the `transItem = ...` line (the reusable item specs block), add a `paramItem` spec:
```matlab
paramItem = fStruct(structOf('name', fText(), 'value', fMatrix()));
```
Then in `s.root = fStruct(structOf( ...`, insert `params` right after `modelName`:
```matlab
s.root = fStruct(structOf( ...
    'irVersion', fText(), ...
    'modelName', fText(), ...
    'params',    fList(paramItem), ...
    'options',   options, ...
    'shocks',    shocks, ...
    'states',    states, ...
    'variables', variables, ...
    'bounds',    fList(boundItem), ...
    'interp',    fList(interpItem), ...
    'model',     model, ...
    'modelInit', opt(modelInit), ...
    'simulate',  simulate, ...
    'hooks',     hooks));
```

- [ ] **Step 4: Update `tSchema` expected root fields**

In `tests/ir/tSchema.m`, the `rootIsStructWithAllSections` test hard-codes the section list. Add `'params'`:
```matlab
            expected = {'irVersion','modelName','params','options','shocks','states', ...
                        'variables','bounds','interp','model','modelInit', ...
                        'simulate','hooks'};
```

- [ ] **Step 5: Add `params` to the `minimalIR` fixture**

In `tests/+gdsgefix/minimalIR.m`, after the `ir.modelName = 'mini';` line, add:
```matlab
ir.params = { struct('name','beta','value',0.95) };
```

- [ ] **Step 6: Add `params` + correct `options` in the HL1996 reference fixture**

In `tests/HeatonLucas1996/ir/buildHL1996IR.m`, after the `ir.modelName = 'HL1996';` line, add:
```matlab
ir.params = { struct('name','beta','value',0.95), ...
              struct('name','gamma','value',1.5), ...
              struct('name','Kb','value',-0.05) };
```
Then replace the `ir.options = struct(...)` block with the corrected defaults (real toolbox
values; `numThreads` stays pinned at 8 for a deterministic JSON golden):
```matlab
ir.options = struct('interpMethod','spline','interpOrder',4,'extrapOrder',2, ...
    'tolEq',1e-6,'numThreads',8,'simuResolve',1,'simuInterp',0, ... % numThreads pinned for golden determinism
    'printFreq',10,'saveFreq',10);
```

- [ ] **Step 7: Regenerate the doc and JSON goldens**

Run (from repo root):
```
matlab -batch "addpath('src'); gdsge.ir.gendoc(fullfile('docs','ir-schema.md'));"
```
```
matlab -batch "addpath('src'); addpath(fullfile('tests','HeatonLucas1996','ir')); ir=buildHL1996IR(); fid=fopen(fullfile('tests','HeatonLucas1996','ir','HL1996.gdsge.json'),'w'); fwrite(fid, gdsge.ir.encode(ir)); fclose(fid);"
```

- [ ] **Step 8: Run the full suite to verify green**

Run: `pwsh -File tests/run.ps1`
Expected: PASS — `tParams` (4) green; `tSchema`, `tIrHL1996` (encoding/round-trip), `tGendoc` (no-drift), and all Phase-0/1 tests still green.

- [ ] **Step 9: Commit**

```bash
git add src/+gdsge/+ir/schema.m tests/ir/tSchema.m tests/ir/tParams.m tests/+gdsgefix/minimalIR.m tests/HeatonLucas1996/ir/buildHL1996IR.m docs/ir-schema.md tests/HeatonLucas1996/ir/HL1996.gdsge.json
git commit -m "feat(ir): add params section; correct HL1996 options to real toolbox defaults"
```

---

## Task 2: `preprocess` — clean text, guard macros

**Files:**
- Create: `src/+gdsge/+parser/preprocess.m`
- Test: `tests/parser/tPreprocess.m`

- [ ] **Step 1: Write the failing test**

`tests/parser/tPreprocess.m`:
```matlab
classdef tPreprocess < matlab.unittest.TestCase
    methods (Test)
        function stripsLineComments(tc)
            out = gdsge.parser.preprocess(sprintf('x = 1; %% hello\ny = 2;'));
            tc.verifyFalse(contains(out, 'hello'));
            tc.verifyTrue(contains(out, 'x = 1;'));
        end
        function joinsContinuations(tc)
            out = gdsge.parser.preprocess(sprintf('x = 1 + ...\n2;'));
            oneLine = regexprep(out, '\s+', ' ');
            tc.verifyTrue(contains(oneLine, '1 + 2'));
        end
        function rewritesDeprecatedKeyword(tc)
            out = gdsge.parser.preprocess('GNDSGE_X = 1;');
            tc.verifyTrue(contains(out, 'GDSGE_X'));
            tc.verifyFalse(contains(out, 'GNDSGE'));
        end
        function macroDirectiveErrors(tc)
            tc.verifyError(@() gdsge.parser.preprocess(sprintf('#define N 3\nx=1;')), ...
                'gdsge:parser:macroUnsupported');
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tPreprocess.m')); disp(table(r))"`
Expected: FAIL — `gdsge.parser.preprocess` undefined.

- [ ] **Step 3: Implement `preprocess`**

`src/+gdsge/+parser/preprocess.m`:
```matlab
function clean = preprocess(rawText)
% PREPROCESS  Clean gmod text into logical lines: strip % comments, join MATLAB
%   line-continuations (...), rewrite deprecated keywords. Raises
%   gdsge:parser:macroUnsupported on any #-macro (full engine is Phase 7).
%   NOTE: the comment rule cuts at the first % on a line (sufficient for the
%   core models); string-aware stripping is a Phase-7 hardening point.
rawLines = regexp(rawText, '\r\n|\r|\n', 'split');

% 1. strip % comments
for i = 1:numel(rawLines)
    rawLines{i} = stripComment(rawLines{i});
end

% 2. join line-continuations
logical = {};
acc = '';
for i = 1:numel(rawLines)
    t = strtrim(rawLines{i});
    if endsWith(t, '...')
        acc = [acc, ' ', t(1:end-3)]; %#ok<AGROW>
    else
        joined = strtrim([acc, ' ', rawLines{i}]);
        logical{end+1} = joined; %#ok<AGROW>
        acc = '';
    end
end
if ~isempty(strtrim(acc)); logical{end+1} = strtrim(acc); end

% 3. macro guard
for i = 1:numel(logical)
    if ~isempty(regexp(logical{i}, '^\s*#', 'once'))
        error('gdsge:parser:macroUnsupported', ...
            'Macro directive not supported until Phase 7: "%s"', strtrim(logical{i}));
    end
end

% 4. deprecated-keyword rewrite (word-boundary)
clean = strjoin(logical, sprintf('\n'));
clean = regexprep(clean, '\<GNDSGE', 'GDSGE');
end

function s = stripComment(line)
idx = strfind(line, '%');
if isempty(idx); s = line; else; s = line(1:idx(1)-1); end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tPreprocess.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 4 tests.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/preprocess.m tests/parser/tPreprocess.m
git commit -m "feat(parser): preprocess stage (comments, continuations, deprecated rewrite, macro guard)"
```

---

## Task 3: `splitBlocks` — declaration lines vs named blocks

**Files:**
- Create: `src/+gdsge/+parser/splitBlocks.m`
- Test: `tests/parser/tSplitBlocks.m`

- [ ] **Step 1: Write the failing test**

`tests/parser/tSplitBlocks.m`:
```matlab
classdef tSplitBlocks < matlab.unittest.TestCase
    methods (Test)
        function separatesModelSimulateAndDecls(tc)
            txt = sprintf(['parameters a;\na = 1;\n', ...
                'model;\n  x = a;\n  equations;\n    x;\n  end;\nend;\n', ...
                'simulate;\n  num_periods = 5;\nend;']);
            out = gdsge.parser.splitBlocks(txt);
            tc.verifyTrue(contains(out.declText, 'parameters a'));
            tc.verifyTrue(contains(out.declText, 'a = 1'));
            tc.verifyFalse(contains(out.declText, 'x = a'));
            tc.verifyTrue(isfield(out.blocks, 'model'));
            tc.verifyTrue(contains(out.blocks.model, 'x = a'));
            tc.verifyTrue(contains(out.blocks.model, 'equations'));   % nested sub-block retained
            tc.verifyTrue(isfield(out.blocks, 'simulate'));
            tc.verifyTrue(contains(out.blocks.simulate, 'num_periods'));
        end
        function unterminatedBlockErrors(tc)
            tc.verifyError(@() gdsge.parser.splitBlocks(sprintf('model;\n  x = 1;')), ...
                'gdsge:parser:unterminatedBlock');
        end
        function strayEndErrors(tc)
            tc.verifyError(@() gdsge.parser.splitBlocks(sprintf('x = 1;\nend;')), ...
                'gdsge:parser:unterminatedBlock');
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tSplitBlocks.m')); disp(table(r))"`
Expected: FAIL — `gdsge.parser.splitBlocks` undefined.

- [ ] **Step 3: Implement `splitBlocks`**

`src/+gdsge/+parser/splitBlocks.m`:
```matlab
function out = splitBlocks(cleanText)
% SPLITBLOCKS  Partition preprocessed gmod into declaration lines and named
%   blocks. Returns struct with:
%     .declText           the non-block lines (joined with \n)
%     .blocks.<name>      raw text of each top-level block (opener/closer removed)
%   Block openers (depth-tracked, so nested equations; ... end; is retained
%   inside model): model simulate model_init pre_model pre_iter post_iter
%   pre_jac_code post_jac_code equations. Closer: end; .
openers = {'model','simulate','model_init','pre_model','pre_iter', ...
           'post_iter','pre_jac_code','post_jac_code','equations'};
lines = regexp(cleanText, '\n', 'split');
declLines = {};
blocks = struct();
stack = {};
curTop = '';
curBody = {};
for i = 1:numel(lines)
    raw = lines{i};
    t = strtrim(raw);
    opener = matchOpener(t, openers);
    isEnd = strcmp(t, 'end;') || strcmp(t, 'end');
    if isempty(stack)
        if ~isempty(opener)
            stack{end+1} = opener; curTop = opener; curBody = {}; %#ok<AGROW>
        elseif isEnd
            error('gdsge:parser:unterminatedBlock', 'Unexpected "end;" at line %d', i);
        else
            declLines{end+1} = raw; %#ok<AGROW>
        end
    else
        if ~isempty(opener)
            stack{end+1} = opener; curBody{end+1} = raw; %#ok<AGROW>
        elseif isEnd
            stack(end) = [];
            if isempty(stack)
                blocks.(curTop) = strjoin(curBody, sprintf('\n'));
                curTop = ''; curBody = {};
            else
                curBody{end+1} = raw; %#ok<AGROW>
            end
        else
            curBody{end+1} = raw; %#ok<AGROW>
        end
    end
end
if ~isempty(stack)
    error('gdsge:parser:unterminatedBlock', 'Block "%s" not closed by "end;"', stack{end});
end
out = struct('declText', strjoin(declLines, sprintf('\n')), 'blocks', blocks);
end

function name = matchOpener(t, openers)
name = '';
tt = regexprep(t, ';\s*$', '');
if ismember(tt, openers)
    name = tt;
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tSplitBlocks.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 3 tests.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/splitBlocks.m tests/parser/tSplitBlocks.m
git commit -m "feat(parser): splitBlocks (declaration lines vs named blocks, depth-tracked)"
```

---

## Task 4: `splitStatements` — bracket-aware statement split

**Files:**
- Create: `src/+gdsge/+parser/splitStatements.m`
- Test: `tests/parser/tSplitStatements.m`

- [ ] **Step 1: Write the failing test**

`tests/parser/tSplitStatements.m`:
```matlab
classdef tSplitStatements < matlab.unittest.TestCase
    methods (Test)
        function splitsOnSemicolonAndNewline(tc)
            s = gdsge.parser.splitStatements(sprintf('x=1; y=2\nz=3;'));
            tc.verifyEqual(numel(s), 3);
            tc.verifyEqual(s{1}, 'x=1');
            tc.verifyEqual(s{3}, 'z=3');
        end
        function keepsMultilineMatrixWhole(tc)
            txt = sprintf('a = 1;\nM = [1 2\n3 4];\nb = 2;');
            s = gdsge.parser.splitStatements(txt);
            tc.verifyEqual(numel(s), 3);
            tc.verifyTrue(contains(s{2}, '1 2'));
            tc.verifyTrue(contains(s{2}, sprintf('\n')));   % row separator preserved
        end
        function matrixStatementEvalsToRightShape(tc)
            txt = sprintf('M = [1 2\n3 4];');
            s = gdsge.parser.splitStatements(txt);
            M = eval([s{1} ';']);
            tc.verifyEqual(size(M), [2 2]);
        end
        function dropsEmptyStatements(tc)
            s = gdsge.parser.splitStatements(sprintf(';;\n\nx=1;;'));
            tc.verifyEqual(numel(s), 1);
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tSplitStatements.m')); disp(table(r))"`
Expected: FAIL — `gdsge.parser.splitStatements` undefined.

- [ ] **Step 3: Implement `splitStatements`**

`src/+gdsge/+parser/splitStatements.m`:
```matlab
function stmts = splitStatements(text)
% SPLITSTATEMENTS  Split MATLAB-ish text into logical statements. Boundaries are
%   ';' or newline at bracket-depth 0; inside ( ) [ ] { } both are ignored, and
%   newlines inside brackets are PRESERVED (they are matrix row separators).
%   Empty statements are dropped. Returns a row cellstr of trimmed statements.
stmts = {};
depth = 0;
cur = '';
nl = sprintf('\n');
cr = sprintf('\r');
for i = 1:numel(text)
    ch = text(i);
    switch ch
        case {'[','(','{'}
            depth = depth + 1; cur(end+1) = ch; %#ok<AGROW>
        case {']',')','}'}
            depth = max(0, depth - 1); cur(end+1) = ch; %#ok<AGROW>
        case ';'
            if depth == 0
                [stmts, cur] = flush(stmts, cur);
            else
                cur(end+1) = ch; %#ok<AGROW>
            end
        case {nl, cr}
            if depth == 0
                [stmts, cur] = flush(stmts, cur);
            else
                cur(end+1) = nl; %#ok<AGROW>
            end
        otherwise
            cur(end+1) = ch; %#ok<AGROW>
    end
end
[stmts, ~] = flush(stmts, cur);
end

function [stmts, cur] = flush(stmts, cur)
t = strtrim(cur);
if ~isempty(t); stmts{end+1} = t; end
cur = '';
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tSplitStatements.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 4 tests.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/splitStatements.m tests/parser/tSplitStatements.m
git commit -m "feat(parser): splitStatements (bracket-aware; multi-line matrices stay whole)"
```

---

## Task 5: `defaultSetupCode` + `evalSetup` — the eval sandbox

**Files:**
- Create: `src/+gdsge/+parser/defaultSetupCode.m`, `src/+gdsge/+parser/evalSetup.m`
- Test: `tests/parser/tEvalSetup.m`

- [ ] **Step 1: Write the failing test**

`tests/parser/tEvalSetup.m`:
```matlab
classdef tEvalSetup < matlab.unittest.TestCase
    methods (Test)
        function capturesVariables(tc)
            script = sprintf('shock_num = 8;\nM = [1 2\n3 4];\nrowsum = sum(M,2);');
            ws = gdsge.parser.evalSetup(script);
            tc.verifyEqual(ws.shock_num, 8);
            tc.verifyEqual(size(ws.M), [2 2]);
            tc.verifyEqual(ws.rowsum, [3;7]);
        end
        function normalizesTransitionRowsToOne(tc)
            script = sprintf(['shock_num = 2;\n', ...
                'shock_trans = [3 1\n1 3];\n', ...
                'shock_trans = shock_trans ./ repmat(sum(shock_trans,2),[1,shock_num]);']);
            ws = gdsge.parser.evalSetup(script);
            tc.verifyEqual(sum(ws.shock_trans, 2), [1;1], 'AbsTol', 1e-12);
        end
        function badLineRaisesLocatedError(tc)
            tc.verifyError(@() gdsge.parser.evalSetup(sprintf('x = 1;\ny = nosuchfn_xyz(3);')), ...
                'gdsge:parser:setupEvalFailed');
        end
        function defaultsResolve(tc)
            ws = gdsge.parser.evalSetup(gdsge.parser.defaultSetupCode());
            tc.verifyEqual(ws.INTERP_ORDER, 4);
            tc.verifyEqual(ws.USE_SPLINE, 1);
            tc.verifyEqual(ws.SIMU_RESOLVE, 1);
            tc.verifyEqual(ws.SIMU_INTERP, 0);
            tc.verifyEqual(ws.TolEq, 1e-6);
            tc.verifyEqual(ws.PrintFreq, 10);
            tc.verifyTrue(ws.NumThreads >= 1);
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tEvalSetup.m')); disp(table(r))"`
Expected: FAIL — `gdsge.parser.evalSetup` / `defaultSetupCode` undefined.

- [ ] **Step 3: Implement `defaultSetupCode`**

`src/+gdsge/+parser/defaultSetupCode.m`:
```matlab
function code = defaultSetupCode()
% DEFAULTSETUPCODE  Parser-owned default flag values — the IR-relevant subset of
%   the old default_mod.nmod + params_template.m. Prepended to a model's setup
%   script before eval, so options resolve and any gmod assignment overrides.
code = strjoin({ ...
    'INTERP_ORDER = 4;', ...
    'EXTRAP_ORDER = 2;', ...
    'USE_SPLINE = 1;', ...
    'USE_ASG = 0;', ...
    'USE_PCHIP = 0;', ...
    'SIMU_INTERP = 0;', ...
    'SIMU_RESOLVE = 1;', ...
    'AsgMaxLevel = 10;', ...
    'AsgThreshold = 1e-2;', ...
    'TolEq = 1e-6;', ...
    'PrintFreq = 10;', ...
    'SaveFreq = 10;', ...
    'SimuPrintFreq = 1000;', ...
    'SimuSaveFreq = Inf;', ...
    'NumThreads = feature(''numcores'');', ...
    'shock_num = 1;', ...
    'shock_trans = 1;' ...
    }, sprintf('\n'));
end
```

- [ ] **Step 4: Implement `evalSetup`**

`src/+gdsge/+parser/evalSetup.m`:
```matlab
function GDSGE_WS = evalSetup(GDSGE_SCRIPT)
% EVALSETUP  Eval a setup script in this isolated function workspace and capture
%   every resulting variable into a struct. On error, raise
%   gdsge:parser:setupEvalFailed echoing the numbered script lines. Internal
%   locals use a GDSGE_ prefix (a reserved prefix in gmod) and are excluded from
%   the captured workspace.
try
    eval(GDSGE_SCRIPT);
catch GDSGE_ME
    GDSGE_LINES = regexp(GDSGE_SCRIPT, '\n', 'split');
    GDSGE_MSG = sprintf('Error evaluating setup code:\n');
    for GDSGE_I = 1:numel(GDSGE_LINES)
        GDSGE_MSG = [GDSGE_MSG, sprintf('%3d\t%s\n', GDSGE_I, GDSGE_LINES{GDSGE_I})]; %#ok<AGROW>
    end
    GDSGE_MSG = [GDSGE_MSG, sprintf('\n%s', GDSGE_ME.message)];
    error('gdsge:parser:setupEvalFailed', '%s', GDSGE_MSG);
end
GDSGE_VARS = who;
GDSGE_WS = struct();
for GDSGE_I = 1:numel(GDSGE_VARS)
    GDSGE_NAME = GDSGE_VARS{GDSGE_I};
    if strncmp(GDSGE_NAME, 'GDSGE_', 6); continue; end
    GDSGE_WS.(GDSGE_NAME) = eval(GDSGE_NAME);
end
end
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tEvalSetup.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 4 tests.

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+parser/defaultSetupCode.m src/+gdsge/+parser/evalSetup.m tests/parser/tEvalSetup.m
git commit -m "feat(parser): eval sandbox (defaultSetupCode + isolated evalSetup with located errors)"
```

---

## Task 6: `parseVarDecls` — pure structural parse of declarations

**Files:**
- Create: `src/+gdsge/+parser/parseVarDecls.m`
- Test: `tests/parser/tParseVarDecls.m`

- [ ] **Step 1: Write the failing test**

`tests/parser/tParseVarDecls.m`:
```matlab
classdef tParseVarDecls < matlab.unittest.TestCase
    methods (Test)
        function parsesNamesLengthsBoundsExprs(tc)
            declText = sprintf([ ...
                'parameters beta;\nbeta = 0.95;\n', ...
                'var_shock z;\nshock_num = 2;\nz = [1 2];\n', ...
                'var_state w1;\nw1 = linspace(0,1,5);\n', ...
                'var_policy c1 w1n[8];\n', ...
                'inbound c1 0 1;\ninbound ps 0 3 adaptive(1.5);\n', ...
                'var_aux a;\nvar_interp f1;\ninitial f1 0.0;\nf1 = c1;\n', ...
                'var_output c1;']);
            d = gdsge.parser.parseVarDecls(declText);
            tc.verifyEqual(d.paramNames, {'beta'});
            tc.verifyEqual(d.shockNames, {'z'});
            tc.verifyEqual(d.stateNames, {'w1'});
            tc.verifyEqual(d.policy{2}.name, 'w1n');
            tc.verifyEqual(d.policy{2}.length, 8);
            tc.verifyEqual(d.bounds{1}.lower, '0');
            tc.verifyEqual(d.bounds{2}.adaptiveFactor, 1.5);
            tc.verifyEqual(d.interpNames, {'f1'});
            tc.verifyEqual(d.interpInitial{1}.expr, '0.0');
            tc.verifyEqual(d.interpUpdate{1}.name, 'f1');
            tc.verifyEqual(d.interpUpdate{1}.expr, 'c1');
            tc.verifyEqual(d.gridText{1}.name, 'w1');
            tc.verifyEqual(d.gridText{1}.expr, 'linspace(0,1,5)');
        end
        function setupStmtsIncludeValuesExcludeInterpUpdate(tc)
            declText = sprintf(['var_state w1;\nw1 = linspace(0,1,5);\n', ...
                'var_interp f1;\ninitial f1 0.0;\nf1 = c1;\nbeta = 0.95;']);
            d = gdsge.parser.parseVarDecls(declText);
            tc.verifyTrue(any(contains(d.setupStmts, 'beta = 0.95')));
            tc.verifyTrue(any(contains(d.setupStmts, 'w1 = linspace')));
            tc.verifyFalse(any(contains(d.setupStmts, 'f1 = c1')));  % interp update not eval'd
        end
        function unknownStatementErrors(tc)
            tc.verifyError(@() gdsge.parser.parseVarDecls('frobnicate the thing'), ...
                'gdsge:parser:unknownDeclaration');
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tParseVarDecls.m')); disp(table(r))"`
Expected: FAIL — `gdsge.parser.parseVarDecls` undefined.

- [ ] **Step 3: Implement `parseVarDecls`**

`src/+gdsge/+parser/parseVarDecls.m`:
```matlab
function decls = parseVarDecls(declText)
% PARSEVARDECLS  Pure structural parse of the declaration region. No eval.
%   Returns a struct of name sets, sized policy/aux items, bounds, interp
%   initial/update expressions, state grid text, and the residual setup
%   statements (param/shock/grid value assignments) to be eval'd later.
stmts = gdsge.parser.splitStatements(declText);

decls = struct('paramNames',{{}}, 'shockNames',{{}}, 'stateNames',{{}}, ...
    'policy',{{}}, 'aux',{{}}, 'interpNames',{{}}, 'outputNames',{{}}, ...
    'tensorNames',{{}}, 'otherNames',{{}}, 'bounds',{{}}, ...
    'interpInitial',{{}}, 'interpUpdate',{{}}, 'gridText',{{}}, 'setupStmts',{{}});

assignIdx = [];

% ---- Pass A: keyword declarations ----------------------------------------
for i = 1:numel(stmts)
    st = stmts{i};
    kw = regexp(st, '^([A-Za-z_]\w*)', 'tokens', 'once');
    key = ''; if ~isempty(kw); key = kw{1}; end
    switch key
        case 'parameters'
            decls.paramNames = [decls.paramNames, splitNames(rest(st,key))];
        case 'var_shock'
            decls.shockNames = [decls.shockNames, splitNames(rest(st,key))];
        case 'var_state'
            decls.stateNames = [decls.stateNames, splitNames(rest(st,key))];
        case 'var_policy'
            decls.policy = [decls.policy, parseSizedNames(rest(st,key))];
        case 'var_aux'
            decls.aux = [decls.aux, parseSizedNames(rest(st,key))];
        case 'var_interp'
            decls.interpNames = [decls.interpNames, splitNames(rest(st,key))];
        case 'var_output'
            decls.outputNames = [decls.outputNames, splitNames(rest(st,key))];
        case 'var_tensor'
            decls.tensorNames = [decls.tensorNames, splitNames(rest(st,key))];
        case 'var_others'
            decls.otherNames = [decls.otherNames, splitNames(rest(st,key))];
        case {'inbound','inbound_init'}
            decls.bounds{end+1} = parseBound(rest(st,key)); %#ok<AGROW>
        case 'initial'
            decls.interpInitial{end+1} = parseInitial(rest(st,key)); %#ok<AGROW>
        case {'var_policy_init','var_aux_init','model_init'}
            % recognized keyword; full handling is Phase 7 (HL1996 has none)
        otherwise
            assignIdx(end+1) = i; %#ok<AGROW>
    end
end

% ---- Pass B: classify assignment statements ------------------------------
for j = 1:numel(assignIdx)
    st = stmts{assignIdx(j)};
    lhs = assignLHS(st);
    if isempty(lhs)
        error('gdsge:parser:unknownDeclaration', 'Unrecognized declaration: "%s"', firstLine(st));
    end
    rhs = strtrim(assignRHS(st));
    if ismember(lhs, decls.interpNames)
        decls.interpUpdate{end+1} = struct('name', lhs, 'expr', rhs); %#ok<AGROW>
    else
        decls.setupStmts{end+1} = st; %#ok<AGROW>
        if ismember(lhs, decls.stateNames)
            decls.gridText{end+1} = struct('name', lhs, 'expr', rhs); %#ok<AGROW>
        end
    end
end
end

% ===== locals =============================================================
function r = rest(st, key)
r = strtrim(st(numel(key)+1:end));
end

function names = splitNames(s)
names = regexp(strtrim(s), '\s+', 'split');
names = names(~cellfun(@isempty, names));
end

function items = parseSizedNames(s)
toks = splitNames(s);
items = cell(1, numel(toks));
for i = 1:numel(toks)
    m = regexp(toks{i}, '^(\w+)(\[(\d+)\])?$', 'tokens', 'once');
    if isempty(m)
        error('gdsge:parser:unknownDeclaration', 'Bad variable token "%s"', toks{i});
    end
    len = 1; if ~isempty(m{3}); len = str2double(m{3}); end
    items{i} = struct('name', m{1}, 'length', len);
end
end

function b = parseBound(s)
toks = splitNames(s);
b = struct('name', toks{1}, 'lower', toks{2}, 'upper', toks{3});
adp = regexp(s, 'adaptive\(\s*([0-9.eE+-]+)\s*\)', 'tokens', 'once');
if ~isempty(adp)
    b.adaptiveFactor = str2double(adp{1});
end
end

function o = parseInitial(s)
[nm, rem] = strtok(strtrim(s));
o = struct('name', nm, 'expr', strtrim(rem));
end

function lhs = assignLHS(st)
m = regexp(firstLine(st), '^\s*([A-Za-z_]\w*)\s*(\([^)]*\)|\[[^\]]*\])?\s*=', 'tokens', 'once');
if isempty(m); lhs = ''; else; lhs = m{1}; end
end

function rhs = assignRHS(st)
idx = strfind(st, '=');
rhs = st(idx(1)+1:end);
end

function l = firstLine(st)
parts = regexp(st, '\n', 'split');
l = parts{1};
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tParseVarDecls.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 3 tests.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/parseVarDecls.m tests/parser/tParseVarDecls.m
git commit -m "feat(parser): parseVarDecls (structural pass: names, bounds, interp exprs, setup split)"
```

---

## Task 7: `parseDeclarations` — eval + assemble IR sections (with slots)

**Files:**
- Create: `src/+gdsge/+parser/parseDeclarations.m`
- Test: `tests/parser/tParseDeclarations.m`

- [ ] **Step 1: Write the failing test**

`tests/parser/tParseDeclarations.m`:
```matlab
classdef tParseDeclarations < matlab.unittest.TestCase
    methods (Test)
        function buildsSectionsWithSlotsAndValues(tc)
            declText = sprintf([ ...
                'parameters beta;\nbeta = 0.95;\n', ...
                'var_shock z;\nshock_num = 2;\nz = [1 2];\n', ...
                'shock_trans = [3 1\n1 3];\n', ...
                'shock_trans = shock_trans ./ repmat(sum(shock_trans,2),[1,shock_num]);\n', ...
                'var_state k;\nk = linspace(0,1,3);\n', ...
                'var_policy c w1n[2];\ninbound c 0 1;\ninbound w1n -1 1;\n', ...
                'var_aux a;\nvar_interp f;\ninitial f 0.0;\nf = c;\n', ...
                'var_output c;']);
            out = gdsge.parser.parseDeclarations(declText);
            tc.verifyEqual(out.params{1}.name, 'beta');
            tc.verifyEqual(out.params{1}.value, 0.95);
            tc.verifyEqual(out.shocks.count, 2);
            tc.verifyEqual(out.shocks.values.z, [1 2]);
            tc.verifyEqual(sum(out.shocks.transitions.shock_trans, 2), [1;1], 'AbsTol', 1e-12);
            tc.verifyEqual(out.states.grids.k, 'linspace(0,1,3)');
            tc.verifyEqual(out.variables.policy{1}.slot, [1 1]);   % c
            tc.verifyEqual(out.variables.policy{2}.slot, [2 3]);   % w1n[2]
            tc.verifyEqual(out.variables.aux{1}.slot, [1 1]);      % a (own slot space)
            tc.verifyEqual(out.interp{1}.name, 'f');
            tc.verifyEqual(out.interp{1}.args, {'k'});
            tc.verifyEqual(out.interp{1}.initialExpr, '0.0');
            tc.verifyEqual(out.interp{1}.updateExpr, 'c');
        end
        function undefinedParameterErrors(tc)
            declText = sprintf('parameters beta;\nvar_state k;\nk = linspace(0,1,3);');
            tc.verifyError(@() gdsge.parser.parseDeclarations(declText), ...
                'gdsge:parser:setupEvalFailed');
        end
        function stateWithoutGridErrors(tc)
            declText = sprintf('var_state k;\nvar_policy c;\ninbound c 0 1;');
            tc.verifyError(@() gdsge.parser.parseDeclarations(declText), ...
                'gdsge:parser:missingGrid');
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tParseDeclarations.m')); disp(table(r))"`
Expected: FAIL — `gdsge.parser.parseDeclarations` undefined.

- [ ] **Step 3: Implement `parseDeclarations`**

`src/+gdsge/+parser/parseDeclarations.m`:
```matlab
function out = parseDeclarations(declText)
% PARSEDECLARATIONS  Declaration region -> IR sections (params, shocks, states,
%   variables with flat slot layout, bounds, interp) plus the eval'd flag
%   workspace .ws (for resolveOptions). Evals param/shock/grid setup in an
%   isolated workspace seeded with parser defaults; interp updates and grids are
%   carried as text, never eval'd.
decls = gdsge.parser.parseVarDecls(declText);

% splitStatements stripped the trailing ';' from each statement; re-terminate
% so eval runs cleanly (no echo) and multi-line matrices end unambiguously.
setupBody = strjoin(decls.setupStmts, sprintf(';\n'));
if ~isempty(setupBody); setupBody = [setupBody, ';']; end
script = [gdsge.parser.defaultSetupCode(), sprintf('\n'), setupBody];
ws = gdsge.parser.evalSetup(script);

% ----- params -------------------------------------------------------------
params = cell(1, numel(decls.paramNames));
for i = 1:numel(decls.paramNames)
    nm = decls.paramNames{i};
    if ~isfield(ws, nm)
        error('gdsge:parser:setupEvalFailed', 'Parameter "%s" not defined', nm);
    end
    params{i} = struct('name', nm, 'value', double(ws.(nm)));
end

% ----- shocks -------------------------------------------------------------
shockVals = struct();
for i = 1:numel(decls.shockNames)
    nm = decls.shockNames{i};
    if ~isfield(ws, nm)
        error('gdsge:parser:setupEvalFailed', 'Shock variable "%s" not defined', nm);
    end
    shockVals.(nm) = double(ws.(nm));
end
shocks = struct('names', {decls.shockNames}, 'count', double(ws.shock_num), ...
    'values', shockVals, 'transitions', struct('shock_trans', double(ws.shock_trans)));

% ----- states -------------------------------------------------------------
grids = struct();
for i = 1:numel(decls.stateNames)
    grids.(decls.stateNames{i}) = gridFor(decls, decls.stateNames{i});
end
states = struct('names', {decls.stateNames}, 'grids', grids);

% ----- variables (+ slots) ------------------------------------------------
variables = struct( ...
    'policy', {assignSlots(decls.policy)}, ...
    'aux',    {assignSlots(decls.aux)}, ...
    'interp', {decls.interpNames}, ...
    'tensor', {decls.tensorNames}, ...
    'output', {decls.outputNames}, ...
    'others', {decls.otherNames});

% ----- interp objects -----------------------------------------------------
interp = cell(1, numel(decls.interpNames));
for i = 1:numel(decls.interpNames)
    nm = decls.interpNames{i};
    interp{i} = struct('name', nm, 'args', {decls.stateNames}, ...
        'initialExpr', lookupExpr(decls.interpInitial, nm), ...
        'updateExpr',  lookupExpr(decls.interpUpdate,  nm));
end

out = struct('params', {params}, 'shocks', shocks, 'states', states, ...
    'variables', variables, 'bounds', {decls.bounds}, 'interp', {interp}, 'ws', ws);
end

% ===== locals =============================================================
function items = assignSlots(rawItems)
items = cell(1, numel(rawItems));
pos = 1;
for i = 1:numel(rawItems)
    len = rawItems{i}.length;
    items{i} = struct('name', rawItems{i}.name, 'length', len, 'slot', [pos, pos+len-1]);
    pos = pos + len;
end
end

function g = gridFor(decls, name)
for i = 1:numel(decls.gridText)
    if strcmp(decls.gridText{i}.name, name); g = decls.gridText{i}.expr; return; end
end
error('gdsge:parser:missingGrid', 'State "%s" has no grid assignment', name);
end

function e = lookupExpr(lst, name)
e = '';
for i = 1:numel(lst)
    if strcmp(lst{i}.name, name); e = lst{i}.expr; return; end
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tParseDeclarations.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 3 tests.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/parseDeclarations.m tests/parser/tParseDeclarations.m
git commit -m "feat(parser): parseDeclarations (eval sandbox + IR sections with slot layout)"
```

---

## Task 8: `parseSimulate` — the simulate block

**Files:**
- Create: `src/+gdsge/+parser/parseSimulate.m`
- Test: `tests/parser/tParseSimulate.m`

- [ ] **Step 1: Write the failing test**

`tests/parser/tParseSimulate.m`:
```matlab
classdef tParseSimulate < matlab.unittest.TestCase
    methods (Test)
        function parsesAllSimulateFields(tc)
            blk = sprintf(['  num_periods = 10000;\n  num_samples = 24;\n', ...
                '  initial w1 0.5;\n  initial shock 1;\n', ...
                '  var_simu c1 c2 ps;\n  w1'' = w1n'';']);
            s = gdsge.parser.parseSimulate(blk);
            tc.verifyEqual(s.numPeriods, 10000);
            tc.verifyEqual(s.numSamples, 24);
            tc.verifyEqual(numel(s.initial), 2);
            tc.verifyEqual(s.initial{1}.var, 'w1');
            tc.verifyEqual(s.initial{1}.value, '0.5');
            tc.verifyEqual(s.initial{2}.var, 'shock');
            tc.verifyEqual(s.varSimu, {'c1','c2','ps'});
            tc.verifyEqual(s.transitions{1}.state, 'w1');
            tc.verifyEqual(s.transitions{1}.expr, 'w1n');
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tParseSimulate.m')); disp(table(r))"`
Expected: FAIL — `gdsge.parser.parseSimulate` undefined.

- [ ] **Step 3: Implement `parseSimulate`**

`src/+gdsge/+parser/parseSimulate.m`:
```matlab
function sim = parseSimulate(blockText)
% PARSESIMULATE  The simulate; ... end; block -> IR simulate section. All exprs
%   are carried as opaque text (no AST). The future-marker quote is stripped
%   from transition state names and right-hand sides.
sim = struct('numPeriods', [], 'numSamples', [], ...
    'initial', {{}}, 'varSimu', {{}}, 'transitions', {{}});
stmts = gdsge.parser.splitStatements(blockText);
for i = 1:numel(stmts)
    st = stmts{i};
    if ~isempty(regexp(st, '^num_periods\s*=', 'once'))
        sim.numPeriods = str2double(strtrim(afterEq(st)));
    elseif ~isempty(regexp(st, '^num_samples\s*=', 'once'))
        sim.numSamples = str2double(strtrim(afterEq(st)));
    elseif ~isempty(regexp(st, '^initial\s+', 'once'))
        m = regexp(st, '^initial\s+(\w+)\s+(.+)$', 'tokens', 'once');
        sim.initial{end+1} = struct('var', m{1}, 'value', strtrim(m{2})); %#ok<AGROW>
    elseif ~isempty(regexp(st, '^var_simu\s+', 'once'))
        toks = regexp(strtrim(st(numel('var_simu')+1:end)), '\s+', 'split');
        sim.varSimu = [sim.varSimu, toks(~cellfun(@isempty, toks))];
    else
        m = regexp(st, '^([A-Za-z_]\w*)''?\s*=\s*(.+)$', 'tokens', 'once');
        if ~isempty(m)
            sim.transitions{end+1} = struct('state', m{1}, ...
                'expr', strtrim(regexprep(m{2}, '''', ''))); %#ok<AGROW>
        end
    end
end
end

function r = afterEq(st)
idx = strfind(st, '=');
r = st(idx(1)+1:end);
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tParseSimulate.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 1 test.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/parseSimulate.m tests/parser/tParseSimulate.m
git commit -m "feat(parser): parseSimulate (num_periods/samples, initial, var_simu, transitions)"
```

---

## Task 9: `resolveOptions` — flag workspace → curated IR options

**Files:**
- Create: `src/+gdsge/+parser/resolveOptions.m`
- Test: `tests/parser/tResolveOptions.m`

- [ ] **Step 1: Write the failing test**

`tests/parser/tResolveOptions.m`:
```matlab
classdef tResolveOptions < matlab.unittest.TestCase
    methods (Test)
        function mapsSplineDefaults(tc)
            ws = gdsge.parser.evalSetup(gdsge.parser.defaultSetupCode());
            o = gdsge.parser.resolveOptions(ws);
            tc.verifyEqual(o.interpMethod, 'spline');
            tc.verifyEqual(o.interpOrder, 4);
            tc.verifyEqual(o.extrapOrder, 2);
            tc.verifyEqual(o.tolEq, 1e-6);
            tc.verifyEqual(o.simuResolve, 1);
            tc.verifyEqual(o.simuInterp, 0);
            tc.verifyEqual(o.printFreq, 10);
            tc.verifyEqual(o.saveFreq, 10);
            tc.verifyTrue(o.numThreads >= 1);
            tc.verifyFalse(isfield(o, 'asgMaxLevel'));   % spline omits ASG fields
        end
        function mapsAsgWithLevels(tc)
            ws = gdsge.parser.evalSetup([gdsge.parser.defaultSetupCode(), ...
                sprintf('\nUSE_SPLINE = 0;\nUSE_ASG = 1;')]);
            o = gdsge.parser.resolveOptions(ws);
            tc.verifyEqual(o.interpMethod, 'asg');
            tc.verifyEqual(o.asgMaxLevel, 10);
            tc.verifyEqual(o.asgThreshold, 1e-2);
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tResolveOptions.m')); disp(table(r))"`
Expected: FAIL — `gdsge.parser.resolveOptions` undefined.

- [ ] **Step 3: Implement `resolveOptions`**

`src/+gdsge/+parser/resolveOptions.m`:
```matlab
function o = resolveOptions(ws)
% RESOLVEOPTIONS  Map the eval'd flag workspace to the curated IR options.
%   ASG level/threshold are emitted only for the asg method (matching the
%   reference IR, which omits them for spline).
if getf(ws, 'USE_ASG', 0)
    method = 'asg';
elseif getf(ws, 'USE_PCHIP', 0)
    method = 'pchip';
else
    method = 'spline';
end
o = struct();
o.interpMethod = method;
o.interpOrder  = getf(ws, 'INTERP_ORDER', 4);
o.extrapOrder  = getf(ws, 'EXTRAP_ORDER', 2);
if strcmp(method, 'asg')
    o.asgMaxLevel  = getf(ws, 'AsgMaxLevel', 10);
    o.asgThreshold = getf(ws, 'AsgThreshold', 1e-2);
end
o.tolEq       = getf(ws, 'TolEq', 1e-6);
o.numThreads  = getf(ws, 'NumThreads', 1);
o.simuResolve = getf(ws, 'SIMU_RESOLVE', 1);
o.simuInterp  = getf(ws, 'SIMU_INTERP', 0);
o.printFreq   = getf(ws, 'PrintFreq', 10);
o.saveFreq    = getf(ws, 'SaveFreq', 10);
end

function v = getf(ws, name, dflt)
if isfield(ws, name); v = double(ws.(name)); else; v = dflt; end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('parser','tResolveOptions.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 2 tests.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/resolveOptions.m tests/parser/tResolveOptions.m
git commit -m "feat(parser): resolveOptions (flag workspace -> curated IR options)"
```

---

## Task 10: `assemblePartialIR` + `parseFrontEnd` + HL1996 golden

**Files:**
- Create: `src/+gdsge/+parser/assemblePartialIR.m`, `src/+gdsge/+parser/parseFrontEnd.m`, `src/+gdsge/+parser/Contents.m`
- Delete: `src/+gdsge/+parser/.gitkeep`
- Test: `tests/HeatonLucas1996/parser/tFrontEndHL1996.m`

- [ ] **Step 1: Write the failing test**

`tests/HeatonLucas1996/parser/tFrontEndHL1996.m`:
```matlab
classdef tFrontEndHL1996 < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addIrFolder(tc)
            here = fileparts(mfilename('fullpath'));            % tests/HeatonLucas1996/parser
            irDir = fullfile(fileparts(here), 'ir');            % tests/HeatonLucas1996/ir
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(irDir));
        end
    end
    methods (Test)
        function partialIrValidatesAndMatchesReference(tc)
            here = fileparts(mfilename('fullpath'));
            gmodPath = fullfile(fileparts(here), 'HL1996.gmod');
            ir = gdsge.parser.parseFrontEnd(fileread(gmodPath), 'HL1996');

            % 1. it validates
            r = gdsge.ir.validate(ir);
            tc.verifyTrue(r.pass, strjoin(r.errors, ' | '));

            % 2. machine-derived numThreads is a positive integer
            tc.verifyGreaterThanOrEqual(ir.options.numThreads, 1);

            % 3. matches the reference, section by section, with the model body
            %    blanked and numThreads neutralized (canonicalize is whole-IR).
            ref = buildHL1996IR();
            ref.model.statements = {};
            ref.model.equations  = {};
            ref.options.numThreads = 0;
            ir.options.numThreads  = 0;

            ca = gdsge.ir.canonicalize(ir);
            cb = gdsge.ir.canonicalize(ref);
            secs = fieldnames(cb);
            for i = 1:numel(secs)
                tc.verifyTrue(isequaln(ca.(secs{i}), cb.(secs{i})), ...
                    sprintf('section differs from reference: %s', secs{i}));
            end
        end
        function modelBodyLeftForPhase3(tc)
            here = fileparts(mfilename('fullpath'));
            gmodPath = fullfile(fileparts(here), 'HL1996.gmod');
            ir = gdsge.parser.parseFrontEnd(fileread(gmodPath), 'HL1996');
            tc.verifyEmpty(ir.model.statements);
            tc.verifyEmpty(ir.model.equations);
        end
        function missingBoundIsReported(tc)
            % a policy var with no inbound must raise a located parser error
            gmod = sprintf([ ...
                'parameters beta;\nbeta = 0.95;\n', ...
                'var_shock z;\nshock_num = 1;\nz = 1;\nshock_trans = 1;\n', ...
                'var_state k;\nk = linspace(0,1,3);\n', ...
                'var_policy c;\n', ...                              % no inbound for c
                'model;\n  equations;\n    c;\n  end;\nend;\n', ...
                'simulate;\n  num_periods = 1;\n  num_samples = 1;\n', ...
                '  initial k 0;\n  initial shock 1;\n  var_simu c;\n  k'' = c;\nend;']);
            tc.verifyError(@() gdsge.parser.parseFrontEnd(gmod, 'tiny'), ...
                'gdsge:parser:missingBound');
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('HeatonLucas1996','parser','tFrontEndHL1996.m')); disp(table(r))"`
Expected: FAIL — `gdsge.parser.parseFrontEnd` / `assemblePartialIR` undefined.

- [ ] **Step 3: Implement `assemblePartialIR`**

`src/+gdsge/+parser/assemblePartialIR.m`:
```matlab
function ir = assemblePartialIR(modelName, decl, sim, options, hooks)
% ASSEMBLEPARTIALIR  Combine front-end sections into a schema-valid partial IR
%   (model statements/equations empty; Phase 3 fills them). Enforces the
%   bounds-completeness semantic gate, then runs gdsge.ir.validate.
ir.irVersion = '1.0.0';
ir.modelName = modelName;
ir.params    = decl.params;
ir.options   = options;
ir.shocks    = decl.shocks;
ir.states    = decl.states;
ir.variables = decl.variables;
ir.bounds    = decl.bounds;
ir.interp    = decl.interp;
ir.model     = struct('statements', {{}}, 'equations', {{}});
ir.simulate  = sim;
ir.hooks     = hooks;

% semantic gate: every policy variable needs an inbound
boundNames = cell(1, numel(decl.bounds));
for i = 1:numel(decl.bounds); boundNames{i} = decl.bounds{i}.name; end
for i = 1:numel(decl.variables.policy)
    pn = decl.variables.policy{i}.name;
    if ~ismember(pn, boundNames)
        error('gdsge:parser:missingBound', 'Policy variable "%s" has no inbound', pn);
    end
end

r = gdsge.ir.validate(ir);
if ~r.pass
    error('gdsge:parser:invalidIR', 'Partial IR failed validation:\n%s', ...
        strjoin(r.errors, sprintf('\n')));
end
end
```

- [ ] **Step 4: Implement `parseFrontEnd`**

`src/+gdsge/+parser/parseFrontEnd.m`:
```matlab
function ir = parseFrontEnd(gmodText, modelName)
% PARSEFRONTEND  gmod text + model name -> validated partial IR. The model;...end;
%   body (statements + equations) is left empty for Phase 3; everything else is
%   populated.
clean = gdsge.parser.preprocess(gmodText);
sb    = gdsge.parser.splitBlocks(clean);

decl    = gdsge.parser.parseDeclarations(sb.declText);
options = gdsge.parser.resolveOptions(decl.ws);

if ~isfield(sb.blocks, 'simulate')
    error('gdsge:parser:missingSimulate', 'No simulate block found');
end
sim = gdsge.parser.parseSimulate(sb.blocks.simulate);

hooks = hooksFromBlocks(sb.blocks);
ir = gdsge.parser.assemblePartialIR(modelName, decl, sim, options, hooks);
end

function h = hooksFromBlocks(blocks)
h = struct('preModel','','preIter','','postIter','', ...
           'preJacCode','','postJacCode','','cxx','');
map = struct('pre_model','preModel', 'pre_iter','preIter', 'post_iter','postIter', ...
             'pre_jac_code','preJacCode', 'post_jac_code','postJacCode');
keys = fieldnames(map);
for i = 1:numel(keys)
    if isfield(blocks, keys{i})
        h.(map.(keys{i})) = blocks.(keys{i});
    end
end
end
```

- [ ] **Step 5: Add the package help and drop the placeholder**

`src/+gdsge/+parser/Contents.m`:
```matlab
% +PARSER  GDSGE gmod front-end: text -> partial IR.
%   parseFrontEnd   - top-level: gmod text + model name -> validated partial IR
%   preprocess      - strip comments, join continuations, rewrite deprecated, guard macros
%   splitBlocks     - separate declaration lines from named blocks
%   splitStatements - bracket-aware split into logical statements
%   defaultSetupCode- parser-owned default flag values
%   evalSetup       - eval a setup script in an isolated workspace; capture vars
%   parseVarDecls   - structural parse of the declaration region (no eval)
%   parseDeclarations - eval + assemble params/shocks/states/variables/bounds/interp
%   parseSimulate   - the simulate block -> IR simulate section
%   resolveOptions  - flag workspace -> curated IR options
%   assemblePartialIR - combine sections into a validated partial IR
```

Delete the placeholder:
```bash
git rm src/+gdsge/+parser/.gitkeep
```

- [ ] **Step 6: Run the golden test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('HeatonLucas1996','parser','tFrontEndHL1996.m')); disp(table(r)); assert(all(~[r.Failed]))"`
Expected: PASS — 3 tests. If a `section differs` failure appears, encode both to compare:
`matlab -batch "addpath('src'); addpath(fullfile('tests','HeatonLucas1996','ir')); ir=gdsge.parser.parseFrontEnd(fileread(fullfile('tests','HeatonLucas1996','HL1996.gmod')),'HL1996'); ref=buildHL1996IR(); ref.model.statements={}; ref.model.equations={}; ir.options.numThreads=0; ref.options.numThreads=0; disp(gdsge.ir.encode(ir)); disp('==== REF ===='); disp(gdsge.ir.encode(ref));"`

- [ ] **Step 7: Run the FULL suite**

Run: `pwsh -File tests/run.ps1`
Expected: PASS — all Phase-0/1/2 tests green; exit 0.

- [ ] **Step 8: Commit**

```bash
git add src/+gdsge/+parser/assemblePartialIR.m src/+gdsge/+parser/parseFrontEnd.m src/+gdsge/+parser/Contents.m tests/HeatonLucas1996/parser/tFrontEndHL1996.m
git rm src/+gdsge/+parser/.gitkeep
git commit -m "feat(parser): assemblePartialIR + parseFrontEnd; HL1996 partial-IR golden green"
```

---

## Task 11: Mark Phase 2 complete

**Files:**
- Modify: `PROGRESS.md`

- [ ] **Step 1: Update the phase tracker**

In `PROGRESS.md`, change the Phase 2 line from `◐` to `☑` and add detail:
```
- ☑ **Phase 2 — Parser front-end** (macros-stub + lexer + block-split + declaration
  parser → partial IR) (done 2026-06-12)
  - ☑ `params` IR section added; HL1996 options corrected to real toolbox defaults
  - ☑ `gdsge.parser`: preprocess, splitBlocks, splitStatements, eval sandbox,
    parseVarDecls, parseDeclarations (+ slot layout), parseSimulate, resolveOptions,
    assemblePartialIR, parseFrontEnd
  - ☑ HL1996 partial-IR golden (section-by-section vs the reference IR) green
```
Change the next line to mark Phase 3 as NEXT:
```
- ☐ **Phase 3 — Model-expression parser** ... → full IR for HL1996) — NEXT
```
Add a Changelog entry at the top of the Changelog list:
```
- 2026-06-12: **Phase 2 complete.** `gdsge.parser` front-end turns HL1996.gmod into a
  schema-valid partial IR (everything but the model body), proven section-by-section
  against the corrected reference IR. Added the `params` IR section; corrected the
  HL1996 options to the real toolbox defaults. Next: Phase 3 (model-expression parser).
```

- [ ] **Step 2: Commit**

```bash
git add PROGRESS.md
git commit -m "docs: mark Phase 2 complete in PROGRESS.md"
```

---

## Self-review notes (coverage against the spec)

- Spec §3 `params` section → Task 1. Spec §4 options reconciliation (incl. fixture
  corrections) → Task 1 (fixture) + Task 9 (`resolveOptions`).
- Spec §5.1 `preprocess` → Task 2 (defaults moved to the eval sandbox, Task 5 — a deliberate
  refinement of §5.1: defaults belong to the workspace the setup script runs in, keeping
  `preprocess` pure text-cleaning).
- Spec §5.2 `splitBlocks` → Task 3. §5 statement tokenizing → Task 4 (`splitStatements`).
- Spec §5.3 declaration classification → Task 6 (`parseVarDecls`). §5.4 eval sandbox →
  Task 5 (`evalSetup`/`defaultSetupCode`) + Task 7 (orchestration). §5.5 slot layout →
  Task 7 (`assignSlots`). §5.6 `parseSimulate` → Task 8. §5.7 `assemblePartialIR` → Task 10.
- Spec §7 testing → unit tests in Tasks 2–9; golden partial-IR test in Task 10.
- Spec §8 error identifiers: `macroUnsupported` (Task 2), `unterminatedBlock` (Task 3),
  `setupEvalFailed` (Tasks 5, 7), `missingBound` (Task 10), `unknownDeclaration` (Task 6).
  `missingGrid` is an additional located error surfaced in Task 7.
- Spec §11 acceptance criteria 1–5 → Task 10 (validate + section match + suite green) and
  Task 1 (doc regen / no-drift).
