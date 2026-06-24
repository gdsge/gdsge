# Phase 7d Macro Engine Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the `preprocess.m` macro guard with a real macro engine covering the full backward-compat surface (`include`, `cinclude`, `#define`, `#strcat_comma`, `#mat`, `#foreach`, `#for`, `#if`).

**Architecture:** One new parser stage, `gdsge.parser.expandMacros(rawText, baseDir) -> struct('text', expandedText, 'cxxIncludes', cellstr)`, called as the **first** step of `parseFrontEnd` (before line-cleaning). Internally it runs ordered passes (each a local sub-function) mirroring the old engine's order. `cinclude` directives are collected into `cxxIncludes`, carried on a new optional IR field `hooks.cxxIncludes`, and joined into the existing `GDSGE_OTHER_INCLUDE` placeholder in the generated C++. Fidelity is **parity where the old engine worked**, with two cleanups (multi-token `#define` values, nested `#for`) and informative `gdsge:parser:macro*` errors where the old engine silently failed.

**Tech Stack:** MATLAB R2025b; `matlab.unittest` headless via `matlab -batch "cd('tests'); run_tests"`; the `+gdsge.parser`, `+gdsge.ir`, `+gdsge.codegen` packages.

**Spec:** `docs/superpowers/specs/2026-06-13-phase7d-macro-engine-design.md`

**Plan refinements vs spec (deliberate, behavior-preserving):**
- Spec §4 "pass 0 moves `GNDSGE→GDSGE`/deprecate out of `preprocess`": the plan **keeps** that rewrite in `preprocess`. It still runs (on the expanded text, just before `splitBlocks`), so the net effect is identical for every realistic input — no macro pass depends on canonical option/keyword spelling. This avoids churning `tPreprocess` and the direct-`preprocess` test paths. The deprecated option-alias rewrites (`InterpOrder` etc.) remain out of 7d (orthogonal to macros; no model uses them).
- Spec §10 "always-present `hooks.cxxIncludes`": the plan makes it **optional, present only when non-empty** (the `modelInit` precedent), so macro-free models' IR JSON snapshots stay byte-identical and need no regeneration.

---

## File Structure

**Create:**
- `src/+gdsge/+parser/expandMacros.m` — the macro engine: a thin ordered pipeline + one local sub-function per pass (`passCinclude`, `passInclude`, `passDefine`, `passStrcatComma`, `passMat`, `passForeach`, `passFor`, `passIf`) and the `evalMacroExpr` / `matchFor` helpers. One file, one responsibility (gmod text -> expanded text + cxxIncludes).
- `tests/parser/tExpandMacros.m` — unit tests for every pass, both cleanups, and the error taxonomy (inline strings + temp files; no solvable model needed).
- `tests/macros/HL1996_macro.gmod` — a macro-ized HL1996 that expands to the same declarations.
- `tests/macros/tMacroEquivIR.m` — IR-equivalence gate: `parseFrontEnd(HL1996_macro)` `isequalIR` `parseFrontEnd(HL1996)`.
- `tests/macros/tCincludeInjection.m` — structural gate: `cinclude` reaches `GDSGE_OTHER_INCLUDE` in the generated `.cpp`.

**Modify:**
- `src/+gdsge/+ir/schema.m:93-99` — add `'cxxIncludes', opt(fList(fText()))` to the `hooks` struct.
- `src/+gdsge/+parser/parseFrontEnd.m` — add optional `baseDir`; call `expandMacros` first; thread `cxxIncludes` into `hooks` (only when non-empty).
- `src/+gdsge/+parser/preprocess.m:31-38` — remove the macro guard.
- `src/+gdsge/+codegen/generateCxx.m:41` — fill `GDSGE_OTHER_INCLUDE` from `ir.hooks.cxxIncludes`.
- `src/+gdsge/+codegen/codegen.m:22` — pass the gmod's directory as `baseDir`.
- `tests/parser/tPreprocess.m:20-23` — remove the `macroDirectiveErrors` test (guard is gone).
- `PROGRESS.md` — mark Phase 7d done; changelog entry.

**Pass order inside `expandMacros` (fixed, regardless of implementation order):**
`passCinclude` → `passInclude` → `passDefine` → `passStrcatComma` → `passMat` → `passForeach` → `passFor` → `passIf`.

---

## Task 1: IR schema — optional `hooks.cxxIncludes`

**Files:**
- Modify: `src/+gdsge/+ir/schema.m:93-99`
- Test: `tests/ir/tHooksCxxIncludes.m` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/ir/tHooksCxxIncludes.m`:

```matlab
classdef tHooksCxxIncludes < matlab.unittest.TestCase
    % Phase 7d: hooks.cxxIncludes is an optional cellstr field that
    % validates and survives JSON round-trip; absence is also valid.
    methods (Test)
        function presentListValidatesAndRoundtrips(tc)
            ir = baseIR();
            ir.hooks.cxxIncludes = {'#include "foo.h"', '#include <bar>'};
            r = gdsge.ir.validate(ir);
            tc.verifyTrue(r.pass, strjoin(r.errors, ' | '));
            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));
        end
        function absentFieldStillValidates(tc)
            ir = baseIR();   % hooks has no cxxIncludes field
            r = gdsge.ir.validate(ir);
            tc.verifyTrue(r.pass, strjoin(r.errors, ' | '));
        end
    end
end

function ir = baseIR()
% Smallest schema-valid IR: reuse the HL1996 reference builder.
here = fileparts(mfilename('fullpath'));                 % tests/ir
irDir = fullfile(fileparts(here), 'HeatonLucas1996', 'ir');
addpath(irDir);
c = onCleanup(@() rmpath(irDir));
ir = buildHL1996IR();
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "cd('tests'); r = runtests('ir/tHooksCxxIncludes.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: FAIL — `presentListValidatesAndRoundtrips` errors because the validator rejects the unknown `cxxIncludes` field.

- [ ] **Step 3: Add the field to the schema**

In `src/+gdsge/+ir/schema.m`, edit the `hooks` block (lines 93-99):

```matlab
hooks = fStruct(structOf( ...
    'preModel',    opt(fText()), ...
    'preIter',     opt(fText()), ...
    'postIter',    opt(fText()), ...
    'preJacCode',  opt(fText()), ...
    'postJacCode', opt(fText()), ...
    'cxx',         opt(fText()), ...
    'cxxIncludes', opt(fList(fText()))));
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "cd('tests'); r = runtests('ir/tHooksCxxIncludes.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: PASS (2 passing).

- [ ] **Step 5: Confirm the doc generator has no drift** (the schema feeds `docs/ir-schema.md`)

Run: `matlab -batch "cd('tests'); r = runtests('ir'); disp(table(r)); exit(any([r.Failed]))"`
Expected: PASS. If a no-drift test fails, regenerate the doc: `matlab -batch "addpath('src'); gdsge.ir.gendoc('docs/ir-schema.md')"`, then re-run.

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+ir/schema.m tests/ir/tHooksCxxIncludes.m docs/ir-schema.md
git commit -m "feat(ir): optional hooks.cxxIncludes field (Phase 7d)"
```

---

## Task 2: Scaffold `expandMacros` (identity) + wire into `parseFrontEnd`

This makes `expandMacros` a no-op pass-through and inserts it into the pipeline, removes the `preprocess` macro guard, and adds the optional `baseDir`. All existing gates must stay green (identity on macro-free models).

**Files:**
- Create: `src/+gdsge/+parser/expandMacros.m`
- Modify: `src/+gdsge/+parser/parseFrontEnd.m`, `src/+gdsge/+parser/preprocess.m:31-38`, `src/+gdsge/+codegen/codegen.m:22`, `tests/parser/tPreprocess.m`
- Test: `tests/parser/tExpandMacros.m` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/parser/tExpandMacros.m`:

```matlab
classdef tExpandMacros < matlab.unittest.TestCase
    methods (Test)
        function identityOnMacroFreeText(tc)
            txt = sprintf('parameters beta;\nbeta = 0.95;  %% note\nvar_state k;\n');
            m = gdsge.parser.expandMacros(txt, pwd);
            tc.verifyEqual(m.text, txt);          % byte-identical: no macros present
            tc.verifyEmpty(m.cxxIncludes);
        end
        function identityOnHL1996Gmod(tc)
            here = fileparts(mfilename('fullpath'));            % tests/parser
            gmod = fileread(fullfile(fileparts(here), 'HeatonLucas1996', 'HL1996.gmod'));
            m = gdsge.parser.expandMacros(gmod, pwd);
            tc.verifyEqual(m.text, gmod);
            tc.verifyEmpty(m.cxxIncludes);
        end
        function returnsStructShape(tc)
            m = gdsge.parser.expandMacros('x=1;', pwd);
            tc.verifyTrue(isstruct(m));
            tc.verifyTrue(isfield(m, 'text') && isfield(m, 'cxxIncludes'));
            tc.verifyTrue(iscell(m.cxxIncludes));
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: FAIL — `gdsge.parser.expandMacros` is undefined.

- [ ] **Step 3: Create the identity scaffold**

Create `src/+gdsge/+parser/expandMacros.m`:

```matlab
function out = expandMacros(rawText, baseDir)
% EXPANDMACROS  Expand gmod preprocessor macros on raw text, BEFORE comment
%   stripping (old-engine order). Returns struct:
%     .text         expanded gmod text
%     .cxxIncludes  cellstr of C++ '#include ...' lines from cinclude directives
%   baseDir is the directory include()/cinclude() file names resolve against
%   (defaults to pwd). Raises gdsge:parser:macro* on malformed directives.
%   Passes, in order: cinclude, include, #define, #strcat_comma, #mat,
%   #foreach, #for, #if. Each pass is a targeted replacement, so text with no
%   macros is returned byte-identical.
if nargin < 2 || isempty(baseDir); baseDir = pwd; end

code = rawText;
cxxIncludes = {};

% (passes are inserted here in subsequent tasks, in the fixed order above)

out = struct('text', code, 'cxxIncludes', {cxxIncludes});
end
```

- [ ] **Step 4: Wire `expandMacros` into `parseFrontEnd` + remove the guard**

Replace the head of `src/+gdsge/+parser/parseFrontEnd.m` (lines 1-3) with:

```matlab
function ir = parseFrontEnd(gmodText, modelName, baseDir)
% PARSEFRONTEND  gmod text + model name -> validated, semantically checked IR.
%   Optional baseDir resolves include()/cinclude() file names (default pwd).
if nargin < 3 || isempty(baseDir); baseDir = pwd; end
macro = gdsge.parser.expandMacros(gmodText, baseDir);
clean = gdsge.parser.preprocess(macro.text);
sb    = gdsge.parser.splitBlocks(clean);
```

Then update the `hooks` assembly near the end of `parseFrontEnd` (the `hooks = hooksFromBlocks(sb.blocks);` line) to thread `cxxIncludes` in only when non-empty:

```matlab
hooks = hooksFromBlocks(sb.blocks);
if ~isempty(macro.cxxIncludes)
    hooks.cxxIncludes = macro.cxxIncludes;
end
ir = gdsge.parser.assemblePartialIR(modelName, decl, sim, options, hooks, model, modelInit);
```

In `src/+gdsge/+parser/preprocess.m`, delete the macro guard (lines 31-38, the `% 3. macro guard` comment block and its `for` loop). Renumber the remaining `% 4.` comment to `% 3.`. The function keeps comment-strip, continuation-join, and the `GNDSGE`/deprecate rewrite.

- [ ] **Step 5: Thread baseDir from the codegen driver**

In `src/+gdsge/+codegen/codegen.m`, change line 22 to pass the gmod's directory:

```matlab
ir = gdsge.parser.parseFrontEnd(fileread(gmodFile), modelName, fileparts(gmodFile));
```

- [ ] **Step 6: Update `tPreprocess` (guard removed)**

In `tests/parser/tPreprocess.m`, delete the `macroDirectiveErrors` test method (lines 20-23). Leave the other four tests unchanged (`preprocess` still strips comments, joins continuations, rewrites `GNDSGE`).

- [ ] **Step 7: Run the affected suites**

Run: `matlab -batch "cd('tests'); r = runtests({'parser/tExpandMacros.m','parser/tPreprocess.m','parser'}); disp(table(r)); exit(any([r.Failed]))"`
Expected: PASS — `tExpandMacros` (3) green; `tPreprocess` green; all front-end gates (`tFrontEndHL1996`, `tFullIRHL1996`, etc.) green because `expandMacros` is identity on macro-free gmods.

- [ ] **Step 8: Commit**

```bash
git add src/+gdsge/+parser/expandMacros.m src/+gdsge/+parser/parseFrontEnd.m \
        src/+gdsge/+parser/preprocess.m src/+gdsge/+codegen/codegen.m \
        tests/parser/tExpandMacros.m tests/parser/tPreprocess.m
git commit -m "feat(parser): expandMacros scaffold wired into parseFrontEnd (Phase 7d)"
```

---

## Task 3: `#define` pass (multi-token value cleanup + error)

**Files:**
- Modify: `src/+gdsge/+parser/expandMacros.m`
- Test: `tests/parser/tExpandMacros.m`

- [ ] **Step 1: Write the failing tests** (append to the `methods (Test)` block)

```matlab
function defineSimpleSubstitution(tc)
    m = gdsge.parser.expandMacros(sprintf('#define N 3\nvar_policy c[N];\n'), pwd);
    tc.verifyTrue(contains(m.text, 'c[3]'));
    tc.verifyFalse(contains(m.text, '#define'));
    tc.verifyFalse(contains(m.text, 'N]'));        % the token N was replaced
end
function defineWordBoundaryOnly(tc)
    % NN must not be touched when N is defined
    m = gdsge.parser.expandMacros(sprintf('#define N 3\nx = NN + N;\n'), pwd);
    tc.verifyTrue(contains(m.text, 'x = NN + 3;'));
end
function defineMultiTokenValue(tc)
    % cleanup: value may contain spaces (old required exactly 3 tokens)
    m = gdsge.parser.expandMacros(sprintf('#define EXPR a + b\ny = EXPR;\n'), pwd);
    tc.verifyTrue(contains(m.text, 'y = a + b;'));
end
function defineStripsTrailingComment(tc)
    m = gdsge.parser.expandMacros(sprintf('#define K 5   %% the count\nz = K;\n'), pwd);
    tc.verifyTrue(contains(m.text, 'z = 5;'));
    tc.verifyFalse(contains(m.text, 'count'));
end
function defineMalformedErrors(tc)
    tc.verifyError(@() gdsge.parser.expandMacros(sprintf('#define\nx=1;'), pwd), ...
        'gdsge:parser:macroDefineMalformed');
    tc.verifyError(@() gdsge.parser.expandMacros(sprintf('#define 9bad 3\nx=1;'), pwd), ...
        'gdsge:parser:macroDefineMalformed');
end
```

- [ ] **Step 2: Run to verify failure**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: FAIL — `#define` lines pass through unexpanded.

- [ ] **Step 3: Implement the pass**

In `expandMacros.m`, add the call in the pipeline (after the `cxxIncludes = {};` line, but note final order — `passDefine` runs after the include passes added later; for now it is the first active pass):

```matlab
code = passDefine(code);
```

Add the sub-function (before the final `end` of the file, alongside future passes):

```matlab
function code = passDefine(code)
% #define NAME <rest-of-line>  ->  word-boundary substitution of NAME by value.
%   Multi-token values allowed (trailing % comment stripped). Processed in
%   source order so a later value may reference an earlier macro name.
BEFORE = '((?<=^)|(?<=\W))';
AFTER  = '((?=$)|(?=\W))';
linePat = '(?<=(^|\n))[ \t]*#define\>[^\n]*';
while true
    [s, e] = regexp(code, linePat, 'start', 'end', 'once');
    if isempty(s); break; end
    line = code(s:e);
    rest = strtrim(regexprep(line, '^[ \t]*#define\>', ''));
    pc = strfind(rest, '%'); if ~isempty(pc); rest = strtrim(rest(1:pc(1)-1)); end
    sp = regexp(rest, '\s', 'once');
    if isempty(sp)
        error('gdsge:parser:macroDefineMalformed', ...
            '#define needs a NAME and a VALUE: "%s"', strtrim(line));
    end
    name  = rest(1:sp-1);
    value = strtrim(rest(sp+1:end));
    if isempty(regexp(name, '^[A-Za-z]\w*$', 'once'))
        error('gdsge:parser:macroDefineMalformed', ...
            'Invalid #define name "%s"', name);
    end
    code = [code(1:s-1), code(e+1:end)];                 % delete the #define line
    value = strrep(value, '$', '$$');                    % literal value in replacement
    code = regexprep(code, [BEFORE name AFTER], value);
end
end
```

- [ ] **Step 4: Run to verify pass**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: PASS (all `tExpandMacros` green).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/expandMacros.m tests/parser/tExpandMacros.m
git commit -m "feat(parser): #define macro with multi-token values (Phase 7d)"
```

---

## Task 4: `include("f")` pass (file splice)

**Files:**
- Modify: `src/+gdsge/+parser/expandMacros.m`
- Test: `tests/parser/tExpandMacros.m`

- [ ] **Step 1: Write the failing tests** (append)

```matlab
function includeSplicesFileContent(tc)
    dir = tempname; mkdir(dir);
    c = onCleanup(@() rmdir(dir, 's'));
    fid = fopen(fullfile(dir, 'frag.gmod'), 'w');
    fprintf(fid, 'var_state k;\nk = linspace(0,1,3);\n'); fclose(fid);
    m = gdsge.parser.expandMacros(sprintf('include("frag.gmod");\nvar_policy c;\n'), dir);
    tc.verifyTrue(contains(m.text, 'var_state k;'));
    tc.verifyTrue(contains(m.text, 'var_policy c;'));
    tc.verifyFalse(contains(m.text, 'include('));
end
function includeMissingFileErrors(tc)
    tc.verifyError(@() gdsge.parser.expandMacros('include("nope.gmod");', tempname), ...
        'gdsge:parser:macroIncludeNotFound');
end
```

- [ ] **Step 2: Run to verify failure**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: FAIL — `include(...)` is not expanded.

- [ ] **Step 3: Implement the pass**

In `expandMacros.m`, add the call **before** `passDefine` (include splices content that later passes process):

```matlab
code = passInclude(code, baseDir);
code = passDefine(code);
```

Add the sub-function:

```matlab
function code = passInclude(code, baseDir)
% include("f")  ->  splice the text of <baseDir>/f in place. Optional trailing ';'.
pat = '(?<=(\n|^))[ \t]*include\(\s*''([^'']*)''\s*\)\s*;?|(?<=(\n|^))[ \t]*include\(\s*"([^"]*)"\s*\)\s*;?';
while true
    [s, e, tok] = regexp(code, pat, 'start', 'end', 'tokens', 'once');
    if isempty(s); break; end
    name = tok{~cellfun(@isempty, tok)};      % the populated quote group
    fpath = name;
    if ~isAbsolute(fpath); fpath = fullfile(baseDir, name); end
    if ~exist(fpath, 'file')
        error('gdsge:parser:macroIncludeNotFound', ...
            'include file not found: %s', fpath);
    end
    content = fileread(fpath);
    code = [code(1:s-1), content, code(e+1:end)];
end
end

function tf = isAbsolute(p)
tf = ~isempty(regexp(p, '^([A-Za-z]:[\\/]|[\\/]{2}|[\\/])', 'once'));
end
```

> Note: gmod include paths use single OR double quotes; the regex handles both and the populated token group is selected. `isAbsolute` covers Windows drive paths and UNC/`/` roots.

- [ ] **Step 4: Run to verify pass**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/expandMacros.m tests/parser/tExpandMacros.m
git commit -m "feat(parser): include() file splice macro (Phase 7d)"
```

---

## Task 5: `cinclude` pass + C++ injection

`cinclude("f")` / `cinclude <f>` records a C++ `#include` line and strips the directive. The collected lines flow through `hooks.cxxIncludes` into `GDSGE_OTHER_INCLUDE`.

**Files:**
- Modify: `src/+gdsge/+parser/expandMacros.m`, `src/+gdsge/+codegen/generateCxx.m:41`
- Test: `tests/parser/tExpandMacros.m`, `tests/macros/tCincludeInjection.m` (create)

- [ ] **Step 1: Write the failing parser tests** (append to `tExpandMacros.m`)

```matlab
function cincludeQuotedCollectsAndStrips(tc)
    m = gdsge.parser.expandMacros(sprintf('cinclude("my_funcs.h");\nvar_state k;\n'), pwd);
    tc.verifyEqual(m.cxxIncludes, {'#include "my_funcs.h"'});
    tc.verifyFalse(contains(m.text, 'cinclude'));
    tc.verifyTrue(contains(m.text, 'var_state k;'));
end
function cincludeAngledCollectsAndStrips(tc)
    m = gdsge.parser.expandMacros(sprintf('cinclude <cmath>;\nx=1;\n'), pwd);
    tc.verifyEqual(m.cxxIncludes, {'#include <cmath>'});
    tc.verifyFalse(contains(m.text, 'cinclude'));
end
function cincludeMultipleInOrder(tc)
    m = gdsge.parser.expandMacros(sprintf('cinclude("a.h");\ncinclude <b>;\n'), pwd);
    tc.verifyEqual(m.cxxIncludes, {'#include "a.h"', '#include <b>'});
end
```

- [ ] **Step 2: Run to verify failure**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: FAIL — `cinclude` not handled.

- [ ] **Step 3: Implement the pass**

In `expandMacros.m`, add the call **first** (before `passInclude`):

```matlab
[code, cxxIncludes] = passCinclude(code, cxxIncludes);
code = passInclude(code, baseDir);
code = passDefine(code);
```

Add the sub-function:

```matlab
function [code, inc] = passCinclude(code, inc)
% cinclude("f") / cinclude <f>  ->  collect '#include "f"' / '#include <f>',
%   strip the directive. Processed top-to-bottom so include order is preserved.
patQ = '(?<=(\n|^))[ \t]*cinclude\(\s*"([^"]*)"\s*\)\s*;?';
patA = '(?<=(\n|^))[ \t]*cinclude\s*<\s*([^>]*?)\s*>\s*;?';
while true
    [sq, eq, tq] = regexp(code, patQ, 'start', 'end', 'tokens', 'once');
    [sa, ea, ta] = regexp(code, patA, 'start', 'end', 'tokens', 'once');
    if isempty(sq) && isempty(sa); break; end
    useQ = ~isempty(sq) && (isempty(sa) || sq < sa);
    if useQ
        inc{end+1} = ['#include "' tq{2} '"']; %#ok<AGROW>
        code = [code(1:sq-1), code(eq+1:end)];
    else
        inc{end+1} = ['#include <' ta{2} '>']; %#ok<AGROW>
        code = [code(1:sa-1), code(ea+1:end)];
    end
end
end
```

> The token group index is `2` because the lookbehind alternation `(\n|^)` is capturing group 1.

- [ ] **Step 4: Wire `cxxIncludes` into the generated C++**

In `src/+gdsge/+codegen/generateCxx.m`, replace line 41 (`'GDSGE_OTHER_INCLUDE', '';`) with a join of the IR field:

```matlab
    'GDSGE_OTHER_INCLUDE',    cxxIncludeText(ir); ...
```

Add a local helper at the bottom of `generateCxx.m` (after the main function's `end`):

```matlab
function s = cxxIncludeText(ir)
% Join hooks.cxxIncludes into a newline-separated block (empty when absent).
s = '';
if isfield(ir.hooks, 'cxxIncludes') && ~isempty(ir.hooks.cxxIncludes)
    s = strjoin(ir.hooks.cxxIncludes, sprintf('\n'));
end
end
```

- [ ] **Step 5: Write the C++ injection structural test**

Create `tests/macros/tCincludeInjection.m`:

```matlab
classdef tCincludeInjection < matlab.unittest.TestCase
    % Phase 7d: a cinclude directive reaches GDSGE_OTHER_INCLUDE in the
    % generated mex_<model>.cpp. Reuses HL1996's known-generatable IR with
    % cxxIncludes injected (no new solvable fixture needed).
    methods (TestClassSetup)
        function addIrFolder(tc)
            here = fileparts(mfilename('fullpath'));            % tests/macros
            irDir = fullfile(fileparts(here), 'HeatonLucas1996', 'ir');
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(irDir));
        end
    end
    methods (Test)
        function injectsIncludeLines(tc)
            ir = buildHL1996IR();
            ir.hooks.cxxIncludes = {'#include "demo_header.h"', '#include <cmath_demo>'};
            outDir = tempname; 
            c = onCleanup(@() rmdir(outDir, 's'));
            files = gdsge.codegen.generateCxx(ir, outDir);
            cpp = fileread(files.cppFile);
            tc.verifyTrue(contains(cpp, '#include "demo_header.h"'));
            tc.verifyTrue(contains(cpp, '#include <cmath_demo>'));
        end
        function noIncludesYieldsByteIdenticalSlot(tc)
            % macro-free model: GDSGE_OTHER_INCLUDE resolves empty (no stray lines)
            ir = buildHL1996IR();
            outDir = tempname; 
            c = onCleanup(@() rmdir(outDir, 's'));
            files = gdsge.codegen.generateCxx(ir, outDir);
            cpp = fileread(files.cppFile);
            tc.verifyFalse(contains(cpp, 'demo_header'));
        end
    end
end
```

- [ ] **Step 6: Run to verify pass**

Run: `matlab -batch "cd('tests'); r = runtests({'parser/tExpandMacros.m','macros/tCincludeInjection.m'}); disp(table(r)); exit(any([r.Failed]))"`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/+gdsge/+parser/expandMacros.m src/+gdsge/+codegen/generateCxx.m \
        tests/parser/tExpandMacros.m tests/macros/tCincludeInjection.m
git commit -m "feat(parser,codegen): cinclude -> GDSGE_OTHER_INCLUDE injection (Phase 7d)"
```

---

## Task 6: `evalMacroExpr` helper + `#strcat_comma` + `#mat`

These three are introduced together because `#strcat_comma` and `#mat` both evaluate MATLAB via the shared `evalMacroExpr` helper.

**Files:**
- Modify: `src/+gdsge/+parser/expandMacros.m`
- Test: `tests/parser/tExpandMacros.m`

- [ ] **Step 1: Write the failing tests** (append)

```matlab
function strcatCommaExpands(tc)
    m = gdsge.parser.expandMacros(sprintf('var_policy #strcat_comma{K, 1:3};\n'), pwd);
    tc.verifyTrue(contains(m.text, 'K(1), K(2), K(3)'));
end
function matInterpolatesNumericResult(tc)
    m = gdsge.parser.expandMacros(sprintf('w = linspace(0,1,#mat{200+1});\n'), pwd);
    tc.verifyTrue(contains(m.text, 'linspace(0,1,201)'));
end
function matCanUseDefinedMacro(tc)
    % #mat is expanded after #define, so it can see substituted values
    m = gdsge.parser.expandMacros(sprintf('#define NP 5\nx = #mat{NP*2};\n'), pwd);
    tc.verifyTrue(contains(m.text, 'x = 10;'));
end
function matEvalErrorIsReported(tc)
    tc.verifyError(@() gdsge.parser.expandMacros(sprintf('x = #mat{nosuchfn(1)};\n'), pwd), ...
        'gdsge:parser:macroEvalFailed');
end
```

- [ ] **Step 2: Run to verify failure**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: FAIL.

- [ ] **Step 3: Implement the helper and two passes**

In `expandMacros.m`, add the calls in the pipeline (after `passDefine`, in order `passStrcatComma` then `passMat`):

```matlab
code = passDefine(code);
code = passStrcatComma(code);
code = passMat(code);
```

Add the sub-functions:

```matlab
function code = passStrcatComma(code)
% #strcat_comma{args, range}  ->  e.g. {K, 1:3} -> "K(1), K(2), K(3)"
%   (faithful port of the old process_strcat_comma expression).
pat = '#strcat_comma\{([^{}]*),([^{}]*)\}';
while true
    [s, e, tok] = regexp(code, pat, 'start', 'end', 'tokens', 'once');
    if isempty(s); break; end
    args = strtrim(tok{1});
    nums = evalMacroExpr(strtrim(tok{2}));
    nums_str = strsplit(num2str(nums), ' ');
    nums_with_bracket = strcat('(', nums_str, ')');
    rep = strip(strjoin(strcat(args, nums_with_bracket, ',')), ',');
    code = [code(1:s-1), rep, code(e+1:end)];
end
end

function code = passMat(code)
% #mat{expr}  ->  char(string(eval(expr))) interpolated in place.
pat = '#mat\{([^\n}]*)\}';
while true
    [s, e, tok] = regexp(code, pat, 'start', 'end', 'tokens', 'once');
    if isempty(s); break; end
    rep = char(string(evalMacroExpr(tok{1})));
    code = [code(1:s-1), rep, code(e+1:end)];
end
end

function val = evalMacroExpr(GDSGE_EXPR)
% EVALMACROEXPR  Eval a macro expression in this isolated workspace and return
%   its value. Only literals/built-ins and already-substituted #define values
%   are in scope (macros expand before parameter/grid eval). Wraps failures as
%   gdsge:parser:macroEvalFailed.
try
    val = eval(GDSGE_EXPR);
catch GDSGE_ME
    error('gdsge:parser:macroEvalFailed', ...
        'Failed to evaluate macro expression "%s": %s', GDSGE_EXPR, GDSGE_ME.message);
end
end
```

> `passMat`'s pattern uses `[^\n}]*` so it stops at the first `}` on the line (sufficient for the documented `#mat` uses and avoids swallowing across braces).

- [ ] **Step 4: Run to verify pass**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/expandMacros.m tests/parser/tExpandMacros.m
git commit -m "feat(parser): #mat and #strcat_comma macros with eval sandbox (Phase 7d)"
```

---

## Task 7: `#foreach … #endfor` pass (recursive)

**Files:**
- Modify: `src/+gdsge/+parser/expandMacros.m`
- Test: `tests/parser/tExpandMacros.m`

- [ ] **Step 1: Write the failing tests** (append)

```matlab
function foreachExpandsValues(tc)
    src = sprintf('#foreach v in ps pb c1\n#v_future = #v;\n#endfor v\n');
    m = gdsge.parser.expandMacros(src, pwd);
    tc.verifyTrue(contains(m.text, 'ps_future = ps;'));
    tc.verifyTrue(contains(m.text, 'pb_future = pb;'));
    tc.verifyTrue(contains(m.text, 'c1_future = c1;'));
    tc.verifyFalse(contains(m.text, '#foreach'));
end
function foreachNestsByDistinctId(tc)
    src = sprintf('#foreach i in 1 2\n#foreach j in a b\nx#i_#j = 0;\n#endfor j\n#endfor i\n');
    m = gdsge.parser.expandMacros(src, pwd);
    tc.verifyTrue(contains(m.text, 'x1_a = 0;'));
    tc.verifyTrue(contains(m.text, 'x2_b = 0;'));
end
```

- [ ] **Step 2: Run to verify failure**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: FAIL.

- [ ] **Step 3: Implement the pass** (faithful recursive port of the old `rec_extract_loop_seg`)

In `expandMacros.m`, add the call after `passMat`:

```matlab
code = passMat(code);
code = passForeach(code);
```

Add the sub-function:

```matlab
function out = passForeach(code)
% #foreach id in v1 v2 ... \n body \n #endfor id  ->  body repeated with #id
%   substituted by each value. Recursive: nested #foreach (distinct ids) work.
allButWord = '[^a-zA-Z0-9_]';
BEFORE = ['((?<=^)|(?<=' allButWord '))'];
AFTER  = ['((?=$)|(?=' allButWord '))'];
code = [newline, code];
rx = '#foreach\s+(?<id>\w+)\s+in\s+([^\n]*)\n(.*?)#endfor\s+\k<id>';
toks   = regexp(code, rx, 'tokens');
starts = regexp(code, rx, 'start');
ends   = regexp(code, rx, 'end');
if isempty(toks)
    out = code(2:end);     % drop the sentinel newline we prepended
    return;
end
out = '';
for i = 1:numel(toks)
    id   = toks{i}{1};
    vals = strsplit(strtrim(toks{i}{2}), ' ');
    body = toks{i}{3};
    sub = '';
    for j = 1:numel(vals)
        sub = [sub, newline, ...
            regexprep(body, [BEFORE '#' id AFTER], strtrim(vals{j}))]; %#ok<AGROW>
    end
    if i == 1
        out = [out, code(1:starts(i)-1), newline]; %#ok<AGROW>
    else
        out = [out, code(ends(i-1)+1:starts(i)-1), newline]; %#ok<AGROW>
    end
    out = [out, passForeach(sub), newline]; %#ok<AGROW>   % recurse into the body
    if i == numel(toks)
        out = [out, code(ends(i)+1:end)]; %#ok<AGROW>
    end
end
end
```

- [ ] **Step 4: Run to verify pass**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/expandMacros.m tests/parser/tExpandMacros.m
git commit -m "feat(parser): recursive #foreach macro (Phase 7d)"
```

---

## Task 8: `#for … #end` pass (recursive, nested-`#for` cleanup + error)

**Files:**
- Modify: `src/+gdsge/+parser/expandMacros.m`
- Test: `tests/parser/tExpandMacros.m`

- [ ] **Step 1: Write the failing tests** (append)

```matlab
function forExpandsRange(tc)
    m = gdsge.parser.expandMacros(sprintf('#for i=1:3\nx#i = #i;\n#end\n'), pwd);
    tc.verifyTrue(contains(m.text, 'x1 = 1;'));
    tc.verifyTrue(contains(m.text, 'x3 = 3;'));
    tc.verifyFalse(contains(m.text, '#for'));
end
function forSubstitutesIdPlusOne(tc)
    m = gdsge.parser.expandMacros(sprintf('#for i=1:2\nk#i_to_#(i+1);\n#end\n'), pwd);
    tc.verifyTrue(contains(m.text, 'k1_to_2;'));
    tc.verifyTrue(contains(m.text, 'k2_to_3;'));
end
function forNestsCleanup(tc)
    % cleanup: nested #for is supported (old errored)
    src = sprintf('#for i=1:2\n#for j=1:2\np#i_#j;\n#end\n#end\n');
    m = gdsge.parser.expandMacros(src, pwd);
    tc.verifyTrue(contains(m.text, 'p1_1;'));
    tc.verifyTrue(contains(m.text, 'p2_2;'));
end
function forUnterminatedErrors(tc)
    tc.verifyError(@() gdsge.parser.expandMacros(sprintf('#for i=1:2\nx#i;\n'), pwd), ...
        'gdsge:parser:macroForUnterminated');
end
```

- [ ] **Step 2: Run to verify failure**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: FAIL.

- [ ] **Step 3: Implement the pass** (balanced matching + recursion)

In `expandMacros.m`, add the call after `passForeach`:

```matlab
code = passForeach(code);
code = passFor(code);
```

Add the sub-functions:

```matlab
function out = passFor(code)
% #for id=range \n body \n #end  ->  body repeated with #id and #(id+1)
%   substituted. Balanced #for/#end matching supports nesting (recursion).
out = '';
pos = 1;
while true
    s = tokenSearch(code, pos, '#for');
    if isempty(s)
        out = [out, code(pos:end)];
        return;
    end
    out = [out, code(pos:s-1)]; %#ok<AGROW>
    [bodyStart, bodyEnd, blockEnd] = matchFor(code, s);
    header = code(s:bodyStart-1);
    body   = code(bodyStart:bodyEnd);
    ht = regexp(header, '#for\s+(\w+)\s*=\s*([^\n]*)', 'tokens', 'once');
    if isempty(ht)
        error('gdsge:parser:macroForMalformed', ...
            'malformed #for header: "%s"', strtrim(header));
    end
    iter = ht{1};
    rangeExpr = strtrim(ht{2});
    pc = strfind(rangeExpr, '%'); 
    if ~isempty(pc); rangeExpr = strtrim(rangeExpr(1:pc(1)-1)); end
    vals = evalMacroExpr(rangeExpr);
    expanded = '';
    for v = vals(:)'
        b = regexprep(body, ['(?<![A-Za-z0-9_])#\(' iter '\+1\)'], int2str(v+1));
        b = regexprep(b, ['(?<![A-Za-z0-9_])#' iter '(?![A-Za-z0-9_])'], int2str(v));
        expanded = [expanded, passFor(b)]; %#ok<AGROW>      % recurse for nesting
    end
    out = [out, expanded]; %#ok<AGROW>
    pos = blockEnd + 1;
end
end

function [bodyStart, bodyEnd, blockEnd] = matchFor(code, s)
% s indexes the '#for' that opens a block. Returns body span and the index of
% the last char of the matching '#end'. Raises macroForUnterminated if none.
nl = regexp(code(s:end), '\n', 'once');
if isempty(nl)
    error('gdsge:parser:macroForUnterminated', 'no body/#end after #for');
end
bodyStart = s + nl;          % first char after the header newline
depth = 1; i = bodyStart; n = numel(code);
while i <= n
    if tokenAt(code, i, '#for'); depth = depth + 1; i = i + 4; continue; end
    if tokenAt(code, i, '#end')
        depth = depth - 1;
        if depth == 0; bodyEnd = i - 1; blockEnd = i + 3; return; end
        i = i + 4; continue;
    end
    i = i + 1;
end
error('gdsge:parser:macroForUnterminated', 'no matching #end for #for');
end

function s = tokenSearch(code, pos, tok)
% First index >= pos where tok appears as a standalone '#'-token.
s = [];
i = pos;
while i <= numel(code)
    if tokenAt(code, i, tok); s = i; return; end
    i = i + 1;
end
end

function tf = tokenAt(code, i, tok)
% True if tok occurs at index i and is not followed by a word char
% (so '#for' does not match inside '#format', '#end' not inside '#endfor').
n = numel(tok);
tf = false;
if i + n - 1 > numel(code); return; end
if ~strcmp(code(i:i+n-1), tok); return; end
if i + n <= numel(code)
    c = code(i+n);
    if ~isempty(regexp(c, '[A-Za-z0-9_]', 'once')); return; end
end
tf = true;
end
```

> `tokenAt` guards against `#end` matching inside `#endfor` (the `#foreach` pass already ran, so no `#endfor` remains, but the guard keeps `passFor` robust if passes are reordered).

- [ ] **Step 4: Run to verify pass**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/expandMacros.m tests/parser/tExpandMacros.m
git commit -m "feat(parser): recursive #for macro with nesting (Phase 7d)"
```

---

## Task 9: `#if … #endif` pass (+ error)

**Files:**
- Modify: `src/+gdsge/+parser/expandMacros.m`
- Test: `tests/parser/tExpandMacros.m`

- [ ] **Step 1: Write the failing tests** (append)

```matlab
function ifTrueKeepsBody(tc)
    m = gdsge.parser.expandMacros(sprintf('#if 1\nvar_policy c;\n#endif\n'), pwd);
    tc.verifyTrue(contains(m.text, 'var_policy c;'));
    tc.verifyFalse(contains(m.text, '#if'));
end
function ifFalseDropsBody(tc)
    m = gdsge.parser.expandMacros(sprintf('#if 0\nvar_policy junk;\n#endif\nvar_state k;\n'), pwd);
    tc.verifyFalse(contains(m.text, 'junk'));
    tc.verifyTrue(contains(m.text, 'var_state k;'));
end
function ifConditionCanUseDefine(tc)
    m = gdsge.parser.expandMacros(sprintf('#define FLAG 1\n#if FLAG\nkeep;\n#endif\n'), pwd);
    tc.verifyTrue(contains(m.text, 'keep;'));
end
function ifUnterminatedErrors(tc)
    tc.verifyError(@() gdsge.parser.expandMacros(sprintf('#if 1\nx=1;\n'), pwd), ...
        'gdsge:parser:macroIfUnterminated');
end
```

- [ ] **Step 2: Run to verify failure**

Run: `matlab -batch "cd('tests'); r = runtests('parser/tExpandMacros.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: FAIL.

- [ ] **Step 3: Implement the pass**

In `expandMacros.m`, add the call after `passFor` (last pass):

```matlab
code = passFor(code);
code = passIf(code);
```

Add the sub-function:

```matlab
function code = passIf(code)
% #if cond \n body \n #endif  ->  keep body iff eval(cond) is truthy, else drop.
%   Parity with old: non-nesting, non-greedy to the first #endif, no #else.
if isempty(regexp(code, '#if\>', 'once')); return; end
if isempty(regexp(code, '#endif\>', 'once'))
    error('gdsge:parser:macroIfUnterminated', '#if without #endif');
end
pat = '#if\s+([^\n]*)\n(.*?)#endif';
while true
    [s, e, tok] = regexp(code, pat, 'start', 'end', 'tokens', 'once');
    if isempty(s); break; end
    keep = evalMacroExpr(strtrim(tok{1}));
    if keep
        rep = tok{2};
    else
        rep = '';
    end
    code = [code(1:s-1), rep, code(e+1:end)];
end
end
```

> The leading `#if`-present / `#endif`-absent guard turns an unterminated `#if` into `macroIfUnterminated` rather than a silent no-op (the regex would simply not match).

- [ ] **Step 4: Run to verify pass + full parser suite**

Run: `matlab -batch "cd('tests'); r = runtests({'parser/tExpandMacros.m','parser'}); disp(table(r)); exit(any([r.Failed]))"`
Expected: PASS — all `tExpandMacros` and all parser gates green.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/expandMacros.m tests/parser/tExpandMacros.m
git commit -m "feat(parser): #if/#endif macro (Phase 7d)"
```

---

## Task 10: IR-equivalence gate — `HL1996_macro.gmod`

A macro-ized HL1996 that expands to the same declarations; parsing it yields an IR `isequalIR` to the plain HL1996's IR. This is the integration differential (no solve, no capture).

**Files:**
- Create: `tests/macros/HL1996_macro.gmod`, `tests/macros/tMacroEquivIR.m`

- [ ] **Step 1: Author the macro-ized fixture**

Create `tests/macros/HL1996_macro.gmod` — identical to `tests/HeatonLucas1996/HL1996.gmod` except for the marked macro uses (each chosen to expand to byte-identical declaration text):

```
% Parameters
parameters beta gamma Kb;
#define BETA 0.95
#define GAMMA 1.5
beta = BETA;  % discount factor
gamma = GAMMA;  % CRRA coefficient
Kb = -0.05;   % borrowing limit in ratio of aggregate output
% Shock variables
var_shock g d eta1;
% Shocks and transition matrix
shock_num = 8;
g = [.9904 1.0470 .9904 1.0470 .9904 1.0470 .9904 1.0470];
d = [.1402 .1437 .1561 .1599 .1402 .1437 .1561 .1599];
eta1 = [.3772 .3772 .3772 .3772 .6228 .6228 .6228 .6228];
shock_trans = [
  0.3932 0.2245 0.0793 0.0453 0.1365 0.0779 0.0275 0.0157
  0.3044 0.3470 0.0425 0.0484 0.1057 0.1205 0.0147 0.0168
  0.0484 0.0425 0.3470 0.3044 0.0168 0.0147 0.1205 0.1057
  0.0453 0.0793 0.2245 0.3932 0.0157 0.0275 0.0779 0.1365
  0.1365 0.0779 0.0275 0.0157 0.3932 0.2245 0.0793 0.0453
  0.1057 0.1205 0.0147 0.0168 0.3044 0.3470 0.0425 0.0484
  0.0168 0.0147 0.1205 0.1057 0.0484 0.0425 0.3470 0.3044
  0.0157 0.0275 0.0779 0.1365 0.0453 0.0793 0.2245 0.3932
  ];
shock_trans = shock_trans ./ repmat(sum(shock_trans,2),[1,shock_num]);
% State variables
var_state w1;  % wealth share
w1 = linspace(-0.05,1.05,#mat{200+1});
% Endogenous variables and bounds
var_policy c1 c2 s1p nb1p nb2p ms1 ms2 mb1 mb2 ps pb w1n[8];
inbound c1 0 1;
inbound c2 0 1;
inbound s1p 0.0 1.0;
inbound nb1p 0.0 1.0;   % nb1p=b1p-Kb
inbound nb2p 0.0 1.0;   
inbound ms1 0 1;        % Multilier for constraints
inbound ms2 0 1;
inbound mb1 0 1;
inbound mb2 0 1;
inbound ps 0 3 adaptive(1.5);
inbound pb 0 3 adaptive(1.5);
inbound w1n -0.5 1.5;
% Extra output variables
var_aux equity_premium;
% Interpolation objects
var_interp ps_future pb_future c1_future c2_future;
initial ps_future 0.0;
initial pb_future 0.0;
initial c1_future w1.*d+eta1;
initial c2_future (1-w1).*d+1-eta1;
#foreach v in ps pb c1 c2
#v_future = #v;
#endfor v
% Variables to be used in simulation if SIMU_RESOLVE=1
var_output c1 c2 ps pb equity_premium w1n;

model;
  % Interpolation
  [psn',pbn',c1n',c2n'] = GDSGE_INTERP_VEC'(w1n');
  % Expectations in Euler Equations
  es1 = GDSGE_EXPECT{g'^(1-gamma)*(c1n'/c1)^(-gamma)*(psn'+d')/ps};
  es2 = GDSGE_EXPECT{g'^(1-gamma)*(c2n'/c2)^(-gamma)*(psn'+d')/ps};
  eb1 = GDSGE_EXPECT{g'^(-gamma)*(c1n'/c1)^(-gamma)/pb};
  eb2 = GDSGE_EXPECT{g'^(-gamma)*(c2n'/c2)^(-gamma)/pb};
  % b transformation
  b1p = nb1p + Kb;  % Transform bond back
  b2p = nb2p + Kb;  
  s2p = 1-s1p;      % Market clear of shares
  % Budget constraint
  budget_1 = w1*(ps+d)+eta1 - c1 - ps*s1p - pb*b1p;
  budget_2 = (1-w1)*(ps+d)+(1-eta1) - c2 - ps*s2p - pb*b2p;
  % Consistency
  w1_consis' = (s1p*(psn'+d') + b1p/g')/(psn'+d') - w1n';
  % Extra output
  equity_premium = GDSGE_EXPECT{(psn'+d')/ps*g'} - 1/pb;
  equations;
    -1+beta*es1+ms1;
    -1+beta*es2+ms2;
    -1+beta*eb1+mb1;
    -1+beta*eb2+mb2;
    ms1*s1p;
    ms2*s2p;
    mb1*nb1p;
    mb2*nb2p;
    b1p+b2p;
    budget_1;
    budget_2;
    w1_consis';
  end;
end;

simulate;
  num_periods = 10000;
  num_samples = 24;
  initial w1 0.5;
  initial shock 1;
  var_simu c1 c2 ps pb equity_premium;
  w1' = w1n';
end;
```

Macro uses exercised: `#define` (BETA/GAMMA → `beta = 0.95;` / `gamma = 1.5;`), `#mat` (`#mat{200+1}` → `201`), `#foreach` (the four `_future = …;` assignments). Each expands to text identical to the plain HL1996 declarations.

- [ ] **Step 2: Write the gate**

Create `tests/macros/tMacroEquivIR.m`:

```matlab
classdef tMacroEquivIR < matlab.unittest.TestCase
    % Phase 7d integration differential: a macro-ized HL1996 expands to the
    % same model. Its IR must equal the plain HL1996 IR (which is already
    % differential-verified end-to-end against the old toolbox).
    methods (Test)
        function macroExpansionMatchesPlainIR(tc)
            here = fileparts(mfilename('fullpath'));               % tests/macros
            root = fileparts(here);                                % tests
            plainGmod = fileread(fullfile(root, 'HeatonLucas1996', 'HL1996.gmod'));
            macroGmod = fileread(fullfile(here, 'HL1996_macro.gmod'));

            irPlain = gdsge.parser.parseFrontEnd(plainGmod, 'HL1996');
            irMacro = gdsge.parser.parseFrontEnd(macroGmod, 'HL1996');

            tc.verifyTrue(gdsge.ir.isequalIR(irMacro, irPlain), ...
                'macro-expanded HL1996 IR differs from the plain HL1996 IR');
        end
        function macroFixtureHasNoCxxIncludes(tc)
            here = fileparts(mfilename('fullpath'));
            macroGmod = fileread(fullfile(here, 'HL1996_macro.gmod'));
            ir = gdsge.parser.parseFrontEnd(macroGmod, 'HL1996');
            tc.verifyFalse(isfield(ir.hooks, 'cxxIncludes'));   % no cinclude used
        end
    end
end
```

- [ ] **Step 3: Run the gate**

Run: `matlab -batch "cd('tests'); r = runtests('macros/tMacroEquivIR.m'); disp(table(r)); exit(any([r.Failed]))"`
Expected: PASS. If `macroExpansionMatchesPlainIR` fails, the failure message identifies a section mismatch — compare `gdsge.ir.canonicalize(irMacro)` vs `...(irPlain)` field by field (as `tFrontEndHL1996` does) to find the diverging declaration, and adjust the fixture so the expansion is byte-identical.

- [ ] **Step 4: Commit**

```bash
git add tests/macros/HL1996_macro.gmod tests/macros/tMacroEquivIR.m
git commit -m "test(macros): IR-equivalence gate for macro-ized HL1996 (Phase 7d)"
```

---

## Task 11: Full-suite regression + PROGRESS.md + branch wrap-up

**Files:**
- Modify: `PROGRESS.md`

- [ ] **Step 1: Run the entire suite**

Run: `matlab -batch "cd('tests'); run_tests"` then inspect exit code (0 = all pass) and `tests/results/junit.xml`.
Expected: all green, including HL1996 / 7a / 7b / 7c end-to-end gates and the new `tExpandMacros` / `tMacroEquivIR` / `tCincludeInjection` / `tHooksCxxIncludes`.

- [ ] **Step 2: Update PROGRESS.md**

In `PROGRESS.md`: change the Phase 7d line from `☐ … — NEXT —` to `☑ … (done 2026-06-13)` with a one-line summary, mark Phase 7e as the new **NEXT**, and add a changelog entry mirroring the 7c entry's style:

```markdown
- 2026-06-13: **Phase 7d complete.** Macro engine replaces the preprocess.m
  guard: new gdsge.parser.expandMacros stage (cinclude/include/#define/
  #strcat_comma/#mat/#foreach/#for/#if) run before line-cleaning; parity where
  the old engine worked, plus two cleanups (multi-token #define values, nested
  #for) and a gdsge:parser:macro* error taxonomy. cinclude flows through a new
  optional IR hooks.cxxIncludes field into the existing GDSGE_OTHER_INCLUDE
  placeholder (no template change). Verified by direct expandMacros unit tests
  (tExpandMacros), a cinclude C++-injection structural gate
  (tCincludeInjection, reusing HL1996's IR), and an IR-equivalence differential
  (tMacroEquivIR: a macro-ized HL1996 parses isequalIR the plain HL1996 IR) —
  zero new golden capture, since macros reduce to an already-covered model. All
  HL1996/7a/7b/7c gates stay green; full suite green. All on branch
  phase7d-macro-engine. Next: Phase 7e (var_tensor).
```

- [ ] **Step 3: Commit**

```bash
git add PROGRESS.md
git commit -m "docs(progress): Phase 7d complete — macro engine"
```

- [ ] **Step 4: Offer the merge** (do not merge unprompted)

Report the green full-suite result and offer to merge `phase7d-macro-engine` into `main` (mirroring the 7a/7b/7c merge pattern) or hand off for review.

---

## Self-Review (completed during planning)

**Spec coverage:** every spec item maps to a task — macro surface §2 → Tasks 3-9; `cinclude` coupling §6 → Tasks 1, 5; cleanups §5 → Tasks 3, 8; eval sandbox §3 → Task 6 (`evalMacroExpr`); testing §7 → Tasks 5 (structural), 10 (IR-equivalence), 2-9 (unit); error taxonomy §8 → Tasks 3, 4, 6, 8, 9; pipeline placement §4 → Task 2. Two spec refinements (deprecate stays in `preprocess`; `cxxIncludes` optional) are documented in the header with justification.

**Placeholder scan:** no TBD/TODO; every code step shows complete, runnable MATLAB and exact run commands.

**Type consistency:** `expandMacros` returns `struct('text', …, 'cxxIncludes', {cellstr})` everywhere; `hooks.cxxIncludes` is the cellstr field consumed by `generateCxx`'s `cxxIncludeText`; pass sub-functions (`passCinclude/passInclude/passDefine/passStrcatComma/passMat/passForeach/passFor/passIf`) and helpers (`evalMacroExpr/matchFor/tokenAt/tokenSearch/isAbsolute`) are referenced with consistent names and signatures across tasks. The fixed pipeline order (cinclude → include → define → strcat → mat → foreach → for → if) is assembled incrementally and ends correct.
