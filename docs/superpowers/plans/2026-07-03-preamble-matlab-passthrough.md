# Preamble MATLAB Passthrough Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restore the old toolbox's ability to put arbitrary MATLAB in the gmod declaration region (multi-output assignments, bare calls, control flow, re-opened declaration blocks), with a typo guard the old toolbox lacked, so `incomp.gmod` runs end-to-end.

**Architecture:** All parser changes live in `src/+gdsge/+parser/parseVarDecls.m`. Pass B stops rejecting statements that aren't `name = rhs` assignments and passes them through as setup text (after a near-miss keyword typo check); Pass C stops rejecting re-opened declaration kinds. The existing plumbing (setup sections → parse-time eval → verbatim replay in generated files) is untouched; the IR schema already allows repeated section kinds.

**Tech Stack:** MATLAB R2025b (`C:\Program Files\MATLAB\R2025b\bin\matlab.exe`, drive with `matlab -batch`), matlab.unittest.

**Spec:** `docs/superpowers/specs/2026-07-03-preamble-matlab-passthrough-design.md`

## Global Constraints

- **MATLAB path golden rule:** each run `addpath`s exactly one source (`src/` here) in its own `matlab -batch` process; never touch the persistent path.
- **One MATLAB process at a time** — each saturates all cores; run sequentially.
- PowerShell 7 (`pwsh`) is not installed; run tests via `matlab -batch "cd('tests'); run_tests"`. `tests/results/junit.xml` + exit code are authoritative (results.tap accumulates stale documents).
- **Backward compatible:** every existing `.gmod` must still run; existing IR snapshots and codegen goldens must be byte-identical after this change.
- TDD: write the failing test first; small, frequent commits.
- Do not commit any file from `D:\macro_local\exchange_rate_dynamics_debug` (private research code). `scratch/` is git-ignored — use it for the end-to-end run.
- `IterRslt_incomp.mat` is irreplaceable (long fine-tuned old-toolbox run). Never overwrite it; never attempt a from-scratch convergence run of incomp.

---

### Task 1: Pass B — permissive residual + probableTypo guard

**Files:**
- Modify: `src/+gdsge/+parser/parseVarDecls.m` (Pass B, lines ~58-78, plus two new local functions)
- Test: `tests/parser/tParseVarDecls.m`

**Interfaces:**
- Consumes: `gdsge.parser.parseVarDecls(declText)` — existing function, existing return struct.
- Produces: unchanged struct shape. New behavior: statements with no single-identifier LHS land in `decls.setupStmts` verbatim instead of erroring; new error id `gdsge:parser:probableTypo` for command-form statements whose first word is edit-distance ≤ 2 from a declaration keyword. Error id `gdsge:parser:unknownDeclaration` remains only for malformed sized-name tokens (`parseSizedNames`).

- [ ] **Step 1: Write the failing tests**

In `tests/parser/tParseVarDecls.m`, **replace** the existing `unknownStatementErrors` method (lines 37-40) with the following methods (keep everything else):

```matlab
        function commandFormFarFromKeywordsPassesThrough(tc)
            % old-toolbox behavior: unrecognized MATLAB passes through as setup text
            d = gdsge.parser.parseVarDecls('frobnicate the thing');
            tc.verifyEqual(d.setupStmts, {'frobnicate the thing'});
        end
        function multiOutputAssignmentPassesThrough(tc)
            declText = sprintf([ ...
                'var_shock z;\nshock_num = 2;\n', ...
                '[zg, wg] = ndgrid([0.9 1.1], 1);\nz = zg(:)'';']);
            d = gdsge.parser.parseVarDecls(declText);
            tc.verifyTrue(any(contains(d.setupStmts, '[zg, wg] = ndgrid')));
            tc.verifyTrue(any(contains(d.setupStmts, 'z = zg(:)')));
        end
        function tildeOutputPassesThrough(tc)
            d = gdsge.parser.parseVarDecls('[~, idx] = max([3 1 2])');
            tc.verifyEqual(d.setupStmts, {'[~, idx] = max([3 1 2])'});
        end
        function bareCallPassesThrough(tc)
            d = gdsge.parser.parseVarDecls(sprintf('rng(0);\na = 1;'));
            tc.verifyEqual(d.setupStmts, {'rng(0)', 'a = 1'});
        end
        function structFieldLHSPassesThrough(tc)
            d = gdsge.parser.parseVarDecls('opts.maxit = 10;');
            tc.verifyEqual(d.setupStmts, {'opts.maxit = 10'});
        end
        function controlFlowPassesThroughInOrder(tc)
            declText = sprintf('kgrid = zeros(1,3);\nfor i = 1:3\nkgrid(i) = i^2;\nend');
            d = gdsge.parser.parseVarDecls(declText);
            tc.verifyEqual(d.setupStmts, ...
                {'kgrid = zeros(1,3)', 'for i = 1:3', 'kgrid(i) = i^2', 'end'});
        end
        function typoKeywordErrors(tc)
            tc.verifyError(@() gdsge.parser.parseVarDecls('paramters beta gamma;'), ...
                'gdsge:parser:probableTypo');
        end
        function capitalizedKeywordErrors(tc)
            % Pass A keyword match is case-sensitive; 'Parameters' is distance 1
            tc.verifyError(@() gdsge.parser.parseVarDecls('Parameters beta;'), ...
                'gdsge:parser:probableTypo');
        end
        function typoMessageNamesSuggestion(tc)
            try
                gdsge.parser.parseVarDecls('var_stat w;');
                tc.verifyFail('expected probableTypo error');
            catch e
                tc.verifyEqual(e.identifier, 'gdsge:parser:probableTypo');
                tc.verifyTrue(contains(e.message, 'var_stat'));
                tc.verifyTrue(contains(e.message, 'var_state'));
            end
        end
```

- [ ] **Step 2: Run the test file to verify the new tests fail**

Run (from repo root `D:\refactor_gdsge`):
```
matlab -batch "addpath('src'); results = runtests('tests/parser/tParseVarDecls.m'); disp(table(results)); assertSuccess(results)"
```
Expected: the nine new tests FAIL (passthrough tests with `gdsge:parser:unknownDeclaration`, typo tests because the error id differs); pre-existing tests still pass. Exit code nonzero.

- [ ] **Step 3: Implement Pass B changes**

In `src/+gdsge/+parser/parseVarDecls.m`, replace the Pass B loop (currently lines 58-78) with:

```matlab
% ---- Pass B: classify residual statements ---------------------------------
% Statements that are not keyword declarations. A single `name = rhs`
% assignment is routed by its LHS (interp update / tensor assign / setup +
% state grid). Any other statement form — multi-output [a,b] = f(...), bare
% function calls, struct-field LHS, control-flow fragments — is MATLAB
% passthrough: it joins setupStmts verbatim, runs in the parse-time eval, and
% is replayed in the generated files (old-toolbox behavior). One hard error:
% a command-form statement whose first word is a near-miss of a declaration
% keyword is almost certainly a typo'd declaration; fail early with a
% suggestion instead of obscurely at eval.
for j = 1:numel(assignIdx)
    st = stmts{assignIdx(j)};
    lhs = assignLHS(st);
    if isempty(lhs)
        checkKeywordTypo(st, DECL_KINDS);
        decls.setupStmts{end+1} = st; %#ok<AGROW>
        continue;
    end
    rhs = strtrim(assignRHS(st));
    if ismember(lhs, decls.interpNames)
        decls.interpUpdate{end+1} = struct('name', lhs, 'expr', rhs); %#ok<AGROW>
    elseif ismember(lhs, decls.tensorNames)
        % var_tensor assignment: capture as opaque text (codegen emits it with
        % ndgrid broadcasting). Plain eval would dimension-mismatch.
        decls.tensorAssign{end+1} = struct('name', lhs, 'expr', rhs); %#ok<AGROW>
    else
        decls.setupStmts{end+1} = st; %#ok<AGROW>
        if ismember(lhs, decls.stateNames)
            decls.gridText{end+1} = struct('name', lhs, 'expr', rhs); %#ok<AGROW>
        end
    end
end
```

Then add two local functions at the end of the file (before or after the existing locals, matching their style):

```matlab
function checkKeywordTypo(st, DECL_KINDS)
% Guard against typo'd declaration keywords. Fires only on command-form
% statements — "<word> <word> ..." with no '=' on the first line, the shape a
% declaration takes — whose first word is within edit distance 2 of a keyword.
% Function-call syntax ("rouwen(...)"), assignments, control-flow headers
% (contain '='), and single-word statements ("end") never match.
l = firstLine(st);
if any(l == '='); return; end
m = regexp(l, '^\s*([A-Za-z_]\w*)\s+[A-Za-z_]', 'tokens', 'once');
if isempty(m); return; end
word = m{1};
KEYWORDS = [DECL_KINDS, {'inbound', 'inbound_init', 'initial'}];
for k = 1:numel(KEYWORDS)
    if levenshtein(word, KEYWORDS{k}) <= 2
        error('gdsge:parser:probableTypo', ...
            ['Unknown keyword "%s" in "%s" — did you mean "%s"? (To pass a ' ...
             'command-form MATLAB statement through to the generated code, ' ...
             'use function syntax: %s(...).)'], ...
            word, strtrim(l), KEYWORDS{k}, word);
    end
end
end

function d = levenshtein(a, b)
% Plain Levenshtein distance. The caller only tests <= 2, so bail out early
% when the length gap alone exceeds that.
la = numel(a); lb = numel(b);
if abs(la - lb) > 2; d = 3; return; end
D = zeros(la + 1, lb + 1);
D(:, 1) = (0:la)';
D(1, :) = 0:lb;
for i = 1:la
    for j = 1:lb
        D(i+1, j+1) = min([D(i, j+1) + 1, D(i+1, j) + 1, D(i, j) + (a(i) ~= b(j))]);
    end
end
d = D(la + 1, lb + 1);
end
```

Also update the file's header comment (lines 2-5): change "and the residual setup statements (param/shock/grid value assignments) to be eval'd later" to "and the residual setup statements (param/shock/grid value assignments plus verbatim MATLAB passthrough) to be eval'd later".

- [ ] **Step 4: Run the test file to verify it passes**

Run:
```
matlab -batch "addpath('src'); results = runtests('tests/parser/tParseVarDecls.m'); disp(table(results)); assertSuccess(results)"
```
Expected: ALL tests pass (new ones and pre-existing ones — including `reopenedDeclKindErrors`, which Task 2 changes but Task 1 must not break). Exit code 0.

- [ ] **Step 5: Commit**

```
git add src/+gdsge/+parser/parseVarDecls.m tests/parser/tParseVarDecls.m
git commit -m "feat(parser): preamble MATLAB passthrough with typo guard

Pass B no longer rejects statements that aren't single-identifier
assignments: multi-output [a,b]=f(...), bare calls, struct-field LHS,
and control flow pass through as setup text, as the old toolbox did.
New gdsge:parser:probableTypo error catches command-form statements
whose first word is a near-miss of a declaration keyword."
```

---

### Task 2: Pass C — allow re-opened declaration blocks

**Files:**
- Modify: `src/+gdsge/+parser/parseVarDecls.m` (Pass C, lines ~80-133)
- Test: `tests/parser/tParseVarDecls.m`

**Interfaces:**
- Consumes: Task 1's parseVarDecls state (Pass C is independent of Pass B's change).
- Produces: `decls.sections` may now contain repeated `kind` values, in source order. Error id `gdsge:parser:reopenedDeclBlock` is removed from the codebase. The section-bodies-concatenate-to-setupStmts contract (verified by existing test `sectionBodiesConcatToSetupText`) still holds.

- [ ] **Step 1: Write the failing tests**

In `tests/parser/tParseVarDecls.m`, **replace** the `reopenedDeclKindErrors` method (the one asserting `gdsge:parser:reopenedDeclBlock`) with:

```matlab
        function reopenedDeclKindAllowed(tc)
            % old-toolbox layouts interleave declaration blocks freely
            d = gdsge.parser.parseVarDecls(sprintf('parameters a;\na = 1;\nparameters b;\nb = 2;'));
            kinds = cellfun(@(s) s.kind, d.sections, 'UniformOutput', false);
            tc.verifyEqual(kinds, {'parameters', 'parameters'});
            tc.verifyEqual(strtrim(d.sections{1}.body), 'a = 1;');
            tc.verifyEqual(strtrim(d.sections{2}.body), 'b = 2;');
        end
        function reopenedKindsPreserveSourceOrder(tc)
            % concatenated section bodies must equal the flat setup text even
            % with a kind re-opened after another block (emitSetup's contract)
            declText = sprintf([ ...
                'parameters a;\na = 1;\n', ...
                'var_state w;\nw = linspace(0,1,3);\n', ...
                'parameters b;\nb = a + 1;\n', ...
                'var_policy c;\ninbound c 0 1;']);
            d = gdsge.parser.parseVarDecls(declText);
            kinds = cellfun(@(s) s.kind, d.sections, 'UniformOutput', false);
            tc.verifyEqual(kinds, {'parameters', 'var_state', 'parameters', 'var_policy'});
            bodies = cellfun(@(s) s.body, d.sections, 'UniformOutput', false);
            bodies = bodies(~cellfun(@isempty, bodies));
            tc.verifyEqual(strjoin(bodies, sprintf('\n')), ...
                [strjoin(d.setupStmts, sprintf(';\n')), ';']);
        end
```

- [ ] **Step 2: Run the test file to verify the new tests fail**

Run:
```
matlab -batch "addpath('src'); results = runtests('tests/parser/tParseVarDecls.m'); disp(table(results)); assertSuccess(results)"
```
Expected: the two new tests FAIL with `gdsge:parser:reopenedDeclBlock`; everything else passes.

- [ ] **Step 3: Implement the Pass C change**

In `src/+gdsge/+parser/parseVarDecls.m` Pass C: delete the `seenKinds` variable and the re-open error. The loop head becomes:

```matlab
sections = {};
curKind = 'global';
curBody = {};
headerOpen = false;
for i = 1:numel(stmts)
    st = stmts{i};
    kw = regexp(st, '^([A-Za-z_]\w*)', 'tokens', 'once');
    key = ''; if ~isempty(kw); key = kw{1}; end
    if ismember(key, DECL_KINDS)
        if headerOpen && strcmp(key, curKind)
            continue;   % contiguous same-kind declaration line — extend header
        end
        % flush the current section (skip an empty leading 'global' preamble).
        % A kind may re-open later (old-toolbox layouts interleave declaration
        % blocks); sections are positional and kinds may repeat in ir.setup.
        if ~(strcmp(curKind, 'global') && isempty(curBody))
            sections{end+1} = struct('kind', curKind, 'body', joinBody(curBody)); %#ok<AGROW>
        end
        curBody = {};
        curKind = key; headerOpen = true;
    else
```

(the `else` branch and everything after the loop are unchanged.)

- [ ] **Step 4: Run the test file to verify it passes**

Run:
```
matlab -batch "addpath('src'); results = runtests('tests/parser/tParseVarDecls.m'); disp(table(results)); assertSuccess(results)"
```
Expected: ALL tests pass, exit code 0.

- [ ] **Step 5: Verify no other code references the removed error id**

Run: `grep -rn "reopenedDeclBlock" src/ tests/ docs/`
Expected: no hits in `src/` or `tests/`. (A hit in `docs/superpowers/specs/` or plans is fine — historical documents.)

- [ ] **Step 6: Commit**

```
git add src/+gdsge/+parser/parseVarDecls.m tests/parser/tParseVarDecls.m
git commit -m "feat(parser): allow re-opened declaration blocks

Old-toolbox gmod files interleave declaration kinds (parameters ...
var_shock ... parameters ...). Sections are positional and ir.setup
may repeat kinds; emitSetup already replays bodies in order."
```

---

### Task 3: Integration — passthrough through eval, front end, and IR round-trip

**Files:**
- Create: `tests/parser/tSetupPassthrough.m`

**Interfaces:**
- Consumes: `gdsge.parser.parseDeclarations(declText)` → struct with `params{i}.value` (numeric), `shocks.values.<name>` (numeric); `gdsge.parser.parseFrontEnd(gmodText, modelName)` → IR with `setup` (list of `{kind, body}`); `gdsge.ir.roundtrip(ir)` / `gdsge.ir.isequalIR(a, b)`.
- Produces: nothing new — this task locks the end-to-end parser behavior with tests only.

- [ ] **Step 1: Write the tests (expected to pass already — this is the integration gate)**

Create `tests/parser/tSetupPassthrough.m`:

```matlab
classdef tSetupPassthrough < matlab.unittest.TestCase
    % Preamble MATLAB passthrough end to end: multi-output assignments and
    % control flow run in parseDeclarations' parse-time eval, land in ir.setup
    % in source order (with re-opened kinds), and survive the JSON round-trip.
    methods (Test)
        function multiOutputAndLoopEvalToValues(tc)
            declText = sprintf([ ...
                'parameters a b;\n', ...
                '[a, b] = deal(2, 3);\n', ...
                'var_shock z;\nshock_num = 2;\n', ...
                'zg = zeros(1, 2);\n', ...
                'for i = 1:2\nzg(i) = i/10;\nend\n', ...
                'z = zg;\n', ...
                'shock_trans = [0.5 0.5; 0.5 0.5];\n', ...
                'var_state k;\nk = linspace(0.5, 1.5, 3);']);
            out = gdsge.parser.parseDeclarations(declText);
            tc.verifyEqual(out.params{1}.value, 2);
            tc.verifyEqual(out.params{2}.value, 3);
            tc.verifyEqual(out.shocks.values.z, [0.1 0.2]);
            tc.verifyEqual(out.shocks.count, 2);
        end
        function frontEndCarriesPassthroughSetup(tc)
            ir = gdsge.parser.parseFrontEnd(gmodText(), 'tinypass');
            kinds = cellfun(@(s) s.kind, ir.setup, 'UniformOutput', false);
            tc.verifyEqual(kinds, {'parameters', 'var_shock', 'var_state', ...
                'parameters', 'var_policy', 'var_interp'});
            allBody = strjoin(cellfun(@(s) s.body, ir.setup, ...
                'UniformOutput', false), sprintf('\n'));
            tc.verifyTrue(contains(allBody, '[zg1, zg2] = ndgrid'));
            tc.verifyEqual(ir.params{2}.value, 2);          % gamma, re-opened block
            tc.verifyEqual(ir.shocks.values.z, [0.9 1.1]);  % via ndgrid passthrough
            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));
        end
    end
end

function txt = gmodText()
txt = sprintf([ ...
    'parameters beta;\n' ...
    'beta = 0.95;\n' ...
    'var_shock z;\n' ...
    'shock_num = 2;\n' ...
    '[zg1, zg2] = ndgrid([0.9 1.1], 1);\n' ...
    'z = zg1(:)'';\n' ...
    'shock_trans = [0.5 0.5; 0.5 0.5];\n' ...
    'var_state k;\n' ...
    'k = linspace(0.5, 1.5, 3);\n' ...
    'parameters gamma;\n' ...
    'gamma = 2;\n' ...
    'var_policy c;\n' ...
    'inbound c 0 2;\n' ...
    'var_interp vf;\n' ...
    'initial vf k;\n' ...
    'vf = q;\n' ...
    'model;\n' ...
    '  q = beta*c + gamma*0 + z*0;\n' ...
    '  equations;\n' ...
    '    c - k;\n' ...
    '  end;\n' ...
    'end;\n' ...
    'simulate;\n' ...
    '  num_periods = 5;\n' ...
    '  num_samples = 1;\n' ...
    '  initial k 1;\n' ...
    '  initial shock 1;\n' ...
    '  var_simu c;\n' ...
    '  k'' = c;\n' ...
    'end;\n']);
end
```

- [ ] **Step 2: Run the test file**

Run:
```
matlab -batch "addpath('src'); results = runtests('tests/parser/tSetupPassthrough.m'); disp(table(results)); assertSuccess(results)"
```
Expected: PASS. If `frontEndCarriesPassthroughSetup` fails on the exact `kinds` list, print `cellfun(@(s) s.kind, ir.setup, 'UniformOutput', false)` interactively, check whether the difference is a genuine ordering bug or an off-by-one in the expected list (e.g. an extra empty trailing section), fix whichever is wrong, and re-run. Do not weaken the body/order assertions.

- [ ] **Step 3: Commit**

```
git add tests/parser/tSetupPassthrough.m
git commit -m "test(parser): integration gate for preamble passthrough

Multi-output + loop preamble evals to IR values; passthrough text
reaches ir.setup in source order with re-opened kinds; JSON round-trip
preserves it."
```

---

### Task 4: Docs, changelog, full regression suite

**Files:**
- Modify: `docs/notes-for-agents.md` (the "gmod language is a MATLAB superset" section, around lines 51-54)
- Modify: `Changelog.md` (add an Unreleased section above `## [0.2.1]`)

**Interfaces:**
- Consumes: Tasks 1-3 merged behavior.
- Produces: documentation only.

- [ ] **Step 1: Update notes-for-agents.md**

In the "gmod language is a MATLAB superset" area of `docs/notes-for-agents.md`, append this paragraph:

```markdown
The declaration region accepts *arbitrary* MATLAB, not just `name = value`
assignments: multi-output calls (`[a,b] = ndgrid(...)`), bare calls (`rng(0)`),
struct-field LHS, and `for`/`if` blocks all pass through verbatim
(`parseVarDecls` Pass B), and declaration kinds may re-open
(`parameters ... var_shock ... parameters ...`, Pass C — `ir.setup` kinds can
repeat). Passthrough code runs twice: once in the parse-time eval and again in
the generated iter/simulate files after the default flags. Two guardrails
replace the old toolbox's blind copying: a command-form statement whose first
word is edit-distance ≤ 2 from a declaration keyword raises
`gdsge:parser:probableTypo`, and broken MATLAB still fails as
`gdsge:parser:setupEvalFailed` with numbered lines. Limitation: state grids and
interp updates are only captured from single `name = rhs` assignments.
```

- [ ] **Step 2: Add the changelog entry**

In `Changelog.md`, insert directly below the intro paragraph (above `## [0.2.1]`):

```markdown
## [Unreleased]

### Added

- The gmod declaration region accepts arbitrary MATLAB again, as the pre-refactor
  toolbox did: multi-output assignments (`[a,b] = ndgrid(...)`), bare function
  calls, struct-field assignments, `for`/`if` blocks, and declaration blocks that
  re-open (`parameters ... var_shock ... parameters ...`). A new guard catches
  misspelled declaration keywords (`paramters`) at parse time with a suggestion
  instead of a confusing downstream error.
```

- [ ] **Step 3: Run the full regression suite**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: exit code 0; check `tests/results/junit.xml` for 0 failures/errors. Existing model IR snapshots and codegen goldens must be untouched (`git status` shows no modified snapshot/golden files — only the files this plan edits).

- [ ] **Step 4: Commit**

```
git add docs/notes-for-agents.md Changelog.md
git commit -m "docs: preamble MATLAB passthrough notes + changelog"
```

---

### Task 5: End-to-end acceptance on incomp.gmod (manual, outside the repo)

**Files:**
- Create: `scratch/incomp_e2e/` (git-ignored working copy — nothing from it is committed)
- Read-only inputs: `D:\macro_local\exchange_rate_dynamics_debug\{incomp.gmod, rouwen.m, IterRslt_incomp.mat}` — never modify or overwrite these.

**Interfaces:**
- Consumes: the merged parser from Tasks 1-2; `gdsge_codegen`, `iter_incomp` from the new toolbox.
- Produces: a written acceptance report in the session (not a repo file): parse/codegen outcome, first-iteration Metric, closeness verdict.

- [ ] **Step 1: Locate the missing helper**

`bivariateNormalProbMidpoints` is called by the gmod preamble but is not in the model folder. Find it:
```
Glob pattern "**/bivariateNormalProbMidpoints.m" under D:\macro_local
```
If not found there, search the user's MATLAB path candidates (e.g. `D:\` one level deep). If it cannot be found, STOP and ask the user where it lives — do not stub it.

- [ ] **Step 2: Stage a scratch copy**

```
mkdir -p scratch/incomp_e2e
cp D:/macro_local/exchange_rate_dynamics_debug/incomp.gmod scratch/incomp_e2e/
cp D:/macro_local/exchange_rate_dynamics_debug/rouwen.m scratch/incomp_e2e/
cp D:/macro_local/exchange_rate_dynamics_debug/IterRslt_incomp.mat scratch/incomp_e2e/
cp <path-found-in-step-1>/bivariateNormalProbMidpoints.m scratch/incomp_e2e/
```

- [ ] **Step 3: Codegen with the new toolbox**

```
matlab -batch "addpath('D:/refactor_gdsge/src'); cd('D:/refactor_gdsge/scratch/incomp_e2e'); gdsge_codegen('incomp')"
```
Expected: parser accepts the file (no `unknownDeclaration` / `reopenedDeclBlock`; at most `setupBlockMismatch` warnings, e.g. for `b_last` assigned before its `var_state` declaration) and C++ compilation succeeds. This is a big model — compilation may take minutes. If codegen fails downstream of the parser, STOP: report the failure to the user and treat it as its own systematic-debugging + test cycle (per spec §6), not an inline patch.

- [ ] **Step 4: Warm-start iteration (first-iteration closeness — NOT a reproduction run)**

Write `scratch/incomp_e2e/run_e2e.m`:

```matlab
warmUp = load('IterRslt_incomp.mat');
options = struct;
options.WarmUp.asg_interp_struct = warmUp.IterRslt.asg_interp_struct;
options.SkipModelInit = 1;
options.PrintFreq = 1;
options.ADAPTIVE_FACTOR = 1.2;
options.AsgMaxLevel = 3;
options.AsgThreshold = 1e-2;
options.AsgOutputMaxLevel = 7;
options.TolEq = 1e-5;
IterRslt = iter_incomp(options);
fprintf('FINAL Metric: %g after %d iters\n', IterRslt.Metric, IterRslt.Iter);
save('IterRslt_new.mat', 'IterRslt');
exit(double(~(isfinite(IterRslt.Metric) && IterRslt.Metric < 1e-4)));
```

Run:
```
matlab -batch "addpath('D:/refactor_gdsge/src'); cd('D:/refactor_gdsge/scratch/incomp_e2e'); run_e2e"
```
Expected: iterations solve cleanly from the warm start and the Metric is small from the first printed iteration (the warm-up is a converged solution; per `main_incomp.m` this is exactly the user's workflow). Success gate: finite Metric < 1e-4. If `iter_incomp` errors or the first Metric is large (order 1), STOP and report — likely the cross-toolbox ASG warm-start path (spec §6), a separate fix cycle.

- [ ] **Step 5: Report**

Report to the user: helper location, codegen wall time, per-iteration Metric trajectory, final Metric/Iter, and any warnings emitted. List any downstream gaps hit (each becomes its own mini design/fix loop). Do not commit anything from `scratch/`.

---

## Self-review notes

- **Spec coverage:** §3.1 Pass B → Task 1; §3.2 Pass C → Task 2; §3.3 backstop → unchanged (existing `tEvalSetup` covers `setupEvalFailed`); §4 semantics → Task 4 docs; §5 tests → Tasks 1-4 (synthetic only, suite regression in Task 4); §5 end-to-end acceptance + §6 risks → Task 5. `GNDSGE_INTERP_VEC` risk from spec §6 is already handled by `preprocess.m:33` (`\<GNDSGE` → `GDSGE` rewrite) — verified during planning, no task needed.
- **Existing tests changed:** `unknownStatementErrors` (replaced in Task 1) and `reopenedDeclKindErrors` (replaced in Task 2) are deliberate behavior flips mandated by the spec; all other `tParseVarDecls` tests must keep passing unmodified.
- **Type consistency:** error ids `gdsge:parser:probableTypo` (new), `gdsge:parser:unknownDeclaration` (kept for sized-name tokens only), `gdsge:parser:reopenedDeclBlock` (removed) — consistent across tasks.
