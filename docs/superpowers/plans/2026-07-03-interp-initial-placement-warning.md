# Exempt interp `initial`/update placement from setupBlockMismatch — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop `gdsge:parser:setupBlockMismatch` warnings from firing on interp `initial X expr` statements and interp-name assignments, which the canonical GDSGE layout forces to appear after `model_init` (positionally inside the `var_aux_init` section).

**Architecture:** One function changes — `expectedSection()` in `src/+gdsge/+parser/parseVarDecls.m` stops claiming `var_interp` as the expected section for those two statement forms, joining `inbound`/`inbound_init` in the not-warned category. Extraction of `interpInitial`/`interpUpdate` (Pass A/B) is untouched; no IR, golden, or codegen output changes.

**Tech Stack:** MATLAB R2025b (`C:\Program Files\MATLAB\R2025b\bin\matlab.exe`), driven headless via `matlab -batch`. Tests: MATLAB unittest (`tests/parser/tParseVarDecls.m`).

**Spec:** `docs/superpowers/specs/2026-07-03-interp-initial-placement-warning-design.md`

## Global Constraints

- Never add a toolbox source to a persistent MATLAB path; each `matlab -batch` process `addpath`s exactly one source (`src/` here).
- PowerShell 7 is not installed — run the suite via `matlab -batch "cd('tests'); run_tests"`, not `tests/run.ps1`.
- Never run two MATLAB processes concurrently.
- `tests/results/junit.xml` + exit code are authoritative (results.tap accumulates stale lines).
- Backward compatible: no IR schema, golden snapshot, or generated-code changes are expected — if any golden diff appears, that's a bug in the change.

---

### Task 1: Parser change — no warning for interp `initial`/update placement

**Files:**
- Modify: `src/+gdsge/+parser/parseVarDecls.m:120-127` (warning branch) and `:211-225` (`expectedSection`)
- Test: `tests/parser/tParseVarDecls.m` (add one test method to the existing class)

**Interfaces:**
- Consumes: `gdsge.parser.parseVarDecls(declText)` — existing; returns struct with `interpInitial` (cell of `struct('name',...,'expr',...)`), `interpUpdate` (same shape), `sections`.
- Produces: no new interfaces — behavior change only (the two statement forms no longer produce `gdsge:parser:setupBlockMismatch`).

- [ ] **Step 1: Write the failing test**

Add this method inside the `methods (Test)` block of `tests/parser/tParseVarDecls.m`, after `warningMessageIsActionable` (keep the file's existing `lastwarn`/`warning('off','all')` idiom):

```matlab
function interpInitialAfterAuxInitDoesNotWarn(tc)
    % canonical GDSGE layout: model_init is excised before this pass, so
    % initial/interp updates land positionally in the var_aux_init section
    declText = sprintf([ ...
        'var_interp Ev;\n', ...
        'var_aux_init E;\n', ...
        'initial Ev E;\nEv = E;']);
    [~] = lastwarn('');
    ws0 = warning('off','all'); restore = onCleanup(@() warning(ws0)); %#ok<NASGU>
    d = gdsge.parser.parseVarDecls(declText);
    [~, id] = lastwarn();
    tc.verifyNotEqual(id, 'gdsge:parser:setupBlockMismatch');
    % extraction must be unaffected by placement
    tc.verifyEqual(d.interpInitial{1}.name, 'Ev');
    tc.verifyEqual(d.interpInitial{1}.expr, 'E');
    tc.verifyEqual(d.interpUpdate{1}.name, 'Ev');
    tc.verifyEqual(d.interpUpdate{1}.expr, 'E');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run (from repo root `D:\refactor_gdsge`):

```
matlab -batch "addpath('src'); results = runtests('tests/parser/tParseVarDecls.m'); disp(table(results)); assert(all([results.Passed]))"
```

Expected: FAIL — `interpInitialAfterAuxInitDoesNotWarn` fails on `verifyNotEqual(id, 'gdsge:parser:setupBlockMismatch')` because both the `initial Ev E` statement and the `Ev = E` assignment currently warn. All other tests pass.

- [ ] **Step 3: Implement the change**

In `src/+gdsge/+parser/parseVarDecls.m`, make two edits.

**Edit 1** — `expectedSection` (currently lines 211–225). Replace the whole function with:

```matlab
function kind = expectedSection(st, key, decls)
% The section a setup/associate statement *belongs* to, by its LHS. '' = no
% expectation (options, helper intermediates, declaration lines, inbound,
% interp initial/updates).
kind = '';
% 'initial' and interp-name assignments are extracted by keyword/LHS name
% regardless of position (never enter a section body), and the canonical
% layout forces them after model_init -> positionally in var_aux_init.
if ismember(key, {'initial','inbound','inbound_init'}); return; end
lhs = assignLHS(st);
if isempty(lhs); return; end
if any(strcmp(lhs, {'shock_trans','shock_num'}));  kind = 'var_shock'; return; end
if ismember(lhs, decls.shockNames);  kind = 'var_shock';  return; end
if ismember(lhs, decls.stateNames);  kind = 'var_state';  return; end
if ismember(lhs, decls.interpNames); return; end
if ismember(lhs, decls.tensorNames); kind = 'var_tensor'; return; end
if ismember(lhs, decls.paramNames);  kind = 'parameters'; return; end
end
```

(The `interpNames` line stays as an explicit early return — a declared interp name must not fall through to the `paramNames` check.)

**Edit 2** — warning branch in Pass C (currently lines 120–127). Every path that now sets a non-empty `kind` requires a non-empty LHS, so the `initial` fallback is dead code. Replace:

```matlab
        if ~isempty(expected) && ~strcmp(expected, curKind)
            lhs = assignLHS(st);
            if isempty(lhs) && strcmp(key, 'initial'); lhs = 'initial'; end
            warning('gdsge:parser:setupBlockMismatch', ...
```

with:

```matlab
        if ~isempty(expected) && ~strcmp(expected, curKind)
            lhs = assignLHS(st);
            warning('gdsge:parser:setupBlockMismatch', ...
```

- [ ] **Step 4: Run the test file to verify it passes**

Run:

```
matlab -batch "addpath('src'); results = runtests('tests/parser/tParseVarDecls.m'); disp(table(results)); assert(all([results.Passed]))"
```

Expected: PASS — all tests including the new one and the existing warning tests (`shockTransAfterStateWarns`, `shockNumBeforeShockWarns`, `warningMessageIsActionable`, `canonicalOrderDoesNotWarn`).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/parseVarDecls.m tests/parser/tParseVarDecls.m
git commit -m "fix(parser): don't warn on interp initial/update placement

The canonical GDSGE layout puts 'initial X expr' and interp-name
assignments after model_init (their RHS references var_aux_init
results), which lands them positionally in the var_aux_init section
once splitBlocks excises model_init. These statements are extracted by
keyword/LHS name regardless of position and never enter a section
body, so placement is functionally inert -- same rationale as the
existing inbound/inbound_init exemption.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01YCcAffmffoVzww6W1CEmSW"
```

---

### Task 2: Amend the sections design doc + full-suite and real-model verification

**Files:**
- Modify: `docs/superpowers/specs/2026-06-16-decl-region-positional-sections-design.md:81-95` (expected-section table + not-warned note)

**Interfaces:**
- Consumes: Task 1's behavior change (needed for the real-model verification step).
- Produces: nothing — documentation and verification only.

- [ ] **Step 1: Update the expected-section table**

In `docs/superpowers/specs/2026-06-16-decl-region-positional-sections-design.md`, delete this table row (currently line 86):

```markdown
| a declared `var_interp` name, or `initial X` | `var_interp`     |
```

and replace the inbound bullet (currently lines 94–95):

```markdown
- `inbound` / `inbound_init` **placement is not warned** — low signal; they are already
  structurally bound to policy/aux bounds regardless of section.
```

with:

```markdown
- `inbound` / `inbound_init` **placement is not warned** — low signal; they are already
  structurally bound to policy/aux bounds regardless of section.
- Interp `initial X` statements and declared-interp-name assignments are likewise **not
  warned** (amended 2026-07-03, see
  `2026-07-03-interp-initial-placement-warning-design.md`): they are extracted by
  keyword/LHS name regardless of position, and the canonical layout forces them after
  `model_init`, which positionally lands them in the `var_aux_init` section.
```

- [ ] **Step 2: Run the full test suite**

Run (from repo root):

```
matlab -batch "cd('tests'); run_tests"
```

Expected: exit code 0, all tests pass (check `tests/results/junit.xml` if in doubt — the TAP file accumulates stale lines). Takes a few minutes. Do not run any other MATLAB process concurrently.

- [ ] **Step 3: Verify on the real model that motivated the change**

Run (from repo root; parses the user's gmod through the front end and asserts no mismatch warning fired):

```
matlab -batch "addpath('D:/refactor_gdsge/src'); cd('D:/macro_local/debug_gdsge_re'); lastwarn(''); txt = fileread('large13calibrationv5utilitiesCCPsREEtweak.gmod'); ir = gdsge.parser.parseFrontEnd(txt, 'large13calibrationv5utilitiesCCPsREEtweak'); [~, id] = lastwarn(); assert(~strcmp(id, 'gdsge:parser:setupBlockMismatch'), 'setupBlockMismatch still fired'); disp('OK: no setupBlockMismatch')"
```

Expected: prints `OK: no setupBlockMismatch`, exit code 0, and no `setupBlockMismatch` warning text anywhere in the output. (Other warnings from this model, if any, are out of scope — only `setupBlockMismatch` matters here.)

- [ ] **Step 4: Commit**

```bash
git add docs/superpowers/specs/2026-06-16-decl-region-positional-sections-design.md
git commit -m "docs(specs): record interp initial/update placement exemption

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01YCcAffmffoVzww6W1CEmSW"
```
