# Preloaded-pyenv mismatch handling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the SymPy backend work (or fail with an actionable message) when other code has already loaded a Python interpreter into the MATLAB session before GDSGE runs.

**Architecture:** `ensurePyenv` gets explicit handling of the `Loaded` pyenv state — probe by importing, terminate-and-repoint when the wrong interpreter is OutOfProcess, warn clearly when it is InProcess — with all side effects behind an injectable `port` struct so every branch is unit-testable without touching real pyenv state. `gdsge_setup_sympy` validates the *session* (not just the venv) at the end, and `resolveBackend`/`generateCxx` messages stop saying "Python not detected" when the real problem is elsewhere.

**Tech Stack:** MATLAB R2025b (`C:\Program Files\MATLAB\R2025b\bin\matlab.exe`) headless via `matlab -batch`; MATLAB unittest; `pyenv`/`terminate` APIs; uv-managed venv at `pyext/.venv`.

**Spec:** `docs/superpowers/specs/2026-07-03-pyenv-preloaded-mismatch-design.md`

## Global Constraints

- Never add a toolbox source to a persistent MATLAB path; each `matlab -batch` process `addpath`s exactly one source (`src/` here).
- Run the suite via `matlab -batch "cd('tests'); run_tests"` (PowerShell 7 is not installed; `tests/run.ps1` is unusable).
- Never run two MATLAB processes concurrently.
- The default (autodiff) path must keep needing no Python: `ensurePyenv` is only ever called on the sympy path and in setup.
- Existing call sites use single-output `ensurePyenv()` (`resolveBackend.m:13`, `generateCxx.m:42`, `tests/+gdsgetest/sympyAvailable.m:5`) and must keep working unchanged in Task 1.
- MATLAB facts this plan relies on (verified by local probes on 2026-07-03): `pyenv` `Version`/`ExecutionMode` are settable only while `Status` is `NotLoaded` or `Terminated`; an InProcess interpreter cannot be unloaded without restarting MATLAB; `terminate(pyenv)` on an OutOfProcess interpreter moves `Status` to `Terminated`, after which re-pointing works.

---

### Task 1: `ensurePyenv` — explicit Loaded-state handling behind an injectable port

**Files:**
- Modify: `src/+gdsge/+codegen/+sympy/ensurePyenv.m` (full rewrite, same public name)
- Create: `tests/codegen/tEnsurePyenv.m`

**Interfaces:**
- Consumes: `gdsge.codegen.sympy.venvPython(pyextDir)` — existing; returns absolute interpreter path or `''`.
- Produces: `[ok, why] = gdsge.codegen.sympy.ensurePyenv(port)`. `ok` logical; `why` char, `''` on success. `port` optional struct with function-handle fields `venvPython` (`@() path-or-''`), `getEnv` (`@() struct/object with Status, ExecutionMode, Executable`), `setEnv` (`@(venvPy) ...`), `terminate` (`@() ...`), `tryImports` (`@() logical`). Omitted → real implementations. Warning id `gdsge:sympy:inProcessMismatch`. Tasks 2 and 3 call the two-output form.

- [ ] **Step 1: Write the failing tests**

Create `tests/codegen/tEnsurePyenv.m`:

```matlab
classdef tEnsurePyenv < matlab.unittest.TestCase
    % ensurePyenv branch coverage via an injected effect port — no real
    % pyenv/Python state is touched. See the 2026-07-03 pyenv-preloaded-
    % mismatch spec: a foreign interpreter already Loaded in the session
    % must be probed, replaced (OutOfProcess) or explained (InProcess).
    methods (TestClassSetup)
        function srcPath(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src')));
        end
    end
    methods (Test)
        function venvMissingFails(tc)
            [port, log] = fakePort('', "NotLoaded", "OutOfProcess", true);
            [ok, why] = gdsge.codegen.sympy.ensurePyenv(port);
            tc.verifyFalse(ok);
            tc.verifyTrue(contains(why, 'gdsge_setup_sympy'));
            tc.verifyEmpty(log('calls'));            % nothing touched
        end
        function notLoadedSetsVenvThenImports(tc)
            [port, log] = fakePort('C:\fake\venv\python.exe', "NotLoaded", "OutOfProcess", true);
            [ok, why] = gdsge.codegen.sympy.ensurePyenv(port);
            tc.verifyTrue(ok);
            tc.verifyEmpty(why);
            tc.verifyEqual(log('calls'), {'setEnv:C:\fake\venv\python.exe', 'tryImports'});
        end
        function terminatedCountsAsSettable(tc)
            [port, log] = fakePort('C:\fake\venv\python.exe', "Terminated", "OutOfProcess", true);
            ok = gdsge.codegen.sympy.ensurePyenv(port);
            tc.verifyTrue(ok);
            tc.verifyEqual(log('calls'), {'setEnv:C:\fake\venv\python.exe', 'tryImports'});
        end
        function loadedWithWorkingImportsUsedAsIs(tc)
            % a foreign-but-SymPy-capable interpreter is acceptable: the
            % bridge only exchanges JSON strings
            [port, log] = fakePort('C:\fake\venv\python.exe', "Loaded", "InProcess", true);
            ok = gdsge.codegen.sympy.ensurePyenv(port);
            tc.verifyTrue(ok);
            tc.verifyEqual(log('calls'), {'tryImports'});   % no setEnv, no terminate
        end
        function loadedOutOfProcessMismatchTerminatesAndRepoints(tc)
            [port, log] = fakePort('C:\fake\venv\python.exe', "Loaded", "OutOfProcess", [false true]);
            [ok, why] = gdsge.codegen.sympy.ensurePyenv(port);
            tc.verifyTrue(ok);
            tc.verifyEmpty(why);
            tc.verifyEqual(log('calls'), ...
                {'tryImports', 'terminate', 'setEnv:C:\fake\venv\python.exe', 'tryImports'});
        end
        function loadedInProcessMismatchWarnsAndFails(tc)
            [port, log] = fakePort('C:\fake\venv\python.exe', "Loaded", "InProcess", false);
            [ok, why] = tc.verifyWarning( ...
                @() gdsge.codegen.sympy.ensurePyenv(port), 'gdsge:sympy:inProcessMismatch');
            tc.verifyFalse(ok);
            tc.verifyTrue(contains(why, 'Restart MATLAB'));
            tc.verifyTrue(contains(why, 'C:\fake\loaded\python3.12'));  % names the culprit
            tc.verifyFalse(any(strcmp(log('calls'), 'terminate')));
            tc.verifyFalse(any(startsWith(log('calls'), 'setEnv')));
        end
        function venvImportFailureExplains(tc)
            [port, ~] = fakePort('C:\fake\venv\python.exe', "NotLoaded", "OutOfProcess", false);
            [ok, why] = gdsge.codegen.sympy.ensurePyenv(port);
            tc.verifyFalse(ok);
            tc.verifyTrue(contains(why, 'gdsge_setup_sympy'));
        end
    end
end

% ---- fake effect port ------------------------------------------------------
function [port, log] = fakePort(venvPy, status, mode, importResults)
% importResults: logical scalar/vector consumed one per tryImports call
% (the last value repeats). log is a containers.Map (handle semantics) so
% the closures below can record calls; key 'calls' is the ordered call list.
log = containers.Map();
log('calls') = {};
log('nImports') = 0;
port = struct( ...
    'venvPython', @() venvPy, ...
    'getEnv',     @() struct('Status', status, 'ExecutionMode', mode, ...
                             'Executable', "C:\fake\loaded\python3.12"), ...
    'setEnv',     @(p) record(log, ['setEnv:' char(p)]), ...
    'terminate',  @() record(log, 'terminate'), ...
    'tryImports', @() nextImport(log, importResults));
end

function record(log, name)
c = log('calls'); c{end+1} = name; log('calls') = c;
end

function tf = nextImport(log, importResults)
record(log, 'tryImports');
n = log('nImports') + 1; log('nImports') = n;
tf = importResults(min(n, numel(importResults)));
end
```

- [ ] **Step 2: Run the tests to verify they fail**

Run (from repo root `D:\refactor_gdsge`):

```
matlab -batch "results = runtests('tests/codegen/tEnsurePyenv.m'); disp(table(results)); assert(all([results.Passed]))"
```

(No `addpath` needed — the class applies a PathFixture.)
Expected: FAIL — current `ensurePyenv` accepts no `port` argument ("Too many input arguments") in every test.

- [ ] **Step 3: Rewrite `ensurePyenv.m`**

Replace the entire contents of `src/+gdsge/+codegen/+sympy/ensurePyenv.m` with:

```matlab
function [ok, why] = ensurePyenv(port)
% ENSUREPYENV  Point MATLAB's pyenv at the uv-managed venv interpreter and make
%   the gdsge_sympy package importable. Idempotent; callable from matlab -batch.
%   ok is true iff `import sympy` and `import gdsge_sympy` succeed; why is ''
%   on success, else a short human-readable reason.
%   The SymPy backend is the ONLY consumer; the default autodiff path never
%   calls this, so a missing/unsynced venv only affects sympy requests.
%
%   Python runs OUT-OF-PROCESS. The bridge only exchanges JSON strings (see
%   callSympy), so the IPC boundary is transparent and results are identical.
%   This isolates MATLAB from a native crash (access violation in MATLAB's
%   error-handling path) observed when the Python interpreter is loaded
%   in-process alongside the OpenMP solver MEX over a long run.
%
%   Other code (startup.m, another toolbox) may have loaded a foreign Python
%   into the session before us; pyenv Version/ExecutionMode are only settable
%   while NotLoaded/Terminated. Strategy when Loaded: probe by importing (any
%   SymPy-capable interpreter works); on failure terminate-and-repoint if
%   OutOfProcess, warn (restart MATLAB) if InProcess.
%
%   port (tests only): struct of effect functions {venvPython,getEnv,setEnv,
%   terminate,tryImports} so unit tests drive every branch without touching
%   the process-global pyenv state.
here = fileparts(mfilename('fullpath'));                 % src/+gdsge/+codegen/+sympy
repoRoot = fileparts(fileparts(fileparts(fileparts(here))));
pyextDir = fullfile(repoRoot, 'pyext');
if nargin < 1
    port = struct( ...
        'venvPython', @() gdsge.codegen.sympy.venvPython(pyextDir), ...
        'getEnv',     @() pyenv, ...
        'setEnv',     @(p) pyenv('Version', p, 'ExecutionMode', 'OutOfProcess'), ...
        'terminate',  @() terminate(pyenv), ...
        'tryImports', @() tryImports(pyextDir));
end

ok = false;
venvPy = port.venvPython();
if isempty(venvPy)
    why = 'the SymPy environment is not set up (run gdsge_setup_sympy)';
    return;
end

try
    pe = port.getEnv();
    if pe.Status ~= "Loaded"                 % NotLoaded or Terminated: settable
        port.setEnv(venvPy);
        ok = port.tryImports();
        why = importsWhy(ok);
        return;
    end

    % Python already loaded (possibly a foreign interpreter loaded by other
    % code). Probe by importing rather than comparing Executable paths — on
    % macOS the venv python is a symlink and pyenv may report a resolved path.
    if port.tryImports()
        ok = true; why = '';
        return;
    end

    if pe.ExecutionMode == "OutOfProcess"
        % Wrong interpreter, but replaceable: terminate and re-point.
        port.terminate();
        port.setEnv(venvPy);
        ok = port.tryImports();
        why = importsWhy(ok);
    else
        % InProcess: cannot be unloaded without restarting MATLAB.
        why = sprintf( ...
            ['a different Python (%s) was already loaded in-process by other ' ...
             'code in this MATLAB session and does not provide SymPy. Restart ' ...
             'MATLAB and re-run; GDSGE will then load its own Python ' ...
             'out-of-process.'], pe.Executable);
        warning('gdsge:sympy:inProcessMismatch', '%s', why);
    end
catch ME
    ok = false;
    why = ME.message;
end
end

function w = importsWhy(ok)
if ok
    w = '';
else
    w = 'the venv interpreter failed to import sympy/gdsge_sympy (re-run gdsge_setup_sympy)';
end
end

function ok = tryImports(pyextDir)
% Make pyext/ importable (the package lives at pyext/gdsge_sympy) and probe.
ok = false;
try
    if count(py.sys.path, pyextDir) == 0
        insert(py.sys.path, int32(0), pyextDir);
    end
    py.importlib.import_module('sympy');
    py.importlib.import_module('gdsge_sympy');
    ok = true;
catch
end
end
```

- [ ] **Step 4: Run the tests to verify they pass**

Run:

```
matlab -batch "results = runtests('tests/codegen/tEnsurePyenv.m'); disp(table(results)); assert(all([results.Passed]))"
```

Expected: PASS (7 tests).

- [ ] **Step 5: Confirm the real single-output call sites still work**

The venv exists on this machine, so the real path should succeed end to end:

```
matlab -batch "addpath('src'); ok = gdsge.codegen.sympy.ensurePyenv(); fprintf('ok=%d\n', ok); assert(ok)"
```

Expected: `ok=1`, exit 0.

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+codegen/+sympy/ensurePyenv.m tests/codegen/tEnsurePyenv.m
git commit -m "fix(sympy): handle a preloaded MATLAB Python in ensurePyenv

pyenv Version/ExecutionMode are only settable while NotLoaded, so a
foreign interpreter already Loaded (by startup.m or another toolbox)
made ensurePyenv silently import sympy in the wrong Python and report
false. Now: probe by importing; terminate-and-repoint when the wrong
interpreter is OutOfProcess; warn actionably (restart MATLAB) when it
is InProcess. Side effects live behind an injectable port so every
branch is unit-tested without touching real pyenv state.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01YCcAffmffoVzww6W1CEmSW"
```

---

### Task 2: `gdsge_setup_sympy` validates the session, not just the venv

**Files:**
- Modify: `src/gdsge_setup_sympy.m:65-72` (the "confirm interpreter + success message" tail)

**Interfaces:**
- Consumes: `[ok, why] = gdsge.codegen.sympy.ensurePyenv()` from Task 1.
- Produces: nothing new — setup output changes only.

- [ ] **Step 1: Replace the setup tail**

In `src/gdsge_setup_sympy.m`, replace:

```matlab
% 3) Confirm the interpreter exists (same resolver ensurePyenv uses).
venvPy = gdsge.codegen.sympy.venvPython(pyextDir);
if isempty(venvPy)
    error('gdsge:setup:noVenv', 'uv sync ran but no interpreter was found.');
end
fprintf(['[gdsge] SymPy backend ready. It is auto-selected when you run a model;\n' ...
         '        force it per-model with UseAutoDiff=0; (UseAutoDiff=1; forces adept).\n']);
end
```

with:

```matlab
% 3) Confirm the interpreter exists (same resolver ensurePyenv uses).
venvPy = gdsge.codegen.sympy.venvPython(pyextDir);
if isempty(venvPy)
    error('gdsge:setup:noVenv', 'uv sync ran but no interpreter was found.');
end

% 4) Confirm THIS session can use it. pyenv may already be loaded with a
% foreign interpreter by other code (startup.m, another toolbox); ensurePyenv
% recovers when it can (OutOfProcess) and explains when it cannot (InProcess).
[ok, why] = gdsge.codegen.sympy.ensurePyenv();
if ok
    fprintf(['[gdsge] SymPy backend ready. It is auto-selected when you run a model;\n' ...
             '        force it per-model with UseAutoDiff=0; (UseAutoDiff=1; forces adept).\n']);
else
    fprintf(['[gdsge] The SymPy environment was created, but this MATLAB session cannot use it:\n' ...
             '        %s\n'], why);
end
end
```

- [ ] **Step 2: Verify the healthy path**

Run:

```
matlab -batch "addpath('src'); gdsge_setup_sympy"
```

Expected: uv sync output, then `[gdsge] SymPy backend ready. ...` (ensurePyenv loads the venv Python OutOfProcess in that session, validating the whole chain), exit 0.

- [ ] **Step 3: Verify the poisoned-session path (mirrors the field report)**

Preload a SymPy-less Python InProcess first, then run setup in the same session:

```
matlab -batch "addpath('src'); pyenv('Version', 'C:\Users\WenlanLUO\AppData\Roaming\uv\python\cpython-3.12.11-windows-x86_64-none\python.exe', 'ExecutionMode', 'InProcess'); py.list; gdsge_setup_sympy"
```

Expected: uv sync output, the `gdsge:sympy:inProcessMismatch` warning, then
`[gdsge] The SymPy environment was created, but this MATLAB session cannot use it:` followed by the restart-MATLAB reason. Exit 0 (it is a notice, not an error).

- [ ] **Step 4: Commit**

```bash
git add src/gdsge_setup_sympy.m
git commit -m "fix(sympy): gdsge_setup_sympy validates the session, not just the venv

Previously setup printed success even when a foreign Python was
already loaded in-process, and the failure only surfaced later as
'Python not detected' at codegen time.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01YCcAffmffoVzww6W1CEmSW"
```

---

### Task 3: Accurate messages + end-to-end validation

**Files:**
- Modify: `src/+gdsge/+codegen/resolveBackend.m:47`
- Modify: `src/+gdsge/+codegen/generateCxx.m:42-45`
- Modify: `tests/codegen/tResolveBackend.m:66-72` (`autoPythonAbsent`)
- Modify: `tests/codegen/tBackendMessage.m:16-20` (`autoAutodiff`)

**Interfaces:**
- Consumes: `[ok, why] = gdsge.codegen.sympy.ensurePyenv()` from Task 1.
- Produces: auto-detect miss reason is exactly `'SymPy Python unavailable (run gdsge_setup_sympy)'`.

- [ ] **Step 1: Write the failing test**

In `tests/codegen/tResolveBackend.m`, extend `autoPythonAbsent` with a reason assertion:

```matlab
        function autoPythonAbsent(tc)
            tc.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'GDSGE_BACKEND', ''));
            d = gdsge.codegen.resolveBackend(tResolveBackend.irNoField(), @() false);
            tc.verifyEqual(d.backend, 'autodiff');
            tc.verifyEqual(d.mode, 'auto');
            tc.verifyEqual(d.reason, 'SymPy Python unavailable (run gdsge_setup_sympy)');
        end
```

And in `tests/codegen/tBackendMessage.m`, update `autoAutodiff` to use the new canonical reason (backendMessage embeds the reason verbatim, so both literals change together):

```matlab
        function autoAutodiff(tc)
            m = gdsge.codegen.backendMessage(struct('backend','autodiff','mode','auto', ...
                'reason','SymPy Python unavailable (run gdsge_setup_sympy)'));
            tc.verifyEqual(m, ['Backend: SymPy Python unavailable (run gdsge_setup_sympy) ' ...
                '— using adept autodiff (auto).']);
        end
```

- [ ] **Step 2: Run the two test files to verify the new assertion fails**

Run:

```
matlab -batch "results = runtests({'tests/codegen/tResolveBackend.m','tests/codegen/tBackendMessage.m'}); disp(table(results)); assert(all([results.Passed]))"
```

Expected: `autoPythonAbsent` FAILS (reason is still `'Python not detected'`). `autoAutodiff` passes already (it only exercises formatting) — that is fine; the failing gate is `autoPythonAbsent`.

- [ ] **Step 3: Update the two source messages**

In `src/+gdsge/+codegen/resolveBackend.m`, replace:

```matlab
    dec = struct('backend', 'autodiff', 'mode', 'auto', 'reason', 'Python not detected');
```

with:

```matlab
    dec = struct('backend', 'autodiff', 'mode', 'auto', ...
        'reason', 'SymPy Python unavailable (run gdsge_setup_sympy)');
```

In `src/+gdsge/+codegen/generateCxx.m`, replace:

```matlab
    if ~gdsge.codegen.sympy.ensurePyenv()
        error('gdsge:codegen:sympyPythonUnavailable', ...
            ['sympy backend requires the uv Python env. Run: uv sync --project pyext']);
    end
```

with:

```matlab
    [pyOk, pyWhy] = gdsge.codegen.sympy.ensurePyenv();
    if ~pyOk
        error('gdsge:codegen:sympyPythonUnavailable', ...
            'sympy backend requires the GDSGE Python environment: %s', pyWhy);
    end
```

- [ ] **Step 4: Run the two test files to verify they pass**

Run:

```
matlab -batch "results = runtests({'tests/codegen/tResolveBackend.m','tests/codegen/tBackendMessage.m'}); disp(table(results)); assert(all([results.Passed]))"
```

Expected: PASS.

- [ ] **Step 5: End-to-end validation of both field-report scenarios**

InProcess poisoning (the macOS user's case — expect warning + adept fallback with the new reason):

```
matlab -batch "addpath('src'); pyenv('Version', 'C:\Users\WenlanLUO\AppData\Roaming\uv\python\cpython-3.12.11-windows-x86_64-none\python.exe', 'ExecutionMode', 'InProcess'); py.list; dec = gdsge.codegen.resolveBackend(struct('model', {{}})); fprintf('backend=%s reason=%s\n', dec.backend, dec.reason)"
```

Expected: the `gdsge:sympy:inProcessMismatch` warning naming the fake-loaded interpreter, then `backend=autodiff reason=SymPy Python unavailable (run gdsge_setup_sympy)`.

OutOfProcess poisoning (expect silent auto-recovery to sympy):

```
matlab -batch "addpath('src'); pyenv('Version', 'C:\Users\WenlanLUO\AppData\Roaming\uv\python\cpython-3.12.11-windows-x86_64-none\python.exe', 'ExecutionMode', 'OutOfProcess'); py.list; dec = gdsge.codegen.resolveBackend(struct('model', {{}})); pe = pyenv; fprintf('backend=%s exe=%s\n', dec.backend, pe.Executable)"
```

Expected: `backend=sympy` and `exe=D:\refactor_gdsge\pyext\.venv\Scripts\python.exe` (the venv interpreter replaced the foreign one).

- [ ] **Step 6: Run the full test suite**

Run:

```
matlab -batch "cd('tests'); run_tests"
```

Expected: exit code 0, all pass (`tests/results/junit.xml` is authoritative). Takes a few minutes; no other MATLAB process running.

- [ ] **Step 7: Commit**

```bash
git add src/+gdsge/+codegen/resolveBackend.m src/+gdsge/+codegen/generateCxx.m tests/codegen/tResolveBackend.m tests/codegen/tBackendMessage.m
git commit -m "fix(codegen): accurate messages when SymPy Python is unavailable

'Python not detected' was misleading when a foreign interpreter was
loaded; the auto-detect reason now says what is missing and how to fix
it, and the explicit-sympy guard error carries ensurePyenv's reason.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01YCcAffmffoVzww6W1CEmSW"
```
