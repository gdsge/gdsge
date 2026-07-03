# SymPy backend vs. a preloaded MATLAB Python — design

**Date:** 2026-07-03
**Status:** Approved.
**Area:** `src/+gdsge/+codegen/+sympy/ensurePyenv.m`, `src/gdsge_setup_sympy.m`,
`src/+gdsge/+codegen/resolveBackend.m`.

## Problem (field report, reproduced locally)

A macOS user ran `gdsge_setup_sympy` successfully, yet codegen reported
`Backend: Python not detected — using adept autodiff (auto)`. Their session showed:

```
pyenv → Status: Loaded, ExecutionMode: InProcess,
        Executable: /opt/homebrew/opt/python@3.12/bin/python3.12
```

Some other code (startup.m, another toolbox — anything touching `py.*`) had loaded a
Python into MATLAB **before** GDSGE ran. The failure chain, confirmed by a local
reproduction (preload a SymPy-less Python InProcess, then call the backend):

1. `ensurePyenv.m` guards `if pe.Status ~= "Loaded"` before pointing `pyenv` at the
   venv — so with Python already loaded it silently skips the re-point. (MATLAB only
   allows changing `Version`/`ExecutionMode` while `NotLoaded`/`Terminated`; an
   InProcess interpreter cannot be unloaded without restarting MATLAB.)
2. It then imports `sympy` in the *foreign* interpreter, which lacks SymPy; the bare
   `catch` swallows the error and returns false.
3. `resolveBackend` reports "Python not detected" — misleading — and
   `gdsge_setup_sympy` never checks the session's pyenv state, so setup prints
   success while the session is already unusable for the venv.

Also verified locally: when the wrong interpreter is loaded **OutOfProcess**,
`terminate(pyenv)` moves Status to `Terminated`, after which re-pointing at the venv
and importing SymPy works.

## Decision (locked)

For a wrong interpreter loaded OutOfProcess: **terminate and re-point** (the user's
in-session `py.*` state is lost; Python transparently restarts on next use). For
InProcess there is no in-session escape — warn clearly and fall back to adept.

## Change 1 — `ensurePyenv`: handle Loaded explicitly, return a reason

New signature: `[ok, why] = ensurePyenv(port)`. `why` is `''` on success, else a
short human-readable reason. `port` is an optional injected struct of effect
functions (for unit tests), defaulting to the real implementations:

- `venvPython` → `@() gdsge.codegen.sympy.venvPython(pyextDir)` (injectable so unit
  tests do not depend on a real venv existing)
- `getEnv`     → `@() pyenv`
- `setEnv`     → `@(venvPy) pyenv('Version', venvPy, 'ExecutionMode', 'OutOfProcess')`
- `terminate`  → `@() terminate(pyenv)`
- `tryImports` → inserts `pyext/` into `py.sys.path` and imports `sympy` +
  `gdsge_sympy` inside try/catch; returns logical

Flow:

1. `venvPython` empty → `ok=false`, `why` = "environment not set up (run gdsge_setup_sympy)".
2. Status ≠ `Loaded` (covers `NotLoaded` and `Terminated`) → `setEnv(venvPy)`,
   then `tryImports` decides `ok` (import failure → `why` = import problem).
3. Status = `Loaded` → `tryImports` **first**, regardless of Executable. If it
   succeeds, `ok=true` — the loaded interpreter has SymPy and the bridge only
   exchanges JSON strings, so any modern interpreter works. (Trying imports first
   also sidesteps path-aliasing issues: on macOS the venv `bin/python` is a symlink
   and `pyenv` may report a resolved path, so string-comparing Executable against
   the venv path is unreliable.)
4. Loaded + imports failed + `ExecutionMode == "OutOfProcess"` → `terminate`,
   `setEnv(venvPy)`, retry `tryImports`.
5. Loaded + imports failed + `ExecutionMode == "InProcess"` → emit
   `warning('gdsge:sympy:inProcessMismatch', ...)`: a different Python
   (`<Executable>`) was already loaded in-process by other code in this session and
   does not provide SymPy; restart MATLAB and re-run — GDSGE will then load its own
   Python out-of-process. `ok=false`, `why` mirrors the warning.

The existing single-output call sites (`resolveBackend`'s `pyAvailableFn`,
`generateCxx`) keep working unchanged (`nargout=1`).

## Change 2 — `gdsge_setup_sympy`: validate the *session*, not just the venv

After the interpreter-exists check, call `[ok, why] = gdsge.codegen.sympy.ensurePyenv()`:

- `ok` → print the current success message (and the venv Python is now loaded
  OutOfProcess, validating the whole path immediately).
- `~ok` → print: environment created, **but this MATLAB session cannot use it**:
  `<why>`. This turns the field-report UX ("setup ready" → "Python not detected")
  into an actionable message at setup time.

## Change 3 — `resolveBackend`: accurate auto-detect reason

The auto-detect miss reason changes from `'Python not detected'` to
`'SymPy Python unavailable (run gdsge_setup_sympy)'` — in the mismatch case the
`inProcessMismatch` warning supplies the detail separately. The `pyAvailableFn`
injection contract (returns logical) is unchanged.

## Testing

- New `tests/codegen/tEnsurePyenv.m` driving `ensurePyenv(port)` with fake ports —
  no real pyenv state touched:
  - venv missing → false, `why` mentions setup.
  - NotLoaded → `setEnv` called with venv path, imports tried, true.
  - Loaded + imports succeed → true; `setEnv`/`terminate` **not** called.
  - Loaded OutOfProcess + first import fails → `terminate` then `setEnv` then
    retry; true when retry succeeds.
  - Loaded InProcess + import fails → `gdsge:sympy:inProcessMismatch` warning,
    false, `terminate` **not** called.
- Existing real-path tests (`tSympyBridge`, `sympyAvailable`) unchanged and passing.
- Manual validation (matlab -batch, mirrors the field report): preload a SymPy-less
  Python InProcess → warning + adept fallback with the new reason; preload it
  OutOfProcess → auto-recovers and selects sympy.

## Risk

Low. The only behavior change for healthy sessions is the imports-first probe when
Python is already loaded (previously it also imported into the loaded interpreter —
same effect). `terminate` fires only when the loaded interpreter is OutOfProcess
*and* provably lacks SymPy. Windows/macOS/Linux all follow the same code path.
