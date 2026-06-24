# Backend auto-selection + informative compile output — design

**Date:** 2026-06-16
**Status:** approved (brainstorming) — pending spec review before planning

## Problem

Two usability gaps in the codegen driver:

1. **Backend choice is static and adept-biased.** The C++ Jacobian backend is fixed at
   parse time: `UseAutoDiff=0` in the gmod bakes `jacobianBackend='sympy'` into the IR,
   anything else falls through to adept. A user with a working Python/SymPy env still gets
   adept unless they know to set `UseAutoDiff=0`, and a SymPy request that hits an
   unsupported construct hard-errors instead of degrading gracefully. There is also no
   feedback about *which* backend ran or *why*.

2. **Compile output is uninformative.** `codegen.m` prints only `Parsing gmod file:`,
   `Compile mex file:`, and `C++ source unchanged, skip compiling.`. The three distinct
   compile steps — kernels, the model solver MEX, the simulate MEX — are unlabeled, and the
   cache-gated kernel compiles print in the middle of the parse line. The user cannot tell
   what stage is running or what is being built.

## Goals

- Make **SymPy the default when Python is detected**, falling back to **adept** when it is
  not — or when SymPy code generation fails — without the user setting anything.
- Keep **explicit user choices authoritative**: an in-gmod `UseAutoDiff` flag is honored or
  errors; it is never silently overridden.
- Print a clear, single line stating which backend is in use and why.
- Replace the compile output with **numbered phase banners** that name each stage (parsing,
  generating MATLAB, generating C++, compiling kernels, compiling solver MEX, compiling
  simulate MEX) and indicate cache skips.
- Keep the **correctness loop Python-free and deterministic** (per `CLAUDE.md`) and keep all
  existing goldens / IR snapshots valid.

## Non-goals

- No change to compile *flags*, the MEX build itself, or the `mex_<model>.cache` mechanism.
- No fallback for **compile-time** failures (a `mex` build error is a real toolchain
  problem and surfaces normally) — fallback covers code-generation failures only.
- No new symbolic capability in the SymPy backend; this is selection + messaging only.
- ASG/pchip backend coverage is unchanged from its current state.

## Decisions (from brainstorming)

| Question | Decision |
|---|---|
| Default when `UseAutoDiff` unset | **Auto** — sympy-if-Python (with fallback), else adept |
| `UseAutoDiff=1` | Force **adept** (even when Python present) |
| `UseAutoDiff=0` | Force **sympy** (error if Python missing — unchanged) |
| Explicit sympy that fails | **Error** — do not silently switch (fallback only in auto mode) |
| Fallback scope | **Codegen-time only** (unsupported syntax); not compile-time |
| Compile output | **Numbered phase banners** naming each stage + cache skips |
| Test determinism | **Global force-adept in the test harness**; sympy gates opt in |
| User override knob | **`GDSGE_BACKEND` env var** (`adept|sympy|auto`) |

## Architecture

The backend decision moves out of the parser and into the **codegen driver**
(`gdsge.codegen.codegen`), reducing `generateCxx` to "emit for the backend I'm told."
Rationale: Python availability is *environmental*, not a model property, and the
auto-fallback needs a codegen-time retry that parse-time resolution cannot express.

```
gdsge_codegen (shim)
  └─ gdsge.codegen.codegen          ← owns phase banners + backend resolution + fallback
       ├─ gdsge.parser.parseFrontEnd
       ├─ gdsge.codegen.generateMatlab
       ├─ gdsge.codegen.resolveBackend(ir)   ← NEW: {backend, mode, reason}
       ├─ gdsge.codegen.generateCxx(ir, dir, struct('backend',…))  ← opts.backend override
       │     (try sympy → catch → adept, only when mode=='auto')
       └─ feval(compile_<model>)   ← cache-gated solver + simulate MEX
```

### 1. Tri-state `resolveOptions`

`src/+gdsge/+parser/resolveOptions.m` currently does:

```matlab
if getf(ws, 'UseAutoDiff', 1) == 0
    o.jacobianBackend = 'sympy';
end
```

Change to a tri-state keyed on **presence**:

```matlab
if isfield(ws, 'UseAutoDiff')
    if double(ws.UseAutoDiff) == 0
        o.jacobianBackend = 'sympy';      % explicit force sympy (unchanged)
    else
        o.jacobianBackend = 'autodiff';   % explicit force adept (NEW)
    end
end
% absent → no field → AUTO (resolved at codegen time)
```

**Backend vocabulary:** the IR `jacobianBackend` enum is `{'autodiff','sympy'}` (see
`gdsge.ir.schema`). The internal backend token used everywhere (`resolveBackend.backend`,
`generateCxx`'s `opts.backend`) is therefore `'autodiff'` or `'sympy'`. Human-facing
messages still say "adept autodiff". The `GDSGE_BACKEND` env var accepts `adept` as a
user-friendly synonym for `autodiff` (both map to the `'autodiff'` token).

**IR meaning:** field *present* ⇒ explicit user choice (authoritative). Field *absent* ⇒
auto. No corpus `.gmod` sets `UseAutoDiff`, so the field stays absent for all of them and
**all 8 IR snapshots remain byte-identical** — no regen.

### 2. `gdsge.codegen.resolveBackend(ir)` — NEW

Pure-ish helper returning a struct `{backend, mode, reason}` where `backend ∈
{'autodiff','sympy'}`, `mode ∈ {'explicit','env','auto'}`, and `reason` is a short string
for the printed message. Precedence:

1. **In-gmod explicit** — if `isfield(ir.options,'jacobianBackend')`: `backend` = that
   value, `mode='explicit'`. Authoritative; env is ignored; no auto-fallback.
2. **`GDSGE_BACKEND` env** — else read `getenv('GDSGE_BACKEND')`:
   - `'adept'` / `'autodiff'` / `'sympy'` → that backend (`adept`→`autodiff`),
     `mode='env'` (explicit pin: honored or errors, no fallback).
   - `'auto'` or `''` (unset) → fall through to (3).
   - any other value → `error('gdsge:codegen:badBackendEnv', …)`.
3. **Auto** — `mode='auto'`. If `gdsge.codegen.sympy.ensurePyenv()` ⇒ `backend='sympy'`,
   `reason='Python detected'`; else `backend='adept'`, `reason='Python not detected'`.

For testability, Python availability is injected: the signature is
`resolveBackend(ir, pyAvailableFn)` with `pyAvailableFn` defaulting to
`@gdsge.codegen.sympy.ensurePyenv`. Tests pass a stub. `getenv` is read directly (tests set
the real env var, since the harness already manipulates it).

### 3. `generateCxx` gains `opts.backend`

`src/+gdsge/+codegen/generateCxx.m` selects the model backend from
`ir.options.jacobianBackend` today. Add an override:

```matlab
if isfield(opts, 'backend')
    useSympy = strcmp(opts.backend, 'sympy');
else
    useSympy = isfield(ir.options,'jacobianBackend') && ...
               strcmp(ir.options.jacobianBackend,'sympy');   % back-compat for direct callers
end
```

All existing guards (sympy+pchip, sympy+cxx-hooks, `ensurePyenv` check) stay. `codegen.m`
always passes an explicit `opts.backend`, so the driver fully controls selection.

### 4. Fallback in `codegen.m`

```matlab
dec = gdsge.codegen.resolveBackend(ir);
printBackend(dec);                              % the "Backend: …" line
try
    files = gdsge.codegen.generateCxx(ir, pwd, struct('backend', dec.backend));
catch err
    if strcmp(dec.mode,'auto') && strcmp(dec.backend,'sympy')
        printBackendFallback(err);              % "SymPy codegen failed (…) — falling back…"
        files = gdsge.codegen.generateCxx(ir, pwd, struct('backend','adept'));
    else
        rethrow(err);                           % explicit/env pin → honor it
    end
end
```

Fallback re-runs `generateCxx` with adept; the kernel-ensure calls inside are cache-gated
and idempotent, so the retry is cheap. The cache (`mex_<model>.cache`) keys on cpp text, so
whichever backend wins, a changed `.cpp` triggers the right recompile.

### 5. Phase banners

A small printer `gdsge.codegen.progress` (or inline helpers in `codegen.m`) emits numbered
banners. Total phase count is computed up front and adapts to the model (ASG vs spline
kernels; whether a simulate MEX exists):

```
[1/5] Parsing gmod (HL1996) ... ok
[2/5] Generating MATLAB (iter_HL1996.m, simulate_HL1996.m) ... ok
[3/5] Generating C++ (mex_HL1996.cpp)
      Backend: Python detected — using SymPy analytic Jacobian (auto; set UseAutoDiff=1 to force adept).
      ... ok
[4/5] Compiling kernels: interp_construct_mex, interp_eval_mex ... up to date
[5/5] Compiling solver MEX (mex_HL1996) ... unchanged — skipped
      Compiling simulate MEX (simulate_HL1996_mex) ...      (only when present)
```

- **Labels** stay honest: *kernels* (the `ensure*Mex` set for the interp method), *solver
  MEX* (`mex_<model>`, the task MEX that drives iter), *simulate MEX*
  (`simulate_<model>_mex`, only when `isSimuMexExpressible`).
- The `ensure*Mex` functions keep their existing per-kernel "Compiling …" lines but the
  driver prints an "up to date" line for the kernel phase when nothing recompiled.
- The backend line prints **every run**, including cache-skip runs, so the active backend is
  always visible.
- Kernel-ensuring stays inside `generateCxx` (it needs `asg.get_mex_constants()` before
  `emitCompile`); the driver wraps it with the phase banner rather than relocating it.
  `generateCxx` returns enough info (resolved backend; which kernels exist for the method)
  for the driver to label phase 4/5 accurately.

## Messaging catalogue

| Situation | Line |
|---|---|
| auto + Python | `Backend: Python detected — using SymPy analytic Jacobian (auto; set UseAutoDiff=1 to force adept).` |
| auto + no Python | `Backend: Python not detected — using adept autodiff (auto).` |
| auto fallback | `Backend: SymPy codegen failed (<short reason>) — falling back to adept autodiff.` |
| explicit adept (gmod) | `Backend: adept autodiff (UseAutoDiff=1).` |
| explicit sympy (gmod) | `Backend: SymPy analytic Jacobian (UseAutoDiff=0).` |
| env pin | `Backend: adept autodiff (GDSGE_BACKEND=adept).` / `Backend: SymPy analytic Jacobian (GDSGE_BACKEND=sympy).` |

`<short reason>` is `err.message` trimmed to its first line.

## Test plan

New/updated tests (TDD — failing test first):

1. **`tResolveBackend`** (unit, no compile): tri-state from IR field; env precedence
   (`GDSGE_BACKEND` honored when no gmod field; gmod field beats env); auto with injected
   Python-available true → sympy, false → adept; bad env value → `badBackendEnv`. Saves and
   restores `GDSGE_BACKEND` around each case.
2. **Auto-fallback gate**: a gmod that SymPy cannot generate (an unsupported construct from
   `docs/deferred-features.md`, or an injected failing emitter) run in auto mode with Python
   available → produces **adept** `mex_<model>.cpp` and prints the fallback notice. Verify
   via `evalc` capture + emitted-cpp inspection.
3. **Explicit-sympy-fails-errors**: same unsupported construct with `UseAutoDiff=0` →
   raises (no fallback).
4. **Messaging/banner assertions** (`evalc`): the numbered banners and the correct
   `Backend:` line appear for representative runs (auto→sympy, forced adept, cache-skip
   second run).
5. **Harness force-adept**: `tests/run_tests.m` sets `setenv('GDSGE_BACKEND','adept')` with
   a cleanup that restores the prior value; a spot-check confirms a default gate still emits
   adept cpp. Existing sympy gates (`UseAutoDiff=0`) are unaffected — gmod field beats env.

Existing functional goldens are unchanged (forced-adept harness). IR snapshots unchanged
(no corpus gmod sets `UseAutoDiff`). The two HL1996 text snapshots are unaffected by
codegen *content* but the **compile_<model>.m** template is untouched, so no codegen-snapshot
regen is expected; if the driver's printed output is captured anywhere, that capture updates.

## Risks / open

- **`evalc` over a real compile is slow.** Keep banner-assertion tests on cache-warm runs or
  small/fake models where possible; reuse existing compiled artifacts.
- **Choosing the "unsupported construct"** for the fallback test: prefer a real one from the
  deferred-features map so the test stays meaningful; fall back to an injected failing
  emitter only if no cheap real construct exists. Decide during planning.
- **Phase-count arithmetic** must match the actual steps (ASG has a different kernel set; no
  simulate MEX for non-`SIMU_INTERP` models). Compute from `ir.options.interpMethod` and the
  `generateCxx` return, not a hard-coded constant.
```
