# Single `addpath src` — self-registering kernels

**Date:** 2026-06-16
**Status:** Approved (design)

## Problem

Running the new toolbox used to require only:

```matlab
addpath('src')
```

The `+gdsge` package and the top-level shims (`gdsge.m`, `gdsge_codegen.m`) resolve from
that single entry. But the flat-named runtime kernels now live in `src/kernels/`
(`interp_construct_mex`, `interp_eval_mex`, `asg_mex`, `gen_discrete_markov_rn.m`,
`essential_blas.dll`, plus `.cpp`/`.cache` build artifacts). MATLAB's `addpath` is **not
recursive**, so this subfolder needs its own path entry, and consumers currently must do:

```matlab
addpath('src'); addpath('src/kernels');
```

Goal: restore `addpath('src')` as the **only** line a user needs, **lazily** — kernels must
resolve once the user actually invokes the toolbox (`gdsge_codegen` / generated `iter_*` /
`simulate_*`), not necessarily from the bare `addpath` itself.

## Key existing facts

- `gdsge.runtime.ensurePath()` already runs first thing in every entry path:
  - `src/+gdsge/+codegen/codegen.m:37` calls it.
  - Every generated `iter_*`/`simulate_*` emits `gdsge.runtime.ensurePath();` at the top
    (see `emitIter.m`, `emitIterAsg.m`, `emitSimulate.m`, `emitSimulateInterp.m`,
    `emitSimulateAsg.m`).
- Today `ensurePath` only fixes the Windows DLL system `PATH` for `essential_blas.dll`, and
  it **assumes** `src/kernels` is already on the MATLAB path — it locates the DLL with
  `which('essential_blas.dll')` and errors `"addpath src/kernels first"` if absent.
- The `ensure*Mex` compile helpers (`ensureInterpEvalMex.m`, `ensureSplineConstructMex.m`,
  `ensureAsgMex.m`) locate kernels via `fullfile(srcRoot, 'kernels')` independently of the
  MATLAB path, and `cd` into that folder to compile. They do **not** depend on the path.

Because the hook already flows through every entry point, it is the natural single place to
make the kernels self-register.

## Chosen approach

Broaden `gdsge.runtime.ensurePath()` from "DLL `PATH` only" to one coherent responsibility:
**make the kernels loadable.**

`src/+gdsge/+runtime/ensurePath.m` does, idempotently:

1. **Derive `kernelsDir` from its own location** — `mfilename('fullpath')` →
   `.../src/+gdsge/+runtime/ensurePath` → up three levels to `.../src` → `fullfile(src,
   'kernels')`. This does not depend on anything already being on the path, which breaks the
   current chicken-and-egg (`which('essential_blas.dll')` needed the folder already added).
2. **Add `kernelsDir` to the MATLAB path** if not already present (membership check against
   `strsplit(path, pathsep)` so it stays idempotent and cheap).
3. **Existing Windows DLL step**, unchanged in effect: append the kernels dir to the system
   `PATH` env var so the MEX kernels can load `essential_blas.dll`. Now guaranteed to find
   it. The misleading `"addpath src/kernels first"` error precondition is removed.

`if ~ispc; return; end` guard: the DLL step stays Windows-only, but the **MATLAB-path
registration must run on all platforms** — restructure so the `addpath` step happens before
(or independent of) the `ispc` early return.

### Why lazy is sufficient

Every public entry point calls `ensurePath()` before any kernel is invoked:

- `gdsge_codegen` → `codegen.m:37` → `ensurePath()` (also compiles kernels via `ensure*Mex`,
  which are path-independent).
- Generated `iter_*` / `simulate_*` emit `ensurePath();` as their first runtime statement,
  ahead of any `interp_construct_mex` / `interp_eval_mex` / `asg_mex` call.

So in a fresh session, `addpath('src')` makes `gdsge.runtime.ensurePath` resolvable, and the
first toolbox call registers the kernels before they are needed — including runtime-only
sessions that load a saved result and call `simulate_*` without re-running codegen.

### Naming

The function now does two related things under one intent ("ensure kernels are usable").
Keep the name **`ensurePath`** to avoid touching the ~6 emit sites + `codegen.m`; update its
doc header to describe the broadened responsibility. (A rename to `ensureKernels` was
considered and declined as unnecessary churn for "refactor a bit".)

## Consumers updated to a single `addpath`

- `README.md:23` — `addpath('src'); addpath('src/kernels');` → `addpath('src');`
- `docs/user-guide.md:12` — same change.
- `tests/run_tests.m:13` — drop the explicit `addpath(.../src/kernels)`; the harness now
  relies on the hook. (`run_tests.m:12` keeps `addpath(.../src)`.)

**Left as-is (dev-only, harmless):** `tests/perf/perf_worker.m`, and the header comments in
`tests/perf/inmex_resolve_bench.m` / `inmex_randomize_bench.m`. Their explicit kernel
`addpath` becomes redundant but does no harm.

## Out of scope / unchanged

- No file moves; `src/kernels/` keeps its flat names and build artifacts.
- `ensure*Mex` compile helpers unchanged.
- Codegen call ordering unchanged.
- `IterRslt`/`SimuRslt` shapes and public API unchanged.

## Testing

1. **New test — single addpath suffices (lazy):** from a path state with **only `src`**
   added (kernels deliberately *not* on the path), invoke an entry point that triggers the
   hook, then assert a kernel resolves — e.g. `~isempty(which('interp_construct_mex'))` and
   `~isempty(which('gen_discrete_markov_rn'))` are true after the call, where they were empty
   before. Restore the path on cleanup (`onCleanup`).
2. **Idempotency:** calling `ensurePath()` twice does not duplicate the path entry.
3. **Regression:** the full existing suite (`tests/run_tests.m`) stays green with the
   `src/kernels` addpath line removed from the harness.

## Risk notes

- Honor the **Golden rule** (`CLAUDE.md`): still `addpath` exactly one source tree per
  process. This change keeps that invariant — `src/kernels` is reached only from within the
  already-added `src` tree, never alongside `base_package/gdsge/source`.
- The path-membership check must use `pathsep` (`;` on Windows) and compare against the live
  `path`, matching the existing style in `ensurePath` (`strsplit(getenv('PATH'), ';')`).
