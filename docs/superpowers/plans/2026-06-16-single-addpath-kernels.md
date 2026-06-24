# Single `addpath src` — self-registering kernels Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `addpath('src')` the only path line a user needs — the flat runtime kernels in `src/kernels/` register themselves lazily through the existing `gdsge.runtime.ensurePath()` hook on first toolbox call.

**Architecture:** `gdsge.runtime.ensurePath()` already runs first in `gdsge_codegen` (`codegen.m:37`) and is emitted at the top of every generated `iter_*`/`simulate_*`. We broaden it from "fix the Windows DLL `PATH`" to "make the kernels loadable": it locates `src/kernels` relative to its own file, `addpath`s it (idempotently), then does the existing DLL step. Consumers drop the redundant second `addpath`. One kernel-only unit test (`tKernels`) that poked kernels directly is made self-sufficient so the harness no longer needs to pre-add the folder.

**Tech Stack:** MATLAB R2025b, `matlab.unittest` suite (`tests/run_tests.m`), POSIX-agnostic path handling (`pathsep`).

---

## File Structure

- **Modify:** `src/+gdsge/+runtime/ensurePath.m` — broadened responsibility (MATLAB-path registration + DLL `PATH`).
- **Modify:** `tests/runtime/tEnsurePath.m` — add the lazy-registration proof test + keep existing DLL tests.
- **Modify:** `tests/kernels/tKernels.m` — add a `TestClassSetup` that calls `ensurePath()` so this direct-call test no longer depends on the harness pre-adding `src/kernels`.
- **Modify:** `tests/run_tests.m` — remove line 13 (`addpath .../src/kernels`); the hook now covers it.
- **Modify:** `README.md:23`, `docs/user-guide.md:12` — single `addpath('src')` in the Quickstart.

Unchanged: the `ensure*Mex` compile helpers (path-independent), codegen ordering, `IterRslt`/`SimuRslt` shapes, public API, and `src/kernels/` file layout. Dev-only `tests/perf/*` addpath lines are intentionally left as-is (harmless redundancy).

---

## Task 1: Broaden `ensurePath` to self-register the kernels

**Files:**
- Modify: `src/+gdsge/+runtime/ensurePath.m`
- Test: `tests/runtime/tEnsurePath.m`

- [ ] **Step 1: Write the failing test**

Add this method inside the existing `methods (Test)` block in `tests/runtime/tEnsurePath.m` (after `idempotent`, before the closing `end` of the block):

```matlab
        function registersKernelsOnMatlabPath(tc)
            % addpath('src') alone must be enough: the hook puts src/kernels
            % on the MATLAB path so the flat kernel names resolve.
            here       = fileparts(mfilename('fullpath'));   % tests/runtime
            root       = fileparts(fileparts(here));         % repo root
            kernelsDir = fullfile(root, 'src', 'kernels');

            origPath = path;
            restore  = onCleanup(@() path(origPath));         %#ok<NASGU>

            if ~isempty(which('interp_construct_mex'))
                rmpath(kernelsDir);
            end
            tc.assertEmpty(which('interp_construct_mex'), ...
                'precondition: src/kernels must be off the MATLAB path');

            gdsge.runtime.ensurePath();

            tc.verifyNotEmpty(which('interp_construct_mex'));
            tc.verifyNotEmpty(which('gen_discrete_markov_rn'));
        end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); addpath(fullfile(fileparts(pwd),'src','kernels')); runtests('runtime/tEnsurePath.m')"`

Expected: `registersKernelsOnMatlabPath` FAILS — after `rmpath`, the current `ensurePath` only touches the system `PATH`, so `which('interp_construct_mex')` is still empty. (`blasDirLandsOnSystemPath` and `idempotent` still pass.)

- [ ] **Step 3: Rewrite `ensurePath.m` with the broadened behavior**

Replace the entire contents of `src/+gdsge/+runtime/ensurePath.m` with:

```matlab
function ensurePath()
% ENSUREPATH  Make the flat runtime kernels in src/kernels loadable.
%   (1) Register src/kernels on the MATLAB path so the flat kernel names
%       (interp_construct_mex, interp_eval_mex, asg_mex, gen_discrete_markov_rn)
%       resolve. The folder is located relative to THIS file, so a single
%       `addpath src` is enough for callers — this hook adds the rest.
%   (2) On Windows, append that directory to the system PATH so the MEX
%       kernels can load essential_blas.dll.
%   Idempotent. Called first by gdsge_codegen and by every generated
%   iter_<model>/simulate_<model>, ahead of any kernel call.
kernelsDir = fullfile( ...
    fileparts(fileparts(fileparts(mfilename('fullpath')))), 'kernels');

if ~any(strcmp(strsplit(path, pathsep), kernelsDir))
    addpath(kernelsDir);
end

if ~ispc; return; end
if ~any(strcmp(strsplit(getenv('PATH'), ';'), kernelsDir))
    setenv('PATH', [getenv('PATH') ';' kernelsDir]);
end
end
```

Note: `mfilename('fullpath')` is `.../src/+gdsge/+runtime/ensurePath` (no extension); three `fileparts` calls strip `ensurePath` → `+runtime` → `+gdsge`, landing on `.../src`, then `/kernels`. `kernelsDir` is the same directory that holds `essential_blas.dll`, so the DLL step no longer needs `which('essential_blas.dll')`.

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); addpath(fullfile(fileparts(pwd),'src','kernels')); runtests('runtime/tEnsurePath.m')"`

Expected: all three tests PASS (`blasDirLandsOnSystemPath`, `idempotent`, `registersKernelsOnMatlabPath`).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+runtime/ensurePath.m tests/runtime/tEnsurePath.m
git commit -m "feat(runtime): ensurePath self-registers src/kernels on the MATLAB path"
```

---

## Task 2: Make `tKernels` self-sufficient and drop the harness kernel `addpath`

`tests/kernels/tKernels.m` calls `gen_discrete_markov_rn` directly with no setup, so it currently relies on `run_tests.m:13`. Make it register kernels itself (mirroring `tConstructSplines`), then remove the harness line.

**Files:**
- Modify: `tests/kernels/tKernels.m`
- Modify: `tests/run_tests.m:13`

- [ ] **Step 1: Add a `TestClassSetup` to `tKernels.m`**

In `tests/kernels/tKernels.m`, insert a setup block immediately after the class-comment lines and before `methods (Test)` (i.e. between line 4 and line 5):

```matlab
    methods (TestClassSetup)
        function registerKernels(~)
            gdsge.runtime.ensurePath();   % put src/kernels on the path
        end
    end
```

The file's `methods (Test)` block and its `markovChainShapeAndRange` test stay unchanged.

- [ ] **Step 2: Verify `tKernels` passes WITHOUT the harness kernel addpath**

Run (note: `src` added, but `src/kernels` deliberately NOT pre-added):
`matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); runtests('kernels/tKernels.m')"`

Expected: `markovChainShapeAndRange` PASSES — the new setup calls `ensurePath`, which adds `src/kernels`, so `gen_discrete_markov_rn` resolves.

- [ ] **Step 3: Remove line 13 from `tests/run_tests.m`**

Delete this line (currently `tests/run_tests.m:13`):

```matlab
addpath(fullfile(fileparts(thisDir), 'src', 'kernels'));   % vendored flat runtime kernels
```

Leave lines 11 and 12 (`addpath(thisDir)` and `addpath(.../src)`) intact.

- [ ] **Step 4: Commit**

```bash
git add tests/kernels/tKernels.m tests/run_tests.m
git commit -m "test: kernel-only test registers its own kernels; drop harness src/kernels addpath"
```

---

## Task 3: Update user-facing docs to a single `addpath`

**Files:**
- Modify: `README.md:23`
- Modify: `docs/user-guide.md:12`

- [ ] **Step 1: Edit `README.md`**

Replace line 23:

```matlab
addpath('src'); addpath('src/kernels');
```

with:

```matlab
addpath('src');
```

- [ ] **Step 2: Edit `docs/user-guide.md`**

Replace line 12 (identical text) the same way:

```matlab
addpath('src');
```

- [ ] **Step 3: Commit**

```bash
git add README.md docs/user-guide.md
git commit -m "docs: addpath('src') alone runs the toolbox"
```

---

## Task 4: Full-suite regression

**Files:** none (verification only)

- [ ] **Step 1: Run the whole suite**

Run: `matlab -batch "cd('tests'); run_tests"`

Expected: exit code 0. Pay attention to the direct-kernel-call tests that no longer get a harness addpath — `tests/kernels/tKernels.m`, `tests/kernels/tInterpEvalMex.m` (has its own `PathFixture`), `tests/runtime/tConstructSplines.m` (calls `ensurePath` in setup) — and the end-to-end tests (`tEndToEnd*`) which fire `ensurePath` through generated `iter_*`/`simulate_*`.

- [ ] **Step 2: Confirm results**

Open `tests/results/junit.xml`; verify zero failures/errors. If any test fails for a missing kernel, it is a direct-call test lacking an `ensurePath`/`PathFixture` setup — add a `TestClassSetup` calling `gdsge.runtime.ensurePath()` to that test class (same one-liner as Task 2 Step 1), then re-run.

- [ ] **Step 3: Commit (only if Step 2 required a fix)**

```bash
git add tests/
git commit -m "test: register kernels in remaining direct-call test setups"
```

---

## Self-Review notes

- **Spec coverage:** ensurePath broadening (Task 1), platform-independent MATLAB-path step before the `ispc` return (Task 1 Step 3), consumer single-`addpath` updates — README/user-guide (Task 3) and `run_tests.m` (Task 2), lazy-proof + idempotency tests (Task 1 / existing `idempotent`), full-suite regression (Task 4). The spec's "left as-is" perf/bench lines are deliberately untouched.
- **Name kept as `ensurePath`** per the approved spec — no emit-site or `codegen.m` churn.
- **Type/name consistency:** `kernelsDir`, `gdsge.runtime.ensurePath()`, `interp_construct_mex`, `gen_discrete_markov_rn` used identically across tasks.
- **Risk handled:** the only test depending on the harness kernel addpath (`tKernels`) is made self-sufficient before line 13 is removed (Task 2 ordering).
