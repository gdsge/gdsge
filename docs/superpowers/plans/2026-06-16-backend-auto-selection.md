# Backend Auto-Selection + Informative Compile Output Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make SymPy the default C++ Jacobian backend when Python is detected (falling back to adept autodiff when it is not, or when SymPy code generation fails), keep explicit user choices authoritative, and replace the codegen driver's terse prints with numbered phase banners that name each compile stage.

**Architecture:** The backend decision moves out of the parser into the codegen driver (`gdsge.codegen.codegen`). A new pure `resolveBackend` returns `{backend, mode, reason}` using precedence *in-gmod `UseAutoDiff` > `GDSGE_BACKEND` env > auto-detect*. `generateCxx` gains an `opts.backend` override so the driver fully controls selection; a thin `generateCxxWithFallback` retries with autodiff only when an *auto*-mode sympy run throws during codegen. The driver prints a one-line backend message and numbered `[i/5]` phase banners.

**Tech Stack:** MATLAB R2025b (`matlab.unittest`), the existing `+gdsge.+codegen`/`+parser` packages, C MEX kernels. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-06-16-backend-auto-selection-design.md`

**Backend vocabulary (important):** the IR `jacobianBackend` enum is `{'autodiff','sympy'}`
(`gdsge.ir.schema`). The internal token used by `resolveBackend.backend` and
`generateCxx`'s `opts.backend` is `'autodiff'` or `'sympy'`. Human-facing messages say
"adept autodiff". `GDSGE_BACKEND` accepts `adept` as a synonym for `autodiff`.

**Test run convention:** from the repo root `D:\refactor_gdsge`, run a single class with:
`matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/<path>.m'); assert(~any([r.Failed]))"`
The full suite is `matlab -batch "cd('tests'); run_tests"` (exit 0 = all pass).

---

### Task 1: Tri-state `UseAutoDiff` in `resolveOptions`

Make backend selection a tri-state keyed on **presence** of `UseAutoDiff`: absent → no IR
field (auto), `=1` → `'autodiff'`, `=0` → `'sympy'`.

**Files:**
- Modify: `src/+gdsge/+parser/resolveOptions.m:38-43`
- Test: `tests/parser/tResolveOptionsBackend.m` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/parser/tResolveOptionsBackend.m`:

```matlab
classdef tResolveOptionsBackend < matlab.unittest.TestCase
    % Tri-state UseAutoDiff -> options.jacobianBackend in the parser.
    methods (TestClassSetup)
        function srcPath(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src')));
        end
    end
    methods (Test)
        function unsetIsAuto(tc)
            o = gdsge.parser.resolveOptions(struct());
            tc.verifyFalse(isfield(o, 'jacobianBackend'), ...
                'absent UseAutoDiff must leave jacobianBackend unset (auto)');
        end
        function oneForcesAutodiff(tc)
            o = gdsge.parser.resolveOptions(struct('UseAutoDiff', 1));
            tc.verifyEqual(o.jacobianBackend, 'autodiff');
        end
        function zeroForcesSympy(tc)
            o = gdsge.parser.resolveOptions(struct('UseAutoDiff', 0));
            tc.verifyEqual(o.jacobianBackend, 'sympy');
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "addpath('src','tests'); r=runtests('tests/parser/tResolveOptionsBackend.m'); assert(~any([r.Failed]))"`
Expected: FAIL — `oneForcesAutodiff` errors (today `UseAutoDiff=1` sets no field).

- [ ] **Step 3: Write minimal implementation**

In `src/+gdsge/+parser/resolveOptions.m`, replace lines 38-43 (the comment + the
`if getf(...)==0` block) with:

```matlab
% Jacobian backend: emit the field ONLY when the gmod sets UseAutoDiff (tri-state).
% Absent -> no field -> AUTO (resolved at codegen time by gdsge.codegen.resolveBackend),
% keeping default-path IR byte-identical to existing goldens. =1 -> autodiff, =0 -> sympy.
if isfield(ws, 'UseAutoDiff')
    if double(ws.UseAutoDiff) == 0
        o.jacobianBackend = 'sympy';
    else
        o.jacobianBackend = 'autodiff';
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "addpath('src','tests'); r=runtests('tests/parser/tResolveOptionsBackend.m'); assert(~any([r.Failed]))"`
Expected: PASS (3/3).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+parser/resolveOptions.m tests/parser/tResolveOptionsBackend.m
git commit -m "feat(parser): tri-state UseAutoDiff -> jacobianBackend (autodiff/sympy/auto)"
```

---

### Task 2: `gdsge.codegen.resolveBackend`

Pure decision helper with precedence in-gmod > env > auto-detect, returning
`{backend, mode, reason}`. Python availability is injected for deterministic tests.

**Files:**
- Create: `src/+gdsge/+codegen/resolveBackend.m`
- Test: `tests/codegen/tResolveBackend.m` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/codegen/tResolveBackend.m`:

```matlab
classdef tResolveBackend < matlab.unittest.TestCase
    % Backend precedence: in-gmod field > GDSGE_BACKEND env > auto-detect.
    methods (TestClassSetup)
        function srcPath(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src')));
        end
    end
    methods (Static)
        function ir = irNoField()
            ir = struct('options', struct('interpMethod', 'spline'));
        end
    end
    methods (Test)
        function explicitSympyWins(tc)
            ir = tResolveBackend.irNoField(); ir.options.jacobianBackend = 'sympy';
            d = gdsge.codegen.resolveBackend(ir, @() false);   % python "absent"
            tc.verifyEqual(d.backend, 'sympy');
            tc.verifyEqual(d.mode, 'explicit');
        end
        function explicitAutodiffWins(tc)
            ir = tResolveBackend.irNoField(); ir.options.jacobianBackend = 'autodiff';
            d = gdsge.codegen.resolveBackend(ir, @() true);
            tc.verifyEqual(d.backend, 'autodiff');
            tc.verifyEqual(d.mode, 'explicit');
        end
        function envAdeptBeatsAutoPython(tc)
            tc.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'GDSGE_BACKEND', 'adept'));
            d = gdsge.codegen.resolveBackend(tResolveBackend.irNoField(), @() true);
            tc.verifyEqual(d.backend, 'autodiff');
            tc.verifyEqual(d.mode, 'env');
        end
        function envSympy(tc)
            tc.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'GDSGE_BACKEND', 'sympy'));
            d = gdsge.codegen.resolveBackend(tResolveBackend.irNoField(), @() false);
            tc.verifyEqual(d.backend, 'sympy');
            tc.verifyEqual(d.mode, 'env');
        end
        function explicitBeatsEnv(tc)
            tc.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'GDSGE_BACKEND', 'sympy'));
            ir = tResolveBackend.irNoField(); ir.options.jacobianBackend = 'autodiff';
            d = gdsge.codegen.resolveBackend(ir, @() false);
            tc.verifyEqual(d.backend, 'autodiff');
            tc.verifyEqual(d.mode, 'explicit');
        end
        function autoPythonPresent(tc)
            tc.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'GDSGE_BACKEND', ''));   % cleared -> auto
            d = gdsge.codegen.resolveBackend(tResolveBackend.irNoField(), @() true);
            tc.verifyEqual(d.backend, 'sympy');
            tc.verifyEqual(d.mode, 'auto');
        end
        function autoPythonAbsent(tc)
            tc.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'GDSGE_BACKEND', ''));
            d = gdsge.codegen.resolveBackend(tResolveBackend.irNoField(), @() false);
            tc.verifyEqual(d.backend, 'autodiff');
            tc.verifyEqual(d.mode, 'auto');
        end
        function badEnvErrors(tc)
            tc.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'GDSGE_BACKEND', 'bogus'));
            tc.verifyError(@() gdsge.codegen.resolveBackend( ...
                tResolveBackend.irNoField(), @() true), 'gdsge:codegen:badBackendEnv');
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "addpath('src','tests'); r=runtests('tests/codegen/tResolveBackend.m'); assert(~any([r.Failed]))"`
Expected: FAIL — `resolveBackend` does not exist.

- [ ] **Step 3: Write minimal implementation**

Create `src/+gdsge/+codegen/resolveBackend.m`:

```matlab
function dec = resolveBackend(ir, pyAvailableFn)
% RESOLVEBACKEND  Decide the C++ Jacobian backend for this codegen run.
%   dec = struct('backend',b,'mode',m,'reason',r) with b in {'autodiff','sympy'},
%   m in {'explicit','env','auto'}. Precedence:
%     1) in-gmod UseAutoDiff (ir.options.jacobianBackend) — authoritative
%     2) GDSGE_BACKEND env (adept|autodiff|sympy|auto)
%     3) auto-detect: Python available => sympy (with codegen fallback), else autodiff
%   pyAvailableFn (default @gdsge.codegen.sympy.ensurePyenv) is injected for tests.
if nargin < 2; pyAvailableFn = @gdsge.codegen.sympy.ensurePyenv; end

% (1) explicit in-gmod choice
if isfield(ir, 'options') && isfield(ir.options, 'jacobianBackend')
    b = ir.options.jacobianBackend;
    dec = struct('backend', b, 'mode', 'explicit', ...
        'reason', sprintf('UseAutoDiff=%d', double(strcmp(b, 'autodiff'))));
    return;
end

% (2) GDSGE_BACKEND env override
env = lower(strtrim(getenv('GDSGE_BACKEND')));
switch env
    case {'adept', 'autodiff'}
        dec = struct('backend', 'autodiff', 'mode', 'env', 'reason', 'GDSGE_BACKEND=adept');
        return;
    case 'sympy'
        dec = struct('backend', 'sympy', 'mode', 'env', 'reason', 'GDSGE_BACKEND=sympy');
        return;
    case {'', 'auto'}
        % fall through to auto-detection
    otherwise
        error('gdsge:codegen:badBackendEnv', ...
            'GDSGE_BACKEND=%s is not recognized (use adept | sympy | auto).', env);
end

% (3) auto-detection
if pyAvailableFn()
    dec = struct('backend', 'sympy', 'mode', 'auto', 'reason', 'Python detected');
else
    dec = struct('backend', 'autodiff', 'mode', 'auto', 'reason', 'Python not detected');
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "addpath('src','tests'); r=runtests('tests/codegen/tResolveBackend.m'); assert(~any([r.Failed]))"`
Expected: PASS (8/8).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/resolveBackend.m tests/codegen/tResolveBackend.m
git commit -m "feat(codegen): resolveBackend decision helper (gmod > env > auto-detect)"
```

---

### Task 3: `gdsge.codegen.backendMessage`

Pure formatter turning a `resolveBackend` decision into the user-facing one-liner.

**Files:**
- Create: `src/+gdsge/+codegen/backendMessage.m`
- Test: `tests/codegen/tBackendMessage.m` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/codegen/tBackendMessage.m`:

```matlab
classdef tBackendMessage < matlab.unittest.TestCase
    methods (TestClassSetup)
        function srcPath(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src')));
        end
    end
    methods (Test)
        function autoSympy(tc)
            m = gdsge.codegen.backendMessage(struct('backend','sympy','mode','auto', ...
                'reason','Python detected'));
            tc.verifyEqual(m, ['Backend: Python detected — using SymPy analytic ' ...
                'Jacobian (auto; set UseAutoDiff=1 to force adept).']);
        end
        function autoAutodiff(tc)
            m = gdsge.codegen.backendMessage(struct('backend','autodiff','mode','auto', ...
                'reason','Python not detected'));
            tc.verifyEqual(m, 'Backend: Python not detected — using adept autodiff (auto).');
        end
        function explicitAutodiff(tc)
            m = gdsge.codegen.backendMessage(struct('backend','autodiff','mode','explicit', ...
                'reason','UseAutoDiff=1'));
            tc.verifyEqual(m, 'Backend: adept autodiff (UseAutoDiff=1).');
        end
        function explicitSympy(tc)
            m = gdsge.codegen.backendMessage(struct('backend','sympy','mode','explicit', ...
                'reason','UseAutoDiff=0'));
            tc.verifyEqual(m, 'Backend: SymPy analytic Jacobian (UseAutoDiff=0).');
        end
        function envAdept(tc)
            m = gdsge.codegen.backendMessage(struct('backend','autodiff','mode','env', ...
                'reason','GDSGE_BACKEND=adept'));
            tc.verifyEqual(m, 'Backend: adept autodiff (GDSGE_BACKEND=adept).');
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "addpath('src','tests'); r=runtests('tests/codegen/tBackendMessage.m'); assert(~any([r.Failed]))"`
Expected: FAIL — `backendMessage` does not exist.

- [ ] **Step 3: Write minimal implementation**

Create `src/+gdsge/+codegen/backendMessage.m`:

```matlab
function msg = backendMessage(dec)
% BACKENDMESSAGE  One-line user-facing description of a resolveBackend decision.
switch dec.backend
    case 'sympy';    name = 'SymPy analytic Jacobian';
    case 'autodiff'; name = 'adept autodiff';
    otherwise; error('gdsge:codegen:badBackend', 'unknown backend %s', dec.backend);
end
switch dec.mode
    case {'explicit', 'env'}
        msg = sprintf('Backend: %s (%s).', name, dec.reason);
    case 'auto'
        if strcmp(dec.backend, 'sympy')
            msg = sprintf(['Backend: %s — using %s (auto; set UseAutoDiff=1 to ' ...
                'force adept).'], dec.reason, name);
        else
            msg = sprintf('Backend: %s — using %s (auto).', dec.reason, name);
        end
    otherwise
        error('gdsge:codegen:badBackendMode', 'unknown mode %s', dec.mode);
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "addpath('src','tests'); r=runtests('tests/codegen/tBackendMessage.m'); assert(~any([r.Failed]))"`
Expected: PASS (5/5).

> Note: the em dash `—` in the strings is a literal U+2014; keep the source file UTF-8.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/backendMessage.m tests/codegen/tBackendMessage.m
git commit -m "feat(codegen): backendMessage formatter for the backend-choice line"
```

---

### Task 4: `generateCxx` honors `opts.backend`

Let the driver force the backend regardless of the IR field. Absent `opts.backend`
keeps the legacy IR-field behavior (back-compat for the many direct callers).

**Files:**
- Modify: `src/+gdsge/+codegen/generateCxx.m:30-32`
- Test: `tests/codegen/tGenerateCxxBackendOverride.m` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/codegen/tGenerateCxxBackendOverride.m`:

```matlab
classdef tGenerateCxxBackendOverride < matlab.unittest.TestCase
    % opts.backend overrides the IR jacobianBackend field. Forcing autodiff over a
    % sympy IR needs no Python, so this stays deterministic and fast.
    methods (TestClassSetup)
        function irFixture(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(here), 'HeatonLucas1996', 'ir')));
        end
    end
    methods (Test)
        function overrideAutodiffBeatsSympyField(tc)
            ir = buildHL1996IR();
            ir.options.jacobianBackend = 'sympy';     % would otherwise need Python
            work = tempname; mkdir(work);
            files = gdsge.codegen.generateCxx(ir, work, struct('backend', 'autodiff'));
            tc.verifyTrue(exist(files.cppFile, 'file') == 2);
            txt = fileread(files.cppFile);
            % autodiff path emits the adept adouble model fn; sympy path does not.
            tc.verifyTrue(contains(txt, 'GDSGE_FUNC_MODEL_1_adouble'), ...
                'override to autodiff should emit the adept model function');
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tGenerateCxxBackendOverride.m'); assert(~any([r.Failed]))"`
Expected: FAIL — today `generateCxx` ignores `opts.backend`, reads the IR field `'sympy'`, and errors `gdsge:codegen:sympyPythonUnavailable` (or emits sympy code).

- [ ] **Step 3: Write minimal implementation**

In `src/+gdsge/+codegen/generateCxx.m`, replace lines 30-32 (the comment +
`useSympy = isfield(...) && strcmp(...)`) with:

```matlab
% --- select the C++ model backend (driver may force via opts.backend; else the
%     IR jacobianBackend field decides, defaulting to autodiff) ---
if isfield(opts, 'backend')
    useSympy = strcmp(opts.backend, 'sympy');
else
    useSympy = isfield(ir.options, 'jacobianBackend') ...
        && strcmp(ir.options.jacobianBackend, 'sympy');
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tGenerateCxxBackendOverride.m'); assert(~any([r.Failed]))"`
Expected: PASS (1/1).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/generateCxx.m tests/codegen/tGenerateCxxBackendOverride.m
git commit -m "feat(codegen): generateCxx opts.backend override (driver-controlled backend)"
```

---

### Task 5: Extract `gdsge.codegen.ensureKernels`

Pull the interp-kernel ensure branch out of `generateCxx` so the driver can announce a
"Compiling kernels" phase and both call sites share one definition.

**Files:**
- Create: `src/+gdsge/+codegen/ensureKernels.m`
- Modify: `src/+gdsge/+codegen/generateCxx.m:15-20`
- Test: `tests/codegen/tEnsureKernels.m` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/codegen/tEnsureKernels.m`:

```matlab
classdef tEnsureKernels < matlab.unittest.TestCase
    % ensureKernels builds the cartesian spline kernels for a spline model.
    methods (TestClassSetup)
        function paths(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(here), 'HeatonLucas1996', 'ir')));
        end
    end
    methods (Test)
        function buildsSplineKernels(tc)
            gdsge.codegen.ensureKernels(buildHL1996IR());   % must not error
            here = fileparts(mfilename('fullpath'));
            kernels = fullfile(fileparts(fileparts(here)), 'src', 'kernels');
            tc.verifyTrue(exist(fullfile(kernels, ['interp_construct_mex.' mexext]), 'file') == 3);
            tc.verifyTrue(exist(fullfile(kernels, ['interp_eval_mex.' mexext]), 'file') == 3);
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tEnsureKernels.m'); assert(~any([r.Failed]))"`
Expected: FAIL — `ensureKernels` does not exist.

- [ ] **Step 3: Write minimal implementation**

Create `src/+gdsge/+codegen/ensureKernels.m`:

```matlab
function ensureKernels(ir)
% ENSUREKERNELS  Compile (cache-gated) the interpolation kernels this model's C++
%   backend needs: asg_mex for ASG models; the fused spline constructor + the generic
%   evaluator for cartesian spline/linear models. Extracted from generateCxx so the
%   codegen driver can announce a "Compiling kernels" phase and reuse the same logic.
if strcmp(ir.options.interpMethod, 'asg')
    gdsge.codegen.ensureAsgMex();
else
    gdsge.codegen.ensureSplineConstructMex();
    gdsge.codegen.ensureInterpEvalMex();
end
end
```

In `src/+gdsge/+codegen/generateCxx.m`, replace the kernel-ensure block (currently lines
15-20: the `if strcmp(...'asg') ... else ... ensureSplineConstructMex/ensureInterpEvalMex ...
end`) with a single call, keeping the `pchip` error immediately above it intact:

```matlab
if strcmp(ir.options.interpMethod, 'pchip')
    error('gdsge:codegen:unsupported', 'interpMethod pchip: Phase 7c');
end
gdsge.codegen.ensureKernels(ir);   % asg_mex, or spline constructor + evaluator
```

- [ ] **Step 4: Run tests to verify pass + no regression**

Run: `matlab -batch "addpath('src','src/kernels','tests'); r=runtests({'tests/codegen/tEnsureKernels.m','tests/codegen/tGenerateCxx.m'}); assert(~any([r.Failed]))"`
Expected: PASS — `tEnsureKernels` (1/1) and the existing `tGenerateCxx` still green.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/ensureKernels.m src/+gdsge/+codegen/generateCxx.m tests/codegen/tEnsureKernels.m
git commit -m "refactor(codegen): extract ensureKernels from generateCxx (shared by driver)"
```

---

### Task 6: `gdsge.codegen.generateCxxWithFallback`

Run `generateCxx` for the resolved backend; in **auto mode only**, on a sympy codegen
failure print the fallback notice and retry with autodiff. Explicit/env pins propagate
the error. `genFn` is injected so the fallback is testable without a real failing model.

**Files:**
- Create: `src/+gdsge/+codegen/generateCxxWithFallback.m`
- Test: `tests/codegen/tBackendFallback.m` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/codegen/tBackendFallback.m`:

```matlab
classdef tBackendFallback < matlab.unittest.TestCase
    % Auto-mode sympy codegen failure falls back to autodiff; explicit/env pins error.
    methods (TestClassSetup)
        function srcPath(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src')));
        end
    end
    methods (Test)
        function autoModeFallsBack(tc)
            ir  = struct('options', struct('interpMethod', 'spline'));
            dec = struct('backend','sympy','mode','auto','reason','Python detected');
            work = tempname; mkdir(work);
            genFn = @failOnSympyGen;
            out = evalc(['[files, fb] = gdsge.codegen.generateCxxWithFallback(' ...
                'ir, work, dec, genFn);']);
            tc.verifyEqual(fb, 'autodiff');
            tc.verifyEqual(files.cppFile, fullfile(work, 'dummy.cpp'));
            tc.verifyTrue(contains(out, 'falling back to adept'));
            tc.verifyTrue(contains(out, 'synthetic sympy failure'));
            tc.verifyFalse(contains(out, 'second line'), 'only the first message line');
        end
        function explicitModeErrors(tc)
            ir  = struct('options', struct('interpMethod', 'spline'));
            dec = struct('backend','sympy','mode','explicit','reason','UseAutoDiff=0');
            work = tempname; mkdir(work); genFn = @failOnSympyGen;
            tc.verifyError(@() gdsge.codegen.generateCxxWithFallback(ir, work, dec, genFn), ...
                'gdsge:test:sympyBoom');
        end
        function envModeErrors(tc)
            ir  = struct('options', struct('interpMethod', 'spline'));
            dec = struct('backend','sympy','mode','env','reason','GDSGE_BACKEND=sympy');
            work = tempname; mkdir(work); genFn = @failOnSympyGen;
            tc.verifyError(@() gdsge.codegen.generateCxxWithFallback(ir, work, dec, genFn), ...
                'gdsge:test:sympyBoom');
        end
    end
end

function files = failOnSympyGen(~, outDir, opts)
% Synthetic generator: throws for sympy, succeeds (dummy files) for autodiff.
if strcmp(opts.backend, 'sympy')
    error('gdsge:test:sympyBoom', 'synthetic sympy failure\nsecond line');
end
files = struct('cppFile', fullfile(outDir, 'dummy.cpp'), 'simuMexFile', '');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "addpath('src','tests'); r=runtests('tests/codegen/tBackendFallback.m'); assert(~any([r.Failed]))"`
Expected: FAIL — `generateCxxWithFallback` does not exist.

- [ ] **Step 3: Write minimal implementation**

Create `src/+gdsge/+codegen/generateCxxWithFallback.m`:

```matlab
function [files, finalBackend] = generateCxxWithFallback(ir, outDir, dec, genFn)
% GENERATECXXWITHFALLBACK  Generate the C++ for the resolved backend, falling back to
%   adept autodiff ONLY when an auto-mode sympy run throws during code generation.
%   Explicit/env-pinned backends are honored (the error propagates). genFn is injected
%   for tests (default @gdsge.codegen.generateCxx).
if nargin < 4; genFn = @gdsge.codegen.generateCxx; end
finalBackend = dec.backend;
try
    files = genFn(ir, outDir, struct('backend', dec.backend));
catch err
    if strcmp(dec.mode, 'auto') && strcmp(dec.backend, 'sympy')
        lines = strsplit(err.message, newline);
        fprintf(['Backend: SymPy codegen failed (%s) — falling back to adept ' ...
            'autodiff.\n'], lines{1});
        finalBackend = 'autodiff';
        files = genFn(ir, outDir, struct('backend', 'autodiff'));
    else
        rethrow(err);
    end
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "addpath('src','tests'); r=runtests('tests/codegen/tBackendFallback.m'); assert(~any([r.Failed]))"`
Expected: PASS (3/3).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/generateCxxWithFallback.m tests/codegen/tBackendFallback.m
git commit -m "feat(codegen): generateCxxWithFallback (auto-mode sympy->autodiff fallback)"
```

---

### Task 7: Name the simulate-MEX compile in `emitCompile`

When a model emits a whole-loop simulate MEX, the generated `compile_<model>.m` should
print a labeled sub-line before that second `mex` call.

**Files:**
- Modify: `src/+gdsge/+codegen/+cxx/emitCompile.m:52-58`
- Test: `tests/codegen/tEmitCompileSimuLabel.m` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/codegen/tEmitCompileSimuLabel.m`:

```matlab
classdef tEmitCompileSimuLabel < matlab.unittest.TestCase
    % The simulate-MEX compile step is labeled in the generated compile_<model>.m.
    methods (TestClassSetup)
        function irFixture(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(here), 'HeatonLucas1996', 'ir')));
        end
    end
    methods (Test)
        function labelsSimuMexWhenPresent(tc)
            txt = gdsge.codegen.cxx.emitCompile(buildHL1996IR(), 'C:/inc', true);
            tc.verifyTrue(contains(txt, 'Compiling simulate MEX (simulate_HL1996_mex)'));
            tc.verifyTrue(contains(txt, 'simu_compile'));   % the existing block remains
        end
        function noLabelWhenAbsent(tc)
            txt = gdsge.codegen.cxx.emitCompile(buildHL1996IR(), 'C:/inc', false);
            tc.verifyFalse(contains(txt, 'Compiling simulate MEX'));
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "addpath('src','tests'); r=runtests('tests/codegen/tEmitCompileSimuLabel.m'); assert(~any([r.Failed]))"`
Expected: FAIL — `labelsSimuMexWhenPresent` finds no "Compiling simulate MEX" line.

- [ ] **Step 3: Write minimal implementation**

In `src/+gdsge/+codegen/+cxx/emitCompile.m`, the `if withSimuMex` branch currently builds
`simuMexCompile` starting with `'simu_cpp = ...'`. Prepend an `fprintf` line so the
generated code announces the step. Replace the assignment's first line:

```matlab
    simuMexCompile = [ ...
        'simu_cpp = fullfile(current_folder,''simulate_' ir.modelName '_mex.cpp'');' newline ...
```

with:

```matlab
    simuMexCompile = [ ...
        'fprintf(''      Compiling simulate MEX (simulate_' ir.modelName '_mex) ...\n'');' newline ...
        'simu_cpp = fullfile(current_folder,''simulate_' ir.modelName '_mex.cpp'');' newline ...
```

(Leave the rest of the `simuMexCompile` concatenation unchanged.)

- [ ] **Step 4: Run tests to verify pass + no snapshot regression**

Run: `matlab -batch "addpath('src','src/kernels','tests'); r=runtests({'tests/codegen/tEmitCompileSimuLabel.m','tests/HeatonLucas1996/codegen/tSnapshotCxxHL1996.m'}); assert(~any([r.Failed]))"`
Expected: PASS — the label test (2/2) and `tSnapshotCxxHL1996` still green (HL1996 is
`SIMU_RESOLVE`, so `withSimuMex=false` and `compile_HL1996.m` is unchanged).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/+cxx/emitCompile.m tests/codegen/tEmitCompileSimuLabel.m
git commit -m "feat(codegen): label the simulate-MEX compile step in compile_<model>.m"
```

---

### Task 8: Wire the driver — phase banners + backend resolution + fallback

Rewrite the print flow in `gdsge.codegen.codegen` to use numbered `[i/5]` phase banners,
resolve the backend, print the backend line every run, and route C++ generation through
`generateCxxWithFallback`. Update `tCodegenDriver` for the new "skipped" wording, pin the
backend with an env fixture, and assert the banners.

**Files:**
- Modify: `src/+gdsge/+codegen/codegen.m:21-42`
- Modify: `tests/codegen/tCodegenDriver.m` (skip wording + env fixture + banner asserts)

- [ ] **Step 1: Update the failing tests**

In `tests/codegen/tCodegenDriver.m`, add a method-level env fixture so the driver-driven
tests resolve to adept deterministically (these pre-seed an adept cache). Insert directly
after the `classdef` line's opening (before `methods (Test)`):

```matlab
    methods (TestMethodSetup)
        function pinAdept(tc)
            tc.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'GDSGE_BACKEND', 'adept'));
        end
    end
```

Then in `skipsCompileAndWritesArtifactsWhenCacheMatches`, replace the skip assertion

```matlab
            tc.verifyTrue(contains(out, 'skip compiling'), ...
                sprintf('expected the skip path, got output:\n%s', out));
```

with assertions for the new banners and skip wording:

```matlab
            tc.verifyTrue(contains(out, 'skipped'), ...
                sprintf('expected the skip path, got output:\n%s', out));
            tc.verifyTrue(contains(out, '[1/5] Parsing gmod (HL1996)'), out);
            tc.verifyTrue(contains(out, '[2/5] Generating MATLAB'), out);
            tc.verifyTrue(contains(out, '[3/5] Compiling kernels'), out);
            tc.verifyTrue(contains(out, '[4/5] Generating C++'), out);
            tc.verifyTrue(contains(out, 'Backend: adept autodiff (GDSGE_BACKEND=adept).'), out);
            tc.verifyTrue(contains(out, '[5/5]'), out);
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tCodegenDriver.m'); assert(~any([r.Failed]))"`
Expected: FAIL — the driver still prints `Parsing gmod file:` / `skip compiling`, so the
new banner/`skipped`/`Backend:` assertions fail.

- [ ] **Step 3: Rewrite the driver print flow**

In `src/+gdsge/+codegen/codegen.m`, replace the block from `fprintf('Parsing gmod file: ');`
(line 21) through the closing `end` of the `if gdsge.codegen.needsCompile(...)` else-branch
(line 42) with:

```matlab
TOTAL = 5;
banner(1, TOTAL, sprintf('Parsing gmod (%s)', modelName));
ir = gdsge.parser.parseFrontEnd(fileread(gmodFile), modelName, fileparts(gmodFile));
gdsge.codegen.writeText(fullfile(pwd, [modelName '.gdsge.json']), gdsge.ir.encode(ir));

banner(2, TOTAL, 'Generating MATLAB (iter + simulate)');
matFiles = gdsge.codegen.generateMatlab(ir, pwd);

banner(3, TOTAL, 'Compiling kernels (cache-gated)');
gdsge.codegen.ensureKernels(ir);

dec = gdsge.codegen.resolveBackend(ir);
banner(4, TOTAL, sprintf('Generating C++ (mex_%s.cpp)', modelName));
fprintf('      %s\n', gdsge.codegen.backendMessage(dec));
files = gdsge.codegen.generateCxxWithFallback(ir, pwd, dec);

cppText   = fileread(files.cppFile);
if isfield(files, 'simuMexFile') && ~isempty(files.simuMexFile)
    cppText = [cppText, fileread(files.simuMexFile)];   % recompile if simulate MEX changed
end
cacheFile = fullfile(pwd, ['mex_' modelName '.cache']);
cppCache  = '';
if exist(cacheFile, 'file'); cppCache = fileread(cacheFile); end
if gdsge.codegen.needsCompile(cppText, cacheFile)
    banner(5, TOTAL, sprintf('Compiling solver MEX (mex_%s)', modelName));
    gdsge.runtime.ensurePath();
    feval(['compile_' modelName]);
    gdsge.codegen.writeText(cacheFile, cppText);
else
    banner(5, TOTAL, sprintf('Solver MEX up to date (mex_%s)', modelName));
    fprintf('      C++ source unchanged — skipped.\n');
end
```

Add a local function at the end of `codegen.m`, after the existing `validateOptions` local
function (a `.m` function file may hold multiple local functions):

```matlab
function banner(idx, total, label)
% One numbered phase header line for the codegen driver's progress output.
fprintf('[%d/%d] %s\n', idx, total, label);
end
```

Note: `cppCache` remains assigned (the `nargout > 1` block below still reads it); `matFiles`
and `files` keep their names, so the existing `if nargout > 1` artifact-exposure block is
untouched.

- [ ] **Step 4: Run tests to verify they pass**

Run: `matlab -batch "addpath('src','src/kernels','tests'); r=runtests('tests/codegen/tCodegenDriver.m'); assert(~any([r.Failed]))"`
Expected: PASS — all `tCodegenDriver` methods (incl. the new banner assertions and the
`skipped` wording).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/codegen.m tests/codegen/tCodegenDriver.m
git commit -m "feat(codegen): numbered phase banners + backend resolution/fallback in driver"
```

---

### Task 9: Force adept in the test harness

Keep the correctness loop Python-free and deterministic (per `CLAUDE.md`): the suite runs
adept by default; the sympy gates opt in via their in-gmod `UseAutoDiff=0` (which beats the
env). Set the env var in `run_tests.m`.

**Files:**
- Modify: `tests/run_tests.m:14` (after the kernels `addpath`)

- [ ] **Step 1: Add the env default**

In `tests/run_tests.m`, immediately after line 13
(`addpath(fullfile(fileparts(thisDir), 'src', 'kernels'));`), insert:

```matlab
% Pin the correctness loop to the adept backend (Python-free, byte-identical goldens).
% Sympy gates opt in with an in-gmod UseAutoDiff=0, which beats this env override.
% Only set a default — respect an explicit value if the caller already exported one.
if isempty(getenv('GDSGE_BACKEND'))
    setenv('GDSGE_BACKEND', 'adept');
end
```

- [ ] **Step 2: Verify a default gate still emits adept and the sympy gate still opts into sympy**

Run the fast no-compile driver tests plus the IR snapshot under the harness path settings:

Run: `matlab -batch "cd('tests'); run_tests"` is the full suite (next task). For a quick
targeted check now, run:
`matlab -batch "setenv('GDSGE_BACKEND','adept'); addpath('src','src/kernels','tests','tests/HeatonLucas1996/ir'); r=runtests({'tests/codegen/tCodegenDriver.m','tests/HeatonLucas1996/codegen/tSnapshotCxxHL1996.m'}); assert(~any([r.Failed]))"`
Expected: PASS — adept cpp snapshot unchanged; driver resolves adept via the env.

- [ ] **Step 3: Commit**

```bash
git add tests/run_tests.m
git commit -m "test(harness): pin suite to adept backend (sympy gates opt in via gmod)"
```

---

### Task 10: Full-suite verification + spec/PROGRESS sync

Confirm the whole suite is green (both backends), then record the change.

**Files:**
- Modify: `PROGRESS.md` (changelog entry)

- [ ] **Step 1: Run the full suite**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: exit code 0; `tests/results/junit.xml` shows 0 failures. The new tests
(`tResolveOptionsBackend`, `tResolveBackend`, `tBackendMessage`,
`tGenerateCxxBackendOverride`, `tEnsureKernels`, `tBackendFallback`,
`tEmitCompileSimuLabel`) pass; existing sympy gates still run sympy (they set
`UseAutoDiff=0`); all adept goldens/snapshots unchanged.

- [ ] **Step 2: If any pre-existing gate regressed, triage before proceeding**

Likely culprits and fixes:
- A gate that asserted the old strings `Parsing gmod file:` / `Compile mex file:` /
  `skip compiling` — update it to the new banners/`skipped` wording (only `tCodegenDriver`
  was found to do so; fix any other).
- A standalone-run gate that resolved sympy because Python is present and `GDSGE_BACKEND`
  was unset — confirm it passes under the suite (harness sets the env). The suite is the
  contract; do not pin every e2e gate.

- [ ] **Step 3: Add the PROGRESS.md changelog entry**

Prepend under the `## Changelog` section of `PROGRESS.md`:

```markdown
- 2026-06-16: **Backend auto-selection + informative compile output.** The C++ Jacobian
  backend is now chosen at codegen time by `gdsge.codegen.resolveBackend` with precedence
  in-gmod `UseAutoDiff` > `GDSGE_BACKEND` env (`adept|sympy|auto`) > auto-detect (Python
  present ⇒ SymPy, else adept). `UseAutoDiff` became tri-state in `resolveOptions` (unset ⇒
  no IR field/auto; =1 ⇒ `'autodiff'`; =0 ⇒ `'sympy'`), so default-path IR stays
  byte-identical (no corpus gmod sets it). In **auto mode only**, a SymPy *codegen* failure
  falls back to adept (`generateCxxWithFallback`); explicit/env pins error instead.
  `generateCxx` gained an `opts.backend` override; the kernel-ensure branch was extracted to
  `gdsge.codegen.ensureKernels`. The driver prints a one-line `Backend: …` message every run
  and numbered `[i/5]` phase banners (parsing → MATLAB → kernels → C++ → solver MEX), and
  `compile_<model>.m` now labels the simulate-MEX step. The test harness pins
  `GDSGE_BACKEND=adept` to keep the correctness loop Python-free; sympy gates opt in via
  `UseAutoDiff=0`. Spec/plan: `docs/superpowers/{specs,plans}/2026-06-16-backend-auto-selection*`.
  Branch `feature/backend-auto-selection`.
```

- [ ] **Step 4: Commit**

```bash
git add PROGRESS.md
git commit -m "docs(progress): record backend auto-selection + informative compile output"
```

---

## Self-Review

**Spec coverage:**
- Tri-state `UseAutoDiff` (spec §1) → Task 1.
- `resolveBackend` precedence + `badBackendEnv` (spec §2) → Task 2.
- Backend-choice messages (spec §3, messaging catalogue) → Task 3 (message) + Task 8 (printed every run) + Task 6 (fallback line).
- `generateCxx` `opts.backend` override (spec §3 arch) → Task 4.
- Codegen-time-only fallback, auto-only, explicit errors (spec §2, §4) → Task 6.
- Phase banners incl. kernels/solver/simulate labels (spec §4, §5) → Task 5 (ensureKernels) + Task 7 (simulate label) + Task 8 (banners).
- Harness force-adept + sympy gates opt in (spec §5, test plan) → Task 9.
- Tests: resolveBackend unit, fallback, explicit-errors, messaging/banner asserts, harness (spec test plan) → Tasks 2,3,6,8,9.

**Placeholder scan:** none — every code/edit step shows the exact content; run commands and expected results are concrete.

**Type/name consistency:** backend token is `'autodiff'`/`'sympy'` throughout (matches the `jacobianBackend` enum); decision struct fields `backend`/`mode`/`reason` are identical across `resolveBackend`, `backendMessage`, and `generateCxxWithFallback`; `GDSGE_BACKEND` and the `gdsge:codegen:badBackendEnv` id are used consistently; `ensureKernels`/`generateCxxWithFallback`/`backendMessage`/`resolveBackend` names match their call sites in `codegen.m`.

**Note on the fallback test (Task 6):** no corpus construct fails sympy-only while
succeeding under adept (pchip and cxx-hooks fail under *both*), so the fallback is verified
with an injected `genFn` rather than a real model — deterministic and Python-free.
```
