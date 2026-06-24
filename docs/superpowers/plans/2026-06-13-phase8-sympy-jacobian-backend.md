# Phase 8 — SymPy Analytic-Jacobian Backend Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an alternative C++ MEX backend that computes the model Jacobian analytically via SymPy + CSE (stack-allocated `double helper_i;`) instead of adept autodiff, for the four spline-interpolation corpus models, with the default autodiff path producing byte-identical output.

**Architecture:** Approach C (spec §5). The MEX outer structure is unchanged (POPN stack locals, OpenMP per-grid-point, CoDoSol). The model function is regenerated to compute residuals **and** the analytic Jacobian in one all-`double` pass. Reductions lower to a single fused shock loop accumulating both value and gradient row (spec §10). MATLAB owns the loop structure and an IR→C++ *gradient registry* (sparse `slot → C++ expr` rows, forward-mode chain rule); Python/SymPy differentiates each body and returns value + partials as CSE'd C++. The chain rule closes through the spline kernel's pure-`double` `search_eval_with_grad_vec_at_array` (`interp'`). Selection is a codegen-time IR option emitted only when non-default, so no goldens regenerate.

**Tech Stack:** MATLAB R2025b (`+gdsge` package, `matlab.unittest`), MATLAB↔Python via `pyenv` (in-process, JSON string in/out), uv-managed Python 3.12 + SymPy (`pyext/`), C++ MEX (MSVC), `pytest`.

---

## Key facts the implementer must know (read before starting)

- **Spec:** `docs/superpowers/specs/2026-06-13-phase8-sympy-jacobian-backend-design.md`. **Branch:** `phase8-sympy-jacobian-backend` (already created).
- **Run MATLAB tests:** `matlab -batch "cd('tests'); run_tests"` (exit 0 = all pass). PowerShell 7 is **not** installed on this machine — do **not** use `tests/run.ps1`. Authoritative result = exit code + `tests/results/junit.xml` (`results.tap` accumulates stale lines across runs — ignore it). MATLAB is at `C:\Program Files\MATLAB\R2025b\bin\matlab.exe`.
- **Run a single MATLAB test class** (faster loop): `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('tFoo'); disp(r); exit(any([r.Failed]))"`. The suite's `run_tests.m` adds `../src` + kernels to the path itself; for single-class runs add `../src` yourself.
- **Run Python tests:** `uv run --project pyext pytest pyext/tests -q`.
- **MATLAB path golden rule** (CLAUDE.md): never add a toolbox source to a persistent path; one source per process. Tests use `PathFixture`/`WorkingFolderFixture`.
- **Slot layout** (the registry's foundation): each `ir.variables.policy{i}` has `.name`, `.length`, `.slot=[start stop]` (1-based, contiguous, covering `1..NX`). Scalar policy `c1` → unknown index `slot(1)-1` in `GDSGE_x`; array policy `w1n` (length 8) → `GDSGE_x[slot(1)-1 .. slot(2)-1]`, element `j` (1-based) → index `slot(1)-1 + (j-1)`. `NX == NUM_EQUATIONS` (square; `emitEquations` asserts).
- **Interp grad layout** (`MatlabInterpEval.h:41`): `search_eval_with_grad_vec_at_array(arrayIdx, xSite, evalResult, grad, cellOfSite)` fills `grad[i_vec + numVec*i_dim] = d(evalResult[i_vec])/d(xSite[i_dim])`. `numVec = #interp vars`, `i_dim` over interp args. HL1996: `xdim=1`, so `grad[interpIdx]` is `d(interp_interpIdx)/d(state)`.
- **The autodiff path is the reference and must stay byte-identical.** Do not edit `emitModel.m`, `model.tpl.cpp`, the adept block of `model_eval.tpl.cpp`, `emitModelBody.m`, `emitExpr.m`. The SymPy backend adds **new** emitters/templates and a **branch** in `generateCxx`.
- **gmod surface flag:** `UseAutoDiff` (numeric, default `1`). `UseAutoDiff=0;` in the gmod selects the SymPy backend. Resolved in `resolveOptions`; the IR field `options.jacobianBackend='sympy'` is emitted **only** when non-default (mirrors how ASG levels are omitted for spline), so existing model IR JSON snapshots stay byte-identical.

---

## File structure

**New files:**
- `src/+gdsge/+codegen/+sympy/ensurePyenv.m` — configure/verify the pyenv bridge.
- `src/+gdsge/+codegen/+sympy/callSympy.m` — JSON-in/JSON-out chokepoint to Python.
- `src/+gdsge/+codegen/+cxx/emitModelSympy.m` — the gradient registry + fused-loop C++ emitter (analytic Jacobian model function).
- `src/+gdsge/+codegen/+cxx/emitInterpSympy.m` — double with-grad spline interp evaluator section.
- `templates/cxx/model_sympy.tpl.cpp` — the all-`double` value+Jacobian model-function template.
- `pyext/gdsge_sympy/__init__.py` — Python entry functions (JSON string in/out).
- `pyext/gdsge_sympy/ast_to_sympy.py` — IR node → SymPy expression.
- `pyext/gdsge_sympy/ccode_helpers.py` — `my_ccode` recipe (hans-derived).
- `pyext/gdsge_sympy/diff_body.py` — differentiate one body w.r.t. its free symbols + CSE.
- `pyext/tests/test_ast_to_sympy.py`, `pyext/tests/test_diff_body.py` — pytest.
- `tests/+gdsgetest/sympyAvailable.m` — test helper wrapping `ensurePyenv` for `assumeTrue`.
- `tests/codegen/tJacobianBackendOption.m` — IR option + guard tests.
- `tests/codegen/tSympyBridge.m` — `ensurePyenv`/`callSympy` round-trip.
- `tests/codegen/tEmitModelSympy.m` — registry/emitter unit + micro-compile FD check.
- `tests/HeatonLucas1996/codegen/tSympyJacobianHL1996.m` — three-way Jacobian cross-check.
- `tests/HeatonLucas1996/codegen/tEndToEndHL1996Sympy.m` — end-to-end golden, sympy backend.
- `tests/Barro_et_al_2017/codegen/tEndToEndSafeAssetsSympy.m`,
  `tests/Mendoza2010/codegen/tEndToEndMendoza2010Sympy.m`,
  `tests/GLSW2020/codegen/tEndToEndGLSWSympy.m` — widen.

**Modified files:**
- `src/+gdsge/+ir/schema.m` — add optional `options.jacobianBackend` enum.
- `src/+gdsge/+ir/validate.m` — guard invariant (sympy ⇒ spline, no cxx/tensor).
- `src/+gdsge/+parser/resolveOptions.m` — map `UseAutoDiff` → `jacobianBackend` (emit only when sympy).
- `src/+gdsge/+codegen/generateCxx.m` — branch on `jacobianBackend`; guards.
- `templates/cxx/task.tpl.cpp` — additive 6th output (analytic Jacobian) for the cross-check, `GDSGE_DEBUG_EVAL_ONLY==2`.
- `templates/cxx/call_fmin.tpl.cpp` — fill `plhs[5]` jac when present in eval mode.
- `PROGRESS.md` — Phase 8 status + changelog.

---

## Task 1: IR option `jacobianBackend` (schema + validator)

**Files:**
- Modify: `src/+gdsge/+ir/schema.m:40-56` (the `options` fStruct)
- Modify: `src/+gdsge/+ir/validate.m:66-81` (`checkOptionInvariants`)
- Test: `tests/codegen/tJacobianBackendOption.m`

- [ ] **Step 1: Write the failing test**

Create `tests/codegen/tJacobianBackendOption.m`:

```matlab
classdef tJacobianBackendOption < matlab.unittest.TestCase
    % PHASE 8: the IR carries an optional jacobianBackend enum; validate
    % accepts 'autodiff'/'sympy'/absent and rejects sympy+asg / sympy+cxx /
    % sympy+tensor. Default-absent keeps existing IR byte-identical.
    methods (TestClassSetup)
        function addSrc(tc)
            here = fileparts(mfilename('fullpath'));        % tests/codegen
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src')));
        end
    end
    methods (Test)
        function acceptsSympyEnum(tc)
            ir = baseSplineIR();
            ir.options.jacobianBackend = 'sympy';
            r = gdsge.ir.validate(ir);
            tc.verifyTrue(r.pass, strjoin(r.errors, newline));
        end
        function acceptsAbsent(tc)
            ir = baseSplineIR();              % no jacobianBackend field
            r = gdsge.ir.validate(ir);
            tc.verifyTrue(r.pass, strjoin(r.errors, newline));
        end
        function rejectsBadEnum(tc)
            ir = baseSplineIR();
            ir.options.jacobianBackend = 'finitediff';
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
        end
        function rejectsSympyWithAsg(tc)
            ir = baseSplineIR();
            ir.options.interpMethod = 'asg';
            ir.options.asgMaxLevel = 10;
            ir.options.jacobianBackend = 'sympy';
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
            tc.verifySubstring(strjoin(r.errors, ' '), 'jacobianBackend');
        end
    end
end

function ir = baseSplineIR()
    % Minimal schema-valid spline IR (one state, one policy, one shock).
    ir = struct();
    ir.irVersion = '1.0.0';
    ir.modelName = 'tiny';
    ir.params = {};
    ir.options = struct('interpMethod','spline','interpOrder',4);
    ir.shocks = struct('names',{{'z'}},'count',1, ...
        'values',struct('z',1), 'transitions',struct('shock_trans',1));
    ir.states = struct('names',{{'k'}}, 'grids',struct('k','linspace(1,3,3)'));
    ir.variables = struct('policy',{{struct('name','c','length',1,'slot',[1 1])}}, ...
        'aux',{{}}, 'interp',{{}}, 'tensor',{{}}, 'output',{{}}, 'others',{{}});
    ir.bounds = {struct('name','c','lower','0','upper','1')};
    ir.interp = {};
    ir.model = struct('statements',{{}}, ...
        'equations',{{struct('expr',gdsge.ir.node.num(0),'primed',0)}});
    ir.simulate = struct('numPeriods',1,'numSamples',1,'initial',{{}}, ...
        'varSimu',{{}},'transitions',{{}});
    ir.hooks = struct();
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); r=runtests('codegen/tJacobianBackendOption'); exit(any([r.Failed]))"`
Expected: FAIL — `acceptsSympyEnum`/`rejectsBadEnum` fail because the schema has no `jacobianBackend` field (it is silently ignored, so `rejectsBadEnum` returns pass=true).

- [ ] **Step 3: Add the schema field**

In `src/+gdsge/+ir/schema.m`, inside the `options = fStruct(structOf(...))` block (after `'interpMethod', fEnum(...)`), add:

```matlab
    'jacobianBackend', opt(fEnum({'autodiff','sympy'})), ...
```

- [ ] **Step 4: Add the validator invariant**

In `src/+gdsge/+ir/validate.m`, in `checkOptionInvariants()`, after the existing asg block, add:

```matlab
        if isfield(o, 'jacobianBackend') && strcmp(o.jacobianBackend, 'sympy')
            if isfield(o,'interpMethod') && ~strcmp(o.interpMethod, 'spline')
                err('options.jacobianBackend', ...
                    'sympy backend supports only interpMethod = spline');
            end
            if isfield(ir,'variables') && isfield(ir.variables,'tensor') ...
                    && iscell(ir.variables.tensor) && ~isempty(ir.variables.tensor)
                err('options.jacobianBackend', ...
                    'sympy backend cannot be used with var_tensor');
            end
            if isfield(ir,'hooks') && isfield(ir.hooks,'cxx') ...
                    && ischar(ir.hooks.cxx) && ~isempty(ir.hooks.cxx)
                err('options.jacobianBackend', ...
                    'sympy backend cannot differentiate cxx blocks');
            end
        end
```

- [ ] **Step 5: Run test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); r=runtests('codegen/tJacobianBackendOption'); exit(any([r.Failed]))"`
Expected: PASS (4/4).

- [ ] **Step 6: Confirm no schema doc drift**

Run the schema-doc no-drift test (regenerates `docs/ir-schema.md`):
Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); r=runtests('ir'); exit(any([r.Failed]))"`
Expected: PASS. If a doc-drift test fails, regenerate per its message (`gdsge.ir.gendoc`) and stage the updated `docs/ir-schema.md`.

- [ ] **Step 7: Commit**

```bash
git add src/+gdsge/+ir/schema.m src/+gdsge/+ir/validate.m tests/codegen/tJacobianBackendOption.m docs/ir-schema.md
git commit -m "feat(ir): optional options.jacobianBackend enum + sympy invariants"
```

---

## Task 2: Resolve `UseAutoDiff` → `jacobianBackend` (parser)

**Files:**
- Modify: `src/+gdsge/+parser/resolveOptions.m`
- Test: add a method to `tests/codegen/tJacobianBackendOption.m`

- [ ] **Step 1: Write the failing test**

Add to `tJacobianBackendOption.m` methods (Test) block:

```matlab
        function resolveOmitsWhenAutodiff(tc)
            o = gdsge.parser.resolveOptions(struct());     % no UseAutoDiff
            tc.verifyFalse(isfield(o, 'jacobianBackend'), ...
                'default must omit jacobianBackend (keeps goldens byte-identical)');
        end
        function resolveSetsSympyWhenFlagOff(tc)
            o = gdsge.parser.resolveOptions(struct('UseAutoDiff', 0));
            tc.verifyEqual(o.jacobianBackend, 'sympy');
        end
        function resolveOmitsWhenFlagOn(tc)
            o = gdsge.parser.resolveOptions(struct('UseAutoDiff', 1));
            tc.verifyFalse(isfield(o, 'jacobianBackend'));
        end
```

- [ ] **Step 2: Run to verify fail**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); r=runtests('codegen/tJacobianBackendOption'); exit(any([r.Failed]))"`
Expected: FAIL — `resolveSetsSympyWhenFlagOff` errors (no such field).

- [ ] **Step 3: Implement**

In `src/+gdsge/+parser/resolveOptions.m`, before the final `end` of the main function (after `o.saveFreq = ...`), add:

```matlab
% Jacobian backend: default autodiff (no Python). Emit the field ONLY when the
% user opts into the SymPy analytic backend, so default-path IR stays
% byte-identical to existing goldens (cf. ASG levels omitted for spline).
if getf(ws, 'UseAutoDiff', 1) == 0
    o.jacobianBackend = 'sympy';
end
```

- [ ] **Step 4: Run to verify pass**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); r=runtests('codegen/tJacobianBackendOption'); exit(any([r.Failed]))"`
Expected: PASS (7/7).

- [ ] **Step 5: Regression — existing front-end goldens unchanged**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); r=runtests('parser'); exit(any([r.Failed]))"`
Expected: PASS (no IR snapshot changed — field is omitted by default).

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+parser/resolveOptions.m tests/codegen/tJacobianBackendOption.m
git commit -m "feat(parser): UseAutoDiff=0 -> options.jacobianBackend='sympy'"
```

---

## Task 3: pyenv bridge — `ensurePyenv` + `callSympy`

**Files:**
- Create: `src/+gdsge/+codegen/+sympy/ensurePyenv.m`
- Create: `src/+gdsge/+codegen/+sympy/callSympy.m`
- Create: `pyext/gdsge_sympy/__init__.py`
- Create: `tests/+gdsgetest/sympyAvailable.m`
- Test: `tests/codegen/tSympyBridge.m`

- [ ] **Step 1: Write a minimal Python entry so the bridge has something to call**

Create `pyext/gdsge_sympy/__init__.py`:

```python
"""SymPy analytic-Jacobian code generation for the GDSGE alt backend.

Entry points take a single JSON string and return a single JSON string, so the
MATLAB<->Python boundary is one char in / one char out (see callSympy.m).
"""
import json


def ping(payload: str) -> str:
    """Round-trip sanity check used by the bridge test."""
    obj = json.loads(payload)
    return json.dumps({"ok": True, "echo": obj})
```

- [ ] **Step 2: Write `ensurePyenv`**

Create `src/+gdsge/+codegen/+sympy/ensurePyenv.m`:

```matlab
function ok = ensurePyenv()
% ENSUREPYENV  Point MATLAB's pyenv at the uv-managed venv interpreter and make
%   the gdsge_sympy package importable. Idempotent; callable from matlab -batch.
%   Returns true iff `import gdsge_sympy` and `import sympy` both succeed.
%   The SymPy backend is the ONLY consumer; the default autodiff path never
%   calls this, so a missing/unsynced venv only affects sympy requests.
ok = false;
try
    here = fileparts(mfilename('fullpath'));                 % src/+gdsge/+codegen/+sympy
    repoRoot = fileparts(fileparts(fileparts(fileparts(here))));
    pyextDir = fullfile(repoRoot, 'pyext');
    venvPy = fullfile(pyextDir, '.venv', 'Scripts', 'python.exe');  % Windows
    if exist(venvPy, 'file') ~= 2
        return;     % env not synced: caller raises the informative guard error
    end
    pe = pyenv;
    if ~strcmp(pe.Executable, venvPy)
        if pe.Status == "Loaded"
            % Cannot switch interpreter once Python is loaded in this process;
            % accept the loaded one if it already has the package (below).
        else
            pyenv('Version', venvPy);
        end
    end
    % Make pyext/ importable (the package lives at pyext/gdsge_sympy).
    if count(py.sys.path, pyextDir) == 0
        insert(py.sys.path, int32(0), pyextDir);
    end
    py.importlib.import_module('sympy');
    py.importlib.import_module('gdsge_sympy');
    ok = true;
catch
    ok = false;
end
end
```

- [ ] **Step 3: Write `callSympy`**

Create `src/+gdsge/+codegen/+sympy/callSympy.m`:

```matlab
function out = callSympy(fnName, requestStruct)
% CALLSYMPY  The single MATLAB->SymPy chokepoint. jsonencode the request,
%   invoke py.gdsge_sympy.<fnName>(jsonChar), jsondecode the returned JSON.
%   Assumes gdsge.codegen.sympy.ensurePyenv() has already returned true.
reqJson = jsonencode(requestStruct);
pyFn = py.getattr(py.importlib.import_module('gdsge_sympy'), fnName);
resJson = char(pyFn(reqJson));
out = jsondecode(resJson);
end
```

- [ ] **Step 4: Write the test helper**

Create `tests/+gdsgetest/sympyAvailable.m`:

```matlab
function tf = sympyAvailable()
% SYMPYAVAILABLE  True iff the uv venv + gdsge_sympy import succeed. SymPy
%   gates call assumeTrue(gdsgetest.sympyAvailable()) so the default
%   (Python-free) correctness loop stays green.
tf = gdsge.codegen.sympy.ensurePyenv();
end
```

- [ ] **Step 5: Write the bridge test**

Create `tests/codegen/tSympyBridge.m`:

```matlab
classdef tSympyBridge < matlab.unittest.TestCase
    % PHASE 8: the pyenv bridge round-trips a struct through JSON to Python and
    % back. Skipped when the uv venv is unsynced (default Python-free machine).
    methods (TestClassSetup)
        function addSrc(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '+gdsgetest')));
            tc.assumeTrue(gdsgetest.sympyAvailable(), ...
                'uv venv not synced (run: uv sync --project pyext)');
        end
    end
    methods (Test)
        function roundTripsStruct(tc)
            out = gdsge.codegen.sympy.callSympy('ping', struct('a', 1, 'b', 'x'));
            tc.verifyTrue(logical(out.ok));
            tc.verifyEqual(out.echo.a, 1);
            tc.verifyEqual(char(out.echo.b), 'x');
        end
    end
end
```

- [ ] **Step 6: Sync the venv and run**

Run: `uv sync --project pyext`
Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('codegen/tSympyBridge'); exit(any([r.Failed]))"`
Expected: PASS (1/1). If the venv cannot be synced in this environment, the test is *skipped* (Incomplete), not failed — verify it reports Incomplete, not Failed.

- [ ] **Step 7: Commit**

```bash
git add src/+gdsge/+codegen/+sympy/ pyext/gdsge_sympy/__init__.py tests/+gdsgetest/sympyAvailable.m tests/codegen/tSympyBridge.m
git commit -m "feat(codegen): pyenv bridge (ensurePyenv + callSympy) + gdsge_sympy.ping"
```

---

## Task 4: Python — IR node → SymPy (`ast_to_sympy`)

**Files:**
- Create: `pyext/gdsge_sympy/ast_to_sympy.py`
- Test: `pyext/tests/test_ast_to_sympy.py`

The IR node JSON shapes (from `schema.m`): `num{kind,value}`, `name{kind,id}`, `primed{kind,id}`, `unop{kind,op,arg}`, `binop{kind,op,lhs,rhs}`, `call{kind,fn,args[]}`. Primed names and plain names both become SymPy symbols; the caller decides their role (see Task 6). Built-in calls map to SymPy functions; `pow` and `^` map to `**`.

- [ ] **Step 1: Write failing pytest**

Create `pyext/tests/test_ast_to_sympy.py`:

```python
import sympy
from gdsge_sympy.ast_to_sympy import to_sympy, symbol_name


def num(v):  return {"kind": "num", "value": v}
def name(i): return {"kind": "name", "id": i}
def primed(i): return {"kind": "primed", "id": i}
def binop(op, a, b): return {"kind": "binop", "op": op, "lhs": a, "rhs": b}
def unop(op, a): return {"kind": "unop", "op": op, "arg": a}
def call(fn, *a): return {"kind": "call", "fn": fn, "args": list(a)}


def test_basic_arithmetic():
    syms = {}
    e = to_sympy(binop("+", name("a"), binop("*", num(2), name("b"))), syms)
    a, b = sympy.symbols("a b")
    assert sympy.simplify(e - (a + 2 * b)) == 0


def test_power_maps_to_pow():
    syms = {}
    e = to_sympy(binop("^", name("x"), num(3)), syms)
    (x,) = sympy.symbols("x"),
    assert sympy.simplify(e - x[0] ** 3) == 0


def test_unary_minus():
    syms = {}
    e = to_sympy(unop("-", name("x")), syms)
    assert sympy.simplify(e + sympy.Symbol("x")) == 0


def test_builtin_call_exp_log():
    syms = {}
    e = to_sympy(call("exp", call("log", name("x"))), syms)
    assert sympy.simplify(e - sympy.Symbol("x")) == 0


def test_primed_distinct_symbol():
    syms = {}
    e = to_sympy(binop("+", name("g"), primed("g")), syms)
    # name g and primed g are DISTINCT symbols (g vs g__prime)
    assert symbol_name(name("g")) != symbol_name(primed("g"))
    assert len(e.free_symbols) == 2
```

- [ ] **Step 2: Run to verify fail**

Run: `uv run --project pyext pytest pyext/tests/test_ast_to_sympy.py -q`
Expected: FAIL — module `gdsge_sympy.ast_to_sympy` does not exist.

- [ ] **Step 3: Implement**

Create `pyext/gdsge_sympy/ast_to_sympy.py`:

```python
"""IR expression node (dict) -> SymPy expression. Inverse of emitExpr.m.

Each `name`/`primed` node becomes a SymPy Symbol; the caller supplies a `syms`
dict that maps the canonical symbol name to a sympy.Symbol so the same name
reuses one symbol. Primed names get a distinct '__prime' suffix so a current
unknown `g` and its shock realization `g'` never collide.
"""
import sympy

# gmod built-ins -> sympy. Anything not here raises (caught upstream).
_BUILTINS = {
    "exp": sympy.exp, "log": sympy.log, "sqrt": sympy.sqrt,
    "abs": sympy.Abs, "sin": sympy.sin, "cos": sympy.cos, "tan": sympy.tan,
    "max": sympy.Max, "min": sympy.Min,
}


def symbol_name(node):
    if node["kind"] == "name":
        return node["id"]
    if node["kind"] == "primed":
        return node["id"] + "__prime"
    raise ValueError("symbol_name expects a name/primed node")


def _sym(node, syms):
    nm = symbol_name(node)
    if nm not in syms:
        syms[nm] = sympy.Symbol(nm, real=True)
    return syms[nm]


def to_sympy(node, syms):
    k = node["kind"]
    if k == "num":
        return sympy.Float(node["value"]) if node["value"] != int(node["value"]) \
            else sympy.Integer(int(node["value"]))
    if k in ("name", "primed"):
        return _sym(node, syms)
    if k == "unop":
        if node["op"] != "-":
            raise ValueError(f"unsupported unop {node['op']}")
        return -to_sympy(node["arg"], syms)
    if k == "binop":
        lhs = to_sympy(node["lhs"], syms)
        rhs = to_sympy(node["rhs"], syms)
        op = node["op"]
        if op == "+": return lhs + rhs
        if op == "-": return lhs - rhs
        if op == "*": return lhs * rhs
        if op == "/": return lhs / rhs
        if op == "^": return lhs ** rhs
        raise ValueError(f"unsupported binop {op}")
    if k == "call":
        fn = node["fn"]
        if fn == "pow":
            a, b = (to_sympy(x, syms) for x in node["args"])
            return a ** b
        if fn in _BUILTINS:
            return _BUILTINS[fn](*[to_sympy(x, syms) for x in node["args"]])
        raise ValueError(f"unsupported call {fn}")
    raise ValueError(f"unknown node kind {k}")
```

- [ ] **Step 4: Run to verify pass**

Run: `uv run --project pyext pytest pyext/tests/test_ast_to_sympy.py -q`
Expected: PASS (5/5).

- [ ] **Step 5: Commit**

```bash
git add pyext/gdsge_sympy/ast_to_sympy.py pyext/tests/test_ast_to_sympy.py
git commit -m "feat(sympy): IR AST node -> SymPy translation"
```

---

## Task 5: Python — `ccode_helpers` + `diff_body`

**Files:**
- Create: `pyext/gdsge_sympy/ccode_helpers.py`
- Create: `pyext/gdsge_sympy/diff_body.py`
- Modify: `pyext/gdsge_sympy/__init__.py` (export `diff_body`)
- Test: `pyext/tests/test_diff_body.py`

**Contract of `diff_body`** (the one thing Python does, spec §5.3): given a body AST + the list of free-symbol names to differentiate by, return a JSON dict:
- `helpers`: list of `[lhs, rhs_ccode]` CSE temporaries (in order; `lhs` like `helper_0`).
- `value`: C++ expr string for the body value (referencing helpers).
- `partials`: map `symbolName -> C++ expr string` for `d(body)/d(symbol)` (referencing helpers; omitted when structurally zero).

The value and all partials are CSE'd **together** (one `cse(...)` call) so they share `helper_i` temporaries — the locality that matters.

- [ ] **Step 1: Write failing pytest**

Create `pyext/tests/test_diff_body.py`:

```python
import json
import math
from gdsge_sympy import diff_body_json


def binop(op, a, b): return {"kind": "binop", "op": op, "lhs": a, "rhs": b}
def name(i): return {"kind": "name", "id": i}
def num(v): return {"kind": "num", "value": v}
def call(fn, *a): return {"kind": "call", "fn": fn, "args": list(a)}


def _eval_cpp(res, env):
    """Tiny evaluator: run helpers then a final expr, in Python, to check ccode
    is numerically correct. Replaces pow(a,b)->a**b and known funcs."""
    import re
    scope = dict(env)
    def fix(s):
        s = re.sub(r"\bpow\(", "__pow(", s)
        return s
    def __pow(a, b): return a ** b
    scope["__pow"] = __pow
    scope["exp"] = math.exp
    scope["log"] = math.log
    for lhs, rhs in res["helpers"]:
        scope[lhs] = eval(fix(rhs), {"__builtins__": {}}, scope)
    return scope, fix


def test_value_and_partial_simple():
    # body = a*b ; d/da = b ; d/db = a
    body = binop("*", name("a"), name("b"))
    res = json.loads(diff_body_json(json.dumps(
        {"body": body, "diffVars": ["a", "b"]})))
    scope, fix = _eval_cpp(res, {"a": 3.0, "b": 5.0})
    assert abs(eval(fix(res["value"]), {"__builtins__": {}}, scope) - 15.0) < 1e-12
    assert abs(eval(fix(res["partials"]["a"]), {"__builtins__": {}}, scope) - 5.0) < 1e-12
    assert abs(eval(fix(res["partials"]["b"]), {"__builtins__": {}}, scope) - 3.0) < 1e-12


def test_zero_partial_omitted():
    # body = a + 7 ; d/db = 0 -> omitted
    body = binop("+", name("a"), num(7))
    res = json.loads(diff_body_json(json.dumps(
        {"body": body, "diffVars": ["a", "b"]})))
    assert "b" not in res["partials"]
    assert "a" in res["partials"]


def test_power_partial():
    # body = a^3 ; d/da = 3*a^2
    body = binop("^", name("a"), num(3))
    res = json.loads(diff_body_json(json.dumps(
        {"body": body, "diffVars": ["a"]})))
    scope, fix = _eval_cpp(res, {"a": 2.0})
    assert abs(eval(fix(res["partials"]["a"]), {"__builtins__": {}}, scope) - 12.0) < 1e-9
```

- [ ] **Step 2: Run to verify fail**

Run: `uv run --project pyext pytest pyext/tests/test_diff_body.py -q`
Expected: FAIL — `diff_body_json` not importable.

- [ ] **Step 3: Implement `ccode_helpers`**

Create `pyext/gdsge_sympy/ccode_helpers.py`:

```python
"""ccode emission tuned for the GDSGE C++ templates (hans recipe, generalized).
pow() is emitted for ** (matches emitExpr.m). Symbol names pass through verbatim
so they bind to the C++ locals MATLAB declares for them."""
import sympy
from sympy.printing.c import C99CodePrinter


class _GdsgePrinter(C99CodePrinter):
    def _print_Pow(self, expr):
        if expr.exp == -1:
            return f"1.0/({self._print(expr.base)})"
        return f"pow({self._print(expr.base)},{self._print(expr.exp)})"


_PRINTER = _GdsgePrinter()


def ccode(expr):
    return _PRINTER.doprint(expr)
```

- [ ] **Step 4: Implement `diff_body`**

Create `pyext/gdsge_sympy/diff_body.py`:

```python
"""Differentiate one body expression w.r.t. a set of free symbols, CSE the value
and all partials together, and emit C++ (helpers + value + partials). Spec §5.3.
"""
import json
import sympy
from sympy.simplify.cse_main import cse
from sympy.utilities.iterables import numbered_symbols
from .ast_to_sympy import to_sympy
from .ccode_helpers import ccode


def diff_body(body_node, diff_var_names):
    syms = {}
    value_expr = to_sympy(body_node, syms)

    # Differentiation targets: only those that actually appear in the body.
    targets = []
    for nm in diff_var_names:
        s = syms.get(nm)
        if s is not None:
            targets.append((nm, s))

    partial_exprs = []
    partial_keys = []
    for nm, s in targets:
        d = sympy.diff(value_expr, s)
        if d != 0:
            partial_keys.append(nm)
            partial_exprs.append(d)

    # One CSE over value + all nonzero partials so they share helpers.
    all_exprs = [value_expr] + partial_exprs
    replacements, reduced = cse(all_exprs, symbols=numbered_symbols("helper_"))

    helpers = [[str(lhs), ccode(rhs)] for lhs, rhs in replacements]
    value_cpp = ccode(reduced[0])
    partials = {nm: ccode(reduced[1 + i]) for i, nm in enumerate(partial_keys)}
    return {"helpers": helpers, "value": value_cpp, "partials": partials}


def diff_body_json(payload: str) -> str:
    req = json.loads(payload)
    return json.dumps(diff_body(req["body"], req["diffVars"]))
```

- [ ] **Step 5: Export from the package**

In `pyext/gdsge_sympy/__init__.py`, add at the bottom:

```python
from .diff_body import diff_body_json  # noqa: E402,F401
```

- [ ] **Step 6: Run to verify pass**

Run: `uv run --project pyext pytest pyext/tests -q`
Expected: PASS (all). 

- [ ] **Step 7: Commit**

```bash
git add pyext/gdsge_sympy/ccode_helpers.py pyext/gdsge_sympy/diff_body.py pyext/gdsge_sympy/__init__.py pyext/tests/test_diff_body.py
git commit -m "feat(sympy): diff_body (value+partials, shared CSE) + ccode helpers"
```

---

## Task 6: Gradient registry + `emitModelSympy` — flat (no reductions) slice

This task builds the MATLAB orchestrator for the *simplest* model shape: unknowns + plain (non-primed) assigns + non-primed equations, no reductions, no interp. It establishes the registry data structure, the chain-rule combiner, and the C++ assembly, verified by a micro model. Reductions (Task 7) and interp (Task 8) extend it.

**Files:**
- Create: `src/+gdsge/+codegen/+cxx/emitModelSympy.m`
- Test: `tests/codegen/tEmitModelSympy.m`

**Registry representation (MATLAB):** a `containers.Map` from name → `gradRow`, where a `gradRow` is a struct array of `{slot:int, expr:char}` entries (sparse; only nonzero slots). Helpers:
- `seedUnknownRow(slotStart0, len, k)` → row `{slot:slotStart0+k, expr:'1'}` (the k-th element).
- `combine(partials, freeSymRows)` → for the new name, `row[m] = Σ_s partials(s) · freeSymRows(s)[m]`, emitted as C++ sums per nonzero slot.

The emitter walks statements, calls `callSympy('diff_body_json', …)` per RHS, and accumulates C++ into a `codeWriter`. The final model function (filled into `model_sympy.tpl.cpp` in Task 9) computes `GDSGE_f[i]` and, guarded by `if (GDSGE_jac)`, `JAC(i,m)`.

- [ ] **Step 1: Write the failing test (registry seed + combine, pure MATLAB)**

Create `tests/codegen/tEmitModelSympy.m`:

```matlab
classdef tEmitModelSympy < matlab.unittest.TestCase
    % PHASE 8: the gradient registry and emitModelSympy. Pure-MATLAB unit tests
    % for the registry algebra; the numeric correctness of the emitted C++ is
    % pinned by the micro-compile FD check (last method, Slow) and by the
    % per-model three-way Jacobian gates.
    methods (TestClassSetup)
        function addSrc(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '+gdsgetest')));
        end
    end
    methods (Test)
        function seedRowIsUnitAtSlot(tc)
            row = gdsge.codegen.cxx.sympymodel.seedUnknownRow(3, 1, 0);
            tc.verifyEqual(numel(row), 1);
            tc.verifyEqual(row(1).slot, 3);
            tc.verifyEqual(row(1).expr, '1');
        end
        function combineAppliesChainRule(tc)
            % new = f(a,b), with d/da, d/db given; a has row {0:'1'},
            % b has row {2:'g_b'}. Expect new row at slots {0,2}.
            partials = containers.Map({'a','b'}, {'pa','pb'});
            rows = containers.Map();
            rows('a') = struct('slot', 0, 'expr', '1');
            rows('b') = struct('slot', 2, 'expr', 'g_b');
            row = gdsge.codegen.cxx.sympymodel.combine(partials, rows);
            slots = sort([row.slot]);
            tc.verifyEqual(slots, [0 2]);
        end
    end
end
```

- [ ] **Step 2: Run to verify fail**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); r=runtests('codegen/tEmitModelSympy'); exit(any([r.Failed]))"`
Expected: FAIL — `emitModelSympy` package/functions do not exist.

- [ ] **Step 3: Implement the registry primitives**

The SymPy emitter is a MATLAB package: **create directory** `src/+gdsge/+codegen/+cxx/+sympymodel/` and put each function in its own file (called as `gdsge.codegen.cxx.sympymodel.<fn>(...)`). For this first slice create the two registry primitives.

`src/+gdsge/+codegen/+cxx/+sympymodel/seedUnknownRow.m`:

```matlab
function row = seedUnknownRow(slotStart0, len, k)
% SEEDUNKNOWNROW  Gradient row of the k-th (0-based) element of an unknown
%   occupying [slotStart0 .. slotStart0+len-1] in GDSGE_x: unit at its own slot.
assert(k >= 0 && k < len);
row = struct('slot', slotStart0 + k, 'expr', '1');
end
```

`src/+gdsge/+codegen/+cxx/+sympymodel/combine.m`:

```matlab
function row = combine(partials, freeSymRows)
% COMBINE  Chain rule: row[m] = sum_s partials(s) * freeSymRows(s)[m].
%   partials: containers.Map symbolName -> C++ partial expr (only nonzero ones).
%   freeSymRows: containers.Map symbolName -> gradRow (struct array slot/expr).
%   Returns a consolidated gradRow (struct array), terms per slot summed as a
%   C++ sum string. Symbols absent from freeSymRows contribute nothing.
acc = containers.Map('KeyType', 'double', 'ValueType', 'char');
keys = partials.keys;
for i = 1:numel(keys)
    s = keys{i};
    if ~isKey(freeSymRows, s); continue; end
    p = partials(s);
    srow = freeSymRows(s);
    for j = 1:numel(srow)
        m = srow(j).slot;
        term = mulExpr(p, srow(j).expr);
        if isKey(acc, m)
            acc(m) = [acc(m) '+' term];
        else
            acc(m) = term;
        end
    end
end
ks = acc.keys;
row = struct('slot', {}, 'expr', {});
for i = 1:numel(ks)
    row(end+1) = struct('slot', ks{i}, 'expr', acc(ks{i})); %#ok<AGROW>
end
end

function s = mulExpr(a, b)
% Drop trivial *1 factors to keep the emitted C++ readable.
if strcmp(a, '1'); s = b; return; end
if strcmp(b, '1'); s = a; return; end
s = ['(' a ')*(' b ')'];
end
```

> All `sympymodel` functions (`seedUnknownRow`, `combine`, and later `generate`, `addReductionLoop`, `diffRHSloop`, `addInterpCall`, `emitTask`) live as one-file-per-function inside `+sympymodel/`. There is **no** top-level `emitModelSympy.m` (a function file and a same-named package cannot coexist).

- [ ] **Step 4: Run to verify pass**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); r=runtests('codegen/tEmitModelSympy'); exit(any([r.Failed]))"`
Expected: PASS (2/2).

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/+cxx/+sympymodel/ tests/codegen/tEmitModelSympy.m
git commit -m "feat(codegen): gradient-registry primitives (seedUnknownRow, combine)"
```

---

## Task 7: `sympymodel.generate` — flat model body → value + Jacobian C++

**Files:**
- Create: `src/+gdsge/+codegen/+cxx/+sympymodel/generate.m`
- Test: extend `tests/codegen/tEmitModelSympy.m`

`generate(ir)` returns a struct with `.argCode` (unknown unpack, reuse a double version of `emitArgUnpack`), `.popCode` (reuse `emitPop`), `.bodyCode` (helpers + value temporaries + JAC assembly), `.numEq`, `.numAux`. For this task, handle only: scalar unknown unpack, non-primed `assign`, non-primed `equation`. The registry seeds rows for each policy unknown; each assign asks Python for partials and `combine`s; each equation asks Python for partials, `combine`s, and emits `GDSGE_f[i]=value; if(GDSGE_jac){ JAC(i,m)=rowexpr; }`.

- [ ] **Step 1: Write the failing test — micro IR, snapshot the structure**

Add to `tEmitModelSympy.m` (methods Test), gated on sympy:

```matlab
        function flatModelEmitsValueAndJac(tc)
            tc.assumeTrue(gdsgetest.sympyAvailable());
            ir = microFlatIR();      % 2 unknowns x1,x2; eqs: x1*x2-3 ; x1+x2-4
            g = gdsge.codegen.cxx.sympymodel.generate(ir);
            tc.verifyEqual(g.numEq, 2);
            % value assignments present
            tc.verifySubstring(g.bodyCode, 'GDSGE_f[0]=');
            tc.verifySubstring(g.bodyCode, 'GDSGE_f[1]=');
            % Jacobian guarded and assigned
            tc.verifySubstring(g.bodyCode, 'if (GDSGE_jac)');
            tc.verifySubstring(g.bodyCode, 'JAC(0,');
            tc.verifySubstring(g.bodyCode, 'JAC(1,');
        end
```

And add the `microFlatIR` local function at the bottom of the test file:

```matlab
function ir = microFlatIR()
    n = @(v) gdsge.ir.node.num(v);
    nm = @(s) gdsge.ir.node.name(s);
    bo = @(op,a,b) gdsge.ir.node.binop(op,a,b);
    ir.irVersion = '1.0.0'; ir.modelName = 'micro'; ir.params = {};
    ir.options = struct('interpMethod','spline','interpOrder',4, ...
        'jacobianBackend','sympy');
    ir.shocks = struct('names',{{'z'}},'count',1,'values',struct('z',1), ...
        'transitions',struct('shock_trans',1));
    ir.states = struct('names',{{'k'}},'grids',struct('k','linspace(1,3,3)'));
    ir.variables = struct( ...
        'policy',{{struct('name','x1','length',1,'slot',[1 1]), ...
                   struct('name','x2','length',1,'slot',[2 2])}}, ...
        'aux',{{}}, 'interp',{{}}, 'tensor',{{}}, 'output',{{}}, 'others',{{}});
    ir.bounds = {struct('name','x1','lower','0','upper','10'), ...
                 struct('name','x2','lower','0','upper','10')};
    ir.interp = {};
    eq1 = struct('expr', bo('-', bo('*', nm('x1'), nm('x2')), n(3)), 'primed', 0);
    eq2 = struct('expr', bo('-', bo('+', nm('x1'), nm('x2')), n(4)), 'primed', 0);
    ir.model = struct('statements',{{}}, 'equations',{{eq1, eq2}});
    ir.simulate = struct('numPeriods',1,'numSamples',1,'initial',{{}}, ...
        'varSimu',{{}},'transitions',{{}});
    ir.hooks = struct();
end
```

- [ ] **Step 2: Run to verify fail**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('codegen/tEmitModelSympy'); exit(any([r.Failed]))"`
Expected: FAIL (or Incomplete if sympy unavailable) — `generate` does not exist.

- [ ] **Step 3: Implement `generate` (flat slice)**

Create `src/+gdsge/+codegen/+cxx/+sympymodel/generate.m`:

```matlab
function g = generate(ir)
% GENERATE  IR -> {argCode, popCode, bodyCode, numEq, numAux} for the SymPy
%   analytic-Jacobian model function (all-double; value + JAC). FLAT slice:
%   scalar-unknown unpack, non-primed assign, non-primed equation. Reductions
%   (addReductionLoop) and interp (addInterpCall) extend this in later tasks.
import gdsge.codegen.cxx.sympymodel.combine
import gdsge.codegen.cxx.sympymodel.seedUnknownRow
w = gdsge.codegen.codeWriter();

% ---- registry: name -> gradRow (struct array slot/expr) -----------------
rows = containers.Map('KeyType','char','ValueType','any');
% Seed unknown rows; also record the C++ symbol each unknown maps to.
% Scalar policy c -> C++ local `c`; array handled in a later task.
unknownSyms = {};
for i = 1:numel(ir.variables.policy)
    p = ir.variables.policy{i};
    assert(p.length == 1, 'array policy is handled in the interp/reduction tasks');
    s0 = p.slot(1) - 1;
    rows(p.name) = seedUnknownRow(s0, 1, 0);
    unknownSyms{end+1} = p.name; %#ok<AGROW>
end

% ---- non-primed assigns: register a new name with a combined row --------
for i = 1:numel(ir.model.statements)
    st = ir.model.statements{i};
    if ~strcmp(st.type, 'assign') || st.primed
        error('gdsge:codegen:sympyUnsupportedStmt', ...
            'flat slice: statement %s (primed=%d) not yet supported', st.type, st.primed);
    end
    d = diffRHS(st.expr, rows.keys);
    emitHelpers(w, d.helpers);
    w.add('double %s = %s;', st.target, d.value);
    rows(st.target) = combine(asMap(d.partials), rows);
end

% ---- equations: GDSGE_f[i] = value; JAC(i,m) = row ----------------------
numEq = 0;
jacLines = {};
for i = 1:numel(ir.model.equations)
    eq = ir.model.equations{i};
    assert(~eq.primed, 'primed equations handled in the reduction task');
    d = diffRHS(eq.expr, rows.keys);
    emitHelpers(w, d.helpers);
    w.add('GDSGE_f[%d]=%s;', numEq, d.value);
    erow = combine(asMap(d.partials), rows);
    for j = 1:numel(erow)
        jacLines{end+1} = sprintf('JAC(%d,%d)=%s;', numEq, erow(j).slot, erow(j).expr); %#ok<AGROW>
    end
    numEq = numEq + 1;
end

w.add('if (GDSGE_jac) {');
w.add('for (int GDSGE_ji=0; GDSGE_ji<%d*%d; ++GDSGE_ji) GDSGE_jac[GDSGE_ji]=0.0;', numEq, numEq);
for i = 1:numel(jacLines)
    w.add('%s', jacLines{i});
end
w.add('}');

g = struct();
g.numEq = numEq;
g.numAux = sum(cellfun(@(a) a.length, ir.variables.aux));
g.bodyCode = w.str();
g.argCode = argUnpackDouble(ir);
g.popCode = gdsge.codegen.cxx.emitPop(ir);
end

% ---- helpers -------------------------------------------------------------
function d = diffRHS(node, freeNames)
% Ask Python for value + partials w.r.t. every currently-registered name.
req = struct('body', node, 'diffVars', {reshape(freeNames,1,[])});
d = gdsge.codegen.sympy.callSympy('diff_body_json', req);
% jsondecode of an empty {} partials -> struct with no fields; normalize:
if ~isfield(d, 'partials'); d.partials = struct(); end
if ~isfield(d, 'helpers');  d.helpers = {}; end
end

function emitHelpers(w, helpers)
for i = 1:numel(helpers)
    h = helpers{i};               % {lhs, rhs}
    w.add('double %s = %s;', h{1}, h{2});
end
end

function m = asMap(partialsStruct)
% jsondecode turns {"a":"..","b":".."} into a struct; convert to a Map.
m = containers.Map('KeyType','char','ValueType','char');
fn = fieldnames(partialsStruct);
for i = 1:numel(fn)
    m(fn{i}) = partialsStruct.(fn{i});
end
end

function txt = argUnpackDouble(ir)
% double unknowns bound by slot (SymPy path has no adouble).
w = gdsge.codegen.codeWriter();
for i = 1:numel(ir.variables.policy)
    p = ir.variables.policy{i};
    if p.length == 1
        w.add('double& %s = GDSGE_x[%d];', p.name, p.slot(1) - 1);
    else
        w.add('double* %s_GDSGE_GRID = GDSGE_x+(%d);', p.name, p.slot(1) - 1);
        w.add('#define %s_GRID(idx) %s_GDSGE_GRID[int(idx)-1]', p.name, p.name);
        w.add('#define %s(idx) %s_GDSGE_GRID[int(idx)-1]', p.name, p.name);
    end
end
txt = w.str();
end
```

> **jsondecode caveat:** `callSympy` returns `helpers` as a cell of 2-cell arrays only if the JSON is a list of 2-lists — `jsondecode([["a","b"]])` yields a 1×2 cell or an N×2 char matrix depending on uniformity. Verify the actual decoded shape in Step 4; if `helpers` decodes as an N×2 cell matrix, change `emitHelpers` to index `helpers{i,1}`/`helpers{i,2}`. Pin this with the assertion in the next step.

- [ ] **Step 4: Run; fix the helpers/partials decode shape**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('codegen/tEmitModelSympy'); exit(any([r.Failed]))"`
Expected: PASS (3/3) when sympy is available. If `flatModelEmitsValueAndJac` errors inside `emitHelpers`/`asMap`, inspect the decoded `d` (add a temporary `disp(d)`), adjust indexing to match (`jsondecode` returns char matrices for uniform string lists), and re-run.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/+cxx/+sympymodel/generate.m tests/codegen/tEmitModelSympy.m
git commit -m "feat(codegen): sympymodel.generate — flat model value + analytic JAC"
```

---

## Task 8: Reductions — fused shock loop (EXPECT)

Extend `generate` to handle `reduction` statements. The body runs in a shock loop over per-shock symbols; primed shock realizations seed **empty** rows (constants), primed future-state arrays are handled in Task 9. For now support EXPECT with a body whose only unknown dependence is through current scalar unknowns (synthetic test); the interp chain rule arrives in Task 9.

**Files:**
- Modify: `src/+gdsge/+codegen/+cxx/+sympymodel/generate.m`
- Create: `src/+gdsge/+codegen/+cxx/+sympymodel/addReductionLoop.m`
- Test: extend `tEmitModelSympy.m`

- [ ] **Step 1: Write the failing test (synthetic EXPECT, FD micro-compile)**

This is the first **numeric** correctness check. Add a Slow test that builds a tiny model with one EXPECT reduction, emits the full model function into a standalone C++ harness, compiles it, and compares the analytic Jacobian to a finite-difference Jacobian. Add to `tEmitModelSympy.m`:

```matlab
    methods (Test, TestTags = {'Slow'})
        function expectJacobianMatchesFiniteDiff(tc)
            tc.assumeTrue(gdsgetest.sympyAvailable());
            % model: 1 unknown x; shock z in {1,2} with prob {0.5,0.5};
            % eq: GDSGE_EXPECT{ z'*x^2 } - 5 = 0  -> residual r(x)=E[z]*x^2-5
            % dr/dx = 2*E[z]*x. Compare emitted JAC vs central difference.
            ir = microExpectIR();
            g = gdsge.codegen.cxx.sympymodel.generate(ir);
            J = gdsgetest.microCompileJac(g, 1, [3.0]);   % NX=1, x=3
            r = gdsgetest.microCompileJacFD(g, 1, [3.0]);
            tc.verifyEqual(J, r, 'AbsTol', 1e-6, 'RelTol', 1e-6);
        end
    end
```

Add `microExpectIR` local at the bottom of the test file:

```matlab
function ir = microExpectIR()
    n = @(v) gdsge.ir.node.num(v);
    nm = @(s) gdsge.ir.node.name(s);
    pr = @(s) gdsge.ir.node.primed(s);
    bo = @(op,a,b) gdsge.ir.node.binop(op,a,b);
    ir.irVersion='1.0.0'; ir.modelName='microexp'; ir.params={};
    ir.options=struct('interpMethod','spline','interpOrder',4,'jacobianBackend','sympy');
    ir.shocks=struct('names',{{'z'}},'count',2,'values',struct('z',[1 2]), ...
        'transitions',struct('shock_trans',[0.5 0.5; 0.5 0.5]));
    ir.states=struct('names',{{'k'}},'grids',struct('k','linspace(1,3,3)'));
    ir.variables=struct('policy',{{struct('name','x','length',1,'slot',[1 1])}}, ...
        'aux',{{}},'interp',{{}},'tensor',{{}},'output',{{}},'others',{{}});
    ir.bounds={struct('name','x','lower','0','upper','10')};
    ir.interp={};
    red = struct('type','reduction','kind','EXPECT','target','es', ...
        'body', bo('*', pr('z'), bo('^', nm('x'), n(2))), 'transRef','shock_trans');
    eq = struct('expr', bo('-', nm('es'), n(5)), 'primed', 0);
    ir.model=struct('statements',{{red}}, 'equations',{{eq}});
    ir.simulate=struct('numPeriods',1,'numSamples',1,'initial',{{}}, ...
        'varSimu',{{}},'transitions',{{}});
    ir.hooks=struct();
end
```

- [ ] **Step 2: Write the micro-compile harness helpers**

Create `tests/+gdsgetest/microCompileJac.m` and `microCompileJacFD.m`. These wrap a tiny standalone (non-MEX) C++ harness so the registry's emitted body can be numerically checked without the full MEX pipeline:

`tests/+gdsgetest/microCompileJac.m`:

```matlab
function J = microCompileJac(g, nx, x)
% MICROCOMPILEJAC  Compile g.bodyCode into a standalone C++ exe that calls the
%   generated body with GDSGE_jac != 0 and prints the NX*NX Jacobian; return it
%   as an nx-by-nx matrix. Used by the synthetic reduction/interp gates.
%   Requires a working MSVC (same toolchain the MEX build uses).
J = gdsgetest.microRun(g, nx, x, true);
end
```

`tests/+gdsgetest/microCompileJacFD.m`:

```matlab
function J = microCompileJacFD(g, nx, x)
% Central-difference Jacobian from the value-only path of the same harness.
h = 1e-6;
J = zeros(nx, nx);
f0 = gdsgetest.microRun(g, nx, x, false);   % residual vector
for m = 1:nx
    xp = x; xp(m) = xp(m) + h;
    xm = x; xm(m) = xm(m) - h;
    fp = gdsgetest.microRun(g, nx, xp, false);
    fm = gdsgetest.microRun(g, nx, xm, false);
    J(:, m) = (fp - fm) / (2*h);
end
end
```

`tests/+gdsgetest/microRun.m` (the shared harness — wraps body in `main`, compiles with the vendored `include/`, runs, parses stdout):

```matlab
function out = microRun(g, nx, x, wantJac)
% MICRORUN  Emit a standalone C++ program that defines the generated model
%   function, calls it at x, and prints either the residual vector (wantJac=0)
%   or the flattened Jacobian (wantJac=1). Compiles with the repo include/.
here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(here));
incl = fullfile(repoRoot, 'include');
work = tempname; mkdir(work);
src = fullfile(work, 'micro.cpp');
% Minimal preamble: pow/exp/log come from <cmath>; the body references
% GDSGE_x/GDSGE_f/GDSGE_jac and the JAC macro + the shock-loop locals the
% generated reduction code emits (shock_trans(...), z_GRID(...), shock).
pre = sprintf([ ...
    '#include <cstdio>\n#include <cmath>\n#include <cstring>\nusing namespace std;\n' ...
    '#define JAC(i_eq,i_var) GDSGE_jac[(i_eq) + %d*(i_var)]\n' ...
    'int main(){\n' ...
    '  double GDSGE_x[%d]={%s};\n' ...
    '  double GDSGE_f[%d]={0};\n' ...
    '  double GDSGE_jac[%d]={0};\n' ...
    '  double* pjac = %s;\n' ...
    '  int GDSGE_NUM_SHOCKS=%d; int shock=1;\n' ...
    '%s\n' ...                      % data setup (shock_trans, z grid) from g.dataSetup
    '%s\n' ...                      % the body, with GDSGE_jac replaced by pjac
    '  for(int i=0;i<%d;++i) printf("%%.17g\\n", GDSGE_f[i]);\n' ...
    '  if(pjac) for(int i=0;i<%d;++i) printf("J %%.17g\\n", GDSGE_jac[i]);\n' ...
    '  return 0;\n}\n'], ...
    nx, nx, strjoin(compose('%.17g', x(:)'), ','), nx, nx*nx, ...
    ternary(wantJac, 'GDSGE_jac', '0'), 0, g.dataSetup, ...
    strrep(g.bodyCode, 'GDSGE_jac', 'pjac'), nx, nx*nx);
gdsge.codegen.writeText(src, pre);
% Compile with MSVC via mex's compiler or system cl. Reuse the project's
% compiler discovery: mex -setup C++ has been run. Use `mex` to compile a
% standalone is awkward; instead use the same cl.exe the MEX build uses.
exe = fullfile(work, 'micro.exe');
cmd = sprintf('cl /nologo /EHsc /I"%s" "%s" /Fe:"%s" /Fo:"%s\\\\"', incl, src, exe, work);
[st, msg] = system(cmd);
assert(st == 0, 'micro compile failed:\n%s', msg);
[st, so] = system(['"' exe '"']);
assert(st == 0, 'micro run failed:\n%s', so);
out = parseMicroStdout(so, nx, wantJac);
end
```

> **Implementation note:** the shock-loop reduction body the registry emits references `shock_trans((shock)+N*(GDSGE_iter-1))`, `z_GRID(GDSGE_iter)`, and `GDSGE_NUM_SHOCKS`. `generate` must therefore also return `g.dataSetup` — C++ that declares `double GDSGE_shock_trans[...]`, the `#define shock_trans(idx)`, the `z_GDSGE_GRID[...]` + `#define z_GRID(idx)` for the synthetic harness (in the real MEX these come from `emitPop`). Add `g.dataSetup` to `generate`'s return, built from the IR shock values/transition for the micro harness only (guard behind a flag so the real path keeps using `emitPop`). Define `parseMicroStdout`, `ternary` as local helpers. If `cl.exe` is not directly on PATH, discover it via `mex.getCompilerConfigurations('C++')` and prepend its `Location`.

- [ ] **Step 3: Run to verify fail**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('codegen/tEmitModelSympy'); exit(any([r.Failed]))"`
Expected: FAIL — `generate` errors on the `reduction` statement (flat slice rejects it).

- [ ] **Step 4: Implement `addReductionLoop` and wire it into `generate`**

Create `src/+gdsge/+codegen/+cxx/+sympymodel/addReductionLoop.m`:

```matlab
function rows = addReductionLoop(w, st, rows, shockCount)
% ADDREDUCTIONLOOP  Emit a fused shock loop for an EXPECT/MIN/MAX/PROD reduction
%   accumulating value (st.target) and gradient row (d<target>[NX]). Registers
%   st.target's row as the sparse set of slots that the loop's dacc touches.
%   Per-shock primed inputs: shock realizations seed empty rows (constants);
%   future-state interp results are added by addInterpCall (Task 9) BEFORE the
%   reduction that consumes them, so their rows are already in `rows`.
import gdsge.codegen.cxx.sympymodel.combine
n = shockCount;
% Differentiate the body w.r.t. all registered names + primed shock names.
% Primed shock realizations (e.g. z') are constants -> not in rows -> drop out.
d = gdsge.codegen.cxx.sympymodel.diffRHSloop(st.body, rows.keys);
% Determine which slots the gradient touches (union over the body's free syms).
brow = combine(d.partialsMap, rows);
touched = unique([brow.slot]);

tgt = st.target;
% init accumulators
switch st.kind
    case 'EXPECT', w.add('double %s = 0.0;', tgt);
    case 'PROD',   w.add('double %s = 1.0;', tgt);
    case 'MIN',    w.add('double %s = 1e20;', tgt);
    case 'MAX',    w.add('double %s = -1e20;', tgt);
end
if ~isempty(touched)
    w.add('double d%s[%d];', tgt, maxSlot(touched)+1);  % sized to max touched+1
    w.add('for(int GDSGE_jz=0; GDSGE_jz<%d; ++GDSGE_jz) d%s[GDSGE_jz]=0.0;', ...
        maxSlot(touched)+1, tgt);
end
w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++) {', n);
% emit helpers + body value/partials inside the loop (per-shock)
for i = 1:numel(d.helpers)
    w.add('double %s = %s;', d.helpers{i}{1}, d.helpers{i}{2});
end
w.add('double GDSGE_body = %s;', d.value);
switch st.kind
    case 'EXPECT'
        w.add('%s += %s((shock)+%d*(GDSGE_iter-1))*GDSGE_body;', tgt, st.transRef, n);
        for s = touched
            w.add('d%s[%d] += %s((shock)+%d*(GDSGE_iter-1))*(%s);', ...
                tgt, s, st.transRef, n, brow(slotIdx(brow,s)).expr);
        end
    case 'PROD'
        for s = touched
            w.add('d%s[%d] = d%s[%d]*GDSGE_body + %s*(%s);', ...
                tgt, s, tgt, s, tgt, brow(slotIdx(brow,s)).expr);
        end
        w.add('%s *= GDSGE_body;', tgt);
    case {'MIN','MAX'}
        cmp = '<'; if strcmp(st.kind,'MAX'); cmp = '>'; end
        w.add('if (GDSGE_body %s %s) {', cmp, tgt);
        w.add('%s = GDSGE_body;', tgt);
        for s = touched
            w.add('d%s[%d] = %s;', tgt, s, brow(slotIdx(brow,s)).expr);
        end
        w.add('}');
end
w.add('}');
% register target's row pointing at the dacc array entries
trow = struct('slot',{},'expr',{});
for s = touched
    trow(end+1) = struct('slot', s, 'expr', sprintf('d%s[%d]', tgt, s)); %#ok<AGROW>
end
rows(tgt) = trow;   %#ok<NASGU>  % containers.Map handle: mutated in place
end

function k = maxSlot(touched), k = max(touched); end
function i = slotIdx(row, s), i = find([row.slot]==s, 1); end
```

Add `diffRHSloop` (a variant of `diffRHS` that also surfaces a `partialsMap` containers.Map) — create `src/+gdsge/+codegen/+cxx/+sympymodel/diffRHSloop.m`:

```matlab
function d = diffRHSloop(node, freeNames)
% DIFFRHSLOOP  Like generate.diffRHS but returns helpers as a cell of {lhs,rhs}
%   pairs and partials as a containers.Map, for the reduction emitter.
req = struct('body', node, 'diffVars', {reshape(freeNames,1,[])});
raw = gdsge.codegen.sympy.callSympy('diff_body_json', req);
d.value = raw.value;
d.helpers = normalizeHelpers(raw);
d.partialsMap = containers.Map('KeyType','char','ValueType','char');
if isfield(raw,'partials')
    fn = fieldnames(raw.partials);
    for i = 1:numel(fn); d.partialsMap(fn{i}) = raw.partials.(fn{i}); end
end
end

function H = normalizeHelpers(raw)
% jsondecode shape-normalize: helpers may decode as {} , Nx2 cell, or 1x2 cell.
H = {};
if ~isfield(raw,'helpers') || isempty(raw.helpers); return; end
h = raw.helpers;
if iscell(h)
    if size(h,2) == 2 && ~iscell(h{1})    % Nx2 cell of chars
        for i = 1:size(h,1); H{end+1} = {h{i,1}, h{i,2}}; end %#ok<AGROW>
    else                                   % cell of 1x2 cells
        for i = 1:numel(h); H{end+1} = {h{i}{1}, h{i}{2}}; end %#ok<AGROW>
    end
end
end
```

In `generate.m`, replace the assign-only statement loop with a dispatch that handles `reduction` via `addReductionLoop` and keeps `assign`/`interpCall` (interp in Task 9). Refactor `diffRHS` to also expose `partials` as a map (reuse `diffRHSloop`). Update the equation loop to allow `primed` equations later (Task 9); for now keep the non-primed assert.

- [ ] **Step 5: Run to verify pass**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('codegen/tEmitModelSympy'); exit(any([r.Failed]))"`
Expected: PASS — `expectJacobianMatchesFiniteDiff` confirms analytic == FD to 1e-6.

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+codegen/+cxx/+sympymodel/ tests/codegen/tEmitModelSympy.m tests/+gdsgetest/microCompileJac.m tests/+gdsgetest/microCompileJacFD.m tests/+gdsgetest/microRun.m
git commit -m "feat(codegen): EXPECT reduction fused loop (value+grad), FD-verified"
```

---

## Task 9: Interp chain rule + primed equations (the HL1996 shape)

Add `interpCall` handling: per-shock interp results `<target>_j = interp(<args>_j)` whose gradient rows chain through `interp'` (the kernel `grad`). Add primed equations (expand across the shock loop), and array-policy unknown unpack (`w1n[8]`). After this task `generate` handles the full HL1996 model body.

**Files:**
- Modify: `src/+gdsge/+codegen/+cxx/+sympymodel/generate.m`
- Create: `src/+gdsge/+codegen/+cxx/+sympymodel/addInterpCall.m`
- Test: extend `tEmitModelSympy.m` with a synthetic interp FD check

The interp result's gradient row (spec §5.2): for `psn_j = interp_psn(w1n_j)` with `w1n_j` the j-th element of array unknown `w1n` (slot `s0..`), `row = {slot(w1n_j): GDSGE_INTERP_GRAD[interpIdx + numInterp*0]}` (single-state ⇒ dim 0). For multi-arg interp, sum a term per arg with its own grad column and the arg's own row.

- [ ] **Step 1: Write the failing synthetic interp test**

Add a Slow test `interpChainRuleMatchesFiniteDiff` building a 1-state model with: array unknown `xn[2]`, an interp `f_future` over the state, an interp call `fn' = GDSGE_INTERP_VEC'(xn')` producing `fn_1,fn_2`, an `EXPECT{ z' * fn' }` reduction, and an equation. The micro harness must provide the interp evaluator + a known cubic spline so `interp'` is analytic and FD-checkable. (Use a linear interp via `interp_lite.h` with a known slope to keep the harness simple, since the chain rule only needs *some* differentiable evaluator with a known derivative buffer.)

```matlab
        function interpChainRuleMatchesFiniteDiff(tc)
            tc.assumeTrue(gdsgetest.sympyAvailable());
            ir = microInterpIR();
            g = gdsge.codegen.cxx.sympymodel.generate(ir);
            x = [1.5; 2.5];
            J  = gdsgetest.microCompileJac(g, 2, x);
            Jf = gdsgetest.microCompileJacFD(g, 2, x);
            tc.verifyEqual(J, Jf, 'AbsTol', 1e-5, 'RelTol', 1e-5);
        end
```

Provide `microInterpIR` (1 state, array policy `xn[2]`, interp `f_future`, interpCall primed `fn`, EXPECT, equation) and extend `microRun` to inject a stub interp evaluator that returns `value = a*site + b` with `grad = a` (a known-slope linear map), so the analytic chain rule and the FD agree. Document the stub clearly in `microRun`.

- [ ] **Step 2: Run to verify fail**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('codegen/tEmitModelSympy'); exit(any([r.Failed]))"`
Expected: FAIL — `interpCall` rejected by `generate`.

- [ ] **Step 3: Implement `addInterpCall` + array-unknown seeding + primed equations**

Create `src/+gdsge/+codegen/+cxx/+sympymodel/addInterpCall.m`:

```matlab
function rows = addInterpCall(w, st, ir, rows, numInterp)
% ADDINTERPCALL  Emit the per-shock interp evaluation (value + grad) and seed
%   each target's per-shock gradient row from interp' and the eval-point row.
%   For GDSGE_INTERP_VEC' the targets map to interp indices 0..K-1; the single
%   arg is a primed array unknown (e.g. xn'), whose element j has its own seed
%   row. Single-state (xdim=1) here; multi-arg sums one term per dim.
%
%   The loop is folded into the reduction loop in the real model (Task 8's
%   GDSGE_iter). addInterpCall runs BEFORE the reduction that consumes its
%   targets, registering rows expressed in terms of the loop index — so it must
%   emit its own evaluation inside the SAME shock loop. Implementation: emit the
%   interp evaluation at the top of the reduction loop body via a deferred
%   snippet returned to generate, OR (simpler) require interpCall + the reduction
%   that uses it to share one loop. Here we emit interpCall results into per-
%   shock arrays <tgt>_GDSGE_GRID[j] and their grad arrays <tgt>_DGRID[j], then
%   the reduction body references <tgt>_GDSGE_GRID[GDSGE_iter-1].
args = st.args;                  % cell of arg nodes; single primed array unknown
assert(numel(args) == 1, 'sympy interp: multi-arg interp is Phase 8 follow-up');
% the arg is primed name of an array unknown; find its slot base
argNode = args{1};
assert(strcmp(argNode.kind,'primed'), 'interp arg must be primed');
argName = argNode.id;
[s0, len] = arrayUnknownSlot(ir, argName);
n = ir.shocks.count;
% evaluate value + grad per shock into arrays
w.add('double %s_VAL[%d]; double %s_DER[%d];', ...
    argName, n, argName, n);   % scratch for the single interp arg's results
for t = 1:numel(st.targets)
    w.add('double %s_GDSGE_GRID[%d]; double %s_DER[%d];', ...
        st.targets{t}, n, st.targets{t}, n);
end
w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++) {', n);
w.add('double GDSGE_site = %s_GDSGE_GRID[GDSGE_iter-1];', argName);  % xn'[j]
w.add('double GDSGE_g[%d];', numInterp);
w.add('GDSGE_INTERP_VEC_double_grad(GDSGE_iter, GDSGE_site, GDSGE_INTERP_RSLT_double, GDSGE_g);');
for t = 1:numel(st.targets)
    w.add('%s_GDSGE_GRID[GDSGE_iter-1]=GDSGE_INTERP_RSLT_double[%d];', st.targets{t}, t-1);
    w.add('%s_DER[GDSGE_iter-1]=GDSGE_g[%d];', st.targets{t}, t-1);
end
w.add('}');
% register each target's per-shock row: at slot (s0 + (j-1)) with expr <tgt>_DER[j-1].
% The reduction loop consumes these by indexing [GDSGE_iter-1]; the registry
% stores the row keyed by the SYMBOL the body sees (the primed target name),
% and combine multiplies by the loop-shock body partial. To express "the slot
% is the j-th element of the array unknown", store a TEMPLATED row whose slot
% and expr depend on GDSGE_iter; addReductionLoop substitutes the loop var.
for t = 1:numel(st.targets)
    rows(st.targets{t}) = struct( ...
        'slot', sprintf('%d+(GDSGE_iter-1)', s0), ...
        'expr', sprintf('%s_DER[GDSGE_iter-1]', st.targets{t}), ...
        'templated', true);
end
end

function [s0, len] = arrayUnknownSlot(ir, name)
s0 = []; len = [];
for i = 1:numel(ir.variables.policy)
    p = ir.variables.policy{i};
    if strcmp(p.name, name); s0 = p.slot(1)-1; len = p.length; return; end
end
error('gdsge:codegen:sympyInterpArg', 'interp arg %s is not an array policy', name);
end
```

> **Design clarification (templated rows):** an interp target's eval point is a *per-shock* unknown (`w1n[j]`), so its gradient row's **slot depends on the loop index** `j`. Represent such rows with `templated=true` and slot/expr strings containing `GDSGE_iter`. In `addReductionLoop`, when combining a body that references templated rows, emit the `dacc[...]` update **inside** the shock loop with the slot expression `s0+(GDSGE_iter-1)` used as the array index: `d<tgt>[s0+(GDSGE_iter-1)] += trans*(partial)*(<tgt>_DER[GDSGE_iter-1])`. Extend `combine`/`addReductionLoop` to handle string slots (templated) by emitting the accumulation inside the loop rather than precomputing `touched`. Non-templated (scalar-unknown) slots keep the integer path. This is the crux of §10; keep the two paths explicit.

Also: seed **array-policy** unknowns in `generate` (currently asserts length==1). For an array unknown, do **not** seed a scalar row; its per-element rows are produced by `addInterpCall`/`seedUnknownRow(s0,len,k)` as needed. Add primed-equation handling mirroring `emitEquations.m` (loop `GDSGE_EQ[base-1+GDSGE_iter]`), with the Jacobian row emitted inside the same loop using templated slots.

- [ ] **Step 4: Run to verify pass**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('codegen/tEmitModelSympy'); exit(any([r.Failed]))"`
Expected: PASS — `interpChainRuleMatchesFiniteDiff` confirms the chain rule.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+codegen/+cxx/+sympymodel/ tests/codegen/tEmitModelSympy.m tests/+gdsgetest/
git commit -m "feat(codegen): interp chain rule (templated per-shock grad rows) + primed eqs"
```

---

## Task 10: Double with-grad interp evaluator + `model_sympy.tpl.cpp` + `generateCxx` branch

Wire `sympymodel.generate` output into a full MEX `.cpp`: a new interp section exposing `GDSGE_INTERP_VEC_double_grad` (calls `search_eval_with_grad_vec_at_array`), the `model_sympy.tpl.cpp` skeleton, and a `generateCxx` branch selected by `jacobianBackend='sympy'`.

**Files:**
- Create: `src/+gdsge/+codegen/+cxx/emitInterpSympy.m`
- Create: `templates/cxx/model_sympy.tpl.cpp`
- Modify: `src/+gdsge/+codegen/generateCxx.m`
- Test: `tests/HeatonLucas1996/codegen/tSympyJacobianHL1996.m` (structure only, fast)

- [ ] **Step 1: Write `emitInterpSympy`**

Create `src/+gdsge/+codegen/+cxx/emitInterpSympy.m`:

```matlab
function frag = emitInterpSympy(ir)
% EMITINTERPSYMPY  Spline interp section for the SymPy backend: a double
%   evaluator that returns BOTH value and gradient (search_eval_with_grad_*).
%   frag.getCode (task scope), frag.threadCode (per-thread scratch + lambda).
import gdsge.codegen.cxx.fillTemplate
import gdsge.codegen.cxx.readTemplate
numInterp = numel(ir.interp);
if numInterp == 0
    frag = struct('getCode','', 'threadCode',''); return;
end
args1 = ir.interp{1}.args; xdim = numel(args1);
assert(all(cellfun(@(it) numel(it.args)==xdim, ir.interp)), ...
    'sympy: all interp vars must share arg dimension');
order = ir.options.interpOrder;
getCode = fillTemplate(readTemplate('interp_spline_construct.tpl.cpp'), ...
    {'VAR_NUM', num2str(xdim); 'NUM_INTERP', num2str(numInterp)});
% per-thread scratch + the with-grad vector evaluator
w = gdsge.codegen.codeWriter();
w.add('int GDSGE_INTERP_CELL[%d] = {0};', xdim);
w.add('double GDSGE_INTERP_RSLT_double[%d] = {0};', numInterp);
w.add('double GDSGE_INTERP_GRAD[%d] = {0};', numInterp*xdim);
% GDSGE_INTERP_VEC_double_grad(shockIdx, site..., out, grad)
sites = strjoin(cellfun(@(v) ['double ' v], args1, 'UniformOutput', false), ',');
varNames = strjoin(args1, ',');
w.add('auto GDSGE_INTERP_VEC_double_grad = [&](int shockIdx, %s, double* GDSGE_out, double* GDSGE_grad){', sites);
w.add('  double xSite[] = {%s};', varNames);
w.add('  GDSGE_CSPLINE_VEC.search_eval_with_grad_vec_at_array(shockIdx-1, xSite, GDSGE_out, GDSGE_grad, GDSGE_INTERP_CELL);');
w.add('};');
frag = struct('getCode', getCode, 'threadCode', w.str());
end
```

> **Reconcile the signature** with `addInterpCall`'s call `GDSGE_INTERP_VEC_double_grad(GDSGE_iter, GDSGE_site, GDSGE_INTERP_RSLT_double, GDSGE_g)`. For HL1996 `xdim=1`, the lambda is `(int shockIdx, double w1, double* out, double* grad)`. Ensure `addInterpCall` passes the single site `GDSGE_site` (the primed eval point) — matches. `GDSGE_g` receives the per-interp gradient; for `xdim=1` it is length `numInterp`. (Multi-state interp is a documented follow-up.)

- [ ] **Step 2: Write the template**

Create `templates/cxx/model_sympy.tpl.cpp`:

```cpp
        // SymPy analytic-Jacobian model function (all double; value + JAC).
        #define JAC(i_eq,i_var) GDSGE_jac[(i_eq) + NUM_EQUATIONS*(i_var)]
        auto GDSGE_FUNC_MODEL_1_double = [&](double* GDSGE_x, double* GDSGE_f, double* GDSGE_jac)
        {
            ARG_CODE;

            MODEL_BODY_CODE;

            AUX_ASSIGN_CODE;
        };
        #undef JAC
        auto GDSGE_OBJ_MODEL_1 = [&](double* GDSGE_x, double* GDSGE_f, double* GDSGE_jac)
        {
            GDSGE_FUNC_MODEL_1_double(GDSGE_x, GDSGE_f, GDSGE_jac);
        };
```

> The `task.tpl.cpp` `CALL_FMIN_CODE` block calls `GDSGE_OBJ_MODEL_NUMBER` and ends by re-evaluating with `GDSGE_FUNC_MODEL_NUMBER_double(&GDSGE_x[0], &GDSGE_EQ[0], 1)`. The SymPy model function ignores a trailing eval-mode arg, but the call_fmin template passes `,1`. Provide a 4th defaulted param `int GDSGE_EVAL=0` on `GDSGE_FUNC_MODEL_1_double` to keep that call valid — OR emit a SymPy-specific call_fmin in Task 11. Choose the defaulted-param route: change the lambda signature to `(double* GDSGE_x, double* GDSGE_f, double* GDSGE_jac=0, int GDSGE_EVAL=0)` so both `GDSGE_OBJ_MODEL_1(x,f,jac)` and the final `..._double(x,EQ,1)` (jac defaults to 0 → value only) compile. Confirm CoDoSol's functor call passes exactly `(x,f,jac)`.

- [ ] **Step 3: Write `emitModelSympyTask` (assembles task.tpl.cpp for the sympy path)**

Create `src/+gdsge/+codegen/+cxx/+sympymodel/emitTask.m` (parallels `emitTask.m`, but uses the sympy model + interp sections, and the existing `call_fmin.tpl.cpp`):

```matlab
function txt = emitTask(ir, taskName)
% EMITTASK (sympy)  Assemble task.tpl.cpp for the analytic-Jacobian backend.
if nargin < 2; taskName = 'task_inf_horizon'; end
import gdsge.codegen.cxx.fillTemplate
import gdsge.codegen.cxx.readTemplate
g = gdsge.codegen.cxx.sympymodel.generate(ir);
interp = gdsge.codegen.cxx.emitInterpSympy(ir);
model = fillTemplate(readTemplate('model_sympy.tpl.cpp'), { ...
    'ARG_CODE',        g.argCode; ...
    'MODEL_BODY_CODE', g.bodyCode; ...
    'AUX_ASSIGN_CODE', auxAssignDouble(ir)});
model = fillTemplate(model, {'NUM_EQUATIONS', num2str(g.numEq)});
model = strrep(model, 'MODEL_NUMBER', 'MODEL_1');
callFMin = fillTemplate(readTemplate('call_fmin.tpl.cpp'), { ...
    'MODEL_CONDITION', '1'; 'NUM_EQUATIONS', num2str(g.numEq)});
callFMin = strrep(callFMin, 'MODEL_NUMBER', 'MODEL_1');
txt = fillTemplate(readTemplate('task.tpl.cpp'), { ...
    'PRE_MODEL_CODE',        ir.hooks.preModel; ...
    'START_LOOP_CODE',       ''; 'FINISH_LOOP_CODE', ''; ...
    'MODEL_CODE',            [model newline callFMin_decl_note()]; ...
    'CALL_FMIN_CODE',        callFMin; ...
    'POP_CODE',              g.popCode; ...
    'INTERP_GET_CODE',       interp.getCode; ...
    'INTERP_GET_THREAD_CODE', interp.threadCode; ...
    'NUM_EQUATIONS',         num2str(g.numEq); ...
    'NUM_AUX',               num2str(g.numAux); ...
    'TASK_NAME',             taskName});
end

function s = callFMin_decl_note(), s = ''; end

function txt = auxAssignDouble(ir)
% aux outputs, plain double (no adept value()).
w = gdsge.codegen.codeWriter();
for i = 1:numel(ir.variables.aux)
    a = ir.variables.aux{i};
    if a.length == 1
        w.add('GDSGE_aux[%d]=%s;', a.slot(1)-1, a.name);
    else
        w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++)', a.length);
        w.add('{ GDSGE_aux[%d+GDSGE_iter]=%s_GDSGE_GRID[GDSGE_iter-1]; }', ...
            a.slot(1)-2, a.name);
    end
end
txt = w.str();
end
```

> Aux that are themselves reduction/interp results are registered as scalar locals by `generate`; for HL1996 the aux (`c1,c2,ps,pb`...) are policy or simple assigns — verify each aux name is in scope at `AUX_ASSIGN_CODE`. If an aux is a per-shock array, it is already a `_GDSGE_GRID` local.

- [ ] **Step 4: Branch `generateCxx`**

In `src/+gdsge/+codegen/generateCxx.m`, after `gdsge.codegen.assertSupportedIR(ir);` and the existing pchip/asg/cxx/tensor guards, add the sympy branch. Replace the single `cpp = fillTemplate(...TASK_INF_HORIZON_CODE..., emitTask(ir))` with:

```matlab
useSympy = isfield(ir.options,'jacobianBackend') && strcmp(ir.options.jacobianBackend,'sympy');
if useSympy
    if ~strcmp(ir.options.interpMethod, 'spline')
        error('gdsge:codegen:sympyInterpUnsupported', ...
            'sympy backend supports only spline interpolation (got %s)', ir.options.interpMethod);
    end
    if ~isempty(ir.hooks.cxx)
        error('gdsge:codegen:cxxUnderSympy', ...
            'sympy backend cannot differentiate cxx blocks');
    end
    if ~gdsge.codegen.sympy.ensurePyenv()
        error('gdsge:codegen:sympyPythonUnavailable', ...
            ['sympy backend requires the uv Python env. Run: ' ...
             'uv sync --project pyext']);
    end
    taskFn = @gdsge.codegen.cxx.sympymodel.emitTask;
else
    taskFn = @gdsge.codegen.cxx.emitTask;
end
```

Then change the task emit lines to use `taskFn`:

```matlab
taskInitCode = '';
if isfield(ir, 'modelInit')
    taskInitCode = taskFn(gdsge.codegen.initView(ir), 'task_init');
end
cpp = fillTemplate(readTemplate('mex.tpl.cpp'), { ...
    'GDSGE_OTHER_INCLUDE',    cxxIncludeText(ir); ...
    'TASK_INIT_CODE',         taskInitCode; ...
    'TASK_INF_HORIZON_CODE',  taskFn(ir)});
```

> HL1996 has no `modelInit`, so `task_init` is not exercised; the sympy `emitTask` is only called for the main task. Keep the `taskFn` indirection so init works once a sympy model with init appears (follow-up).

- [ ] **Step 5: Write a fast structure test**

Create `tests/HeatonLucas1996/codegen/tSympyJacobianHL1996.m` with a fast (non-compiling) method first:

```matlab
classdef tSympyJacobianHL1996 < matlab.unittest.TestCase
    % PHASE 8: HL1996 with jacobianBackend='sympy' generates a .cpp whose model
    % function computes value + analytic JAC (fast structural check), and
    % (Slow) the analytic Jacobian agrees with adept and finite differences.
    properties
        Work
    end
    methods (TestClassSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', '..', 'src')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', 'ir')));               % buildHL1996IR
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', '+gdsgetest')));
            tc.assumeTrue(gdsgetest.sympyAvailable());
        end
    end
    methods (Test)
        function generatesAnalyticModel(tc)
            ir = buildHL1996IR();
            ir.options.jacobianBackend = 'sympy';
            tc.Work = tempname; mkdir(tc.Work);
            gdsge.codegen.generateCxx(ir, tc.Work);
            cpp = fileread(fullfile(tc.Work, 'mex_HL1996.cpp'));
            tc.verifySubstring(cpp, 'GDSGE_FUNC_MODEL_1_double');
            tc.verifySubstring(cpp, 'JAC(');
            tc.verifySubstring(cpp, 'search_eval_with_grad_vec_at_array');
            tc.verifyEmpty(regexp(cpp, '\badouble\b', 'once'), ...
                'sympy model must be adouble-free');
        end
    end
end
```

- [ ] **Step 6: Run to verify pass**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('HeatonLucas1996/codegen/tSympyJacobianHL1996'); exit(any([r.Failed]))"`
Expected: PASS (1/1) — the `.cpp` generates with the analytic structure, no adouble. Iterate on `generate`/emitters until HL1996 emits without error.

- [ ] **Step 7: Commit**

```bash
git add src/+gdsge/+codegen/+cxx/emitInterpSympy.m templates/cxx/model_sympy.tpl.cpp src/+gdsge/+codegen/+cxx/+sympymodel/emitTask.m src/+gdsge/+codegen/generateCxx.m tests/HeatonLucas1996/codegen/tSympyJacobianHL1996.m
git commit -m "feat(codegen): generateCxx sympy branch + model_sympy template + double interp eval"
```

---

## Task 11: Compile defines for the SymPy MEX + 6th-output Jacobian hook

The SymPy MEX needs no adept Jacobian but still includes adept headers (the templates reference `value()` etc. only on the autodiff path — verify the sympy `.cpp` compiles as-is). The cross-check needs the analytic Jacobian out of the MEX: add an **additive 6th output** filled in eval mode.

**Files:**
- Modify: `templates/cxx/task.tpl.cpp` (additive `plhs[5]` jac, `GDSGE_DEBUG_EVAL_ONLY==2`)
- Modify: `templates/cxx/call_fmin.tpl.cpp` (fill jac when `GDSGE_DEBUG_EVAL_ONLY==2`)
- Modify: `src/+gdsge/+codegen/+cxx/emitCompile.m` (no change unless adept omitted — see note)
- Test: covered by Task 12's compile.

- [ ] **Step 1: Add the optional Jacobian output to `task.tpl.cpp`**

In `templates/cxx/task.tpl.cpp`, after the `plhs[4]` block (line ~42), add (additive — only when the caller requests ≥6 outputs):

```cpp
  // Optional analytic-Jacobian debug output (Phase 8 cross-check): per problem,
  // NUM_EQUATIONS*NUM_EQUATIONS column-major. Filled only when nlhs>=6 and
  // GDSGE_DEBUG_EVAL_ONLY==2 (the eval path writes it). Existing 5-output
  // callers are unaffected.
  double* GDSGE_JAC_OUT = 0;
  if (nlhs >= 6) {
    plhs[5] = mxCreateDoubleMatrix(NUM_EQUATIONS*NUM_EQUATIONS, GDSGE_NPROB, mxREAL);
    GDSGE_JAC_OUT = mxGetPr(plhs[5]);
  }
```

> This block is shared by both backends (it is in the common task template). The adept path can also fill it (it already computes a Jacobian via the functor), giving the adept witness for the cross-check. Verify the adept `.cpp` still compiles (the new lines are plain C++). Because the autodiff `.cpp` is regenerated from the same template, **regenerate and re-run the existing HL1996 autodiff gates** to confirm byte-compatible behavior (the cache key changes since the template changed — expect one recompile, then identical results).

- [ ] **Step 2: Fill the Jacobian in `call_fmin.tpl.cpp`**

In `templates/cxx/call_fmin.tpl.cpp`, change the final eval block to also write the Jacobian when requested:

```cpp
  // Evaluate at the solution
  double* GDSGE_x = &GDSGE_SOL(1, GDSGE_I);
  #if MAXDIM>MAX_STACK_DIM
  vector<double> GDSGE_EQ(NUM_EQUATIONS);
  #else
  double GDSGE_EQ[NUM_EQUATIONS];
  #endif

  if (GDSGE_DEBUG_EVAL_ONLY==2 && GDSGE_JAC_OUT) {
    GDSGE_OBJ_MODEL_NUMBER(&GDSGE_x[0], &GDSGE_EQ[0], &GDSGE_JAC_OUT[(GDSGE_I-1)*NUM_EQUATIONS*NUM_EQUATIONS]);
  } else {
    GDSGE_FUNC_MODEL_NUMBER_double(&GDSGE_x[0], &GDSGE_EQ[0], 1);
  }
}
```

> `GDSGE_OBJ_MODEL_NUMBER` is the solver functor `(x,f,jac)` for both backends (adept: records; sympy: analytic). When `jac != 0` both fill the Jacobian. This is the single line that exposes both witnesses through one code path. Confirm the adept `GDSGE_OBJ_MODEL_NUMBER` accepts a non-null jac in eval context (it does — it is the same functor CoDoSol calls).

- [ ] **Step 3: Decide adept inclusion for the sympy MEX**

The sympy `.cpp` includes `adept_extension.h` via `mex.tpl.cpp` and declares `Stack _stack;` in `task.tpl.cpp`. That is harmless (the sympy model never records), and keeps `mex.tpl.cpp` shared/untouched. **Do not** strip adept now (YAGNI; parity bar is correctness, not compile-time). Leave `emitCompile.m` unchanged for spline sympy models. Note this decision in the commit message.

- [ ] **Step 4: Commit**

```bash
git add templates/cxx/task.tpl.cpp templates/cxx/call_fmin.tpl.cpp
git commit -m "feat(codegen): additive 6th-output analytic Jacobian (GDSGE_DEBUG_EVAL_ONLY==2)"
```

---

## Task 12: HL1996 three-way Jacobian cross-check (the decisive gate)

**Files:**
- Modify: `tests/HeatonLucas1996/codegen/tSympyJacobianHL1996.m` (add the Slow cross-check)
- Reuse: `tests/+gdsgetest/buildHL1996MexInputs.m`, `tests/HeatonLucas1996/oldmex` or the autodiff MEX.

- [ ] **Step 1: Write the failing Slow cross-check**

Add to `tSympyJacobianHL1996.m`:

```matlab
    methods (Test, TestTags = {'Slow'})
        function analyticJacMatchesAdeptAndFD(tc)
            % Build BOTH MEX (adept + sympy) for HL1996, evaluate the Jacobian
            % at the converged golden solution via the 6th output (eval mode 2),
            % and assert sympy == adept == central-difference, per grid point.
            ir = buildHL1996IR();
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            gdsge.runtime.ensurePath();

            % --- adept MEX ---
            adir = fullfile(work, 'ad'); mkdir(adir);
            gdsge.codegen.generateCxx(ir, adir);
            buildMex(tc, adir, 'HL1996');
            copyfile(fullfile(adir, ['mex_HL1996.' mexext]), ...
                fullfile(work, ['mex_HL1996_ad.' mexext]));

            % --- sympy MEX ---
            sir = ir; sir.options.jacobianBackend = 'sympy';
            sdir = fullfile(work, 'sy'); mkdir(sdir);
            gdsge.codegen.generateCxx(sir, sdir);
            buildMex(tc, sdir, 'HL1996');
            copyfile(fullfile(sdir, ['mex_HL1996.' mexext]), ...
                fullfile(work, ['mex_HL1996_sy.' mexext]));

            oldCd = pwd; cd(work); tc.addTeardown(@() cd(oldCd));
            S = gdsgetest.buildHL1996MexInputs();
            % evaluate the analytic Jacobian (6 outputs, eval mode 2)
            Jad = evalJac(@mex_HL1996_ad, S);
            Jsy = evalJac(@mex_HL1996_sy, S);
            Jfd = fdJac(@mex_HL1996_ad, S);          % central diff via residuals

            % compare on a sample of grid points (all 1608 is fine; sample for speed)
            cols = 1:50:size(Jad, 2);
            r = gdsgetest.compareNumericClose(Jsy(:,cols), Jad(:,cols), 1e-6, 1e-8);
            tc.verifyTrue(r.pass, ['sympy vs adept: ' strjoin(r.failures, newline)]);
            r = gdsgetest.compareNumericClose(Jsy(:,cols), Jfd(:,cols), 1e-4, 1e-6);
            tc.verifyTrue(r.pass, ['sympy vs FD: ' strjoin(r.failures, newline)]);
        end
    end
```

Add local functions in the test file:

```matlab
function buildMex(tc, dir, model)
    oldCd = pwd; cd(dir); c = onCleanup(@() cd(oldCd));
    gdsge.runtime.ensurePath();
    feval(['compile_' model]);
    tc.assertTrue(exist(fullfile(dir, ['mex_' model '.' mexext]), 'file') == 3, ...
        'MEX did not compile');
end

function J = evalJac(mexFn, S)
    % 6th output: NUM_EQ^2 x NPROB, eval mode 2 (analytic Jacobian at sol0).
    cfg = S.cfg; cfg.debugEvalOnly = 2;
    [~,~,~,~,~,jac] = callMex6(mexFn, S, cfg);
    J = jac;
end

function J = fdJac(mexFn, S)
    % Per-problem central-difference Jacobian from residuals (eval mode 1).
    nx = size(S.sol0, 1); np = size(S.sol0, 2); h = 1e-6;
    J = zeros(nx*nx, np);
    cfg = S.cfg; cfg.debugEvalOnly = 1;
    for m = 1:nx
        Sp = S; Sp.sol0 = S.sol0; Sp.sol0(m,:) = Sp.sol0(m,:) + h;
        Sm = S; Sm.sol0 = S.sol0; Sm.sol0(m,:) = Sm.sol0(m,:) - h;
        [~,~,~,ep] = gdsge.runtime.solveProblems(mexFn, Sp.sol0, S.lb, S.ub, S.data, S.f0, S.aux0, S.eqval0, cfg);
        [~,~,~,em] = gdsge.runtime.solveProblems(mexFn, Sm.sol0, S.lb, S.ub, S.data, S.f0, S.aux0, S.eqval0, cfg);
        d = (ep - em) / (2*h);            % nx x np residual derivative wrt x_m
        for ieq = 1:nx
            J((ieq-1) + nx*(m-1) + 1, :) = d(ieq, :);
        end
    end
end

function varargout = callMex6(mexFn, S, cfg)
    % Mirror solveProblems' caller-workspace contract for a single 6-output
    % eval call (no resolve loop). The MEX reads these locals via mexGetVariable.
    TolFun=cfg.tolFun; TolSol=cfg.tolSol; SolMaxIter=cfg.solMaxIter; %#ok<NASGU>
    NumThreads=cfg.numThreads; GDSGE_DEBUG_EVAL_ONLY=cfg.debugEvalOnly; %#ok<NASGU>
    UseBroyden=cfg.useBroyden; FiniteDiffDelta=cfg.finiteDiffDelta; %#ok<NASGU>
    GDSGE_USE_BROYDEN_NOW=cfg.useBroydenNow; MEX_TASK_NAME=cfg.taskName; %#ok<NASGU>
    MEX_TASK_INIT=0; MEX_TASK_INF_HORIZON=1; GDSGE_SPLINE_VEC=cfg.splineVec; %#ok<NASGU>
    for iPP=1:numel(cfg.ppNames); assert(isvarname(cfg.ppNames{iPP})); %#ok<*NASGU>
        eval([cfg.ppNames{iPP} ' = cfg.ppCell{iPP};']); end
    skip = zeros(1, size(S.sol0,2));
    [varargout{1:6}] = mexFn(S.sol0, S.lb, S.ub, S.data, skip, S.f0, S.aux0, S.eqval0);
end
```

> `callMex6` duplicates the MEX caller-workspace contract from `solveProblems.m` because the 6-output eval call must stay in one function body (the MEX reads via `mexGetVariable('caller',...)`). Keep it local to the test. If a future refactor adds a `solveProblems` jac-return mode, this can be replaced.

- [ ] **Step 2: Run to verify it fails first (then drives the fix loop)**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('HeatonLucas1996/codegen/tSympyJacobianHL1996'); exit(any([r.Failed]))"`
Expected: initially FAIL (compile errors in the sympy `.cpp`, or Jacobian mismatch). This is the integration gate that drives debugging of Tasks 6–11. Fix `generate`/emitters until `sympy == adept == FD`.

- [ ] **Step 3: Iterate to green**

Debug systematically (use superpowers:systematic-debugging): if compile fails, read the sympy `.cpp` model function; if a specific JAC column is wrong, the registry chain rule for that unknown is wrong — check templated rows (interp-of-`w1n[j]`) and the EXPECT loop accumulation. The FD witness localizes which equation/unknown.

- [ ] **Step 4: Commit**

```bash
git add tests/HeatonLucas1996/codegen/tSympyJacobianHL1996.m
git commit -m "test(codegen): HL1996 three-way Jacobian cross-check (sympy=adept=FD)"
```

---

## Task 13: HL1996 end-to-end golden (sympy backend)

**Files:**
- Create: `tests/HeatonLucas1996/codegen/tEndToEndHL1996Sympy.m`

- [ ] **Step 1: Write the gate**

Model it on `tEndToEndHL1996.m`, but force the sympy backend and reuse the SAME golden (both backends solve the same equations to the same `TolEq`):

```matlab
classdef tEndToEndHL1996Sympy < matlab.unittest.TestCase
    % PHASE 8 GATE: HL1996 with the SymPy analytic-Jacobian backend converges to
    % the committed golden (Iter/Metric/policy/aux/interp + simulate shock path).
    properties (Constant)
        RelTol = 1e-4; AbsTol = 1e-4;
    end
    methods (TestClassSetup)
        function gate(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', '..', 'src')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', '+gdsgetest')));
            tc.assumeTrue(gdsgetest.sympyAvailable());
        end
    end
    methods (Test, TestTags = {'Slow'})
        function sympyBackendMatchesGolden(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            % Write a UseAutoDiff=0 variant of the gmod so the public surface
            % selects the sympy backend.
            gmod = fileread(fullfile(modelDir, 'HL1996.gmod'));
            gmod = ['UseAutoDiff=0;' newline gmod];
            fid = fopen(fullfile(work, 'HL1996.gmod'), 'w'); fwrite(fid, gmod); fclose(fid);

            ir = gdsge_codegen('HL1996');
            tc.assertEqual(ir.options.jacobianBackend, 'sympy');
            tc.assertTrue(exist(fullfile(work, ['mex_HL1996.' mexext]), 'file') == 3);

            opts = struct('SaveFreq', inf, 'NoSave', 1);
            IterRslt = iter_HL1996(opts);
            golden = load(fullfile(modelDir, 'golden', 'IterRslt.mat'));
            G = golden.IterRslt;
            tc.verifyLessThan(IterRslt.Metric, 1e-6);
            tc.verifyGreaterThan(IterRslt.Iter, 100);
            tc.verifyLessThan(IterRslt.Iter, 400);
            r = gdsgetest.compareNumericClose(IterRslt.var_policy, G.var_policy, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_aux, G.var_aux, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));
            r = gdsgetest.compareNumericClose(IterRslt.var_interp, G.var_interp, tc.RelTol, tc.AbsTol);
            tc.verifyTrue(r.pass, strjoin(r.failures, newline));

            rng(0823);
            SimuRslt = simulate_HL1996(IterRslt);
            gs = load(fullfile(modelDir, 'golden', 'SimuRslt.mat')); GS = gs.SimuRslt;
            tc.verifyEqual(SimuRslt.shock, GS.shock);
            flds = {'w1','c1','c2','ps','pb','equity_premium'};
            for i = 1:numel(flds)
                a = SimuRslt.(flds{i}); b = GS.(flds{i});
                tc.verifyLessThan(abs(mean(a(:)) - mean(b(:))), 5e-3, flds{i});
                tc.verifyLessThan(abs(std(a(:)) - std(b(:))), 5e-3, flds{i});
            end
        end
    end
end
```

> `iter_HL1996`/`simulate_HL1996` are generated identically for both backends (MATLAB backend untouched); only `mex_HL1996` differs. The `UseAutoDiff=0;` line is a plain MATLAB assignment in the gmod setup — `resolveOptions` reads it from the eval'd workspace.

- [ ] **Step 2: Run to verify pass**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('HeatonLucas1996/codegen/tEndToEndHL1996Sympy'); exit(any([r.Failed]))"`
Expected: PASS — the sympy MEX drives policy iteration to the golden (Iter≈209, Metric<1e-6).

- [ ] **Step 3: Commit**

```bash
git add tests/HeatonLucas1996/codegen/tEndToEndHL1996Sympy.m
git commit -m "test(codegen): HL1996 end-to-end gate on the sympy backend"
```

---

## Task 14: Synthetic MIN/MAX/PROD unit tests

The corpus gates EXPECT (all 4) + MIN (safe_assets). PROD/MAX share the lowering; pin them with synthetic FD checks so the lowering is exercised even though no gated model uses them.

**Files:**
- Modify: `tests/codegen/tEmitModelSympy.m`

- [ ] **Step 1: Add MIN, MAX, PROD synthetic FD checks**

For each kind, build a micro IR like `microExpectIR` but with the reduction kind swapped, and assert `microCompileJac == microCompileJacFD`. For MIN/MAX, choose the body and `x` so the argmin/argmax is unique at the test point (away from the kink), e.g. body `= z' * x` with distinct `z'` per shock so the extremum is at a single shock. Add three Slow methods mirroring `expectJacobianMatchesFiniteDiff`.

```matlab
        function minJacobianMatchesFiniteDiff(tc)
            tc.assumeTrue(gdsgetest.sympyAvailable());
            ir = microReductionIR('MIN');     % body z'*x, x>0 -> argmin = min z'
            J  = gdsgetest.microCompileJac(ir2g(ir), 1, [3.0]);
            Jf = gdsgetest.microCompileJacFD(ir2g(ir), 1, [3.0]);
            tc.verifyEqual(J, Jf, 'AbsTol', 1e-6, 'RelTol', 1e-6);
        end
        % maxJacobianMatchesFiniteDiff, prodJacobianMatchesFiniteDiff: identical
        % shape with 'MAX' / 'PROD'.
```

Add `microReductionIR(kind)` (parameterize `microExpectIR`'s reduction `kind`) and `ir2g = @(ir) gdsge.codegen.cxx.sympymodel.generate(ir)` as a local. Pick shock values `z=[1 2]` (PROD), `z=[1 2]` with `x>0` (MIN→shock 1, MAX→shock 2) so the extremum branch is stable under the ±h perturbation.

- [ ] **Step 2: Run to verify pass**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('codegen/tEmitModelSympy'); exit(any([r.Failed]))"`
Expected: PASS — MIN/MAX/PROD analytic Jacobians match FD.

- [ ] **Step 3: Commit**

```bash
git add tests/codegen/tEmitModelSympy.m
git commit -m "test(codegen): synthetic MIN/MAX/PROD reduction Jacobian FD checks"
```

---

## Task 15: Widen — safe_assets, Mendoza2010, GLSW end-to-end (sympy)

**Files:**
- Create: `tests/Barro_et_al_2017/codegen/tEndToEndSafeAssetsSympy.m`
- Create: `tests/Mendoza2010/codegen/tEndToEndMendoza2010Sympy.m`
- Create: `tests/GLSW2020/codegen/tEndToEndGLSWSympy.m`

Each mirrors its existing autodiff `tEndToEnd*` gate (copy structure, fields, tolerances, seeds), prepending `UseAutoDiff=0;` to the gmod and asserting `ir.options.jacobianBackend=='sympy'`. safe_assets exercises the **MIN** reduction; Mendoza is the **large-NX** stress (verify the `vector<double>` fallback path); GLSW exercises its EXPECT shape.

- [ ] **Step 1: safe_assets gate**

Copy `tests/Barro_et_al_2017/codegen/tEndToEndSafeAssets.m` → `...Sympy.m`, rename the class, add the `assumeTrue(gdsgetest.sympyAvailable())` setup, prepend `UseAutoDiff=0;` to the copied gmod, assert the backend, and keep the same golden comparisons. Read the original first to copy its exact field list/tolerances (do not invent them).

- [ ] **Step 2: Run safe_assets**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('Barro_et_al_2017/codegen/tEndToEndSafeAssetsSympy'); exit(any([r.Failed]))"`
Expected: PASS. If the MIN reduction's gradient causes slow/failed convergence, revisit `addReductionLoop` MIN branch (extremal-`j` subgradient) — the cross-check style (add a `tSympyJacobianSafeAssets` three-way check) localizes it.

- [ ] **Step 3: Commit safe_assets**

```bash
git add tests/Barro_et_al_2017/codegen/tEndToEndSafeAssetsSympy.m
git commit -m "test(codegen): safe_assets end-to-end on sympy backend (MIN reduction)"
```

- [ ] **Step 4: Mendoza2010 gate**

Copy `tests/Mendoza2010/codegen/tEndToEndMendoza2010.m` → `...Sympy.m` likewise. Mendoza is the large-NX case: confirm the generated sympy `.cpp` uses the `#if MAXDIM>MAX_STACK_DIM` fallback for `GDSGE_EQ` and that the `d<name>[NX]` registry arrays fit (if a model's NX is large enough to stress the stack, switch those arrays to the `vector<double>` fallback in `generate` under the same `MAXDIM>MAX_STACK_DIM` guard — add this only if the build/run shows stack issues; log it).

- [ ] **Step 5: Run + commit Mendoza**

Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('Mendoza2010/codegen/tEndToEndMendoza2010Sympy'); exit(any([r.Failed]))"`
Expected: PASS.

```bash
git add tests/Mendoza2010/codegen/tEndToEndMendoza2010Sympy.m src/+gdsge/+codegen/+cxx/+sympymodel/
git commit -m "test(codegen): Mendoza2010 end-to-end on sympy backend (large-NX)"
```

- [ ] **Step 6: GLSW gate + run + commit**

Copy `tests/GLSW2020/codegen/tEndToEndGLSW.m` → `...Sympy.m` likewise.
Run: `matlab -batch "cd('tests'); addpath(genpath('../src')); addpath('+gdsgetest'); r=runtests('GLSW2020/codegen/tEndToEndGLSWSympy'); exit(any([r.Failed]))"`
Expected: PASS.

```bash
git add tests/GLSW2020/codegen/tEndToEndGLSWSympy.m
git commit -m "test(codegen): GLSW end-to-end on sympy backend"
```

---

## Task 16: Full suite, autodiff non-regression, and PROGRESS update

**Files:**
- Modify: `PROGRESS.md`

- [ ] **Step 1: Run the full MATLAB suite**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: exit 0. Inspect `tests/results/junit.xml`. Sympy gates show as **Incomplete (skipped)** if the venv is unsynced, **Passed** if synced; nothing **Failed**.

- [ ] **Step 2: Confirm autodiff non-regression**

The `task.tpl.cpp`/`call_fmin.tpl.cpp` changes (Task 11) touched the shared templates. Confirm the existing autodiff gates still match goldens: `tEndToEndHL1996`, `tEndToEndSafeAssets`, `tEndToEndMendoza2010`, `tEndToEndGLSW`, `tMexEquivalenceHL1996`, and the ASG gates (`tEndToEndCaoKS2016`, `tEndToEndBianchi2011`) all pass. These are included in the full-suite run; verify none failed in `junit.xml`.

- [ ] **Step 3: Run the Python suite**

Run: `uv run --project pyext pytest pyext/tests -q`
Expected: PASS (all).

- [ ] **Step 4: Update PROGRESS.md**

Mark Phase 8 done. Edit the phase line:

```markdown
- ☑ **Phase 8 — SymPy analytic-Jacobian backend** (done 2026-06-13) — Approach C:
  MATLAB owns the fused shock-loop structure + a gradient registry (sparse slot→C++
  rows, forward-mode chain rule); SymPy (pyenv bridge, JSON in/out) differentiates
  each body and returns value+partials as shared-CSE C++. Chain rule closes through
  the spline `search_eval_with_grad_vec_at_array` (double). Selected by
  `UseAutoDiff=0` → `options.jacobianBackend='sympy'` (emitted only when non-default;
  zero golden regen). Gated three-way (sympy=adept=finite-diff Jacobian) + end-to-end
  golden on HL1996/safe_assets/Mendoza2010/GLSW; MIN/MAX/PROD synthetic FD checks.
  cxx/asg/pchip/var_tensor under sympy → clear errors. ASG → Phase 8b (deferred).
```

Add a changelog entry at the top of the Changelog section summarizing the work, branch `phase8-sympy-jacobian-backend`, suite count, and "Next: Phase 8b (ASG) or Phase 9 (polish)".

- [ ] **Step 5: Commit**

```bash
git add PROGRESS.md
git commit -m "docs(progress): Phase 8 complete (SymPy analytic-Jacobian backend)"
```

---

## Self-review checklist (run before handing off for execution)

- **Spec coverage:** §4 selection+guards → Tasks 1,2,10; §5.1 bridge → Task 3; §5.2/5.3 registry+contract → Tasks 4–9; §6 emitted C++ (template, loop shapes, storage, interp eval) → Tasks 8–11; §7 testing (Python unit, MATLAB unit, three-way cross-check, end-to-end, availability gating) → Tasks 4–9,12,13,15,16; §8 deferrals (ASG/pchip/cxx/tensor, PROD/MAX synthetic) → Tasks 1,10,14; §9 risks (stack fallback, pyenv gating, MIN kink) → Tasks 3,9,14,15. ✓
- **Non-regression of autodiff (byte-identical):** guaranteed by emitting `jacobianBackend` only when non-default (Task 2) and by leaving the autodiff emitters/templates untouched except the additive 6th-output (Task 11, re-verified Task 16). ✓
- **MIN/MAX kink:** Task 14 places test points away from the kink (documented). ✓
- **Known follow-ups (documented, not gated):** multi-state interp args in `emitInterpSympy`/`addInterpCall` (single-state today); ASG (Phase 8b); `task_init` for a sympy model with `model_init` (the `taskFn` indirection is in place but untested — no spline corpus model with init uses sympy). These are noted inline, not silent.
```
