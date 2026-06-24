# Phase 8b (ASG slice) — SymPy backend for ASG — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make ASG models run on the SymPy analytic-Jacobian backend (`UseAutoDiff=0`), matching their autodiff goldens, so the ASG path reaches backend parity with the spline path.

**Architecture:** The ASG interpolant already computes value+gradient in pure `double` (`asg_adouble.h::eval_vec_adept`). Add a pure-double class method `eval_vec_with_grad` (mirrors `eval_vec_adouble` minus the adept injection), emit a `GDSGE_INTERP_VEC_double_grad` lambda backed by it (same name/signature the spline SymPy path emits, so the whole `+sympymodel` registry/chain-rule layer is untouched), and lift the `generateCxx` guard from "spline only" to "spline or asg". Validate with a SymPy↔adept↔FD Jacobian cross-check and three ASG end-to-end gates (CaoKS2016, Bianchi2011_asg, and a new corpus model CaoNie2016_asg that also exercises ASG × Phase-9b conditional regions).

**Tech Stack:** MATLAB R2025b (`matlab.unittest`), C++ MEX (adept autodiff + the vendored `asg_mex`/`asg_adouble.h`), SymPy via the uv pyenv. Headless run: `matlab -batch`. Spec: `docs/superpowers/specs/2026-06-15-phase8b-sympy-asg-design.md`. Branch: `phase8b-sympy-asg`.

**Conventions used below**
- **Full suite:** `matlab -batch "cd('tests'); run_tests"` (exit 0 = all pass; report `tests/results/junit.xml`). Never run two MATLAB processes at once (each saturates all cores).
- **Single test (fast iteration)** — define once, reused in steps as *RUN(<relpath>)*:
  ```
  matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile('..','src')); addpath(fullfile('..','src','kernels')); r=runtests('<relpath>'); disp(table([r.Passed]',[r.Failed]','VariableNames',{'Passed','Failed'},'RowNames',{r.Name})); assert(all([r.Passed]),'TEST FAILED');"
  ```
  `<relpath>` is relative to `tests/`, no `.m` (e.g. `codegen/tEmitInterpSympyAsg`).
- SymPy tests must guard with `tc.assumeTrue(gdsgetest.sympyAvailable())` so the default Python-free loop stays green.

---

## Task 1: SymPy-ASG coverage probe (investigation — go/no-go, no code)

Confirm the three ASG models use only SymPy-supported constructs **before** writing code. The SymPy backend was validated on the four spline models; this checks the ASG corpus introduces no new construct.

**Files:** read-only — `tests/CaoKS2016/ir/CaoKS2016.gdsge.json`, `tests/Bianchi2011_asg/ir/bianchi2011.gdsge.json`, `tests/CaoNie2016/CaoNie2016.gmod` (the cartesian twin of the new model).

- [ ] **Step 1: Decode and inspect the two existing ASG IRs**

Run:
```
matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile('..','src')); for f={'CaoKS2016/ir/CaoKS2016.gdsge.json','Bianchi2011_asg/ir/bianchi2011.gdsge.json'}, ir=gdsge.ir.decode(fileread(f{1})); fprintf('--- %s: interpMethod=%s, #regions=%d, #interp=%d\n', f{1}, ir.options.interpMethod, numel(ir.model.regions), numel(ir.interp)); end"
```
Expected: both `interpMethod=asg`; `#regions=1` (no conditional regions); `#interp>=1`.

- [ ] **Step 2: Confirm reduction kinds and equation kinds are SymPy-supported**

The SymPy backend handles reductions `EXPECT/MIN/MAX/PROD`, single/multi-state `GDSGE_INTERP_VEC'`, named scalar interp, primed assigns/equations, constant-indexed per-shock access, and (since Phase 9b) conditional regions + `if/else` equations. Inspect each region's `equations` kinds and any reduction nodes; confirm none is outside this set. CaoNie2016_asg's body is **identical** to the cartesian `tests/CaoNie2016/CaoNie2016.gmod`, which already passes SymPy (`tEndToEndCaoNie2016Sympy`), so it needs no separate construct check — only its ASG+region *combination* is new (Tasks 7–10).

- [ ] **Step 3: Record findings**

If everything is in the supported set, write one line in the PR/commit notes ("coverage probe: CaoKS2016/Bianchi2011 use only EXPECT/INTERP_VEC/named-interp; no conditional regions; CaoNie2016_asg body == cartesian CaoNie2016 (already SymPy-green)") and proceed. **If any construct is unsupported, stop and revise the spec** — do not code around it silently.

---

## Task 2: `eval_vec_with_grad` (C++) + the SymPy ASG interp emitter

Add the pure-double ASG value+grad method and make `emitInterpSympy` emit an ASG-backed `GDSGE_INTERP_VEC_double_grad` lambda. The emitter is unit-tested here; the C++ method's correctness is gated by compilation in Task 4 and the cross-check in Task 6.

**Files:**
- Modify: `include/asg_adouble.h` (add a method to `AsgInterpArrayAdoubleEvaluator`, after `eval_vec_adouble`, ~line 623)
- Modify: `src/+gdsge/+codegen/+cxx/emitInterpSympy.m`
- Test: `tests/codegen/tEmitInterpSympyAsg.m` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/codegen/tEmitInterpSympyAsg.m`:
```matlab
classdef tEmitInterpSympyAsg < matlab.unittest.TestCase
    % Phase 8b: the SymPy interp emitter, given an ASG IR, emits the ASG
    % construct + a GDSGE_INTERP_VEC_double_grad lambda backed by
    % eval_vec_with_grad, with the native->spline gradient transpose.
    methods (TestClassSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src')));
        end
    end
    methods (Test)
        function asgEmitsDoubleGradLambda(tc)
            here = fileparts(mfilename('fullpath'));
            ir = gdsge.ir.decode(fileread(fullfile(here, '..', ...
                'CaoKS2016', 'ir', 'CaoKS2016.gdsge.json')));
            ir.options.jacobianBackend = 'sympy';
            frag = gdsge.codegen.cxx.emitInterpSympy(ir);
            % getCode: the ASG construct (class-handle evaluator)
            tc.verifySubstring(frag.getCode, 'convertMat2Ptr<AsgInterpArray>');
            tc.verifySubstring(frag.getCode, 'GDSGE_CPP_ASG');
            % threadCode: the shared lambda name + the ASG double-grad call + transpose
            tc.verifySubstring(frag.threadCode, 'GDSGE_INTERP_VEC_double_grad');
            tc.verifySubstring(frag.threadCode, 'eval_vec_with_grad');
            tc.verifySubstring(frag.threadCode, 'GDSGE_ASG_GRAD_NATIVE');
            % must NOT fall through to the spline kernel call
            tc.verifyEmpty(strfind(frag.threadCode, 'search_eval_with_grad_vec_at_array'), ...
                'ASG path must not emit the spline kernel call');
        end
        function splineStillEmitsKernelCall(tc)
            % regression: spline path unchanged
            here = fileparts(mfilename('fullpath'));
            ir = gdsge.ir.decode(fileread(fullfile(here, '..', ...
                'HeatonLucas1996', 'ir', 'HL1996.gdsge.json')));
            ir.options.jacobianBackend = 'sympy';
            frag = gdsge.codegen.cxx.emitInterpSympy(ir);
            tc.verifySubstring(frag.threadCode, 'search_eval_with_grad_vec_at_array');
            tc.verifyEmpty(strfind(frag.threadCode, 'eval_vec_with_grad'));
        end
    end
end
```

- [ ] **Step 2: Run the test, verify it fails**

Run *RUN(`codegen/tEmitInterpSympyAsg`)*.
Expected: FAIL — the current `emitInterpSympy` ignores `interpMethod` and emits the spline kernel call for the ASG IR (so `asgEmitsDoubleGradLambda` fails on the `convertMat2Ptr` / `eval_vec_with_grad` substrings).

- [ ] **Step 3: Add the C++ method to `include/asg_adouble.h`**

Insert into `class AsgInterpArrayAdoubleEvaluator`, immediately after the closing `}` of `eval_vec_adouble` (just before the class-closing `};` near line 623):
```cpp
        // Pure-double value + gradient (w.r.t. the real site). Mirrors
        // eval_vec_adouble's evalGradFlag==1 branch but stays in double and
        // does NOT inject into an adept tape — for the SymPy analytic backend.
        // gradient native layout: gradient[i_vec*numDim + i_dim].
        void eval_vec_with_grad(int i_array, double* site, double* cell, double* ratio,
            double* slope, double* rslt, double* gradient)
        {
            double siteScaled[ASG_MAX_DIM];
            for (int i_dim = 0; i_dim < numDim; i_dim++)
                siteScaled[i_dim] = (site[i_dim] - stateMin[i_dim]) / stateRange[i_dim];
            AsgInterp& interp = interps[i_array];
            search_adept_with_slope(siteScaled, cell, ratio, slope, currentLevel, numDim);
            eval_vec_adept(numVec, cell, ratio, slope, currentLevel, numDim,
                interp.info, interp.setOfLevelCombinations, rslt, gradient);
            for (int i_vec = 0; i_vec < numVec; i_vec++)
                for (int i_dim = 0; i_dim < numDim; i_dim++)
                    gradient[i_vec*numDim + i_dim] /= stateRange[i_dim];
        }
```

- [ ] **Step 4: Add the ASG branch to `src/+gdsge/+codegen/+cxx/emitInterpSympy.m`**

After the `xdim`/assert block (after current line 17) and before the `getCode = fillTemplate(...interp_spline_construct...)` line, factor the site param strings up and branch. Replace the body from line 19 onward so the function reads:
```matlab
siteParams = strjoin(arrayfun(@(d) sprintf('double GDSGE_site%d', d-1), 1:xdim, ...
    'UniformOutput', false), ', ');
siteList = strjoin(arrayfun(@(d) sprintf('GDSGE_site%d', d-1), 1:xdim, ...
    'UniformOutput', false), ',');

if strcmp(ir.options.interpMethod, 'asg')
    frag = asgFrag(xdim, numInterp, siteParams, siteList);
    return;
end

getCode = fillTemplate(readTemplate('interp_spline_construct.tpl.cpp'), ...
    {'VAR_NUM', num2str(xdim); 'NUM_INTERP', num2str(numInterp)});
w = gdsge.codegen.codeWriter();
w.add('int GDSGE_INTERP_CELL[%d] = {0};', xdim);
w.add('auto GDSGE_INTERP_VEC_double_grad = [&GDSGE_CSPLINE_VEC,&GDSGE_INTERP_CELL](int shockIdx, %s, double* GDSGE_out, double* GDSGE_grad){', siteParams);
w.add('  double xSite[] = {%s};', siteList);
w.add('  GDSGE_CSPLINE_VEC.search_eval_with_grad_vec_at_array(shockIdx-1, xSite, GDSGE_out, GDSGE_grad, GDSGE_INTERP_CELL);');
w.add('};');
frag = struct('getCode', getCode, 'threadCode', w.str());
end

function frag = asgFrag(xdim, numInterp, siteParams, siteList)
% SymPy ASG interp: reuse the autodiff ASG construct (declares GDSGE_CPP_ASG),
% emit a double-only per-thread scratch block + the shared-name double-grad
% lambda backed by eval_vec_with_grad, transposing native [v*DIM+d] -> the
% spline layout [v + NVEC*d] that addInterpCall consumes. ASG_MAX_LEVEL stays
% a compile-time define; DIM/NVEC are substituted numerically.
import gdsge.codegen.cxx.readTemplate
getCode = readTemplate('interp_asg_construct.tpl.cpp');
w = gdsge.codegen.codeWriter();
w.add('double GDSGE_ASG_CELL[(ASG_MAX_LEVEL+2)*%d] = {0};', xdim);
w.add('double GDSGE_ASG_RATIO[(ASG_MAX_LEVEL+2)*%d] = {0};', xdim);
w.add('double GDSGE_ASG_SLOPE[(ASG_MAX_LEVEL+2)*%d] = {0};', xdim);
w.add('double GDSGE_ASG_GRAD_NATIVE[%d*%d] = {0};', numInterp, xdim);
w.add('auto GDSGE_INTERP_VEC_double_grad = [&GDSGE_CPP_ASG,&GDSGE_ASG_CELL,&GDSGE_ASG_RATIO,&GDSGE_ASG_SLOPE,&GDSGE_ASG_GRAD_NATIVE](int shockIdx, %s, double* GDSGE_out, double* GDSGE_grad){', siteParams);
w.add('  double xSite[] = {%s};', siteList);
w.add('  GDSGE_CPP_ASG.eval_vec_with_grad(shockIdx-1, xSite, GDSGE_ASG_CELL, GDSGE_ASG_RATIO, GDSGE_ASG_SLOPE, GDSGE_out, GDSGE_ASG_GRAD_NATIVE);');
w.add('  for(int GDSGE_v=0; GDSGE_v<%d; GDSGE_v++) for(int GDSGE_d=0; GDSGE_d<%d; GDSGE_d++)', numInterp, xdim);
w.add('    GDSGE_grad[GDSGE_v + %d*GDSGE_d] = GDSGE_ASG_GRAD_NATIVE[GDSGE_v*%d + GDSGE_d];', numInterp, xdim);
w.add('};');
frag = struct('getCode', getCode, 'threadCode', w.str());
end
```
(Keep the existing function header/comment and the `numInterp==0 → empty frag` early return at the top unchanged.)

- [ ] **Step 5: Run the test, verify it passes**

Run *RUN(`codegen/tEmitInterpSympyAsg`)*.
Expected: PASS (both methods).

- [ ] **Step 6: Commit**

```
git add include/asg_adouble.h src/+gdsge/+codegen/+cxx/emitInterpSympy.m tests/codegen/tEmitInterpSympyAsg.m
git commit -m "feat(sympy-asg): eval_vec_with_grad + ASG branch in emitInterpSympy"
```

---

## Task 3: Lift the `generateCxx` SymPy guard + structural .cpp gate

Allow `interpMethod=='asg'` under the SymPy backend (keep pchip + cxx-hook rejections). Verify the generated `.cpp` is structurally correct without a full compile (fast).

**Files:**
- Modify: `src/+gdsge/+codegen/generateCxx.m:32-36`
- Test: `tests/codegen/tGenerateCxxAsgSympy.m` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/codegen/tGenerateCxxAsgSympy.m`:
```matlab
classdef tGenerateCxxAsgSympy < matlab.unittest.TestCase
    % Phase 8b: generateCxx accepts asg under jacobianBackend='sympy' and emits
    % an analytic-Jacobian ASG model .cpp (no MEX build — content checks only).
    methods (TestClassSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', 'src', 'kernels')));  % asg.get_mex_constants
        end
    end
    methods (Test)
        function asgSympyGenerates(tc)
            here = fileparts(mfilename('fullpath'));
            ir = gdsge.ir.decode(fileread(fullfile(here, '..', ...
                'CaoKS2016', 'ir', 'CaoKS2016.gdsge.json')));
            ir.options.jacobianBackend = 'sympy';
            tc.assumeTrue(gdsgetest.sympyAvailable());
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            gdsge.codegen.generateCxx(ir, work);
            cpp = fileread(fullfile(work, 'mex_CaoKS2016.cpp'));
            tc.verifySubstring(cpp, 'GDSGE_INTERP_VEC_double_grad');
            tc.verifySubstring(cpp, 'eval_vec_with_grad');
            tc.verifySubstring(cpp, 'JAC(');
            tc.verifyEmpty(strfind(cpp, 'GDSGE_FUNC_MODEL_1_adouble'), ...
                'sympy model must not emit the adouble model lambda');
        end
        function pchipStillRejected(tc)
            here = fileparts(mfilename('fullpath'));
            ir = gdsge.ir.decode(fileread(fullfile(here, '..', ...
                'CaoKS2016', 'ir', 'CaoKS2016.gdsge.json')));
            ir.options.interpMethod = 'pchip';
            ir.options.jacobianBackend = 'sympy';
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            tc.verifyError(@() gdsge.codegen.generateCxx(ir, work), ...
                'gdsge:codegen:unsupported');   % pchip is rejected before the sympy guard
        end
    end
end
```

- [ ] **Step 2: Run the test, verify it fails**

Run *RUN(`codegen/tGenerateCxxAsgSympy`)*.
Expected: `asgSympyGenerates` FAILS with `gdsge:codegen:sympyInterpUnsupported` (current guard rejects asg). `pchipStillRejected` should already pass.

- [ ] **Step 3: Lift the guard in `src/+gdsge/+codegen/generateCxx.m`**

Replace lines 32-36:
```matlab
if useSympy
    if ~strcmp(ir.options.interpMethod, 'spline')
        error('gdsge:codegen:sympyInterpUnsupported', ...
            'sympy backend supports only spline interpolation (got %s)', ir.options.interpMethod);
    end
```
with:
```matlab
if useSympy
    if ~ismember(ir.options.interpMethod, {'spline', 'asg'})
        error('gdsge:codegen:sympyInterpUnsupported', ...
            'sympy backend supports spline or asg interpolation (got %s)', ir.options.interpMethod);
    end
```
(The pchip `gdsge:codegen:unsupported` error at lines 12-14 runs *before* this block, so pchip still fails fast on either backend.)

- [ ] **Step 4: Run the test, verify it passes**

Run *RUN(`codegen/tGenerateCxxAsgSympy`)*.
Expected: PASS (both methods).

- [ ] **Step 5: Commit**

```
git add src/+gdsge/+codegen/generateCxx.m tests/codegen/tGenerateCxxAsgSympy.m
git commit -m "feat(sympy-asg): accept interpMethod=asg under the sympy backend"
```

---

## Task 4: CaoKS2016 SymPy end-to-end gate (first functional proof)

Drive the public API end-to-end on the SymPy backend and match the existing autodiff golden. This is the first build+converge of a SymPy ASG MEX — it gates the `eval_vec_with_grad` C++ method.

**Files:**
- Test: `tests/CaoKS2016/codegen/tEndToEndCaoKS2016Sympy.m` (create)
- Reference (read): `tests/CaoKS2016/codegen/tEndToEndCaoKS2016.m` (the autodiff gate to mirror), `tests/CaoNie2016/codegen/tEndToEndCaoNie2016Sympy.m` (the `UseAutoDiff=0` prepend pattern)

- [ ] **Step 1: Write the failing test**

Create `tests/CaoKS2016/codegen/tEndToEndCaoKS2016Sympy.m` — the autodiff gate's structure with the SymPy prepend and the same golden:
```matlab
classdef tEndToEndCaoKS2016Sympy < matlab.unittest.TestCase
    % Phase 8b GATE (CaoKS2016, sympy backend): ASG under UseAutoDiff=0
    % converges to the same golden as the autodiff path. Evaluation-based ASG
    % interpolant comparison (tiny drift changes which grids refine). Slow.
    properties (Constant)
        RelTol = 1e-3; AbsTol = 1e-3;
    end
    methods (TestClassSetup)
        function gate(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', '..', 'src')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', '..', 'src', 'kernels')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..')));     % tests/ -> gdsgetest.*
            tc.assumeTrue(gdsgetest.sympyAvailable());
        end
    end
    methods (Test, TestTags = {'Slow'})
        function sympyBackendMatchesGolden(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            gmod = fileread(fullfile(modelDir, 'CaoKS2016.gmod'));
            gmod = ['UseAutoDiff=0;' newline gmod];
            fid = fopen(fullfile(work, 'CaoKS2016.gmod'), 'w'); fwrite(fid, gmod); fclose(fid);

            ir = gdsge_codegen('CaoKS2016');
            tc.assertEqual(ir.options.jacobianBackend, 'sympy');
            tc.assertTrue(exist(fullfile(work, ['mex_CaoKS2016.' mexext]), 'file') == 3, ...
                'sympy MEX did not compile');

            IterOptions.PrintFreq = 10;
            IterOptions.SaveFreq = 100;
            IterRslt = iter_CaoKS2016(IterOptions);

            G = load(fullfile(modelDir, 'golden', 'IterRslt.mat')).IterRslt;
            tc.verifyLessThan(IterRslt.Metric, 1e-4);
            tc.verifyLessThan(abs(IterRslt.Iter - G.Iter), 0.2*G.Iter + 20, ...
                sprintf('Iter=%d vs golden %d', IterRslt.Iter, G.Iter));
            assertInterpClose(tc, IterRslt.asg_interp_struct, G.asg_interp_struct, ...
                tc.RelTol, tc.AbsTol);
            assertInterpClose(tc, IterRslt.asg_output_struct, G.asg_output_struct, ...
                tc.RelTol, tc.AbsTol);
        end
    end
end

function assertInterpClose(tc, newStruct, goldStruct, relTol, absTol)
A = asg.construct_from_struct(newStruct);
B = asg.construct_from_struct(goldStruct);
grids = B.get_grids_info();
for j = 1:numel(grids)
    g = grids{j};
    if isempty(g); continue; end
    idx = j*ones(1, size(g, 2));
    va = A.eval_vec(idx, g);
    vb = B.eval_vec(idx, g);
    r = gdsgetest.compareNumericClose(va, vb, relTol, absTol);
    tc.verifyTrue(r.pass, sprintf('interp eval mismatch, shock %d: %s', j, ...
        strjoin(r.failures, newline)));
end
end
```

- [ ] **Step 2: Run the test, verify it fails (or errors) before the implementation is correct**

Run *RUN(`CaoKS2016/codegen/tEndToEndCaoKS2016Sympy`)*.
Expected (if run before Tasks 2–3 are merged): error at `gdsge_codegen` (sympyInterpUnsupported) or compile failure. After Tasks 2–3, this is the real gate — it should compile and converge. If the Jacobian transpose/scaling is wrong, expect non-convergence (`Metric` not < 1e-4) or a wrong `Iter` — that is the signal to recheck §5 of the spec.

- [ ] **Step 3: (Implementation already in Tasks 2–3.) Diagnose only if it fails.**

If it does not converge: the likely cause is the §5 gradient contract. Re-verify in `asgFrag` (Task 2 Step 4) that the transpose is `GDSGE_grad[v + NVEC*d] = native[v*DIM+d]` and that `eval_vec_with_grad` divides by `stateRange`. Use Task 6's cross-check to localize. Do not loosen tolerances to force a pass.

- [ ] **Step 4: Run the test, verify it passes**

Run *RUN(`CaoKS2016/codegen/tEndToEndCaoKS2016Sympy`)*.
Expected: PASS (Metric < 1e-4, Iter within band, interpolants match golden within 1e-3).

- [ ] **Step 5: Commit**

```
git add tests/CaoKS2016/codegen/tEndToEndCaoKS2016Sympy.m
git commit -m "test(sympy-asg): CaoKS2016 end-to-end gate on the sympy backend"
```

---

## Task 5: Bianchi2011_asg SymPy end-to-end gate

Same pattern for the two-stage ASG model. Mirror the autodiff gate `tests/Bianchi2011_asg/codegen/tEndToEndBianchi2011.m` (read it first — it drives a **two-stage** solve with WarmUp/SkipModelInit/MaxIter staging) and add the `UseAutoDiff=0;` prepend + `assumeTrue(sympyAvailable)`, comparing to the same golden.

**Files:**
- Test: `tests/Bianchi2011_asg/codegen/tEndToEndBianchi2011Sympy.m` (create)
- Reference (read): `tests/Bianchi2011_asg/codegen/tEndToEndBianchi2011.m`

- [ ] **Step 1: Read the autodiff gate to copy its exact two-stage driver**

Run: `sed -n '1,200p' tests/Bianchi2011_asg/codegen/tEndToEndBianchi2011.m` (or open it). Note the exact stage-1/stage-2 option structs, golden field names, and tolerances — copy them verbatim into the new test (do not re-derive the staging).

- [ ] **Step 2: Write the failing test**

Create `tests/Bianchi2011_asg/codegen/tEndToEndBianchi2011Sympy.m` with the same class skeleton as Task 4 (the `TestClassSetup` path fixtures + `assumeTrue(sympyAvailable)`), but:
- prepend `UseAutoDiff=0;` to the gmod exactly as in Task 4 Step 1,
- assert `ir.options.jacobianBackend == 'sympy'`,
- run the **identical two-stage driver and golden assertions copied from `tEndToEndBianchi2011.m`** (same `IterOptions`/WarmUp staging, same `assertInterpClose`-style evaluation-based comparison, same tolerances).

- [ ] **Step 3: Run the test, verify it fails before / passes after**

Run *RUN(`Bianchi2011_asg/codegen/tEndToEndBianchi2011Sympy`)*.
Expected after Tasks 2–3: PASS (both stages converge, interpolants match the golden within the autodiff gate's tolerances).

- [ ] **Step 4: Commit**

```
git add tests/Bianchi2011_asg/codegen/tEndToEndBianchi2011Sympy.m
git commit -m "test(sympy-asg): Bianchi2011_asg two-stage end-to-end gate on the sympy backend"
```

---

## Task 6: SymPy↔adept↔FD Jacobian cross-check (ASG)

The strongest correctness gate: build both MEX, evaluate the analytic Jacobian (6th output, `GDSGE_DEBUG_EVAL_ONLY==2`) at a converged solution on the same grid problems, and assert **sympy == adept == finite-difference**. Mirrors `tSympyJacobianHL1996.analyticJacMatchesAdeptAndFD` but with an ASG handle. The MEX caller-workspace contract is taken verbatim from `gdsge.runtime.solveProblemsAsg` (the ASG path uses `GDSGE_ASG_HANDLE`, sets `GDSGE_SPLINE_VEC=[]`).

**Files:**
- Test: `tests/CaoKS2016/codegen/tSympyJacobianCaoKS2016.m` (create)
- Reference (read): `tests/HeatonLucas1996/codegen/tSympyJacobianHL1996.m` (the spline cross-check to mirror), `src/+gdsge/+runtime/solveProblemsAsg.m:28-47` (the ASG caller contract), and the generated `iter_CaoKS2016.m` (the data/bounds assembly — available from Task 4's compile)

- [ ] **Step 1: Write the failing test (harness + a single converged batch)**

Create `tests/CaoKS2016/codegen/tSympyJacobianCaoKS2016.m`. The harness builds inputs from a converged **autodiff** run (its `asg_interp_struct` is the value-fn interp; its `sol_asg_interp_struct` gives the solution at the refined grids), then calls both MEX with 6 outputs:
```matlab
classdef tSympyJacobianCaoKS2016 < matlab.unittest.TestCase
    % Phase 8b cross-check: ASG analytic Jacobian (sympy) == adept == finite
    % difference, at a converged solution, per grid problem. Slow.
    methods (TestClassSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', '..', 'src')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..', '..', 'src', 'kernels')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(here, '..', '..')));
            tc.assumeTrue(gdsgetest.sympyAvailable());
        end
    end
    methods (Test, TestTags = {'Slow'})
        function analyticJacMatchesAdeptAndFD(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'CaoKS2016.gmod'), work);
            gdsge.runtime.ensurePath();

            % --- build BOTH MEX (autodiff + sympy) ---
            ir = gdsge_codegen('CaoKS2016');                 % autodiff build in work
            copyfile(fullfile(work, ['mex_CaoKS2016.' mexext]), ...
                     fullfile(work, ['mex_CaoKS2016_ad.' mexext]));
            sir = ir; sir.options.jacobianBackend = 'sympy';
            sdir = fullfile(work, 'sy'); mkdir(sdir);
            gdsge.codegen.generateCxx(sir, sdir);
            oldCd = pwd; cd(sdir); compile_CaoKS2016(); cd(oldCd);
            tc.assertTrue(exist(fullfile(sdir, ['mex_CaoKS2016.' mexext]),'file')==3, 'sympy MEX');
            copyfile(fullfile(sdir, ['mex_CaoKS2016.' mexext]), ...
                     fullfile(work, ['mex_CaoKS2016_sy.' mexext]));

            % --- a converged run supplies the asg handle + grids + solution ---
            cd(work); tc.addTeardown(@() cd(oldCd));
            IterRslt = iter_CaoKS2016(struct('PrintFreq',10,'SaveFreq',inf,'NoSave',1));
            S = gdsgetest.buildAsgMexInputs(IterRslt, ir);   % helper, Step 2

            Jad = evalJac(@mex_CaoKS2016_ad, S);
            Jsy = evalJac(@mex_CaoKS2016_sy, S);
            Jfd = fdJac(@mex_CaoKS2016_ad, S);

            r = gdsgetest.compareNumericClose(Jsy, Jad, 1e-6, 1e-8);
            tc.verifyTrue(r.pass, ['sympy vs adept: ' strjoin(r.failures, newline)]);
            r = gdsgetest.compareNumericClose(Jsy, Jfd, 1e-4, 1e-6);
            tc.verifyTrue(r.pass, ['sympy vs finite-diff: ' strjoin(r.failures, newline)]);
        end
    end
end

function J = evalJac(mexFn, S)
cfg = S.cfg; cfg.debugEvalOnly = 2;
J = callMex6Asg(mexFn, S, cfg);
end

function J = fdJac(mexFn, S)
nx = size(S.sol0,1); np = size(S.sol0,2); h = 1e-6;
J = zeros(nx*nx, np);
for var = 1:nx
    sp = S.sol0; sp(var,:) = sp(var,:) + h; ep = residAsg(mexFn, S, sp);
    sm = S.sol0; sm(var,:) = sm(var,:) - h; em = residAsg(mexFn, S, sm);
    d = (ep - em) / (2*h);
    for eq = 1:nx
        J((eq-1) + nx*(var-1) + 1, :) = d(eq,:);
    end
end
end

function eqval = residAsg(mexFn, S, sol)
cfg = S.cfg; cfg.debugEvalOnly = 1;
[~,~,~,eqval] = callMexAsg(mexFn, S, sol, cfg);   % 4th output = eqval
end
```
Add two local caller helpers `callMex6Asg` (6-output, eval mode 2) and `callMexAsg` (eval mode 1) that set the **ASG caller-workspace contract** copied verbatim from `solveProblemsAsg.m:28-42` (`TolFun/TolSol/SolMaxIter/NumThreads/GDSGE_DEBUG_EVAL_ONLY/UseBroyden/FiniteDiffDelta/GDSGE_USE_BROYDEN_NOW=0/MEX_TASK_NAME=cfg.taskName/MEX_TASK_INIT=0/MEX_TASK_INF_HORIZON=1/GDSGE_SPLINE_VEC=[]/GDSGE_ASG_HANDLE=S.cfg.asgHandle`), then call `mexFn(sol, S.lb, S.ub, S.data, skip, S.f0, S.aux0, S.eqval0)` with the right nlhs (`skip=zeros(1,np)` so all problems evaluate). Model these two helpers on `callMex6` in `tSympyJacobianHL1996.m` (same idea, ASG contract).

- [ ] **Step 2: Write the input-builder helper `gdsgetest.buildAsgMexInputs`**

Create `tests/+gdsgetest/buildAsgMexInputs.m`. It returns `S` with fields `sol0, lb, ub, data, f0, aux0, eqval0, cfg` for the converged grid problems. Build it by reusing the converged interpolants and **mirroring the data/bounds assembly in the generated `iter_CaoKS2016.m`** (open it from the Task-4 work dir; locate the block that builds `GDSGE_DATA`, the bounds `GDSGE_LB/GDSGE_UB`, and `GDSGE_ASG_HANDLE`):
```matlab
function S = buildAsgMexInputs(IterRslt, ir)
% Build one converged ASG MEX batch for the Jacobian cross-check.
A    = asg.construct_from_struct(IterRslt.asg_interp_struct);      % value-fn interp
Asol = asg.construct_from_struct(IterRslt.sol_asg_interp_struct);  % solution interp
grids = A.get_grids_info();                 % {1..numArray}, each [numDim x np_j]
% Stack all shock arrays' grid points into one batch, tagging each with its shock idx.
arrayIdx = []; pts = [];
for j = 1:numel(grids)
    g = grids{j};
    if isempty(g); continue; end
    arrayIdx = [arrayIdx, j*ones(1, size(g,2))]; %#ok<AGROW>
    pts = [pts, g]; %#ok<AGROW>
end
np = size(pts, 2);
S.sol0 = Asol.eval_vec(arrayIdx, pts);      % [numSol x np] converged solution = warm start
% data column = [params (dataLayout order); state coords]. MIRROR the generated
% iter_CaoKS2016.m GDSGE_DATA assembly: params from ir setup replicated across
% np columns, then pts (the state coords) stacked beneath, in state order.
S.data = gdsgetest.asgDataPack(ir, pts);    % helper below, or inline the iter file's pack
% bounds: from the model inbound; MIRROR the generated iter's GDSGE_LB/GDSGE_UB
% (adaptive bounds widen per point — for a Jacobian-at-a-point check, the static
% inbound bounds are sufficient since we evaluate, not solve).
[S.lb, S.ub] = gdsgetest.asgBounds(ir, S.sol0);
S.f0    = 1e20 * ones(1, np);
S.aux0  = zeros(sum(cellfun(@(a) a.length, ir.variables.aux)) + 1, np);
S.eqval0 = zeros(size(S.sol0));
S.cfg = struct('tolFun',1e-10,'tolSol',1e-6,'solMaxIter',200,'numThreads',1, ...
    'useBroyden',0,'finiteDiffDelta',1e-6,'useBroydenNow',0, ...
    'taskName','inf_horizon','asgHandle',A.objectHandle);
end
```
`asgDataPack` and `asgBounds` are thin helpers that reproduce the exact `GDSGE_DATA`/`GDSGE_LB`/`GDSGE_UB` construction found in the generated `iter_CaoKS2016.m`. Copy those few lines directly (param vector order, state stacking, bound expansion) — do not re-derive the layout.

- [ ] **Step 3: Run the test, verify it fails first**

Run *RUN(`CaoKS2016/codegen/tSympyJacobianCaoKS2016`)*.
Expected before the helpers are correct: an error (missing helper) or a Jacobian mismatch. Iterate on `buildAsgMexInputs` until `sympy == adept` (a clean signal the harness inputs are valid — both backends consume the identical contract).

- [ ] **Step 4: Run the test, verify it passes**

Run *RUN(`CaoKS2016/codegen/tSympyJacobianCaoKS2016`)*.
Expected: PASS — `sympy == adept` (rtol 1e-6) and `sympy == finite-diff` (rtol 1e-4), per grid problem.

> **Fallback (if the batch reconstruction proves brittle):** the three end-to-end gates (Tasks 4, 5, 10) already exercise this exact Jacobian through real convergence to the goldens. If `buildAsgMexInputs` cannot be made robust in a reasonable budget, descope this task to a follow-up and rely on the end-to-end gates for functional correctness — note the descope in `PROGRESS.md`. Do **not** weaken the cross-check tolerances to force a pass.

- [ ] **Step 5: Commit**

```
git add tests/CaoKS2016/codegen/tSympyJacobianCaoKS2016.m tests/+gdsgetest/buildAsgMexInputs.m tests/+gdsgetest/asgDataPack.m tests/+gdsgetest/asgBounds.m
git commit -m "test(sympy-asg): SymPy<->adept<->FD Jacobian cross-check on CaoKS2016 (ASG)"
```

---

## Task 7: Create the CaoNie2016_asg corpus model (`.gmod`)

The capstone model: the proven cartesian `CaoNie2016` body with the interpolation axis flipped to ASG. First model combining ASG × conditional regions/equations.

**Files:**
- Create: `tests/CaoNie2016_asg/CaoNie2016.gmod`
- Reference (read): `tests/CaoNie2016/CaoNie2016.gmod` (the cartesian twin), `tests/CaoKS2016/CaoKS2016.gmod:10-14` (canonical new-pipeline ASG option names)

- [ ] **Step 1: Copy the cartesian gmod and flip the interpolation options**

Copy `tests/CaoNie2016/CaoNie2016.gmod` → `tests/CaoNie2016_asg/CaoNie2016.gmod`, then edit the options block (lines 13-16 of the cartesian file) from:
```
USE_SPLINE = 1;
USE_ASG = 0;
AsgThreshhold = 1e-4;
InterpOrder = 4;
```
to (use the **canonical** new-pipeline option names whitelisted in `resolveOptions.m:16-21` — `AsgThreshhold` with the double-h typo is NOT whitelisted and is only silently ignored for spline):
```
USE_SPLINE = 0;
USE_ASG = 1;
AsgMaxLevel = 10;
AsgThreshold = 1e-4;
```
Leave the entire model body, `var_others shock_trans;`, the empty `var_tensor;`, the two `model(...)` regions, the `if/else` equations, and the `simulate` block byte-for-byte unchanged. Keep `TolEq`/`PrintFreq`/`SaveFreq`/`NumThreads` as in the cartesian file.

- [ ] **Step 2: Sanity-parse on the new pipeline (no compile)**

Run:
```
matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile('..','src')); addpath(fullfile('..','src','kernels')); w=tempname; mkdir(w); copyfile('CaoNie2016_asg/CaoNie2016.gmod', w); cd(w); ir=gdsge.parser.parseFrontEnd('CaoNie2016.gmod'); fprintf('interpMethod=%s #regions=%d #interp=%d #aux=%d\n', ir.options.interpMethod, numel(ir.model.regions), numel(ir.interp), numel(ir.variables.aux));"
```
Expected: `interpMethod=asg #regions=2 #interp=3 ...` with no error (the empty `var_tensor;` and the conditional regions parse). If `AsgThreshold`/`AsgMaxLevel` raise an option-whitelist error, fix the names to match `resolveOptions.m`.

- [ ] **Step 3: Commit**

```
git add tests/CaoNie2016_asg/CaoNie2016.gmod
git commit -m "test(sympy-asg): add CaoNie2016_asg corpus gmod (ASG variant of CaoNie2016)"
```

---

## Task 8: Capture the CaoNie2016_asg golden (old toolbox)

Reproducible golden capture from the vendored old toolbox, mirroring `tests/golden/capture_CaoNie2016.m` / `capture_CaoKS2016.m`. One MATLAB process; never concurrent.

**Files:**
- Create: `tests/golden/capture_CaoNie2016_asg.m`
- Create (output): `tests/CaoNie2016_asg/golden/IterRslt.mat`, `tests/CaoNie2016_asg/golden/SimuRslt.mat`
- Reference (read): `tests/golden/capture_CaoNie2016.m`

- [ ] **Step 1: Write the capture script**

Create `tests/golden/capture_CaoNie2016_asg.m` modeled on `capture_CaoNie2016.m`, pointing at `tests/CaoNie2016_asg`:
```matlab
function capture_CaoNie2016_asg()
% Golden capture for CaoNie2016_asg from the OLD toolbox. Same model body as
% the cartesian CaoNie2016, ASG interpolation. Reduced seeded simulate 6x1000.
here     = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(here));
oldSrc   = fullfile(repoRoot, 'base_package', 'gdsge', 'source');
modelDir = fullfile(repoRoot, 'tests', 'CaoNie2016_asg');
goldenDir = fullfile(modelDir, 'golden');
if ~exist(goldenDir, 'dir'); mkdir(goldenDir); end

work = tempname; mkdir(work);
copyfile(fullfile(modelDir, 'CaoNie2016.gmod'), work);

oldPath = path; restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(oldSrc);
oldCd = pwd; cdRestore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(work);

t0 = tic;
gdsge_codegen('CaoNie2016');
IterOptions.PrintFreq = 10;
IterOptions.SaveFreq = 100;
IterRslt = iter_CaoNie2016(IterOptions);
fprintf('iter wall-clock: %.1fs (Iter=%d Metric=%g)\n', toc(t0), IterRslt.Iter, IterRslt.Metric);

rng(0823);
SimuRslt = simulate_CaoNie2016(IterRslt, struct('num_samples', 6, 'num_periods', 1000));

save(fullfile(goldenDir, 'IterRslt.mat'), 'IterRslt', '-v7');
save(fullfile(goldenDir, 'SimuRslt.mat'), 'SimuRslt', '-v7');
fprintf('Golden captured: Iter=%d Metric=%g\n', IterRslt.Iter, IterRslt.Metric);
end
```

- [ ] **Step 2: Run the capture**

Run: `matlab -batch "addpath('tests/golden'); capture_CaoNie2016_asg"`
Expected: console shows convergence (the base_package reference is `IterRslt_CaoNie2016_1077.mat` → expect Iter≈1077, Metric<1e-6) and writes both `.mat` files. If the **old toolbox** rejects `var_others` together with ASG, drop the `var_others shock_trans;` line from the gmod (Task 7) and re-run — the cartesian golden carries `var_others.shock_trans`, but if the ASG old-toolbox path can't, the new-pipeline gate simply won't compare that field (note it).

- [ ] **Step 3: Verify golden size and commit**

Run: `ls -la tests/CaoNie2016_asg/golden/` — confirm both files exist and are reasonably sized (≤ a few MB; if `IterRslt.mat` is large, trim uncompared bulk fields per `tests/golden/trim_mendoza_golden.m` precedent before committing).
```
git add tests/golden/capture_CaoNie2016_asg.m tests/CaoNie2016_asg/golden/IterRslt.mat tests/CaoNie2016_asg/golden/SimuRslt.mat
git commit -m "test(sympy-asg): capture CaoNie2016_asg golden from the old toolbox"
```

---

## Task 9: CaoNie2016_asg autodiff end-to-end gate + IR snapshot

Bring the new model up on the **autodiff** ASG path first — this is where any ASG × conditional-region wiring bug surfaces (independent of SymPy). Also locks the IR JSON snapshot.

**Files:**
- Create: `tests/CaoNie2016_asg/codegen/tEndToEndCaoNie2016_asg.m`
- Create: `tests/CaoNie2016_asg/ir/CaoNie2016_asg.gdsge.json` (generated, committed)
- Reference (read): `tests/CaoKS2016/codegen/tEndToEndCaoKS2016.m` (ASG evaluation-based comparison + driver), `tests/CaoNie2016/codegen/tEndToEndCaoNie2016.m` (the cartesian twin's golden fields: `var_others`, simulate `X`)

- [ ] **Step 1: Write the failing test (autodiff)**

Create `tests/CaoNie2016_asg/codegen/tEndToEndCaoNie2016_asg.m` modeled on `tEndToEndCaoKS2016.m` (artifacts exist, IR round-trips, second-run compile-skip, iter converges, **evaluation-based** ASG interpolant comparison via `assertInterpClose`), plus the CaoNie2016-specific bits from the cartesian twin: pin `rng(0823)` before iter (multi-root KKT restart determinism — see the safe_assets RNG-restart memory), compare `IterRslt.var_others.shock_trans` if present, and a seeded `simulate` comparing `SimuRslt.shock` + `SimuRslt.X`. Use `RelTol=AbsTol=1e-3` for the ASG eval comparison. Tag `Slow`.

- [ ] **Step 2: Run the test, verify it fails first**

Run *RUN(`CaoNie2016_asg/codegen/tEndToEndCaoNie2016_asg`)*.
Expected: FAIL initially — `tests/CaoNie2016_asg/ir/` snapshot does not exist yet, and/or the ASG×region intersection surfaces a wiring bug (e.g. the `model(X==0)` region or an `if/else` equation mis-emitted under ASG). Debug per `superpowers:systematic-debugging` — the model body is proven on cartesian, so suspect the ASG driver (`emitIterAsg`) region handling or the ASG interp section, not the parser.

- [ ] **Step 3: Generate and commit the IR snapshot**

After the model compiles, write the IR JSON snapshot (the public API already emits `CaoNie2016.gdsge.json` into the work dir; copy it to the committed location, renamed):
```
matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile('..','src')); addpath(fullfile('..','src','kernels')); w=tempname; mkdir(w); copyfile('CaoNie2016_asg/CaoNie2016.gmod', w); cd(w); gdsge_codegen('CaoNie2016'); copyfile('CaoNie2016.gdsge.json', fullfile('$REPO','tests','CaoNie2016_asg','ir','CaoNie2016_asg.gdsge.json'));"
```
(Replace `$REPO` with the absolute repo path. Create `tests/CaoNie2016_asg/ir/` first.) Have the test assert round-trip: `decoded = gdsge.ir.decode(fileread(<work>/CaoNie2016.gdsge.json)); verifyTrue(isequalIR(decoded, ir))`.

- [ ] **Step 4: Run the test, verify it passes**

Run *RUN(`CaoNie2016_asg/codegen/tEndToEndCaoNie2016_asg`)*.
Expected: PASS — converges (Iter≈golden, Metric<1e-6), interpolants match golden within 1e-3, simulate shock+X match.

- [ ] **Step 5: Commit**

```
git add tests/CaoNie2016_asg/codegen/tEndToEndCaoNie2016_asg.m tests/CaoNie2016_asg/ir/CaoNie2016_asg.gdsge.json
git commit -m "test(sympy-asg): CaoNie2016_asg autodiff end-to-end gate + IR snapshot"
```

---

## Task 10: CaoNie2016_asg SymPy end-to-end gate (the capstone)

The full intersection: ASG × conditional regions/equations on the SymPy backend.

**Files:**
- Create: `tests/CaoNie2016_asg/codegen/tEndToEndCaoNie2016_asgSympy.m`
- Reference (read): the Task 9 autodiff gate + `tests/CaoNie2016/codegen/tEndToEndCaoNie2016Sympy.m` (the `UseAutoDiff=0` prepend)

- [ ] **Step 1: Write the failing test**

Create `tests/CaoNie2016_asg/codegen/tEndToEndCaoNie2016_asgSympy.m`: the Task 9 test's body with the SymPy class skeleton (path fixtures + `assumeTrue(sympyAvailable)`), the `UseAutoDiff=0;` gmod prepend, `assertEqual(ir.options.jacobianBackend,'sympy')`, the same `rng(0823)` pin, and the same golden assertions (evaluation-based ASG comparison + simulate). Tag `Slow`.

- [ ] **Step 2: Run the test, verify it fails before / passes after**

Run *RUN(`CaoNie2016_asg/codegen/tEndToEndCaoNie2016_asgSympy`)*.
Expected after Tasks 2–3 + 9: PASS. If the SymPy build fails to converge while autodiff (Task 9) converged, the suspect is the SymPy per-region/per-branch guarded Jacobian *combined with* the ASG interp grad — use Task 6's cross-check harness pointed at CaoNie2016_asg to localize.

- [ ] **Step 3: Commit**

```
git add tests/CaoNie2016_asg/codegen/tEndToEndCaoNie2016_asgSympy.m
git commit -m "test(sympy-asg): CaoNie2016_asg SymPy end-to-end gate (capstone)"
```

---

## Task 11: Docs, deferred-map update, and full-suite gate

**Files:**
- Modify: `docs/deferred-features.md` (the `sympyInterpUnsupported` row)
- Modify: `PROGRESS.md` (mark Phase 8b ASG slice done; changelog entry)

- [ ] **Step 1: Update the deferred-features map**

In `docs/deferred-features.md`, edit the `gdsge:codegen:sympyInterpUnsupported` row (currently "ASG or pchip") to cover **pchip only**, and update the "Open deferred sub-phases" Phase 8b line to "pchip only (ASG landed 2026-06-15)":
```
| `gdsge:codegen:sympyInterpUnsupported` | `+codegen/generateCxx.m` | SymPy backend (`UseAutoDiff=0`) with the **pchip** interpolant | Use the autodiff backend (ASG now supported); pchip is a separate gap. |
```

- [ ] **Step 2: Update PROGRESS.md**

Flip the Phase 8b status (the deferred bullet at PROGRESS.md:165-167) to note the ASG half is done, and add a `## Changelog` entry dated 2026-06-15 summarizing: `eval_vec_with_grad` + SymPy ASG interp emitter + guard lift; gates CaoKS2016/Bianchi2011_asg/CaoNie2016_asg (new corpus model, ASG × conditional regions, both backends) + Jacobian cross-check; pchip remains deferred. Reference the spec + this plan + branch `phase8b-sympy-asg`. Convert "Phase 8b (deferred)" to "Phase 8b — ASG done; pchip deferred".

- [ ] **Step 3: Run the full suite**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: exit 0; `tests/results/junit.xml` shows all pass (the prior 491 + the new SymPy-ASG and CaoNie2016_asg tests; SymPy tests run only if `sympyAvailable`). Confirm no autodiff golden regressed (the existing ASG autodiff gates `tEndToEndCaoKS2016`/`tEndToEndBianchi2011` must stay green — the `asg_adouble.h` addition is additive).

- [ ] **Step 4: Commit**

```
git add docs/deferred-features.md PROGRESS.md
git commit -m "docs(sympy-asg): mark Phase 8b ASG slice done; pchip remains the only sympy-interp deferral"
```

---

## Self-review (completed)

**Spec coverage:**
- §3 enabling fact / §4.1 C++ method → Task 2 (Step 3).
- §4.2 SymPy ASG emitter → Task 2 (Step 4) + unit test.
- §4.3 guard lift → Task 3.
- §5 gradient contract (scaling + transpose) → encoded in Task 2 Step 3/4; pinned by Task 6.
- §6 CaoNie2016_asg (gmod, golden, autodiff-first, IR snapshot, SymPy) → Tasks 7–10.
- §7 verification: cross-check → Task 6; three end-to-end gates → Tasks 4, 5, 10; Python gating → every SymPy test's `assumeTrue`; structural snapshot → Task 3.
- §10 pchip stays deferred → Task 3 (guard keeps pchip rejected) + Task 11 (docs).

**Placeholder scan:** the only non-literal step is Task 6's data/bounds pack, which is explicitly "copy the exact lines from the generated `iter_CaoKS2016.m`" (a concrete, available source) with a stated fallback — not an open-ended TODO.

**Type/name consistency:** `eval_vec_with_grad` (C++ method), `GDSGE_INTERP_VEC_double_grad` (lambda, shared with the spline path), `GDSGE_ASG_GRAD_NATIVE` (scratch), `asgFrag` (local fn), `buildAsgMexInputs`/`asgDataPack`/`asgBounds` (test helpers), `A.objectHandle` (asg handle) — used consistently across tasks. Gradient layouts: native `[v*DIM+d]`, spline-consumed `[v + NVEC*d]` — consistent in spec §5 and Task 2.
