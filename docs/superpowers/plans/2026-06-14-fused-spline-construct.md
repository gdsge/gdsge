# Fused In-MEX Spline Construction — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the per-iteration cartesian spline/linear rebuild (`myppual` per var + `convert_to_interp_eval_array` permute/concat) with one minimal, self-contained `interp_construct_mex` call that constructs all interpolants and writes `GDSGE_SPLINE_VEC` directly in eval order — bit-exact with today.

**Architecture:** A new construction-only MEX (`src/kernels/interp_construct_mex.cpp`) copies the verbatim construction driver from the owner's `myppual_mex.cpp` (`mkl_start` / `mkl_start_with_extrap` + helpers, over `mkl_dummy_interp.h`'s `construct_cubic_spline_notaknot` / `construct_linear_interp`), then (a) returns per-variable natural-order pp structs identical to `myppual(pp)`, and (b) packs the coefficients directly into the eval-order layout that `convert_to_interp_eval_array` produces. `gdsge.runtime.constructSplines` calls it once (with all `Values` passed via a struct). The C++ solver MEX, the simulate path (keeps the prebuilt `myppual_mex`), ASG, and pchip are untouched.

**Tech Stack:** C++ MEX (MSVC, OpenMP, `/fp:precise`), MATLAB R2025b, `matlab.unittest`.

**Reference spec:** `docs/superpowers/specs/2026-06-14-fused-spline-construct-design.md`.

**Source of the verbatim driver:** `C:\Dropbox\Finite_Agents_General_Toolkit\update_mac\gdsge\source\myppual_mex.cpp` (same machine). Key functions and lines: `fullConstruction` (368–423), `fullConstructionWithExtrap` (298–367), `mkl_start` (2647–2739), `mkl_start_with_extrap` (2476–2646), `intmemcpy`/`doublememcpy`/`vectorProd` (2740–2760). The MATLAB MKL-construct branch it mirrors: `src/kernels/myppual.m` lines 513–646 (esp. 590–642: extrap break-extension, `mkl_start[_with_extrap]`, the returned `MKLpp` struct shape).

> **Bit-exactness rule (applies to every C++ step):** never alter the copied driver's arithmetic or loop order. Compile with MSVC default `/fp:precise` (never `/fp:fast`). The Task 2/3 driver tests are the oracle: coefs must equal `myppual`/`convert_to_interp_eval_array` to the **last bit** (`isequal`, not `verifyEqual` with tolerance).

---

## Task 0: Branch & baseline

**Files:** none (already on branch `fused-spline-construct`).

- [ ] **Step 1: Confirm branch and clean tree**

Run: `git status -sb`
Expected: `## fused-spline-construct`, no unstaged changes except this plan.

- [ ] **Step 2: Record the green baseline (kernel suite)**

Run: `matlab -batch "cd('tests'); results = runtests('tKernels'); disp(table(results)); assert(all([results.Passed]))"`
Expected: PASS. (Sanity that the kernels load before we add one.)

---

## Task 1: Compile gate `ensureSplineConstructMex`

Mirror `ensureAsgMex` so the new kernel auto-compiles (cache-gated) and is invoked at codegen entry.

**Files:**
- Create: `src/+gdsge/+codegen/ensureSplineConstructMex.m`
- Test: `tests/codegen/tEnsureSplineConstructMex.m`

- [ ] **Step 1: Write the failing test**

```matlab
classdef tEnsureSplineConstructMex < matlab.unittest.TestCase
    methods (Test)
        function compilesAndIsCallable(tc)
            gdsge.codegen.ensureSplineConstructMex();
            tc.verifyEqual(exist('interp_construct_mex', 'file'), 3, ...
                'interp_construct_mex MEX not on path after ensure');
        end
        function secondCallNoOps(tc)
            gdsge.codegen.ensureSplineConstructMex();   % warm
            t = tic; gdsge.codegen.ensureSplineConstructMex(); dt = toc(t);
            tc.verifyLessThan(dt, 2.0, 'cache-gated re-call should be instant');
        end
    end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `matlab -batch "cd('tests'); r = runtests('tEnsureSplineConstructMex'); assert(all([r.Passed]))"`
Expected: FAIL — `Unrecognized function 'gdsge.codegen.ensureSplineConstructMex'` (the .cpp does not exist yet either; Task 2 adds it). This task lands the .m; the test goes green only after Task 2's .cpp exists. **Mark this test `@Disabled`-equivalent by leaving it failing until Task 2** — note in the commit. (Alternatively run it at the end of Task 2.)

- [ ] **Step 3: Write `ensureSplineConstructMex.m`** (copy of `ensureAsgMex` with names + `myppual` compile flags)

```matlab
function ensureSplineConstructMex()
% ENSURESPLINECONSTRUCTMEX  Compile src/kernels/interp_construct_mex.cpp when
%   needed (hash cache interp_construct_mex.cache beside it, same mechanism as
%   ensureAsgMex / mex_<model>.cache). The kernel is the fused spline/linear
%   constructor consumed by gdsge.runtime.constructSplines. Flags mirror the
%   owner's compile_myppual.m (/openmp, -DUSE_OMP) plus /O2; /fp:precise is the
%   MSVC default and MUST stay (bit-exact with the prebuilt myppual_mex).
here      = fileparts(mfilename('fullpath'));   % src/+gdsge/+codegen
srcRoot   = fileparts(fileparts(here));         % src
repoRoot  = fileparts(srcRoot);
kernels   = fullfile(srcRoot, 'kernels');
cppFile   = fullfile(kernels, 'interp_construct_mex.cpp');
cacheFile = fullfile(kernels, 'interp_construct_mex.cache');
mexFile   = fullfile(kernels, ['interp_construct_mex.' mexext]);
cppText   = fileread(cppFile);
if exist(mexFile, 'file') == 3 && ~gdsge.codegen.needsCompile(cppText, cacheFile)
    return;
end
fprintf('Compiling interp_construct_mex (cache-gated):\n');
includeDir = fullfile(repoRoot, 'include');
oldCd = pwd; restore = onCleanup(@() cd(oldCd)); %#ok<NASGU>
cd(kernels);
mex(['-I' includeDir], '-DUSE_OMP', 'interp_construct_mex.cpp', ...
    'OPTIMFLAGS=/O2 /DNDEBUG', ...
    'COMPFLAGS=$COMPFLAGS /wd4267 /wd4068 /wd4091 /openmp');
gdsge.codegen.writeText(cacheFile, cppText);
end
```

- [ ] **Step 4: Commit**

```bash
git add src/+gdsge/+codegen/ensureSplineConstructMex.m tests/codegen/tEnsureSplineConstructMex.m
git commit -m "feat(codegen): ensureSplineConstructMex compile gate (mirrors ensureAsgMex)"
```

---

## Task 2: The construct-only MEX — per-variable pp output

Build `interp_construct_mex.cpp` with the verbatim driver and a first output: the per-variable natural-order pp cell, **bit-identical** to `myppual(pp)`. (Eval-order packing is Task 3.)

**Files:**
- Create: `src/kernels/interp_construct_mex.cpp`
- Test: `tests/kernels/tInterpConstructMex.m`

- [ ] **Step 1: Write the failing test (ppCell parity)**

```matlab
classdef tInterpConstructMex < matlab.unittest.TestCase
    % Bit-exact parity of the fused constructor vs myppual + convert.
    methods (TestClassSetup)
        function ensureBuilt(tc)
            gdsge.codegen.ensureSplineConstructMex();
            root = fileparts(which('essential_blas.dll'));   % myppual_mex needs it
            if ~isempty(root) && ~any(strcmp(strsplit(getenv('PATH'),';'),root))
                setenv('PATH', [getenv('PATH') ';' root]);
            end
        end
    end
    methods (Static)
        function pp = mkpp(values, breaks, order, extrap)
            % Legacy reference pp via myppual (one interp var). values is
            % [numArray, gridDims...]; breaks is a 1xD cell of grids.
            sizeState = num2cell(cellfun(@numel, breaks));
            pp = struct('form','MKL','breaks',{breaks}, ...
                'Values',reshape(values,[],sizeState{:}),'coefs',[], ...
                'order',order*ones(1,numel(breaks)),'Method',[], ...
                'ExtrapolationOrder',extrap,'thread',1,'orient','curvefit');
            pp = myppual(pp);
        end
        function vals = build(breaks, numArray)
            sz = [numArray, cellfun(@numel, breaks)];
            vals = reshape(1:prod(sz), sz);           % deterministic data
        end
    end
    methods (Test)
        function ppCellMatchesMyppual_1D_cubic(tc)
            breaks = {linspace(0,1,7)};
            v = tInterpConstructMex.build(breaks, 3);  % numArray=3
            ref = tInterpConstructMex.mkpp(v, breaks, 4, 4);
            S.a = reshape(v, [], numel(breaks{1}));    % first dim = numArray
            ppCell = interp_construct_mex(breaks, S, int32(4*ones(1,1)), int32([]), int32(1));
            tc.verifyTrue(isequal(ppCell{1}.coefs(:), ref.coefs(:)), '1D cubic coefs not bit-identical');
        end
        function ppCellMatchesMyppual_2D_cubic_extrap(tc)
            breaks = {linspace(0,1,6), linspace(-1,2,5)};
            v = tInterpConstructMex.build(breaks, 4);
            ref = tInterpConstructMex.mkpp(v, breaks, 4, 2);  % extrapOrder=2
            S.a = reshape(v, [], numel(breaks{1}), numel(breaks{2}));
            ppCell = interp_construct_mex(breaks, S, int32(4*ones(1,2)), int32(2*ones(1,2)), int32(1));
            tc.verifyTrue(isequal(ppCell{1}.coefs(:), ref.coefs(:)), '2D cubic+extrap coefs not bit-identical');
        end
        function ppCellMatchesMyppual_linear(tc)
            breaks = {linspace(0,1,8)};
            v = tInterpConstructMex.build(breaks, 2);
            ref = tInterpConstructMex.mkpp(v, breaks, 2, []);
            S.a = reshape(v, [], numel(breaks{1}));
            ppCell = interp_construct_mex(breaks, S, int32(2), int32([]), int32(1));
            tc.verifyTrue(isequal(ppCell{1}.coefs(:), ref.coefs(:)), 'linear coefs not bit-identical');
        end
    end
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `matlab -batch "cd('tests'); r = runtests('tInterpConstructMex'); assert(all([r.Passed]))"`
Expected: FAIL — `interp_construct_mex.cpp` does not exist, `ensureSplineConstructMex` errors on missing source.

- [ ] **Step 3: Create `interp_construct_mex.cpp` — header, driver, gateway (ppCell only)**

Create the file with this exact skeleton. **Copy the five driver functions verbatim** from the reference `myppual_mex.cpp` at the lines listed in the plan header (do not edit their bodies). The short helpers are reproduced inline below; paste `mkl_start` (2647–2739) and `mkl_start_with_extrap` (2476–2646) verbatim where indicated.

```cpp
// interp_construct_mex.cpp — minimal fused spline/linear constructor.
// Construction driver copied verbatim from myppual_mex.cpp (do not edit bodies).
// [out: ppCell, splineVec] = interp_construct_mex(breaks, values, order, extrapOrder, numThreads)
#include "mex.h"
#include "mkl_dummy_interp.h"   // construct_cubic_spline_notaknot / construct_linear_interp / mkl_domatcopy / dfd*
#include <string.h>
#ifdef USE_OMP
#include <omp.h>
#endif

#define MAXDIM 10
#define MAX(a,b) ((a>b) ? a : b)
#define MIN(a,b) ((a<b) ? a : b)

// --- verbatim helpers (myppual_mex.cpp 2740-2760) ---
void intmemcpy(int* des, int* src, int n) { memcpy(des, src, sizeof(int) * n); }
inline void doublememcpy(double* des, double* src, int n) { memcpy(des, src, sizeof(double) * n); }
inline int vectorProd(int* x, int n) { int r = x[0]; for (int i=1;i<n;++i) r*=x[i]; return r; }

// --- PASTE VERBATIM: mkl_start_with_extrap (myppual_mex.cpp 2476-2646) ---
// --- PASTE VERBATIM: mkl_start            (myppual_mex.cpp 2647-2739) ---

// Construct one interp variable's natural-order coefs into `coeff`.
// `valuesMKL` is the value array already permuted to MKL orientation
// (grid dims reversed, array dim last); see buildValuesMKL below. Mirrors
// fullConstruction / fullConstructionWithExtrap (myppual_mex.cpp 298-423).
static void constructOne(int xDim, int* xPts, double** xGrid, int yDim,
        double* valuesMKL, int* s_order, int* extrap_order, double* coeff, int coeffn)
{
    if (extrap_order == NULL) {
        mkl_start(xDim, xPts, xGrid, yDim, valuesMKL, s_order, coeff, coeffn);
    } else {
        mkl_start_with_extrap(xDim, xPts, xGrid, yDim, valuesMKL, s_order, coeff, coeffn, extrap_order);
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    const mxArray* breaks_  = prhs[0];
    const mxArray* values_  = prhs[1];   // struct, one field per interp var
    const mxArray* order_   = prhs[2];   // int32 [1 x xDim]
    const mxArray* extrap_  = prhs[3];   // int32 [1 x xDim] or empty
    const mxArray* nthr_    = prhs[4];   // int32 scalar

    int xDim    = (int) mxGetNumberOfElements(breaks_);
    int numVec  = (int) mxGetNumberOfFields(values_);
    int* s_order = (int*) mxGetData(order_);
    int  numThreads = *((int*) mxGetData(nthr_));
    bool hasExtrap = (mxGetNumberOfElements(extrap_) > 0);
    int* extrap_order = hasExtrap ? (int*) mxGetData(extrap_) : NULL;
#ifdef USE_OMP
    omp_set_num_threads(numThreads > 0 ? numThreads : 1);
#endif

    // --- grids + (cubic) extrap break-extension. Mirror myppual.m:626-632
    //     (extend each cubic break by one knot at each end, pieces += 2) and
    //     fullConstructionWithExtrap:313-321 (construction consumes the interior
    //     grid: xGridConstruct = extendedGrid+1, xPtsConstruct = extendedPts-2).
    double* xGridExt[MAXDIM];     // breaks carried by the pp structs (extended if cubic+extrap)
    int     xPtsExt[MAXDIM];      // numel of each xGridExt
    double* xGridCon[MAXDIM];     // grid passed to mkl_start[_with_extrap]
    int     xPtsCon[MAXDIM];
    int     pieces[MAXDIM];       // post-extension pieces per dim (xPtsExt-1)
    double* extBuf[MAXDIM];       // owned extended-knot buffers (freed at end)
    for (int d = 0; d < xDim; ++d) {
        const mxArray* g = mxGetCell(breaks_, d);
        double* gp = mxGetPr(g);
        int     gn = (int) mxGetNumberOfElements(g);
        if (hasExtrap && s_order[d] == 4) {
            extBuf[d] = (double*) mxMalloc(sizeof(double) * (gn + 2));
            extBuf[d][0] = gp[0] - 1.0;
            memcpy(extBuf[d] + 1, gp, sizeof(double) * gn);
            extBuf[d][gn + 1] = gp[gn - 1] + 1.0;
            xGridExt[d] = extBuf[d];  xPtsExt[d] = gn + 2;
            xGridCon[d] = gp;         xPtsCon[d] = gn;          // interior == original grid
        } else {
            extBuf[d] = NULL;
            xGridExt[d] = gp;         xPtsExt[d] = gn;
            xGridCon[d] = gp;         xPtsCon[d] = gn;
        }
        pieces[d] = xPtsExt[d] - 1;
    }

    plhs[0] = mxCreateCellMatrix(1, numVec);   // ppCell
    double** natCoefsPerVar = (double**) mxMalloc(sizeof(double*) * numVec);  // for Task 3
    int numArray = 0;
    for (int iv = 0; iv < numVec; ++iv) {
        const mxArray* field = mxGetFieldByNumber(values_, 0, iv);  // [numArray, g0pts, ..., g(D-1)pts]
        numArray = (int) mxGetM(field);                            // array dim is first

        // Values_MKL = permute(values, [xDim+1:-1:1]) as an mxArray (kept alive as
        // the pp's "Values" field). Natural col-major dims (fastest first) are
        // [numArray, g0, ..., g(D-1)]; reversed result dims are [g(D-1),...,g0,numArray].
        // (mkl_start uses the ORIGINAL grid points for construction → valuesMKL is on
        // xPtsCon, not extended.)
        double* src = mxGetPr(field);
        int natDims[MAXDIM + 1];           // fastest first
        natDims[0] = numArray;
        for (int d = 0; d < xDim; ++d) natDims[d + 1] = xPtsCon[d];
        mwSize revDims[MAXDIM + 1];         // reversed, fastest first
        for (int d = 0; d < xDim; ++d) revDims[d] = xPtsCon[xDim - 1 - d];
        revDims[xDim] = numArray;
        mxArray* valuesMKLArr = mxCreateNumericArray(xDim + 1, revDims, mxDOUBLE_CLASS, mxREAL);
        double* valuesMKL = mxGetPr(valuesMKLArr);
        int N = numArray; for (int d = 0; d < xDim; ++d) N *= xPtsCon[d];
        for (int s = 0; s < N; ++s) {
            int t = s, idx[MAXDIM + 1];
            for (int k = 0; k <= xDim; ++k) { idx[k] = t % natDims[k]; t /= natDims[k]; }
            int dst = 0, stride = 1;        // dest fastest = idx[xDim] (g(D-1)); numArray slowest
            for (int k = xDim; k >= 1; --k) { dst += idx[k] * stride; stride *= natDims[k]; }
            dst += idx[0] * stride;
            valuesMKL[dst] = src[s];
        }

        // coefs size (mirror fullConstruction:401-410 / extrap variant:340-356)
        mwSize coeffns[MAXDIM + 1];
        coeffns[xDim] = numArray;          // fortran order: array dim last
        int coeffn = numArray;
        for (int d = 0; d < xDim; ++d) {
            coeffns[xDim - 1 - d] = (mwSize) (pieces[d] * s_order[d]);   // pieces already post-extension
            coeffn *= (int) coeffns[xDim - 1 - d];
        }
        mxArray* coefsArr = mxCreateNumericArray(xDim + 1, coeffns, mxDOUBLE_CLASS, mxREAL);
        double* coeff = mxGetPr(coefsArr);
        constructOne(xDim, xPtsCon, xGridCon, numArray, valuesMKL, s_order,
                     hasExtrap ? extrap_order : NULL, coeff, coeffn);
        natCoefsPerVar[iv] = coeff;        // alias into the pp's coefs (lives in plhs[0])

        // Build pp struct matching myppual(pp)'s MKLpp field set (myppual.m:636-638)
        // — IterRslt.pp is a frozen-shape contract, so carry ALL fields.
        const char* fns[] = {"form","breaks","Values","coefs","pieces","order",
                             "extrap_order","dim","Method","orient","thread"};
        mxArray* pp = mxCreateStructMatrix(1, 1, 11, fns);
        mxSetField(pp, 0, "form",   mxCreateString("MKLpp"));
        mxArray* bcell = mxCreateCellMatrix(1, xDim);
        for (int d = 0; d < xDim; ++d) {
            mxArray* gd = mxCreateDoubleMatrix(1, xPtsExt[d], mxREAL);
            memcpy(mxGetPr(gd), xGridExt[d], sizeof(double) * xPtsExt[d]);
            mxSetCell(bcell, d, gd);
        }
        mxSetField(pp, 0, "breaks", bcell);
        mxSetField(pp, 0, "Values", valuesMKLArr);
        mxSetField(pp, 0, "coefs",  coefsArr);
        mxArray* piecesArr = mxCreateNumericMatrix(1, xDim, mxINT32_CLASS, mxREAL);
        mxArray* orderArr  = mxCreateNumericMatrix(1, xDim, mxINT32_CLASS, mxREAL);
        int* pp_pieces = (int*) mxGetData(piecesArr);
        int* pp_order  = (int*) mxGetData(orderArr);
        for (int d = 0; d < xDim; ++d) { pp_pieces[d] = pieces[d]; pp_order[d] = s_order[d]; }
        mxSetField(pp, 0, "pieces", piecesArr);
        mxSetField(pp, 0, "order",  orderArr);
        if (hasExtrap) {
            mxArray* exArr = mxCreateNumericMatrix(1, xDim, mxINT32_CLASS, mxREAL);
            memcpy(mxGetData(exArr), extrap_order, sizeof(int) * xDim);
            mxSetField(pp, 0, "extrap_order", exArr);
        } else {
            mxSetField(pp, 0, "extrap_order", mxCreateDoubleMatrix(0, 0, mxREAL));
        }
        mxArray* dimArr = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        *((int*) mxGetData(dimArr)) = numArray;
        mxSetField(pp, 0, "dim", dimArr);
        mxSetField(pp, 0, "Method", mxCreateString("not-a-knot"));
        mxSetField(pp, 0, "orient", mxCreateString("MKLC"));
        mxArray* thrArr = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        *((int*) mxGetData(thrArr)) = numThreads;
        mxSetField(pp, 0, "thread", thrArr);
        mxSetCell(plhs[0], iv, pp);
    }
    // Task 3 reads numArray / pieces / natCoefsPerVar here to build plhs[1].
    for (int d = 0; d < xDim; ++d) if (extBuf[d]) mxFree(extBuf[d]);
    mxFree(natCoefsPerVar);
}
```

> The driver test (Step 1) is the bit-exact oracle: if the reverse-permute, the extrap break-extension, or `coeffns` is off, `isequal(ppCell{1}.coefs(:), ref.coefs(:))` fails on the 1-D/2-D/extrap cases. Iterate against it; do not relax to a tolerance.
>
> The 11-field set mirrors `myppual`'s `MKLpp` output so `IterRslt.pp` keeps its frozen shape. The solver's named-scalar interp (`MatlabPp`) reads only `breaks`/`coefs`/`order`/`pieces`/`dim`; the rest exist for `IterRslt.pp` parity. If `myppual`'s field *order* or exact types differ from this list, match `myppual.m:636-638` (the Task 6 `IterRslt.pp` `fieldnames` parity check pins it).

- [ ] **Step 4: Build and iterate against the test until green**

Run: `matlab -batch "cd('tests'); gdsge.codegen.ensureSplineConstructMex; r = runtests('tInterpConstructMex'); assert(all([r.Passed]))"`
Expected: PASS — all three `ppCellMatchesMyppual_*` cases bit-identical. Fix orientation/extrap until `isequal` holds. (`tEnsureSplineConstructMex` now also passes — run it too.)

- [ ] **Step 5: Commit**

```bash
git add src/kernels/interp_construct_mex.cpp src/kernels/interp_construct_mex.mexw64 src/kernels/interp_construct_mex.cache tests/kernels/tInterpConstructMex.m
git commit -m "feat(kernels): interp_construct_mex — per-var construction, bit-exact vs myppual"
```

---

## Task 3: Eval-order packing — `GDSGE_SPLINE_VEC` second output

Add the second output: coefficients written directly in the eval order that `convert_to_interp_eval_array` produces, plus the scalar/grid metadata fields.

**Files:**
- Modify: `src/kernels/interp_construct_mex.cpp`
- Modify: `tests/kernels/tInterpConstructMex.m`

- [ ] **Step 1: Add the failing parity test (splineVec vs convert_to_interp_eval_array)**

Append these tests to `tInterpConstructMex.m`:

```matlab
        function splineVecMatchesConvert_helper(tc, breaks, order, extrap, numArray, numVars)
            sizeState = num2cell(cellfun(@numel, breaks));
            ppCellRef = cell(1, numVars); S = struct();
            for k = 1:numVars
                v = tInterpConstructMex.build(breaks, numArray) + 10*k;
                ppCellRef{k} = tInterpConstructMex.mkpp(v, breaks, order, extrap);
                S.(sprintf('v%d', k)) = reshape(v, [], sizeState{:});
            end
            refVec = convert_to_interp_eval_array(ppCellRef);
            [~, vec] = interp_construct_mex(breaks, S, int32(order*ones(1,numel(breaks))), ...
                int32(extrap), int32(1));
            tc.verifyTrue(isequal(vec.coefs(:), refVec.coefs(:)), 'coefs not bit-identical');
            for f = ["fullVecEvalCoefsLength","singleVecEvalCoefsLength","xDim","dim","arrayOffset"]
                tc.verifyEqual(double(vec.(f)), double(refVec.(f)), sprintf('%s differs', f));
            end
            tc.verifyTrue(isequal(int32(vec.order(:)),  int32(refVec.order(:))));
            tc.verifyTrue(isequal(int32(vec.pieces(:)), int32(refVec.pieces(:))));
        end
        function splineVec_1D(tc)
            tc.splineVecMatchesConvert_helper({linspace(0,1,7)}, 4, 4, 3, 2);
        end
        function splineVec_2D_extrap(tc)
            tc.splineVecMatchesConvert_helper({linspace(0,1,6),linspace(-1,2,5)}, 4, 2, 4, 3);
        end
        function splineVec_3D(tc)
            tc.splineVecMatchesConvert_helper({linspace(0,1,5),linspace(0,1,4),linspace(0,1,6)}, 4, 4, 2, 2);
        end
        function splineVec_linear(tc)
            tc.splineVecMatchesConvert_helper({linspace(0,1,8)}, 2, [], 2, 2);
        end
```

(Change the helper to a `methods (Test, ...)`-callable instance method, or inline it — keep it a `methods` helper invoked by the four test methods.)

- [ ] **Step 2: Run to verify it fails**

Run: `matlab -batch "cd('tests'); r = runtests('tInterpConstructMex'); assert(all([r.Passed]))"`
Expected: FAIL — `interp_construct_mex` returns only 1 output (`plhs[1]` unassigned) → "One or more output arguments not assigned".

- [ ] **Step 3: Implement the eval-order pack in the MEX**

Add, after the per-var construction loop, code that writes `plhs[1]`. The eval-order layout (derived from `convert_to_interp_eval_array.m`) is, fastest dimension first:

```
[ i_vec , order_1..order_D , pieces_1..pieces_D , i_array ]
```

so `evalIdx = i_vec + numVec * base`, where `base` is the mixed-radix index of `(o_1..o_D, p_1..p_D, a)` with radices `[order_1..order_D, pieces_1..pieces_D, numArray]`. Each per-var natural coefs buffer (length `L = prod(order)*prod(pieces)*numArray`) is laid out (fastest first) as `[order_D, pieces_D, ..., order_1, pieces_1, numArray]`. Build the natural→base index map once, then scatter every variable:

```cpp
// metadata (== convert_to_interp_eval_array)
int prodOrder = 1, prodPieces = 1;
for (int d = 0; d < xDim; ++d) { prodOrder *= s_order[d]; prodPieces *= pieces[d]; }
int fullVecEvalCoefsLength   = numVec * prodOrder;
int singleVecEvalCoefsLength = prodOrder;                       // == fullVec/numVec (convert divides by dim)
int arrayOffset              = numVec * prodOrder * prodPieces;
int L                        = prodOrder * prodPieces * numArray;   // per-var coefs length

// natIdx -> base map (independent of i_vec)
int* baseMap = (int*) mxMalloc(sizeof(int) * L);
for (int natIdx = 0; natIdx < L; ++natIdx) {
    int t = natIdx;
    int o[MAXDIM], p[MAXDIM], a;
    // decompose natIdx, fastest first: order_D, pieces_D, ..., order_1, pieces_1, numArray
    for (int d = xDim - 1; d >= 0; --d) {
        o[d] = t % s_order[d]; t /= s_order[d];
        p[d] = t % pieces[d];  t /= pieces[d];
    }
    a = t;   // remaining is numArray index
    // recompose base, fastest first: order_1..order_D, then pieces_1..pieces_D, then numArray
    int base = 0, stride = 1;
    for (int d = 0; d < xDim; ++d) { base += o[d] * stride; stride *= s_order[d]; }
    for (int d = 0; d < xDim; ++d) { base += p[d] * stride; stride *= pieces[d]; }
    base += a * stride;
    baseMap[natIdx] = base;
}

mwSize evalDims[2] = { (mwSize)(arrayOffset * numArray), 1 };   // flat column; C++ reads mxGetPr flat, tests compare (:)
mxArray* coefsOut = mxCreateNumericArray(2, evalDims, mxDOUBLE_CLASS, mxREAL);
double* evalCoefs = mxGetPr(coefsOut);
for (int iv = 0; iv < numVec; ++iv) {
    double* nat = natCoefsPerVar[iv];   // the per-var natural coefs from Task 2
    for (int natIdx = 0; natIdx < L; ++natIdx)
        evalCoefs[iv + numVec * baseMap[natIdx]] = nat[natIdx];
}
mxFree(baseMap);
```

Then assemble `plhs[1]` as a struct with fields exactly matching `convert_to_interp_eval_array`'s output **and the types `MatlabInterpEval.h` expects**:
- **`double` scalars** (read via `mxGetPr`): `fullVecEvalCoefsLength`, `singleVecEvalCoefsLength` (= `prodOrder`, **not** `fullVecEvalCoefsLength`), `xDim`, `dim` (= `numVec`), `arrayOffset`.
- **`int32` rows** (read via `mxGetData` with `memcpy ... sizeof(int)*xDim` — **must be int32, not double**): `order` (= `s_order`), `pieces` (post-extension), `xPts` (= `pieces + 1`).
- `coefs` (= `coefsOut`, double), `breaks` (the input cell of double grids — extended under extrap, matching the pp structs).

(`convert_to_interp_eval_array` gets these types for free because it copies `pp.order`/`pp.pieces`, which `myppual` already stored as `int32`; reproduce them as `int32` explicitly here.)

> Note on `pieces` under extrap: `convert_to_interp_eval_array` reads `pp.pieces`, which for cubic+extrap already includes the `+2` per dim. Use the **post-extension** `pieces` here (same value the Task 2 pp structs carry), so coefs length and `arrayOffset` match.

- [ ] **Step 4: Run to verify green**

Run: `matlab -batch "cd('tests'); gdsge.codegen.ensureSplineConstructMex; r = runtests('tInterpConstructMex'); assert(all([r.Passed]))"`
Expected: PASS — 1D/2D/3D/linear, with/without extrap, multi-var, `numArray>1` all bit-identical to `convert_to_interp_eval_array`.

- [ ] **Step 5: Commit**

```bash
git add src/kernels/interp_construct_mex.cpp src/kernels/interp_construct_mex.mexw64 src/kernels/interp_construct_mex.cache tests/kernels/tInterpConstructMex.m
git commit -m "feat(kernels): interp_construct_mex eval-order pack — bit-exact vs convert_to_interp_eval_array"
```

---

## Task 4: Wire `constructSplines` to the fused MEX (with internal legacy switch)

**Files:**
- Modify: `src/+gdsge/+runtime/constructSplines.m`
- Create: `src/+gdsge/+runtime/useFusedConstruct.m` (the internal switch)
- Test: `tests/runtime/tConstructSplinesEquivalence.m`

- [ ] **Step 1: Write the failing equivalence test (fused vs legacy, bit-identical)**

```matlab
classdef tConstructSplinesEquivalence < matlab.unittest.TestCase
    methods (TestClassSetup)
        function setup(tc)
            gdsge.codegen.ensureSplineConstructMex();
            root = fileparts(which('essential_blas.dll'));
            if ~isempty(root) && ~any(strcmp(strsplit(getenv('PATH'),';'),root))
                setenv('PATH', [getenv('PATH') ';' root]);
            end
            tc.addTeardown(@() gdsge.runtime.useFusedConstruct(true));  % restore default
        end
    end
    methods (Test)
        function fusedEqualsLegacy_2D(tc)
            breaks = {linspace(0,1,6), linspace(-1,2,5)};
            sizeState = num2cell(cellfun(@numel, breaks));
            valueCell = {reshape(1:4*30,4,30)+0.5, reshape(2:4*30+1,4,30)-0.25};
            gdsge.runtime.useFusedConstruct(false);
            [ppL, vecL] = gdsge.runtime.constructSplines(valueCell, breaks, sizeState, 4, 2, 1);
            gdsge.runtime.useFusedConstruct(true);
            [ppF, vecF] = gdsge.runtime.constructSplines(valueCell, breaks, sizeState, 4, 2, 1);
            tc.verifyTrue(isequal(vecL.coefs(:), vecF.coefs(:)), 'splineVec.coefs diverged');
            for k = 1:numel(ppL)
                tc.verifyTrue(isequal(ppL{k}.coefs(:), ppF{k}.coefs(:)), 'ppCell coefs diverged');
            end
        end
    end
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `matlab -batch "cd('tests'); r = runtests('tConstructSplinesEquivalence'); assert(all([r.Passed]))"`
Expected: FAIL — `gdsge.runtime.useFusedConstruct` undefined.

- [ ] **Step 3: Add the switch helper**

```matlab
function out = useFusedConstruct(set)
% USEFUSEDCONSTRUCT  Internal toggle for gdsge.runtime.constructSplines:
%   true (default) = fused interp_construct_mex; false = legacy myppual +
%   convert_to_interp_eval_array. Production stays true; the A/B gates flip it.
persistent val
if isempty(val); val = true; end
if nargin > 0; val = logical(set); end
out = val;
end
```

- [ ] **Step 4: Rewrite `constructSplines.m`**

```matlab
function [ppCell, splineVec] = constructSplines(valueCell, breaks, sizeState, interpOrder, extrapOrder, numThreads)
% CONSTRUCTSPLINES  Build interpolants for each interp variable plus the
%   vectorized evaluation array consumed by the MEX (GDSGE_SPLINE_VEC).
%   Default path: one fused interp_construct_mex call (constructs all vars and
%   writes GDSGE_SPLINE_VEC directly in eval order). Legacy path (per-var
%   myppual + convert_to_interp_eval_array) stays reachable via
%   gdsge.runtime.useFusedConstruct(false) for the A/B differential gates.
orderVec = interpOrder * ones(1, numel(sizeState));
if interpOrder == 4
    extrapVec = extrapOrder * ones(1, numel(sizeState));
else
    extrapVec = [];
end
nv = numel(valueCell);
if nv == 0
    ppCell = {}; splineVec = struct(); return;
end
if gdsge.runtime.useFusedConstruct()
    values = struct();
    for i = 1:nv
        values.(sprintf('GDSGE_V%d', i)) = reshape(valueCell{i}, [], sizeState{:});
    end
    [ppCell, splineVec] = interp_construct_mex(breaks, values, ...
        int32(orderVec), int32(extrapVec), int32(numThreads));
else
    ppCell = cell(1, nv);
    for i = 1:nv
        pp = struct('form','MKL','breaks',{breaks}, ...
            'Values',reshape(valueCell{i},[],sizeState{:}),'coefs',[],'order',orderVec, ...
            'Method',[],'ExtrapolationOrder',extrapVec,'thread',numThreads, ...
            'orient','curvefit');
        ppCell{i} = myppual(pp);
    end
    splineVec = convert_to_interp_eval_array(ppCell);
end
end
```

- [ ] **Step 5: Run to verify green**

Run: `matlab -batch "cd('tests'); r = runtests('tConstructSplinesEquivalence'); assert(all([r.Passed]))"`
Expected: PASS — fused and legacy bit-identical.

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+runtime/constructSplines.m src/+gdsge/+runtime/useFusedConstruct.m tests/runtime/tConstructSplinesEquivalence.m
git commit -m "feat(runtime): constructSplines uses fused interp_construct_mex (legacy behind useFusedConstruct switch)"
```

---

## Task 5: Invoke the compile gate at codegen entry

So a fresh checkout (or a changed kernel) builds `interp_construct_mex` before `iter_<model>` runs, exactly as `ensureAsgMex`.

**Files:**
- Modify: `src/+gdsge/+codegen/generateCxx.m:16`
- Test: `tests/codegen/tEnsureSplineConstructMex.m` (extend)

- [ ] **Step 1: Add the failing assertion**

Append to `tEnsureSplineConstructMex.m` a check that `generateCxx`'s source text invokes the gate (cheap, non-brittle — the functional end-to-end coverage is Tasks 6/7):

```matlab
        function generateCxxInvokesEnsure(tc)
            src = fileread(which('gdsge.codegen.generateCxx'));
            tc.verifyTrue(contains(src, 'ensureSplineConstructMex'), ...
                'generateCxx must call ensureSplineConstructMex');
            tc.verifyEqual(exist('interp_construct_mex','file'), 3);
        end
```

- [ ] **Step 2: Wire it in `generateCxx.m`** (beside the existing `ensureAsgMex`)

At `src/+gdsge/+codegen/generateCxx.m:16`, after:
```matlab
    gdsge.codegen.ensureAsgMex();   % emitCompile reads asg.get_mex_constants()
```
add:
```matlab
    gdsge.codegen.ensureSplineConstructMex();  % fused spline constructor for constructSplines
```

- [ ] **Step 3: Run to verify green**

Run: `matlab -batch "cd('tests'); r = runtests('tEnsureSplineConstructMex'); assert(all([r.Passed]))"`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/+gdsge/+codegen/generateCxx.m tests/codegen/tEnsureSplineConstructMex.m
git commit -m "feat(codegen): generateCxx ensures interp_construct_mex (parity with ensureAsgMex)"
```

---

## Task 6: A/B end-to-end bit-exact gates (HL1996 + safe_assets)

Prove the fused constructor yields a **bit-identical** `IterRslt` vs the legacy path, through the real generated `iter_<model>`, on a pure-vectorized model (HL1996) and a named-scalar-interp + multi-root model (safe_assets).

**Files:**
- Create: `tests/HeatonLucas1996/codegen/tFusedConstructHL1996.m`
- Create: `tests/Barro_et_al_2017/codegen/tFusedConstructSafeAssets.m`

- [ ] **Step 1: Write the HL1996 A/B gate** (mirrors `tInMexResolveHL1996.m`; toggles the constructor, pins `rng`)

```matlab
classdef tFusedConstructHL1996 < matlab.unittest.TestCase
    % A/B differential: fused interp_construct_mex vs legacy myppual+convert
    % must produce a bit-identical IterRslt.
    methods (Test, TestTags = {'Slow'})
        function fusedEqualsLegacy(tc)
            here = fileparts(mfilename('fullpath'));
            modelDir = fileparts(here);
            work = tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture).Folder;
            copyfile(fullfile(modelDir, 'HL1996.gmod'), work);
            tc.addTeardown(@() gdsge.runtime.useFusedConstruct(true));

            gdsge_codegen('HL1996');
            opts = struct('SaveFreq', inf, 'NoSave', 1);

            gdsge.runtime.useFusedConstruct(true);  rng(0823); Rf = iter_HL1996(opts);
            gdsge.runtime.useFusedConstruct(false); rng(0823); Rl = iter_HL1996(opts);

            tc.verifyEqual(Rf.Iter, Rl.Iter, 'iteration count diverged');
            tc.verifyEqual(Rf.Metric, Rl.Metric, 'final metric diverged');
            flds = {'var_policy','var_aux','var_interp'};
            for i = 1:numel(flds)
                a = Rf.(flds{i}); b = Rl.(flds{i}); fn = fieldnames(a);
                for k = 1:numel(fn)
                    tc.verifyTrue(isequal(a.(fn{k}), b.(fn{k})), ...
                        sprintf('%s.%s not bit-identical (fused vs legacy)', flds{i}, fn{k}));
                end
            end
            % IterRslt.pp is a frozen-shape contract: the fused pp struct must
            % carry the same field set as the legacy myppual MKLpp.
            ppn = fieldnames(Rf.pp);
            for k = 1:numel(ppn)
                tc.verifyEqual(fieldnames(Rf.pp.(ppn{k})), fieldnames(Rl.pp.(ppn{k})), ...
                    sprintf('IterRslt.pp.%s field set changed (frozen shape)', ppn{k}));
            end
        end
    end
end
```

- [ ] **Step 2: Write the safe_assets A/B gate** (same shape; model = `safe_assets`, `iter_safe_assets`, gmod under `tests/Barro_et_al_2017`)

Copy the HL1996 gate, replacing the model name/gmod/`iter_*` with `safe_assets` (confirm the gmod filename in `tests/Barro_et_al_2017` — it is the file the existing `tInMexResolveSafeAssets.m` copies). Keep `rng(0823)` before each run (safe_assets' multi-root restart is RNG-fragile; the in-MEX-resolve gate pins the same).

- [ ] **Step 3: Run both gates**

Run: `matlab -batch "cd('tests'); r = runtests({'tFusedConstructHL1996','tFusedConstructSafeAssets'}); disp(table(r)); assert(all([r.Passed]))"`
Expected: PASS — both `Iter`/`Metric`/all fields bit-identical.

- [ ] **Step 4: Commit**

```bash
git add tests/HeatonLucas1996/codegen/tFusedConstructHL1996.m tests/Barro_et_al_2017/codegen/tFusedConstructSafeAssets.m
git commit -m "test: A/B gates — fused constructor == legacy, bit-identical (HL1996 + safe_assets)"
```

---

## Task 7: Full-suite regression + perf bench + docs + changelog

**Files:**
- Modify: `docs/architecture.md`, `PROGRESS.md`
- Create: `tests/perf/fused_construct_bench.m`

- [ ] **Step 1: Run the full test suite (all corpus models, both backends)**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: exit 0; `tests/results/junit.xml` all green (prior count + the new kernel/runtime/codegen/A-B tests). Every existing end-to-end golden (HL1996, safe_assets, Mendoza2010, GLSW, Cao2011EZ, CaoNie2016; autodiff + SymPy) stays green — the constructor is internal and bit-exact.

- [ ] **Step 2: Perf micro-bench (rng-pinned; construct-time + iter-time, fused vs legacy)**

Create `tests/perf/fused_construct_bench.m` modeled on `tests/perf/inmex_resolve_bench.m`: for GLSW and safe_assets, run `iter_<model>` under `useFusedConstruct(true)` and `(false)` with `rng` pinned, assert identical `Iter`/`Metric`, and print wall-clock for each. (One MATLAB process per source per the golden-rule; reuse the `tests/perf` harness.)

Run: `matlab -batch "cd('tests/perf'); fused_construct_bench"`
Expected: identical Iter/Metric; fused ≤ legacy wall-clock (expected larger win on GLSW per the Phase-10 perf report's ~30% spline-rebuild share).

- [ ] **Step 3: Update `docs/architecture.md`**

Add `interp_construct_mex` to the kernel inventory: minimal, self-contained fused spline/linear constructor; construction-only (legacy `myppual_mex` retained, prebuilt, for the simulate/`output_interp` eval path); cache-gated compile via `ensureSplineConstructMex` (flags from `compile_myppual.m` + `/O2`, `/fp:precise`); consumed by `gdsge.runtime.constructSplines`; bit-exact with the legacy path.

- [ ] **Step 4: Add the `PROGRESS.md` changelog entry**

Prepend under `## Changelog` a `2026-06-14: **Fused in-MEX spline construction**` entry: what changed (one fused construct MEX replaces per-var `myppual` + `convert_to_interp_eval_array` permute/concat on the cartesian spline/linear iteration path; construction-only minimal kernel, prebuilt `myppual_mex` kept for simulate eval); bit-exact (driver test vs `myppual`/`convert`; A/B end-to-end gates vs legacy, rng-pinned); perf result from Step 2; scope boundary (ASG/pchip/simulate unchanged, solver MEX untouched); branch `fused-spline-construct`; spec/plan paths.

- [ ] **Step 5: Commit**

```bash
git add docs/architecture.md PROGRESS.md tests/perf/fused_construct_bench.m
git commit -m "docs(perf): fused spline construct — perf bench, architecture, PROGRESS changelog"
```

- [ ] **Step 6: Final verification**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: exit 0, full suite green. Report the count from `tests/results/junit.xml`.

---

## Self-review notes (coverage map)

- Spec §2 goal 1 (one fused call, struct input, eval-order direct, both forms) → Tasks 2–4.
- Spec §2 goal 2 (bit-exact, no new goldens) → Tasks 2/3 driver tests + Task 6 A/B gates + Task 7 full suite.
- Spec §2 goal 3 (solver MEX untouched) → constructSplines returns the same two contracts; no template edits (verified by Task 6/7 goldens).
- Spec §2 goal 4 (minimal self-contained kernel; prebuilt myppual_mex kept) → Task 2 file; Task 7 docs.
- Spec §6 eval-order layout → Task 3 (normative parity vs `convert_to_interp_eval_array`).
- Spec §7 MATLAB integration (constructSplines, ensure gate, codegen entry) → Tasks 4/5.
- Spec §9 edge cases (extrap, linear, empty list) → Task 2/3 tests (extrap, linear) + `nv==0` guard in Task 4 constructSplines.
- Spec §10 testing → Tasks 2/3/6/7.
- Spec §8 compile flags (`/fp:precise`, `/openmp`, `-DUSE_OMP`) → Task 1.
