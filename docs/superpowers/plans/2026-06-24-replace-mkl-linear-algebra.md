# Replace MKL Linear Solver with Eigen — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Profile MKL vs Eigen dense linear-algebra variants in the CoDoSol nonlinear solver, then replace the MKL path (`dgetrf`/`dgetrs`/`dgemv` + `essential_blas.lib`) with the winning Eigen variant so the toolbox has one dependency-free dense path.

**Architecture:** A self-contained micro-benchmark MEX (`bench_linalg`) ranks five solver variants on synthetic 23×23 systems in a single process (apples-to-apples). The current MKL build supplies an end-to-end baseline. The winner is implemented in `include/codosol.h` by collapsing the `#ifdef USE_MKL` branches into one Eigen path, the build template drops `-DUSE_MKL` and the `essential_blas.lib` link, goldens are regenerated, and the full test suite plus an end-to-end re-time confirm correctness and the ≤5% rule.

**Tech Stack:** C++ (MEX), Eigen (vendored under `include/Eigen`), Intel MKL via `include/essential_blas.lib` (being removed), MATLAB R2025b test harness.

## Global Constraints

- **MATLAB:** R2025b at `C:\Program Files\MATLAB\R2025b\bin\matlab.exe`; drive headless with `matlab -batch "<expr>"`. Exit 0 = success.
- **Run tests with:** `matlab -batch "cd('tests'); run_tests"`. Authoritative report: `tests/results/junit.xml` + exit code (NOT `results.tap` — it accumulates stale lines).
- **One MATLAB process at a time** — each saturates all cores via OpenMP; concurrent runs distort timing. Run sequentially.
- **PowerShell 7 (`pwsh`) is NOT installed** — do not use `tests/run.ps1`; use `matlab -batch`.
- **Golden MATLAB-path rule:** never add a toolbox source to a persistent path; each run `addpath`s exactly one source in its own `matlab -batch` process.
- **Backward compatible:** every existing `.gmod` must still run; `IterRslt`/`SimuRslt` shapes frozen. Default path stays MATLAB + C MEX, no Python.
- **Decision rule:** adopt Eigen if the winner is within **≤5% of MKL on total HL1996 end-to-end iter time** AND golden outputs match. Otherwise keep MKL and report.
- **Scratch is git-ignored:** throwaway benchmark sources live under `scratch/bench_mkl/`. Committed evidence goes into the spec/results doc.
- HL1996: `NUM_EQUATIONS = MAXDIM = 23`. `jac` is a contiguous **n×n column-major buffer, leading dim = n** (NOT padded to MAXDIM) — so a compile-time `Matrix<double,MAXDIM,MAXDIM>` map is correct only when `n == MAXDIM`.

---

### Task 1: Build the micro-benchmark MEX and record per-solve timings

**Files:**
- Create: `scratch/bench_mkl/bench_linalg.cpp`
- Create: `scratch/bench_mkl/build_and_run.m`

**Interfaces:**
- Produces: a MEX `bench_linalg(n, nMat, nReps)` returning a struct with fields `v0_mkl_lu, v1_eigen_qr, v2_eigen_lu, v3_eigen_lu_fixed, v4_eigen_lu_nowrite, gemv_mkl, gemv_eigen` (ns per operation), consumed by the decision in Task 3.

**Note on synthetic data:** LU/QR solve *timing* is essentially value-independent (same flop count regardless of matrix entries; pivoting branches are a negligible fraction). So synthetic well-conditioned 23×23 matrices match real Newton jacobians for *timing*. Real-jacobian *correctness* is covered by the golden suite in Task 6, not here.

- [ ] **Step 1: Write the benchmark source**

Create `scratch/bench_mkl/bench_linalg.cpp`:

```cpp
// bench_linalg.cpp — micro-benchmark of CoDoSol's per-solve linear algebra.
// Compares the MKL path (dgetrf/dgetrs, dgemv) against Eigen variants on
// synthetic n x n systems, all within one process for apples-to-apples timing.
#include "mex.h"
#include <chrono>
#include <vector>
#include <random>
#include <cstring>
#include <Eigen/Dense>
#include "essential_blas.h"

using clk = std::chrono::high_resolution_clock;
using Eigen::MatrixXd; using Eigen::VectorXd; using Eigen::Map;

#ifndef MAXDIM
#define MAXDIM 23
#endif

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    int n     = (nrhs>0)? (int)mxGetScalar(prhs[0]) : MAXDIM;
    int nMat  = (nrhs>1)? (int)mxGetScalar(prhs[1]) : 256;
    int nReps = (nrhs>2)? (int)mxGetScalar(prhs[2]) : 4000;

    // batch of well-conditioned matrices + rhs (fixed seed -> reproducible)
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    std::vector<std::vector<double>> A(nMat, std::vector<double>(n*n));
    std::vector<std::vector<double>> b(nMat, std::vector<double>(n));
    for (int k=0;k<nMat;k++){
        for (int j=0;j<n;j++){
            for (int i=0;i<n;i++) A[k][i+j*n] = U(rng);
            A[k][j+j*n] += n;            // diagonal dominance -> well conditioned
            b[k][j] = U(rng);
        }
    }

    double checksum = 0;   // defeats dead-code elimination
    auto bench = [&](auto solveOne)->double{
        auto t0 = clk::now();
        for (int r=0;r<nReps;r++)
            for (int k=0;k<nMat;k++)
                checksum += solveOne(A[k].data(), b[k].data());
        auto t1 = clk::now();
        return std::chrono::duration<double,std::nano>(t1-t0).count()
               / (double)((long)nReps*nMat);   // ns per solve
    };

    std::vector<double> tmp(n*n), sol(n);
    std::vector<int> ipiv(n);

    double v0 = bench([&](const double* a, const double* rhs){     // MKL LU
        std::memcpy(tmp.data(), a, sizeof(double)*n*n);
        std::memcpy(sol.data(), rhs, sizeof(double)*n);
        _dgetrf(n,n,tmp.data(),n,ipiv.data());
        _dgetrs('N',n,1,tmp.data(),n,ipiv.data(),sol.data(),n);
        return sol[0];
    });
    double v1 = bench([&](const double* a, const double* rhs){     // Eigen colPiv QR
        Map<const MatrixXd> J(a,n,n); Map<const VectorXd> f(rhs,n);
        VectorXd s = J.colPivHouseholderQr().solve(f);
        return s[0];
    });
    double v2 = bench([&](const double* a, const double* rhs){     // Eigen partialPiv LU
        Map<const MatrixXd> J(a,n,n); Map<const VectorXd> f(rhs,n);
        VectorXd s = J.partialPivLu().solve(f);
        return s[0];
    });
    double v3 = -1;                                                // Eigen LU, fixed size
    if (n==MAXDIM){
        typedef Eigen::Matrix<double,MAXDIM,MAXDIM> MatN;
        typedef Eigen::Matrix<double,MAXDIM,1> VecN;
        v3 = bench([&](const double* a, const double* rhs){
            Map<const MatN> J(a); Map<const VecN> f(rhs);
            VecN s = J.partialPivLu().solve(f);
            return s[0];
        });
    }
    double v4 = bench([&](const double* a, const double* rhs){     // LU, write to caller buf
        Map<const MatrixXd> J(a,n,n); Map<const VectorXd> f(rhs,n);
        Map<VectorXd> s(sol.data(),n);
        s = J.partialPivLu().solve(f);
        return s[0];
    });

    auto benchGemv = [&](auto one)->double{
        auto t0=clk::now();
        for(int r=0;r<nReps;r++) for(int k=0;k<nMat;k++) checksum+=one(A[k].data(),b[k].data());
        auto t1=clk::now();
        return std::chrono::duration<double,std::nano>(t1-t0).count()/(double)((long)nReps*nMat);
    };
    double g_mkl = benchGemv([&](const double* a,const double* x){
        std::vector<double> y(n);
        _dgemv('T',n,n,1.0,a,n,x,1,0.0,y.data(),1);
        return y[0];
    });
    double g_eig = benchGemv([&](const double* a,const double* x){
        Map<const MatrixXd> J(a,n,n); Map<const VectorXd> X(x,n);
        VectorXd y = J.transpose()*X;
        return y[0];
    });

    const char* fields[] = {"n","nMat","nReps","v0_mkl_lu","v1_eigen_qr",
        "v2_eigen_lu","v3_eigen_lu_fixed","v4_eigen_lu_nowrite",
        "gemv_mkl","gemv_eigen","checksum"};
    plhs[0] = mxCreateStructMatrix(1,1,11,fields);
    auto setf=[&](const char* f,double v){ mxSetField(plhs[0],0,f,mxCreateDoubleScalar(v)); };
    setf("n",n); setf("nMat",nMat); setf("nReps",nReps);
    setf("v0_mkl_lu",v0); setf("v1_eigen_qr",v1); setf("v2_eigen_lu",v2);
    setf("v3_eigen_lu_fixed",v3); setf("v4_eigen_lu_nowrite",v4);
    setf("gemv_mkl",g_mkl); setf("gemv_eigen",g_eig); setf("checksum",checksum);
}
```

- [ ] **Step 2: Write the build+run driver**

Create `scratch/bench_mkl/build_and_run.m`:

```matlab
function build_and_run()
% Builds bench_linalg with the SAME flags production uses (USE_MKL + essential_blas.lib)
% so V0's MKL linkage matches the real kernel, then runs and prints results.
here = fileparts(mfilename('fullpath'));
inc  = fullfile(here, '..', '..', 'include');
src  = fullfile(here, 'bench_linalg.cpp');
lib  = fullfile(inc, 'essential_blas.lib');
cmd = ['mex -DMAXDIM=23 -DEIGEN_DONT_PARALLELIZE -DUSE_MKL ' ...
       'OPTIMFLAGS="/O2 /DNDEBUG" COMPFLAGS="$COMPFLAGS /openmp" ' ...
       sprintf('"%s" "%s" -I"%s" -outdir "%s"', src, lib, inc, here)];
eval(cmd);
addpath(here);
% warm up, then take the median of several runs
R = bench_linalg(23, 256, 4000);
for k = 1:5, R = bench_linalg(23, 256, 4000); end
fprintf('\n=== bench_linalg  n=%d  nMat=%d  nReps=%d ===\n', R.n, R.nMat, R.nReps);
fprintf('V0 MKL LU            : %8.1f ns/solve\n', R.v0_mkl_lu);
fprintf('V1 Eigen colPiv QR   : %8.1f ns/solve  (%.2fx V0)\n', R.v1_eigen_qr,  R.v1_eigen_qr/R.v0_mkl_lu);
fprintf('V2 Eigen partialPiv  : %8.1f ns/solve  (%.2fx V0)\n', R.v2_eigen_lu,  R.v2_eigen_lu/R.v0_mkl_lu);
fprintf('V3 Eigen LU fixed    : %8.1f ns/solve  (%.2fx V0)\n', R.v3_eigen_lu_fixed, R.v3_eigen_lu_fixed/R.v0_mkl_lu);
fprintf('V4 Eigen LU nowrite  : %8.1f ns/solve  (%.2fx V0)\n', R.v4_eigen_lu_nowrite, R.v4_eigen_lu_nowrite/R.v0_mkl_lu);
fprintf('gemv MKL / Eigen     : %8.1f / %8.1f ns\n', R.gemv_mkl, R.gemv_eigen);
end
```

- [ ] **Step 3: Build and run**

Run: `matlab -batch "cd('scratch/bench_mkl'); build_and_run"`
Expected: MEX compiles without error; a printed table with V0..V4 and gemv timings. V1 (QR) is expected to be the slowest Eigen variant; V2 the closest to V0.

- [ ] **Step 4: Record results**

Paste the printed table into the spec under a new `## Results` section (the spec file `docs/superpowers/specs/2026-06-24-replace-mkl-linear-algebra-design.md`). Do not edit `codosol.h` yet.

- [ ] **Step 5: Commit**

```bash
git add scratch/bench_mkl/bench_linalg.cpp scratch/bench_mkl/build_and_run.m docs/superpowers/specs/2026-06-24-replace-mkl-linear-algebra-design.md
git commit -m "perf(codosol): micro-benchmark MKL vs Eigen dense solve variants"
```
(Note: `scratch/` is git-ignored, so the two `scratch/bench_mkl/*` paths may not stage — if `git add` reports they're ignored, that is expected; commit only the spec change. Force-add with `git add -f` only if you want the bench sources tracked.)

---

### Task 2: Capture the MKL end-to-end baseline for HL1996

This MUST run on the **current code** (MKL path still active) — the baseline disappears once Task 4 deletes the MKL branch.

**Files:**
- Create: `scratch/bench_mkl/time_hl1996.m`

**Interfaces:**
- Produces: a baseline median iter time (seconds) for HL1996, recorded in the spec `## Results`, consumed by the ≤5% decision in Task 7.

- [ ] **Step 1: Write the timing script**

Create `scratch/bench_mkl/time_hl1996.m`:

```matlab
function t = time_hl1996()
% Codegen + compile + time iter_HL1996 to convergence. Mirrors tEndToEndHL1996.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot,'src'));
work = tempname; mkdir(work);
copyfile(fullfile(repoRoot,'tests','HeatonLucas1996','HL1996.gmod'), work);
oldp = pwd; cd(work); cleanup = onCleanup(@() cd(oldp));
gdsge_codegen('HL1996');                 % compiles the MEX
opts = struct('SaveFreq', inf, 'NoSave', 1);
iter_HL1996(opts);                        % warm-up run (discard)
ts = zeros(1,3);
for k = 1:3
    tic; R = iter_HL1996(opts); ts(k) = toc;
end
t = median(ts);
fprintf('\nHL1996 iter: median %.3f s over 3 runs [%.3f %.3f %.3f], Iter=%d Metric=%.2e\n', ...
    t, ts(1), ts(2), ts(3), R.Iter, R.Metric);
end
```

- [ ] **Step 2: Run it on current (MKL) code**

Run: `matlab -batch "cd('scratch/bench_mkl'); time_hl1996"`
Expected: prints a median iter time; `Metric` < 1e-6 and `Iter` ≈ 209 (matching the golden), confirming a real converged run.

- [ ] **Step 3: Record the MKL baseline**

Add the median time to the spec `## Results` as "MKL end-to-end baseline".

- [ ] **Step 4: Commit**

```bash
git add docs/superpowers/specs/2026-06-24-replace-mkl-linear-algebra-design.md
git commit -m "perf(HL1996): record MKL end-to-end iter baseline"
```

---

### Task 3: Decision checkpoint — pick the winning variant

No code change. This is an analysis gate; its deliverable is a one-paragraph decision recorded in the spec.

- [ ] **Step 1: Evaluate the micro-benchmark table**

Confirm the core hypothesis (V1 QR → V2 LU recovers most of the gap). Choose the production variant:
- If V3 (fixed size) is meaningfully faster than V2 **and** within the same ballpark of practicality → production uses V3 **gated on `n == MAXDIM`** with a V2 fallback.
- Otherwise → production uses **V2** (dynamic `partialPivLu`), the simplest faithful match to MKL.
- If V4 ≈ V2 (expected), prefer V2 for clarity.

- [ ] **Step 2: Record the decision**

Write the chosen variant + the one-line rationale into the spec `## Results`. Reference the exact ns/solve numbers.

- [ ] **Step 3: Commit**

```bash
git add docs/superpowers/specs/2026-06-24-replace-mkl-linear-algebra-design.md
git commit -m "docs(spec): record solver-variant decision from micro-benchmark"
```

---

### Task 4: Implement the winning variant in `codosol.h`

**Files:**
- Modify: `include/codosol.h` (include block lines 15-25; the 7 `#ifdef`/`#elif defined(USE_MKL)` sites at lines ~364, 390, 415, 486, 499, 588, 621)

**Interfaces:**
- Consumes: the variant chosen in Task 3.
- Produces: a `codosol.h` with two compile paths only — `USE_SPARSE_JACOBIAN` and a single dense-Eigen path. No `USE_MKL` symbol remains.

The reference replacement below assumes **V2** (dynamic `partialPivLu`). If Task 3 chose V3, additionally gate the Newton-step solve on `n == MAXDIM` (see Step 4).

- [ ] **Step 1: Simplify the include block**

In `include/codosol.h` replace lines 15-21:

```cpp
#ifdef USE_SPARSE_JACOBIAN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#elif !defined(USE_MKL)
#include <Eigen/Dense>
#endif
```

with:

```cpp
#ifdef USE_SPARSE_JACOBIAN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#else
#include <Eigen/Dense>
#endif
```

- [ ] **Step 2: Collapse each `#elif defined(USE_MKL)` site**

At each of the 6 three-way sites (grad `dgemv 'T'`; jdgrad `dgemv 'N'`; the Newton step; `aa += jac*pciv`; `bb = jac*seg`; `fxpjactp += jac*p`), delete the `#elif defined(USE_MKL)` MKL branch entirely and keep the existing `#else` Eigen branch, turning each block from three-way into:

```cpp
#ifdef USE_SPARSE_JACOBIAN
    ...sparse...
#else
    ...existing Eigen expression...
#endif
```

Example — the grad computation (currently lines 364-380) becomes:

```cpp
            #ifdef USE_SPARSE_JACOBIAN
            Eigen::Map<Eigen::MatrixXd> jac_dense(jac, n, n);
            Eigen::SparseMatrix<double> jac_sparse;
            jac_sparse = jac_dense.sparseView();
            jac_sparse.makeCompressed();
            Eigen::Map<Eigen::MatrixXd> fx_dense(fx, n, 1);
            Eigen::Map<Eigen::MatrixXd> grad_dense(grad, n, 1);
            grad_dense = jac_sparse.transpose() * fx_dense;
            #else
            Eigen::Map<Eigen::MatrixXd> jac_dense(jac, n, n);
            Eigen::Map<Eigen::MatrixXd> fx_dense(fx, n, 1);
            Eigen::Map<Eigen::MatrixXd> grad_dense(grad, n, 1);
            grad_dense = jac_dense.transpose() * fx_dense;
            #endif
```

- [ ] **Step 3: Switch the Newton step from QR to LU**

At the Newton-step site, the `#else` branch currently reads:

```cpp
            Eigen::Map<Eigen::MatrixXd> sn_dense(sn, n, 1);
            sn_dense = jac_dense.colPivHouseholderQr().solve(fx_dense);
```

Change `colPivHouseholderQr()` to `partialPivLu()`:

```cpp
            Eigen::Map<Eigen::MatrixXd> sn_dense(sn, n, 1);
            sn_dense = jac_dense.partialPivLu().solve(fx_dense);
```

- [ ] **Step 4: (Only if Task 3 chose V3) gate the fixed-size path**

If V3 won, replace the Newton-step `#else` body with:

```cpp
            Eigen::Map<Eigen::MatrixXd> sn_dense(sn, n, 1);
            if (n == MAXDIM) {
                Eigen::Map<Eigen::Matrix<double,MAXDIM,MAXDIM>> jacN(jac);
                Eigen::Map<Eigen::Matrix<double,MAXDIM,1>> fxN(fx);
                Eigen::Map<Eigen::Matrix<double,MAXDIM,1>> snN(sn);
                snN = jacN.partialPivLu().solve(fxN);
            } else {
                sn_dense = jac_dense.partialPivLu().solve(fx_dense);
            }
```

- [ ] **Step 5: Collapse the Broyden-update site**

The last site (currently lines 621-627) is two-way (`#ifdef USE_MKL` / `#else`). Delete the MKL branch, keep the Eigen `#else` body unconditionally:

```cpp
            Eigen::Map<Eigen::MatrixXd> jacdx_dense(jacdx, n, 1);
            jacdx_dense = jac_dense * Eigen::Map<Eigen::MatrixXd>(deltax, n, 1);
```

- [ ] **Step 6: Verify no `USE_MKL` symbol remains**

Run: `grep -n "USE_MKL" include/codosol.h`
Expected: no output.

- [ ] **Step 7: Commit**

```bash
git add include/codosol.h
git commit -m "perf(codosol): use Eigen partialPivLu for dense Newton step; drop MKL path"
```

---

### Task 5: Drop `-DUSE_MKL` and the `essential_blas.lib` link from the build template

**Files:**
- Modify: `templates/cxx/compile.tpl.m:36,39,49`

**Interfaces:**
- Consumes: the MKL-free `codosol.h` from Task 4.
- Produces: a compile template that no longer defines `USE_MKL` or links `essential_blas.lib`.

- [ ] **Step 1: Remove `-DUSE_MKL` from the MSVC and Intel flag strings**

In `templates/cxx/compile.tpl.m`, the MSVC case (line 36) currently starts:

```matlab
            flag2 = ' -DUSE_MKL OPTIMFLAGS="/O2 /DNDEBUG" COMPFLAGS="$COMPFLAGS /wd4267 /wd4068 /wd4091 /diagnostics:caret /openmp /Z7"';
```

Remove the ` -DUSE_MKL` token:

```matlab
            flag2 = ' OPTIMFLAGS="/O2 /DNDEBUG" COMPFLAGS="$COMPFLAGS /wd4267 /wd4068 /wd4091 /diagnostics:caret /openmp /Z7"';
```

Do the same for the Intel `otherwise` case (line 39): delete the leading ` -DUSE_MKL`.

- [ ] **Step 2: Stop linking `essential_blas.lib`**

In the `switch compiler_name` for `link_to_lib` (lines 45-52), the MSVC case is:

```matlab
        case 'MSVC'
            link_to_lib = sprintf(' "%s/essential_blas.lib"',include_folder);
```

Change it to link nothing:

```matlab
        case 'MSVC'
            link_to_lib = '';
```

- [ ] **Step 3: Sanity-check the template still references its placeholders**

Run: `grep -n "USE_MKL\|essential_blas" templates/cxx/compile.tpl.m`
Expected: no output.

- [ ] **Step 4: Commit**

```bash
git add templates/cxx/compile.tpl.m
git commit -m "build: stop defining USE_MKL / linking essential_blas.lib"
```

---

### Task 6: Regenerate goldens and run the full correctness gate

**Files:**
- Modify (regenerated): `tests/HeatonLucas1996/codegen/golden/compile_HL1996_golden.txt` and any other codegen snapshots that capture the compile string
- Use: `tests/HeatonLucas1996/codegen/regen_snapshots.m`

**Interfaces:**
- Consumes: the changes from Tasks 4-5.
- Produces: regenerated goldens reflecting the new compile string; a green full suite.

- [ ] **Step 1: Inspect which goldens capture the compile flags**

Run: `grep -rln "USE_MKL\|essential_blas" tests/`
Expected: at least `tests/HeatonLucas1996/codegen/golden/compile_HL1996_golden.txt` (the `oldmex/` reference files are a frozen record of the OLD toolbox and must NOT be regenerated).

- [ ] **Step 2: Regenerate the HL1996 codegen snapshots**

Run: `matlab -batch "cd('tests/HeatonLucas1996/codegen'); regen_snapshots"`
Expected: `compile_HL1996_golden.txt` updates so its flag line no longer contains `-DUSE_MKL`, and the golden no longer links `essential_blas.lib`. Confirm with `grep -n "USE_MKL\|essential_blas" tests/HeatonLucas1996/codegen/golden/compile_HL1996_golden.txt` → no output.

- [ ] **Step 3: Run the full test suite**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: exit code 0. Check `tests/results/junit.xml` for zero failures. The end-to-end golden gates (`tEndToEndHL1996`, the safe_assets / Bianchi / CaoNie gates) must pass — `partialPivLu` is the same algorithm as MKL's `dgetrf`, so policy goldens should match within their existing tolerances. If a model that previously relied on colPiv-QR's rank-deficiency tolerance fails, STOP and apply systematic-debugging before proceeding.

- [ ] **Step 4: Commit**

```bash
git add tests/HeatonLucas1996/codegen/golden/
git commit -m "test: regenerate HL1996 compile golden without USE_MKL"
```

---

### Task 7: End-to-end re-time, apply the decision rule, dispose `essential_blas.h`

**Files:**
- Modify: `docs/superpowers/specs/2026-06-24-replace-mkl-linear-algebra-design.md` (`## Results`)
- Modify or delete: `include/essential_blas.h`

**Interfaces:**
- Consumes: the MKL baseline (Task 2), the new Eigen build (Tasks 4-6).
- Produces: a final pass/fail against the ≤5% rule and a disposition for the dead header.

- [ ] **Step 1: Re-time HL1996 on the Eigen build**

Run: `matlab -batch "cd('scratch/bench_mkl'); time_hl1996"`
Expected: a converged run (`Metric` < 1e-6, `Iter` ≈ 209) with a median iter time. Record it in the spec `## Results` as "Eigen end-to-end".

- [ ] **Step 2: Apply the decision rule**

Compute `(eigen_time - mkl_time) / mkl_time`. If ≤ 5% AND the Task 6 suite was green → decision is **adopt Eigen** (already implemented). If > 5% → record the gap; flag to the user that the dependency-removal benefit must be weighed against the measured slowdown (per the spec's fallback). Write the verdict into the spec `## Results`.

- [ ] **Step 3: Dispose of the dead header**

`include/essential_blas.h` is now unreferenced (Task 4 removed its only includer). Confirm: `grep -rln "essential_blas.h" include/ templates/ src/ tests/` → expect no production references (the `scratch/bench_mkl` bench may still include it; that's fine and git-ignored). Delete the header:

```bash
git rm include/essential_blas.h
```

Leave `include/essential_blas.lib` in place (harmless; removing the binary is out of scope and may be referenced by archived builds).

- [ ] **Step 4: Final full-suite confirmation**

Run: `matlab -batch "cd('tests'); run_tests"`
Expected: exit 0, zero failures in `tests/results/junit.xml` — confirming the header deletion broke nothing.

- [ ] **Step 5: Commit**

```bash
git add include/essential_blas.h docs/superpowers/specs/2026-06-24-replace-mkl-linear-algebra-design.md
git commit -m "perf(codosol): drop MKL — record end-to-end verdict, remove dead essential_blas.h"
```

---

## Self-Review

**Spec coverage:**
- Two-tier methodology → Task 1 (micro-bench) + Tasks 2/7 (end-to-end). ✓
- Variant matrix V0-V4 → Task 1 benchmark covers all five + the gemv line item. ✓
- Decision rule (≤5% + golden match) → Task 3 (variant pick) + Task 7 Step 2 (final rule). ✓
- Implementation in `codosol.h` (collapse ifdefs, QR→LU, V3 guard) → Task 4. ✓
- `compile.tpl.m` drops `-DUSE_MKL` + link → Task 5. ✓
- `essential_blas.h` disposition → Task 7 Step 3. ✓
- Regenerate compile goldens → Task 6. ✓
- Full-corpus correctness gate → Task 6 Step 3 + Task 7 Step 4. ✓
- Buffer-layout fact / `n==MAXDIM` guard → Task 4 Step 4. ✓
- Sparse path untouched → every Task-4 edit preserves the `#ifdef USE_SPARSE_JACOBIAN` branch. ✓

**Placeholder scan:** No TBD/TODO; all code shown in full; all commands have expected output. ✓

**Type consistency:** `bench_linalg` struct field names match between `bench_linalg.cpp` and `build_and_run.m`. `jac_dense`/`fx_dense`/`sn_dense` names match the existing `codosol.h` identifiers. `partialPivLu()` used consistently. ✓
