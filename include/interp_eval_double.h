#pragma once
// interp_eval_double.h — double-only tensor-spline evaluator for a Matlab pp
// struct (breaks/coefs/order/dim), shared by interp_eval_mex.cpp and the
// generated simulate_<model>_mex.cpp.
//
// hunt() and receval2/receval4() are copied VERBATIM from include/interp_lite.h
// (the same arithmetic the model MEX's MatlabPp uses) so results are numerically
// identical — but without interp_lite.h's `#include "adept_source.h"`, since the
// simulation path needs no autodiff. Uniform interpolation order only (every
// breaks dim shares the order), matching the stacked uniform-order layout that
// interp_construct_mex builds for output_interp / the resolve warm-start interp.
#include "mex.h"

namespace gdsge_eval {

#ifndef GDSGE_EVAL_MAXDIM
#define GDSGE_EVAL_MAXDIM 10
#endif

// --- verbatim from interp_lite.h (1-based grid convention) ---
inline void hunt(double xx[], int n, double x, int* jlo)
{
    int jm, jhi, inc;
    if (*jlo <= 0 || *jlo > n) { *jlo = 0; jhi = n + 1; }
    else {
        inc = 1;
        if (x >= xx[*jlo]) {
            if (*jlo == n) return;
            jhi = (*jlo) + 1;
            while (x >= xx[jhi]) {
                *jlo = jhi; inc += inc; jhi = (*jlo) + inc;
                if (jhi > n) { jhi = n + 1; break; }
            }
        }
        else {
            if (*jlo == 1) { *jlo = 0; return; }
            jhi = (*jlo)--;
            while (x < xx[*jlo]) {
                jhi = (*jlo); inc <<= 1;
                if (inc >= jhi) { *jlo = 0; break; }
                else *jlo = jhi - inc;
            }
        }
    }
    while (jhi - (*jlo) != 1) {
        jm = (jhi + (*jlo)) >> 1;
        if (x >= xx[jm]) *jlo = jm; else jhi = jm;
    }
    if (x == xx[n]) *jlo = n - 1;
    if (x == xx[1]) *jlo = 1;
}

inline double receval4(double* coefs, int* coefsSize, int xDim, int idim,
                       double* xsite, int* cellOfSite, int shift)
{
    double r;
    if (idim == xDim - 1) {
        double* pCoefs = coefs + shift + (*cellOfSite + 1) * 4;
        r = *(--pCoefs); r *= *xsite;
        r += *(--pCoefs); r *= *xsite;
        r += *(--pCoefs); r *= *xsite;
        r += *(--pCoefs);
    }
    else {
        shift += (*cellOfSite + 1) * 4;
        r = receval4(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
        r *= *xsite;
        r += receval4(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
        r *= *xsite;
        r += receval4(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
        r *= *xsite;
        r += receval4(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
    }
    return r;
}

inline double receval2(double* coefs, int* coefsSize, int xDim, int idim,
                       double* xsite, int* cellOfSite, int shift)
{
    double r;
    if (idim == xDim - 1) {
        double* pCoefs = coefs + shift + (*cellOfSite + 1) * 2;
        r = *(--pCoefs); r *= *xsite;
        r += *(--pCoefs);
    }
    else {
        shift += (*cellOfSite + 1) * 2;
        r = receval2(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
        r *= *xsite;
        r += receval2(coefs, coefsSize, xDim, idim + 1, xsite + 1, cellOfSite + 1, (--shift) * coefsSize[idim + 2]);
    }
    return r;
}

// pp evaluator: read a Matlab pp (interp_construct_mex / myppual output), then
// search once per site and nosearch_eval each stacked vector index.
template <int XDIM, int ORDER>
struct PpEval {
    double* xGrid[XDIM];
    int xPts[XDIM];
    int yDim;
    double* coefs;
    int coefsSize[XDIM + 1];

    PpEval(const mxArray* pp)
    {
        const mxArray* ppBreaks = mxGetField(pp, 0, "breaks");
        coefs = mxGetPr(mxGetField(pp, 0, "coefs"));
        yDim  = *((int*)mxGetData(mxGetField(pp, 0, "dim")));
        for (int j = 0; j < XDIM; ++j) {
            const mxArray* b = mxGetCell(ppBreaks, j);
            xGrid[j] = mxGetPr(b);
            xPts[j]  = (int)mxGetNumberOfElements(b);
        }
        coefsSize[0] = yDim;
        for (int i = 0; i < XDIM; ++i) coefsSize[i + 1] = (xPts[i] - 1) * ORDER;
    }

    void search(double* xSite, int* cellOfSite, double* xSiteToLeft)
    {
        for (int j = 0; j < XDIM; ++j) {
            hunt(xGrid[j], xPts[j] - 2, xSite[j], cellOfSite + j);
            xSiteToLeft[j] = xSite[j] - xGrid[j][cellOfSite[j]];
        }
    }

    double nosearch_eval(double* xSiteToLeft, int* cellOfSite, int vecIdx)
    {
        if (ORDER == 4)
            return receval4(coefs, coefsSize, XDIM, 0, xSiteToLeft, cellOfSite, vecIdx * coefsSize[1]);
        else
            return receval2(coefs, coefsSize, XDIM, 0, xSiteToLeft, cellOfSite, vecIdx * coefsSize[1]);
    }
};

} // namespace gdsge_eval
