// interp_der_mex.cpp — analytic partial derivative of stacked pp splines.
//
//   ders = interp_der_mex(int32 numThreads, pp, sites, int32 idx, int32 derDim)
//
//   pp     : struct with breaks (1xD cell), coefs, order (int32 1xD, uniform),
//            dim (int32 numArray) — as produced by interp_construct_mex.
//   sites  : D x N matrix of evaluation points (column = one site).
//   idx    : K x 1 or K x N int32, 0-based stacked-vector indices to evaluate.
//   derDim : scalar int32, 1-based state dimension to differentiate along.
//   ders   : K x N, ders(k,n) = d/d(x_derDim) of pp vector idx(k[,n]) at
//            sites(:,n).
//
// The GDSGE_DER (global var_aux) evaluation core: same dispatch and layout as
// interp_eval_mex, but the leaf Horner along derDim is the derivative Horner
// (recder4/recder2 in include/interp_eval_double.h). Uniform order 2 or 4.
#include "mex.h"
#include "interp_eval_double.h"
#ifdef USE_OMP
#include <omp.h>
#endif

#define MAXDIM 10

using gdsge_eval::PpEval;

template <int XDIM, int ORDER>
static void run(const mxArray* pp, const double* sites, int N,
                const int* idx, int K, bool perSite, int derDim,
                double* out, int nThreads)
{
    PpEval<XDIM, ORDER> ev(pp);
#ifdef USE_OMP
    omp_set_num_threads(nThreads);
    #pragma omp parallel for
#endif
    for (int n = 0; n < N; ++n)
    {
        double xSite[MAXDIM], xLeft[MAXDIM];
        int cell[MAXDIM];
        for (int d = 0; d < XDIM; ++d) xSite[d] = sites[(size_t)n * XDIM + d];
        ev.search(xSite, cell, xLeft);
        const int* col = perSite ? idx + (size_t)n * K : idx;
        for (int k = 0; k < K; ++k)
            out[(size_t)n * K + k] = ev.nosearch_eval_der(xLeft, cell, col[k], derDim);
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs != 5) mexErrMsgIdAndTxt("gdsge:interpDer:nargin",
        "interp_der_mex(numThreads, pp, sites, idx, derDim)");
    int nThreads = *((int*)mxGetData(prhs[0]));
    const mxArray* pp = prhs[1];
    const double* sites = mxGetPr(prhs[2]);
    int XDIM = (int)mxGetM(prhs[2]);
    int N    = (int)mxGetN(prhs[2]);
    const int* idx = (const int*)mxGetData(prhs[3]);
    int K = (int)mxGetM(prhs[3]);
    bool perSite = (mxGetN(prhs[3]) == (mwSize)N) && N > 1;
    int derDim = *((int*)mxGetData(prhs[4])) - 1;   // 1-based from MATLAB
    if (derDim < 0 || derDim >= XDIM)
        mexErrMsgIdAndTxt("gdsge:interpDer:derDim",
            "derDim %d out of range for state dim %d", derDim + 1, XDIM);

    int order0 = *((int*)mxGetData(mxGetField(pp, 0, "order")));

    plhs[0] = mxCreateDoubleMatrix(K, N, mxREAL);
    double* out = mxGetPr(plhs[0]);

    #define DISPATCH(D) \
        if (XDIM == D) { if (order0 == 4) run<D,4>(pp,sites,N,idx,K,perSite,derDim,out,nThreads); \
                         else              run<D,2>(pp,sites,N,idx,K,perSite,derDim,out,nThreads); return; }
    DISPATCH(1) DISPATCH(2) DISPATCH(3) DISPATCH(4) DISPATCH(5) DISPATCH(6) DISPATCH(7) DISPATCH(8)
    #undef DISPATCH
    mexErrMsgIdAndTxt("gdsge:interpDer:dim", "unsupported state dim %d", XDIM);
}
