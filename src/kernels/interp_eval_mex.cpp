// interp_eval_mex.cpp — generic uniform-order evaluator for stacked pp structs.
//
//   vals = interp_eval_mex(int32 numThreads, pp, sites, int32 idx)
//
//   pp     : struct with breaks (1xD cell), coefs, order (int32 1xD, uniform),
//            dim (int32 numArray) — as produced by interp_construct_mex.
//   sites  : D x N matrix of evaluation points (column = one site).
//   idx    : K x 1 or K x N int32, 0-based stacked-vector indices to evaluate.
//            K x 1  -> same K vectors at every site;
//            K x N  -> per-site vector indices (one column per site).
//   vals   : K x N, vals(k,n) = pp vector idx(k[,n]) evaluated at sites(:,n).
//
// Uniform order only (every breaks dim shares order[0]); supports order 2 or 4.
// Evaluation core is include/interp_eval_double.h (double-only copy of the
// model MEX's MatlabPp arithmetic — same results, no adept dependency).
#include "mex.h"
#include "interp_eval_double.h"
#ifdef USE_OMP
#include <omp.h>
#endif

#define MAXDIM 10

using gdsge_eval::PpEval;

template <int XDIM, int ORDER>
static void run(const mxArray* pp, const double* sites, int N,
                const int* idx, int K, bool perSite, double* out, int nThreads)
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
            out[(size_t)n * K + k] = ev.nosearch_eval(xLeft, cell, col[k]);
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs != 4) mexErrMsgIdAndTxt("gdsge:interpEval:nargin",
        "interp_eval_mex(numThreads, pp, sites, idx)");
    int nThreads = *((int*)mxGetData(prhs[0]));
    const mxArray* pp = prhs[1];
    const double* sites = mxGetPr(prhs[2]);
    int XDIM = (int)mxGetM(prhs[2]);
    int N    = (int)mxGetN(prhs[2]);
    const int* idx = (const int*)mxGetData(prhs[3]);
    int K = (int)mxGetM(prhs[3]);
    bool perSite = (mxGetN(prhs[3]) == (mwSize)N) && N > 1;

    int order0 = *((int*)mxGetData(mxGetField(pp, 0, "order")));

    plhs[0] = mxCreateDoubleMatrix(K, N, mxREAL);
    double* out = mxGetPr(plhs[0]);

    #define DISPATCH(D) \
        if (XDIM == D) { if (order0 == 4) run<D,4>(pp,sites,N,idx,K,perSite,out,nThreads); \
                         else              run<D,2>(pp,sites,N,idx,K,perSite,out,nThreads); return; }
    DISPATCH(1) DISPATCH(2) DISPATCH(3) DISPATCH(4) DISPATCH(5) DISPATCH(6) DISPATCH(7) DISPATCH(8)
    #undef DISPATCH
    mexErrMsgIdAndTxt("gdsge:interpEval:dim", "unsupported state dim %d", XDIM);
}
