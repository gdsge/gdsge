// interp_construct_mex.cpp — minimal fused spline/linear constructor.
//
// [ppCell, splineVec] = interp_construct_mex(breaks, values, order, extrapOrder, numThreads)
//
//   breaks      : 1xD cell of grid row vectors (the shared state grid)
//   values      : struct, one field per interp variable; each field is that
//                 variable's value array [numArray, gridDims...] (no concat copy)
//   order       : int32 [1 x D]   (4 = cubic not-a-knot, 2 = linear)
//   extrapOrder : int32 [1 x D] or [] (cubic extrapolation order; [] = none)
//   numThreads  : int32 scalar
//
//   ppCell      : 1xnumVec cell of natural-order pp structs (== myppual(pp))
//   splineVec   : GDSGE_SPLINE_VEC (eval-order coefs == convert_to_interp_eval_array)
//
// The construction driver (mkl_start / mkl_start_with_extrap + helpers) is copied
// VERBATIM from the owner's myppual_mex.cpp — do NOT edit those bodies (bit-exact).
// Compile: mex interp_construct_mex.cpp -DUSE_OMP -Iinclude COMPFLAGS="$COMPFLAGS /openmp"
#include "mex.h"
#include "mkl_dummy_interp.h"   // construct_*_spline / mkl_domatcopy / dfd* / DF_* / MKL_INT / DFTaskPtr
#include <string.h>
#ifdef USE_OMP
#include <omp.h>
#endif

#define MAXDIM 10
#define MAX(a,b) ((a>b) ? a : b)
#define MIN(a,b) ((a<b) ? a : b)

// ===================== verbatim helpers (myppual_mex.cpp 2740-2760) =====================
void intmemcpy(int* des, int* src, int n) { memcpy(des, src, sizeof(int) * n); }
inline void doublememcpy(double* des, double* src, int n) { memcpy(des, src, sizeof(double) * n); }
inline int vectorProd(int* x, int n) { int result = x[0]; for (int i = 1; i < n; ++i) result *= x[i]; return result; }

// ===================== verbatim: mkl_start_with_extrap (myppual_mex.cpp 2476-2646) =====================
void mkl_start_with_extrap(int xDim, int* xPts, double** xGrid, int valDim, double* y, int* s_order, double* coeff, int coeffn,
	int* extrap_order)
{
	// copy from Wenlan's NdInterp::start(), modified without using Matrix class
	// extend grid and copy coefficients for extrapolation
	double* lastCoeff = (double*)malloc(sizeof(double) * coeffn);
	int lastCoeffns[MAXDIM];
	int lastCoeffn;

	double* extendedCoeff = (double*)malloc(sizeof(double) * coeffn);

	// for the first coordinate, the lastCoeff is just y
	int initialCoeffSize[MAXDIM + 1];
	initialCoeffSize[0] = valDim;
	for (int i = 0; i < xDim; ++i) {
		initialCoeffSize[i + 1] = xPts[i];
	}

	int coeffns[MAXDIM];
	intmemcpy(coeffns, initialCoeffSize, xDim + 1);
	coeffn = vectorProd(coeffns, xDim + 1);
	doublememcpy(coeff, y, coeffn);

	// the following is a mimic of csape
	MKL_INT s_type[MAXDIM];
	MKL_INT bc_type[MAXDIM];
	for (int i = 0; i < xDim; ++i) {
		switch (s_order[i]) {
		case 4:
			s_type[i] = DF_PP_NATURAL;
			bc_type[i] = DF_BC_NOT_A_KNOT;
			break;
		case 2:
			s_type[i] = DF_PP_DEFAULT;
			bc_type[i] = DF_NO_BC;
			break;
		}
	}

	// Let me do Mkl Interp at low level to boost performance
	DFTaskPtr interp0;
	int currentVecSize = 1; // not used
	for (int i = xDim - 1; i >= 0; --i) {
		// carry out coordinatewise interpolation at coordinate i
		// the interpolation is always w.r.t. the last coordinate, w.r.t. to which y are stored adjacently
		// calculate new size of vector function, i.e. the total size except the last dimension
		currentVecSize = vectorProd(coeffns, xDim);

		// interp0 = new MklInterp(interpType, xGrid[i], xPts[i], lastCoeff->pt(), currentVecSize, 'n');
		// direct coefficients of MklInterp to the array in NdInterp

		// compute the new size of coeffs
		// the last dimension reduces by 1, and times by s_order
		coeffns[xDim] = coeffns[xDim] - 1;
		coeffns[xDim] = coeffns[xDim] * s_order[i];
		// update last Coeff
		intmemcpy(lastCoeffns, coeffns, xDim + 1);
		lastCoeffn = vectorProd(lastCoeffns, xDim + 1);

		// compute coefficients;
		dfdNewTask1D(&interp0, xPts[i], xGrid[i], DF_NO_HINT, currentVecSize, coeff, DF_NO_HINT);
		dfdEditPPSpline1D(interp0, s_order[i], s_type[i], bc_type[i], 0, DF_NO_IC, 0, lastCoeff, DF_NO_HINT);
		dfdConstruct1D(interp0, DF_PP_SPLINE, DF_METHOD_STD);
		dfDeleteTask(&interp0);

		// At this point lastCoeff has memory order (C-order) (vec, grid, order)
		// Now I'm going to plug in the extrapolation coefficient here
		// This is done by expanding the grid below and above, for EACH vec
		// I have to do this in a new array (pre-allocated)
		int vecLength = lastCoeffn / lastCoeffns[xDim];
		int blockLength = lastCoeffns[xDim];
        if (s_order[i] == 4)
        {
            // only do extrapolation for cubic splines
            if (extrap_order[i] == 3)
            {
                #pragma omp parallel for
                for (int i_vec=0; i_vec < vecLength; ++i_vec)
                {
                    int ptr = (blockLength + 2 * s_order[i]) * i_vec;
                    memcpy(extendedCoeff + ptr + s_order[i], lastCoeff+i_vec*blockLength, sizeof(double)*blockLength);
                    // Copy coefficient of the piece and the last piece
                    {
                        double* cn = extendedCoeff + ptr;
                        double* c = cn + s_order[i];
                        cn[3] = 0;
                        cn[2] = c[2];
                        cn[1] = c[1] - 2*cn[2];
                        cn[0] = c[0] - cn[1] - cn[2];
                    }

                    {
                        double* cn = extendedCoeff + ptr + blockLength + s_order[i];
                        double* c = cn - s_order[i];
                        double delta = xGrid[i][xPts[i]-1] - xGrid[i][xPts[i]-2];
                        cn[0] = c[0] + delta*c[1] + delta*delta*c[2] + delta*delta*delta*c[3];
                        cn[1] = c[1] + 2*delta*c[2] + 3*delta*delta*c[3];
                        cn[2] = c[2] + 3*delta*c[3];
                        cn[3] = 0;
                    }
                }
            }
            else if (extrap_order[i]==2)
            {
                for (int i_vec=0; i_vec < vecLength; ++i_vec)
                {
                    int ptr = (blockLength + 2 * s_order[i]) * i_vec;
                    memcpy(extendedCoeff + ptr + s_order[i], lastCoeff+i_vec*blockLength, sizeof(double)*blockLength);
                    // Copy coefficient of the piece and the last piece
                    {
                        double* cn = extendedCoeff + ptr;
                        double* c = cn + s_order[i];
                        cn[3] = 0;
                        cn[2] = 0;
                        cn[1] = c[1];
                        cn[0] = c[0] - cn[1];
                    }

                    {
                        double* cn = extendedCoeff + ptr + blockLength + s_order[i];
                        double* c = cn - s_order[i];
                        double delta = xGrid[i][xPts[i]-1] - xGrid[i][xPts[i]-2];
                        cn[0] = c[0] + delta*c[1] + delta*delta*c[2] + delta*delta*delta*c[3];
                        cn[1] = c[1] + 2*delta*c[2] + 3*delta*delta*c[3];
                        cn[2] = 0;
                        cn[3] = 0;
                    }

                    // shift pointer forward
                    ptr += (blockLength + 2*s_order[i]);
                }
            }
            // Update coefs
            // Update the coefs size after extrapolation
            lastCoeffns[xDim] += 2*s_order[i];
            lastCoeffn = vectorProd(lastCoeffns, xDim + 1);
            memcpy(lastCoeff, extendedCoeff, sizeof(double)*lastCoeffn);
        }

		// forming new interpolation problem, by permuting the last dimension to the first one
		// The last dimension is ((xPts-1)*order)
		if (xDim > 1) {
			int vecShift = lastCoeffn / valDim;

			double* des = coeff;
			double* src = lastCoeff;
			int transRow = vecShift / lastCoeffns[xDim];
			int transCol = lastCoeffns[xDim];
			#pragma omp parallel for
			for (int i = 0; i < valDim; ++i) {
				mkl_domatcopy('R', 'T', transRow, transCol, 1.0, src + i * vecShift, transCol, des + i * vecShift, transRow);
			}
			// swap dimension
			int temp = lastCoeffns[xDim];
			for (int i = xDim; i >= 2; --i) {
				lastCoeffns[i] = lastCoeffns[i - 1];
			}

			lastCoeffns[1] = temp;
			intmemcpy(coeffns, lastCoeffns, xDim + 1);
			coeffn = vectorProd(coeffns, xDim + 1);
		}
		else {
			doublememcpy(coeff, lastCoeff, lastCoeffn);
		}
	}

	free(lastCoeff);
	free(extendedCoeff);
}

// ===================== verbatim: mkl_start (myppual_mex.cpp 2647-2739) =====================
void mkl_start(int xDim, int* xPts, double** xGrid, int valDim, double* y, int* s_order, double* coeff, int coeffn)
{
	// copy from Wenlan's NdInterp::start(), modified without using Matrix class
	double* lastCoeff = (double*)malloc(sizeof(double) * coeffn);
	int lastCoeffns[MAXDIM];
	int lastCoeffn;

	// for the first coordinate, the lastCoeff is just y
	int initialCoeffSize[MAXDIM + 1];
	initialCoeffSize[0] = valDim;
	for (int i = 0; i < xDim; ++i) {
		initialCoeffSize[i + 1] = xPts[i];
	}

	int coeffns[MAXDIM];
	intmemcpy(coeffns, initialCoeffSize, xDim + 1);
	coeffn = vectorProd(coeffns, xDim + 1);
	doublememcpy(coeff, y, coeffn);

	// the following is a mimic of csape
	MKL_INT s_type[MAXDIM];
	MKL_INT bc_type[MAXDIM];
	for (int i = 0; i < xDim; ++i) {
		switch (s_order[i]) {
		case 4:
			s_type[i] = DF_PP_NATURAL;
			bc_type[i] = DF_BC_NOT_A_KNOT;
			break;
		case 2:
			s_type[i] = DF_PP_DEFAULT;
			bc_type[i] = DF_NO_BC;
			break;
		}
	}

	// Let me do Mkl Interp at low level to boost performance
	DFTaskPtr interp0;
	int currentVecSize = 1; // not used
	for (int i = xDim - 1; i >= 0; --i) {
		// carry out coordinatewise interpolation at coordinate i
		// the interpolation is always w.r.t. the last coordinate, w.r.t. to which y are stored adjacently
		// calculate new size of vector function, i.e. the total size except the last dimension
		currentVecSize = vectorProd(coeffns, xDim);

		// interp0 = new MklInterp(interpType, xGrid[i], xPts[i], lastCoeff->pt(), currentVecSize, 'n');
		// direct coefficients of MklInterp to the array in NdInterp

		// compute the new size of coeffs
		// the last dimension reduces by 1, and times by s_order
		coeffns[xDim] = coeffns[xDim] - 1;
		coeffns[xDim] = coeffns[xDim] * s_order[i];
		// update last Coeff
		intmemcpy(lastCoeffns, coeffns, xDim + 1);
		lastCoeffn = vectorProd(lastCoeffns, xDim + 1);

		// compute coefficients
		dfdNewTask1D(&interp0, xPts[i], xGrid[i], DF_NO_HINT, currentVecSize, coeff, DF_NO_HINT);
		dfdEditPPSpline1D(interp0, s_order[i], s_type[i], bc_type[i], 0, DF_NO_IC, 0, lastCoeff, DF_NO_HINT);
		dfdConstruct1D(interp0, DF_PP_SPLINE, DF_METHOD_STD);
		dfDeleteTask(&interp0);

		// forming new interpolation problem, by permuting the last dimension to the first one
		if (xDim > 1) {
			int vecShift = lastCoeffn / valDim;

			double* des = coeff;
			double* src = lastCoeff;
			int transRow = vecShift / lastCoeffns[xDim];
			int transCol = lastCoeffns[xDim];
			for (int i = 0; i < valDim; ++i) {
				mkl_domatcopy('R', 'T', transRow, transCol, 1.0, src, transCol, des, transRow);
				des += vecShift;
				src += vecShift;
			}
			// swap dimension
			int temp = lastCoeffns[xDim];
			for (int i = xDim; i >= 2; --i) {
				lastCoeffns[i] = lastCoeffns[i - 1];
			}

			lastCoeffns[1] = temp;
			intmemcpy(coeffns, lastCoeffns, xDim + 1);
			coeffn = vectorProd(coeffns, xDim + 1);
		}
		else {
			doublememcpy(coeff, lastCoeff, lastCoeffn);
		}
	}

	free(lastCoeff);
}

// Construct one interp variable's natural-order coefs into `coeff`. `valuesMKL`
// is the value array already permuted to MKL orientation. Mirrors
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
    if (hasExtrap && (int) mxGetNumberOfElements(extrap_) != xDim)
        mexErrMsgIdAndTxt("gdsge:kernels:interpConstructExtrapLen",
            "interp_construct_mex: extrapOrder must have one entry per state dimension");
#ifdef USE_OMP
    omp_set_num_threads(numThreads > 0 ? numThreads : 1);
#endif

    // --- grids + (cubic) extrap break-extension. Mirror myppual.m:626-632
    //     (extend each cubic break by one knot at each end, pieces += 2) and
    //     fullConstructionWithExtrap:313-321 (construction consumes the interior grid).
    // Extra knots are attached ONLY when the extrapolant is LOWER order than the
    // cubic spline: extrap_order 2 (linear) or 3 (quadratic). extrap_order==4
    // (cubic) needs none — the boundary cubic extrapolates naturally, so reuse the
    // plain construction (the old myppual errors on extrap_order=4; we handle it
    // cleanly as plain cubic, no extra knots).
    double* xGridExt[MAXDIM];     // breaks carried by the pp structs (extended if cubic+flatten-extrap)
    int     xPtsExt[MAXDIM];
    double* xGridCon[MAXDIM];     // grid passed to mkl_start[_with_extrap]
    int     xPtsCon[MAXDIM];
    int     pieces[MAXDIM];       // post-extension pieces per dim (xPtsExt-1)
    double* extBuf[MAXDIM];       // owned extended-knot buffers
    bool    anyExtend = false;
    for (int d = 0; d < xDim; ++d) {
        const mxArray* g = mxGetCell(breaks_, d);
        double* gp = mxGetPr(g);
        int     gn = (int) mxGetNumberOfElements(g);
        bool extendThis = hasExtrap && s_order[d] == 4 &&
                          (extrap_order[d] == 2 || extrap_order[d] == 3);
        if (extendThis) {
            extBuf[d] = (double*) mxMalloc(sizeof(double) * (gn + 2));
            extBuf[d][0] = gp[0] - 1.0;
            memcpy(extBuf[d] + 1, gp, sizeof(double) * gn);
            extBuf[d][gn + 1] = gp[gn - 1] + 1.0;
            xGridExt[d] = extBuf[d];  xPtsExt[d] = gn + 2;
            xGridCon[d] = gp;         xPtsCon[d] = gn;
            anyExtend = true;
        } else {
            extBuf[d] = NULL;
            xGridExt[d] = gp;         xPtsExt[d] = gn;
            xGridCon[d] = gp;         xPtsCon[d] = gn;
        }
        pieces[d] = xPtsExt[d] - 1;
    }
    // Mixed flatten/cubic extrapolation across dims is not supported (constructSplines
    // always passes a uniform extrap order); guard against silent garbage.
    if (anyExtend) {
        for (int d = 0; d < xDim; ++d)
            if (s_order[d] == 4 && !(extrap_order[d] == 2 || extrap_order[d] == 3))
                mexErrMsgIdAndTxt("gdsge:kernels:interpConstructMixedExtrap",
                    "interp_construct_mex: mixed flatten/cubic extrap orders are unsupported");
    }

    plhs[0] = mxCreateCellMatrix(1, numVec);   // ppCell
    double** natCoefsPerVar = (double**) mxMalloc(sizeof(double*) * numVec);
    int numArray = 0;
    for (int iv = 0; iv < numVec; ++iv) {
        const mxArray* field = mxGetFieldByNumber(values_, 0, iv);   // [numArray, g0pts, ...]
        numArray = (int) mxGetM(field);                             // array dim is first

        // Values_MKL = permute(values, [xDim+1:-1:1]) as an mxArray (kept alive as
        // the pp's "Values"). natural fastest-first dims [numArray, g0, ..., g(D-1)];
        // reversed result dims [g(D-1), ..., g0, numArray]. Construction uses the
        // ORIGINAL grid points (xPtsCon), not the extended ones.
        double* src = mxGetPr(field);
        int natDims[MAXDIM + 1];                 // fastest first
        natDims[0] = numArray;
        for (int d = 0; d < xDim; ++d) natDims[d + 1] = xPtsCon[d];
        mwSize revDims[MAXDIM + 1];              // reversed, fastest first
        for (int d = 0; d < xDim; ++d) revDims[d] = xPtsCon[xDim - 1 - d];
        revDims[xDim] = numArray;
        mxArray* valuesMKLArr = mxCreateNumericArray(xDim + 1, revDims, mxDOUBLE_CLASS, mxREAL);
        double* valuesMKL = mxGetPr(valuesMKLArr);
        int N = numArray; for (int d = 0; d < xDim; ++d) N *= xPtsCon[d];
        for (int s = 0; s < N; ++s) {
            int t = s, idx[MAXDIM + 1];
            for (int k = 0; k <= xDim; ++k) { idx[k] = t % natDims[k]; t /= natDims[k]; }
            int dst = 0, stride = 1;             // dest fastest = idx[xDim]; numArray slowest
            for (int k = xDim; k >= 1; --k) { dst += idx[k] * stride; stride *= natDims[k]; }
            dst += idx[0] * stride;
            valuesMKL[dst] = src[s];
        }

        // coefs size (mirror fullConstruction:401-410 / extrap variant:340-356)
        mwSize coeffns[MAXDIM + 1];
        coeffns[xDim] = numArray;                // fortran order: array dim last
        int coeffn = numArray;
        for (int d = 0; d < xDim; ++d) {
            coeffns[xDim - 1 - d] = (mwSize) (pieces[d] * s_order[d]);   // pieces post-extension
            coeffn *= (int) coeffns[xDim - 1 - d];
        }
        mxArray* coefsArr = mxCreateNumericArray(xDim + 1, coeffns, mxDOUBLE_CLASS, mxREAL);
        double* coeff = mxGetPr(coefsArr);
        constructOne(xDim, xPtsCon, xGridCon, numArray, valuesMKL, s_order,
                     anyExtend ? extrap_order : NULL, coeff, coeffn);
        natCoefsPerVar[iv] = coeff;              // alias into the pp's coefs (lives in plhs[0])

        // pp struct == myppual(pp)'s MKLpp field set (myppual.m:636-638): frozen IterRslt.pp shape.
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

    // ===================== eval-order pack -> GDSGE_SPLINE_VEC (plhs[1]) =====================
    if (nlhs > 1) {
        int prodOrder = 1, prodPieces = 1;
        for (int d = 0; d < xDim; ++d) { prodOrder *= s_order[d]; prodPieces *= pieces[d]; }
        int fullVecEvalCoefsLength   = numVec * prodOrder;
        int singleVecEvalCoefsLength = prodOrder;                 // == fullVec/numVec (convert divides by dim)
        int arrayOffset              = numVec * prodOrder * prodPieces;
        int L                        = prodOrder * prodPieces * numArray;   // per-var coefs length

        // Eval-base stride (in numVec-block units, i.e. eval index = iv + numVec*base)
        // for each natural dimension. order dims: evalStrideO[d]=prod(order[0..d-1]);
        // piece dims: evalStrideP[d]=prodOrder*prod(pieces[0..d-1]); array dim: prodOrder*prodPieces.
        int evalStrideO[MAXDIM], evalStrideP[MAXDIM];
        int so = 1;
        for (int d = 0; d < xDim; ++d) { evalStrideO[d] = so; so *= s_order[d]; }   // so -> prodOrder
        int sp = so;
        for (int d = 0; d < xDim; ++d) { evalStrideP[d] = sp; sp *= pieces[d]; }     // sp -> prodOrder*prodPieces
        // Odometer digits in NATURAL memory order (fastest first): o_D, p_D,
        // o_{D-1}, p_{D-1}, ..., o_1, p_1, a — matching the natural coefs layout.
        int digitSize[2*MAXDIM + 1], digitEvalStride[2*MAXDIM + 1];
        int j = 0;
        for (int d = xDim - 1; d >= 0; --d) {
            digitSize[j] = s_order[d]; digitEvalStride[j] = evalStrideO[d]; ++j;   // o_d
            digitSize[j] = pieces[d];  digitEvalStride[j] = evalStrideP[d]; ++j;   // p_d
        }
        digitSize[j] = numArray; digitEvalStride[j] = sp;                          // a
        int nDigits = 2*xDim + 1;

        mwSize evalDims[2] = { (mwSize)((mwSize)arrayOffset * numArray), 1 };  // flat column
        mxArray* coefsOut = mxCreateNumericArray(2, evalDims, mxDOUBLE_CLASS, mxREAL);
        double* evalCoefs = mxGetPr(coefsOut);

        // Single pass: walk the natural buffer sequentially (cache-friendly per-var
        // reads), accumulate the eval base via an odometer (no per-element
        // modulo/divide), and write each numVec burst contiguously — one scattered
        // write-line per natIdx instead of numVec full scattered sweeps.
        int digits[2*MAXDIM + 1] = {0};
        int base = 0;
        for (int natIdx = 0; natIdx < L; ++natIdx) {
            double* dst = evalCoefs + (size_t)numVec * base;
            for (int iv = 0; iv < numVec; ++iv) dst[iv] = natCoefsPerVar[iv][natIdx];
            for (int jd = 0; jd < nDigits; ++jd) {
                base += digitEvalStride[jd];
                if (++digits[jd] < digitSize[jd]) break;
                digits[jd] = 0;
                base -= digitSize[jd] * digitEvalStride[jd];
            }
        }

        const char* sf[] = {"coefs","fullVecEvalCoefsLength","singleVecEvalCoefsLength",
                            "xDim","order","pieces","xPts","breaks","dim","arrayOffset"};
        mxArray* sv = mxCreateStructMatrix(1, 1, 10, sf);
        mxSetField(sv, 0, "coefs", coefsOut);
        mxSetField(sv, 0, "fullVecEvalCoefsLength",   mxCreateDoubleScalar((double) fullVecEvalCoefsLength));
        mxSetField(sv, 0, "singleVecEvalCoefsLength", mxCreateDoubleScalar((double) singleVecEvalCoefsLength));
        mxSetField(sv, 0, "xDim", mxCreateDoubleScalar((double) xDim));
        mxArray* ordD = mxCreateNumericMatrix(1, xDim, mxINT32_CLASS, mxREAL);  // C++ reads via mxGetData as int
        mxArray* pieD = mxCreateNumericMatrix(1, xDim, mxINT32_CLASS, mxREAL);
        mxArray* xptD = mxCreateNumericMatrix(1, xDim, mxINT32_CLASS, mxREAL);
        int* ordP = (int*) mxGetData(ordD);
        int* pieP = (int*) mxGetData(pieD);
        int* xptP = (int*) mxGetData(xptD);
        for (int d = 0; d < xDim; ++d) { ordP[d] = s_order[d]; pieP[d] = pieces[d]; xptP[d] = pieces[d] + 1; }
        mxSetField(sv, 0, "order",  ordD);
        mxSetField(sv, 0, "pieces", pieD);
        mxSetField(sv, 0, "xPts",   xptD);
        mxArray* bcellSV = mxCreateCellMatrix(1, xDim);
        for (int d = 0; d < xDim; ++d) {
            mxArray* gd = mxCreateDoubleMatrix(1, xPtsExt[d], mxREAL);
            memcpy(mxGetPr(gd), xGridExt[d], sizeof(double) * xPtsExt[d]);
            mxSetCell(bcellSV, d, gd);
        }
        mxSetField(sv, 0, "breaks", bcellSV);
        mxSetField(sv, 0, "dim", mxCreateDoubleScalar((double) numVec));
        mxSetField(sv, 0, "arrayOffset", mxCreateDoubleScalar((double) arrayOffset));
        plhs[1] = sv;
    }

    for (int d = 0; d < xDim; ++d) if (extBuf[d]) mxFree(extBuf[d]);
    mxFree(natCoefsPerVar);
}
