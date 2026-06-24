/* 
Convert Matlab InterpEval to C Struct
@Author: Wenlan Luo
@Date: 8/5/2019
*/

#include "mex.h"
#include "adept_source.h"
#include "InterpEval.h"
#include <string.h>

using adept::adouble;

//pow
template <class T>
constexpr T pow2(T exponent) {
	return (T(1) << exponent);
}

template <int INTERP_EVAL_XDIM, int INTERP_EVAL_NUMVEC>
class MatlabInterpEval {
public:
	int xDim;
	int numVec;
	int arrayOffset;
	int fullVecEvalCoefsLength;
	int singleVecEvalCoefsLength;
	double* xGrid[INTERP_EVAL_XDIM];
	int order[INTERP_EVAL_XDIM];
	int pieces[INTERP_EVAL_XDIM];
	int xPts[INTERP_EVAL_XDIM];
	double* coefs;


	inline void search_eval_vec_at_array(int arrayIdx, double* xSite, double* evalResult, int* cellOfSite)
	{
		InterpEval::search_eval_vec<INTERP_EVAL_XDIM, INTERP_EVAL_NUMVEC>(numVec, xSite, xDim, xGrid, xPts, coefs + arrayIdx*arrayOffset, fullVecEvalCoefsLength,
			order, pieces, evalResult, cellOfSite);
	}

	inline void search_eval_with_grad_vec_at_array(int arrayIdx, double* xSite, double* evalResult, double* grad, int* cellOfSite)
	{
		InterpEval::search_eval_with_grad_vec<INTERP_EVAL_XDIM, INTERP_EVAL_NUMVEC>(xSite, xGrid, xPts, coefs + arrayIdx*arrayOffset, fullVecEvalCoefsLength,
			order, pieces, evalResult, grad, cellOfSite);
	}

	inline void search_eval_vec_at_array_adouble(int arrayIdx, adouble* xSiteAdouble, adouble* evalResultAdouble, double* grad, int* cellOfSite, int evalGradFlag)
	{
		double xSite[INTERP_EVAL_XDIM];
		double evalResult[INTERP_EVAL_NUMVEC];

		for (int i_dim = 0; i_dim < xDim ; i_dim++)
		{
			xSite[i_dim] = value(xSiteAdouble[i_dim]);
		}
		if (evalGradFlag == 1)
		{
			search_eval_with_grad_vec_at_array(arrayIdx, xSite, evalResult, grad, cellOfSite);
		}
		else
		{
			search_eval_vec_at_array(arrayIdx, xSite, evalResult, cellOfSite);
		}
		for (int i_vec = 0; i_vec < numVec ; i_vec++)
		{
			evalResultAdouble[i_vec].set_value(evalResult[i_vec]);
			evalResultAdouble[i_vec].add_derivative_dependence(xSiteAdouble, grad + i_vec, xDim, numVec);
		}
	}

	inline void search_eval_vec(double* xSite, double* evalResult)
	{
		search_eval_vec_at_array(0, xSite, evalResult);
	}

	inline void search_eval_with_grad_vec(double* xSite, double* evalResult, double* grad)
	{
		search_eval_with_grad_vec_at_array(0, xSite, evalResult, grad);
	}


	MatlabInterpEval(const mxArray* pp)
	{
#define GET_INT_STRUCT(var, varinstruct) var = (int) *(mxGetPr(mxGetField(pp, 0, #varinstruct))) 
		GET_INT_STRUCT(xDim, xDim);
		GET_INT_STRUCT(numVec, dim);
		GET_INT_STRUCT(arrayOffset, arrayOffset);
		GET_INT_STRUCT(fullVecEvalCoefsLength, fullVecEvalCoefsLength);
		GET_INT_STRUCT(singleVecEvalCoefsLength, fullVecEvalCoefsLength);
#undef GET_INT_STRUCT 

#define GET_DBL_PTR_STRUCT(var, varinstruct) var = mxGetPr(mxGetField(pp, 0, #varinstruct))
		GET_DBL_PTR_STRUCT(coefs, coefs);
#undef GET_DBL_PTR_STRUCT

#define GET_LOCAL_INT_ARRAY_STRUCT(var, varinstruct) memcpy(var, mxGetData(mxGetField(pp,0,#varinstruct)), sizeof(int)*xDim)
		GET_LOCAL_INT_ARRAY_STRUCT(order, order);
		GET_LOCAL_INT_ARRAY_STRUCT(pieces, pieces);
		GET_LOCAL_INT_ARRAY_STRUCT(xPts, xPts);
#undef GET_LOCAL_INT_ARRAY_STRUCT

		const mxArray* ppBreaks = mxGetField(pp, 0, "breaks");
		for (int i_dim = 0; i_dim < xDim; i_dim++)
		{
			xGrid[i_dim] = mxGetPr(mxGetCell(ppBreaks, i_dim));
		}
	}
};

