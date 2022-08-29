#pragma once 

#include "asg.h"
#include "mex.h"
#include "assert.h"

using namespace AdaptiveSparseGridInterp;

#define GET_INT_FROM_STRUCT(structname, var) \
	const mxArray* mx_##var = mxGetField(structname,0,#var); \
	if (!mx_##var || !mxIsDouble(mx_##var)) \
		mexErrMsgTxt("Not double: "#var); \
	int var = (int)(*mxGetPr(mx_##var));

#define GET_PR_FROM_STRUCT(structname, var) \
	const mxArray* mx_##var = mxGetField(structname,0,#var); \
	if (!mx_##var || !mxIsDouble(mx_##var)) \
		mexErrMsgTxt("Not double: "#var); \
	double* var = mxGetPr(mx_##var);


template<size_t SIZE, class T> inline size_t array_size(T (&arr)[SIZE]) {
    return SIZE;
}

AsgInterp convert_matlab_struct_to_asg(const mxArray* mx_interp)
{
	if (!mxIsStruct)
		mexErrMsgTxt("Not struct.");

	GET_INT_FROM_STRUCT(mx_interp, nDim);
	GET_INT_FROM_STRUCT(mx_interp, nVec);
	GET_INT_FROM_STRUCT(mx_interp, maxLevel);

	GET_PR_FROM_STRUCT(mx_interp, grid);
	GET_PR_FROM_STRUCT(mx_interp, surplus);

	AsgInterp interp(nDim, nVec);

	interp.currentLevel = maxLevel;

	int sizeGrid = mxGetN(mx_grid);

	double* p_grid = mxGetPr(mx_grid);
	double* p_surplus = mxGetPr(mx_surplus);

	for (int i = 0; i < sizeGrid ; i++)
	{
		inarray grid = {};
		outarray surplusAtGrid = {};

		memcpy(grid.data(), p_grid, sizeof(double)*nDim);
		memcpy(surplusAtGrid.data(), p_surplus, sizeof(double)*nVec);
		interp.surplus[grid] = surplusAtGrid;
		
		p_grid += nDim;
		p_surplus += nVec;
	}

	return interp;
}

mxArray* convert_asg_to_matlab_struct(AsgInterp& interp)
{
	const char* fieldNames[] = {
		"nDim",
		"nVec",
		"maxLevel",
		"grid",
		"surplus",
	};

	mxArray* mx_interp = mxCreateStructMatrix(1, 1, array_size(fieldNames), fieldNames);

	int nDim = interp.nDim;
	int nVec = interp.nVec;

	mxArray* mx_nDim = mxCreateDoubleScalar((double)nDim);
	mxArray* mx_nVec = mxCreateDoubleScalar((double)nVec);
	mxArray* mx_maxLevel = mxCreateDoubleScalar((double)interp.currentLevel);

	// Create output array
	mxArray* mx_grid = mxCreateDoubleMatrix(nDim, interp.surplus.size(), mxREAL);
	mxArray* mx_surplus = mxCreateDoubleMatrix(nVec, interp.surplus.size(), mxREAL);

	// Assign grid and surplus to output
	double* p_grid = mxGetPr(mx_grid);
	double* p_surplus = mxGetPr(mx_surplus);
	for (std::map<inarray, outarray>::iterator it = interp.surplus.begin(); it != interp.surplus.end(); ++it)
	{
		memcpy(p_grid, it->first.data(), nDim * sizeof(double));
		memcpy(p_surplus, it->second.data(), nVec*sizeof(double));

		p_grid += nDim;
		p_surplus += nVec;
	}

	// Assign to structure
	int ptr = 0;

	mxSetFieldByNumber(mx_interp, 0, ptr, mx_nDim);
	ptr++;

	mxSetFieldByNumber(mx_interp, 0, ptr, mx_nVec);
	ptr++;

	mxSetFieldByNumber(mx_interp, 0, ptr, mx_maxLevel);
	ptr++;

	mxSetFieldByNumber(mx_interp, 0, ptr, mx_grid);
	ptr++;

	mxSetFieldByNumber(mx_interp, 0, ptr, mx_surplus);
	ptr++;

	assert(ptr == array_size(fieldNames));

	return mx_interp;
}
