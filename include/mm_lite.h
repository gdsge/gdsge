#pragma once 

// #define MM_SAFE

template <class T> class Vector {
public:
	int pts;
	T* data;
	Vector(int _pts) :pts(_pts) {}
	Vector(int _pts, T* _data) :pts(_pts), data(_data - 1) {}
	Vector(T* _data) : pts(0), data(_data - 1) {}
	T& operator()  (int idx)
	{
		// 1 based
		return data[idx];
	}

	const T& operator() (int idx) const
	{
		// 1 based
		return data[idx];
	}
};

// Column major matrix
template <class T> class Matrix {
public:
	int stride;
	T* data;
	Matrix(T* data_, int stride_) : stride(stride_), data(data_) {}

	T& operator() (int m, int n)
	{
		// 1 based 
		return data[(m - 1) + (n - 1)*stride];
	}

	const T& operator() (int m, int n) const
	{
		// 1 based
		return data[(m - 1) + (n - 1)*stride];
	}
};

#ifdef MM_SAFE
#define GET_DMAT0_VIEW(var) const mxArray* __##var = mexGetVariablePtr("caller",#var); \
	if(__##var==0) mexErrMsgTxt("Variable doesn't exist: "#var); \
	if (!mxIsDouble(__##var)) mexErrMsgTxt("Not double: "#var); \
	double* _##var = mxGetPr(__##var)
#else
#define GET_DMAT0_VIEW(var) \
	const mxArray* __##var = mexGetVariablePtr("caller",#var); \
	double* _##var = mxGetPr(__##var)
#endif

#define GET_DMAT0_VIEW_FROM_MX(var, mx) const mxArray* __##var=mx; \
	if(__##var==0) mexErrMsgTxt("Variable doesn't exist: "#var); \
	if (!mxIsDouble(__##var)) mexErrMsgTxt("Not double: "#var); \
	double* _##var = mxGetPr(__##var)


#define GET_DMAT0_VIEW_RENAME(var,newvar) const mxArray* __##newvar = mexGetVariablePtr("caller",#var); \
	if(__##newvar==0) mexErrMsgTxt("Variable doesn't exist: "#var); \
	if (!mxIsDouble(__##newvar)) mexErrMsgTxt("Not double: "#var); \
	double* _##newvar = mxGetPr(__##newvar)

#define GET_DV_VIEW(var) GET_DMAT0_VIEW(var); \
	Vector<double> var(_##var)

#define GET_DV_VIEW_RENAME(var,newvar) GET_DMAT0_VIEW_RENAME(var,newvar); \
	Vector<double> newvar(_##newvar)

#define GET_DM_VIEW(var) GET_DMAT0_VIEW(var); \
	Matrix<double> var(_##var,*(mxGetDimensions(__##var)));

#define GET_DV_VIEW_FROM_MX(var, mx) GET_DMAT0_VIEW_FROM_MX(var, mx); \
	Vector<double> var(_##var)

#define GET_DM_VIEW_FROM_MX(var, mx) GET_DMAT0_VIEW_FROM_MX(var, mx); \
	Matrix<double> var(_##var,*(mxGetDimensions(__##var)));

#define GET_MX_ARRAY(var) mxArray* __##var = mexGetVariable("caller",#var);

#ifdef MM_SAFE
#define GET_INT(var) mxArray* __##var = mexGetVariable("caller",#var); \
	if(__##var==0) mexErrMsgTxt("Variable doesn't exist: "#var); \
	if (!mxIsDouble(__##var)) mexErrMsgTxt("Not double: "#var); \
	int var = (int)*mxGetPr(__##var);
#else
#define GET_INT(var) mxArray* __##var = mexGetVariable("caller",#var); \
	int var = (int)*mxGetPr(__##var);
#endif

#ifdef MM_SAFE
#define GET_DBL(var) mxArray* __##var = mexGetVariable("caller",#var); \
	if (__##var == 0) mexErrMsgTxt("Variable doesn't exist: "#var); \
	if (!mxIsDouble(__##var)) mexErrMsgTxt("Not double: "#var); \
	double var = *mxGetPr(__##var);
#else
#define GET_DBL(var) mxArray* __##var = mexGetVariable("caller",#var); \
	double var = *mxGetPr(__##var);
#endif

#define GET_STRUCT(structname)\
	const mxArray* structname = mexGetVariablePtr("caller", #structname); \
	if (structname==0) mexErrMsgTxt("Struct doesn't exist: "#structname); \
	if (!mxIsStruct(structname)) mexErrMsgTxt("Not struct: "#structname);

#define GET_STRUCT_VAR(structname,var) \
	const mxArray* __##structname##var##_struct = mexGetVariablePtr("caller", #structname); \
	if (__##structname##var##_struct==0) mexErrMsgTxt("Struct doesn't exist: "#structname); \
	if (!mxIsStruct(__##structname##var##_struct)) mexErrMsgTxt("Not struct: "#structname); \
	const mxArray* __##structname##var##_var = mxGetField(__##structname##var##_struct, 0, #var); \
	if (__##structname##var##_var==0) mexErrMsgTxt("Var "#var" doesn't exist in struct "#structname);

#define GET_DV_VIEW_STRUCT(structname, var) \
	GET_STRUCT_VAR(structname,var); \
	if (!mxIsDouble(__##structname##var##_var)) mexErrMsgTxt("Not double: "#var); \
	double* _##structname##var = mxGetPr(__##structname##var##_var); \
	Vector<double> ##structname##_##var(_##structname##var)

#define GET_DM_VIEW_STRUCT(structname, var) \
	GET_STRUCT_VAR(structname,var); \
	if (!mxIsDouble(__##structname##var##_var)) mexErrMsgTxt("Not double: "#var); \
	double* _##structname##var = mxGetPr(__##structname##var##_var); \
	Matrix<double> ##structname##_##var(_##structname##var, *(mxGetDimensions(__##structname##var##_var)));
