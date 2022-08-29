// This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE

#ifndef _HAS_CONDITIONAL_EXPLICIT
#if defined(__CUDACC__)
#define _HAS_CONDITIONAL_EXPLICIT 0 // TRANSITION
#elif defined(__INTEL_COMPILER)
#define _HAS_CONDITIONAL_EXPLICIT 0
#elif defined(__EDG__)
#define _HAS_CONDITIONAL_EXPLICIT 1
#elif defined(__clang__)
#define _HAS_CONDITIONAL_EXPLICIT 0 // TRANSITION, Clang 9
#else // vvv C1XX vvv
#define _HAS_CONDITIONAL_EXPLICIT 1
#endif // ^^^ C1XX ^^^
#endif // _HAS_CONDITIONAL_EXPLICIT

#define ADEPT_REMOVE_NULL_STATEMENTS 1

#define MAX_STACK_DIM 10000

#include "mex.h"
#include "mm_lite.h"
#include "interp_lite.h"
#include "codosol.h"
#include "MatlabPp.h"
#include "MatlabInterpEval.h"
#include "cmath"

#include<vector>
using std::vector;

#ifndef NO_OMP
#include "omp.h"
#endif

#include "adept_extension.h"
#include "wl_math.h"
#ifdef USE_ASG
  #include "class_handle.hpp"
  #include "asg_adouble.h"
  using namespace AdaptiveSparseGridInterp;
#endif

GDSGE_OTHER_INCLUDE

using adept::adouble;
using adept::Stack;

#ifdef HAS_INIT
void task_init(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
#endif

void task_inf_horizon(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  GET_INT(MEX_TASK_NAME);
  GET_INT(MEX_TASK_INIT);
  GET_INT(MEX_TASK_INF_HORIZON);

#ifdef HAS_INIT  
  if (MEX_TASK_NAME==MEX_TASK_INIT)
  {
    task_init(nlhs, plhs, nrhs, prhs);
    return;
  }
#endif
  
  if (MEX_TASK_NAME==MEX_TASK_INF_HORIZON)
  {
    task_inf_horizon(nlhs, plhs, nrhs, prhs);
    return;
  }
  
  mexErrMsgTxt("No mex task executed.");
}

#ifdef HAS_INIT
TASK_INIT_CODE
#endif

TASK_INF_HORIZON_CODE
