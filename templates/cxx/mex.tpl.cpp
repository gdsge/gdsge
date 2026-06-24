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
#include "adept_extension.h"
#include "mm_lite.h"
#include "codosol.h"
#include "interp_lite.h"
#include "MatlabPp.h"
#include "MatlabInterpEval.h"
#include "cmath"
#include <cstdint>

#include<vector>
using std::vector;

#ifndef NO_OMP
#include "omp.h"
#endif

#include "wl_math.h"
#ifdef USE_ASG
  #include "class_handle.hpp"
  #include "asg_adouble.h"
  using namespace AdaptiveSparseGridInterp;
#endif

GDSGE_OTHER_INCLUDE

using adept::adouble;
using adept::Stack;

// ---------------------------------------------------------------------------
// Counter-based RNG for in-MEX randomized restart (UseMexRandomize).
// Every draw is a pure function of (seed, salt, point, trial, component), so
// results are independent of OpenMP thread count/scheduling and resumable
// across batches. splitmix64 mixing -> uniform double in [0,1).
// ---------------------------------------------------------------------------
static inline uint64_t gdsge_splitmix64(uint64_t x)
{
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

static inline double gdsge_rng_uniform(uint64_t seed, uint64_t salt, uint64_t point,
                                       uint64_t trial, uint64_t comp)
{
  uint64_t s = seed;
  s = gdsge_splitmix64(s ^ (salt  * 0x9E3779B97F4A7C15ULL + 0x1ULL));
  s = gdsge_splitmix64(s ^ (point * 0xD1B54A32D192ED03ULL + 0x2ULL));
  s = gdsge_splitmix64(s ^ (trial * 0x9E3779B97F4A7C15ULL + 0x5ULL));
  s = gdsge_splitmix64(s ^ (comp  * 0xD1B54A32D192ED03ULL + 0x7ULL));
  // 53-bit mantissa -> [0,1)
  return (double)(s >> 11) * (1.0 / 9007199254740992.0);
}

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
