// This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE
void TASK_NAME(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  GET_DBL(TolFun);
  GET_DBL(TolSol);
  GET_INT(SolMaxIter);
  GET_INT(NumThreads);
  GET_INT(GDSGE_DEBUG_EVAL_ONLY);
  GET_INT(UseBroyden);
  GET_DBL(FiniteDiffDelta);
  GET_INT(GDSGE_USE_BROYDEN_NOW);

  omp_set_num_threads(NumThreads);
  
  // 
  int GDSGE_NPROB = mxGetN(prhs[3]);
  int GDSGE_NUM_SOL0 = mxGetN(prhs[0]);

  // Input
  GET_DM_VIEW_FROM_MX(GDSGE_SOL0,prhs[0]);
  GET_DM_VIEW_FROM_MX(GDSGE_LB,prhs[1]);
  GET_DM_VIEW_FROM_MX(GDSGE_UB,prhs[2]);
  GET_DM_VIEW_FROM_MX(GDSGE_DATA,prhs[3]);
  GET_DV_VIEW_FROM_MX(GDSGE_SKIP,prhs[4]);
  GET_DV_VIEW_FROM_MX(GDSGE_F0,prhs[5]);
  GET_DM_VIEW_FROM_MX(GDSGE_AUX0,prhs[6]);
  GET_DM_VIEW_FROM_MX(GDSGE_EQVAL0,prhs[7]);
  
  // Output
  plhs[0] = mxCreateDoubleMatrix(NUM_EQUATIONS,GDSGE_NPROB,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,GDSGE_NPROB,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(NUM_AUX,GDSGE_NPROB,mxREAL);
  plhs[3] = mxCreateDoubleMatrix(NUM_EQUATIONS,GDSGE_NPROB,mxREAL);
  plhs[4] = mxCreateDoubleMatrix(1,GDSGE_NPROB,mxREAL);
  
  GET_DM_VIEW_FROM_MX(GDSGE_SOL,plhs[0]);
  GET_DV_VIEW_FROM_MX(GDSGE_F,plhs[1]);
  GET_DM_VIEW_FROM_MX(GDSGE_AUX,plhs[2]);
  GET_DM_VIEW_FROM_MX(GDSGE_EQVAL,plhs[3]);
  GET_DV_VIEW_FROM_MX(GDSGE_OPT_INFO,plhs[4]);
  
  // Input and output
  INTERP_GET_CODE;
  
#pragma omp parallel for schedule(dynamic)
	for (int GDSGE_I = 1; GDSGE_I <= GDSGE_NPROB; GDSGE_I++)
	{
    GDSGE_F(GDSGE_I) = GDSGE_F0(GDSGE_I);
    if (GDSGE_SKIP(GDSGE_I) == 1)
    { 
      memcpy(&GDSGE_SOL(1,GDSGE_I), &GDSGE_SOL0(1,GDSGE_I), sizeof(double)*NUM_EQUATIONS);
      memcpy(&GDSGE_EQVAL(1,GDSGE_I), &GDSGE_EQVAL0(1,GDSGE_I), sizeof(double)*NUM_EQUATIONS);
      memcpy(&GDSGE_AUX(1,GDSGE_I), &GDSGE_AUX0(1,GDSGE_I), sizeof(double)*NUM_AUX);
			continue;
    }
    
    //Initiate adept
    Stack _stack;
	int GDSGE_EVAL_GRAD_FLAG = 0;
  
    INTERP_GET_THREAD_CODE;

		double* GDSGE_data = &GDSGE_DATA(1, GDSGE_I);
		double* GDSGE_aux = &GDSGE_AUX(1, GDSGE_I);
		double* GDSGE_eqval = &GDSGE_EQVAL(1, GDSGE_I);

	int GDSGE_POP_I = 0;
#define POPN(var) double var = GDSGE_data[GDSGE_POP_I++]
#define POPNARRAY(var, length) double* var = &GDSGE_data[GDSGE_POP_I]; GDSGE_POP_I+=length
#define POPNARRAY_BASE1(var, length) double* var = &GDSGE_data[GDSGE_POP_I-1]; GDSGE_POP_I+=length

		POP_CODE;
        
        PRE_MODEL_CODE;
    
    MODEL_CODE;
    
    CALL_FMIN_CODE;
	}
}