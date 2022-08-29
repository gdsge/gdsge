// This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE
if (MODEL_CONDITION) {
  if (GDSGE_DEBUG_EVAL_ONLY==0)
  {
    // Copy sol0 to sol
    memcpy(&GDSGE_SOL(1,GDSGE_I), &GDSGE_SOL0(1,GDSGE_I), sizeof(double)*NUM_EQUATIONS);

    GDSGE_F(GDSGE_I) = CoDoSol::solve(NUM_EQUATIONS, &GDSGE_SOL(1, GDSGE_I), &GDSGE_LB(1, GDSGE_I), &GDSGE_UB(1, GDSGE_I), UseBroyden, TolFun, 0, SolMaxIter, GDSGE_OBJ_MODEL_NUMBER, &GDSGE_OPT_INFO(GDSGE_I));
  }
  else if (GDSGE_DEBUG_EVAL_ONLY==1)
  {
    memcpy(&GDSGE_SOL(1,GDSGE_I), &GDSGE_SOL0(1,GDSGE_I), sizeof(double)*NUM_EQUATIONS);
  }
  
  // Evaluate at the solution
  double* GDSGE_x = &GDSGE_SOL(1, GDSGE_I);
  #if MAXDIM>MAX_STACK_DIM
  vector<double> GDSGE_EQ(NUM_EQUATIONS);
  #else
  double GDSGE_EQ[NUM_EQUATIONS];
  #endif
  
  GDSGE_FUNC_MODEL_NUMBER_double(&GDSGE_x[0], &GDSGE_EQ[0], 1);
}