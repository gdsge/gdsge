// This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE
// SymPy backend: same CoDoSol drive as call_fmin.tpl.cpp, but the model
// function is all-double (x,f,jac) — no adept eval flag. Value-only eval omits
// the jac pointer; the 6th-output cross-check fills it via the OBJ functor.
if (MODEL_CONDITION) {
  if (GDSGE_DEBUG_EVAL_ONLY==0)
  {
    // Warm start already in GDSGE_SOL (placed by the caller); solve in place.
    GDSGE_F(GDSGE_I) = CoDoSol::solve(NUM_EQUATIONS, &GDSGE_SOL(1, GDSGE_I), &GDSGE_LB(1, GDSGE_I), &GDSGE_UB(1, GDSGE_I), UseBroyden, TolFun, 0, SolMaxIter, GDSGE_OBJ_MODEL_NUMBER, &GDSGE_OPT_INFO(GDSGE_I));
  }

  // Evaluate at the solution
  double* GDSGE_x = &GDSGE_SOL(1, GDSGE_I);
  #if MAXDIM>MAX_STACK_DIM
  vector<double> GDSGE_EQ(NUM_EQUATIONS);
  #else
  double GDSGE_EQ[NUM_EQUATIONS];
  #endif

  if (GDSGE_DEBUG_EVAL_ONLY==2 && GDSGE_JAC_OUT)
  {
    GDSGE_OBJ_MODEL_NUMBER(&GDSGE_x[0], &GDSGE_EQ[0], &GDSGE_JAC_OUT[(GDSGE_I-1)*NUM_EQUATIONS*NUM_EQUATIONS]);
  }
  else
  {
    GDSGE_FUNC_MODEL_NUMBER_double(&GDSGE_x[0], &GDSGE_EQ[0]);
  }
}
