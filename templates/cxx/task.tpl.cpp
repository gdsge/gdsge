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

  #ifndef NO_OMP
  omp_set_num_threads(NumThreads);
  #endif

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

  // Optional analytic-Jacobian debug output (Phase 8 cross-check): per problem,
  // NUM_EQUATIONS*NUM_EQUATIONS column-major, filled only when nlhs>=6 and
  // GDSGE_DEBUG_EVAL_ONLY==2. Existing 5-output callers are unaffected.
  double* GDSGE_JAC_OUT = 0;
  if (nlhs >= 6) {
    plhs[5] = mxCreateDoubleMatrix(NUM_EQUATIONS*NUM_EQUATIONS, GDSGE_NPROB, mxREAL);
    GDSGE_JAC_OUT = mxGetPr(plhs[5]);
  }

  // Input and output
  INTERP_GET_CODE;

  // Optional in-MEX nearby-warmup resolve. Read the per-dimension problem
  // strides from the caller workspace; absent/empty => no resolve (Pass 0 only),
  // which keeps the init task, ASG, and the MATLAB-sweep fallback on today's path.
  int GDSGE_NUM_DIMS = 0;
  double* GDSGE_STRIDES = 0;
  {
    const mxArray* GDSGE_STRIDES_MX = mexGetVariablePtr("caller", "GDSGE_PROBLEM_STRIDES");
    if (GDSGE_STRIDES_MX != 0 && mxIsDouble(GDSGE_STRIDES_MX) && mxGetNumberOfElements(GDSGE_STRIDES_MX) > 0) {
      GDSGE_NUM_DIMS = (int) mxGetNumberOfElements(GDSGE_STRIDES_MX);
      GDSGE_STRIDES = mxGetPr(GDSGE_STRIDES_MX);
    }
  }

  // Optional in-MEX randomized restart (UseMexRandomize). Absent/<=0 => today's
  // behavior (no in-MEX randomize; MATLAB drives restarts). Read from caller.
  int      GDSGE_RAND_BATCH = 0;
  uint64_t GDSGE_RAND_SEED = 0, GDSGE_RAND_SALT = 0, GDSGE_RAND_TRIAL_OFFSET = 0;
  {
    const mxArray* p = mexGetVariablePtr("caller", "GDSGE_RANDOMIZE_BATCH");
    if (p != 0 && mxIsNumeric(p) && mxGetNumberOfElements(p) > 0)
      GDSGE_RAND_BATCH = (int) mxGetScalar(p);
    if (GDSGE_RAND_BATCH > 0) {
      const mxArray* ps = mexGetVariablePtr("caller", "GDSGE_RANDOM_SEED");
      if (ps != 0 && mxGetNumberOfElements(ps) > 0) GDSGE_RAND_SEED = (uint64_t) mxGetScalar(ps);
      const mxArray* pa = mexGetVariablePtr("caller", "GDSGE_RANDOMIZE_SALT");
      if (pa != 0 && mxGetNumberOfElements(pa) > 0) GDSGE_RAND_SALT = (uint64_t) mxGetScalar(pa);
      const mxArray* pt = mexGetVariablePtr("caller", "GDSGE_RANDOMIZE_TRIAL_OFFSET");
      if (pt != 0 && mxGetNumberOfElements(pt) > 0) GDSGE_RAND_TRIAL_OFFSET = (uint64_t) mxGetScalar(pt);
    }
  }

  // Per-grid-point solve, reusable across the initial pass and each resolve pass.
  // CONTRACT: the warm-start guess must already be in GDSGE_SOL(:,GDSGE_I) before
  // this runs. Pass 0 copies it from GDSGE_SOL0; a resolve pass copies it from a
  // converged neighbor. The solve/eval operate in place on GDSGE_SOL.
  auto GDSGE_solve_one = [&](int GDSGE_I)
  {
    START_LOOP_CODE;

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

    FINISH_LOOP_CODE;
#undef POPN
#undef POPNARRAY
#undef POPNARRAY_BASE1
  };

  // Pass 0: solve every non-skipped point from its supplied guess.
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
    memcpy(&GDSGE_SOL(1,GDSGE_I), &GDSGE_SOL0(1,GDSGE_I), sizeof(double)*NUM_EQUATIONS);
    GDSGE_solve_one(GDSGE_I);
  }

  // In-MEX nearby-warmup resolve (deterministic snapshot semantics; bit-exact
  // with solveProblems.m's MATLAB sweep). Reusable: called by the standalone
  // resolve sweep and by each batched randomize trial below.
  auto GDSGE_neighbor_fixpoint = [&]()
  {
    if (GDSGE_NUM_DIMS <= 0) return;
    vector<char> GDSGE_SOLVED0(GDSGE_NPROB + 1);
    auto GDSGE_count_unconv = [&]() -> int {
      int c = 0;
      for (int i = 1; i <= GDSGE_NPROB; i++) if (!(GDSGE_F(i) <= TolSol)) c++;
      return c;
    };
    int GDSGE_BEFORE = GDSGE_count_unconv();
    int GDSGE_AFTER = 0x7fffffff;
    while (GDSGE_AFTER != GDSGE_BEFORE)
    {
      GDSGE_BEFORE = GDSGE_count_unconv();
      for (int GDSGE_DIM = 0; GDSGE_DIM < GDSGE_NUM_DIMS; GDSGE_DIM++)
      {
        int GDSGE_STRIDE = (int) GDSGE_STRIDES[GDSGE_DIM];

        // Lower neighbor: snapshot the converged set, then fused copy + solve.
        for (int i = 1; i <= GDSGE_NPROB; i++)
          GDSGE_SOLVED0[i] = (GDSGE_F(i) <= TolSol) ? 1 : 0;
#pragma omp parallel for schedule(dynamic)
        for (int GDSGE_I = 1; GDSGE_I <= GDSGE_NPROB; GDSGE_I++)
        {
          int GDSGE_J = GDSGE_I - GDSGE_STRIDE;
          if (!GDSGE_SOLVED0[GDSGE_I] && GDSGE_J >= 1 && GDSGE_SOLVED0[GDSGE_J])
          {
            memcpy(&GDSGE_SOL(1,GDSGE_I), &GDSGE_SOL(1,GDSGE_J), sizeof(double)*NUM_EQUATIONS);
            GDSGE_solve_one(GDSGE_I);
          }
        }

        // Upper neighbor.
        for (int i = 1; i <= GDSGE_NPROB; i++)
          GDSGE_SOLVED0[i] = (GDSGE_F(i) <= TolSol) ? 1 : 0;
#pragma omp parallel for schedule(dynamic)
        for (int GDSGE_I = 1; GDSGE_I <= GDSGE_NPROB; GDSGE_I++)
        {
          int GDSGE_J = GDSGE_I + GDSGE_STRIDE;
          if (!GDSGE_SOLVED0[GDSGE_I] && GDSGE_J <= GDSGE_NPROB && GDSGE_SOLVED0[GDSGE_J])
          {
            memcpy(&GDSGE_SOL(1,GDSGE_I), &GDSGE_SOL(1,GDSGE_J), sizeof(double)*NUM_EQUATIONS);
            GDSGE_solve_one(GDSGE_I);
          }
        }
      }
      GDSGE_AFTER = GDSGE_count_unconv();
    }
  };

  // Dispatch: batched in-MEX randomize (UseMexRandomize) takes precedence; else
  // the standalone neighbor sweep (UseMexResolve); else nothing (today's Pass-0-only).
  if (GDSGE_RAND_BATCH > 0)
  {
    auto GDSGE_count_left = [&]() -> int {
      int c = 0;
      for (int i = 1; i <= GDSGE_NPROB; i++) if (!(GDSGE_F(i) <= TolSol)) c++;
      return c;
    };
    // Gated neighbor sweep, mirroring solveProblems.m's numNeedResolvedAfter gate:
    // run the (expensive, full-grid) neighbor fixpoint only when the unconverged
    // count changed since the last sweep — i.e. when a random restart actually
    // converged new points to propagate. Running it on every trial (even when a
    // restart made no progress) launches thousands of no-op OpenMP grid sweeps and
    // is the dominant overhead on restart-heavy models. Skipping a no-op sweep is
    // result-preserving (the converged set is unchanged), so this stays bit-exact.
    int GDSGE_NAFTER = -1;   // -1 forces the first sweep
    int GDSGE_TRIALS_USED = 0;   // random-restart rounds actually run this call
    for (int GDSGE_TRIAL = 0; GDSGE_TRIAL < GDSGE_RAND_BATCH; GDSGE_TRIAL++)
    {
      int GDSGE_NLEFT = GDSGE_count_left();
      if (GDSGE_NLEFT == 0) break;

      if (GDSGE_NAFTER != GDSGE_NLEFT)
      {
        GDSGE_neighbor_fixpoint();          // propagate newly converged points
        GDSGE_NAFTER = GDSGE_count_left();
        if (GDSGE_NAFTER == 0) break;
      }

      uint64_t GDSGE_T = GDSGE_RAND_TRIAL_OFFSET + (uint64_t) GDSGE_TRIAL;
#pragma omp parallel for schedule(dynamic)
      for (int GDSGE_I = 1; GDSGE_I <= GDSGE_NPROB; GDSGE_I++)
      {
        if (GDSGE_F(GDSGE_I) <= TolSol) continue;   // keep converged points
        for (int GDSGE_C = 1; GDSGE_C <= NUM_EQUATIONS; GDSGE_C++)
        {
          double u = gdsge_rng_uniform(GDSGE_RAND_SEED, GDSGE_RAND_SALT,
                                       (uint64_t)(GDSGE_I - 1), GDSGE_T, (uint64_t)(GDSGE_C - 1));
          GDSGE_SOL(GDSGE_C, GDSGE_I) =
              GDSGE_LB(GDSGE_C, GDSGE_I) + u * (GDSGE_UB(GDSGE_C, GDSGE_I) - GDSGE_LB(GDSGE_C, GDSGE_I));
        }
        GDSGE_solve_one(GDSGE_I);
      }
      GDSGE_TRIALS_USED++;
    }
    // Report the actual restart-round count back to the caller (solveProblems),
    // which pre-seeds this var to the batch size; it advances minorIter and the
    // resolve-round print by the real count, not the batch ceiling.
    {
      mxArray* GDSGE_TU = mxCreateDoubleScalar((double) GDSGE_TRIALS_USED);
      mexPutVariable("caller", "GDSGE_RANDOMIZE_TRIALS_USED", GDSGE_TU);
      mxDestroyArray(GDSGE_TU);
    }
  }
  else if (GDSGE_NUM_DIMS > 0)
  {
    GDSGE_neighbor_fixpoint();   // standalone UseMexResolve sweep (unchanged behavior)
  }
}
