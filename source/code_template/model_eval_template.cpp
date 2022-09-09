        auto GDSGE_OBJ_MODEL_NUMBER = [&](double* GDSGE_x, double* GDSGE_f, double* GDSGE_jac)
		{
            if (GDSGE_jac == 0)
            {
                GDSGE_FUNC_MODEL_NUMBER_double(&GDSGE_x[0], &GDSGE_f[0]);
            }
            else
            {
                #ifdef USE_FINITE_DIFF
                    GDSGE_FUNC_MODEL_NUMBER_double(&GDSGE_x[0], &GDSGE_f[0]);
                    #if MAXDIM>MAX_STACK_DIM
                    vector<double> GDSGE_x_new(NUM_EQUATIONS);
                    vector<double> GDSGE_f_new(NUM_EQUATIONS);
                    #else
                    double GDSGE_x_new[NUM_EQUATIONS];
                    double GDSGE_f_new[NUM_EQUATIONS];
                    #endif
                    memcpy(&GDSGE_x_new[0], &GDSGE_x[0], sizeof(double)*NUM_EQUATIONS);
                    for (int i_var = 0; i_var < NUM_EQUATIONS; ++i_var)
                    {
                        // Perturb
                        GDSGE_x_new[i_var] = GDSGE_x[i_var] + FiniteDiffDelta;
                        GDSGE_FUNC_MODEL_NUMBER_double(&GDSGE_x_new[0], &GDSGE_f_new[0]);
                        
                        // Jacobian, eq goes first in memory
                        #define JAC(i_eq,i_var) GDSGE_jac[(i_eq) + NUM_EQUATIONS*(i_var)] 
                        for (int i_eq = 0; i_eq < NUM_EQUATIONS; ++i_eq)
                        {
                            JAC(i_eq,i_var) = (GDSGE_f_new[i_eq] - GDSGE_f[i_eq]) / FiniteDiffDelta;
                        }
                        #undef JAC
                        
                        // Revert
                        GDSGE_x_new[i_var] = GDSGE_x[i_var];
                    }
                #else
                    // Copy to adouble
                    #if MAXDIM>MAX_STACK_DIM
                    vector<adouble> GDSGE_x_adept(NUM_EQUATIONS);
                    vector<adouble> GDSGE_EQ(NUM_EQUATIONS);
                    #else
                    adouble GDSGE_x_adept[NUM_EQUATIONS];
                    adouble GDSGE_EQ[NUM_EQUATIONS];
                    #endif
                    for (int j = 0; j < NUM_EQUATIONS; ++j)
                    {
                        GDSGE_x_adept[j] = GDSGE_x[j];
                    }
                    _stack.new_recording();
                    _stack.set_max_jacobian_threads(1);
                    GDSGE_EVAL_GRAD_FLAG = 1;

                    // Evaluate functions
                    GDSGE_FUNC_MODEL_NUMBER_adouble(&GDSGE_x_adept[0], &GDSGE_EQ[0]);
                    if (GDSGE_f) {
                        for (int i = 0; i < NUM_EQUATIONS ; i++)
                        {
                            GDSGE_f[i] = GDSGE_EQ[i].value();
                        }
                    }
                    // Calculate jacobian
                    
                    PRE_JAC_CODE
                    _stack.independent(&GDSGE_x_adept[0], NUM_EQUATIONS);
                    _stack.dependent(&GDSGE_EQ[0], NUM_EQUATIONS);
                    _stack.jacobian_forward_vec(GDSGE_jac);
                    POST_JAC_CODE
                #endif
            }
		};