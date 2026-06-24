        // SymPy analytic-Jacobian model function (all double; value + JAC).
        #define JAC(i_eq,i_var) GDSGE_jac[(i_eq) + NUM_EQUATIONS*(i_var)]
        auto GDSGE_FUNC_MODEL_NUMBER_double = [&](double* GDSGE_x, double* GDSGE_f, double* GDSGE_jac=0, int GDSGE_EVAL=0)
        {
            ARG_CODE

            MODEL_BODY_CODE;

            AUX_ASSIGN_CODE;

            // mirror residuals into GDSGE_eqval (autodiff parity)
            for (int GDSGE_ieq=0; GDSGE_ieq<NUM_EQUATIONS; ++GDSGE_ieq) GDSGE_eqval[GDSGE_ieq]=GDSGE_f[GDSGE_ieq];
        };
        #undef JAC
        auto GDSGE_OBJ_MODEL_NUMBER = [&](double* GDSGE_x, double* GDSGE_f, double* GDSGE_jac)
        {
            GDSGE_FUNC_MODEL_NUMBER_double(GDSGE_x, GDSGE_f, GDSGE_jac);
        };
