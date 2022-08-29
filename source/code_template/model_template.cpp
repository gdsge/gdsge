		auto GDSGE_FUNC_MODEL_NUMBER_adouble = [&](adouble* GDSGE_x, adouble* GDSGE_EQ, int GDSGE_EVAL=0)
		{
			ARG_CODE;

			DECLARE_CODE;

			MODEL_CODE;

			HEADER_AUX_ASSIGN_CODE;

			EQUATION_CODE;

			for (int i = 0; i < NUM_EQUATIONS; ++i)
			{
				GDSGE_eqval[i] = value(GDSGE_EQ[i]);
			}
		};