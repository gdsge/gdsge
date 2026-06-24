int GDSGE_INTERP_CELL[GDSGE_INTERP_XDIM] = { 0 };
adouble GDSGE_INTERP_RSLT_adouble[GDSGE_NUM_INTERP] = {0};
double GDSGE_INTERP_RSLT_double[GDSGE_NUM_INTERP] = {0};
bool INTERP_VEC_FLAG[GDSGE_NUM_INTERP] = {false};
double GDSGE_INTERP_GRAD[GDSGE_NUM_INTERP*GDSGE_INTERP_XDIM] = { 0 };

// Vector evaluation
auto GDSGE_INTERP_VEC_adouble = [&GDSGE_CSPLINE_VEC,&GDSGE_INTERP_RSLT_adouble,&GDSGE_INTERP_GRAD,&GDSGE_INTERP_CELL,&GDSGE_EVAL_GRAD_FLAG] (int shockIdx, ADOUBLE_VAR_NAME)
{
  // Search once
  adouble xSite[] = {VAR_NAME};

  GDSGE_CSPLINE_VEC.search_eval_vec_at_array_adouble(shockIdx - 1, xSite, GDSGE_INTERP_RSLT_adouble, GDSGE_INTERP_GRAD, GDSGE_INTERP_CELL, GDSGE_EVAL_GRAD_FLAG);
};

// Vector evaluation
auto GDSGE_INTERP_VEC_double = [&GDSGE_CSPLINE_VEC,&GDSGE_INTERP_RSLT_double,&GDSGE_INTERP_CELL,&GDSGE_EVAL_GRAD_FLAG] (int shockIdx, DOUBLE_VAR_NAME)
{
  // Search once
  double xSite[] = {VAR_NAME};

  GDSGE_CSPLINE_VEC.search_eval_vec_at_array(shockIdx - 1, xSite, GDSGE_INTERP_RSLT_double, GDSGE_INTERP_CELL);
};
