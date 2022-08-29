GET_MX_ARRAY(GDSGE_ASG_HANDLE);
AsgInterpArray* P_GDSGE_CPP_ASG = convertMat2Ptr<AsgInterpArray>(__GDSGE_ASG_HANDLE);
AsgInterpArrayAdoubleEvaluator GDSGE_CPP_ASG(P_GDSGE_CPP_ASG);