GDSGE_SKIP(:) = 0;
[GDSGE_SOL,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL,GDSGE_OPT_INFO] = mex_MODEL_NAME(GDSGE_SOL,GDSGE_LB,GDSGE_UB,GDSGE_DATA,GDSGE_SKIP,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL);

% Randomzie for nonconvert point
GDSGE_MinorIter = 0;
while ((max(isnan(GDSGE_F)) || max(GDSGE_F(:))>TolSol) && GDSGE_MinorIter<MaxMinorIter)
    % Use randomize as initial guess
    GDSGE_X0Rand = rand(size(GDSGE_SOL)) .* (GDSGE_UB-GDSGE_LB) + GDSGE_LB;
    NeedResolved = (GDSGE_F>TolSol) | isnan(GDSGE_F);
    GDSGE_SOL(:,NeedResolved) = GDSGE_X0Rand(:,NeedResolved);
    GDSGE_SKIP(:) = 0;
    GDSGE_SKIP(~NeedResolved) = 1;
    
    [GDSGE_SOL,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL,GDSGE_OPT_INFO] = mex_MODEL_NAME(GDSGE_SOL,GDSGE_LB,GDSGE_UB,GDSGE_DATA,GDSGE_SKIP,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL);
    
    GDSGE_MinorIter = GDSGE_MinorIter+1;
end