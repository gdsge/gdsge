GDSGE_SKIP(:) = 0;
PRE_CALL_MEX
[GDSGE_SOL,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL,GDSGE_OPT_INFO] = mex_MODEL_NAME(GDSGE_SOL,GDSGE_LB,GDSGE_UB,GDSGE_DATA,GDSGE_SKIP,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL);
POST_CALL_MEX

% Use current asg as warmup for non-convergent point
if GDSGE_ASG_INTERP_NEW.get_current_level>=0
    NeedResolved = (GDSGE_F>TolSol) | isnan(GDSGE_F);
    GDSGE_SOL0 = GDSGE_SOL_ASG_INTERP_NEW.eval_vec(GDSGE_evalArrayIdx,GDSGE_evalGrids);
    GDSGE_SOL(:,NeedResolved) = GDSGE_SOL0(:,NeedResolved);
    
    GDSGE_SKIP(:) = 0;
    GDSGE_SKIP(~NeedResolved) = 1;
    
    if UseAdaptiveBound==1
        INBOUND_ADAPTIVE_TIGHT_CODE
    end
    PRE_CALL_MEX
    [GDSGE_SOL,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL,GDSGE_OPT_INFO] = mex_MODEL_NAME(GDSGE_SOL,GDSGE_LB,GDSGE_UB,GDSGE_DATA,GDSGE_SKIP,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL);
    POST_CALL_MEX
end

% Restore bound using last-iteration GDSGE_SOL
if GDSGE_SOL_ASG_INTERP.get_current_level>=0
    NeedResolved = (GDSGE_F>TolSol) | isnan(GDSGE_F);
    
    % Interp to get the warmup
    GDSGE_SOL0 = GDSGE_SOL_ASG_INTERP.eval_vec(GDSGE_evalArrayIdx,GDSGE_evalGrids);
    GDSGE_SOL(:,NeedResolved) = GDSGE_SOL0(:,NeedResolved);
    
    if UseAdaptiveBound==1
        INBOUND_ADAPTIVE_TIGHT_CODE
    end
end

% Randomzie for nonconvert point
GDSGE_MinorIter = 0;
while ((max(isnan(GDSGE_F)) || max(GDSGE_F(:))>TolSol) && GDSGE_MinorIter<MaxMinorIter)
    % Use randomize as initial guess
    GDSGE_X0Rand = rand(size(GDSGE_SOL)) .* (GDSGE_UB-GDSGE_LB) + GDSGE_LB;
    NeedResolved = (GDSGE_F>TolSol) | isnan(GDSGE_F);
    GDSGE_SOL(:,NeedResolved) = GDSGE_X0Rand(:,NeedResolved);
    GDSGE_SKIP(:) = 0;
    GDSGE_SKIP(~NeedResolved) = 1;
    
    PRE_CALL_MEX
    [GDSGE_SOL,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL,GDSGE_OPT_INFO] = mex_MODEL_NAME(GDSGE_SOL,GDSGE_LB,GDSGE_UB,GDSGE_DATA,GDSGE_SKIP,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL);
    POST_CALL_MEX
    
    if UseAdaptiveBoundInSol==1
        % Tentatively adjust the bound
        GDSGE_LB_OLD = GDSGE_LB;
        GDSGE_UB_OLD = GDSGE_UB;
        
        
        INBOUND_ADAPTIVE_TIGHT_CODE
        
        % Hitting lower bound
        GDSGE_SOL_hitting_lower_bound = abs(GDSGE_SOL - GDSGE_LB_OLD) < 1e-8;
        GDSGE_SOL_hitting_upper_bound = abs(GDSGE_SOL - GDSGE_UB_OLD) < 1e-8;
        
        % Adjust for those hitting lower bound or upper bound
        GDSGE_LB(~GDSGE_SOL_hitting_lower_bound) = GDSGE_LB_OLD(~GDSGE_SOL_hitting_lower_bound);
        GDSGE_UB(~GDSGE_SOL_hitting_upper_bound) = GDSGE_UB_OLD(~GDSGE_SOL_hitting_upper_bound);
        
        GDSGE_MinorIter = GDSGE_MinorIter+1;
    end
end
GDSGE_solved = (GDSGE_F<TolSol) & ~isnan(GDSGE_F);