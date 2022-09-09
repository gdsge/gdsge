GDSGE_F_MODEL_ID(:) = 1e20;
GDSGE_SKIP_MODEL_ID(:) = 0;
if exist('GDSGE_Iter','var')>0 && GDSGE_Iter>1
    GDSGE_USE_BROYDEN_NOW = GDSGE_USE_BROYDEN;
end
PRE_CALL_MEX
[GDSGE_SOL_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_AUX_MODEL_ID,GDSGE_EQVAL_MODEL_ID,GDSGE_OPT_INFO_MODEL_ID] = mex_MODEL_NAME(GDSGE_SOL_MODEL_ID,GDSGE_LB_MODEL_ID,GDSGE_UB_MODEL_ID,GDSGE_DATA_MODEL_ID,GDSGE_SKIP_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_AUX_MODEL_ID,GDSGE_EQVAL_MODEL_ID);
POST_CALL_MEX
NeedResolved = (GDSGE_F_MODEL_ID>TolSol) | isnan(GDSGE_F_MODEL_ID);
GDSGE_USE_BROYDEN_NOW = 0;
% Randomzie for nonconvert point
GDSGE_MinorIter = 0;
numNeedResolvedAfter = inf;
while ((max(isnan(GDSGE_F_MODEL_ID)) || max(GDSGE_F_MODEL_ID(:))>TolSol) && GDSGE_MinorIter<MaxMinorIter)
    % Repeatedly use the nearest point as initial guess
    NeedResolved = (GDSGE_F_MODEL_ID>TolSol) | isnan(GDSGE_F_MODEL_ID);
    numNeedResolved = sum(NeedResolved);
    while numNeedResolvedAfter ~= numNeedResolved
        NeedResolved = (GDSGE_F_MODEL_ID>TolSol) | isnan(GDSGE_F_MODEL_ID);
        numNeedResolved = sum(NeedResolved);
        % Use the nearest point as initial guess
        for i_dim = 1:length(GDSGE_SIZE_MODEL_ID)
            stride = prod(GDSGE_SIZE_MODEL_ID(1:i_dim-1));
            
            NeedResolved = (GDSGE_F_MODEL_ID>TolSol) | isnan(GDSGE_F_MODEL_ID);
            GDSGE_SKIP_MODEL_ID(:) = 1;
            for GDSGE_i = 1:numel(GDSGE_F_MODEL_ID)
                if NeedResolved(GDSGE_i) && GDSGE_i-stride>=1
                    GDSGE_idx = GDSGE_i-stride;
                    if ~NeedResolved(GDSGE_idx)
                        GDSGE_SOL_MODEL_ID(:,GDSGE_i) = GDSGE_SOL_MODEL_ID(:,GDSGE_idx);
                        GDSGE_SKIP_MODEL_ID(GDSGE_i) = 0;
                    end
                end
            end
            PRE_CALL_MEX
            [GDSGE_SOL_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_AUX_MODEL_ID,GDSGE_EQVAL_MODEL_ID,GDSGE_OPT_INFO_MODEL_ID] = mex_MODEL_NAME(GDSGE_SOL_MODEL_ID,GDSGE_LB_MODEL_ID,GDSGE_UB_MODEL_ID,GDSGE_DATA_MODEL_ID,GDSGE_SKIP_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_AUX_MODEL_ID,GDSGE_EQVAL_MODEL_ID);
            POST_CALL_MEX
            
            NeedResolved = (GDSGE_F_MODEL_ID>TolSol) | isnan(GDSGE_F_MODEL_ID);
            GDSGE_SKIP_MODEL_ID(:) = 1;
            for GDSGE_i = 1:numel(GDSGE_F_MODEL_ID)
                if NeedResolved(GDSGE_i) && GDSGE_i+stride<=numel(GDSGE_F_MODEL_ID)
                    GDSGE_idx = GDSGE_i+stride;
                    if ~NeedResolved(GDSGE_idx)
                        GDSGE_SOL_MODEL_ID(:,GDSGE_i) = GDSGE_SOL_MODEL_ID(:,GDSGE_idx);
                        GDSGE_SKIP_MODEL_ID(GDSGE_i) = 0;
                    end
                end
            end
            PRE_CALL_MEX
            [GDSGE_SOL_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_AUX_MODEL_ID,GDSGE_EQVAL_MODEL_ID,GDSGE_OPT_INFO_MODEL_ID] = mex_MODEL_NAME(GDSGE_SOL_MODEL_ID,GDSGE_LB_MODEL_ID,GDSGE_UB_MODEL_ID,GDSGE_DATA_MODEL_ID,GDSGE_SKIP_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_AUX_MODEL_ID,GDSGE_EQVAL_MODEL_ID);
            POST_CALL_MEX
        end
        
        NeedResolvedAfter = (GDSGE_F_MODEL_ID>TolSol) | isnan(GDSGE_F_MODEL_ID);
        numNeedResolvedAfter = sum(NeedResolvedAfter);
    end
    
    % Use randomize as initial guess
    GDSGE_X0Rand = rand(size(GDSGE_SOL_MODEL_ID)) .* (GDSGE_UB_MODEL_ID-GDSGE_LB_MODEL_ID) + GDSGE_LB_MODEL_ID;
    NeedResolved = (GDSGE_F_MODEL_ID>TolSol) | isnan(GDSGE_F_MODEL_ID);
    GDSGE_SOL_MODEL_ID(:,NeedResolved) = GDSGE_X0Rand(:,NeedResolved);
    GDSGE_SKIP_MODEL_ID(:) = 0;
    GDSGE_SKIP_MODEL_ID(~NeedResolved) = 1;
    
    PRE_CALL_MEX
    [GDSGE_SOL_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_AUX_MODEL_ID,GDSGE_EQVAL_MODEL_ID,GDSGE_OPT_INFO_MODEL_ID] = mex_MODEL_NAME(GDSGE_SOL_MODEL_ID,GDSGE_LB_MODEL_ID,GDSGE_UB_MODEL_ID,GDSGE_DATA_MODEL_ID,GDSGE_SKIP_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_AUX_MODEL_ID,GDSGE_EQVAL_MODEL_ID);
    POST_CALL_MEX
    
    if UseAdaptiveBoundInSol==1 && exist('GDSGE_Iter','var')>0
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
    
    GDSGE_MinorIter = GDSGE_MinorIter+1;
end