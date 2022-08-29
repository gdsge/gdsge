GDSGE_NPROB = size(GDSGE_evalGrids,2);
GDSGE_SIZE = [1,GDSGE_NPROB];

% Solve the problem
POLICY_INBOUND_CODE

% Interp the warmup using last period solution
GDSGE_SOL0 = (GDSGE_LB + GDSGE_UB)/2;
GDSGE_SOL = GDSGE_SOL0;
if GDSGE_SOL_ASG_INTERP.get_current_level>=0
    % Interp to get the warmup
    GDSGE_SOL0 = GDSGE_SOL_ASG_INTERP.eval_vec(GDSGE_evalArrayIdx,GDSGE_evalGrids);
    GDSGE_SOL = GDSGE_SOL0;
    if UseAdaptiveBound==1
        INBOUND_ADAPTIVE_TIGHT_CODE
    end
end

GDSGE_EQVAL = 1e20*ones(GDSGE_MAXDIM,GDSGE_NPROB);
GDSGE_F = 1e20*ones(1,GDSGE_NPROB);
GDSGE_AUX = zeros(GDSGE_NUM_AUX,GDSGE_NPROB);
GDSGE_SKIP = zeros(1,GDSGE_NPROB);
GDSGE_DATA = zeros(GDSGE_MAXDATA,GDSGE_NPROB);

GDSGE_DATA(:) = [repmat([shock_num;PARAMS_SEMI_COLON;shock_trans(:);SHOCK_GRID_SEMI_COLON],1,GDSGE_NPROB);GDSGE_evalArrayIdx;GDSGE_evalGrids];

SOLVE_PROBLEM_CODE

POST_SOL_CODE

POLICY_ASSIGN_CODE

AUX_ASSIGN_CODE