if ~SkipModelInit
    MEX_TASK_NAME = MEX_TASK_INIT;
    
    %% Construct interpolation
    GDSGE_ASG_INTERP = asg({STATE_COMMA},GDSGE_NUM_INTERP,shock_num);
    
    while GDSGE_ASG_INTERP.get_current_level < AsgMaxLevel
        if GDSGE_ASG_INTERP.get_current_level<AsgMinLevel
            [GDSGE_evalArrayIdx,GDSGE_evalGrids,GDSGE_evalGridsLength] = GDSGE_ASG_INTERP.get_eval_grids(0.0);
        else
            [GDSGE_evalArrayIdx,GDSGE_evalGrids,GDSGE_evalGridsLength] = GDSGE_ASG_INTERP.get_eval_grids(AsgThreshold);
        end
        if isempty(GDSGE_evalGrids)
            break;
        end
        
        GDSGE_NPROB = size(GDSGE_evalGrids,2);
        GDSGE_SIZE = [1,GDSGE_NPROB];
        
        % Solve the problem
        POLICY_INIT_INBOUND_CODE
        
        POLICY_INIT_INITIALIZE_CODE
        
        GDSGE_EQVAL = 1e20*ones(GDSGE_MAXDIM_INIT,GDSGE_NPROB);
        GDSGE_F = 1e20*ones(1,GDSGE_NPROB);
        GDSGE_SOL = zeros(GDSGE_MAXDIM_INIT,GDSGE_NPROB);
        GDSGE_SOL(:) = GDSGE_X0;
        GDSGE_AUX = zeros(GDSGE_NUM_AUX_INIT,GDSGE_NPROB);
        GDSGE_SKIP = zeros(1,GDSGE_NPROB);
        GDSGE_DATA = zeros(GDSGE_MAXDATA,GDSGE_NPROB);
        
        GDSGE_DATA(:) = [repmat([shock_num;PARAMS_SEMI_COLON;shock_trans(:);SHOCK_GRID_SEMI_COLON],1,GDSGE_NPROB);GDSGE_evalArrayIdx;GDSGE_evalGrids];
        
        SOLVE_PROBLEM_CODE
        
        POLICY_INIT_ASSIGN_CODE
        
        AUX_INIT_ASSIGN_CODE
        
        % Map variables to initial interp
        INTERP_INITIALIZE_CODE
        GDSGE_evalRslts = [INTERP_VAR_SEMI_COLON];
        GDSGE_ASG_INTERP.push_eval_results(GDSGE_evalRslts);
    end
end