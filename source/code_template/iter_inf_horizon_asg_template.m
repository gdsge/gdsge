MEX_TASK_NAME = MEX_TASK_INF_HORIZON;

GDSGE_SOL_ASG_INTERP = asg({STATE_COMMA},GDSGE_NUM_SOL,shock_num);

GDSGE_Metric = 1;
GDSGE_Iter = 0;

if nargin>=1 && isfield(GDSGE_OPTIONS,'WarmUp')
    if isfield(GDSGE_OPTIONS.WarmUp,'asg_interp_struct')
        GDSGE_ASG_INTERP = asg.construct_from_struct(GDSGE_OPTIONS.WarmUp.asg_interp_struct);
    end
    if isfield(GDSGE_OPTIONS.WarmUp,'Iter')
        GDSGE_Iter = GDSGE_OPTIONS.WarmUp.Iter;
    end
    if isfield(GDSGE_OPTIONS.WarmUp,'sol_asg_interp_struct')
        GDSGE_SOL_ASG_INTERP = asg.construct_from_struct(GDSGE_OPTIONS.WarmUp.sol_asg_interp_struct);
    end
end

stopFlag = false;
tic;
while(~stopFlag)
    GDSGE_Iter = GDSGE_Iter+1;
    
    PRE_ITER_CODE
    
    % Construction interpolation
    GDSGE_ASG_HANDLE = GDSGE_ASG_INTERP.objectHandle;
    GDSGE_ASG_INTERP_NEW = asg({STATE_COMMA},GDSGE_NUM_INTERP,shock_num);
    GDSGE_SOL_ASG_INTERP_NEW = asg({STATE_COMMA},GDSGE_NUM_SOL,shock_num);
    GDSGE_ASG_STORE_evalArrayIdx = cell(0);
    GDSGE_ASG_STORE_evalGridsUnscaled = cell(0);
    GDSGE_ASG_STORE_output = cell(0);
    
    while GDSGE_ASG_INTERP_NEW.get_current_level < AsgMaxLevel
        if GDSGE_ASG_FIX_GRID==1
            [GDSGE_TEMP_grids, GDSGE_TEMP_surplus, GDSGE_TEMP_levels, GDSGE_TEMP_unscaledGrids] = GDSGE_ASG_INTERP.get_grids_info_at_level(GDSGE_ASG_INTERP_NEW.get_current_level+1);
            GDSGE_evalArrayIdx = [];
            for GDSGE_I_ARRAY=1:GDSGE_OPTIONS.WarmUp.asg_interp_struct.numArray
                GDSGE_evalArrayIdx = [GDSGE_evalArrayIdx,GDSGE_I_ARRAY*ones(1,size(GDSGE_TEMP_grids{GDSGE_I_ARRAY},2))];
            end
            GDSGE_evalGrids = cat(2,GDSGE_TEMP_grids{:});
            GDSGE_evalGridsUnscaled = cat(2,GDSGE_TEMP_unscaledGrids{:});
        else
            if GDSGE_ASG_INTERP_NEW.get_current_level<AsgMinLevel
                [GDSGE_evalArrayIdx,GDSGE_evalGrids,GDSGE_evalGridsLength,GDSGE_evalGridsUnscaled] = GDSGE_ASG_INTERP_NEW.get_eval_grids(0.0);
            else
                [GDSGE_evalArrayIdx,GDSGE_evalGrids,GDSGE_evalGridsLength,GDSGE_evalGridsUnscaled] = GDSGE_ASG_INTERP_NEW.get_eval_grids(AsgThreshold);
            end
        end
        if isempty(GDSGE_evalGrids)
            break;
        end
        
        ASG_PROPOSE_GRIDS_AND_SOLVE
        
        INTERP_ASG_ASSIGN_CODE
        
        GDSGE_evalRslts = [INTERP_VAR_SEMI_COLON];
        
        GDSGE_SOL_ASG_INTERP_NEW.push_eval_results_at_grids(GDSGE_evalArrayIdx, GDSGE_evalGridsUnscaled, GDSGE_SOL, GDSGE_SOL_ASG_INTERP_NEW.get_current_level);
        if GDSGE_ASG_FIX_GRID==1
            GDSGE_ASG_INTERP_NEW.push_eval_results_at_grids(GDSGE_evalArrayIdx, GDSGE_evalGridsUnscaled, GDSGE_evalRslts, GDSGE_ASG_INTERP_NEW.get_current_level);
        else
            GDSGE_ASG_INTERP_NEW.push_eval_results(GDSGE_evalRslts);
        end
        GDSGE_ASG_STORE_evalArrayIdx = [GDSGE_ASG_STORE_evalArrayIdx,GDSGE_evalArrayIdx];
        GDSGE_ASG_STORE_evalGridsUnscaled = [GDSGE_ASG_STORE_evalGridsUnscaled,GDSGE_evalGridsUnscaled];
        GDSGE_ASG_STORE_output = [GDSGE_ASG_STORE_output,[OUTPUT_VAR_SEMI_COLON]];
    end
    PRE_UPDATE_CODE
    
    % Compute Metric
    [GDSGE_Metric,GDSGE_MetricVec] = asg.compute_inf_metric(GDSGE_ASG_INTERP_NEW, GDSGE_ASG_INTERP);
    
    % Update
    GDSGE_ASG_INTERP = GDSGE_ASG_INTERP_NEW;
    GDSGE_SOL_ASG_INTERP = GDSGE_SOL_ASG_INTERP_NEW;
    
    POST_ITER_CODE
    
    stopFlag = GDSGE_Metric<TolEq || GDSGE_Iter>=MaxIter;
    
    if ( mod(GDSGE_Iter,PrintFreq)==0 || stopFlag == true )
      fprintf(['Iter:%d, Metric:%g, maxF:%g\n'],GDSGE_Iter,GDSGE_Metric,max(GDSGE_F));
      toc;
      tic;
    end
    
    if ( mod(GDSGE_Iter,SaveFreq)==0 || stopFlag == true )
        % Construct output
        % GDSGE_ASG_OUTPUT = asg({STATE_COMMA},GDSGE_NUM_OUTPUT,shock_num);
        % OUTPUT_xxx_CONSTRUCT_CODE
        
        % Solve the problem and get output variables
        GDSGE_ASG_HANDLE = GDSGE_ASG_INTERP.objectHandle;
        GDSGE_ASG_INTERP_NEW = asg({STATE_COMMA},GDSGE_NUM_OUTPUT,shock_num);
        
        OUTPUT_CONSTRUCT_CODE

        IterRslt.Metric = GDSGE_Metric;
        IterRslt.MetricVec = GDSGE_MetricVec;
        IterRslt.Iter = GDSGE_Iter;
        IterRslt.shock_num = shock_num;
        IterRslt.shock_trans = shock_trans;
        IterRslt.var_shock = v2struct(RSLT_SHOCK);
        IterRslt.var_state = v2struct(RSLT_STATE);
        IterRslt.params = v2struct(RSLT_PARAMS);
        IterRslt.asg_interp_struct = GDSGE_ASG_INTERP.convert_to_struct();
        IterRslt.sol_asg_interp_struct = GDSGE_SOL_ASG_INTERP.convert_to_struct();
        IterRslt.var_others = v2struct(RSLT_VAR_OTHERS);

        if ~NoSave
            if IterSaveAll
                save(['IterRslt_MODEL_NAME_' num2str(GDSGE_Iter) '.mat']);
            else
                save(['IterRslt_MODEL_NAME_' num2str(GDSGE_Iter) '.mat'],'IterRslt');
            end
        end
    end
end