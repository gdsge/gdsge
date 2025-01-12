while GDSGE_ASG_INTERP_NEW.get_current_level < AsgOutputMaxLevel
    % Fix grids
    [GDSGE_TEMP_grids, GDSGE_TEMP_surplus, GDSGE_TEMP_levels, GDSGE_TEMP_unscaledGrids] = GDSGE_ASG_INTERP.get_grids_info_at_level(GDSGE_ASG_INTERP_NEW.get_current_level+1);
    GDSGE_evalArrayIdx = [];
    for GDSGE_I_ARRAY=1:shock_num
        GDSGE_evalArrayIdx = [GDSGE_evalArrayIdx,GDSGE_I_ARRAY*ones(1,size(GDSGE_TEMP_grids{GDSGE_I_ARRAY},2))];
    end
    GDSGE_evalGrids = cat(2,GDSGE_TEMP_grids{:});
    GDSGE_evalGridsUnscaled = cat(2,GDSGE_TEMP_unscaledGrids{:});
    if isempty(GDSGE_evalGrids)
        break;
    end
    
    ASG_PROPOSE_GRIDS_AND_SOLVE
    
    GDSGE_evalRslts = [OUTPUT_VAR_SEMI_COLON];
    % Fix grids
    GDSGE_ASG_INTERP_NEW.push_eval_results_at_grids(GDSGE_evalArrayIdx(GDSGE_solved), GDSGE_evalGridsUnscaled(:, GDSGE_solved), GDSGE_evalRslts(:, GDSGE_solved), GDSGE_ASG_INTERP_NEW.get_current_level);
end

OUTPUT_INDEX_ASSIGN_CODE
IterRslt.output_var_index = output_var_index;
IterRslt.asg_output_struct = GDSGE_ASG_INTERP_NEW.convert_to_struct();
