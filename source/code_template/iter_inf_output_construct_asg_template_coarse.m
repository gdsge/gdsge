for i_level=1:length(GDSGE_ASG_STORE_evalArrayIdx)
    GDSGE_ASG_OUTPUT.push_eval_results_at_grids(GDSGE_ASG_STORE_evalArrayIdx{i_level}, ...
        GDSGE_ASG_STORE_evalGridsUnscaled{i_level}, ...
        GDSGE_ASG_STORE_output{i_level}, GDSGE_ASG_OUTPUT.get_current_level);
end

OUTPUT_INDEX_ASSIGN_CODE
IterRslt.output_var_index = output_var_index;
IterRslt.asg_output_struct = GDSGE_ASG_OUTPUT.convert_to_struct();
