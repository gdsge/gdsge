% This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE
function SimuRslt = simulate_MODEL_NAME(IterRslt,GDSGE_OPTIONS)
PARAMS_PRESET_CODE

SIMU_PRE_CODE

GEN_SHOCK_START_PERIOD=1;

v2struct(IterRslt.var_state);

%% Overwrite parameters
if nargin>=2
    v2struct(GDSGE_OPTIONS)
end

%% Reconstruct interpolations
SIMULATE_INTERP_CONSTRUCT_CODE

shock_trans = IterRslt.shock_trans;
output_var_index = IterRslt.output_var_index;

%% Construct space
SimuRslt.shock = ones(num_samples,num_periods+1);
SIMU_CONSTRUCT_CODE

%% Initiate state
SIMULATE_INIT_CODE

if nargin>1 && isfield(GDSGE_OPTIONS,'init')
    SIMULATE_INIT_OVERWRITE_CODE
end

if any(SimuRslt.shock(:,1)>shock_num)
    error('initial shock exceeds shock_num');
end

%% Generate random number
SimuRslt.shock(:,GEN_SHOCK_START_PERIOD:end) = gen_discrete_markov_rn(shock_trans,num_samples,length(GEN_SHOCK_START_PERIOD:num_periods+1),...
    SimuRslt.shock(:,GEN_SHOCK_START_PERIOD));

%% Apply the interpolation for simulation
shock_num = size(shock_trans,1);
GDSGE_SHOCK_VAR_INDEX_BASE = ([0:num_samples-1]')*shock_num;
for GDSGE_t=1:num_periods
    % constrain states inbound
    if EnforceSimuStateInbound==1
        SIMU_STATE_INBOUND_CODE
    end
    
    % Intepolate values
    SIMULATE_INTERP_EVAL_CODE
    
    % Assign to output variables from evaluation results
    ASSIGN_INTERP_TO_OUTPUT_CODE
    
    % Move states forward
    % Linear index
    GDSGE_SHOCK_VAR_LINEAR_INDEX = SimuRslt.shock(:,GDSGE_t+1) + GDSGE_SHOCK_VAR_INDEX_BASE;
    SIMU_POST_ASSIGN_CODE
    
    SIMU_POST_ITER_CODE
    
    % Print something
    if mod(GDSGE_t,SimuPrintFreq)==0
        fprintf('Periods: %d\n', GDSGE_t);
        SimuRsltNames = fieldnames(SimuRslt);
        for GDSGE_field = 1:length(SimuRsltNames)
            fprintf('%8s', SimuRsltNames{GDSGE_field});
        end
        fprintf('\n');
        for GDSGE_field = 1:length(SimuRsltNames)
            fprintf('%8.4g', SimuRslt.(SimuRsltNames{GDSGE_field})(1,GDSGE_t));
        end
        fprintf('\n');
    end
    
    if mod(GDSGE_t,SimuSaveFreq)==0
        save(['SimuRslt_MODEL_NAME_' num2str(GDSGE_t) '.mat'], 'SimuRslt');
    end
end

end