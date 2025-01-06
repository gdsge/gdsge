% This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE
function SimuRslt = simulate_MODEL_NAME(IterRslt,GDSGE_OPTIONS)
%% Add path
if ispc
    BLAS_FILE = 'essential_blas.dll';
    PATH_DELIMITER = ';';
    
    GDSGE_TOOLBOX_ROOT = fileparts(which(BLAS_FILE));
    if ~any(strcmp(strsplit(getenv('PATH'),PATH_DELIMITER),GDSGE_TOOLBOX_ROOT))
        setenv('PATH',[getenv('PATH'),PATH_DELIMITER,GDSGE_TOOLBOX_ROOT]);
    end
    
    clear BLAS_FILE PATH_DELIMITER GDSGE_TOOLBOX_ROOT
elseif ismac
end

%% Simulate code starts here
PARAMS_PRESET_CODE

SIMU_PRE_CODE

GEN_SHOCK_START_PERIOD = 1;

shock_trans = IterRslt.shock_trans;
v2struct(IterRslt.params);
v2struct(IterRslt.var_shock);
v2struct(IterRslt.var_state);

if nargin>=2
    v2struct(GDSGE_OPTIONS)
end

% Reconstruct asg interpolations
GDSGE_ASG_INTERP = asg.construct_from_struct(IterRslt.asg_interp_struct);
GDSGE_ASG_HANDLE = GDSGE_ASG_INTERP.objectHandle;

GDSGE_NPROB = num_samples;

SimuRslt.shock = ones(num_samples,num_periods+1);
SIMU_CONSTRUCT_CODE

SIMULATE_INIT_CODE

if nargin>1 && isfield(GDSGE_OPTIONS,'init')
    SIMULATE_INIT_OVERWRITE_CODE
end

if any(SimuRslt.shock(:,1)>shock_num)
    error('initial shock exceeds shock_num');
end

POLICY_INBOUND_CODE

GDSGE_SKIP_MODEL_ID = zeros(1,GDSGE_NPROB);
GDSGE_EQVAL_MODEL_ID = 1e20*ones(GDSGE_MAXDIM,GDSGE_NPROB);
GDSGE_F_MODEL_ID = 1e20*ones(1,GDSGE_NPROB);
GDSGE_SOL_MODEL_ID = zeros(GDSGE_MAXDIM,GDSGE_NPROB);
GDSGE_data0_MODEL_ID = repmat([shock_num;PARAMS_SEMI_COLON;shock_trans(:);SHOCK_SEMI_COLON],1,GDSGE_NPROB);
GDSGE_AUX_MODEL_ID = zeros(GDSGE_NUM_AUX,GDSGE_NPROB);

%% Generate random number
SimuRslt.shock(:,GEN_SHOCK_START_PERIOD:end) = gen_discrete_markov_rn(shock_trans,num_samples,length(GEN_SHOCK_START_PERIOD:num_periods+1),...
    SimuRslt.shock(:,GEN_SHOCK_START_PERIOD));

MEX_TASK_NAME = MEX_TASK_INF_HORIZON;
GDSGE_SHOCK_VAR_INDEX_BASE = ([0:num_samples-1]')*shock_num;

%% Extract interpolation for sol
GDSGE_SOL_ASG_INTERP = asg.construct_from_struct(IterRslt.sol_asg_interp_struct);

tic;
for GDSGE_t=1:num_periods
    % Reuse the init inbound and tensor code
    % Map grid to current variable
    shock = SimuRslt.shock(:,GDSGE_t);
    SIMU_PRE_ASSIGN_CODE
    
    GDSGE_SOL = GDSGE_SOL_ASG_INTERP.eval_vec(SimuRslt.shock(:,GDSGE_t)',[SIMU_RSLT_STATE_SEMI_COLON]);
    
    if UseAdaptiveBoundInSol==1
        INBOUND_ADAPTIVE_TIGHT_CODE
    end

    % Construct tensor variable
    SIMU_TENSOR
    
    % Reconstruct data
    GDSGE_DATA_MODEL_ID = [GDSGE_data0_MODEL_ID;shock(:)';STATE_SEMI_COLON TENSOR_SEMI_COLON];
    
    % Solve problems
    GDSGE_SKIP_MODEL_ID = zeros(1,GDSGE_NPROB);
    [GDSGE_SOL_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_AUX_MODEL_ID,GDSGE_EQVAL_MODEL_ID,GDSGE_OPT_INFO] = mex_MODEL_NAME(GDSGE_SOL_MODEL_ID,GDSGE_LB_MODEL_ID,GDSGE_UB_MODEL_ID,GDSGE_DATA_MODEL_ID,GDSGE_SKIP_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_AUX_MODEL_ID,GDSGE_EQVAL_MODEL_ID);
    while (max(isnan(GDSGE_F_MODEL_ID)) || max(GDSGE_F_MODEL_ID(:))>TolSol)
        GDSGE_X0Rand = rand(size(GDSGE_SOL_MODEL_ID)) .* (GDSGE_UB_MODEL_ID-GDSGE_LB_MODEL_ID) + GDSGE_LB_MODEL_ID;
        NeedResolved = (GDSGE_F_MODEL_ID>TolSol) | isnan(GDSGE_F_MODEL_ID);
        GDSGE_SOL_MODEL_ID(:,NeedResolved) = GDSGE_X0Rand(:,NeedResolved);
        GDSGE_SKIP_MODEL_ID(~NeedResolved) = 1;
        
        [GDSGE_SOL_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_AUX_MODEL_ID,GDSGE_EQVAL_MODEL_ID,GDSGE_OPT_INFO] = mex_MODEL_NAME(GDSGE_SOL_MODEL_ID,GDSGE_LB_MODEL_ID,GDSGE_UB_MODEL_ID,GDSGE_DATA_MODEL_ID,GDSGE_SKIP_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_AUX_MODEL_ID,GDSGE_EQVAL_MODEL_ID);
    
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
        end
    end
    
    POLICY_ASSIGN_CODE
    
    AUX_ASSIGN_CODE
    
    GDSGE_SHOCK_VAR_LINEAR_INDEX = SimuRslt.shock(:,GDSGE_t+1) + GDSGE_SHOCK_VAR_INDEX_BASE;
    SIMU_POST_ASSIGN_CODE
    
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
        toc;
        tic;
    end
    
    if mod(GDSGE_t,SimuSaveFreq)==0
        save(['SimuRslt_MODEL_NAME_' num2str(GDSGE_t) '.mat'], 'SimuRslt');
    end
end
end