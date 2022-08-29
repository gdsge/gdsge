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
    if exist('./essential_blas.dylib','file') == 0
        copyfile(which('essential_blas.dylib'),'./');
    end
end

%% Simulate code starts here
PARAMS_PRESET_CODE

v2struct(IterRslt.pp);

SIMU_PRE_CODE

GEN_SHOCK_START_PERIOD = 1;

shock_trans = IterRslt.shock_trans;
v2struct(IterRslt.params);
v2struct(IterRslt.var_shock);
v2struct(IterRslt.var_state);

if nargin>=2
    v2struct(GDSGE_OPTIONS)
end

%% Construct interpolation for solutions
if shock_num>1
    GDSGE_PP=struct('form','MKL','breaks',{{1:shock_num,STATE_COMMA}},...
        'Values',reshape(IterRslt.GDSGE_PROB.GDSGE_SOL, [numel(IterRslt.GDSGE_PROB.GDSGE_SOL)/prod(IterRslt.GDSGE_PROB.GDSGE_SIZE),IterRslt.GDSGE_PROB.GDSGE_SIZE]),...
        'coefs',[],'order',[2 OutputInterpOrder*ones(1,length({STATE_COMMA}))],'Method',[],...
        'ExtrapolationOrder',[],'thread',NumThreads, ...
        'orient','curvefit');
else
    GDSGE_PP=struct('form','MKL','breaks',{{STATE_COMMA}},...
        'Values',reshape(IterRslt.GDSGE_PROB.GDSGE_SOL, [numel(IterRslt.GDSGE_PROB.GDSGE_SOL)/prod(IterRslt.GDSGE_PROB.GDSGE_SIZE),IterRslt.GDSGE_PROB.GDSGE_SIZE(2:end)]),...
        'coefs',[],'order',[OutputInterpOrder*ones(1,length({STATE_COMMA}))],'Method',[],...
        'ExtrapolationOrder',[],'thread',NumThreads, ...
        'orient','curvefit');
end
GDSGE_PP=myppual(GDSGE_PP);

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

% Use the largest bound in IterSol
if UseAdaptiveBound==1
GDSGE_LB = repmat(min(IterRslt.GDSGE_PROB.GDSGE_LB,[],2),[1,GDSGE_NPROB]);
GDSGE_UB = repmat(max(IterRslt.GDSGE_PROB.GDSGE_UB,[],2),[1,GDSGE_NPROB]);
end

GDSGE_SKIP_MODEL_ID = zeros(1,GDSGE_NPROB);
GDSGE_EQVAL_MODEL_ID = 1e20*ones(GDSGE_MAXDIM,GDSGE_NPROB);
GDSGE_F_MODEL_ID = 1e20*ones(1,GDSGE_NPROB);
GDSGE_SOL_MODEL_ID = zeros(GDSGE_MAXDIM,GDSGE_NPROB);
GDSGE_AUX_MODEL_ID = zeros(GDSGE_NUM_AUX,GDSGE_NPROB);

%% Generate random number
SimuRslt.shock(:,GEN_SHOCK_START_PERIOD:end) = gen_discrete_markov_rn(shock_trans,num_samples,length(GEN_SHOCK_START_PERIOD:num_periods+1),...
    SimuRslt.shock(:,GEN_SHOCK_START_PERIOD));

MEX_TASK_NAME = MEX_TASK_INF_HORIZON;
GDSGE_SHOCK_VAR_INDEX_BASE = ([0:num_samples-1]')*shock_num;

tic;
for GDSGE_t=1:num_periods
    % Reuse the init inbound and tensor code
    % Map grid to current variable
    shock = SimuRslt.shock(:,GDSGE_t);
    SIMU_PRE_ASSIGN_CODE
    
    GDSGE_data0_MODEL_ID = repmat([shock_num;PARAMS_SEMI_COLON;shock_trans(:);SHOCK_SEMI_COLON],1,GDSGE_NPROB);
    
    SIMU_PRE_ITER_CODE
    
    % Use interpolation as initial values
%     %{
    if GDSGE_t>0
        if shock_num>1
            GDSGE_SOL_MODEL_ID = myppual_mex(int32(NumThreads),GDSGE_PP.breaks,GDSGE_PP.coefs,...
                int32(GDSGE_PP.pieces),int32(GDSGE_PP.order),int32(GDSGE_PP.dim),'not-a-knot',[SimuRslt.shock(:,GDSGE_t)';SIMU_RSLT_STATE_SEMI_COLON],[],[],[]);
        else
            GDSGE_SOL_MODEL_ID = myppual_mex(int32(NumThreads),GDSGE_PP.breaks,GDSGE_PP.coefs,...
                int32(GDSGE_PP.pieces),int32(GDSGE_PP.order),int32(GDSGE_PP.dim),'not-a-knot',[SIMU_RSLT_STATE_SEMI_COLON],[],[],[]);
        end
    end
    %}
    
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
    
    SIMU_POST_ITER_CODE
    
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