TENSOR_CONSTRUCT_CODE

TENSOR_ASSIGN_CODE

PREPARE_SPACE_CODE

if ~( nargin>=1 && isfield(GDSGE_OPTIONS,'WarmUp') && isfield(GDSGE_OPTIONS.WarmUp,'var_interp') )
    INTERP_INITIALIZE_CODE
    
    CONSTRUCT_SPLINE_CODE
end

GDSGE_Metric = 1;
GDSGE_Iter = 0;
IS_WARMUP_LOOP = 0;

if nargin>=1 && isfield(GDSGE_OPTIONS,'WarmUp')
    if isfield(GDSGE_OPTIONS.WarmUp,'var_interp')
        v2struct(GDSGE_OPTIONS.WarmUp.var_interp);
        GDSGE_TEMP = v2struct(RSLT_STATE,GDSGE_SIZE_STATE);
        if isfield(GDSGE_OPTIONS.WarmUp,'GDSGE_PROB')
        GDSGE_SIZE_STATE = num2cell(GDSGE_OPTIONS.WarmUp.GDSGE_PROB.GDSGE_SIZE(2:end));
        end
        if isfield(GDSGE_OPTIONS.WarmUp,'var_state')
        v2struct(GDSGE_OPTIONS.WarmUp.var_state);    
        end
        CONSTRUCT_SPLINE_CODE
        v2struct(GDSGE_TEMP);
        IS_WARMUP_LOOP = 1;
    end
    if isfield(GDSGE_OPTIONS.WarmUp,'Iter')
        GDSGE_Iter = GDSGE_OPTIONS.WarmUp.Iter;
    end
    if isfield(GDSGE_OPTIONS.WarmUp,'GDSGE_SOL')
        if ~isequal(size(GDSGE_OPTIONS.WarmUp.GDSGE_SOL,1),size(GDSGE_SOL,1))
            error('Wrong size of GDSGE_SOL in WarmUp');
        end
        GDSGE_SOL = GDSGE_OPTIONS.WarmUp.GDSGE_SOL;
    end
    if REUSE_WARMUP_SOL==1 && isfield(GDSGE_OPTIONS.WarmUp,'GDSGE_PROB') && size(GDSGE_OPTIONS.WarmUp.GDSGE_PROB.GDSGE_SOL,1)==size(GDSGE_SOL,1)
        if INTERP_WARMUP_SOL==1 && shock_num>=4
        % Interpolate SOL, LB, and UB
        GDSGE_TEMP = v2struct(RSLT_STATE,GDSGE_SIZE_STATE);
        GDSGE_SIZE_STATE = num2cell(GDSGE_OPTIONS.WarmUp.GDSGE_PROB.GDSGE_SIZE);
        v2struct(GDSGE_OPTIONS.WarmUp.var_state);
        GDSGE_SOL_interp=struct('form','MKL','breaks',{{[1:shock_num],RSLT_STATE}},'Values',reshape(GDSGE_OPTIONS.WarmUp.GDSGE_PROB.GDSGE_SOL,[],GDSGE_SIZE_STATE{:}),'coefs',[],'order',[2*ones(1,length(GDSGE_SIZE_STATE))],'Method',[],'ExtrapolationOrder',[],'thread',NumThreads,'orient','curvefit');
        GDSGE_SOL_interp=myppual(GDSGE_SOL_interp);
        GDSGE_LB_interp=struct('form','MKL','breaks',{{[1:shock_num],RSLT_STATE}},'Values',reshape(GDSGE_OPTIONS.WarmUp.GDSGE_PROB.GDSGE_LB,[],GDSGE_SIZE_STATE{:}),'coefs',[],'order',[2*ones(1,length(GDSGE_SIZE_STATE))],'Method',[],'ExtrapolationOrder',[],'thread',NumThreads,'orient','curvefit');
        GDSGE_LB_interp=myppual(GDSGE_LB_interp);
        GDSGE_UB_interp=struct('form','MKL','breaks',{{[1:shock_num],RSLT_STATE}},'Values',reshape(GDSGE_OPTIONS.WarmUp.GDSGE_PROB.GDSGE_UB,[],GDSGE_SIZE_STATE{:}),'coefs',[],'order',[2*ones(1,length(GDSGE_SIZE_STATE))],'Method',[],'ExtrapolationOrder',[],'thread',NumThreads,'orient','curvefit');
        GDSGE_UB_interp=myppual(GDSGE_UB_interp);
        
        v2struct(GDSGE_TEMP);
        GDSGE_SOL = reshape(myppual(GDSGE_SOL_interp,[GDSGE_TENSOR_shockIdx(:)';STATE_SEMI_COLON]),size(GDSGE_SOL));
        GDSGE_LB_NEW = reshape(myppual(GDSGE_LB_interp,[GDSGE_TENSOR_shockIdx(:)';STATE_SEMI_COLON]),size(GDSGE_LB));
        GDSGE_UB_NEW = reshape(myppual(GDSGE_UB_interp,[GDSGE_TENSOR_shockIdx(:)';STATE_SEMI_COLON]),size(GDSGE_UB));
        
        WARMUP_BOUND_ADAPTIVE_CODE
        else
            if isfield(GDSGE_OPTIONS.WarmUp.GDSGE_PROB,'GDSGE_SOL')
            GDSGE_SOL = GDSGE_OPTIONS.WarmUp.GDSGE_PROB.GDSGE_SOL;
            end
            
            if isfield(GDSGE_OPTIONS.WarmUp.GDSGE_PROB,'GDSGE_LB')
            GDSGE_LB = GDSGE_OPTIONS.WarmUp.GDSGE_PROB.GDSGE_LB;
            end
            
            if isfield(GDSGE_OPTIONS.WarmUp.GDSGE_PROB,'GDSGE_UB')
            GDSGE_UB = GDSGE_OPTIONS.WarmUp.GDSGE_PROB.GDSGE_UB;
            end
        end
    end
end

stopFlag = false;
tic;
while(~stopFlag)
    GDSGE_Iter = GDSGE_Iter+1;
    
    PRE_ITER_CODE

    SOLVE_AND_ASSIGN_CODE
    
    POST_SOL_CODE
    
    RESHAPE_SUB_CODE
    
    INTERP_NEW_ASSIGN_CODE
    
    % Compute Metric
    try
        METRIC_CODE
    catch
        GDSGE_Metric = nan;
    end
    
    GDSGE_Metric0 = GDSGE_Metric;
    
    if IS_WARMUP_LOOP==1
%         GDSGE_Metric = inf;
        IS_WARMUP_LOOP = 0;
    end
    
    % Update
    INTERP_ASSIGN_CODE
    
    INBOUND_ADAPTIVE_CODE
    
    POST_ITER_CODE
    
    stopFlag = (isempty(GDSGE_Metric) || GDSGE_Metric<TolEq) || GDSGE_Iter>=MaxIter;
    
    CONSTRUCT_SPLINE_CODE
    
    if ( ~NoPrint && (mod(GDSGE_Iter,PrintFreq)==0 || stopFlag == true) )
      fprintf(['Iter:%d, Metric:%g, maxF:%g\n'],GDSGE_Iter,GDSGE_Metric,max(GDSGE_F_MODEL_ID));
      toc;
      tic;
    end
    
    if ( (mod(GDSGE_Iter,SaveFreq)==0 || stopFlag == true) )
        RESHAPE_CODE
        
        if CONSTRUCT_OUTPUT==1
        OUTPUT_CONSTRUCT_CODE
        end
        
        if CONSTRUCT_OUTPUT==1
        IterRslt.shock_num = shock_num;
        IterRslt.shock_trans = shock_trans;
        IterRslt.params = v2struct(RSLT_PARAMS);
        IterRslt.var_shock = v2struct(RSLT_SHOCK);
        IterRslt.var_policy = get_scalar(v2struct(RSLT_POLICY),{[1:shock_num],RSLT_STATE});
        IterRslt.var_aux = get_scalar(v2struct(RSLT_AUX),{[1:shock_num],RSLT_STATE});
        IterRslt.var_tensor = v2struct(RSLT_TENSOR);
        IterRslt.pp = v2struct(RSLT_PP,GDSGE_SPLINE_VEC,GDSGE_EMPTY);
        end
        IterRslt.Metric0 = GDSGE_Metric0;
        IterRslt.Metric = GDSGE_Metric;
        IterRslt.Iter = GDSGE_Iter;
        IterRslt.var_state = v2struct(RSLT_STATE);
        IterRslt.var_interp = v2struct(RSLT_INTERP);
        IterRslt.GDSGE_PROB = v2struct(GDSGE_LB,GDSGE_UB,GDSGE_SOL_MODEL_ID,GDSGE_F_MODEL_ID,GDSGE_SIZE_MODEL_ID);
        IterRslt.var_others = v2struct(RSLT_VAR_OTHERS);
        IterRslt.NeedResolved = NeedResolved;

        if ~NoSave
            if IterSaveAll
                save(['IterRslt_MODEL_NAME_' num2str(GDSGE_Iter) '.mat']);
            else
                save(['IterRslt_MODEL_NAME_' num2str(GDSGE_Iter) '.mat'],'IterRslt');
            end
        end
    end
end
