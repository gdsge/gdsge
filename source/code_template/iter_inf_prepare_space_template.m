GDSGE_SIZE_MODEL_ID = size(GDSGE_TENSOR_shockIdx_MODEL_ID);
GDSGE_SIZE_STATE_MODEL_ID = num2cell(GDSGE_SIZE_MODEL_ID(2:end));

POLICY_INBOUND_CODE

GDSGE_EQVAL_MODEL_ID = 1e20*ones(GDSGE_MAXDIM,GDSGE_NPROB_MODEL_ID);
GDSGE_F_MODEL_ID = 1e20*ones(1,GDSGE_NPROB_MODEL_ID);
GDSGE_SOL_MODEL_ID = zeros(GDSGE_MAXDIM,GDSGE_NPROB_MODEL_ID);
GDSGE_X0_MODEL_ID = rand(size(GDSGE_SOL_MODEL_ID)) .* (GDSGE_UB_MODEL_ID-GDSGE_LB_MODEL_ID) + GDSGE_LB_MODEL_ID;
GDSGE_SOL_MODEL_ID(:) = GDSGE_X0_MODEL_ID;
GDSGE_AUX_MODEL_ID = zeros(GDSGE_NUM_AUX,GDSGE_NPROB_MODEL_ID);
GDSGE_SKIP_MODEL_ID = zeros(1,GDSGE_NPROB_MODEL_ID);
GDSGE_DATA_MODEL_ID = zeros(GDSGE_MAXDATA,GDSGE_NPROB_MODEL_ID);