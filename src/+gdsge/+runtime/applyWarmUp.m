function [sol, lb, ub] = applyWarmUp(warmUp, sol, lb, ub, spec)
% APPLYWARMUP  Transfer a previous run's solution/bounds into this run.
%   Behavior-parity port of the old generated WarmUp block (plain GDSGE_SOL
%   copy + the REUSE_WARMUP_SOL interp/direct transfer; see
%   scratch/old_generated/iter_HL1996.m lines 257-304).
%   spec fields: shockNum, stateNames (decl order), evalPoints
%     ([shockIdx(:)'; state tensors...]), adaptiveSlots, reuseSol, interpSol,
%     numThreads.
if isfield(warmUp, 'GDSGE_SOL')
    if ~isequal(size(warmUp.GDSGE_SOL, 1), size(sol, 1))
        error('gdsge:runtime:warmUpSize', 'Wrong size of GDSGE_SOL in WarmUp');
    end
    sol = warmUp.GDSGE_SOL;
end
if ~(spec.reuseSol == 1 && isfield(warmUp, 'GDSGE_PROB') ...
        && size(warmUp.GDSGE_PROB.GDSGE_SOL, 1) == size(sol, 1))
    return;
end
if spec.interpSol == 1
    if spec.shockNum < 2
        error('gdsge:runtime:warmUpShock', ...
            'WarmUp solution can only be applied with shock_num>=2. Please set REUSE_WARMUP_SOL=0.');
    end
    sizeFull = num2cell(warmUp.GDSGE_PROB.GDSGE_SIZE);
    breaks = cell(1, 1 + numel(spec.stateNames));
    breaks{1} = 1:spec.shockNum;
    for s = 1:numel(spec.stateNames)
        breaks{1+s} = warmUp.var_state.(spec.stateNames{s});
    end
    solI = mkInterp(warmUp.GDSGE_PROB.GDSGE_SOL, breaks, sizeFull, spec.numThreads);
    lbI  = mkInterp(warmUp.GDSGE_PROB.GDSGE_LB,  breaks, sizeFull, spec.numThreads);
    ubI  = mkInterp(warmUp.GDSGE_PROB.GDSGE_UB,  breaks, sizeFull, spec.numThreads);
    sol   = reshape(evalAll(solI, spec.evalPoints, spec.numThreads), size(sol));
    lbNew = reshape(evalAll(lbI,  spec.evalPoints, spec.numThreads), size(lb));
    ubNew = reshape(evalAll(ubI,  spec.evalPoints, spec.numThreads), size(ub));
    lb(spec.adaptiveSlots, :) = lbNew(spec.adaptiveSlots, :);
    ub(spec.adaptiveSlots, :) = ubNew(spec.adaptiveSlots, :);
else
    if isfield(warmUp.GDSGE_PROB, 'GDSGE_SOL'); sol = warmUp.GDSGE_PROB.GDSGE_SOL; end
    if isfield(warmUp.GDSGE_PROB, 'GDSGE_LB');  lb  = warmUp.GDSGE_PROB.GDSGE_LB;  end
    if isfield(warmUp.GDSGE_PROB, 'GDSGE_UB');  ub  = warmUp.GDSGE_PROB.GDSGE_UB;  end
end
end

function pp = mkInterp(vals, breaks, sizeFull, numThreads)
% order-2 (linear) transfer interp over {1:shock_num, state grids...} — parity.
% Uniform order 2 across all dims (shock is a real spline dim here), so the
% generic interp_construct_mex/interp_eval_mex kernels evaluate it directly.
S = struct('GDSGE_V1', reshape(vals, [], sizeFull{:}));
[ppCell, ~] = interp_construct_mex(breaks, S, ...
    int32(2*ones(1, numel(sizeFull))), int32([]), int32(numThreads));
pp = ppCell{1};
end

function v = evalAll(pp, sites, numThreads)
% Evaluate every stacked vector of pp at each site column (sites include the
% shock index as the first coordinate, matching breaks{1}=1:shock_num).
nVec = double(pp.dim);
v = interp_eval_mex(int32(numThreads), pp, sites, int32((0:nVec-1)'));
end
