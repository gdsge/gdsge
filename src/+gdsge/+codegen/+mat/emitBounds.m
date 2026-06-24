function frag = emitBounds(ir, mode)
% EMITBOUNDS  Bound-related fragments from ir.bounds + policy slot layout.
%   .init            GDSGE_LB/UB allocation + per-variable fills (iter)
%   .initNoTensor    same, omitting fills that reference a var_tensor (simulate
%                    pre-loop; tensors don't exist there yet)
%   .tensorInit      only the tensor-referencing fills, stripped of the
%                    GDSGE_TENSOR_ prefix (simulate per-period, after recompute)
%   .adaptiveWiden   in-loop widening      (old: min/max against SOL*factor)
%   .adaptiveTight   UseAdaptiveBoundInSol (old: SOL/factor .. SOL*factor)
%   .adaptiveSlots   row vector of adaptive slot indices (warmup LB/UB
%                    overwrite happens inside gdsge.runtime.applyWarmUp)
%   .maxDim          number of solution rows
%   mode (optional, default 'tensor'): selects how state/shock names are
%   rewritten in state-dependent bound expressions. 'tensor' uses the
%   ndgrid tensors (GDSGE_TENSOR_<name>(:)'); 'asg' uses the proposed
%   per-problem eval grids (GDSGE_evalGrids(k,:) for states,
%   <shock>(GDSGE_evalArrayIdx) for shocks).
if nargin < 2; mode = 'tensor'; end
pol = ir.variables.policy;
slotOf = struct();
maxDim = 0;
for i = 1:numel(pol)
    slotOf.(pol{i}.name) = pol{i}.slot;
    maxDim = max(maxDim, pol{i}.slot(2));
end
tnames = cellfun(@(t) t.name, ir.variables.tensor, 'UniformOutput', false);
sNames = [ir.states.names, ir.shocks.names];
if strcmp(mode, 'asg')
    % ASG: per-problem columns are the proposed eval grids, not tensors.
    % (ASG + var_tensor is rejected, so no tensor names enter here.)
    names = sNames;
    reps = {};
    for k = 1:numel(ir.states.names)
        reps{end+1} = sprintf('GDSGE_evalGrids(%d,:)', k); %#ok<AGROW>
    end
    for k = 1:numel(ir.shocks.names)
        reps{end+1} = sprintf('%s(GDSGE_evalArrayIdx)', ir.shocks.names{k}); %#ok<AGROW>
    end
else
    names = [sNames, tnames];
    reps = cellfun(@(n) ['GDSGE_TENSOR_' n '(:)'''], names, 'UniformOutput', false);
end

wInit = gdsge.codegen.codeWriter();
wNoT  = gdsge.codegen.codeWriter();
wTensor = gdsge.codegen.codeWriter();
allocLB = sprintf('GDSGE_LB = zeros(%d,GDSGE_NPROB);', maxDim);
allocUB = sprintf('GDSGE_UB = 1e3*ones(%d,GDSGE_NPROB);', maxDim);
wInit.add('%s', allocLB); wInit.add('%s', allocUB);
wNoT.add('%s', allocLB);  wNoT.add('%s', allocUB);
wWiden = gdsge.codegen.codeWriter();
wTight = gdsge.codegen.codeWriter();
adaptiveSlots = [];
for i = 1:numel(ir.bounds)
    b = ir.bounds{i};
    s = slotOf.(b.name);
    lo = gdsge.codegen.mat.rewriteNames(b.lower, names, reps);
    hi = gdsge.codegen.mat.rewriteNames(b.upper, names, reps);
    addClassified(wInit, wNoT, wTensor, ...
        sprintf('GDSGE_LB(%d:%d,:)=%s;', s(1), s(2), lo), refsTensor(b.lower, tnames));
    addClassified(wInit, wNoT, wTensor, ...
        sprintf('GDSGE_UB(%d:%d,:)=%s;', s(1), s(2), hi), refsTensor(b.upper, tnames));
    if isfield(b, 'adaptiveFactor') && ~isempty(b.adaptiveFactor)
        fct = mat2str(b.adaptiveFactor, 17);
        wWiden.add('GDSGE_LB(%d:%d,:)=min(GDSGE_LB(%d:%d,:),GDSGE_SOL(%d:%d,:)*%s);', ...
            s(1), s(2), s(1), s(2), s(1), s(2), fct);
        wWiden.add('GDSGE_UB(%d:%d,:)=max(GDSGE_UB(%d:%d,:),GDSGE_SOL(%d:%d,:)*%s);', ...
            s(1), s(2), s(1), s(2), s(1), s(2), fct);
        wTight.add('GDSGE_LB(%d:%d,:)=GDSGE_SOL(%d:%d,:)/%s;', s(1), s(2), s(1), s(2), fct);
        wTight.add('GDSGE_UB(%d:%d,:)=GDSGE_SOL(%d:%d,:)*%s;', s(1), s(2), s(1), s(2), fct);
        adaptiveSlots = [adaptiveSlots, s(1):s(2)]; %#ok<AGROW>
    end
end
frag = struct('init', wInit.str(), 'initNoTensor', wNoT.str(), ...
    'tensorInit', wTensor.str(), 'adaptiveWiden', wWiden.str(), ...
    'adaptiveTight', wTight.str(), ...
    'adaptiveSlots', adaptiveSlots, 'maxDim', maxDim);
end

% A bound fill referencing a var_tensor goes to wInit (iter, full tensor form)
% AND wTensor (simulate, stripped of GDSGE_TENSOR_); everything else goes to
% wInit AND wNoT (so wInit stays byte-identical to the pre-Phase-9a output).
function addClassified(wInit, wNoT, wTensor, line, isTensor)
wInit.add('%s', line);
if isTensor
    wTensor.add('%s', strrep(line, 'GDSGE_TENSOR_', ''));
else
    wNoT.add('%s', line);
end
end

function tf = refsTensor(exprText, tnames)
tf = false;
for i = 1:numel(tnames)
    if ~isempty(regexp(exprText, ['\<' tnames{i} '\>'], 'once')); tf = true; return; end
end
end
