function frag = emitBounds(ir, mode)
% EMITBOUNDS  Bound-related fragments from ir.bounds + policy slot layout.
%   .init            GDSGE_LB/UB allocation + per-variable fills (iter)
%   .initConstant    same, omitting per-point fills — those referencing a
%                    state, shock or var_tensor (simulate pre-loop; the iter
%                    GDSGE_TENSOR_* grids don't exist there)
%   .perPeriod       only the per-point fills, in simulate form: state s ->
%                    s(:)' (the per-period sample column), shock z -> its
%                    per-sample value, tensor t -> t(:)' (simulate per-period,
%                    after the state assignment and tensor recompute)
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
pointNames = [sNames, tnames];
% simulate per-period form: problems are samples, states are the per-period
% columns, shocks are looked up by each sample's shock index
simuReps = [cellfun(@(n) [n '(:)'''], ir.states.names, 'UniformOutput', false), ...
    cellfun(@(n) ['reshape(IterRslt.var_shock.' n '(shock),1,[])'], ir.shocks.names, 'UniformOutput', false), ...
    cellfun(@(n) [n '(:)'''], tnames, 'UniformOutput', false)];
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
wConst = gdsge.codegen.codeWriter();
wPoint = gdsge.codegen.codeWriter();
allocLB = sprintf('GDSGE_LB = zeros(%d,GDSGE_NPROB);', maxDim);
allocUB = sprintf('GDSGE_UB = 1e3*ones(%d,GDSGE_NPROB);', maxDim);
wInit.add('%s', allocLB); wInit.add('%s', allocUB);
wConst.add('%s', allocLB);  wConst.add('%s', allocUB);
wWiden = gdsge.codegen.codeWriter();
wTight = gdsge.codegen.codeWriter();
adaptiveSlots = [];
for i = 1:numel(ir.bounds)
    b = ir.bounds{i};
    s = slotOf.(b.name);
    addFill(wInit, wConst, wPoint, sprintf('GDSGE_LB(%d:%d,:)=', s(1), s(2)), ...
        b.lower, names, reps, pointNames, simuReps);
    addFill(wInit, wConst, wPoint, sprintf('GDSGE_UB(%d:%d,:)=', s(1), s(2)), ...
        b.upper, names, reps, pointNames, simuReps);
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
frag = struct('init', wInit.str(), 'initConstant', wConst.str(), ...
    'perPeriod', wPoint.str(), 'adaptiveWiden', wWiden.str(), ...
    'adaptiveTight', wTight.str(), ...
    'adaptiveSlots', adaptiveSlots, 'maxDim', maxDim);
end

% Every bound fill goes to wInit (iter, mode-rewritten form, byte-identical to
% the pre-split output). A per-point fill (references a state, shock or
% var_tensor) also goes to wPoint in simulate form; a constant fill to wConst.
function addFill(wInit, wConst, wPoint, lhs, exprText, names, reps, pointNames, simuReps)
wInit.add('%s', [lhs gdsge.codegen.mat.rewriteNames(exprText, names, reps) ';']);
if refsAny(exprText, pointNames)
    wPoint.add('%s', [lhs gdsge.codegen.mat.rewriteNames(exprText, pointNames, simuReps) ';']);
else
    wConst.add('%s', [lhs exprText ';']);
end
end

function tf = refsAny(exprText, nameList)
tf = false;
for i = 1:numel(nameList)
    if ~isempty(regexp(exprText, ['(?<![A-Za-z0-9_])' regexptranslate('escape', nameList{i}) '(?![A-Za-z0-9_])'], 'once'))
        tf = true; return;
    end
end
end
