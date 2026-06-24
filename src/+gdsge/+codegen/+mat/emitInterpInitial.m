function txt = emitInterpInitial(ir)
% EMITINTERPINITIAL  Interp-variable preallocation + initial expressions,
%   with state/shock references rewritten to GDSGE_TENSOR_<name>(:)'.
tnames = cellfun(@(t) t.name, ir.variables.tensor, 'UniformOutput', false);
names = [ir.states.names, ir.shocks.names, tnames];
reps = cellfun(@(n) ['GDSGE_TENSOR_' n '(:)'''], names, 'UniformOutput', false);
w = gdsge.codegen.codeWriter();
for i = 1:numel(ir.interp)
    w.add('%s=zeros(GDSGE_SIZE);', ir.interp{i}.name);
end
for i = 1:numel(ir.interp)
    it = ir.interp{i};
    w.add('%s(:)=%s;', it.name, ...
        gdsge.codegen.mat.rewriteNames(it.initialExpr, names, reps));
end
txt = w.str();
end
