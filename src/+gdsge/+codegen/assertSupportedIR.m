function assertSupportedIR(ir)
% ASSERTSUPPORTEDIR  Defense-in-depth invariant: refuse IR features the
%   refactored pipeline does not yet generate code for. Called at the top of
%   the IR-consuming generators (generateMatlab/generateCxx) so a hand-authored
%   or hand-edited IR cannot silently slip an unsupported combination past codegen.
%   Phase 9a: MATLAB-side var_tensor (feeding inbound/initial) is supported; a
%   tensor used in the model body is rejected by analyzeModel. ASG would need
%   the deferred C++-pop path, so var_tensor + ASG is rejected here.
if ~isempty(ir.variables.tensor) && strcmp(ir.options.interpMethod, 'asg')
    names = cellfun(@(t) t.name, ir.variables.tensor, 'UniformOutput', false);
    error('gdsge:codegen:varTensorAsgUnsupported', ...
        ['var_tensor is not supported with ASG (tensor(s): %s); the C++-pop ' ...
         'path ASG would require is deferred.'], strjoin(names, ', '));
end
% Global var_aux (GDSGE_DER) needs spline construction + analytic gradient
% evaluation on the MATLAB side; asg_mex exposes no gradient command yet.
if isfield(ir.variables, 'auxGlobal') && ~isempty(ir.variables.auxGlobal) ...
        && ~strcmp(ir.options.interpMethod, 'spline')
    names = cellfun(@(g) g.name, ir.variables.auxGlobal, 'UniformOutput', false);
    error('gdsge:codegen:globalAuxAsgUnsupported', ...
        ['Global var_aux (%s) require interpMethod = spline; ASG support needs ' ...
         'an asg_mex gradient command (deferred).'], strjoin(names, ', '));
end
end
