function tf = isSimuMexExpressible(ir)
% ISSIMUMEXEXPRESSIBLE  true if a SIMU_INTERP spline model's simulate block can
%   run in simulate_<model>_mex: var_simu entries resolve to output vars, and
%   every transition RHS is a single output-var name (primed or unprimed).
%   Otherwise the model uses the per-period interp_eval_mex fallback loop.
tf = false;
if strcmp(ir.options.interpMethod, 'asg'); return; end
if ~(isfield(ir.options,'simuInterp') && ir.options.simuInterp == 1); return; end
outNames = ir.variables.output;
isOut = @(n) any(strcmp(n, outNames));
for k = 1:numel(ir.simulate.varSimu)
    if ~isOut(ir.simulate.varSimu{k}); return; end
end
for k = 1:numel(ir.simulate.transitions)
    e = strtrim(ir.simulate.transitions{k}.expr);
    if isempty(regexp(e, '^[A-Za-z_]\w*$', 'once')) || ~isOut(e); return; end
end
tf = true;
end
