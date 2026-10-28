function tf = isSimuMexExpressible(ir)
% ISSIMUMEXEXPRESSIBLE  true if a spline model's simulate block can run in
%   simulate_<model>_mex: var_simu entries resolve to output vars, and every
%   transition RHS is a single output-var name (primed or unprimed), and no
%   pre/post_simulate_iter hook is declared. Otherwise
%   the SIMU_INTERP variant uses the per-period interp_eval_mex fallback loop.
%   Deliberately independent of the SIMU_RESOLVE/SIMU_INTERP toggle (ADR-0002):
%   the simulate MEX is generated for every expressible model, so the compile
%   cache is invariant to the toggle.
tf = false;
if strcmp(ir.options.interpMethod, 'asg'); return; end
% a pre/post_simulate_iter hook is MATLAB run every period: MATLAB loop only
if ~isempty(gdsge.codegen.mat.simuHook(ir, 'pre')) || ...
        ~isempty(gdsge.codegen.mat.simuHook(ir, 'post')); return; end
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
