function v = regionView(ir, k)
% REGIONVIEW  Adapt one model region into an ir-like struct whose .model has the
%   flat .statements/.equations/.condition the per-body C++ emitters read. The
%   drivers (emitModel, sympy emitTask, analyzeModel) loop regions through this,
%   so emitModelBody/emitEquations/modelVars/generate stay region-agnostic.
if nargin < 2; k = 1; end
v = ir;
r = ir.model.regions{k};
v.model = struct('statements', {r.statements}, 'equations', {r.equations}, ...
                 'condition', r.condition);
end
