function txt = emitTask(ir, taskName)
% EMITTASK (sympy)  Assemble task.tpl.cpp for the analytic-Jacobian backend:
%   sympymodel.generate body + the double with-grad interp evaluator + the
%   shared CoDoSol call block. All-double; no adept.
if nargin < 2; taskName = 'task_inf_horizon'; end
import gdsge.codegen.cxx.fillTemplate
import gdsge.codegen.cxx.readTemplate
interp = gdsge.codegen.cxx.emitInterpSympy(ir);

regions = ir.model.regions;
modelParts = cell(1, numel(regions));
callParts  = cell(1, numel(regions));
numEqs = zeros(1, numel(regions));
g = [];
for k = 1:numel(regions)
    v = gdsge.ir.regionView(ir, k);
    g = gdsge.codegen.cxx.sympymodel.generate(v);
    numEqs(k) = g.numEq;
    mk = fillTemplate(readTemplate('model_sympy.tpl.cpp'), { ...
        'ARG_CODE',        g.argCode; ...
        'MODEL_BODY_CODE', g.bodyCode; ...
        'AUX_ASSIGN_CODE', auxAssignDouble(ir)});
    mk = fillTemplate(mk, {'NUM_EQUATIONS', num2str(g.numEq)});
    modelParts{k} = strrep(mk, 'MODEL_NUMBER', sprintf('MODEL_%d', k));
    cfm = fillTemplate(readTemplate('call_fmin_sympy.tpl.cpp'), { ...
        'MODEL_CONDITION', gdsge.codegen.cxx.lowerCondition(regions{k}.condition, ir); ...
        'NUM_EQUATIONS', num2str(g.numEq)});
    callParts{k} = strrep(cfm, 'MODEL_NUMBER', sprintf('MODEL_%d', k));
end
assert(all(numEqs == numEqs(1)), 'sympy: regions have unequal equation counts (%s)', mat2str(numEqs));
% g.numAux / g.popCode / g.argCode are region-invariant (same unknowns/aux); use
% the last region's g for them below.
model = strjoin(modelParts, newline);
callFMin = strjoin(callParts, newline);

txt = fillTemplate(readTemplate('task.tpl.cpp'), { ...
    'PRE_MODEL_CODE',        ir.hooks.preModel; ...
    'START_LOOP_CODE',       ''; ...
    'FINISH_LOOP_CODE',      ''; ...
    'MODEL_CODE',            model; ...
    'CALL_FMIN_CODE',        callFMin; ...
    'POP_CODE',              g.popCode; ...
    'INTERP_GET_CODE',       interp.getCode; ...
    'INTERP_GET_THREAD_CODE', interp.threadCode; ...
    'NUM_EQUATIONS',         num2str(g.numEq); ...
    'NUM_AUX',               num2str(g.numAux); ...
    'TASK_NAME',             taskName});
end

function txt = auxAssignDouble(ir)
% aux outputs, plain double (no adept value()).
w = gdsge.codegen.codeWriter();
for i = 1:numel(ir.variables.aux)
    a = ir.variables.aux{i};
    if a.length == 1
        w.add('GDSGE_aux[%d]=%s;', a.slot(1) - 1, a.name);
    else
        w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++)', a.length);
        w.add('{ GDSGE_aux[%d+GDSGE_iter]=%s_GDSGE_GRID[GDSGE_iter-1]; }', ...
            a.slot(1) - 2, a.name);
    end
end
txt = w.str();
end
