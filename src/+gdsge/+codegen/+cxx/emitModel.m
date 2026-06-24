function [modelCode, callFMinCode, numEq] = emitModel(ir)
% EMITMODEL  Assemble the model lambdas (adouble + double + OBJ wrapper) and
%   the CoDoSol call block. The double lambda is derived from the adouble
%   text by word-ish-boundary adouble->double replacement, which also
%   retargets the _adouble interp suffixes (old-generator trick; duplicate
%   #defines across the two lambdas are identical, so the preprocessor
%   accepts the redefinitions).
import gdsge.codegen.cxx.fillTemplate
import gdsge.codegen.cxx.readTemplate
regions = ir.model.regions;
modelParts = cell(1, numel(regions));
callParts  = cell(1, numel(regions));
numEqs = zeros(1, numel(regions));
% suffix interp-ish names with _adouble ((?<=_) included: old parity)
suffixNames = [cellfun(@(s) s.name, ir.interp, 'UniformOutput', false), ...
    {'GDSGE_INTERP_VEC', 'GDSGE_INTERP_RSLT'}];
for k = 1:numel(regions)
    v = gdsge.ir.regionView(ir, k);
    body = gdsge.codegen.cxx.emitModelBody(v);
    for i = 1:numel(suffixNames)
        body = regexprep(body, ...
            ['((?<=^)|(?<=\W)|(?<=_))' suffixNames{i} '(?=\W)'], ...
            [suffixNames{i} '_adouble']);
    end
    [eqCode, numEqs(k)] = gdsge.codegen.cxx.emitEquations(v);
    lam = fillTemplate(readTemplate('model.tpl.cpp'), { ...
        'ARG_CODE',               gdsge.codegen.cxx.emitArgUnpack(ir); ...
        'DECLARE_CODE',           gdsge.codegen.cxx.emitDeclare(v); ...
        'MODEL_CODE',             body; ...
        'HEADER_AUX_ASSIGN_CODE', gdsge.codegen.cxx.emitAux(ir); ...
        'EQUATION_CODE',          eqCode});
    lamDouble = regexprep(lam, '((?<=^)|(?<=\W)|(?<=_))adouble(?=\W)', 'double');
    evalCode = fillTemplate(readTemplate('model_eval.tpl.cpp'), { ...
        'PRE_JAC_CODE',  ir.hooks.preJacCode; ...
        'POST_JAC_CODE', ir.hooks.postJacCode});
    part = [lam newline lamDouble newline evalCode];
    part = fillTemplate(part, {'NUM_EQUATIONS', num2str(numEqs(k))});
    % MODEL_NUMBER lives only in underscore-compound tokens (e.g.
    % GDSGE_FUNC_MODEL_NUMBER_adouble) — fillTemplate's word-boundary regex
    % cannot see it; strrep MODEL_NUMBER -> MODEL_k preserves the prefix.
    assert(contains(part, 'MODEL_NUMBER'), ...
        'MODEL_NUMBER not found in model code — template mismatch');
    modelParts{k} = strrep(part, 'MODEL_NUMBER', sprintf('MODEL_%d', k));

    cfm = fillTemplate(readTemplate('call_fmin.tpl.cpp'), { ...
        'MODEL_CONDITION', gdsge.codegen.cxx.lowerCondition(regions{k}.condition, ir); ...
        'NUM_EQUATIONS',   num2str(numEqs(k))});
    assert(contains(cfm, 'MODEL_NUMBER'), ...
        'MODEL_NUMBER not found in call_fmin code — template mismatch');
    callParts{k} = strrep(cfm, 'MODEL_NUMBER', sprintf('MODEL_%d', k));
end
assert(all(numEqs == numEqs(1)), ...
    'model regions must have equal equation counts (got %s)', mat2str(numEqs));
numEq = numEqs(1);
modelCode = strjoin(modelParts, newline);
callFMinCode = strjoin(callParts, newline);
end
