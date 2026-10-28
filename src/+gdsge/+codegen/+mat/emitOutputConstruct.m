function txt = emitOutputConstruct(ir)
% EMITOUTPUTCONSTRUCT  Build IterRslt.output_interp (+ output_var_index) from
%   the var_output values currently in the iter workspace. Used by the final
%   IterRslt assembly (emitResultIter, under CONSTRUCT_OUTPUT) and by the old
%   OUTPUT_CONSTRUCT_CODE placeholder inside a pre/post_iter hook.
lenOf = struct();
for i = 1:numel(ir.variables.policy)
    lenOf.(ir.variables.policy{i}.name) = ir.variables.policy{i}.length;
end
for i = 1:numel(ir.variables.aux)
    lenOf.(ir.variables.aux{i}.name) = ir.variables.aux{i}.length;
end
auxGlobal = {};
if isfield(ir.variables, 'auxGlobal'); auxGlobal = ir.variables.auxGlobal; end
for i = 1:numel(auxGlobal)
    lenOf.(auxGlobal{i}.name) = 1;   % globals are always length 1
end
stateList = strjoin(ir.states.names, ',');
gridCell = ['{' stateList '}'];
w = gdsge.codegen.codeWriter();
parts = cell(1, numel(ir.variables.output));
for i = 1:numel(ir.variables.output)
    n = ir.variables.output{i};
    parts{i} = sprintf('reshape(%s,%d,[])', n, lenOf.(n));
end
w.add('outputVarStack = cat(1,%s);', strjoin(parts, ','));
w.add('GDSGE_NCOMP = size(outputVarStack,1);');
% Stacked uniform-order layout: shock is the stacked vector index (not a spline
% dim). For realized shock s, output row r is vector index (s-1)*nComp + r.
% Values: [shock_num*nComp, stateDims], comp fastest. Built (and evaluated in
% simulate) via interp_construct_mex / interp_eval_mex — myppual is retired.
w.add('if shock_num>1');
w.add('    outputVarStack = reshape(outputVarStack, GDSGE_NCOMP, shock_num, GDSGE_SIZE_STATE{:});');
w.add('    GDSGE_OUTPUT_VALUES = reshape(outputVarStack, GDSGE_NCOMP*shock_num, GDSGE_SIZE_STATE{:});');
w.add('else');
w.add('    GDSGE_OUTPUT_VALUES = reshape(outputVarStack, GDSGE_NCOMP, GDSGE_SIZE_STATE{:});');
w.add('end');
w.add('[GDSGE_OUTPUT_PP_CELL, ~] = interp_construct_mex({%s}, struct(''GDSGE_V1'', GDSGE_OUTPUT_VALUES), ...', stateList);
w.add('    int32(OutputInterpOrder*ones(1,length(%s))), int32([]), int32(NumThreads));', gridCell);
w.add('IterRslt.output_interp = GDSGE_OUTPUT_PP_CELL{1};');
w.add('IterRslt.output_interp.GDSGE_NCOMP = GDSGE_NCOMP;');
w.add('output_var_index=struct();');
idx = 1;
for i = 1:numel(ir.variables.output)
    n = ir.variables.output{i};
    w.add('output_var_index.%s=%d:%d;', n, idx, idx + lenOf.(n) - 1);
    idx = idx + lenOf.(n);
end
w.add('IterRslt.output_var_index = output_var_index;');
txt = w.str();
end
