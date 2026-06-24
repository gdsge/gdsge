function txt = emitResultIter(ir)
% EMITRESULTITER  Output construction + explicit IterRslt assembly.
%   Frozen result shape: same fields as the old generator, including the
%   GDSGE_EMPTY placeholder fields v2struct used to create (params, var_aux,
%   var_tensor, pp, var_others). No v2struct anywhere.
lenOf = struct();
for i = 1:numel(ir.variables.policy)
    lenOf.(ir.variables.policy{i}.name) = ir.variables.policy{i}.length;
end
for i = 1:numel(ir.variables.aux)
    lenOf.(ir.variables.aux{i}.name) = ir.variables.aux{i}.length;
end
stateList = strjoin(ir.states.names, ',');            % 'w1'
gridCell = ['{' stateList '}'];                        % '{w1}'

w = gdsge.codegen.codeWriter();
w.add('if CONSTRUCT_OUTPUT==1');
w.add('%% Construct output variables');
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
w.blank();
w.add('IterRslt.shock_num = shock_num;');
w.add('IterRslt.shock_trans = shock_trans;');
for i = 1:numel(ir.params)
    w.add('IterRslt.params.%s = %s;', ir.params{i}.name, ir.params{i}.name);
end
w.add('IterRslt.params.GDSGE_EMPTY = GDSGE_EMPTY;');
for i = 1:numel(ir.shocks.names)
    w.add('IterRslt.var_shock.%s = %s;', ir.shocks.names{i}, ir.shocks.names{i});
end
w.add('GDSGE_VAR_POLICY = struct();');
for i = 1:numel(ir.variables.policy)
    n = ir.variables.policy{i}.name;
    w.add('GDSGE_VAR_POLICY.%s = %s;', n, n);
end
w.add('IterRslt.var_policy = GDSGE_VAR_POLICY;');
w.add('GDSGE_VAR_AUX = struct();');
for i = 1:numel(ir.variables.aux)
    n = ir.variables.aux{i}.name;
    w.add('GDSGE_VAR_AUX.%s = %s;', n, n);
end
w.add('GDSGE_VAR_AUX.GDSGE_EMPTY = GDSGE_EMPTY;');
w.add('IterRslt.var_aux = GDSGE_VAR_AUX;');
w.add('IterRslt.var_tensor.GDSGE_EMPTY = GDSGE_EMPTY;');
for i = 1:numel(ir.variables.interp)
    n = ir.variables.interp{i};
    w.add('IterRslt.pp.GDSGE_PP_%s = GDSGE_PP_%s;', n, n);
end
w.add('IterRslt.pp.GDSGE_SPLINE_VEC = GDSGE_SPLINE_VEC;');
w.add('IterRslt.pp.GDSGE_EMPTY = GDSGE_EMPTY;');
w.add('end');
w.add('IterRslt.Metric0 = GDSGE_Metric0;');
w.add('IterRslt.Metric = GDSGE_Metric;');
w.add('IterRslt.Iter = GDSGE_Iter;');
for i = 1:numel(ir.states.names)
    w.add('IterRslt.var_state.%s = %s;', ir.states.names{i}, ir.states.names{i});
end
for i = 1:numel(ir.variables.interp)
    n = ir.variables.interp{i};
    w.add('IterRslt.var_interp.%s = %s;', n, n);
end
w.add('IterRslt.GDSGE_PROB.GDSGE_LB = GDSGE_LB;');
w.add('IterRslt.GDSGE_PROB.GDSGE_UB = GDSGE_UB;');
w.add('IterRslt.GDSGE_PROB.GDSGE_SOL = GDSGE_SOL;');
w.add('IterRslt.GDSGE_PROB.GDSGE_F = GDSGE_F;');
w.add('IterRslt.GDSGE_PROB.GDSGE_SIZE = GDSGE_SIZE;');
for i = 1:numel(ir.variables.others)
    nm = ir.variables.others{i};
    w.add('IterRslt.var_others.%s = %s;', nm, nm);
end
w.add('IterRslt.var_others.GDSGE_EMPTY = GDSGE_EMPTY;');
w.add('IterRslt.NeedResolved = NeedResolved;');
txt = w.str();
end
