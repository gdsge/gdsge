function L = dataLayout(ir)
% DATALAYOUT  Ordered descriptor of the packed GDSGE_DATA vector — THE single
%   source of truth shared by the MATLAB packer (gdsge.codegen.mat.emitDataPack)
%   and the C++ POP sequence (gdsge.codegen.cxx.emitPop). Order: shock count,
%   params (decl order), transition matrices (field order), shock values
%   (decl order) — constant per problem — then shockIdx and state values
%   (per-problem rows).
e = entry('shock_num', 'shockCount', 1, false);
for i = 1:numel(ir.params)
    e(end+1) = entry(ir.params{i}.name, 'param', numel(ir.params{i}.value), false); %#ok<AGROW>
end
tn = fieldnames(ir.shocks.transitions);
for i = 1:numel(tn)
    e(end+1) = entry(tn{i}, 'trans', numel(ir.shocks.transitions.(tn{i})), false); %#ok<AGROW>
end
for i = 1:numel(ir.shocks.names)
    nm = ir.shocks.names{i};
    e(end+1) = entry(nm, 'shockVals', numel(ir.shocks.values.(nm)), false); %#ok<AGROW>
end
e(end+1) = entry('shockIdx', 'shockIdx', 1, true);
for i = 1:numel(ir.states.names)
    % States are always scalar per-problem grid values (each state contributes
    % exactly one row to GDSGE_DATA). If a vector-valued state is ever added,
    % revisit this entry and the tensor packing in emitDataPack.
    e(end+1) = entry(ir.states.names{i}, 'state', 1, true); %#ok<AGROW>
end
L = struct('entries', {e}, 'nData', sum([e.length]));
end

function s = entry(name, role, len, perProb)
s = struct('name', name, 'role', role, 'length', len, 'perProb', perProb);
end
