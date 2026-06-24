function rows = addInterpCall(w, st, ir, rows, perShockBases)
% ADDINTERPCALL  Emit per-shock spline evaluation (value + per-dim derivative)
%   and register each interp result's chain row. Handles GDSGE_INTERP_VEC' (all
%   interp results) and named scalar interp calls (one named result), over any
%   state dimension. The interp result's gradient row is the chain rule
%   sum_dim interp_grad[dim] * d(arg_dim)/d(unknowns), assembled via combine.
import gdsge.codegen.cxx.sympymodel.combine
n = ir.shocks.count;
numInterp = numel(ir.interp);
interpNames = cellfun(@(it) it.name, ir.interp, 'UniformOutput', false);
xdim = numel(st.args);

% target t -> interp index in GDSGE_SPLINE_VEC
isVec = strcmp(st.interpRef, 'GDSGE_INTERP_VEC');
if isVec
    assert(numel(st.targets) == numInterp, ...
        'gdsge:codegen:sympyInterpTargets', ...
        'GDSGE_INTERP_VEC call has %d targets but %d interp vars', numel(st.targets), numInterp);
    idxFixed = 0;
else
    assert(numel(st.targets) == 1, ...
        'gdsge:codegen:sympyInterpTargets', 'named interp %s must have one target', st.interpRef);
    ni = find(strcmp(interpNames, st.interpRef), 1);
    assert(~isempty(ni), 'gdsge:codegen:sympyInterpRef', 'unknown interp %s', st.interpRef);
    idxFixed = ni - 1;
end

% per-arg site value expression + registry symbol/row
siteExprs = cell(1, xdim); argSyms = cell(1, xdim); argRows = cell(1, xdim);
for d = 1:xdim
    [siteExprs{d}, argSyms{d}] = siteInfo(st.args{d}, ir, perShockBases);
    assert(isKey(rows, argSyms{d}), ...
        'gdsge:codegen:sympyInterpArg', 'interp arg symbol %s is not registered', argSyms{d});
    argRows{d} = rows(argSyms{d});
end

% value arrays are pre-declared in generate (a name may also be a primed assign);
% here we declare only the per-dim derivative arrays.
for t = 1:numel(st.targets)
    for d = 1:xdim
        w.add('double %s_DER%d[%d];', st.targets{t}, d-1, n);
    end
end
w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++) {', n);
w.add('double GDSGE_rslt[%d]; double GDSGE_grad[%d];', numInterp, numInterp*xdim);
w.add('GDSGE_INTERP_VEC_double_grad(GDSGE_iter, %s, GDSGE_rslt, GDSGE_grad);', strjoin(siteExprs, ', '));
for t = 1:numel(st.targets)
    if isVec; ix = t - 1; else; ix = idxFixed; end
    w.add('%s_GDSGE_GRID[GDSGE_iter-1]=GDSGE_rslt[%d];', st.targets{t}, ix);
    for d = 1:xdim
        w.add('%s_DER%d[GDSGE_iter-1]=GDSGE_grad[%d];', st.targets{t}, d-1, ix + numInterp*(d-1));
    end
end
w.add('}');

% register each target's chain row: sum_dim DER<dim> * (arg_dim's row)
for t = 1:numel(st.targets)
    partials = containers.Map('KeyType', 'char', 'ValueType', 'char');
    argrowMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for d = 1:xdim
        partials(argSyms{d}) = sprintf('%s_DER%d[GDSGE_iter-1]', st.targets{t}, d-1);
        argrowMap(argSyms{d}) = argRows{d};
    end
    rows([st.targets{t} '__prime']) = combine(partials, argrowMap);
end
end

% --------------------------------------------------------------------------
function [expr, sym] = siteInfo(node, ir, perShockBases)
% C++ value expression of an interp eval coordinate + its registry symbol name.
if strcmp(node.kind, 'primed')
    sym = [node.id '__prime'];
    base = node.id;
    if ismember(base, ir.shocks.names)
        expr = sprintf('%s_GRID(GDSGE_iter)', base);
    elseif ismember(base, perShockBases)
        expr = sprintf('%s_GDSGE_GRID[GDSGE_iter-1]', base);
    else
        error('gdsge:codegen:sympyInterpArg', ...
            'primed interp arg %s is neither a shock nor a per-shock array', sym);
    end
elseif strcmp(node.kind, 'name')
    sym = node.id;
    if ismember(node.id, perShockBases)
        expr = sprintf('%s_GDSGE_GRID[GDSGE_iter-1]', node.id);
    else
        expr = node.id;   % scalar in scope (current unknown / assign)
    end
else
    error('gdsge:codegen:sympyInterpArg', ...
        'interp arg must be a name or primed name (Phase 8)');
end
end
