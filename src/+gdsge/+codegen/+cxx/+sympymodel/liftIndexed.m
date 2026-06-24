function node = liftIndexed(node, rows, w, perShockBases)
% LIFTINDEXED  Replace constant-indexed accesses to a per-shock array —
%   call{fn,[const]} where fn has a per-shock <fn>_GDSGE_GRID (interp/primed-
%   assign target or array policy) — with a lifted scalar symbol fn_GDSGE_AT_k
%   bound to fn_GDSGE_GRID[k-1], and register its gradient row = fn's per-shock
%   row specialized to shock k. Used for safe_assets `Re_n(2)` (Re_n' is a
%   per-shock assign indexed in scalar scope). rows doubles as the dedup set.
switch node.kind
    case 'call'
        if ismember(node.fn, perShockBases)
            assert(numel(node.args) == 1 && strcmp(node.args{1}.kind, 'num'), ...
                'gdsge:codegen:sympyIndexedAccess', ...
                'indexed access %s(...) requires a constant index', node.fn);
            k = node.args{1}.value;
            assert(k == floor(k), 'gdsge:codegen:sympyIndexedAccess', ...
                '%s index must be an integer', node.fn);
            sym = sprintf('%s_GDSGE_AT_%d', node.fn, k);
            if ~isKey(rows, sym)
                srcKey = node.fn;
                if ~isKey(rows, srcKey) && isKey(rows, [node.fn '__prime'])
                    srcKey = [node.fn '__prime'];   % interp targets register under name'
                end
                assert(isKey(rows, srcKey), 'gdsge:codegen:sympyIndexedAccess', ...
                    'indexed access %s(%d) before %s is computed', node.fn, k, node.fn);
                w.add('double %s = %s_GDSGE_GRID[%d];', sym, node.fn, k-1);
                rows(sym) = specializeRow(rows(srcKey), k);
            end
            node = gdsge.ir.node.name(sym);
        else
            for i = 1:numel(node.args)
                node.args{i} = gdsge.codegen.cxx.sympymodel.liftIndexed( ...
                    node.args{i}, rows, w, perShockBases);
            end
        end
    case 'binop'
        node.lhs = gdsge.codegen.cxx.sympymodel.liftIndexed(node.lhs, rows, w, perShockBases);
        node.rhs = gdsge.codegen.cxx.sympymodel.liftIndexed(node.rhs, rows, w, perShockBases);
    case 'unop'
        node.arg = gdsge.codegen.cxx.sympymodel.liftIndexed(node.arg, rows, w, perShockBases);
    % num, name, primed: leaves
end
end

function out = specializeRow(row, k)
% Resolve a per-shock row to a fixed shock k: per-shock array exprs index [k-1];
% templated (diagonal) slots resolve to base+(k-1), fixed.
out = struct('slot', {}, 'expr', {}, 'templated', {});
for i = 1:numel(row)
    e = row(i);
    expr = strrep(e.expr, 'GDSGE_iter-1', num2str(k-1));
    if isfield(e, 'templated') && e.templated
        out(end+1) = struct('slot', e.slot + (k-1), 'expr', expr, 'templated', false); %#ok<AGROW>
    else
        out(end+1) = struct('slot', e.slot, 'expr', expr, 'templated', false); %#ok<AGROW>
    end
end
end
