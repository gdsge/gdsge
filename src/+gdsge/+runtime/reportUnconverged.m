function msg = reportUnconverged(needResolved, f, probSize, stateNames, stateGrids, topK)
% REPORTUNCONVERGED  Human-readable report of unconverged grid points.
%   needResolved : logical row over problems (column-linear over probSize)
%   probSize     : [shock_num, state grid sizes...]
%   stateGrids   : cell of grid vectors aligned with stateNames
%   Returns the report text; the caller decides warning vs fprintf.
idx = find(needResolved);
if isempty(idx)
    msg = 'All problems converged.';
    return;
end
[~, order] = sort(f(idx), 'descend');
idx = idx(order);
k = min(topK, numel(idx));
lines = cell(1, k+1);
lines{1} = sprintf('%d of %d problems unconverged; worst %d:', ...
    nnz(needResolved), numel(needResolved), k);
subs = cell(1, numel(probSize));
for i = 1:k
    [subs{:}] = ind2sub(probSize, idx(i));
    parts = sprintf('shock=%d', subs{1});
    for s = 1:numel(stateNames)
        parts = sprintf('%s, %s=%g', parts, stateNames{s}, stateGrids{s}(subs{s+1}));
    end
    lines{i+1} = sprintf('  [%s] residual=%g', parts, f(idx(i)));
end
msg = strjoin(lines, newline);
end
