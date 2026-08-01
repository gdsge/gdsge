function msg = reportUnconvergedSamples(needResolved, f, shockVec, stateNames, stateVals, topK)
% REPORTUNCONVERGEDSAMPLES  Human-readable report of unconverged simulate samples.
%   The simulate path solves one problem per sample, so a problem's coordinates
%   are that sample's own current shock and state values - not grid subscripts,
%   which is reportUnconverged's job. Everything here is indexed by the linear
%   problem (= sample) index.
%
%   needResolved : logical row over samples
%   f            : residual per sample
%   shockVec     : per-sample shock index
%   stateNames   : cell of state names
%   stateVals    : cell of per-sample state vectors, aligned with stateNames
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
for i = 1:k
    j = idx(i);
    parts = sprintf('sample=%d, shock=%d', j, shockVec(j));
    for s = 1:numel(stateNames)
        parts = sprintf('%s, %s=%g', parts, stateNames{s}, stateVals{s}(j));
    end
    lines{i+1} = sprintf('  [%s] residual=%g', parts, f(j));
end
msg = strjoin(lines, newline);
end
