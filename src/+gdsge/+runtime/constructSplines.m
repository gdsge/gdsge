function [ppCell, splineVec] = constructSplines(valueCell, breaks, sizeState, interpOrder, extrapOrder, numThreads)
% CONSTRUCTSPLINES  Build interpolants + GDSGE_SPLINE_VEC via the fused
%   interp_construct_mex kernel: it constructs all interp variables and writes
%   the eval-order array (GDSGE_SPLINE_VEC) in one call (no per-var dispatch,
%   no permute/concat).
%   valueCell : cell of value arrays (one per interp variable)
%   breaks    : cell of state grid vectors, e.g. {w1}
%   sizeState : num2cell of the state-space size, e.g. {201}
%   interpOrder, extrapOrder, numThreads: from the generated file's options.
orderVec = interpOrder * ones(1, numel(sizeState));
if interpOrder == 4
    extrapVec = extrapOrder * ones(1, numel(sizeState));
else
    extrapVec = [];
end
nv = numel(valueCell);
if nv == 0
    ppCell = {}; splineVec = struct(); return;
end
values = struct();
for i = 1:nv
    values.(sprintf('GDSGE_V%d', i)) = reshape(valueCell{i}, [], sizeState{:});
end
[ppCell, splineVec] = interp_construct_mex(breaks, values, ...
    int32(orderVec), int32(extrapVec), int32(numThreads));
end
