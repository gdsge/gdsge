function metric = computeMetric(newCell, oldCell)
% COMPUTEMETRIC  Sup-norm distance across interp variables.
%   NaN on any error (size mismatch etc.) — parity with the old generated
%   try/catch -> NaN. Assumes at least one interp variable.
try
    m = -inf;
    for i = 1:numel(newCell)
        m = max(m, max(abs(newCell{i}(:) - oldCell{i}(:))));
    end
    metric = m;
catch
    metric = nan;
end
end
