function d = constructAuxGlobalDer(values, breaks, sizeState, interpOrder, extrapOrder, wrtDim, numThreads)
% CONSTRUCTAUXGLOBALDER  Global var_aux construction (GDSGE_DER operator).
%   Fit the model-specification spline (same interp/extrap orders as
%   var_interp) to VALUES = [shock_num, prod(sizeState)] and return its
%   analytic partial derivative along state dim WRTDIM (1-based), evaluated at
%   every grid point — same [shock_num, prod(sizeState)] shape, ready for
%   var_interp update expressions and IterRslt assembly.
[ppCell, ~] = gdsge.runtime.constructSplines({values}, breaks, sizeState, ...
    interpOrder, extrapOrder, numThreads);
nDim = numel(breaks);
gridTensors = cell(1, nDim);
[gridTensors{:}] = ndgrid(breaks{:});
nSites = numel(gridTensors{1});
sites = zeros(nDim, nSites);
for j = 1:nDim
    sites(j, :) = gridTensors{j}(:)';
end
shockNum = size(values, 1);
idx = int32((0:shockNum-1)');
d = interp_der_mex(int32(numThreads), ppCell{1}, sites, idx, int32(wrtDim));
end
