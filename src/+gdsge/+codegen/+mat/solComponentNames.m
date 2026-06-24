function names = solComponentNames(ir)
% SOLCOMPONENTNAMES  Per-row labels for the GDSGE_SOL/LB/UB solution stack.
%   Returns a 1×maxDim cell aligned with the policy slot layout: a scalar
%   policy var contributes its name; an array var of length L contributes
%   name(1)..name(L). Ordered by slot so it matches GDSGE_SOL's row order.
pol = ir.variables.policy;
maxDim = 0;
for i = 1:numel(pol); maxDim = max(maxDim, pol{i}.slot(2)); end
names = repmat({''}, 1, maxDim);
for i = 1:numel(pol)
    p = pol{i};
    lo = p.slot(1); hi = p.slot(2);
    if p.length == 1
        names{lo} = p.name;
    else
        for k = lo:hi
            names{k} = sprintf('%s(%d)', p.name, k - lo + 1);
        end
    end
end
end
