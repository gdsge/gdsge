function code = simuHook(ir, which)
% SIMUHOOK  The gmod's pre_simulate_iter / post_simulate_iter MATLAB text
%   (WHICH = 'pre' | 'post'), '' when not declared. Spliced verbatim into
%   every simulate variant's period loop (old SIMU_PRE/POST_ITER_CODE).
f = [which 'SimulateIter'];
code = '';
if isfield(ir, 'hooks') && isfield(ir.hooks, f); code = ir.hooks.(f); end
end
