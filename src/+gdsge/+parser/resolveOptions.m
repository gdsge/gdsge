function o = resolveOptions(ws)
% RESOLVEOPTIONS  Map the eval'd flag workspace to the curated IR options.
%   ASG level/threshold are emitted only for the asg method (matching the
%   reference IR, which omits them for spline).
if getf(ws, 'USE_ASG', 0)
    method = 'asg';
elseif getf(ws, 'USE_PCHIP', 0)
    method = 'pchip';
else
    method = 'spline';
end
o = struct();
o.interpMethod = method;
o.interpOrder  = getf(ws, 'INTERP_ORDER', 4);
o.extrapOrder  = getf(ws, 'EXTRAP_ORDER', 2);
if strcmp(method, 'asg')
    o.asgMinLevel        = getf(ws, 'AsgMinLevel', 4);
    o.asgMaxLevel        = getf(ws, 'AsgMaxLevel', 10);
    o.asgThreshold       = getf(ws, 'AsgThreshold', 1e-2);
    o.asgOutputMaxLevel  = getf(ws, 'AsgOutputMaxLevel', 10);
    o.asgOutputThreshold = getf(ws, 'AsgOutputThreshold', 1e-2);
end
o.tolEq       = getf(ws, 'TolEq', 1e-6);
o.numThreads  = getf(ws, 'NumThreads', 0);   % 0 = dynamic (feature('numcores') at runtime)
o.simuResolve   = getf(ws, 'SIMU_RESOLVE', 1);
o.simuInterp    = getf(ws, 'SIMU_INTERP', 0);
% spec §7: the two simulate modes are mutually exclusive; exactly one active.
if (o.simuResolve ~= 0) == (o.simuInterp ~= 0)
    error('gdsge:parser:simuModeConflict', ...
        ['SIMU_RESOLVE (=%g) and SIMU_INTERP (=%g) are mutually exclusive: ' ...
         'exactly one must be nonzero.'], o.simuResolve, o.simuInterp);
end
o.simuPrintFreq = getf(ws, 'SimuPrintFreq', 1000);
o.simuSaveFreq  = getf(ws, 'SimuSaveFreq', Inf);
o.printFreq     = getf(ws, 'PrintFreq', 10);
o.resolvePrintFreq = getf(ws, 'ResolvePrintFreq', 100);
o.saveFreq      = getf(ws, 'SaveFreq', 10);
% Jacobian backend: emit the field ONLY when the gmod sets UseAutoDiff (tri-state).
% Absent -> no field -> AUTO (resolved at codegen time by gdsge.codegen.resolveBackend),
% keeping default-path IR byte-identical to existing goldens. =1 -> autodiff, =0 -> sympy.
if isfield(ws, 'UseAutoDiff')
    if double(ws.UseAutoDiff) == 0
        o.jacobianBackend = 'sympy';
    else
        o.jacobianBackend = 'autodiff';
    end
end
end

function v = getf(ws, name, dflt)
if isfield(ws, name); v = double(ws.(name)); else; v = dflt; end
end
