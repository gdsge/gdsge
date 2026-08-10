function files = generateMatlab(ir, outDir)
% GENERATEMATLAB  IR -> iter_<model>.m / simulate_<model>.m written to outDir.
%   interpMethod 'asg' -> emitIterAsg + (emitSimulateInterp when
%   options.simuInterp==1, else emitSimulateAsg). Otherwise the spline path:
%   emitIter + (emitSimulateInterp when options.simuInterp==1, else emitSimulate).
if nargin < 2; outDir = pwd; end
gdsge.codegen.assertSupportedIR(ir);
files = struct();
files.iterFile = fullfile(outDir, ['iter_' ir.modelName '.m']);
simuInterp = isfield(ir.options, 'simuInterp') && ir.options.simuInterp == 1;
if simuInterp && ~strcmp(ir.options.interpMethod, 'asg') ...
        && ~gdsge.codegen.cxx.isSimuMexExpressible(ir)
    warning('gdsge:codegen:simuInterpFallback', ...
        ['SIMU_INTERP: the simulate block is not MEX-expressible (a var_simu ' ...
         'entry outside var_output, or a non-name transition RHS); simulate_%s ' ...
         'uses the slower per-period MATLAB interpolation loop.'], ir.modelName);
end
if strcmp(ir.options.interpMethod, 'asg')
    gdsge.codegen.writeText(files.iterFile, gdsge.codegen.mat.emitIterAsg(ir));
    if simuInterp
        simTxt = gdsge.codegen.mat.emitSimulateInterp(ir);
    else
        simTxt = gdsge.codegen.mat.emitSimulateAsg(ir);
    end
else
    gdsge.codegen.writeText(files.iterFile, gdsge.codegen.mat.emitIter(ir));
    if simuInterp
        simTxt = gdsge.codegen.mat.emitSimulateInterp(ir);
    else
        simTxt = gdsge.codegen.mat.emitSimulate(ir);
    end
end
files.simulateFile = fullfile(outDir, ['simulate_' ir.modelName '.m']);
gdsge.codegen.writeText(files.simulateFile, simTxt);
end
