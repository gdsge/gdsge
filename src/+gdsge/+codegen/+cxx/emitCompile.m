function txt = emitCompile(ir, includeDir, withSimuMex)
% EMITCOMPILE  compile_<model>.m from the cleaned template.
%
%   withSimuMex (default false): also compile simulate_<model>_mex.cpp (the
%   whole-loop SIMU_INTERP simulate MEX) with the same OpenMP flags. It needs
%   only -DMAXDIM/-DINTERP_ORDER (interp_eval_double.h is header-only).
%
%   extraDef: empty when ir.modelInit is absent; '-DHAS_INIT' when it is
%   present.  GDSGE_MAXDIM = max(main policy total, init policy total) + 4,
%   matching the old generator's parity.  GDSGE_SOLVE_DIM = exact solver
%   dimension (no +4) = max(main policy total, init policy total); this sizes
%   solver working-storage to exactly n (no +4 slack).  Under interpMethod='asg'
%   the ASG defines (-DUSE_ASG -DASG_MAX_DIM/NVEC/LEVEL, from asg.get_mex_constants)
%   are appended.  Finite-diff and sparse defines are not yet wired (Phase 7c work).
%
%   MODEL_NAME appears in compound tokens (compile_MODEL_NAME, mex_MODEL_NAME.cpp)
%   so it is substituted with strrep (boundary-less); all other keys are
%   standalone and are filled via fillTemplate.
import gdsge.codegen.cxx.fillTemplate
import gdsge.codegen.cxx.readTemplate
if nargin < 3; withSimuMex = false; end
pol = ir.variables.policy;
solveDim = sum(cellfun(@(p) p.length, pol));
maxDim = solveDim + 4;
extraDef = '';
if isfield(ir, 'modelInit')
    initDim = sum(cellfun(@(p) p.length, ir.modelInit.variables.policyInit)) + 4;
    maxDim = max(maxDim, initDim);     % old: max(num_policy_total, num_policy_init_total)+4
    solveDim = max(solveDim, sum(cellfun(@(p) p.length, ir.modelInit.variables.policyInit)));
    extraDef = '-DHAS_INIT';
end
if strcmp(ir.options.interpMethod, 'asg')
    [asgMaxDim, asgMaxNvec, asgMaxLevel] = asg.get_mex_constants();
    if numel(ir.states.names) > asgMaxDim
        error('gdsge:codegen:asgDim', ...
            'model has %d states but asg_mex was built with ASG_MAX_DIM=%d', ...
            numel(ir.states.names), asgMaxDim);
    end
    if ir.options.asgMaxLevel > asgMaxLevel
        error('gdsge:codegen:asgLevel', ...
            'AsgMaxLevel=%d exceeds asg_mex ASG_MAX_LEVEL=%d', ...
            ir.options.asgMaxLevel, asgMaxLevel);
    end
    extraDef = strtrim([extraDef ' -DUSE_ASG' ...
        ' -DASG_MAX_DIM='   num2str(asgMaxDim) ...
        ' -DASG_MAX_NVEC='  num2str(asgMaxNvec) ...
        ' -DASG_MAX_LEVEL=' num2str(asgMaxLevel)]);
end
txt = readTemplate('compile.tpl.m');
% MODEL_NAME is underscore-adjacent in compile_MODEL_NAME and mex_MODEL_NAME.cpp
% so fillTemplate can't see it — use strrep.
assert(contains(txt, 'MODEL_NAME'), ...
    'MODEL_NAME not found in compile template — template mismatch');
txt = strrep(txt, 'MODEL_NAME', ir.modelName);
% Second mex call for the whole-loop simulate MEX (numbers baked, no placeholders).
if withSimuMex
    simuMexCompile = [ ...
        'fprintf(''      Compiling simulate MEX (simulate_' ir.modelName '_mex) ...\n'');' newline ...
        'simu_cpp = fullfile(current_folder,''simulate_' ir.modelName '_mex.cpp'');' newline ...
        'simu_compile = [mexCommand '' -DMAXDIM=' num2str(maxDim) ...
        ' -DINTERP_ORDER=' num2str(ir.options.interpOrder) ''' flag2 flag3 ' ...
        'sprintf('' "%s"'',simu_cpp) link_to_lib '' -outdir "'' current_folder ''"'' '' -I"'' include_folder ''"''];' newline ...
        'eval(simu_compile);'];
else
    simuMexCompile = '';
end
% Remaining keys are standalone (appear at word boundaries or adjacent to = / ' / space).
txt = fillTemplate(txt, { ...
    'INCLUDE_FOLDER',      strrep(includeDir, '\', '/'); ...
    'GDSGE_MAXDIM',        num2str(maxDim); ...
    'GDSGE_SOLVE_DIM',     num2str(solveDim); ...
    'GDSGE_INTERP_ORDER',  num2str(ir.options.interpOrder); ...
    'EXTRA_DEF',           extraDef; ...
    'SIMU_MEX_COMPILE',    simuMexCompile});
end
