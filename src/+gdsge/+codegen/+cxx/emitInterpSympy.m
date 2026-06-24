function frag = emitInterpSympy(ir)
% EMITINTERPSYMPY  Spline interp section for the SymPy backend: a double
%   evaluator returning value AND gradient (search_eval_with_grad_vec_at_array),
%   over GDSGE_SPLINE_VEC (which bundles ALL interp vars, named or VEC). Works
%   for any state dimension xdim; addInterpCall picks indices/dims per call.
%   frag.getCode (task scope: constructs GDSGE_CSPLINE_VEC), frag.threadCode
%   (per-thread scratch + the GDSGE_INTERP_VEC_double_grad lambda).
import gdsge.codegen.cxx.fillTemplate
import gdsge.codegen.cxx.readTemplate
numInterp = numel(ir.interp);
if numInterp == 0
    frag = struct('getCode', '', 'threadCode', '');
    return;
end
xdim = numel(ir.interp{1}.args);
assert(all(cellfun(@(it) numel(it.args) == xdim, ir.interp)), ...
    'gdsge:codegen:sympyMixedXdim', 'all interp vars must share arg dimension');

siteParams = strjoin(arrayfun(@(d) sprintf('double GDSGE_site%d', d-1), 1:xdim, ...
    'UniformOutput', false), ', ');
siteList = strjoin(arrayfun(@(d) sprintf('GDSGE_site%d', d-1), 1:xdim, ...
    'UniformOutput', false), ',');

if strcmp(ir.options.interpMethod, 'asg')
    frag = asgFrag(xdim, numInterp, siteParams, siteList);
    return;
end

getCode = fillTemplate(readTemplate('interp_spline_construct.tpl.cpp'), ...
    {'VAR_NUM', num2str(xdim); 'NUM_INTERP', num2str(numInterp)});
w = gdsge.codegen.codeWriter();
w.add('int GDSGE_INTERP_CELL[%d] = {0};', xdim);
w.add('auto GDSGE_INTERP_VEC_double_grad = [&GDSGE_CSPLINE_VEC,&GDSGE_INTERP_CELL](int shockIdx, %s, double* GDSGE_out, double* GDSGE_grad){', siteParams);
w.add('  double xSite[] = {%s};', siteList);
w.add('  GDSGE_CSPLINE_VEC.search_eval_with_grad_vec_at_array(shockIdx-1, xSite, GDSGE_out, GDSGE_grad, GDSGE_INTERP_CELL);');
w.add('};');
frag = struct('getCode', getCode, 'threadCode', w.str());
end

function frag = asgFrag(xdim, numInterp, siteParams, siteList)
% SymPy ASG interp: reuse the autodiff ASG construct (declares GDSGE_CPP_ASG),
% emit a double-only per-thread scratch block + the shared-name double-grad
% lambda backed by eval_vec_with_grad, transposing the native [v*DIM+d]
% gradient to the spline layout [v + NVEC*d] that addInterpCall consumes.
% ASG_MAX_LEVEL stays a compile-time define; DIM/NVEC are substituted here.
import gdsge.codegen.cxx.readTemplate
getCode = readTemplate('interp_asg_construct.tpl.cpp');
w = gdsge.codegen.codeWriter();
w.add('double GDSGE_ASG_CELL[(ASG_MAX_LEVEL+2)*%d] = {0};', xdim);
w.add('double GDSGE_ASG_RATIO[(ASG_MAX_LEVEL+2)*%d] = {0};', xdim);
w.add('double GDSGE_ASG_SLOPE[(ASG_MAX_LEVEL+2)*%d] = {0};', xdim);
w.add('double GDSGE_ASG_GRAD_NATIVE[%d*%d] = {0};', numInterp, xdim);
w.add('auto GDSGE_INTERP_VEC_double_grad = [&GDSGE_CPP_ASG,&GDSGE_ASG_CELL,&GDSGE_ASG_RATIO,&GDSGE_ASG_SLOPE,&GDSGE_ASG_GRAD_NATIVE](int shockIdx, %s, double* GDSGE_out, double* GDSGE_grad){', siteParams);
w.add('  double xSite[] = {%s};', siteList);
w.add('  GDSGE_CPP_ASG.eval_vec_with_grad(shockIdx-1, xSite, GDSGE_ASG_CELL, GDSGE_ASG_RATIO, GDSGE_ASG_SLOPE, GDSGE_out, GDSGE_ASG_GRAD_NATIVE);');
w.add('  for(int GDSGE_v=0; GDSGE_v<%d; GDSGE_v++) for(int GDSGE_d=0; GDSGE_d<%d; GDSGE_d++)', numInterp, xdim);
w.add('    GDSGE_grad[GDSGE_v + %d*GDSGE_d] = GDSGE_ASG_GRAD_NATIVE[GDSGE_v*%d + GDSGE_d];', numInterp, xdim);
w.add('};');
frag = struct('getCode', getCode, 'threadCode', w.str());
end
