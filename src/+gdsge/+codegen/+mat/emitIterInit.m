function txt = emitIterInit(ir)
% EMITITERINIT  The model-init (last-period) solve segment of iter_<model>.m.
%   Returns '' when the model has no model_init. Old iter_init_template.m
%   parity: gated ONLY by SkipModelInit (runs even under WarmUp; the warm-up
%   interpolants supersede its output later), allocates the init problem over
%   the SAME state tensor, solves task_init, and leaves the init policy/aux
%   values in the workspace as rows (consumed by the interp initial exprs).
%   Runs BEFORE the main bounds/space section, which then re-allocates every
%   GDSGE_* array it touches.
if ~isfield(ir, 'modelInit'); txt = ''; return; end
iv     = gdsge.codegen.initView(ir);
bounds = gdsge.codegen.mat.emitBounds(iv);
pack   = gdsge.codegen.mat.emitDataPack(ir);   % layout identical to the main task
m      = ir.modelName;
nAux   = 0;
for i = 1:numel(iv.variables.aux); nAux = max(nAux, iv.variables.aux{i}.slot(2)); end
stateList     = strjoin(ir.states.names, ',');
stateNameCell = ['{''' strjoin(ir.states.names, ''',''') '''}'];
stateGridCell = ['{' stateList '}'];

w = gdsge.codegen.codeWriter();
w.add('%%%% Model-init solve (last-period problem; old semantics: SkipModelInit only)');
w.add('if ~SkipModelInit');
w.in();
w.addRaw(indent4(bounds.init));
w.add('GDSGE_EQVAL = 1e20*ones(%d,GDSGE_NPROB);', bounds.maxDim);
w.add('GDSGE_F = 1e20*ones(1,GDSGE_NPROB);');
w.add('GDSGE_SOL = zeros(%d,GDSGE_NPROB);', bounds.maxDim);
w.add('GDSGE_X0 = rand(size(GDSGE_SOL)) .* (GDSGE_UB-GDSGE_LB) + GDSGE_LB;');
w.add('GDSGE_SOL(:) = GDSGE_X0;');
w.add('GDSGE_AUX = zeros(%d,GDSGE_NPROB);', max(nAux, 1));
w.add('GDSGE_DATA = zeros(%d,GDSGE_NPROB);', pack.maxData);
w.add('%s', pack.iterPack);
for cfgLine = gdsge.codegen.mat.emitCfgCommon(); w.add('%s', cfgLine{1}); end
w.add('GDSGE_CFG.useBroydenNow = 0;');
w.add('GDSGE_CFG.taskName = MEX_TASK_INIT;');
w.add('GDSGE_CFG.splineVec = [];');
w.add('GDSGE_CFG.ppNames = {}; GDSGE_CFG.ppCell = {};');
w.add('GDSGE_CFG.maxMinorIter = MaxMinorIter;');
w.add('GDSGE_CFG.probSize = GDSGE_SIZE;');
w.add('GDSGE_CFG.useNearestNeighbor = true;');
w.add('GDSGE_CFG.verboseRetry = ~NoPrint;');
w.add('GDSGE_CFG.resolvePrintFreq = ResolvePrintFreq;');
w.add('GDSGE_CFG.majorIter = 0;');   % model_init runs before the iteration loop
w.add('GDSGE_CFG.adaptInSol = [];');
w.add('[GDSGE_SOL,GDSGE_F,GDSGE_AUX,GDSGE_EQVAL,GDSGE_OPT_INFO,GDSGE_DIAG] = gdsge.runtime.solveProblems(@mex_%s, GDSGE_SOL, GDSGE_LB, GDSGE_UB, GDSGE_DATA, GDSGE_F, GDSGE_AUX, GDSGE_EQVAL, GDSGE_CFG); %%#ok<ASGLU>', m);
w.add('if any(GDSGE_DIAG.needResolved)');
w.add('    warning(''gdsge:runtime:unconverged'', ''model_init: %%s'', gdsge.runtime.reportUnconverged(GDSGE_DIAG.needResolved, GDSGE_F, GDSGE_SIZE, %s, %s, 5));', ...
    stateNameCell, stateGridCell);
w.add('end');
for i = 1:numel(iv.variables.policy)
    p = iv.variables.policy{i};
    w.add('%s = GDSGE_SOL(%d:%d,:);', p.name, p.slot(1), p.slot(2));
end
for i = 1:numel(iv.variables.aux)
    a = iv.variables.aux{i};
    w.add('%s = GDSGE_AUX(%d:%d,:);', a.name, a.slot(1), a.slot(2));
end
w.out();
w.add('end');
txt = w.str();
end

function s = indent4(s)
s = gdsge.codegen.mat.indentBy(s, '    ');
end
