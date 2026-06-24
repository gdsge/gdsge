function irv = initView(ir)
% INITVIEW  A main-model-shaped view of ir.modelInit so every existing emitter
%   (emitModel/emitModelBody/emitArgUnpack/emitDeclare/emitEquations/emitAux/
%   emitPop/emitInterp) is reused verbatim for the init task.
%   The init problem has no interpolants and no hook code, and its GDSGE_DATA
%   layout is identical to the main task's (dataLayout ignores policy/aux).
%
%   Invariant: emitters reachable from emitTask (the cxx task path) must NOT
%   read irv.variables.output, irv.variables.others, or irv.variables.tensor —
%   those fields are copied from the MAIN model and their values are
%   meaningless in the init context.  (tensor is separately refused by
%   generateCxx before initView is called.)
%
%   NOTE on emitBounds: emitBounds is listed in the reused-emitter set above
%   because it consumes this view on the MATLAB side (the iter-init segment),
%   not in the cxx task path.
irv = ir;
irv.variables.policy = ir.modelInit.variables.policyInit;
irv.variables.aux    = ir.modelInit.variables.auxInit;
irv.variables.interp = {};
irv.interp = {};
irv.bounds = ir.modelInit.bounds;
irv.model  = struct('regions', {{ struct('condition', '', ...
                    'statements', {ir.modelInit.statements}, ...
                    'equations',  {ir.modelInit.equations}) }});
irv.hooks  = struct('preModel','','preIter','','postIter','', ...
                    'preJacCode','','postJacCode','','cxx','');
irv = rmfield(irv, 'modelInit');
end
