function b = modelBody(ir)
% MODELBODY  The {statements, equations} the per-body C++ emitters consume.
%   Accepts either a regionView (ir.model already has .statements/.equations)
%   or a full IR (ir.model has .regions). A full IR defaults to region 1 — the
%   single-region case used by callers/tests that hand a whole IR to a per-body
%   emitter directly. Multi-region drivers always pass a regionView, so they
%   pick the right region.
if isfield(ir.model, 'statements')
    b = ir.model;
else
    b = ir.model.regions{1};
end
end
