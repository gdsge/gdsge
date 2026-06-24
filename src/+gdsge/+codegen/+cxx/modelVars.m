function V = modelVars(ir)
% MODELVARS  Classify model-body variables for the C++ emitters.
%   scalarTemps: unprimed assign targets + reduction targets (declared
%       `adouble x(0)`), sorted unique.
%   futureVars: primed interpCall targets + primed assign targets (declared
%       as per-shock `_GDSGE_GRID` arrays), sorted unique.
%   ctx/ctxLoop: ready-made emitExpr contexts (outside / inside shock loop).
scalarTemps = {};
futureVars = {};
body = gdsge.codegen.cxx.modelBody(ir);
for i = 1:numel(body.statements)
    st = body.statements{i};
    switch st.type
        case 'assign'
            if st.primed
                futureVars{end+1} = st.target; %#ok<AGROW>
            else
                scalarTemps{end+1} = st.target; %#ok<AGROW>
            end
        case 'reduction'
            scalarTemps{end+1} = st.target; %#ok<AGROW>
        case 'interpCall'
            if st.primed
                futureVars = [futureVars, st.targets]; %#ok<AGROW>
            end
    end
end
V = struct();
V.scalarTemps = reshape(unique(scalarTemps), 1, []);
V.futureVars  = reshape(unique(futureVars),  1, []);
V.ctx = struct('shockNames', {ir.shocks.names}, ...
               'futureVars', {V.futureVars}, 'inLoop', false);
V.ctxLoop = V.ctx; V.ctxLoop.inLoop = true;
end
