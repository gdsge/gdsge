function V = modelVars(ir)
% MODELVARS  Classify model-body variables for the C++ emitters.
%   scalarTemps: unprimed assign targets + reduction targets + unprimed
%       GDSGE_INTERP_VEC targets (declared `adouble x(0)`), sorted unique.
%   futureVars: primed interpCall targets + primed assign targets (declared
%       as per-shock `_GDSGE_GRID` arrays), sorted unique.
%   ctx/ctxLoop: ready-made emitExpr contexts (outside / inside shock loop).
body = gdsge.codegen.cxx.modelBody(ir);
[scalarTemps, futureVars] = collect(body.statements);
V = struct();
V.scalarTemps = reshape(unique(scalarTemps), 1, []);
V.futureVars  = reshape(unique(futureVars),  1, []);
V.ctx = struct('shockNames', {ir.shocks.names}, ...
               'futureVars', {V.futureVars}, 'inLoop', false);
V.ctxLoop = V.ctx; V.ctxLoop.inLoop = true;
end

function [scalarTemps, futureVars] = collect(stmts)
% Recurses into conditional statements (ADR-0003): a branch-assigned scalar is
% still declared once at the top of the model lambda, so the C++ compiles
% whichever branch runs and the value survives the join.
scalarTemps = {};
futureVars = {};
for i = 1:numel(stmts)
    st = stmts{i};
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
            else
                scalarTemps = [scalarTemps, st.targets]; %#ok<AGROW>
            end
        case 'conditional'
            for c = 1:numel(st.cases)
                [s2, f2] = collect(st.cases{c}.statements);
                scalarTemps = [scalarTemps, s2]; %#ok<AGROW>
                futureVars  = [futureVars,  f2]; %#ok<AGROW>
            end
    end
end
end
