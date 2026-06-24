function [txt, numEq] = emitEquations(ir)
% EMITEQUATIONS  EQUATION_CODE: residual assignments GDSGE_EQ[i] = expr.
%   primed equations expand across the shock loop (old indexing:
%   GDSGE_EQ[<base-1>+GDSGE_iter], GDSGE_iter 1..shock_num). Asserts the
%   square system against the policy slot total.
V = gdsge.codegen.cxx.modelVars(ir);
n = ir.shocks.count;
w = gdsge.codegen.codeWriter();
numEq = 0;
body = gdsge.codegen.cxx.modelBody(ir);
for i = 1:numel(body.equations)
    eq = body.equations{i};
    if isfield(eq, 'kind') && strcmp(eq.kind, 'conditional')
        numEq = emitConditional(w, eq, ir, V, n, numEq);
    else
        numEq = emitPlain(w, eq, V, n, numEq);
    end
end
pol = ir.variables.policy;
numPolicyTotal = sum(cellfun(@(p) p.length, pol));
if numEq ~= numPolicyTotal
    error('gdsge:codegen:notSquare', ...
        '%d equations vs %d policy unknowns', numEq, numPolicyTotal);
end
txt = w.str();
end

% A single residual slot (or a primed expansion over the shock loop).
function numEq = emitPlain(w, eq, V, n, numEq)
import gdsge.codegen.cxx.emitExpr
if eq.primed
    w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++)', n);
    w.add('{');
    % 0-based indexing with a 1-based loop var: next slot is numEq, so
    % emit GDSGE_EQ[numEq-1+GDSGE_iter] (HL1996: GDSGE_EQ[10+GDSGE_iter])
    w.add('GDSGE_EQ[%d+GDSGE_iter]=%s;', numEq - 1, emitExpr(eq.expr, V.ctxLoop));
    w.add('}');
    numEq = numEq + n;
else
    w.add('GDSGE_EQ[%d]=%s;', numEq, emitExpr(eq.expr, V.ctx));
    numEq = numEq + 1;
end
end

% A runtime branch on a shock guard; each case fills the SAME slot range, so
% the cases must have equal slot arity (CaoNie2016: if A==AGood; Xp(1); else;
% consis1; end).
function numEq = emitConditional(w, eq, ir, V, n, numEq)
base = numEq; arity = [];
for c = 1:numel(eq.cases)
    cse = eq.cases{c};
    if c == 1
        w.add('if (%s) {', gdsge.codegen.cxx.lowerCondition(cse.cond, ir));
    elseif isempty(strtrim(cse.cond))
        w.add('} else {');
    else
        w.add('} else if (%s) {', gdsge.codegen.cxx.lowerCondition(cse.cond, ir));
    end
    nq = base;
    for j = 1:numel(cse.equations)
        nq = emitPlain(w, cse.equations{j}, V, n, nq);
    end
    if isempty(arity); arity = nq - base;
    else
        assert(nq - base == arity, ...
            'conditional equation cases have unequal slot counts (%d vs %d)', nq - base, arity);
    end
end
w.add('}');
numEq = base + arity;
end
