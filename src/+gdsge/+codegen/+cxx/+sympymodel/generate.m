function g = generate(ir)
% GENERATE  IR -> {argCode, popCode, bodyCode, numEq, numAux} for the SymPy
%   analytic-Jacobian model function (all-double; value + JAC). Walks the model
%   statements maintaining a gradient registry (name -> gradRow), then emits the
%   residual values and the analytic Jacobian. Handles: scalar/array unknown
%   unpack, interp calls (chain rule via interp'), EXPECT/MIN/MAX/PROD
%   reductions, scalar + primed assigns, scalar + primed equations.
import gdsge.codegen.cxx.sympymodel.seedUnknownRow
import gdsge.codegen.cxx.sympymodel.combine
import gdsge.codegen.cxx.sympymodel.diffRHS
import gdsge.codegen.cxx.sympymodel.addReductionLoop
import gdsge.codegen.cxx.sympymodel.addInterpCall
import gdsge.codegen.cxx.sympymodel.addPrimedAssign
import gdsge.codegen.cxx.sympymodel.loadLoopLocals

n = ir.shocks.count;
NX = sum(cellfun(@(p) p.length, ir.variables.policy));
body = gdsge.codegen.cxx.modelBody(ir);
w = gdsge.codegen.codeWriter();

% per-shock array bases: interp/primed-assign targets (modelVars) + array policies
futureVars = gdsge.codegen.cxx.modelVars(ir).futureVars;
arrayPolicies = {};
for i = 1:numel(ir.variables.policy)
    if ir.variables.policy{i}.length > 1
        arrayPolicies{end+1} = ir.variables.policy{i}.name; %#ok<AGROW>
    end
end
perShockBases = [futureVars, arrayPolicies];

% ---- registry: name -> gradRow (struct array of mkEntry) ----------------
rows = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:numel(ir.variables.policy)
    p = ir.variables.policy{i};
    s0 = p.slot(1) - 1;
    if p.length == 1
        rows(p.name) = seedUnknownRow(s0, 1, 0);           % fixed seed
    else
        % array policy = per-shock future state, referenced as <name>' :
        % d(<name>'_j)/d(unknowns) = 1 at the j-th array slot (templated).
        rows([p.name '__prime']) = gdsge.codegen.cxx.sympymodel.mkEntry(s0, '1', true);
    end
end

% Declare every per-shock value array ONCE (a name may be both an interp target
% and a primed-assign target, e.g. Mendoza flow_future); emitters only write.
for i = 1:numel(futureVars)
    w.add('double %s_GDSGE_GRID[%d];', futureVars{i}, n);
end

% ---- statements ---------------------------------------------------------
% Track scalar-assign targets already declared, so a reassignment (e.g.
% CaoKS2016's `rhs1 = lambda1; ... rhs1 = rhs1 + EXPECT{...}`) emits a bare
% `name = ...` instead of redeclaring `double name` (C2374). The registry row
% is overwritten either way (correct forward-mode accumulation).
declaredAssigns = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for i = 1:numel(body.statements)
    st = body.statements{i};
    switch st.type
        case 'interpCall'
            assert(st.primed, 'sympy: unprimed interp call unsupported');
            rows = addInterpCall(w, st, ir, rows, perShockBases);
        case 'reduction'
            rows = addReductionLoop(w, st, ir, rows, perShockBases, NX, sprintf('hr%d_', i));
        case 'assign'
            if st.primed
                rows = addPrimedAssign(w, st, ir, rows, perShockBases, sprintf('hp%d_', i));
            else
                expr = gdsge.codegen.cxx.sympymodel.liftIndexed(st.expr, rows, w, perShockBases);
                d = diffRHS(expr, rows.keys, sprintf('ha%d_', i));
                emitHelpers(w, d.helpers);
                if isKey(declaredAssigns, st.target)
                    w.add('%s = %s;', st.target, d.value);          % reassignment
                else
                    w.add('double %s = %s;', st.target, d.value);   % first declaration
                    declaredAssigns(st.target) = true;
                end
                rows(st.target) = combine(d.partialsMap, rows);
            end
        otherwise
            error('gdsge:codegen:sympyUnsupportedStmt', 'statement %s', st.type);
    end
end

% ---- equations: values, then (guarded) Jacobian ------------------------
w.add('if (GDSGE_jac) { for(int GDSGE_ji=0; GDSGE_ji<%d*%d; ++GDSGE_ji) GDSGE_jac[GDSGE_ji]=0.0; }', NX, NX);
numEq = 0;
for i = 1:numel(body.equations)
    eq = body.equations{i};
    if isfield(eq, 'kind') && strcmp(eq.kind, 'conditional')
        % Phase 9b: each case fills the SAME slot range under a shock guard;
        % value + Jacobian rows are emitted per branch (CaoNie2016 if/else).
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
                nq = addPlainEq(w, cse.equations{j}, ir, rows, perShockBases, n, nq, ...
                    sprintf('%dc%dq%d', i, c, j));
            end
            if isempty(arity); arity = nq - base;
            else; assert(nq - base == arity, 'sympy conditional cases unequal slot counts'); end
        end
        w.add('}');
        numEq = base + arity;
    else
        numEq = addPlainEq(w, eq, ir, rows, perShockBases, n, numEq, sprintf('%d', i));
    end
end

g = struct();
g.numEq = numEq;
if isempty(ir.variables.aux)
    g.numAux = 0;
else
    g.numAux = sum(cellfun(@(a) a.length, ir.variables.aux));
end
g.bodyCode = w.str();
g.argCode = argUnpackDouble(ir);
g.popCode = gdsge.codegen.cxx.emitPop(ir);
assert(numEq == NX, 'sympy: %d equations vs %d unknowns (non-square)', numEq, NX);
end

% ---- helpers -------------------------------------------------------------
function emitHelpers(w, helpers)
for i = 1:numel(helpers)
    w.add('double %s = %s;', helpers(i).lhs, helpers(i).rhs);
end
end

% One residual slot: value + (guarded) Jacobian rows. `tag` makes generated
% helper names unique across the (possibly branched) equation set.
function numEq = addPlainEq(w, eq, ir, rows, perShockBases, n, numEq, tag)
import gdsge.codegen.cxx.sympymodel.diffRHS
import gdsge.codegen.cxx.sympymodel.combine
import gdsge.codegen.cxx.sympymodel.loadLoopLocals
eqExpr = gdsge.codegen.cxx.sympymodel.liftIndexed(eq.expr, rows, w, perShockBases);
if eq.primed
    d = diffRHS(eqExpr, rows.keys, sprintf('hpe%s_', tag));
    erow = combine(d.partialsMap, rows);
    w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++) {', n);
    loadLoopLocals(w, d.freeSymbols, ir.shocks.names, perShockBases);
    emitHelpers(w, d.helpers);
    w.add('GDSGE_f[%d+GDSGE_iter]=%s;', numEq - 1, d.value);
    w.add('if (GDSGE_jac) {');
    for k = 1:numel(erow)
        if erow(k).templated
            w.add('JAC(%d+GDSGE_iter,%d+(GDSGE_iter-1))=%s;', ...
                numEq - 1, erow(k).slot, erow(k).expr);
        else
            w.add('JAC(%d+GDSGE_iter,%d)=%s;', numEq - 1, erow(k).slot, erow(k).expr);
        end
    end
    w.add('}');
    w.add('}');
    numEq = numEq + n;
else
    d = diffRHS(eqExpr, rows.keys, sprintf('he%s_', tag));
    emitHelpers(w, d.helpers);
    w.add('GDSGE_f[%d]=%s;', numEq, d.value);
    erow = combine(d.partialsMap, rows);
    w.add('if (GDSGE_jac) {');
    for k = 1:numel(erow)
        assert(~erow(k).templated, ...
            'scalar equation %s has a per-shock gradient slot (modeling error)', tag);
        w.add('JAC(%d,%d)=%s;', numEq, erow(k).slot, erow(k).expr);
    end
    w.add('}');
    numEq = numEq + 1;
end
end

function txt = argUnpackDouble(ir)
% double unknowns bound by slot (SymPy path has no adouble).
w = gdsge.codegen.codeWriter();
for i = 1:numel(ir.variables.policy)
    p = ir.variables.policy{i};
    if p.length == 1
        w.add('double& %s = GDSGE_x[%d];', p.name, p.slot(1) - 1);
    else
        w.add('double* %s_GDSGE_GRID = GDSGE_x+(%d);', p.name, p.slot(1) - 1);
        w.add('#define %s_GRID(idx) %s_GDSGE_GRID[int(idx)-1]', p.name, p.name);
        w.add('#define %s(idx) %s_GDSGE_GRID[int(idx)-1]', p.name, p.name);
    end
end
txt = w.str();
end
