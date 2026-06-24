function rows = addPrimedAssign(w, st, ir, rows, perShockBases, prefix)
% ADDPRIMEDASSIGN  Emit a per-shock loop computing a primed assign's value array
%   <tgt>_GDSGE_GRID[j] and, per gradient slot it touches, a per-shock gradient
%   array (<tgt>_G_<slot> fixed, <tgt>_Gt_<base> templated). Registers a
%   per-shock row whose exprs read those arrays at [GDSGE_iter-1] — consumed by
%   the primed equation that references this assign inside its own shock loop.
import gdsge.codegen.cxx.sympymodel.combine
n = ir.shocks.count;
expr = gdsge.codegen.cxx.sympymodel.liftIndexed(st.expr, rows, w, perShockBases);
d = gdsge.codegen.cxx.sympymodel.diffRHS(expr, rows.keys, prefix);
brow = combine(d.partialsMap, rows);
tgt = st.target;

% value array <tgt>_GDSGE_GRID is pre-declared in generate; declare grad arrays.
for k = 1:numel(brow)
    w.add('double %s[%d];', gradArr(tgt, brow(k)), n);
end
w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++) {', n);
gdsge.codegen.cxx.sympymodel.loadLoopLocals(w, d.freeSymbols, ir.shocks.names, perShockBases);
for i = 1:numel(d.helpers)
    w.add('double %s = %s;', d.helpers(i).lhs, d.helpers(i).rhs);
end
w.add('%s_GDSGE_GRID[GDSGE_iter-1] = %s;', tgt, d.value);
for k = 1:numel(brow)
    w.add('%s[GDSGE_iter-1] = %s;', gradArr(tgt, brow(k)), brow(k).expr);
end
w.add('}');

prow = struct('slot', {}, 'expr', {}, 'templated', {});
for k = 1:numel(brow)
    prow(end+1) = struct('slot', brow(k).slot, ...
        'expr', sprintf('%s[GDSGE_iter-1]', gradArr(tgt, brow(k))), ...
        'templated', brow(k).templated); %#ok<AGROW>
end
% register under both the plain and primed names: a primed assign may be
% consumed as <name> (in a primed equation, HL1996 w1_consis) or as <name>'
% (in another per-shock body, Mendoza flow_future).
rows(tgt) = prow;
rows([tgt '__prime']) = prow;   %#ok<NASGU> containers.Map handle: mutated in place
end

function nm = gradArr(tgt, entry)
if isfield(entry, 'templated') && entry.templated
    nm = sprintf('%s_Gt_%d', tgt, entry.slot);
else
    nm = sprintf('%s_G_%d', tgt, entry.slot);
end
end
