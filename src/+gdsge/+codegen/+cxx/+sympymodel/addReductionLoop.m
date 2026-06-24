function rows = addReductionLoop(w, st, ir, rows, perShockBases, NX, prefix)
% ADDREDUCTIONLOOP  Emit a fused shock loop for an EXPECT/MIN/MAX/PROD reduction.
%   Accumulates the value (scalar st.target) and its gradient into d<target>[NX].
%   The body's gradient (combine of the registry) may have templated (per-shock)
%   entries from interp results / future-state unknowns; the loop writes those
%   into the per-shock diagonal slot base+(GDSGE_iter-1). The reduction COLLAPSES
%   per-shock gradients into a scalar row (fixed integer slots) which it
%   registers for st.target.
import gdsge.codegen.cxx.sympymodel.combine
n = ir.shocks.count;
body = gdsge.codegen.cxx.sympymodel.liftIndexed(st.body, rows, w, perShockBases);
d = gdsge.codegen.cxx.sympymodel.diffRHS(body, rows.keys, prefix);
brow = combine(d.partialsMap, rows);
tgt = st.target;

switch st.kind
    case 'EXPECT', w.add('double %s = 0.0;', tgt);
    case 'PROD',   w.add('double %s = 1.0;', tgt);
    case 'MIN',    w.add('double %s = 1e20;', tgt);
    case 'MAX',    w.add('double %s = -1e20;', tgt);
    otherwise
        error('gdsge:codegen:sympyReductionKind', 'reduction kind %s', st.kind);
end
hasTempl = any([brow.templated]);
if any(strcmp(st.kind, {'PROD'})) && hasTempl
    error('gdsge:codegen:sympyProdTemplated', ...
        'PROD over per-shock (interp/future-state) bodies is not supported');
end
w.add('double d%s[%d];', tgt, NX);
w.add('for(int GDSGE_jz=0; GDSGE_jz<%d; ++GDSGE_jz) d%s[GDSGE_jz]=0.0;', NX, tgt);

w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++) {', n);
gdsge.codegen.cxx.sympymodel.loadLoopLocals(w, d.freeSymbols, ir.shocks.names, perShockBases);
for i = 1:numel(d.helpers)
    w.add('double %s = %s;', d.helpers(i).lhs, d.helpers(i).rhs);
end
w.add('double GDSGE_body = %s;', d.value);
trans = sprintf('%s((shock)+%d*(GDSGE_iter-1))', st.transRef, n);
switch st.kind
    case 'EXPECT'
        w.add('%s += %s*GDSGE_body;', tgt, trans);
        for k = 1:numel(brow)
            w.add('d%s[%s] += %s*(%s);', tgt, slotIdx(brow(k)), trans, brow(k).expr);
        end
    case 'PROD'
        for k = 1:numel(brow)
            w.add('d%s[%d] = d%s[%d]*GDSGE_body + %s*(%s);', ...
                tgt, brow(k).slot, tgt, brow(k).slot, tgt, brow(k).expr);
        end
        w.add('%s *= GDSGE_body;', tgt);
    case {'MIN','MAX'}
        cmp = '<'; if strcmp(st.kind, 'MAX'); cmp = '>'; end
        w.add('if (GDSGE_body %s %s) {', cmp, tgt);
        w.add('%s = GDSGE_body;', tgt);
        % gradient of the extremum = winning shock's body gradient; clear stale.
        w.add('for(int GDSGE_jz=0; GDSGE_jz<%d; ++GDSGE_jz) d%s[GDSGE_jz]=0.0;', NX, tgt);
        for k = 1:numel(brow)
            w.add('d%s[%s] = %s;', tgt, slotIdx(brow(k)), brow(k).expr);
        end
        w.add('}');
end
w.add('}');

% Register the COLLAPSED scalar row over the UNIQUE set of touched slots: a
% templated entry enumerates to its n diagonal slots base..base+n-1, a fixed
% entry contributes its slot. d<tgt>[slot] already accumulates ALL (templated +
% fixed) contributions to that slot, so each slot is registered ONCE — never
% twice (the MIN-shift pattern can otherwise alias a templated and a fixed entry
% onto the same slot, which a duplicate row would double via combine).
touchedSlots = [];
for k = 1:numel(brow)
    if brow(k).templated
        touchedSlots = [touchedSlots, brow(k).slot + (0:n-1)]; %#ok<AGROW>
    else
        touchedSlots(end+1) = brow(k).slot; %#ok<AGROW>
    end
end
touchedSlots = unique(touchedSlots);
trow = struct('slot', {}, 'expr', {}, 'templated', {});
for k = 1:numel(touchedSlots)
    s = touchedSlots(k);
    trow(end+1) = struct('slot', s, ...
        'expr', sprintf('d%s[%d]', tgt, s), 'templated', false); %#ok<AGROW>
end
rows(tgt) = trow;   %#ok<NASGU> containers.Map handle: mutated in place
end

function s = slotIdx(entry)
% C++ index into the dacc array for an entry: fixed slot, or per-shock diagonal.
if isfield(entry, 'templated') && entry.templated
    s = sprintf('%d+(GDSGE_iter-1)', entry.slot);
else
    s = sprintf('%d', entry.slot);
end
end
