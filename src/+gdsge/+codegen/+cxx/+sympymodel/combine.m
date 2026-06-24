function row = combine(partials, freeSymRows)
% COMBINE  Chain rule: row[m] = sum_s partials(s) * freeSymRows(s)[m].
%   partials: containers.Map symbolName -> C++ partial expr (nonzero only).
%   freeSymRows: containers.Map symbolName -> gradRow (struct array of mkEntry).
%   Returns a consolidated gradRow. Entries are keyed by (templated, baseSlot):
%   fixed and templated slots never merge (structurally disjoint — scalar
%   unknowns vs array-unknown elements), and two contributions to the same slot
%   are summed. Symbols absent from freeSymRows contribute nothing (constants).
acc = containers.Map('KeyType', 'char', 'ValueType', 'any');
keys = partials.keys;
for i = 1:numel(keys)
    s = keys{i};
    if ~isKey(freeSymRows, s); continue; end
    p = partials(s);
    srow = freeSymRows(s);
    for j = 1:numel(srow)
        e = srow(j);
        t = isfield(e, 'templated') && e.templated;
        key = sprintf('%d|%d', t, e.slot);
        term = mulExpr(p, e.expr);
        if isKey(acc, key)
            cur = acc(key);
            cur.expr = [cur.expr '+' term];
            acc(key) = cur;
        else
            acc(key) = struct('slot', e.slot, 'expr', term, 'templated', t);
        end
    end
end
ks = acc.keys;
row = struct('slot', {}, 'expr', {}, 'templated', {});
for i = 1:numel(ks)
    row(end+1) = acc(ks{i}); %#ok<AGROW>
end
end

function s = mulExpr(a, b)
% Drop trivial *1 factors to keep the emitted C++ readable.
if strcmp(a, '1'); s = b; return; end
if strcmp(b, '1'); s = a; return; end
s = ['(' a ')*(' b ')'];
end
