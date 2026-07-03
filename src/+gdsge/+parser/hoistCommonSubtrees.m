function statements = hoistCommonSubtrees(statements)
% HOISTCOMMONSUBTREES  Region-wide common-subexpression hoisting.
%   Any expression subtree of >= THRESH nodes is computed once into a
%   generated GDSGE_CSE_<n> local (or reuses an earlier statement whose whole
%   RHS is that subtree) and referenced by name afterwards. Values are
%   unchanged: subtrees are pure, hoisting preserves the operation order
%   inside each subtree, and a memo entry dies as soon as any variable it
%   (transitively) depends on is reassigned.
%
%   Why: MSVC gives every adept expression temporary its own stack slot with
%   no cross-statement reuse, so a CCP-style model that repeats a 13-deep max
%   chain 15x per line produced a 2.3 MB eval frame against vcomp's fixed
%   1 MB worker stacks (matlab.exe's PE stack reserve). Hoisting shrinks the
%   frame and the autodiff tape together. Two tiers: any subtree of >=
%   THRESH_BIG nodes hoists unconditionally; a subtree of >= THRESH_SMALL
%   nodes hoists only when its canonical form occurs more than once in the
%   region (counted up front), so hand-sized one-off expressions stay inline.
THRESH_BIG = 24;
THRESH_SMALL = 8;
counts = tally(statements, THRESH_SMALL);
memo = struct('key',{},'refName',{},'primed',{},'names',{});
out = {};
counter = 0;
for i = 1:numel(statements)
    s = statements{i};
    pre = {};
    rootMemo = [];
    switch s.type
        case 'assign'
            [s.expr, pre, memo, counter] = rewrite(s.expr, true, memo, counter, ...
                THRESH_BIG, THRESH_SMALL, counts);
            % whole-RHS memo: a later identical subtree (or whole RHS) reuses
            % this target instead of recomputing. Self-referencing RHS
            % (x = f(..x..)) is never memoized: the cached value would be
            % computed from the pre-assignment x.
            [key, cnt, names, prm] = canon(s.expr); %#ok<ASGLU>
            if cnt >= THRESH_SMALL
                idx = lookup(memo, key);
                if ~isempty(idx)
                    s.expr = refNode(memo(idx));
                elseif ~ismember(s.target, names)
                    rootMemo = struct('key',key, 'refName',s.target, ...
                        'primed',s.primed ~= 0, ...
                        'names',{expandDeps(memo, names)});
                end
            end
        case 'reduction'
            [s.body, pre, memo, counter] = rewrite(s.body, true, memo, counter, ...
                THRESH_BIG, THRESH_SMALL, counts);
    end
    out = [out, pre, {s}]; %#ok<AGROW>
    memo = invalidate(memo, targetsOf(s));
    if ~isempty(rootMemo)
        memo(end+1) = rootMemo; %#ok<AGROW>
    end
end
statements = out;
end

% ---------------------------------------------------------------------------
function [node, pre, memo, counter] = rewrite(node, isRoot, memo, counter, BIG, SMALL, counts)
% Post-order: rewrite children, then hoist/reuse this subtree if it is big,
% or medium-sized and repeated somewhere in the region (counts, keyed on the
% pre-rewrite canonical — identical originals rewrite identically). The
% statement root itself is never hoisted (its statement already names it).
pre = {};
[key0, cnt0] = canon(node);           % pre-rewrite form, for the repeat tally
if cnt0 < SMALL                       % too small to ever hoist: skip descent
    return                            % (children are smaller still)
end
switch node.kind
    case 'unop'
        [node.arg, pre, memo, counter] = rewrite(node.arg, false, memo, counter, BIG, SMALL, counts);
    case 'binop'
        [node.lhs, p1, memo, counter] = rewrite(node.lhs, false, memo, counter, BIG, SMALL, counts);
        [node.rhs, p2, memo, counter] = rewrite(node.rhs, false, memo, counter, BIG, SMALL, counts);
        pre = [p1, p2];
    case 'call'
        for k = 1:numel(node.args)
            [node.args{k}, pk, memo, counter] = rewrite(node.args{k}, false, memo, counter, BIG, SMALL, counts);
            pre = [pre, pk]; %#ok<AGROW>
        end
    case 'index'
        [node.base, p1, memo, counter] = rewrite(node.base, false, memo, counter, BIG, SMALL, counts);
        pre = p1;
        for k = 1:numel(node.args)
            [node.args{k}, pk, memo, counter] = rewrite(node.args{k}, false, memo, counter, BIG, SMALL, counts);
            pre = [pre, pk]; %#ok<AGROW>
        end
end
if isRoot
    return
end
[key, cnt, names, prm] = canon(node);
if ~(cnt >= BIG || (cnt >= SMALL && counts.isKey(key0) && counts(key0) > 1))
    return
end
idx = lookup(memo, key);
if ~isempty(idx)
    node = refNode(memo(idx));
    return
end
counter = counter + 1;
genName = sprintf('GDSGE_CSE_%d', counter);
pre{end+1} = struct('type','assign', 'target',genName, 'primed',double(prm), 'expr',node);
memo(end+1) = struct('key',key, 'refName',genName, 'primed',prm, ...
    'names',{expandDeps(memo, names)});
node = refNode(memo(end));
end

% ---------------------------------------------------------------------------
function counts = tally(statements, SMALL)
% Occurrence count of every subtree (>= SMALL nodes) canonical form across
% the region's assign exprs and reduction bodies, on the original ASTs.
counts = containers.Map('KeyType','char','ValueType','double');
for i = 1:numel(statements)
    s = statements{i};
    switch s.type
        case 'assign';    tallyNode(s.expr);
        case 'reduction'; tallyNode(s.body);
    end
end
    function cnt = tallyNode(node)
        kids = gdsge.ir.node.children(node);
        cnt = 1;
        for k = 1:numel(kids)
            cnt = cnt + tallyNode(kids{k});
        end
        if cnt >= SMALL
            key = canon(node);
            if counts.isKey(key)
                counts(key) = counts(key) + 1;
            else
                counts(key) = 1;
            end
        end
    end
end

% ---------------------------------------------------------------------------
function [key, cnt, names, prm] = canon(node)
% One walk: canonical serialization, node count, referenced identifiers
% (variables and call targets, for reassignment invalidation), primed flag.
switch node.kind
    case 'num'
        key = sprintf('#%.17g', node.value); cnt = 1; names = {}; prm = false;
    case 'name'
        key = ['n:' node.id]; cnt = 1; names = {node.id}; prm = false;
    case 'primed'
        key = ['p:' node.id]; cnt = 1; names = {node.id}; prm = true;
    case 'unop'
        [k1, c1, names, prm] = canon(node.arg);
        key = ['u' node.op '(' k1 ')']; cnt = c1 + 1;
    case 'binop'
        [k1, c1, n1, p1] = canon(node.lhs);
        [k2, c2, n2, p2] = canon(node.rhs);
        key = ['b' node.op '(' k1 ',' k2 ')'];
        cnt = c1 + c2 + 1; names = [n1, n2]; prm = p1 || p2;
    case 'call'
        cnt = 1; names = {node.fn}; prm = false; parts = cell(1, numel(node.args));
        for k = 1:numel(node.args)
            [parts{k}, ck, nk, pk] = canon(node.args{k});
            cnt = cnt + ck; names = [names, nk]; prm = prm || pk; %#ok<AGROW>
        end
        key = ['c:' node.fn '(' strjoin(parts, ',') ')'];
    case 'index'
        [kb, cb, nb, pb] = canon(node.base);
        cnt = cb + 1; names = nb; prm = pb; parts = cell(1, numel(node.args));
        for k = 1:numel(node.args)
            [parts{k}, ck, nk, pk] = canon(node.args{k});
            cnt = cnt + ck; names = [names, nk]; prm = prm || pk; %#ok<AGROW>
        end
        key = ['i:' kb '(' strjoin(parts, ',') ')'];
    otherwise
        error('gdsge:parser:badExpr', 'unknown node kind "%s" in CSE pass', node.kind);
end
end

% ---------------------------------------------------------------------------
function idx = lookup(memo, key)
idx = find(strcmp({memo.key}, key), 1);
end

function node = refNode(entry)
if entry.primed
    node = gdsge.ir.node.primed(entry.refName);
else
    node = gdsge.ir.node.name(entry.refName);
end
end

function names = expandDeps(memo, names)
% Union in the dependency sets of any referenced memo targets so that
% invalidation is transitive (chain -> CSE_j -> b3: reassigning b3 must kill
% entries built on CSE_j too).
for j = 1:numel(memo)
    if ismember(memo(j).refName, names)
        names = [names, memo(j).names]; %#ok<AGROW>
    end
end
names = unique(names);
end

function t = targetsOf(s)
switch s.type
    case {'assign','reduction'}
        t = {s.target};
    case 'interpCall'
        t = s.targets;
    otherwise
        t = {};
end
end

function memo = invalidate(memo, targets)
% Kill entries that depend on any reassigned variable — including entries
% whose cached value IS that variable (refName).
if isempty(memo) || isempty(targets); return; end
keep = arrayfun(@(e) ~any(ismember(targets, e.names)) ...
    && ~any(strcmp(e.refName, targets)), memo);
memo = memo(keep);
end
