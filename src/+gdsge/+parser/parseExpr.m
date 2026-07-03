function [ast, pos] = parseExpr(toks, pos)
% PARSEEXPR  Precedence-climbing expression parser over TOKENIZE's stream.
%   [AST, POS] = PARSEEXPR(TOKS[, POS]) parses the longest expression starting
%   at POS (default 1) and returns the AST plus the position of the first
%   unconsumed token (callers decide what may legally follow: eof, ',', ')').
%   MATLAB semantics: every binary op is left-associative, including ^;
%   unary +/-/~ bind between power and multiplicative; power accepts a unary
%   right operand (2^-3). Elementwise ops normalize to scalar ('.*' -> '*').
%   Row-vector literals [a b c] follow MATLAB's matrix-row rule (a depth-0
%   comma or whitespace gap separates elements; "[a -b]" is two elements,
%   "[a - b]" is one) and exist only transiently: ops and single-argument
%   functions broadcast elementwise, sum/prod/max/min collapse to scalar
%   (nested 2-arg max/min calls — the old toolbox's str2sym+ccode expansion),
%   and a vector that survives to the top of the expression is an error, so
%   the IR stays scalar.
if nargin < 2; pos = 1; end
[ast, pos] = parseBinary(toks, pos, 1);
if strcmp(ast.kind, 'vec')
    error('gdsge:parser:badExpr', ...
        'a vector literal must be reduced to a scalar by sum/prod/max/min (line %d, col %d)', ...
        ast.line, ast.col);
end
end

function [ast, pos] = parseBinary(toks, pos, level)
levels = { {'||'}, {'&&'}, {'|'}, {'&'}, ...
           {'==','~=','<','>','<=','>='}, ...
           {'+','-'}, ...
           {'*','/','.*','./'} };
if level > numel(levels)
    [ast, pos] = parseUnary(toks, pos);
    return
end
[ast, pos] = parseBinary(toks, pos, level + 1);
while strcmp(toks(pos).type, 'op') && any(strcmp(toks(pos).text, levels{level}))
    op = normOp(toks(pos).text);
    opTok = toks(pos);
    [rhs, pos] = parseBinary(toks, pos + 1, level + 1);
    ast = mkBinop(op, ast, rhs, opTok);
end
end

function [ast, pos] = parseUnary(toks, pos)
t = toks(pos);
if strcmp(t.type, 'op') && any(strcmp(t.text, {'+','-','~'}))
    [arg, pos] = parseUnary(toks, pos + 1);
    ast = mkUnop(t.text, arg);
else
    [ast, pos] = parsePower(toks, pos);
end
end

function [ast, pos] = parsePower(toks, pos)
[ast, pos] = parsePostfix(toks, pos);
while strcmp(toks(pos).type, 'op') && any(strcmp(toks(pos).text, {'^','.^'}))
    opTok = toks(pos);
    [rhs, pos] = parsePowerOperand(toks, pos + 1);
    ast = mkBinop('^', ast, rhs, opTok);
end
end

function [ast, pos] = parsePowerOperand(toks, pos)
% MATLAB's level-2 "power with unary minus": 2^-3 is legal.
t = toks(pos);
if strcmp(t.type, 'op') && any(strcmp(t.text, {'+','-','~'}))
    [arg, pos] = parsePowerOperand(toks, pos + 1);
    ast = mkUnop(t.text, arg);
else
    [ast, pos] = parsePostfix(toks, pos);
end
end

function [ast, pos] = parsePostfix(toks, pos)
[ast, pos] = parsePrimary(toks, pos);
while true
    t = toks(pos);
    if strcmp(t.type, 'prime')
        if ~strcmp(ast.kind, 'name')
            error('gdsge:parser:badExpr', ...
                'the future marker '' can only follow a variable name (line %d, col %d)', ...
                t.line, t.col);
        end
        ast = gdsge.ir.node.primed(ast.id);
        pos = pos + 1;
    elseif strcmp(t.type, 'punct') && strcmp(t.text, '(') && strcmp(ast.kind, 'name')
        fn = ast.id;
        [args, pos] = parseArgs(toks, pos + 1);
        ast = mkCall(fn, args, t);
    elseif strcmp(t.type, 'punct') && strcmp(t.text, '(') && strcmp(ast.kind, 'primed')
        error('gdsge:parser:badExpr', ...
            'a named-interp call (%s''(...)) must be the entire right-hand side of its own statement (line %d, col %d)', ...
            ast.id, t.line, t.col);
    else
        break
    end
end
end

function [ast, pos] = parsePrimary(toks, pos)
t = toks(pos);
if strcmp(t.type, 'number')
    ast = gdsge.ir.node.num(str2double(t.text));
    pos = pos + 1;
elseif strcmp(t.type, 'name')
    ast = gdsge.ir.node.name(t.text);
    pos = pos + 1;
elseif strcmp(t.type, 'punct') && strcmp(t.text, '(')
    [ast, pos] = parseBinary(toks, pos + 1, 1);
    if ~(strcmp(toks(pos).type, 'punct') && strcmp(toks(pos).text, ')'))
        error('gdsge:parser:badExpr', 'expected ")" at line %d, col %d', ...
            toks(pos).line, toks(pos).col);
    end
    pos = pos + 1;
elseif strcmp(t.type, 'punct') && strcmp(t.text, '[')
    [ast, pos] = parseVector(toks, pos);
else
    error('gdsge:parser:badExpr', 'unexpected token "%s" at line %d, col %d', ...
        t.text, t.line, t.col);
end
end

function [args, pos] = parseArgs(toks, pos)
args = {};
if strcmp(toks(pos).type, 'punct') && strcmp(toks(pos).text, ')')
    pos = pos + 1; return
end
while true
    [a, pos] = parseBinary(toks, pos, 1);
    args{end+1} = a; %#ok<AGROW>
    t = toks(pos);
    if strcmp(t.type, 'punct') && strcmp(t.text, ',')
        pos = pos + 1;
    elseif strcmp(t.type, 'punct') && strcmp(t.text, ')')
        pos = pos + 1; return
    else
        error('gdsge:parser:badExpr', 'expected "," or ")" at line %d, col %d', ...
            t.line, t.col);
    end
end
end

function op = normOp(op)
switch op
    case '.*'; op = '*';
    case './'; op = '/';
    case '.^'; op = '^';
end
end

% ===== row-vector literals ==================================================
% A 'vec' node (kind/elems/line/col) never reaches the IR: it is built here,
% broadcast by mkBinop/mkUnop/mkCall, and must collapse before parseExpr
% returns.

function [ast, pos] = parseVector(toks, openPos)
% openPos points at '['. Split the bracketed span into element token slices
% (depth-0 commas and MATLAB's whitespace rule), parse each as a full
% expression, and splice nested vectors (concatenation).
closePos = matchBracket(toks, openPos);
inner = openPos + 1 : closePos - 1;
if isempty(inner)
    error('gdsge:parser:badExpr', 'empty vector literal at line %d, col %d', ...
        toks(openPos).line, toks(openPos).col);
end
% element start indices within the token stream
starts = inner(1); depth = 0;
isComma = false(1, 0);
for i = inner
    isPunct = strcmp(toks(i).type, 'punct');
    if depth == 0
        if isPunct && strcmp(toks(i).text, ',')
            if i == inner(end)
                error('gdsge:parser:badExpr', ...
                    'trailing "," in vector literal at line %d, col %d', ...
                    toks(i).line, toks(i).col);
            end
            starts(end+1) = i + 1; %#ok<AGROW>
            isComma(numel(starts)-1) = true;
        elseif i > inner(1) && ~any(starts == i) ...
                && gapBefore(toks, i) && endsOperand(toks(i-1)) && startsElement(toks, i)
            starts(end+1) = i; %#ok<AGROW>
        end
    end
    if isPunct
        if any(strcmp(toks(i).text, {'(','[','{'}))
            depth = depth + 1;
        elseif any(strcmp(toks(i).text, {')',']','}'}))
            depth = depth - 1;
        end
    end
end
elems = {};
for k = 1:numel(starts)
    if k < numel(starts)
        stop = starts(k+1) - 1;
        if numel(isComma) >= k && isComma(k); stop = stop - 1; end   % drop the ','
    else
        stop = closePos - 1;
    end
    slice = toks(starts(k):stop);
    if isempty(slice)
        error('gdsge:parser:badExpr', 'empty vector element at line %d, col %d', ...
            toks(starts(k)).line, toks(starts(k)).col);
    end
    slice(end+1) = struct('type','eof','text','','line',slice(end).line, ...
        'col',slice(end).col + numel(slice(end).text)); %#ok<AGROW>
    [e, q] = parseBinary(slice, 1, 1);
    if ~strcmp(slice(q).type, 'eof')
        error('gdsge:parser:badExpr', ...
            'unexpected "%s" in vector literal at line %d, col %d', ...
            slice(q).text, slice(q).line, slice(q).col);
    end
    if strcmp(e.kind, 'vec')
        elems = [elems, e.elems]; %#ok<AGROW>  concatenation: [[a b] c]
    else
        elems{end+1} = e; %#ok<AGROW>
    end
end
ast = struct('kind','vec', 'elems',{elems}, ...
    'line',toks(openPos).line, 'col',toks(openPos).col);
pos = closePos + 1;
end

function closePos = matchBracket(toks, openPos)
depth = 0;
for i = openPos:numel(toks)
    if strcmp(toks(i).type, 'punct')
        if strcmp(toks(i).text, '['); depth = depth + 1;
        elseif strcmp(toks(i).text, ']')
            depth = depth - 1;
            if depth == 0; closePos = i; return; end
        end
    end
end
error('gdsge:parser:badExpr', 'unmatched "[" at line %d, col %d', ...
    toks(openPos).line, toks(openPos).col);
end

function tf = gapBefore(toks, i)
% whitespace (or a line break) between token i-1 and token i
tf = toks(i).line ~= toks(i-1).line ...
    || toks(i).col > toks(i-1).col + numel(toks(i-1).text);
end

function tf = endsOperand(t)
tf = any(strcmp(t.type, {'name','number','prime'})) ...
    || (strcmp(t.type, 'punct') && any(strcmp(t.text, {')',']'})));
end

function tf = startsElement(toks, i)
% MATLAB matrix-row rule: after a gap, a new element starts at an operand
% token, or at a sign glued (no gap) to the operand it negates ("[a -b]").
t = toks(i);
if any(strcmp(t.type, {'name','number'})) ...
        || (strcmp(t.type, 'punct') && any(strcmp(t.text, {'(','['})))
    tf = true;
elseif strcmp(t.type, 'op') && any(strcmp(t.text, {'+','-','~'})) ...
        && i < numel(toks) && ~gapBefore(toks, i+1)
    nt = toks(i+1);
    tf = any(strcmp(nt.type, {'name','number'})) ...
        || (strcmp(nt.type, 'punct') && any(strcmp(nt.text, {'(','['})));
else
    tf = false;
end
end

% ===== vector broadcasting ==================================================
function tf = isVec(n)
tf = strcmp(n.kind, 'vec');
end

function tf = isZeroNum(n)
tf = strcmp(n.kind, 'num') && n.value == 0;
end

function ast = mkBinop(op, l, r, opTok)
lv = isVec(l); rv = isVec(r);
if ~lv && ~rv
    % Fold 0*x -> 0 (the old toolbox's str2sym simplification did this;
    % models rely on it to reference undeclared names, e.g. "+0*m'").
    if strcmp(op, '*') && (isZeroNum(l) || isZeroNum(r))
        ast = gdsge.ir.node.num(0);
        return
    end
    ast = gdsge.ir.node.binop(op, l, r);
    return
end
if lv && rv
    nl = numel(l.elems); nr = numel(r.elems);
    if nl == nr
        elems = cellfun(@(a,b) gdsge.ir.node.binop(op, a, b), ...
            l.elems, r.elems, 'UniformOutput', false);
    elseif nl == 1
        elems = cellfun(@(b) gdsge.ir.node.binop(op, l.elems{1}, b), ...
            r.elems, 'UniformOutput', false);
    elseif nr == 1
        elems = cellfun(@(a) gdsge.ir.node.binop(op, a, r.elems{1}), ...
            l.elems, 'UniformOutput', false);
    else
        error('gdsge:parser:badExpr', ...
            'vector length mismatch (%d vs %d) at "%s" (line %d, col %d)', ...
            nl, nr, opTok.text, opTok.line, opTok.col);
    end
    ast = l; ast.elems = elems;
elseif lv
    ast = l;
    ast.elems = cellfun(@(a) gdsge.ir.node.binop(op, a, r), ...
        l.elems, 'UniformOutput', false);
else
    ast = r;
    ast.elems = cellfun(@(b) gdsge.ir.node.binop(op, l, b), ...
        r.elems, 'UniformOutput', false);
end
end

function ast = mkUnop(op, a)
if isVec(a)
    ast = a;
    ast.elems = cellfun(@(e) gdsge.ir.node.unop(op, e), ...
        a.elems, 'UniformOutput', false);
else
    ast = gdsge.ir.node.unop(op, a);
end
end

function ast = mkCall(fn, args, t)
vecMask = cellfun(@(a) isVec(a), args);
if ~any(vecMask)
    ast = gdsge.ir.node.call(fn, args);
    return
end
if numel(args) ~= 1
    error('gdsge:parser:badExpr', ...
        'a vector literal argument to %s must be the only argument (line %d, col %d)', ...
        fn, t.line, t.col);
end
v = args{1};
switch fn
    case {'max','min'}   % collapse to nested 2-arg calls (old ccode expansion)
        ast = v.elems{1};
        for k = 2:numel(v.elems)
            ast = gdsge.ir.node.call(fn, {ast, v.elems{k}});
        end
    case 'sum'
        ast = foldBinop('+', v.elems);
    case 'prod'
        ast = foldBinop('*', v.elems);
    otherwise            % elementwise map: exp([...]), log([...]), ...
        ast = v;
        ast.elems = cellfun(@(e) gdsge.ir.node.call(fn, {e}), ...
            v.elems, 'UniformOutput', false);
end
end

function ast = foldBinop(op, elems)
ast = elems{1};
for k = 2:numel(elems)
    ast = gdsge.ir.node.binop(op, ast, elems{k});
end
end
