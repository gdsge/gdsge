function [ast, pos] = parseExpr(toks, pos)
% PARSEEXPR  Precedence-climbing expression parser over TOKENIZE's stream.
%   [AST, POS] = PARSEEXPR(TOKS[, POS]) parses the longest expression starting
%   at POS (default 1) and returns the AST plus the position of the first
%   unconsumed token (callers decide what may legally follow: eof, ',', ')').
%   MATLAB semantics: every binary op is left-associative, including ^;
%   unary +/-/~ bind between power and multiplicative; power accepts a unary
%   right operand (2^-3). Elementwise ops normalize to scalar ('.*' -> '*').
if nargin < 2; pos = 1; end
[ast, pos] = parseBinary(toks, pos, 1);
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
    [rhs, pos] = parseBinary(toks, pos + 1, level + 1);
    ast = gdsge.ir.node.binop(op, ast, rhs);
end
end

function [ast, pos] = parseUnary(toks, pos)
t = toks(pos);
if strcmp(t.type, 'op') && any(strcmp(t.text, {'+','-','~'}))
    [arg, pos] = parseUnary(toks, pos + 1);
    ast = gdsge.ir.node.unop(t.text, arg);
else
    [ast, pos] = parsePower(toks, pos);
end
end

function [ast, pos] = parsePower(toks, pos)
[ast, pos] = parsePostfix(toks, pos);
while strcmp(toks(pos).type, 'op') && any(strcmp(toks(pos).text, {'^','.^'}))
    [rhs, pos] = parsePowerOperand(toks, pos + 1);
    ast = gdsge.ir.node.binop('^', ast, rhs);
end
end

function [ast, pos] = parsePowerOperand(toks, pos)
% MATLAB's level-2 "power with unary minus": 2^-3 is legal.
t = toks(pos);
if strcmp(t.type, 'op') && any(strcmp(t.text, {'+','-','~'}))
    [arg, pos] = parsePowerOperand(toks, pos + 1);
    ast = gdsge.ir.node.unop(t.text, arg);
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
        ast = gdsge.ir.node.call(fn, args);
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
