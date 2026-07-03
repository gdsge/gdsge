function model = parseModel(modelText)
% PARSEMODEL  model-block text -> ir.model struct:
%   .statements  ordered cell of assign / interpCall / reduction statements
%   interpCall forms: [a',b',...] = GDSGE_INTERP_VEC'(args) (one target per
%   declared var_interp) and the named single-target form xn' = name'(args)
%   (interpRef = the named var_interp).
%   .equations   cell of struct('expr',ast,'primed',0|1)
%   Reductions: GDSGE_<EXPECT|MIN|MAX|PROD>{ body [| transName] }, default
%   transition 'shock_trans'. A statement whose whole RHS is one reduction
%   becomes that reduction (target = LHS); an inline reduction is hoisted to a
%   generated GDSGE_<KIND>_<n> statement inserted just before its host. At most
%   one reduction per statement; reductions cannot nest. Statement order is
%   source order. Line numbers in errors are relative to the model block.

% --- 1. split off the nested equations;...end; region (line-wise) ----------
lines = regexp(modelText, '\n', 'split');
trimmed = cellfun(@strtrim, lines, 'UniformOutput', false);
eqOpen = find(strcmp(trimmed, 'equations;'));
if numel(eqOpen) ~= 1
    error('gdsge:parser:missingEquations', ...
        'model block must contain exactly one "equations;...end;" block (found %d)', ...
        numel(eqOpen));
end
eqClose = find(strcmp(trimmed, 'end;'));
eqClose = eqClose(eqClose > eqOpen);
if isempty(eqClose)
    error('gdsge:parser:missingEquations', '"equations;" block not closed by "end;"');
end
eqClose = eqClose(1);
if any(~cellfun(@isempty, trimmed(eqClose+1:end)))
    error('gdsge:parser:badStatement', ...
        'unexpected content after the equations block in the model body');
end
bodyText = strjoin(lines(1:eqOpen-1), sprintf('\n'));
eqText   = strjoin(lines(eqOpen+1:eqClose-1), sprintf('\n'));

% --- 2. model statements ----------------------------------------------------
[bodyStmts, bodyLines] = gdsge.parser.splitStatements(bodyText);
statements = {};
counter = 0;
for i = 1:numel(bodyStmts)
    % bodyLines(i) seeds tokenize's lineOffset so tokens carry block-relative lines
    toks = gdsge.parser.tokenize(bodyStmts{i}, bodyLines(i));
    [stmtList, counter] = parseStatement(toks, counter);
    statements = [statements, stmtList]; %#ok<AGROW>
end
statements = gdsge.parser.hoistCommonSubtrees(statements);

% --- 3. equations (line-aware: if/elseif/else/end group into conditionals) ---
equations = parseEquationRegion(eqText, eqOpen);

model = struct('statements', {statements}, 'equations', {equations});
end

% ===== equation region (plain + if/else conditional) =======================
function entries = parseEquationRegion(eqText, eqLineOffset)
% Walk the equations region line-wise. `if/elseif/else/end` (no trailing ';')
% are control lines that group plain `<expr>;` equations into a 'conditional'
% entry (CaoNie2016: if A==AGood; Xp(1); else; consis1; end). Plain equations
% are accumulated and split by ';' (so a multi-line equation still parses).
% One level of nesting only; a nested `if` errors.
lines = regexp(eqText, '\n', 'split');
n = numel(lines);
entries = {};
i = 1;
while i <= n
    t = strtrim(lines{i});
    if isempty(t); i = i + 1; continue; end
    if isIfLine(t)
        cases = { struct('cond', ifCondOf(t), 'lines', {{}}) };
        i = i + 1; closed = false;
        while i <= n
            tt = strtrim(lines{i});
            if isIfLine(tt)
                error('gdsge:parser:nestedConditionalEquation', ...
                    'nested "if" inside a conditional equation is not supported (line %d)', ...
                    eqLineOffset + i);
            elseif ~isempty(regexp(tt, '^elseif\W', 'once'))
                cases{end+1} = struct('cond', strtrim(regexprep(tt, '^elseif', '')), 'lines', {{}}); %#ok<AGROW>
                i = i + 1;
            elseif strcmp(tt, 'else')
                cases{end+1} = struct('cond', '', 'lines', {{}}); %#ok<AGROW>
                i = i + 1;
            elseif strcmp(tt, 'end') || strcmp(tt, 'end;')
                i = i + 1; closed = true; break;
            else
                cases{end}.lines{end+1} = lines{i}; %#ok<AGROW>
                i = i + 1;
            end
        end
        if ~closed
            error('gdsge:parser:unterminatedConditional', ...
                'conditional equation "if" not closed by "end" (line %d)', eqLineOffset + i);
        end
        caseEntries = cell(1, numel(cases));
        for c = 1:numel(cases)
            caseEntries{c} = struct('cond', cases{c}.cond, ...
                'equations', {parsePlainEquations(strjoin(cases{c}.lines, sprintf('\n')), eqLineOffset)});
        end
        entries{end+1} = struct('kind','conditional','cases', {caseEntries}); %#ok<AGROW>
    else
        buf = {lines{i}}; i = i + 1;
        while i <= n
            tt = strtrim(lines{i});
            if isIfLine(tt) || ~isempty(regexp(tt, '^elseif\W', 'once')) ...
                    || strcmp(tt, 'else') || strcmp(tt, 'end') || strcmp(tt, 'end;')
                break;
            end
            buf{end+1} = lines{i}; i = i + 1; %#ok<AGROW>
        end
        plain = parsePlainEquations(strjoin(buf, sprintf('\n')), eqLineOffset);
        entries = [entries, plain]; %#ok<AGROW>
    end
end
end

function tf = isIfLine(t)
tf = ~isempty(regexp(t, '^if\W', 'once')) || strcmp(t, 'if');
end

function c = ifCondOf(t)
c = strtrim(regexprep(t, '^if', ''));
end

function eqs = parsePlainEquations(text, lineOffset)
[stmts, slines] = gdsge.parser.splitStatements(text);
eqs = cell(1, numel(stmts));
for q = 1:numel(stmts)
    toks = gdsge.parser.tokenize(stmts{q}, lineOffset + slines(q));
    e0 = parseEquation(toks);
    eqs{q} = struct('kind','plain','expr', e0.expr, 'primed', e0.primed);
end
end

% ===== statement classification ============================================
function [stmtList, counter] = parseStatement(toks, counter)
if strcmp(toks(1).type, 'punct') && strcmp(toks(1).text, '[')
    stmtList = { parseInterpCall(toks) };
    return
end
if ~strcmp(toks(1).type, 'name')
    error('gdsge:parser:badStatement', ...
        'model statement must start with "[" or a variable name (line %d)', toks(1).line);
end
target = toks(1).text;
p = 2; primed = 0;
if strcmp(toks(p).type, 'prime'); primed = 1; p = p + 1; end
if ~(strcmp(toks(p).type, 'punct') && strcmp(toks(p).text, '='))
    error('gdsge:parser:badStatement', ...
        'expected "=" after "%s" in model statement (line %d)', target, toks(1).line);
end
rest = toks(p+1:end);   % keeps the eof sentinel

% Hoist embedded named-interp calls (name'(args) used as a sub-expression,
% not the whole RHS) to generated GDSGE_INTERP_<n>' interpCall statements, so
% the host references a primed name (Cao2011EZ: qp' = qFuture'(Wp')+d').
[icStmts, rest, counter] = hoistEmbeddedInterp(rest, counter);

% named-interp call: the entire RHS is interpName'(args)
if numel(rest) >= 4 && strcmp(rest(1).type, 'name') ...
        && strcmp(rest(2).type, 'prime') ...
        && strcmp(rest(3).type, 'punct') && strcmp(rest(3).text, '(')
    stmtList = [icStmts, { parseNamedInterpCall(target, primed, rest) }];
    return
end

[redIdx, kinds] = findReductions(rest);
if isempty(redIdx)
    [expr, q] = gdsge.parser.parseExpr(rest, 1);
    expectEof(rest, q);
    stmtList = [icStmts, { struct('type','assign','target',target,'primed',primed,'expr',expr) }];
    return
end

% whole-RHS single reduction: it becomes the statement itself (fast path —
% keeps existing models' IR snapshots byte-identical: no hoist, no counter).
if numel(redIdx) == 1
    closeIdx = matchBrace(rest, redIdx(1) + 1);
    assertNoNested(rest, redIdx(1), closeIdx);
    if redIdx(1) == 1 && strcmp(rest(closeIdx+1).type, 'eof')
        if primed
            error('gdsge:parser:badStatement', ...
                'a primed target cannot take a reduction directly (line %d)', toks(1).line);
        end
        [bodyAst, transRef] = parseReductionInner(rest(redIdx(1)+2:closeIdx-1), rest(redIdx(1)).line);
        stmtList = [icStmts, { struct('type','reduction','kind',kinds{1},'target',target, ...
                            'body',bodyAst,'transRef',transRef) }];
        return
    end
end

% general case: hoist every reduction (left-to-right) to its own generated
% GDSGE_<KIND>_<n> reduction statement, substituting its name into the host.
genStmts = {};
host = rest;
while true
    [hostRed, hostKinds] = findReductions(host);
    if isempty(hostRed); break; end
    ri = hostRed(1);
    ci = matchBrace(host, ri + 1);
    assertNoNested(host, ri, ci);
    [bodyAst, transRef] = parseReductionInner(host(ri+2:ci-1), host(ri).line);
    counter = counter + 1;
    genName = sprintf('GDSGE_%s_%d', hostKinds{1}, counter);
    genStmts{end+1} = struct('type','reduction','kind',hostKinds{1}, ...
        'target',genName,'body',bodyAst,'transRef',transRef); %#ok<AGROW>
    nameTok = host(ri); nameTok.type = 'name'; nameTok.text = genName;
    host = [host(1:ri-1), nameTok, host(ci+1:end)];
end
[expr, q] = gdsge.parser.parseExpr(host, 1);
expectEof(host, q);
stmtList = [icStmts, genStmts, { struct('type','assign','target',target, ...
                               'primed',primed,'expr',expr) }];
end

% Reject a reduction nested inside another reduction's braces.
function assertNoNested(toks, openRed, closeIdx)
inner = findReductions(toks(openRed+1:closeIdx-1));
if ~isempty(inner)
    error('gdsge:parser:nestedReduction', ...
        'reductions cannot nest (line %d); use an intermediate variable', ...
        toks(openRed).line);
end
end

% ===== embedded named-interp hoisting ======================================
function [icStmts, rest, counter] = hoistEmbeddedInterp(rest, counter)
% Repeatedly hoist a named-interp call name'(args) that appears embedded in a
% larger RHS (not the whole RHS, not inside reduction braces) into a generated
% GDSGE_INTERP_<n>' = name'(args) interpCall, replacing the span with the
% primed token GDSGE_INTERP_<n>'.
icStmts = {};
while true
    pos = findEmbeddedInterp(rest);
    if pos == 0; break; end
    openPar = pos + 2;
    closePar = matchParen(rest, openPar);
    interpName = rest(pos).text;
    argToks = [rest(openPar+1:closePar-1), eofTok(rest(closePar))];
    args = {}; q = 1;
    if ~strcmp(argToks(1).type, 'eof')
        while true
            [a, q] = gdsge.parser.parseExpr(argToks, q);
            args{end+1} = a; %#ok<AGROW>
            if strcmp(argToks(q).type,'punct') && strcmp(argToks(q).text,','); q = q + 1;
            elseif strcmp(argToks(q).type,'eof'); break;
            else
                error('gdsge:parser:badInterpCall', ...
                    'bad interp args (line %d, col %d)', argToks(q).line, argToks(q).col);
            end
        end
    end
    counter = counter + 1;
    genName = sprintf('GDSGE_INTERP_%d', counter);
    icStmts{end+1} = struct('type','interpCall','targets',{{genName}}, ...
        'primed',1,'args',{args},'interpRef',interpName); %#ok<AGROW>
    nameTok = rest(pos);  nameTok.type = 'name';  nameTok.text = genName;
    primeTok = rest(pos); primeTok.type = 'prime'; primeTok.text = '''';
    rest = [rest(1:pos-1), nameTok, primeTok, rest(closePar+1:end)];
end
end

function pos = findEmbeddedInterp(rest)
% Index of a name,prime,'(' triple at brace-depth 0 that is NOT the whole RHS
% and not the multi-target GDSGE_INTERP_VEC; 0 if none.
pos = 0; depth = 0;
for i = 1:numel(rest)-2
    if strcmp(rest(i).type, 'punct')
        if strcmp(rest(i).text, '{'); depth = depth + 1; continue;
        elseif strcmp(rest(i).text, '}'); depth = depth - 1; continue; end
    end
    if depth == 0 && strcmp(rest(i).type, 'name') ...
            && strcmp(rest(i+1).type, 'prime') ...
            && strcmp(rest(i+2).type, 'punct') && strcmp(rest(i+2).text, '(') ...
            && ~strcmp(rest(i).text, 'GDSGE_INTERP_VEC')
        cp = matchParen(rest, i+2);
        isWholeRhs = (i == 1) && strcmp(rest(cp+1).type, 'eof');
        if ~isWholeRhs; pos = i; return; end
    end
end
end

function cp = matchParen(rest, openIdx)
% openIdx points at '('; return the index of its matching ')'.
depth = 0;
for i = openIdx:numel(rest)
    if strcmp(rest(i).type, 'punct') && strcmp(rest(i).text, '('); depth = depth + 1;
    elseif strcmp(rest(i).type, 'punct') && strcmp(rest(i).text, ')')
        depth = depth - 1;
        if depth == 0; cp = i; return; end
    end
end
error('gdsge:parser:badInterpCall', 'unmatched "(" in interp call (line %d)', rest(openIdx).line);
end

% ===== interp call ==========================================================
function s = parseInterpCall(toks)
% [a',b',...] = GDSGE_INTERP_VEC'(arg1, arg2, ...)
p = 2;   % after '['
targets = {};
while true
    if ~strcmp(toks(p).type, 'name')
        error('gdsge:parser:badInterpCall', ...
            'expected a variable name in the interp target list (line %d, col %d)', ...
            toks(p).line, toks(p).col);
    end
    tname = toks(p).text; p = p + 1;
    if ~strcmp(toks(p).type, 'prime')
        error('gdsge:parser:badInterpCall', ...
            'interp targets must be primed, e.g. [psn'',...] (line %d, col %d)', ...
            toks(p).line, toks(p).col);
    end
    p = p + 1;
    targets{end+1} = tname; %#ok<AGROW>
    if strcmp(toks(p).type, 'punct') && strcmp(toks(p).text, ',')
        p = p + 1;
    elseif strcmp(toks(p).type, 'punct') && strcmp(toks(p).text, ']')
        p = p + 1; break
    else
        error('gdsge:parser:badInterpCall', 'expected "," or "]" (line %d, col %d)', ...
            toks(p).line, toks(p).col);
    end
end
if ~(strcmp(toks(p).type, 'punct') && strcmp(toks(p).text, '='))
    error('gdsge:parser:badInterpCall', 'expected "=" (line %d, col %d)', ...
        toks(p).line, toks(p).col);
end
p = p + 1;
if ~(strcmp(toks(p).type, 'name') && strcmp(toks(p).text, 'GDSGE_INTERP_VEC'))
    error('gdsge:parser:badInterpCall', ...
        'multi-target assignment must call GDSGE_INTERP_VEC (line %d, col %d)', ...
        toks(p).line, toks(p).col);
end
p = p + 1;
if ~strcmp(toks(p).type, 'prime')
    error('gdsge:parser:badInterpCall', ...
        'unprimed GDSGE_INTERP_VEC is not supported until Phase 7b (line %d, col %d)', ...
        toks(p).line, toks(p).col);
end
p = p + 1;
if ~(strcmp(toks(p).type, 'punct') && strcmp(toks(p).text, '('))
    error('gdsge:parser:badInterpCall', 'expected "(" (line %d, col %d)', ...
        toks(p).line, toks(p).col);
end
p = p + 1;
args = {};
while true
    [a, p] = gdsge.parser.parseExpr(toks, p);
    args{end+1} = a; %#ok<AGROW>
    if strcmp(toks(p).type, 'punct') && strcmp(toks(p).text, ',')
        p = p + 1;
    elseif strcmp(toks(p).type, 'punct') && strcmp(toks(p).text, ')')
        p = p + 1; break
    else
        error('gdsge:parser:badInterpCall', 'expected "," or ")" (line %d, col %d)', ...
            toks(p).line, toks(p).col);
    end
end
expectEof(toks, p);
s = struct('type','interpCall', 'targets',{targets}, 'primed',1, ...
           'args',{args}, 'interpRef','GDSGE_INTERP_VEC');
end

% ===== named-interp call ====================================================
function s = parseNamedInterpCall(target, primed, rest)
% target' = interpName'(arg1, ...): single-target named-interp evaluation.
% The call must be the entire RHS and the target must be primed (the result
% is per-shock, like an interpCall target).
if ~primed
    error('gdsge:parser:badInterpCall', ...
        'the target of a named-interp call must be primed, e.g. xn'' = %s''(...) (line %d)', ...
        rest(1).text, rest(1).line);
end
p = 4;   % after name, prime, '('
args = {};
while true
    [a, p] = gdsge.parser.parseExpr(rest, p);
    args{end+1} = a; %#ok<AGROW>
    if strcmp(rest(p).type, 'punct') && strcmp(rest(p).text, ',')
        p = p + 1;
    elseif strcmp(rest(p).type, 'punct') && strcmp(rest(p).text, ')')
        p = p + 1; break
    else
        error('gdsge:parser:badInterpCall', 'expected "," or ")" (line %d, col %d)', ...
            rest(p).line, rest(p).col);
    end
end
if ~strcmp(rest(p).type, 'eof')
    error('gdsge:parser:badInterpCall', ...
        'a named-interp call must be the entire right-hand side of its statement; use an intermediate variable (line %d, col %d)', ...
        rest(p).line, rest(p).col);
end
s = struct('type','interpCall', 'targets',{{target}}, 'primed',1, ...
           'args',{args}, 'interpRef',rest(1).text);
end

% ===== equations ============================================================
function e = parseEquation(toks)
% a bare "name'" is a future-indexed equation: the prime becomes a flag
if numel(toks) == 3 && strcmp(toks(1).type, 'name') && strcmp(toks(2).type, 'prime')
    e = struct('expr', gdsge.ir.node.name(toks(1).text), 'primed', 1);
    return
end
[ast, pos] = gdsge.parser.parseExpr(toks, 1);
expectEof(toks, pos);
e = struct('expr', ast, 'primed', 0);
end

% ===== reduction helpers ====================================================
function [idx, kinds] = findReductions(toks)
names = {'GDSGE_EXPECT','GDSGE_MIN','GDSGE_MAX','GDSGE_PROD'};
idx = []; kinds = {};
for i = 1:numel(toks)-1
    if strcmp(toks(i).type, 'name') && ismember(toks(i).text, names) ...
            && strcmp(toks(i+1).type, 'punct') && strcmp(toks(i+1).text, '{')
        idx(end+1) = i; %#ok<AGROW>
        kinds{end+1} = extractAfter(toks(i).text, 'GDSGE_'); %#ok<AGROW>
    end
end
end

function closeIdx = matchBrace(toks, openIdx)
% openIdx points at '{'; return the index of its matching '}'
depth = 0;
for i = openIdx:numel(toks)
    if strcmp(toks(i).type, 'punct')
        if strcmp(toks(i).text, '{')
            depth = depth + 1;
        elseif strcmp(toks(i).text, '}')
            depth = depth - 1;
            if depth == 0; closeIdx = i; return; end
        end
    end
end
error('gdsge:parser:unterminatedReduction', ...
    'unmatched "{" in reduction at line %d, col %d', ...
    toks(openIdx).line, toks(openIdx).col);
end

function [bodyAst, transRef] = parseReductionInner(inner, line)
% inner: tokens strictly between { and } (no eof). A depth-0 "|" splits
% body | transitionName; the default transition is shock_trans.
transRef = 'shock_trans';
depth = 0; pipeAt = 0;
for i = 1:numel(inner)
    if strcmp(inner(i).type, 'punct') && any(strcmp(inner(i).text, {'(','[','{'}))
        depth = depth + 1;
    elseif strcmp(inner(i).type, 'punct') && any(strcmp(inner(i).text, {')',']','}'}))
        depth = depth - 1;
    elseif strcmp(inner(i).type, 'op') && strcmp(inner(i).text, '|') && depth == 0
        pipeAt = i; break
    end
end
if pipeAt > 0
    transToks = inner(pipeAt+1:end);
    if numel(transToks) ~= 1 || ~strcmp(transToks(1).type, 'name')
        error('gdsge:parser:badStatement', ...
            'a reduction transition must be a single matrix name (line %d)', line);
    end
    transRef = transToks(1).text;
    inner = inner(1:pipeAt-1);
end
if isempty(inner)
    error('gdsge:parser:badStatement', 'empty reduction body (line %d)', line);
end
inner = [inner, eofTok(inner(end))];
[bodyAst, q] = gdsge.parser.parseExpr(inner, 1);
expectEof(inner, q);
end

% ===== small shared helpers =================================================
function expectEof(toks, pos)
if ~strcmp(toks(pos).type, 'eof')
    error('gdsge:parser:badExpr', 'unexpected "%s" at line %d, col %d', ...
        toks(pos).text, toks(pos).line, toks(pos).col);
end
end

function t = eofTok(prev)
t = struct('type','eof','text','','line',prev.line,'col',prev.col + numel(prev.text));
end
