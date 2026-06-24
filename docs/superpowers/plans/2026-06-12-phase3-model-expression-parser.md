# Phase 3 — Model-Expression Parser Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Parse the `model;…equations;…end;…end;` body of a gmod file into typed IR statements + equations so `gdsge.parser.parseFrontEnd` produces the **complete IR for HL1996**, matching the hand-authored reference (`buildHL1996IR`) via `gdsge.ir.isequalIR`.

**Architecture:** Four new single-purpose stages in `src/+gdsge/+parser/`: `tokenize` (text → token stream with positions) → `parseExpr` (precedence-climbing recursive descent → the 7 schema AST node kinds) → `parseModel` (statement classification + inline-reduction hoisting → `ir.model`) → `analyzeModel` (name resolution, square-system, interp arity, transRef checks). `parseFrontEnd` wires them in; the Phase-1 schema/validator/node API are untouched.

**Tech Stack:** MATLAB R2025b (`C:\Program Files\MATLAB\R2025b\bin\matlab.exe`), `matlab.unittest` via `tests/run.ps1`. No Python.

**Approved spec:** `docs/superpowers/specs/2026-06-12-phase3-model-expression-parser-design.md`

---

## Conventions used by every task

- **Run a single test class** (from the repo root `D:\refactor_gdsge`, PowerShell):

  ```powershell
  & 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe' -batch "addpath('src'); r = runtests('tests/parser/tTokenize.m'); exit(double(any([r.Failed])))"
  ```

  Exit code 0 = pass, 1 = failure. Swap the file path per task.

- **Run the full suite:** `pwsh -File tests/run.ps1` (exit 0 = all green).

- **MATLAB path policy:** never `addpath` persistently. The commands above add `src/` only inside a one-shot `matlab -batch` process — that is the project rule.

- **Token struct** (produced by `tokenize`, consumed by everything else):
  `struct('type',T,'text',S,'line',L,'col',C)` where `type` ∈
  `'name' | 'number' | 'op' | 'punct' | 'prime' | 'eof'`. The stream always ends
  with one `'eof'` sentinel — parsers peek freely without bounds checks.

- **Statement structs** (must match `gdsge.ir.schema` exactly — same field names and order as the fixture builders in `buildHL1996IR.m`):
  - `struct('type','assign','target',T,'primed',P,'expr',AST)`
  - `struct('type','interpCall','targets',{C},'primed',1,'args',{C},'interpRef','GDSGE_INTERP_VEC')`
  - `struct('type','reduction','kind',K,'target',T,'body',AST,'transRef',R)` with K ∈ EXPECT/MIN/MAX/PROD
  - equations: `struct('expr',AST,'primed',P)`

### Task 0: branch

- [ ] **Step 1: Create the phase branch**

```powershell
git checkout -b phase3-expression-parser
```

---

### Task 1: `tokenize.m` — model-body tokenizer

**Files:**
- Create: `src/+gdsge/+parser/tokenize.m`
- Test: `tests/parser/tTokenize.m`

- [ ] **Step 1: Write the failing test**

Create `tests/parser/tTokenize.m`:

```matlab
classdef tTokenize < matlab.unittest.TestCase
    methods (Test)
        function namesNumbersOps(tc)
            toks = gdsge.parser.tokenize('es1 = beta*c1 + 1.5e-2;');
            texts = {toks.text};
            types = {toks.type};
            tc.verifyEqual(texts(1:8), {'es1','=','beta','*','c1','+','1.5e-2',';'});
            tc.verifyEqual(types{1}, 'name');
            tc.verifyEqual(types{2}, 'punct');   % '=' is punctuation, not an expression op
            tc.verifyEqual(types{4}, 'op');
            tc.verifyEqual(types{7}, 'number');
            tc.verifyEqual(types{end}, 'eof');
        end
        function elementwiseAndComparisonOps(tc)
            toks = gdsge.parser.tokenize('a.*b ./ c .^ d == e ~= f <= g && h || i');
            ops = {toks(strcmp({toks.type},'op')).text};
            tc.verifyEqual(ops, {'.*','./','.^','==','~=','<=','&&','||'});
        end
        function primesAndReductionBraces(tc)
            toks = gdsge.parser.tokenize('GDSGE_EXPECT{g''^2 | my_trans}');
            tc.verifyEqual(toks(1).type, 'name');
            tc.verifyEqual(toks(1).text, 'GDSGE_EXPECT');
            tc.verifyEqual(toks(2).text, '{');
            tc.verifyEqual(toks(4).type, 'prime');
            tc.verifyEqual(sum(strcmp({toks.text}, '|')), 1);
            tc.verifyEqual(toks(end-1).text, '}');
        end
        function numberForms(tc)
            forms = {'1','0.95','1e-3','1.5e+2','.5','3.'};
            for i = 1:numel(forms)
                toks = gdsge.parser.tokenize(forms{i});
                tc.verifyEqual(toks(1).type, 'number', forms{i});
                tc.verifyEqual(toks(1).text, forms{i});
                tc.verifyEqual(toks(2).type, 'eof', forms{i});
            end
        end
        function dotBeforeStarIsElementwiseNotNumber(tc)
            toks = gdsge.parser.tokenize('2.*x');
            tc.verifyEqual({toks(1:3).text}, {'2','.*','x'});
        end
        function lineColPositionsWithOffset(tc)
            toks = gdsge.parser.tokenize(sprintf('a +\nbb*c'), 10);
            tc.verifyEqual([toks(1).line toks(1).col], [10 1]);   % a
            tc.verifyEqual([toks(3).line toks(3).col], [11 1]);   % bb
            tc.verifyEqual([toks(5).line toks(5).col], [11 4]);   % c
        end
        function badCharacterRaises(tc)
            tc.verifyError(@() gdsge.parser.tokenize('a $ b'), 'gdsge:parser:badToken');
        end
    end
end
```

- [ ] **Step 2: Run it to make sure it fails**

```powershell
& 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe' -batch "addpath('src'); r = runtests('tests/parser/tTokenize.m'); exit(double(any([r.Failed])))"
```

Expected: FAIL — `gdsge.parser.tokenize` undefined; exit code 1.

- [ ] **Step 3: Implement `tokenize.m`**

Create `src/+gdsge/+parser/tokenize.m`:

```matlab
function toks = tokenize(text, lineOffset)
% TOKENIZE  Model-body text -> token struct array.
%   TOKS = TOKENIZE(TEXT[, LINEOFFSET]) returns a struct array with fields
%   type/text/line/col, terminated by one 'eof' sentinel. Types:
%     'name' | 'number' | 'op' | 'punct' | 'prime' | 'eof'
%   '=' is 'punct' (assignment is statement-level, not an expression op).
%   Elementwise ops ('.*','./','.^') are kept distinct here; parseExpr
%   normalizes them to scalar forms. In gmod model bodies a quote is always
%   the future marker (never transpose / string), so no ambiguity exists.
%   LINEOFFSET (default 1) makes lines block-relative when callers tokenize
%   one statement at a time.
if nargin < 2; lineOffset = 1; end
toks = struct('type',{},'text',{},'line',{},'col',{});
nl = sprintf('\n'); cr = sprintf('\r'); tab = sprintf('\t');
twoCharOps = {'==','~=','<=','>=','&&','||','.*','./','.^'};
oneCharOps = '+-*/^<>&|~';
punct = '()[]{},;=';
line = lineOffset; col = 1;
i = 1; n = numel(text);
while i <= n
    ch = text(i);
    if ch == nl
        line = line + 1; col = 1; i = i + 1;
    elseif ch == ' ' || ch == tab || ch == cr
        col = col + 1; i = i + 1;
    elseif isletter(ch)
        j = i;
        while j <= n && (isletter(text(j)) || isDigitChar(text(j)) || text(j) == '_')
            j = j + 1;
        end
        toks = push(toks, 'name', text(i:j-1), line, col);
        col = col + (j - i); i = j;
    elseif isDigitChar(ch) || (ch == '.' && i < n && isDigitChar(text(i+1)))
        j = i;
        while j <= n && isDigitChar(text(j)); j = j + 1; end
        % a '.' belongs to the number unless it starts an elementwise op
        if j <= n && text(j) == '.' && ~(j < n && any(text(j+1) == '*/^'))
            j = j + 1;
            while j <= n && isDigitChar(text(j)); j = j + 1; end
        end
        if j <= n && (text(j) == 'e' || text(j) == 'E')
            k = j + 1;
            if k <= n && (text(k) == '+' || text(k) == '-'); k = k + 1; end
            if k <= n && isDigitChar(text(k))
                j = k;
                while j <= n && isDigitChar(text(j)); j = j + 1; end
            end
        end
        toks = push(toks, 'number', text(i:j-1), line, col);
        col = col + (j - i); i = j;
    elseif i < n && ismember(text(i:i+1), twoCharOps)
        toks = push(toks, 'op', text(i:i+1), line, col);
        col = col + 2; i = i + 2;
    elseif ch == ''''
        toks = push(toks, 'prime', '''', line, col);
        col = col + 1; i = i + 1;
    elseif any(ch == punct)
        toks = push(toks, 'punct', ch, line, col);
        col = col + 1; i = i + 1;
    elseif any(ch == oneCharOps)
        toks = push(toks, 'op', ch, line, col);
        col = col + 1; i = i + 1;
    else
        error('gdsge:parser:badToken', ...
            'Unexpected character "%s" at line %d, col %d', ch, line, col);
    end
end
toks = push(toks, 'eof', '', line, col);
end

function toks = push(toks, type, text, line, col)
toks(end+1) = struct('type',type,'text',text,'line',line,'col',col);
end

function tf = isDigitChar(c)
tf = c >= '0' && c <= '9';
end
```

- [ ] **Step 4: Run the test and make sure it passes**

Same command as Step 2. Expected: PASS, exit code 0.

- [ ] **Step 5: Commit**

```powershell
git add src/+gdsge/+parser/tokenize.m tests/parser/tTokenize.m
git commit -m "feat(parser): tokenize - model-body tokenizer with positions"
```

---

### Task 2: `parseExpr.m` — precedence-climbing expression parser

**Files:**
- Create: `src/+gdsge/+parser/parseExpr.m`
- Test: `tests/parser/tParseExpr.m`

MATLAB precedence facts the implementation must honor (from the spec):
all binary ops left-associative **including `^`**; unary `+ - ~` bind between
power and multiplicative (`-2^2 = -4`, `-a*b = (-a)*b`); power accepts a unary
right operand (`2^-3`); elementwise ops normalize to scalar in the AST.

- [ ] **Step 1: Write the failing test**

Create `tests/parser/tParseExpr.m`:

```matlab
classdef tParseExpr < matlab.unittest.TestCase
    methods (Test)
        function valueMatchesMatlabEval(tc)
            % property test: our AST evaluates exactly like MATLAB on scalars
            a = 2; b = 3; c = 5; %#ok<NASGU> % consumed via eval below
            env = struct('a', 2, 'b', 3, 'c', 5);
            exprs = {'a+b*c','(a+b)*c','-a^b','a^-b','a^b^2','a-b-c','a/b/c', ...
                     '-a*b','a*-b','2^-3','-2^2','a+b-c*a/b^c', ...
                     'a==b','a<=b','a&&b||c==c','a~=b','a>b'};
            for i = 1:numel(exprs)
                ast = pe(exprs{i});
                tc.verifyEqual(evalAst(ast, env), double(eval(exprs{i})), ...
                    sprintf('mismatch for "%s"', exprs{i}), 'AbsTol', 1e-12);
            end
        end
        function powerIsLeftAssociative(tc)
            N = @gdsge.ir.node.num; B = @gdsge.ir.node.binop;
            tc.verifyEqual(pe('2^3^2'), B('^', B('^', N(2), N(3)), N(2)));
        end
        function unaryBindsBetweenPowerAndMultiplicative(tc)
            NM = @gdsge.ir.node.name; B = @gdsge.ir.node.binop; U = @gdsge.ir.node.unop;
            tc.verifyEqual(pe('-a*b'), B('*', U('-', NM('a')), NM('b')));
            tc.verifyEqual(pe('-a^b'), U('-', B('^', NM('a'), NM('b'))));
            tc.verifyEqual(pe('a^-b'), B('^', NM('a'), U('-', NM('b'))));
        end
        function elementwiseOpsNormalize(tc)
            NM = @gdsge.ir.node.name; B = @gdsge.ir.node.binop;
            tc.verifyEqual(pe('a.*b'), B('*', NM('a'), NM('b')));
            tc.verifyEqual(pe('a./b'), B('/', NM('a'), NM('b')));
            tc.verifyEqual(pe('a.^b'), B('^', NM('a'), NM('b')));
        end
        function primesCallsParens(tc)
            NM = @gdsge.ir.node.name; P = @gdsge.ir.node.primed;
            B = @gdsge.ir.node.binop; C = @gdsge.ir.node.call;
            tc.verifyEqual(pe('g'''), P('g'));
            tc.verifyEqual(pe('exp(a+b)'), C('exp', {B('+', NM('a'), NM('b'))}));
            tc.verifyEqual(pe('f(a,b)'), C('f', {NM('a'), NM('b')}));
        end
        function fixtureEs1Body(tc)
            % the exact es1 reduction body from buildHL1996IR
            N = @gdsge.ir.node.num; NM = @gdsge.ir.node.name; P = @gdsge.ir.node.primed;
            B = @gdsge.ir.node.binop; U = @gdsge.ir.node.unop;
            gPow = B('^', P('g'), B('-', N(1), NM('gamma')));
            psnD = B('+', P('psn'), P('d'));
            want = B('/', B('*', B('*', gPow, ...
                B('^', B('/', P('c1n'), NM('c1')), U('-', NM('gamma')))), psnD), NM('ps'));
            tc.verifyEqual(pe('g''^(1-gamma)*(c1n''/c1)^(-gamma)*(psn''+d'')/ps'), want);
        end
        function fixtureBudget1(tc)
            N = @gdsge.ir.node.num; NM = @gdsge.ir.node.name; B = @gdsge.ir.node.binop;
            want = B('-', B('-', B('-', ...
                B('+', B('*', NM('w1'), B('+', NM('ps'), NM('d'))), NM('eta1')), ...
                NM('c1')), B('*', NM('ps'), NM('s1p'))), B('*', NM('pb'), NM('b1p')));
            tc.verifyEqual(pe('w1*(ps+d)+eta1 - c1 - ps*s1p - pb*b1p'), want);
        end
        function fixtureW1Consis(tc)
            NM = @gdsge.ir.node.name; P = @gdsge.ir.node.primed; B = @gdsge.ir.node.binop;
            psnD = B('+', P('psn'), P('d'));
            want = B('-', B('/', B('+', B('*', NM('s1p'), psnD), ...
                B('/', NM('b1p'), P('g'))), psnD), P('w1n'));
            tc.verifyEqual(pe('(s1p*(psn''+d'') + b1p/g'')/(psn''+d'') - w1n'''), want);
        end
        function badExprRaises(tc)
            tc.verifyError(@() pe('a+'), 'gdsge:parser:badExpr');
            tc.verifyError(@() pe('(a'), 'gdsge:parser:badExpr');
            tc.verifyError(@() pe('a + {b}'), 'gdsge:parser:badExpr');
        end
        function primedCallIsPhase7(tc)
            tc.verifyError(@() pe('ps_future''(x)'), 'gdsge:parser:badExpr');
        end
    end
end

% ----- local helpers -------------------------------------------------------
function ast = pe(txt)
toks = gdsge.parser.tokenize(txt);
[ast, pos] = gdsge.parser.parseExpr(toks);
if ~strcmp(toks(pos).type, 'eof')
    error('gdsge:parser:badExpr', 'expression not fully consumed at "%s"', toks(pos).text);
end
end

function v = evalAst(n, env)
switch n.kind
    case 'num';  v = n.value;
    case 'name'; v = env.(n.id);
    case 'unop'
        a = evalAst(n.arg, env);
        switch n.op
            case '-'; v = -a;
            case '+'; v = a;
            case '~'; v = ~a;
        end
    case 'binop'
        l = evalAst(n.lhs, env); r = evalAst(n.rhs, env);
        switch n.op
            case '+';  v = l + r;   case '-';  v = l - r;
            case '*';  v = l * r;   case '/';  v = l / r;
            case '^';  v = l ^ r;
            case '=='; v = l == r;  case '~='; v = l ~= r;
            case '<';  v = l < r;   case '>';  v = l > r;
            case '<='; v = l <= r;  case '>='; v = l >= r;
            case '&&'; v = l && r;  case '||'; v = l || r;
            case '&';  v = l & r;   case '|';  v = l | r;
        end
end
v = double(v);
end
```

- [ ] **Step 2: Run it to make sure it fails**

```powershell
& 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe' -batch "addpath('src'); r = runtests('tests/parser/tParseExpr.m'); exit(double(any([r.Failed])))"
```

Expected: FAIL — `gdsge.parser.parseExpr` undefined; exit code 1.

- [ ] **Step 3: Implement `parseExpr.m`**

Create `src/+gdsge/+parser/parseExpr.m`:

```matlab
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
            'calling a primed object (named-interp call) is not supported until Phase 7 (line %d, col %d)', ...
            t.line, t.col);
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
```

- [ ] **Step 4: Run the test and make sure it passes**

Same command as Step 2. Expected: PASS, exit code 0.

- [ ] **Step 5: Run the tokenizer test too (no regression)**

```powershell
& 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe' -batch "addpath('src'); r = runtests('tests/parser'); exit(double(any([r.Failed])))"
```

Expected: PASS (all classes under `tests/parser`).

- [ ] **Step 6: Commit**

```powershell
git add src/+gdsge/+parser/parseExpr.m tests/parser/tParseExpr.m
git commit -m "feat(parser): parseExpr - precedence-climbing expression parser (MATLAB semantics)"
```

---

### Task 3: `splitStatements` gains statement start lines

`parseModel` needs block-relative line numbers for error messages, but
`splitStatements` discards them. Add an optional second output — existing
callers are unaffected.

**Files:**
- Modify: `src/+gdsge/+parser/splitStatements.m`
- Test: `tests/parser/tSplitStatements.m` (add one method)

- [ ] **Step 1: Add the failing test method**

In `tests/parser/tSplitStatements.m`, add inside `methods (Test)`:

```matlab
        function reportsStatementStartLines(tc)
            txt = sprintf('a = 1;\n\nM = [1 2\n3 4];\nb = 2; c = 3;');
            [s, lines] = gdsge.parser.splitStatements(txt);
            tc.verifyEqual(numel(s), 4);
            tc.verifyEqual(lines, [1 3 5 5]);   % M starts line 3; b and c share line 5
        end
```

- [ ] **Step 2: Run it to make sure it fails**

```powershell
& 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe' -batch "addpath('src'); r = runtests('tests/parser/tSplitStatements.m'); exit(double(any([r.Failed])))"
```

Expected: FAIL — "Too many output arguments"; exit code 1.

- [ ] **Step 3: Implement**

Replace the full contents of `src/+gdsge/+parser/splitStatements.m` with:

```matlab
function [stmts, startLines] = splitStatements(text)
% SPLITSTATEMENTS  Split MATLAB-ish text into logical statements. Boundaries are
%   ';' or newline at bracket-depth 0; inside ( ) [ ] { } both are ignored, and
%   newlines inside brackets are PRESERVED (they are matrix row separators).
%   Empty statements are dropped. Returns a row cellstr of trimmed statements
%   and, optionally, each statement's 1-based start line within TEXT.
stmts = {};
startLines = [];
depth = 0;
cur = '';
curStart = 1;
line = 1;
nl = sprintf('\n');
cr = sprintf('\r');
for i = 1:numel(text)
    ch = text(i);
    if isempty(strtrim(cur)) && ~isspace(ch)
        curStart = line;   % first visible char of a fresh statement
    end
    switch ch
        case {'[','(','{'}
            depth = depth + 1; cur(end+1) = ch; %#ok<AGROW>
        case {']',')','}'}
            depth = max(0, depth - 1); cur(end+1) = ch; %#ok<AGROW>
        case ';'
            if depth == 0
                [stmts, startLines, cur] = flush(stmts, startLines, cur, curStart);
            else
                cur(end+1) = ch; %#ok<AGROW>
            end
        case {nl, cr}
            if ch == nl; line = line + 1; end
            if depth == 0
                [stmts, startLines, cur] = flush(stmts, startLines, cur, curStart);
            else
                cur(end+1) = nl; %#ok<AGROW>
            end
        otherwise
            cur(end+1) = ch; %#ok<AGROW>
    end
end
[stmts, startLines, ~] = flush(stmts, startLines, cur, curStart);
end

function [stmts, startLines, cur] = flush(stmts, startLines, cur, curStart)
t = strtrim(cur);
if ~isempty(t)
    stmts{end+1} = t;
    startLines(end+1) = curStart;
end
cur = '';
end
```

- [ ] **Step 4: Run the test class — all five methods must pass**

Same command as Step 2. Expected: PASS, exit code 0.

- [ ] **Step 5: Commit**

```powershell
git add src/+gdsge/+parser/splitStatements.m tests/parser/tSplitStatements.m
git commit -m "feat(parser): splitStatements optionally reports statement start lines"
```

---

### Task 4: `parseModel.m` — statement classification + reduction hoisting

**Files:**
- Create: `src/+gdsge/+parser/parseModel.m`
- Test: `tests/parser/tParseModel.m`

- [ ] **Step 1: Write the failing test**

Create `tests/parser/tParseModel.m`:

```matlab
classdef tParseModel < matlab.unittest.TestCase
    methods (Test)
        function plainAssign(tc)
            m = gdsge.parser.parseModel(body('b1p = nb1p + Kb;'));
            tc.verifyEqual(numel(m.statements), 1);
            s = m.statements{1};
            B = @gdsge.ir.node.binop; NM = @gdsge.ir.node.name;
            tc.verifyEqual(s, struct('type','assign','target','b1p','primed',0, ...
                'expr', B('+', NM('nb1p'), NM('Kb'))));
        end
        function primedAssign(tc)
            m = gdsge.parser.parseModel(body('w1_consis'' = s1p - w1n'';'));
            s = m.statements{1};
            tc.verifyEqual(s.type, 'assign');
            tc.verifyEqual(s.target, 'w1_consis');
            tc.verifyEqual(s.primed, 1);
        end
        function wholeRhsReductionBecomesReductionStmt(tc)
            m = gdsge.parser.parseModel(body('es1 = GDSGE_EXPECT{g''^2};'));
            s = m.statements{1};
            N = @gdsge.ir.node.num; P = @gdsge.ir.node.primed; B = @gdsge.ir.node.binop;
            tc.verifyEqual(s, struct('type','reduction','kind','EXPECT', ...
                'target','es1','body',B('^', P('g'), N(2)),'transRef','shock_trans'));
        end
        function reductionWithCustomTransition(tc)
            m = gdsge.parser.parseModel(body('lo = GDSGE_MIN{c1 | my_trans};'));
            s = m.statements{1};
            tc.verifyEqual(s.kind, 'MIN');
            tc.verifyEqual(s.transRef, 'my_trans');
        end
        function inlineReductionIsHoisted(tc)
            m = gdsge.parser.parseModel(body('equity_premium = GDSGE_EXPECT{g''} - 1/pb;'));
            tc.verifyEqual(numel(m.statements), 2);
            red = m.statements{1}; host = m.statements{2};
            N = @gdsge.ir.node.num; NM = @gdsge.ir.node.name;
            P = @gdsge.ir.node.primed; B = @gdsge.ir.node.binop;
            tc.verifyEqual(red, struct('type','reduction','kind','EXPECT', ...
                'target','GDSGE_EXPECT_1','body',P('g'),'transRef','shock_trans'));
            tc.verifyEqual(host, struct('type','assign','target','equity_premium', ...
                'primed',0,'expr',B('-', NM('GDSGE_EXPECT_1'), B('/', N(1), NM('pb')))));
        end
        function hoistCounterIsGlobalAcrossKinds(tc)
            m = gdsge.parser.parseModel(body(sprintf( ...
                'x = GDSGE_EXPECT{a} + 1;\ny = GDSGE_MAX{b} + 2;')));
            tc.verifyEqual(m.statements{1}.target, 'GDSGE_EXPECT_1');
            tc.verifyEqual(m.statements{3}.target, 'GDSGE_MAX_2');
        end
        function interpCallForm(tc)
            m = gdsge.parser.parseModel(body('[psn'',pbn''] = GDSGE_INTERP_VEC''(w1n'');'));
            s = m.statements{1};
            P = @gdsge.ir.node.primed;
            tc.verifyEqual(s, struct('type','interpCall','targets',{{'psn','pbn'}}, ...
                'primed',1,'args',{{P('w1n')}},'interpRef','GDSGE_INTERP_VEC'));
        end
        function equationsPlainAndPrimed(tc)
            txt = sprintf('a = 1;\nequations;\n  -1+beta*es1+ms1;\n  w1_consis'';\nend;');
            m = gdsge.parser.parseModel(txt);
            tc.verifyEqual(numel(m.equations), 2);
            tc.verifyEqual(m.equations{1}.primed, 0);
            tc.verifyEqual(m.equations{2}, ...
                struct('expr', gdsge.ir.node.name('w1_consis'), 'primed', 1));
        end
        function errorOnMultipleReductions(tc)
            tc.verifyError(@() gdsge.parser.parseModel( ...
                body('x = GDSGE_EXPECT{a} + GDSGE_EXPECT{b};')), ...
                'gdsge:parser:multipleReductions');
        end
        function errorOnNestedReduction(tc)
            tc.verifyError(@() gdsge.parser.parseModel( ...
                body('x = GDSGE_EXPECT{GDSGE_MAX{a}};')), ...
                'gdsge:parser:nestedReduction');
        end
        function errorOnUnterminatedReduction(tc)
            tc.verifyError(@() gdsge.parser.parseModel( ...
                body('x = GDSGE_EXPECT{a;')), ...
                'gdsge:parser:unterminatedReduction');
        end
        function errorOnPrimedReductionTarget(tc)
            tc.verifyError(@() gdsge.parser.parseModel( ...
                body('x'' = GDSGE_EXPECT{a};')), ...
                'gdsge:parser:badStatement');
        end
        function errorOnMissingEquationsBlock(tc)
            tc.verifyError(@() gdsge.parser.parseModel('x = 1;'), ...
                'gdsge:parser:missingEquations');
        end
        function errorOnContentAfterEquations(tc)
            txt = sprintf('a = 1;\nequations;\n  x;\nend;\nb = 2;');
            tc.verifyError(@() gdsge.parser.parseModel(txt), ...
                'gdsge:parser:badStatement');
        end
    end
end

% ----- local helper --------------------------------------------------------
function txt = body(stmts)
% wrap model statements with a minimal equations block (parseModel requires one)
txt = sprintf('%s\nequations;\n  x;\nend;', stmts);
end
```

- [ ] **Step 2: Run it to make sure it fails**

```powershell
& 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe' -batch "addpath('src'); r = runtests('tests/parser/tParseModel.m'); exit(double(any([r.Failed])))"
```

Expected: FAIL — `gdsge.parser.parseModel` undefined; exit code 1.

- [ ] **Step 3: Implement `parseModel.m`**

Create `src/+gdsge/+parser/parseModel.m`:

```matlab
function model = parseModel(modelText)
% PARSEMODEL  model-block text -> ir.model struct:
%   .statements  ordered cell of assign / interpCall / reduction statements
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
    toks = gdsge.parser.tokenize(bodyStmts{i}, bodyLines(i));
    [stmtList, counter] = parseStatement(toks, counter);
    statements = [statements, stmtList]; %#ok<AGROW>
end

% --- 3. equations -----------------------------------------------------------
[eqStmts, eqLines] = gdsge.parser.splitStatements(eqText);
equations = cell(1, numel(eqStmts));
for i = 1:numel(eqStmts)
    toks = gdsge.parser.tokenize(eqStmts{i}, eqLines(i) + eqOpen);
    equations{i} = parseEquation(toks);
end

model = struct('statements', {statements}, 'equations', {equations});
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

[redIdx, kinds] = findReductions(rest);
if isempty(redIdx)
    [expr, q] = gdsge.parser.parseExpr(rest, 1);
    expectEof(rest, q);
    stmtList = { struct('type','assign','target',target,'primed',primed,'expr',expr) };
    return
end

closeIdx = matchBrace(rest, redIdx(1) + 1);
nested = redIdx(redIdx > redIdx(1) & redIdx < closeIdx);
if ~isempty(nested)
    error('gdsge:parser:nestedReduction', ...
        'reductions cannot nest (line %d); use an intermediate variable', ...
        rest(nested(1)).line);
end
if numel(redIdx) > 1
    error('gdsge:parser:multipleReductions', ...
        'only one GDSGE_EXPECT/MIN/MAX/PROD per statement (line %d); use an intermediate variable', ...
        rest(redIdx(2)).line);
end
kind = kinds{1};
[bodyAst, transRef] = parseReductionInner(rest(redIdx(1)+2:closeIdx-1), rest(redIdx(1)).line);

if redIdx(1) == 1 && strcmp(rest(closeIdx+1).type, 'eof')
    % the whole RHS is the reduction: it becomes the statement itself
    if primed
        error('gdsge:parser:badStatement', ...
            'a primed target cannot take a reduction directly (line %d)', toks(1).line);
    end
    stmtList = { struct('type','reduction','kind',kind,'target',target, ...
                        'body',bodyAst,'transRef',transRef) };
    return
end

% inline: hoist under a generated name placed just before the host statement
counter = counter + 1;
genName = sprintf('GDSGE_%s_%d', kind, counter);
nameTok = rest(redIdx(1)); nameTok.type = 'name'; nameTok.text = genName;
host = [rest(1:redIdx(1)-1), nameTok, rest(closeIdx+1:end)];
[expr, q] = gdsge.parser.parseExpr(host, 1);
expectEof(host, q);
stmtList = { struct('type','reduction','kind',kind,'target',genName, ...
                    'body',bodyAst,'transRef',transRef), ...
             struct('type','assign','target',target,'primed',primed,'expr',expr) };
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
        'unprimed GDSGE_INTERP_VEC is not supported until Phase 7 (line %d, col %d)', ...
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
```

- [ ] **Step 4: Run the test and make sure it passes**

Same command as Step 2. Expected: PASS, exit code 0.

- [ ] **Step 5: Commit**

```powershell
git add src/+gdsge/+parser/parseModel.m tests/parser/tParseModel.m
git commit -m "feat(parser): parseModel - statement classification + reduction hoisting"
```

---

### Task 5: fixture rename `__ep_expect` → `GDSGE_EXPECT_1` + regenerate JSON golden

The hand-authored reference IR used `__ep_expect` for the hoisted reduction temp.
The parser generates `GDSGE_EXPECT_1` (deterministic; `GDSGE_` is the reserved
namespace; `__`-prefixed identifiers are reserved in C++). Rename the fixture
and regenerate the JSON golden — a pure rename, per the approved spec §6.

**Files:**
- Modify: `tests/HeatonLucas1996/ir/buildHL1996IR.m` (2 lines)
- Regenerate: `tests/HeatonLucas1996/ir/HL1996.gdsge.json`

- [ ] **Step 1: Edit the fixture**

In `tests/HeatonLucas1996/ir/buildHL1996IR.m` replace:

```matlab
    stReduce('EXPECT','__ep_expect', epBody), ... % __ep_expect: parser-hoisted reduction temp (generated name, not a policy var)
    stAssign('equity_premium', 0, B('-', NM('__ep_expect'), B('/', N(1), NM('pb')))) };
```

with:

```matlab
    stReduce('EXPECT','GDSGE_EXPECT_1', epBody), ... % GDSGE_EXPECT_1: parser-hoisted reduction temp (kind + global counter)
    stAssign('equity_premium', 0, B('-', NM('GDSGE_EXPECT_1'), B('/', N(1), NM('pb')))) };
```

- [ ] **Step 2: Confirm the golden is now stale (test fails)**

```powershell
& 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe' -batch "addpath('src'); r = runtests('tests/HeatonLucas1996/ir/tIrHL1996.m'); exit(double(any([r.Failed])))"
```

Expected: FAIL — `encodingMatchesGolden` reports "HL1996.gdsge.json is stale"; exit code 1.

- [ ] **Step 3: Regenerate the JSON golden**

```powershell
& 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe' -batch "addpath('src'); addpath('tests/HeatonLucas1996/ir'); txt = gdsge.ir.encode(buildHL1996IR()); fid = fopen('tests/HeatonLucas1996/ir/HL1996.gdsge.json','w'); fwrite(fid, txt); fclose(fid); disp('golden regenerated')"
```

Expected: prints `golden regenerated`.

- [ ] **Step 4: Verify the IR fixture suite is green again**

Same command as Step 2. Expected: PASS, exit code 0.

- [ ] **Step 5: Verify no other `__ep_expect` references remain in code**

```powershell
git grep -n "__ep_expect" -- src tests
```

Expected: no matches (git grep exits 1 on zero hits — that is the pass condition;
docs/ still mentions the old name when describing this very rename, which is fine).

- [ ] **Step 6: Commit**

```powershell
git add tests/HeatonLucas1996/ir/buildHL1996IR.m tests/HeatonLucas1996/ir/HL1996.gdsge.json
git commit -m "test(ir): rename hoisted temp __ep_expect -> GDSGE_EXPECT_1; regen JSON golden"
```

---

### Task 6: `analyzeModel.m` — semantic analysis

**Files:**
- Create: `src/+gdsge/+parser/analyzeModel.m`
- Test: `tests/parser/tAnalyzeModel.m`

- [ ] **Step 1: Write the failing test**

Create `tests/parser/tAnalyzeModel.m`:

```matlab
classdef tAnalyzeModel < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addIrFixtureFolder(tc)
            here = fileparts(mfilename('fullpath'));            % tests/parser
            irDir = fullfile(fileparts(here), 'HeatonLucas1996', 'ir');
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(irDir));
        end
    end
    methods (Test)
        function referenceIrPasses(tc)
            tc.verifyWarningFree(@() gdsge.parser.analyzeModel(buildHL1996IR()));
        end
        function unknownNameInEquation(tc)
            ir = buildHL1996IR();
            ir.model.equations{1}.expr = gdsge.ir.node.name('zzz_undefined');
            tc.verifyError(@() gdsge.parser.analyzeModel(ir), 'gdsge:parser:unknownName');
        end
        function unknownNameInStatement(tc)
            ir = buildHL1996IR();
            ir.model.statements{6}.expr = gdsge.ir.node.name('zzz_undefined');
            tc.verifyError(@() gdsge.parser.analyzeModel(ir), 'gdsge:parser:unknownName');
        end
        function paramCannotBePrimed(tc)
            ir = buildHL1996IR();
            % statements{2} is the es1 reduction; beta is a param, not primeable
            ir.model.statements{2}.body = gdsge.ir.node.primed('beta');
            tc.verifyError(@() gdsge.parser.analyzeModel(ir), 'gdsge:parser:unknownName');
        end
        function localUsedBeforeDefinitionFails(tc)
            ir = buildHL1996IR();
            % statements{6} (b1p = nb1p + Kb) referencing budget_1, defined later
            ir.model.statements{6}.expr = gdsge.ir.node.name('budget_1');
            tc.verifyError(@() gdsge.parser.analyzeModel(ir), 'gdsge:parser:unknownName');
        end
        function nonSquareSystem(tc)
            ir = buildHL1996IR();
            ir.model.equations(end) = [];   % drop the primed w1_consis: 11 vs 19
            tc.verifyError(@() gdsge.parser.analyzeModel(ir), 'gdsge:parser:notSquare');
        end
        function interpArityMismatch(tc)
            ir = buildHL1996IR();
            ir.model.statements{1}.targets(end) = [];   % 3 targets, 4 var_interp
            tc.verifyError(@() gdsge.parser.analyzeModel(ir), 'gdsge:parser:interpArity');
        end
        function badTransRef(tc)
            ir = buildHL1996IR();
            ir.model.statements{2}.transRef = 'no_such_trans';
            tc.verifyError(@() gdsge.parser.analyzeModel(ir), 'gdsge:parser:badTransRef');
        end
    end
end
```

- [ ] **Step 2: Run it to make sure it fails**

```powershell
& 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe' -batch "addpath('src'); r = runtests('tests/parser/tAnalyzeModel.m'); exit(double(any([r.Failed])))"
```

Expected: FAIL — `gdsge.parser.analyzeModel` undefined; exit code 1.

- [ ] **Step 3: Implement `analyzeModel.m`**

Create `src/+gdsge/+parser/analyzeModel.m`:

```matlab
function analyzeModel(ir)
% ANALYZEMODEL  Model-level semantic checks on a fully assembled IR. The IR's
%   *shape* is gdsge.ir.validate's job; this checks what the validator cannot:
%   1. every name leaf resolves (param/state/shock/policy/aux/interp-output/
%      reduction target/earlier local/builtin); primed names are shocks,
%      interp outputs, or policy variables
%   2. the system is square after expanding primed equations by shock count
%   3. the interp call has one target per declared var_interp
%   4. every reduction transRef names a declared transition matrix
%   Call-node function names (exp, log, ...) are not variable references and
%   are not resolved here.
builtins = {'GDSGE_Iter','TASK','shock'};
known = builtins;
for i = 1:numel(ir.params); known{end+1} = ir.params{i}.name; end %#ok<AGROW>
known = [known, reshape(ir.states.names, 1, []), reshape(ir.shocks.names, 1, [])];
policyNames = cell(1, numel(ir.variables.policy));
for i = 1:numel(ir.variables.policy); policyNames{i} = ir.variables.policy{i}.name; end
auxNames = cell(1, numel(ir.variables.aux));
for i = 1:numel(ir.variables.aux); auxNames{i} = ir.variables.aux{i}.name; end
known = [known, policyNames, auxNames];
primeable = [reshape(ir.shocks.names, 1, []), policyNames];  % + interp outputs below

stmts = ir.model.statements;
for i = 1:numel(stmts)
    s = stmts{i};
    switch s.type
        case 'assign'
            checkNames(s.expr, known, primeable, sprintf('statement %d (assign %s)', i, s.target));
            known{end+1} = s.target; %#ok<AGROW>
        case 'reduction'
            if ~isfield(ir.shocks.transitions, s.transRef)
                error('gdsge:parser:badTransRef', ...
                    'statement %d: reduction transition "%s" is not a declared transition matrix', ...
                    i, s.transRef);
            end
            checkNames(s.body, known, primeable, sprintf('statement %d (reduction %s)', i, s.target));
            known{end+1} = s.target; %#ok<AGROW>
        case 'interpCall'
            nDecl = numel(ir.variables.interp);
            if numel(s.targets) ~= nDecl
                error('gdsge:parser:interpArity', ...
                    'statement %d: GDSGE_INTERP_VEC has %d targets but %d var_interp are declared', ...
                    i, numel(s.targets), nDecl);
            end
            for a = 1:numel(s.args)
                checkNames(s.args{a}, known, primeable, sprintf('statement %d (interp call)', i));
            end
            known = [known, reshape(s.targets, 1, [])]; %#ok<AGROW>
            primeable = [primeable, reshape(s.targets, 1, [])]; %#ok<AGROW>
    end
end

for i = 1:numel(ir.model.equations)
    checkNames(ir.model.equations{i}.expr, known, primeable, sprintf('equation %d', i));
end

% square system after prime expansion
nUnknown = 0;
for i = 1:numel(ir.variables.policy)
    nUnknown = nUnknown + ir.variables.policy{i}.length;
end
nEq = 0;
for i = 1:numel(ir.model.equations)
    if ir.model.equations{i}.primed
        nEq = nEq + ir.shocks.count;
    else
        nEq = nEq + 1;
    end
end
if nUnknown ~= nEq
    error('gdsge:parser:notSquare', ...
        'system is not square: %d unknowns vs %d equations (primed equations count shock_num = %d rows)', ...
        nUnknown, nEq, ir.shocks.count);
end
end

function checkNames(ast, known, primeable, where)
if strcmp(ast.kind, 'name') && ~ismember(ast.id, known)
    error('gdsge:parser:unknownName', '%s: unknown name "%s"', where, ast.id);
end
if strcmp(ast.kind, 'primed') && ~ismember(ast.id, primeable)
    error('gdsge:parser:unknownName', ...
        '%s: "%s''" is not a shock, interp output, or policy variable', where, ast.id);
end
kids = gdsge.ir.node.children(ast);
for i = 1:numel(kids)
    checkNames(kids{i}, known, primeable, where);
end
end
```

- [ ] **Step 4: Run the test and make sure it passes**

Same command as Step 2. Expected: PASS, exit code 0.

- [ ] **Step 5: Commit**

```powershell
git add src/+gdsge/+parser/analyzeModel.m tests/parser/tAnalyzeModel.m
git commit -m "feat(parser): analyzeModel - name resolution, square system, interp arity, transRef"
```

---

### Task 7: wire it into `parseFrontEnd` — full IR for HL1996 (the phase gate)

**Files:**
- Modify: `src/+gdsge/+parser/parseFrontEnd.m`
- Modify: `src/+gdsge/+parser/assemblePartialIR.m`
- Modify: `tests/HeatonLucas1996/parser/tFrontEndHL1996.m`
- Create: `tests/HeatonLucas1996/parser/tFullIRHL1996.m`

- [ ] **Step 1: Write the failing phase-gate test**

Create `tests/HeatonLucas1996/parser/tFullIRHL1996.m`:

```matlab
classdef tFullIRHL1996 < matlab.unittest.TestCase
    % Phase 3 gate: the parser must reproduce the hand-authored reference IR.
    methods (TestClassSetup)
        function addIrFolder(tc)
            here = fileparts(mfilename('fullpath'));            % tests/HeatonLucas1996/parser
            irDir = fullfile(fileparts(here), 'ir');            % tests/HeatonLucas1996/ir
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(irDir));
        end
    end
    methods (Test)
        function parserReproducesReferenceIR(tc)
            ir = parseHL1996();
            ref = buildHL1996IR();
            ref.options.numThreads = 0;   % machine-derived; neutralize
            ir.options.numThreads  = 0;
            tc.verifyTrue(gdsge.ir.isequalIR(ir, ref), ...
                'parsed HL1996 IR differs from the hand-authored reference');
        end
        function parsedIrSurvivesJsonRoundtrip(tc)
            ir = parseHL1996();
            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));
        end
    end
end

function ir = parseHL1996()
here = fileparts(mfilename('fullpath'));
gmodPath = fullfile(fileparts(here), 'HL1996.gmod');
ir = gdsge.parser.parseFrontEnd(fileread(gmodPath), 'HL1996');
end
```

- [ ] **Step 2: Run it to make sure it fails**

```powershell
& 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe' -batch "addpath('src'); r = runtests('tests/HeatonLucas1996/parser/tFullIRHL1996.m'); exit(double(any([r.Failed])))"
```

Expected: FAIL — `parserReproducesReferenceIR` fails (parsed `ir.model` is still empty); exit code 1.

- [ ] **Step 3: Update `assemblePartialIR.m` to accept the model section**

Replace the full contents of `src/+gdsge/+parser/assemblePartialIR.m` with:

```matlab
function ir = assemblePartialIR(modelName, decl, sim, options, hooks, model)
% ASSEMBLEPARTIALIR  Combine front-end sections into a schema-valid IR.
%   MODEL (statements + equations from parseModel) is optional; omitting it
%   yields the Phase-2 partial IR (empty model body). Enforces the
%   bounds-completeness semantic gate, then runs gdsge.ir.validate.
if nargin < 6
    model = struct('statements', {{}}, 'equations', {{}});
end
ir.irVersion = '1.0.0';
ir.modelName = modelName;
ir.params    = decl.params;
ir.options   = options;
ir.shocks    = decl.shocks;
ir.states    = decl.states;
ir.variables = decl.variables;
ir.bounds    = decl.bounds;
ir.interp    = decl.interp;
ir.model     = model;
ir.simulate  = sim;
ir.hooks     = hooks;

% semantic gate: every policy variable needs an inbound
boundNames = cell(1, numel(decl.bounds));
for i = 1:numel(decl.bounds); boundNames{i} = decl.bounds{i}.name; end
for i = 1:numel(decl.variables.policy)
    pn = decl.variables.policy{i}.name;
    if ~ismember(pn, boundNames)
        error('gdsge:parser:missingBound', 'Policy variable "%s" has no inbound', pn);
    end
end

r = gdsge.ir.validate(ir);
if ~r.pass
    error('gdsge:parser:invalidIR', 'IR failed validation:\n%s', ...
        strjoin(r.errors, sprintf('\n')));
end
end
```

- [ ] **Step 4: Update `parseFrontEnd.m`**

Replace the full contents of `src/+gdsge/+parser/parseFrontEnd.m` with:

```matlab
function ir = parseFrontEnd(gmodText, modelName)
% PARSEFRONTEND  gmod text + model name -> validated, semantically checked IR.
clean = gdsge.parser.preprocess(gmodText);
sb    = gdsge.parser.splitBlocks(clean);

decl    = gdsge.parser.parseDeclarations(sb.declText);
options = gdsge.parser.resolveOptions(decl.ws);

if ~isfield(sb.blocks, 'simulate')
    error('gdsge:parser:missingSimulate', 'No simulate block found');
end
sim = gdsge.parser.parseSimulate(sb.blocks.simulate);

if ~isfield(sb.blocks, 'model')
    error('gdsge:parser:missingModel', 'No model block found');
end
model = gdsge.parser.parseModel(sb.blocks.model);

hooks = hooksFromBlocks(sb.blocks);
ir = gdsge.parser.assemblePartialIR(modelName, decl, sim, options, hooks, model);
gdsge.parser.analyzeModel(ir);
end

function h = hooksFromBlocks(blocks)
h = struct('preModel','','preIter','','postIter','', ...
           'preJacCode','','postJacCode','','cxx','');
map = struct('pre_model','preModel', 'pre_iter','preIter', 'post_iter','postIter', ...
             'pre_jac_code','preJacCode', 'post_jac_code','postJacCode');
keys = fieldnames(map);
for i = 1:numel(keys)
    if isfield(blocks, keys{i})
        h.(map.(keys{i})) = blocks.(keys{i});
    end
end
end
```

- [ ] **Step 5: Update `tFrontEndHL1996.m` (the model body is no longer empty)**

Replace the full contents of `tests/HeatonLucas1996/parser/tFrontEndHL1996.m` with:

```matlab
classdef tFrontEndHL1996 < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addIrFolder(tc)
            here = fileparts(mfilename('fullpath'));            % tests/HeatonLucas1996/parser
            irDir = fullfile(fileparts(here), 'ir');            % tests/HeatonLucas1996/ir
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(irDir));
        end
    end
    methods (Test)
        function irValidatesAndMatchesReference(tc)
            here = fileparts(mfilename('fullpath'));
            gmodPath = fullfile(fileparts(here), 'HL1996.gmod');
            ir = gdsge.parser.parseFrontEnd(fileread(gmodPath), 'HL1996');

            % 1. it validates
            r = gdsge.ir.validate(ir);
            tc.verifyTrue(r.pass, strjoin(r.errors, ' | '));

            % 2. machine-derived numThreads is a positive integer
            tc.verifyGreaterThanOrEqual(ir.options.numThreads, 1);

            % 3. matches the reference, section by section, numThreads
            %    neutralized (canonicalize is whole-IR).
            ref = buildHL1996IR();
            ref.options.numThreads = 0;
            ir.options.numThreads  = 0;

            ca = gdsge.ir.canonicalize(ir);
            cb = gdsge.ir.canonicalize(ref);
            secs = fieldnames(cb);
            for i = 1:numel(secs)
                tc.verifyTrue(isequaln(ca.(secs{i}), cb.(secs{i})), ...
                    sprintf('section differs from reference: %s', secs{i}));
            end
        end
        function missingBoundIsReported(tc)
            % a policy var with no inbound must raise a located parser error
            gmod = sprintf([ ...
                'parameters beta;\nbeta = 0.95;\n', ...
                'var_shock z;\nshock_num = 1;\nz = 1;\nshock_trans = 1;\n', ...
                'var_state k;\nk = linspace(0,1,3);\n', ...
                'var_policy c;\n', ...                              % no inbound for c
                'model;\n  equations;\n    c;\n  end;\nend;\n', ...
                'simulate;\n  num_periods = 1;\n  num_samples = 1;\n', ...
                '  initial k 0;\n  initial shock 1;\n  var_simu c;\n  k'' = c;\nend;']);
            tc.verifyError(@() gdsge.parser.parseFrontEnd(gmod, 'tiny'), ...
                'gdsge:parser:missingBound');
        end
    end
end
```

(Changes: the `ref.model` blanking lines and the `modelBodyLeftForPhase3` method
are gone — the model section now participates in the section-by-section match.)

- [ ] **Step 6: Run the phase gate and the front-end tests**

```powershell
& 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe' -batch "addpath('src'); r = runtests('tests/HeatonLucas1996/parser'); exit(double(any([r.Failed])))"
```

Expected: PASS — both `tFullIRHL1996` and `tFrontEndHL1996` green, exit code 0.
If `parserReproducesReferenceIR` fails, debug by diffing section `model`
statement-by-statement against `buildHL1996IR` (order: interpCall, es1, es2,
eb1, eb2, b1p, b2p, s2p, budget_1, budget_2, w1_consis, GDSGE_EXPECT_1,
equity_premium).

- [ ] **Step 7: Commit**

```powershell
git add src/+gdsge/+parser/parseFrontEnd.m src/+gdsge/+parser/assemblePartialIR.m tests/HeatonLucas1996/parser/tFrontEndHL1996.m tests/HeatonLucas1996/parser/tFullIRHL1996.m
git commit -m "feat(parser): full IR for HL1996 - model body wired into parseFrontEnd"
```

---

### Task 8: full suite, Contents.m, PROGRESS.md

**Files:**
- Modify: `src/+gdsge/+parser/Contents.m`
- Modify: `PROGRESS.md`

- [ ] **Step 1: Run the complete suite**

```powershell
pwsh -File tests/run.ps1
```

Expected: exit 0, all tests pass (Phase 0/1/2 suites + the new parser tests).
If anything fails, fix before proceeding — do not mark the phase complete with
red tests.

- [ ] **Step 2: Update `src/+gdsge/+parser/Contents.m`**

Replace the full contents with:

```matlab
% +PARSER  GDSGE gmod front-end: text -> validated IR.
%   parseFrontEnd   - top-level: gmod text + model name -> validated full IR
%   preprocess      - strip comments, join continuations, rewrite deprecated, guard macros
%   splitBlocks     - separate declaration lines from named blocks
%   splitStatements - bracket-aware split into logical statements (+ start lines)
%   defaultSetupCode- parser-owned default flag values
%   evalSetup       - eval a setup script in an isolated workspace; capture vars
%   parseVarDecls   - structural parse of the declaration region (no eval)
%   parseDeclarations - eval + assemble params/shocks/states/variables/bounds/interp
%   parseSimulate   - the simulate block -> IR simulate section
%   resolveOptions  - flag workspace -> curated IR options
%   tokenize        - model-body text -> token stream (name/number/op/punct/prime)
%   parseExpr       - precedence-climbing expression parser -> AST nodes
%   parseModel      - model block -> typed statements + equations (reduction hoisting)
%   analyzeModel    - semantic checks: names, square system, interp arity, transRefs
%   assemblePartialIR - combine sections (incl. model) into a validated IR
```

- [ ] **Step 3: Update `PROGRESS.md`**

In the Phases list, change the Phase 3 line from:

```markdown
- ◐ **Phase 3 — Model-expression parser** (tokenizer + recursive-descent AST; reductions, primes, interp; semantic analysis → full IR for HL1996) — NEXT
- ☐ **Phase 4 — MATLAB codegen** (IR → `iter_HL1996.m`/`simulate_HL1996.m`, no v2struct, error reporting)
```

to:

```markdown
- ☑ **Phase 3 — Model-expression parser** (done 2026-06-12)
  - ☑ `tokenize` (positions, primes, reduction braces) + `parseExpr` (full MATLAB
    precedence: `^` left-assoc, unary between power and multiplicative; all 7 node kinds)
  - ☑ `parseModel`: assign / primed assign / `GDSGE_INTERP_VEC'` / all 4 reductions
    (+ `| trans` pipe syntax); inline reductions hoisted as `GDSGE_<KIND>_<n>`
  - ☑ `analyzeModel`: name resolution, square system, interp arity, transRef
  - ☑ Phase gate green: `parseFrontEnd(HL1996.gmod)` `isequalIR` the reference IR
  - ☑ Fixture temp renamed `__ep_expect` → `GDSGE_EXPECT_1` (+ JSON golden regen)
- ☐ **Phase 4 — MATLAB codegen** (IR → `iter_HL1996.m`/`simulate_HL1996.m`, no v2struct, error reporting) — NEXT
```

Add at the top of the Changelog section:

```markdown
- 2026-06-12: **Phase 3 complete.** Model-expression parser: `tokenize` →
  `parseExpr` (precedence climbing, MATLAB semantics) → `parseModel` (statement
  classification, inline-reduction hoisting to `GDSGE_<KIND>_<n>`) →
  `analyzeModel` (names, square system, interp arity, transRef). `parseFrontEnd`
  now emits the complete HL1996 IR, equal to the hand-authored reference
  (`tFullIRHL1996` gate + JSON round-trip). Fixture hoisted-temp renamed to
  `GDSGE_EXPECT_1`. All on branch `phase3-expression-parser`. Next: Phase 4
  (MATLAB codegen).
```

- [ ] **Step 4: Commit**

```powershell
git add src/+gdsge/+parser/Contents.m PROGRESS.md
git commit -m "docs: mark Phase 3 complete in PROGRESS.md"
```

- [ ] **Step 5: Finish the branch**

Use the superpowers:finishing-a-development-branch skill to decide merge/PR/cleanup
(previous phases merged their branch into `main`).

---

## Spec coverage map

- Spec §4 tokenizer → Task 1. §5 grammar/precedence/normalization → Task 2.
- §6 statement layer: classification, hoisting (`GDSGE_<KIND>_<n>`), pipe transRef,
  equations split, error cases → Tasks 3–4; fixture rename → Task 5.
- §7 semantic analysis (4 checks; bounds stay in `assemblePartialIR`) → Task 6.
- §3 module layout + `parseFrontEnd`/`assemblePartialIR` edits → Task 7.
- §8 testing: tTokenize/tParseExpr/tParseModel/tAnalyzeModel/tFullIRHL1996 →
  Tasks 1, 2, 4, 6, 7; full-suite gate → Task 8.
- §2 out-of-scope items produce explicit `Phase 7` errors (primed-call,
  unprimed `GDSGE_INTERP_VEC`) — Tasks 2 and 4.
