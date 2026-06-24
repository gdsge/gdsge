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
        lexeme = text(i:j-1);
        % Phase 9a: name_GRID(idx) normalizes to name(idx) — the old toolbox's
        % internal per-shock accessor; keep the AST free of the _GRID name so
        % both backends treat qp_GRID(1) exactly like qp(1).
        if numel(lexeme) > 5 && strcmp(lexeme(end-4:end), '_GRID') ...
                && j <= n && text(j) == '('
            lexeme = lexeme(1:end-5);
        end
        toks = push(toks, 'name', lexeme, line, col);
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
