function out = expandMacros(rawText, baseDir)
% EXPANDMACROS  Expand gmod preprocessor macros on raw text, BEFORE comment
%   stripping (old-engine order). Returns struct:
%     .text         expanded gmod text
%     .cxxIncludes  cellstr of C++ '#include ...' lines from cinclude directives
%   baseDir is the directory include()/cinclude() file names resolve against
%   (defaults to pwd). (Error-raising on malformed directives is added as each
%   pass is implemented.)
%   Passes to be inserted here, in order: cinclude, include, #define,
%   #strcat_comma, #mat, #foreach, #for, #if.
%   Each pass is a targeted replacement, so text with no macros is returned
%   byte-identical.
if nargin < 2 || isempty(baseDir); baseDir = pwd; end

code = rawText;
cxxIncludes = {};

% (passes are inserted here in subsequent tasks, in the fixed order above)
[code, cxxIncludes] = passCinclude(code, cxxIncludes);
code = passInclude(code, baseDir);
code = passDefine(code);
code = passStrcatComma(code);
code = passMat(code);
code = passForeach(code);
code = passFor(code);
code = passIf(code);

out = struct('text', code, 'cxxIncludes', {cxxIncludes});
end

function code = passInclude(code, baseDir)
% include("f") / include('f')  ->  splice the text of <baseDir>/f in place.
%   Optional trailing ';'. Single or double quotes. Resolves relative names
%   against baseDir (the top-level gmod directory), NOT relative to the
%   including file's own directory.
%   Anchored to line start so it never matches cinclude(...).
%
%   Pattern: group 1 (['"]) captures the opening quote; group 2 ([^'"\n]*)
%   captures the filename (no quote chars, no newline — self-documenting and
%   cannot cross lines); \1 requires the matching closing quote.
%   tok{1}=quote char, tok{2}=filename. Works for both quote styles.
%
%   Note: the while-loop re-scans spliced text, so a self-referential include
%   would loop forever. This is an accepted limitation — no corpus model uses
%   recursive includes, and cycle detection is not implemented.
pat = '(?<=\n|^)[ \t]*include\(\s*([''"])([^''"\n]*)\1\s*\)\s*;?';
while true
    [s, e, tok] = regexp(code, pat, 'start', 'end', 'tokens', 'once');
    if isempty(s); break; end
    name = tok{2};                       % tok{1}=quote char, tok{2}=filename
    fpath = name;
    if ~isAbsolute(fpath); fpath = fullfile(baseDir, name); end
    if ~exist(fpath, 'file')
        error('gdsge:parser:macroIncludeNotFound', ...
            'include file not found: %s', fpath);
    end
    content = fileread(fpath);
    code = [code(1:s-1), content, code(e+1:end)];
end
end

function tf = isAbsolute(p)
tf = ~isempty(regexp(p, '^([A-Za-z]:[\\/]|[\\/]{2}|[\\/])', 'once'));
end

function [code, inc] = passCinclude(code, inc)
% cinclude("f") / cinclude <f>  ->  collect '#include "f"' / '#include <f>'
%   and strip the directive. Top-to-bottom so include order is preserved.
%   Non-capturing lookbehind; the name is tok{1}.
%
%   Accepted limitations (no corpus model uses cinclude; malformed names
%   surface at C++ compile time, not here):
%     - The quoted form ("f") is NOT whitespace-trimmed — leading/trailing
%       spaces in the filename are taken verbatim.
%     - An empty filename (cinclude("") or cinclude <>) is not rejected.
patQ = '(?<=\n|^)[ \t]*cinclude\(\s*"([^"\n]*)"\s*\)\s*;?';
patA = '(?<=\n|^)[ \t]*cinclude\s*<\s*([^>\n]*?)\s*>\s*;?';
while true
    [sq, eq, tq] = regexp(code, patQ, 'start', 'end', 'tokens', 'once');
    [sa, ea, ta] = regexp(code, patA, 'start', 'end', 'tokens', 'once');
    if isempty(sq) && isempty(sa); break; end
    useQ = ~isempty(sq) && (isempty(sa) || sq < sa);
    if useQ
        inc{end+1} = ['#include "' tq{1} '"']; %#ok<AGROW>
        code = [code(1:sq-1), code(eq+1:end)];
    else
        inc{end+1} = ['#include <' ta{1} '>']; %#ok<AGROW>
        code = [code(1:sa-1), code(ea+1:end)];
    end
end
end

function code = passDefine(code)
% #define NAME <rest-of-line>  ->  word-boundary substitution of NAME by value.
%   Multi-token values allowed (trailing % comment stripped). Processed in
%   source order so a later value may reference an earlier macro name.
%   Names must start with a letter; \< and \> treat '_' as a word char so
%   _foo would also match — the '^[A-Za-z]\w*$' check intentionally rejects it.
%
% DIVERGENCE (deliberate improvements over old engine):
%   1. Line-anchored: a commented-out '% #define ...' is NOT processed (the
%      lookbehind requires the line to start with optional whitespace then '#').
%   2. Redefining a name errors immediately (one-at-a-time scan raises
%      macroDefineMalformed on a duplicate) — old engine silently overwrote.
linePat = '(?<=(^|\n))[ \t]*#define\>[^\n]*';
while true
    [s, e] = regexp(code, linePat, 'start', 'end', 'once');
    if isempty(s); break; end
    line = code(s:e);
    rest = strtrim(regexprep(line, '^[ \t]*#define\>', ''));
    pc = strfind(rest, '%'); if ~isempty(pc); rest = strtrim(rest(1:pc(1)-1)); end
    sp = regexp(rest, '\s', 'once');
    if isempty(sp)
        % rest is either empty (no name) or a bare name with no value
        if isempty(rest)
            error('gdsge:parser:macroDefineMalformed', ...
                '#define needs a NAME and a VALUE: "%s"', strtrim(line));
        else
            error('gdsge:parser:macroDefineMalformed', ...
                '#define %s requires a value', rest);
        end
    end
    name  = rest(1:sp-1);
    value = strtrim(rest(sp+1:end));
    if isempty(regexp(name, '^[A-Za-z]\w*$', 'once'))
        error('gdsge:parser:macroDefineMalformed', ...
            'Invalid #define name "%s"', name);
    end
    code = [code(1:s-1), code(e+1:end)];                 % delete the #define line
    % Escape value for use as a literal regexprep replacement string:
    % backslashes first (so we don't double-escape the \ we add for $),
    % then dollars (regexprep interprets $0/$1/... as backreferences).
    value = strrep(value, '\', '\\');
    value = strrep(value, '$', '\$');
    code = regexprep(code, ['\<' name '\>'], value);
end
end

function code = passStrcatComma(code)
% #strcat_comma{args, range}  ->  e.g. {K, 1:3} -> "K(1), K(2), K(3)"
%   (faithful port of the old process_strcat_comma expression).
%   Pattern forbids newlines inside braces (consistent with passMat).
pat = '#strcat_comma\{([^\n{}]*),([^\n{}]*)\}';
while true
    [s, e, tok] = regexp(code, pat, 'start', 'end', 'tokens', 'once');
    if isempty(s); break; end
    args = strtrim(tok{1});
    nums = evalMacroExpr(strtrim(tok{2}));
    nums_str = strsplit(num2str(nums), ' ');
    nums_with_bracket = strcat('(', nums_str, ')');
    % strjoin default delimiter is a space, producing "K(1), K(2), K(3),";
    % strip(...,',') removes the trailing comma left by strcat(...,',').
    rep = strip(strjoin(strcat(args, nums_with_bracket, ',')), ',');
    code = [code(1:s-1), rep, code(e+1:end)];
end
end

function code = passMat(code)
% #mat{expr}  ->  evaluate expr and interpolate the result as text. The result
%   must be a char row (string — the old engine's primary use) or a scalar
%   (number/logical/string); a multi-element numeric/string array is an error.
%
% DIVERGENCE: the pattern '[^\n}]*' captures up to the FIRST '}' on the line
%   (deliberate restriction — narrower than the old engine's greedy-to-last-'}'.
%   A #mat{} expression containing a literal '}' (cell/struct literal) is not
%   supported. No corpus model uses #mat with such a construct.
pat = '#mat\{([^\n}]*)\}';
while true
    [s, e, tok] = regexp(code, pat, 'start', 'end', 'tokens', 'once');
    if isempty(s); break; end
    val = evalMacroExpr(tok{1});
    if ischar(val) && (isrow(val) || isempty(val))
        rep = val;
    elseif isscalar(val)
        rep = char(string(val));
    else
        error('gdsge:parser:macroEvalFailed', ...
            '#mat{%s} must evaluate to a scalar or a string; got a %s of size %s', ...
            tok{1}, class(val), mat2str(size(val)));
    end
    code = [code(1:s-1), rep, code(e+1:end)];
end
end

function out = passForeach(code)
% #foreach id in v1 v2 ... \n body \n #endfor id  ->  body repeated with #id
%   substituted by each value. Recursive: nested #foreach (distinct ids) work.
%   A #id token is substituted only at non-alphanumeric boundaries on BOTH sides
%   (underscore counts as a boundary, so p_#i and #v_future substitute correctly),
%   matching the old engine exactly. An alphanumeric character immediately before
%   '#' blocks substitution (e.g. K#i is left unchanged).
allButWord = '[^a-zA-Z0-9]';
BEFORE = ['((?<=^)|(?<=' allButWord '))'];
AFTER  = ['((?=$)|(?=' allButWord '))'];
code = [newline, code];
rx = '#foreach\s+(?<id>\w+)\s+in\s+([^\n]*)\n(.*?)#endfor\s+\k<id>';
toks   = regexp(code, rx, 'tokens', 'dotall');
starts = regexp(code, rx, 'start',  'dotall');
ends   = regexp(code, rx, 'end',    'dotall');
if isempty(toks)
    if ~isempty(regexp(code, '#foreach\>', 'once'))
        error('gdsge:parser:macroForeachUnterminated', ...
            '#foreach with no matching #endfor');
    end
    out = code(2:end);     % drop the sentinel newline we prepended
    return;
end
out = '';
for i = 1:numel(toks)
    id   = toks{i}{1};
    vals = strsplit(strtrim(toks{i}{2}), ' ');
    body = toks{i}{3};
    sub = '';
    for j = 1:numel(vals)
        % Escape value for literal use in regexprep replacement string, mirroring
        % passDefine: backslash first (avoids double-escaping the \ added for $),
        % then dollar ($0/$1/... are backreferences in regexprep).
        % Note: 'dotall' on the regexp calls above is harmless belt-and-suspenders
        % since the lazy .*? already crosses newlines; the escaping here is what
        % prevents misinterpretation of foreach values (identifiers/numbers) by
        % regexprep even though current values are benign.
        GDSGE_val = strrep(strtrim(vals{j}), '\', '\\');
        GDSGE_val = strrep(GDSGE_val, '$', '\$');
        sub = [sub, newline, ...
            regexprep(body, [BEFORE '#' id AFTER], GDSGE_val)]; %#ok<AGROW>
    end
    if i == 1
        out = [out, code(1:starts(i)-1), newline]; %#ok<AGROW>
    else
        out = [out, code(ends(i-1)+1:starts(i)-1), newline]; %#ok<AGROW>
    end
    out = [out, passForeach(sub), newline]; %#ok<AGROW>   % recurse into the body
    if i == numel(toks)
        out = [out, code(ends(i)+1:end)]; %#ok<AGROW>
    end
end
end

function out = passFor(code)
% #for id=range \n body \n #end  ->  body repeated with #id and #(id+1)
%   substituted. Balanced #for/#end matching supports nesting (recursion) —
%   a cleanup over the old engine which errored on nested #for.
%   #(id+1) is substituted first via strrep (the literal ')' prevents a
%   collision with e.g. #(i+10)); then #id is substituted via regexprep with a
%   TRAILING alphanumeric boundary (?![A-Za-z0-9]) (underscore excluded) so
%   nested prefix-named iterators (outer 'i', inner 'i2') don't collide —
%   '#i2' is skipped because '2' is alphanumeric, while '#i_' is matched because
%   '_' is not in the class. There is deliberately NO leading boundary so
%   'x#i' -> 'x1' matches the old engine (the forExpandsRange test relies on this).
out = '';
pos = 1;
while true
    s = tokenSearch(code, pos, '#for');
    if isempty(s)
        out = [out, code(pos:end)];
        return;
    end
    out = [out, code(pos:s-1)]; %#ok<AGROW>
    [bodyStart, bodyEnd, blockEnd] = matchFor(code, s);
    header = code(s:bodyStart-1);
    body   = code(bodyStart:bodyEnd);
    ht = regexp(header, '#for\s+(\w+)\s*=\s*([^\n]*)', 'tokens', 'once');
    if isempty(ht)
        error('gdsge:parser:macroForMalformed', ...
            'malformed #for header: "%s"', strtrim(header));
    end
    iter = ht{1};
    rangeExpr = strtrim(ht{2});
    pc = strfind(rangeExpr, '%');
    if ~isempty(pc); rangeExpr = strtrim(rangeExpr(1:pc(1)-1)); end
    vals = evalMacroExpr(rangeExpr);
    expanded = '';
    for v = vals(:)'
        b = strrep(body, ['#(' iter '+1)'], int2str(v+1));
        b = regexprep(b, ['#' iter '(?![A-Za-z0-9])'], int2str(v));
        expanded = [expanded, passFor(b)]; %#ok<AGROW>      % recurse for nesting
    end
    out = [out, expanded]; %#ok<AGROW>
    pos = blockEnd + 1;
end
end

function [bodyStart, bodyEnd, blockEnd] = matchFor(code, s)
% s indexes the '#for' that opens a block. Returns the body span and the index
% of the last char of the matching '#end'. Raises macroForUnterminated if none.
nl = regexp(code(s:end), '\n', 'once');
if isempty(nl)
    error('gdsge:parser:macroForUnterminated', 'no body/#end after #for');
end
bodyStart = s + nl;          % first char after the header newline
depth = 1; i = bodyStart; n = numel(code);
while i <= n
    if tokenAt(code, i, '#for'); depth = depth + 1; i = i + 4; continue; end
    if tokenAt(code, i, '#end')
        depth = depth - 1;
        if depth == 0; bodyEnd = i - 1; blockEnd = i + 3; return; end
        i = i + 4; continue;
    end
    i = i + 1;
end
error('gdsge:parser:macroForUnterminated', 'no matching #end for #for');
end

function s = tokenSearch(code, pos, tok)
% First index >= pos where tok appears as a standalone '#'-token.
s = [];
i = pos;
while i <= numel(code)
    if tokenAt(code, i, tok); s = i; return; end
    i = i + 1;
end
end

function tf = tokenAt(code, i, tok)
% True if tok occurs at index i and is not immediately followed by a word char
% (so '#for' does not match inside '#format', '#end' not inside '#endfor').
n = numel(tok);
tf = false;
if i + n - 1 > numel(code); return; end
if ~strcmp(code(i:i+n-1), tok); return; end
if i + n <= numel(code)
    c = code(i+n);
    if ~isempty(regexp(c, '[A-Za-z0-9_]', 'once')); return; end
end
tf = true;
end

function code = passIf(code)
% #if cond \n body \n #endif  ->  keep body iff eval(cond) is truthy, else drop.
%   Parity with old: non-nesting, non-greedy to the first #endif, no #else.
%
% DIVERGENCE (old-engine parity, documented here for maintainers):
%   1. #if is applied LAST, so eval-bearing macros (#define/#mat/#for) nested
%      inside a FALSE '#if 0' block are NOT conditionally suppressed — they
%      already ran in earlier passes. This matches old-engine ordering.
%   2. A dangling '#if' (no matching '#endif') raises macroIfUnterminated; but
%      a dangling '#endif' (e.g. from an unsupported nested '#if') is left for
%      the downstream parser — old-engine parity; nesting is not supported.
if isempty(regexp(code, '#if\>', 'once')); return; end
if isempty(regexp(code, '#endif\>', 'once'))
    error('gdsge:parser:macroIfUnterminated', '#if without #endif');
end
% Note: MATLAB's lazy (.*?) crosses newlines by default, so no 'dotall' option
% is needed (unlike passForeach which sets it as belt-and-suspenders).
% The \> word-boundary on #endif prevents '#endifx' in a body from matching.
pat = '#if\s+([^\n]*)\n(.*?)#endif\>';
while true
    [s, e, tok] = regexp(code, pat, 'start', 'end', 'tokens', 'once');
    if isempty(s); break; end
    keep = evalMacroExpr(strtrim(tok{1}));
    if keep
        rep = tok{2};
    else
        rep = '';
    end
    code = [code(1:s-1), rep, code(e+1:end)];
end
if ~isempty(regexp(code, '#if\>', 'once'))
    error('gdsge:parser:macroIfUnterminated', 'unmatched #if (no matching #endif)');
end
end

function val = evalMacroExpr(GDSGE_EXPR)
% EVALMACROEXPR  Eval a macro expression in this isolated workspace and return
%   its value. Only literals/built-ins and already-substituted #define values
%   are in scope (macros expand before parameter/grid eval). Wraps failures as
%   gdsge:parser:macroEvalFailed.
try
    val = eval(GDSGE_EXPR);
catch GDSGE_ME
    error('gdsge:parser:macroEvalFailed', ...
        'Failed to evaluate macro expression "%s": %s', GDSGE_EXPR, GDSGE_ME.message);
end
end
