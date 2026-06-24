function [stmts, startLines] = splitStatements(text)
% SPLITSTATEMENTS  Split MATLAB-ish text into logical statements. Boundaries are
%   ';' or newline at bracket-depth 0; inside ( ) [ ] { } both are ignored, and
%   newlines inside brackets are PRESERVED (they are matrix row separators).
%   Empty statements are dropped. Returns a row cellstr of trimmed statements
%   and, optionally, each statement's 1-based start line within TEXT.
%   Line numbers assume \n-normalized input (as produced by preprocess).
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
