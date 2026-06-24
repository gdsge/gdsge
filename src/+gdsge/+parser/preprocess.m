function clean = preprocess(rawText)
% PREPROCESS  Clean gmod text into logical lines: strip % comments, join MATLAB
%   line-continuations (...), rewrite deprecated keywords. Macro directives (#)
%   are handled by expandMacros before this stage and will not appear here.
%   NOTE: the comment rule cuts at the first % on a line (sufficient for the
%   core models); string-aware stripping is a Phase-7c hardening point.
rawLines = regexp(rawText, '\r\n|\r|\n', 'split');

% 1. strip % comments
for i = 1:numel(rawLines)
    rawLines{i} = stripComment(rawLines{i});
end

% 2. join line-continuations
logical = {};
acc = '';
for i = 1:numel(rawLines)
    t = strtrim(rawLines{i});
    if endsWith(t, '...')
        acc = [acc, ' ', t(1:end-3)]; %#ok<AGROW>
    else
        logical{end+1} = strtrim([acc, ' ', t]); %#ok<AGROW>
        acc = '';
    end
end
accTrim = strtrim(acc);
if ~isempty(accTrim)
    logical{end+1} = accTrim;
end

% 3. deprecated-keyword rewrite (word-boundary)
clean = strjoin(logical, sprintf('\n'));
clean = regexprep(clean, '\<GNDSGE', 'GDSGE');
end

function s = stripComment(line)
idx = strfind(line, '%');
if isempty(idx); s = line; else; s = line(1:idx(1)-1); end
end
