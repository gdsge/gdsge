function s = indentBy(s, pad)
% INDENTBY  Prepend pad to every non-empty line in s.
%   Shared by emitIter, emitSimulate, emitIterInit.
if isempty(s); return; end
lines = strsplit(s, newline, 'CollapseDelimiters', false);
for i = 1:numel(lines)
    if ~isempty(lines{i}); lines{i} = [pad lines{i}]; end
end
s = strjoin(lines, newline);
end
