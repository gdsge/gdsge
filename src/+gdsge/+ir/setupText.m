function txt = setupText(ir)
% SETUPTEXT  Concatenate ir.setup section bodies into one eval-able string,
%   in source (list) order. Empty bodies are dropped so the result has no
%   blank-line padding. Returns '' for an IR with no setup sections.
txt = '';
if ~isfield(ir, 'setup') || ~iscell(ir.setup); return; end
bodies = {};
for i = 1:numel(ir.setup)
    s = ir.setup{i};
    if isstruct(s) && isfield(s, 'body') && ~isempty(strtrim(s.body))
        bodies{end+1} = strtrim(s.body); %#ok<AGROW>
    end
end
txt = strjoin(bodies, sprintf('\n'));
end
