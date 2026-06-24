function unpackOptions(opts, validNames)
% UNPACKOPTIONS  Assign whitelisted option fields into the caller workspace.
%   Replaces v2struct(GDSGE_OPTIONS): explicit, checked, errors on unknown
%   fields instead of silently creating/overwriting workspace variables.
if isempty(opts); return; end
if ~isstruct(opts)
    error('gdsge:options:notAStruct', 'GDSGE_OPTIONS must be a struct.');
end
fn = fieldnames(opts);
unknown = setdiff(fn, validNames);
if ~isempty(unknown)
    error('gdsge:options:unknownField', ...
        'Unknown option field(s): %s\nValid fields: %s', ...
        strjoin(unknown', ', '), strjoin(sort(validNames(:))', ', '));
end
for i = 1:numel(fn)
    assignin('caller', fn{i}, opts.(fn{i}));
end
end
