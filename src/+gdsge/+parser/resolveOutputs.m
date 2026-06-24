function variables = resolveOutputs(variables, varSimu, stateNames)
% RESOLVEOUTPUTS  Old-generator output promotion (gdsge_parser.m:619-685):
%   1) error if a declared var_output is a state;
%   2) var_simu names that are neither states nor existing outputs append to
%      variables.output, in var_simu order;
%   3) output names not found in policy or aux append to variables.aux with
%      length 1 (model-body locals the MEX must export), extending the slot
%      layout. Prints the same informational lines as the old parser.
for i = 1:numel(variables.output)
    if ismember(variables.output{i}, stateNames)
        error('gdsge:parser:stateInOutput', ...
            'state %s should not be in var_output', variables.output{i});
    end
end
added = {};
for i = 1:numel(varSimu)
    n = varSimu{i};
    if ismember(n, stateNames) || ismember(n, variables.output); continue; end
    variables.output{end+1} = n;
    added{end+1} = n; %#ok<AGROW>
end
if ~isempty(added)
    fprintf('The following var_simulate are added to var_output: %s\n', strjoin(added, ' '));
end

polNames = cellfun(@(p) p.name, variables.policy, 'UniformOutput', false);
auxNames = cellfun(@(a) a.name, variables.aux, 'UniformOutput', false);
pos = 1;
for i = 1:numel(variables.aux); pos = max(pos, variables.aux{i}.slot(2) + 1); end
added = {};
for i = 1:numel(variables.output)
    n = variables.output{i};
    if ismember(n, polNames) || ismember(n, auxNames); continue; end
    variables.aux{end+1} = struct('name', n, 'length', 1, 'slot', [pos, pos]);
    auxNames{end+1} = n; %#ok<AGROW>
    pos = pos + 1;
    added{end+1} = n; %#ok<AGROW>
end
if ~isempty(added)
    fprintf('The following var_output are added to var_aux: %s\n', strjoin(added, ' '));
end
end
