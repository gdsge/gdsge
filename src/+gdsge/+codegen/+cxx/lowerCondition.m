function cpp = lowerCondition(cond, ir)
% LOWERCONDITION  Translate a model-region / if-else-equation guard (opaque gmod
%   text on states/shocks/params) to a C++ boolean. Every identifier must be a
%   state, shock, or parameter — each is a C++ local or #define in model scope.
%   '~=' -> '!='; other comparison/logical operators are already valid C++.
%   Empty condition -> always true. Mirrors the old gdsge_parser MODEL_CONDITION
%   / replace_grid_variable handling (the C++ locals are named like the gmod
%   identifiers, so the translation is near-identity).
cond = strtrim(cond);
if isempty(cond); cpp = '1'; return; end
known = [reshape(ir.states.names, 1, []), reshape(ir.shocks.names, 1, [])];
for i = 1:numel(ir.params); known{end+1} = ir.params{i}.name; end %#ok<AGROW>
ids = regexp(cond, '[A-Za-z_]\w*', 'match');
for i = 1:numel(ids)
    if ~ismember(ids{i}, known)
        error('gdsge:codegen:badCondition', ...
            'condition "%s" references unknown identifier "%s" (allow: states, shocks, params)', ...
            cond, ids{i});
    end
end
cpp = strrep(cond, '~=', '!=');
end
