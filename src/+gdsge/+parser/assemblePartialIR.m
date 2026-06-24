function ir = assemblePartialIR(modelName, decl, sim, options, hooks, model, modelInit)
% ASSEMBLEPARTIALIR  Combine front-end sections into a schema-valid IR.
%   MODEL (statements + equations from parseModel) is optional; omitting it
%   yields the Phase-2 partial IR (empty model body). Enforces the
%   bounds-completeness semantic gate, then runs gdsge.ir.validate.
if nargin < 6 || isempty(model)
    model = struct('regions', {{ struct('condition','', ...
        'statements', {{}}, 'equations', {{}}) }});
end
if nargin < 7; modelInit = []; end
ir.irVersion = '2.0.0';
ir.modelName = modelName;
ir.params    = decl.params;
ir.setup     = decl.setupBlocks;
ir.setupNames = decl.setupNames;
ir.options   = options;
ir.shocks    = decl.shocks;
ir.states    = decl.states;
ir.variables = decl.variables;
ir.bounds    = decl.bounds;
ir.interp    = decl.interp;
ir.model     = normalizeModel(model);
ir.simulate  = sim;
ir.hooks     = hooks;
if ~isempty(modelInit)
    modelInit.equations = tagPlain(modelInit.equations);
    ir.modelInit = modelInit;
end

% semantic gate: every policy variable needs an inbound
boundNames = cell(1, numel(decl.bounds));
for i = 1:numel(decl.bounds); boundNames{i} = decl.bounds{i}.name; end
for i = 1:numel(decl.variables.policy)
    pn = decl.variables.policy{i}.name;
    if ~ismember(pn, boundNames)
        error('gdsge:parser:missingBound', 'Policy variable "%s" has no inbound', pn);
    end
end

% semantic gate: every init policy variable needs an inbound_init
if ~isempty(modelInit)
    bn = cellfun(@(b) b.name, modelInit.bounds, 'UniformOutput', false);
    for i = 1:numel(modelInit.variables.policyInit)
        pn = modelInit.variables.policyInit{i}.name;
        if ~ismember(pn, bn)
            error('gdsge:parser:missingBound', ...
                'Init policy variable "%s" has no inbound_init', pn);
        end
    end
end

r = gdsge.ir.validate(ir);
if ~r.pass
    error('gdsge:parser:invalidIR', 'IR failed validation:\n%s', ...
        strjoin(r.errors, sprintf('\n')));
end
end

% ===== model normalization (Phase 9b) =====================================
function m = normalizeModel(model)
% Accept a single-body struct {statements,equations[,condition]} OR an already
% region-shaped struct {regions}. Single body -> one region (condition ''). In
% both cases every equation is tagged 'plain' if it is not already kinded, so
% pre-9b parser output and the 9b region/if-else output both normalize cleanly.
if isfield(model, 'regions')
    regions = model.regions;
    for k = 1:numel(regions)
        regions{k}.equations = tagPlain(regions{k}.equations);
    end
    m = struct('regions', {regions});
else
    cond = '';
    if isfield(model, 'condition'); cond = model.condition; end
    m = struct('regions', {{ struct('condition', cond, ...
        'statements', {model.statements}, ...
        'equations',  {tagPlain(model.equations)}) }});
end
end

function eqs = tagPlain(eqs)
% Idempotent: wrap an un-kinded {expr,primed} equation as a 'plain' record.
for i = 1:numel(eqs)
    if ~isfield(eqs{i}, 'kind')
        eqs{i} = struct('kind','plain','expr', eqs{i}.expr, 'primed', eqs{i}.primed);
    end
end
end
