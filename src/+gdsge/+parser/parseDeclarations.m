function out = parseDeclarations(declText)
% PARSEDECLARATIONS  Declaration region -> IR sections (params, shocks, states,
%   variables with flat slot layout, bounds, interp) plus the eval'd flag
%   workspace .ws (for resolveOptions). Evals param/shock/grid setup in an
%   isolated workspace seeded with parser defaults; interp updates and grids are
%   carried as text, never eval'd.
decls = gdsge.parser.parseVarDecls(declText);

% Phase 9a: var_tensor is now a supported MATLAB-side construct (its
% assignments feed inbound bounds + initial interp values, broadcast over the
% state×shock grid). A tensor *used in the model body* still needs the deferred
% C++-pop path — analyzeModel rejects that case (gdsge:parser:varTensorInBodyUnsupported).

% splitStatements stripped the trailing ';' from each statement; re-terminate
% so eval runs cleanly (no echo) and multi-line matrices end unambiguously.
setupBody = strjoin(decls.setupStmts, sprintf(';\n'));
if ~isempty(setupBody); setupBody = [setupBody, ';']; end
script = [gdsge.parser.defaultSetupCode(), sprintf('\n'), setupBody];
ws = gdsge.parser.evalSetup(script);

% adaptive(<factor>) may name a setup variable (e.g. ADAPTIVE_FACTOR);
% resolve every non-numeric factor against the eval'd workspace.
decls.bounds     = resolveAdaptive(decls.bounds, ws);
decls.boundsInit = resolveAdaptive(decls.boundsInit, ws);

% ----- params -------------------------------------------------------------
params = cell(1, numel(decls.paramNames));
for i = 1:numel(decls.paramNames)
    nm = decls.paramNames{i};
    if ~isfield(ws, nm)
        error('gdsge:parser:setupEvalFailed', 'Parameter "%s" not defined', nm);
    end
    params{i} = struct('name', nm, 'value', double(ws.(nm)));
end

% ----- shocks -------------------------------------------------------------
shockVals = struct();
for i = 1:numel(decls.shockNames)
    nm = decls.shockNames{i};
    if ~isfield(ws, nm)
        error('gdsge:parser:setupEvalFailed', 'Shock variable "%s" not defined', nm);
    end
    shockVals.(nm) = double(ws.(nm));
end
shocks = struct('names', {decls.shockNames}, 'count', double(ws.shock_num), ...
    'values', shockVals, 'transitions', struct('shock_trans', double(ws.shock_trans)));

% ----- states -------------------------------------------------------------
grids = struct();
for i = 1:numel(decls.stateNames)
    grids.(decls.stateNames{i}) = gridFor(decls, decls.stateNames{i});
end
states = struct('names', {decls.stateNames}, 'grids', grids);

% ----- variables (+ slots) ------------------------------------------------
% Every declared var_tensor must carry an assignment (the broadcast RHS).
% Aux declared with a definition statement in the var_aux section are *global*
% aux (constructed post-iteration from the whole solved function); they carry
% no GDSGE_AUX slot, so ordinary aux slots stay contiguous over the remainder.
tensorAssign = orderTensorAssign(decls.tensorNames, decls.tensorAssign);
auxGlobal = resolveAuxGlobal(decls);
globalNames = cellfun(@(g) g.name, auxGlobal, 'UniformOutput', false);
ordinaryAux = decls.aux(~cellfun(@(a) ismember(a.name, globalNames), decls.aux));
variables = struct( ...
    'policy', {assignSlots(decls.policy)}, ...
    'aux',    {assignSlots(ordinaryAux)}, ...
    'auxGlobal', {auxGlobal}, ...
    'interp', {decls.interpNames}, ...
    'tensor', {tensorAssign}, ...
    'output', {decls.outputNames}, ...
    'others', {decls.otherNames});

% ----- interp objects -----------------------------------------------------
interp = cell(1, numel(decls.interpNames));
for i = 1:numel(decls.interpNames)
    nm = decls.interpNames{i};
    interp{i} = struct('name', nm, 'args', {decls.stateNames}, ...
        'initialExpr', lookupExpr(decls.interpInitial, nm), ...
        'updateExpr',  lookupExpr(decls.interpUpdate,  nm));
end

out = struct('params', {params}, 'shocks', shocks, 'states', states, ...
    'variables', variables, 'bounds', {decls.bounds}, 'interp', {interp}, ...
    'ws', ws, 'setupText', setupBody, 'setupBlocks', {decls.sections}, ...
    'setupNames', {reshape(setdiff(fieldnames(ws), {'ans'}, 'stable'), 1, [])}, ...
    'boundsInit', {decls.boundsInit}, ...
    'initVariables', struct( ...
        'policyInit', {assignSlots(decls.policyInit)}, ...
        'auxInit',    {assignSlots(decls.auxInit)}));
end

% ===== locals =============================================================
function bounds = resolveAdaptive(bounds, ws)
for i = 1:numel(bounds)
    if isfield(bounds{i}, 'adaptiveFactor') && ~isnumeric(bounds{i}.adaptiveFactor)
        txt = bounds{i}.adaptiveFactor;
        if ~isfield(ws, txt)
            error('gdsge:parser:badBound', ...
                'adaptive factor "%s" is neither a number nor a setup variable', txt);
        end
        bounds{i}.adaptiveFactor = double(ws.(txt));
    end
end
end

function items = assignSlots(rawItems)
items = cell(1, numel(rawItems));
pos = 1;
for i = 1:numel(rawItems)
    len = rawItems{i}.length;
    items{i} = struct('name', rawItems{i}.name, 'length', len, 'slot', [pos, pos+len-1]);
    pos = pos + len;
end
end

function g = gridFor(decls, name)
for i = 1:numel(decls.gridText)
    if strcmp(decls.gridText{i}.name, name); g = decls.gridText{i}.expr; return; end
end
error('gdsge:parser:missingGrid', 'State "%s" has no grid assignment', name);
end

function e = lookupExpr(lst, name)
e = '';
for i = 1:numel(lst)
    if strcmp(lst{i}.name, name); e = lst{i}.expr; return; end
end
end

function items = resolveAuxGlobal(decls)
% Validate captured global var_aux definitions and build the IR items
% {name, op, of, wrt} in definition order. GDSGE_DER(v, s): v must be a
% var_policy, an ordinary var_aux, or an EARLIER-defined global var_aux
% (chaining); s must be a var_state.
defs = decls.auxGlobalDefs;
items = cell(1, numel(defs));
policyNames = cellfun(@(p) p.name, decls.policy, 'UniformOutput', false);
auxNames    = cellfun(@(a) a.name, decls.aux,    'UniformOutput', false);
defNames    = cellfun(@(d) d.name, defs, 'UniformOutput', false);
ordinaryNames = setdiff(auxNames, defNames, 'stable');
earlier = {};
for i = 1:numel(defs)
    d = defs{i};
    if ismember(d.name, cellfun(@(g) g.name, items(1:i-1), 'UniformOutput', false))
        error('gdsge:parser:globalAuxDuplicateDef', ...
            'Global var_aux "%s" is defined more than once', d.name);
    end
    m = regexp(d.expr, '^GDSGE_DER\s*\(\s*(\w+)\s*,\s*(\w+)\s*\)$', 'tokens', 'once');
    if isempty(m)
        error('gdsge:parser:globalAuxBadOperator', ...
            ['Definition of var_aux "%s" must be a global-aux operator, e.g. ' ...
             '"%s = GDSGE_DER(<var>, <state>);" — got "%s = %s;". (An ordinary ' ...
             'var_aux is assigned in the model block instead.)'], ...
            d.name, d.name, d.name, d.expr);
    end
    of = m{1}; wrt = m{2};
    if lengthOf(decls.aux, d.name) ~= 1
        error('gdsge:parser:globalAuxVectorUnsupported', ...
            'Global var_aux "%s" must be declared with length 1 (no [N] size)', d.name);
    end
    srcLen = max(lengthOf(decls.policy, of), lengthOf(decls.aux, of));
    if srcLen > 1
        error('gdsge:parser:globalAuxVectorUnsupported', ...
            'GDSGE_DER source "%s" (for "%s") must be a length-1 variable', of, d.name);
    end
    if ~ismember(wrt, decls.stateNames)
        error('gdsge:parser:globalAuxBadState', ...
            'GDSGE_DER(%s, %s) for "%s": "%s" is not a declared var_state', ...
            of, wrt, d.name, wrt);
    end
    if ~ismember(of, [policyNames, ordinaryNames, earlier])
        error('gdsge:parser:globalAuxUnknownSource', ...
            ['GDSGE_DER(%s, %s) for "%s": "%s" is not a var_policy, an ordinary ' ...
             'var_aux, or an earlier-defined global var_aux'], of, wrt, d.name, of);
    end
    items{i} = struct('name', d.name, 'op', 'der', 'of', of, 'wrt', wrt);
    earlier{end+1} = d.name; %#ok<AGROW>
end
end

function n = lengthOf(sizedItems, name)
% Declared length of NAME in a parseSizedNames list; 0 when absent.
n = 0;
for i = 1:numel(sizedItems)
    if strcmp(sizedItems{i}.name, name); n = sizedItems{i}.length; return; end
end
end

function items = orderTensorAssign(names, assigns)
% Build the IR tensor list as {name, expr} items in declaration order; every
% declared var_tensor must have an assignment.
items = cell(1, numel(names));
for i = 1:numel(names)
    nm = names{i};
    expr = '';
    for j = 1:numel(assigns)
        if strcmp(assigns{j}.name, nm); expr = assigns{j}.expr; break; end
    end
    if isempty(expr)
        error('gdsge:parser:tensorNoAssign', ...
            'var_tensor "%s" has no assignment expression', nm);
    end
    items{i} = struct('name', nm, 'expr', expr);
end
end
