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
tensorAssign = orderTensorAssign(decls.tensorNames, decls.tensorAssign);
variables = struct( ...
    'policy', {assignSlots(decls.policy)}, ...
    'aux',    {assignSlots(decls.aux)}, ...
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
