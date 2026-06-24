function analyzeModel(ir)
% ANALYZEMODEL  Model-level semantic checks on a fully assembled IR. The IR's
%   *shape* is gdsge.ir.validate's job; this checks what the validator cannot:
%   1. every name leaf resolves (param/state/shock/policy/aux/interp-output/
%      reduction target/earlier local/builtin); primed names are shocks,
%      interp outputs, policy variables, or earlier primed assign targets
%   2. the system is square after expanding primed equations by shock count
%   3. GDSGE_INTERP_VEC has one target per declared var_interp; a named
%      interp call references a declared var_interp
%   4. every reduction transRef names a declared transition matrix
%   Call-node function names (exp, log, ...) are not variable references and
%   are not resolved here.
builtins = {'GDSGE_Iter','TASK','shock'};
baseKnown = builtins;
for i = 1:numel(ir.params); baseKnown{end+1} = ir.params{i}.name; end %#ok<AGROW>
baseKnown = [baseKnown, reshape(ir.states.names, 1, []), reshape(ir.shocks.names, 1, [])];
policyNames = cell(1, numel(ir.variables.policy));
for i = 1:numel(ir.variables.policy); policyNames{i} = ir.variables.policy{i}.name; end
auxNames = cell(1, numel(ir.variables.aux));
for i = 1:numel(ir.variables.aux); auxNames{i} = ir.variables.aux{i}.name; end
baseKnown = [baseKnown, policyNames, auxNames];
basePrimeable = [reshape(ir.shocks.names, 1, []), policyNames];  % + interp outputs below

tensorNames = cellfun(@(t) t.name, ir.variables.tensor, 'UniformOutput', false);
stateAllowed = [reshape(ir.states.names, 1, []), paramNames(ir)];
shockAllowed = [reshape(ir.shocks.names, 1, []), paramNames(ir)];

nUnknown = 0;
for i = 1:numel(ir.variables.policy)
    nUnknown = nUnknown + ir.variables.policy{i}.length;
end

% Each region is an independent square system selected per grid point by its
% state condition; name resolution + square check run per region.
for kR = 1:numel(ir.model.regions)
    R = ir.model.regions{kR};
    validateCondition(R.condition, stateAllowed, sprintf('model region %d', kR));

    % Phase 9a: a tensor referenced in the model body needs the deferred
    % C++-pop path — reject it clearly (before name resolution).
    if ~isempty(tensorNames)
        used = {};
        for i = 1:numel(R.statements); used = [used, bodyNames(R.statements{i})]; end %#ok<AGROW>
        for i = 1:numel(R.equations);  used = [used, eqNames(R.equations{i})];    end %#ok<AGROW>
        hit = intersect(tensorNames, used);
        if ~isempty(hit)
            error('gdsge:parser:varTensorInBodyUnsupported', ...
                ['var_tensor in the model body is not yet supported (deferred ' ...
                 'C++-pop path). Tensor(s) used in body: %s. Tensors may feed ' ...
                 'only inbound bounds and initial interp values.'], strjoin(hit, ', '));
        end
    end

    known = baseKnown; primeable = basePrimeable;
    for i = 1:numel(R.statements)
        s = R.statements{i};
        switch s.type
            case 'assign'
                checkNames(s.expr, known, primeable, ...
                    sprintf('region %d statement %d (assign %s)', kR, i, s.target));
                known{end+1} = s.target; %#ok<AGROW>
                if s.primed; primeable{end+1} = s.target; end %#ok<AGROW>
            case 'reduction'
                if ~(isfield(ir.shocks.transitions, s.transRef) || isSquareParam(ir, s.transRef))
                    error('gdsge:parser:badTransRef', ...
                        'region %d statement %d: reduction transition "%s" is not a declared transition matrix', ...
                        kR, i, s.transRef);
                end
                checkNames(s.body, known, primeable, ...
                    sprintf('region %d statement %d (reduction %s)', kR, i, s.target));
                known{end+1} = s.target; %#ok<AGROW>
            case 'interpCall'
                if strcmp(s.interpRef, 'GDSGE_INTERP_VEC')
                    nDecl = numel(ir.variables.interp);
                    if numel(s.targets) ~= nDecl
                        error('gdsge:parser:interpArity', ...
                            'region %d statement %d: GDSGE_INTERP_VEC has %d targets but %d var_interp are declared', ...
                            kR, i, numel(s.targets), nDecl);
                    end
                elseif ~ismember(s.interpRef, ir.variables.interp)
                    error('gdsge:parser:unknownInterp', ...
                        'region %d statement %d: "%s" is not a declared var_interp', kR, i, s.interpRef);
                else
                    % Named-interp arity cross-check: call arg count must match
                    % the declared interp's state-arg list.
                    j = find(cellfun(@(x) strcmp(x.name, s.interpRef), ir.interp), 1);
                    if ~isempty(j)
                        nExpected = numel(ir.interp{j}.args);
                        nGot      = numel(s.args);
                        if nGot ~= nExpected
                            error('gdsge:parser:interpArity', ...
                                'region %d statement %d: interp "%s" expects %d arg(s) but got %d', ...
                                kR, i, s.interpRef, nExpected, nGot);
                        end
                    end
                end
                for a = 1:numel(s.args)
                    checkNames(s.args{a}, known, primeable, ...
                        sprintf('region %d statement %d (interp call)', kR, i));
                end
                known = [known, reshape(s.targets, 1, [])]; %#ok<AGROW>
                primeable = [primeable, reshape(s.targets, 1, [])]; %#ok<AGROW>
        end
    end

    nEq = 0;
    for e = 1:numel(R.equations)
        eq = R.equations{e};
        if strcmp(eq.kind, 'conditional')
            for c = 1:numel(eq.cases)
                validateCondition(eq.cases{c}.cond, shockAllowed, ...
                    sprintf('region %d equation %d case %d', kR, e, c));
                for jj = 1:numel(eq.cases{c}.equations)
                    checkNames(eq.cases{c}.equations{jj}.expr, known, primeable, ...
                        sprintf('region %d equation %d case %d', kR, e, c));
                end
            end
        else
            checkNames(eq.expr, known, primeable, sprintf('region %d equation %d', kR, e));
        end
        nEq = nEq + slotCount(eq, ir.shocks.count);
    end
    if nUnknown ~= nEq
        error('gdsge:parser:notSquare', ...
            'region %d is not square: %d unknowns vs %d equations (primed = shock_num=%d rows)', ...
            kR, nUnknown, nEq, ir.shocks.count);
    end
end

if isfield(ir, 'modelInit')
    analyzeInit(ir);
end
end

% ---- equation slot counting (Phase 9b) -----------------------------------
function n = slotCount(eq, shockNum)
if strcmp(eq.kind, 'conditional')
    a = caseArity(eq.cases{1}, shockNum);
    for c = 2:numel(eq.cases)
        if caseArity(eq.cases{c}, shockNum) ~= a
            error('gdsge:parser:conditionalArity', ...
                'conditional equation cases have unequal slot counts');
        end
    end
    n = a;
else
    n = 1; if eq.primed; n = shockNum; end
end
end

function a = caseArity(cse, shockNum)
a = 0;
for j = 1:numel(cse.equations)
    a = a + 1;
    if cse.equations{j}.primed; a = a + shockNum - 1; end
end
end

function validateCondition(cond, allowed, where)
cond = strtrim(cond);
if isempty(cond); return; end
ids = regexp(cond, '[A-Za-z_]\w*', 'match');
for i = 1:numel(ids)
    if ~ismember(ids{i}, allowed)
        error('gdsge:parser:badCondition', ...
            '%s: condition "%s" references "%s", not an allowed identifier', where, cond, ids{i});
    end
end
end

function ns = paramNames(ir)
ns = cellfun(@(p) p.name, ir.params, 'UniformOutput', false);
end

function out = eqNames(eq)
% All `name` leaves referenced by a plain or conditional equation entry.
out = {};
if strcmp(eq.kind, 'conditional')
    for c = 1:numel(eq.cases)
        for j = 1:numel(eq.cases{c}.equations)
            out = [out, namesIn(eq.cases{c}.equations{j}.expr)]; %#ok<AGROW>
        end
    end
else
    out = namesIn(eq.expr);
end
end

function analyzeInit(ir)
% Same discipline as the main pass, with the init variable tables. The init
% problem runs before any interpolant exists, so interp calls are an error;
% only shocks are primeable (no policy priming in a static problem).
known = {'GDSGE_Iter','TASK','shock'};
for i = 1:numel(ir.params); known{end+1} = ir.params{i}.name; end %#ok<AGROW>
known = [known, reshape(ir.states.names, 1, []), reshape(ir.shocks.names, 1, [])];
polNames = cellfun(@(p) p.name, ir.modelInit.variables.policyInit, 'UniformOutput', false);
auxNames = cellfun(@(a) a.name, ir.modelInit.variables.auxInit, 'UniformOutput', false);
known = [known, reshape(polNames, 1, []), reshape(auxNames, 1, [])];
primeable = reshape(ir.shocks.names, 1, []);

stmts = ir.modelInit.statements;
for i = 1:numel(stmts)
    s = stmts{i};
    switch s.type
        case 'assign'
            checkNames(s.expr, known, primeable, ...
                sprintf('model_init statement %d (assign %s)', i, s.target));
            known{end+1} = s.target; %#ok<AGROW>
        case 'reduction'
            if ~(isfield(ir.shocks.transitions, s.transRef) || isSquareParam(ir, s.transRef))
                error('gdsge:parser:badTransRef', ...
                    'model_init statement %d: reduction transition "%s" is not declared', ...
                    i, s.transRef);
            end
            checkNames(s.body, known, primeable, ...
                sprintf('model_init statement %d (reduction %s)', i, s.target));
            known{end+1} = s.target; %#ok<AGROW>
        case 'interpCall'
            error('gdsge:parser:interpInInit', ...
                'model_init statement %d: GDSGE_INTERP_VEC is not available in model_init', i);
    end
end
for i = 1:numel(ir.modelInit.equations)
    checkNames(ir.modelInit.equations{i}.expr, known, primeable, ...
        sprintf('model_init equation %d', i));
end
nUnknown = 0;
for i = 1:numel(ir.modelInit.variables.policyInit)
    nUnknown = nUnknown + ir.modelInit.variables.policyInit{i}.length;
end
nEq = 0;
for i = 1:numel(ir.modelInit.equations)
    if ir.modelInit.equations{i}.primed
        nEq = nEq + ir.shocks.count;
    else
        nEq = nEq + 1;
    end
end
if nUnknown ~= nEq
    error('gdsge:parser:notSquare', ...
        'model_init system is not square: %d unknowns vs %d equations', nUnknown, nEq);
end
end

function checkNames(ast, known, primeable, where)
if strcmp(ast.kind, 'name') && ~ismember(ast.id, known)
    error('gdsge:parser:unknownName', '%s: unknown name "%s"', where, ast.id);
end
if strcmp(ast.kind, 'primed') && ~ismember(ast.id, primeable)
    error('gdsge:parser:unknownName', ...
        '%s: "%s''" is not a shock, interp output, policy variable, or primed local', where, ast.id);
end
kids = gdsge.ir.node.children(ast);
for i = 1:numel(kids)
    checkNames(kids{i}, known, primeable, where);
end
end

function out = bodyNames(s)
% All `name` leaves referenced by a model-body statement or equation.
if isfield(s, 'type')
    switch s.type
        case 'assign';    out = namesIn(s.expr);
        case 'reduction'; out = namesIn(s.body);
        case 'interpCall'
            out = {};
            for a = 1:numel(s.args); out = [out, namesIn(s.args{a})]; end %#ok<AGROW>
        otherwise;        out = {};
    end
else
    out = namesIn(s.expr);   % equation item: struct('expr',ast,'primed',..)
end
end

function out = namesIn(ast)
out = {};
if isempty(ast) || ~isstruct(ast) || ~isfield(ast, 'kind'); return; end
if strcmp(ast.kind, 'name'); out = {ast.id}; end
kids = gdsge.ir.node.children(ast);
for i = 1:numel(kids); out = [out, namesIn(kids{i})]; end %#ok<AGROW>
end

function tf = isSquareParam(ir, name)
% A parameter whose value is shock_num x shock_num is a legal reduction
% transition (CaoKS2016: GDSGE_EXPECT{ ... | shock_trans2 }).
tf = false;
n = ir.shocks.count;
for i = 1:numel(ir.params)
    if strcmp(ir.params{i}.name, name) && numel(ir.params{i}.value) == n*n
        tf = true;
        return;
    end
end
end
