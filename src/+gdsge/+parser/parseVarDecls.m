function decls = parseVarDecls(declText)
% PARSEVARDECLS  Pure structural parse of the declaration region. No eval.
%   Returns a struct of name sets, sized policy/aux items, bounds, interp
%   initial/update expressions, state grid text, and the residual setup
%   statements (param/shock/grid value assignments plus verbatim MATLAB
%   passthrough) to be eval'd later.
stmts = gdsge.parser.splitStatements(declText);

decls = struct('paramNames',{{}}, 'shockNames',{{}}, 'stateNames',{{}}, ...
    'policy',{{}}, 'aux',{{}}, 'interpNames',{{}}, 'outputNames',{{}}, ...
    'tensorNames',{{}}, 'tensorAssign',{{}}, 'otherNames',{{}}, 'bounds',{{}}, ...
    'interpInitial',{{}}, 'interpUpdate',{{}}, 'gridText',{{}}, 'setupStmts',{{}}, ...
    'policyInit',{{}}, 'auxInit',{{}}, 'boundsInit',{{}});

DECL_KINDS = {'parameters','var_shock','var_state','var_policy','var_aux', ...
    'var_interp','var_tensor','var_output','var_others','var_policy_init','var_aux_init'};

assignIdx = [];

% ---- Pass A: keyword declarations ----------------------------------------
for i = 1:numel(stmts)
    st = stmts{i};
    kw = regexp(st, '^([A-Za-z_]\w*)', 'tokens', 'once');
    key = ''; if ~isempty(kw); key = kw{1}; end
    switch key
        case 'parameters'
            decls.paramNames = [decls.paramNames, splitNames(rest(st,key))];
        case 'var_shock'
            decls.shockNames = [decls.shockNames, splitNames(rest(st,key))];
        case 'var_state'
            decls.stateNames = [decls.stateNames, splitNames(rest(st,key))];
        case 'var_policy'
            decls.policy = [decls.policy, parseSizedNames(rest(st,key))];
        case 'var_aux'
            decls.aux = [decls.aux, parseSizedNames(rest(st,key))];
        case 'var_interp'
            decls.interpNames = [decls.interpNames, splitNames(rest(st,key))];
        case 'var_output'
            decls.outputNames = [decls.outputNames, splitNames(rest(st,key))];
        case 'var_tensor'
            decls.tensorNames = [decls.tensorNames, splitNames(rest(st,key))];
        case 'var_others'
            decls.otherNames = [decls.otherNames, splitNames(rest(st,key))];
        case 'inbound'
            decls.bounds{end+1} = parseBound(rest(st,key)); %#ok<AGROW>
        case 'inbound_init'
            decls.boundsInit{end+1} = parseBound(rest(st,key)); %#ok<AGROW>
        case 'var_policy_init'
            decls.policyInit = [decls.policyInit, parseSizedNames(rest(st,key))];
        case 'var_aux_init'
            decls.auxInit = [decls.auxInit, parseSizedNames(rest(st,key))];
        case 'initial'
            decls.interpInitial{end+1} = parseInitial(rest(st,key)); %#ok<AGROW>
        otherwise
            assignIdx(end+1) = i; %#ok<AGROW>
    end
end

% ---- Pass B: classify residual statements ---------------------------------
% Statements that are not keyword declarations. A single `name = rhs`
% assignment is routed by its LHS (interp update / tensor assign / setup +
% state grid). Any other statement form — multi-output [a,b] = f(...), bare
% function calls, struct-field LHS, control-flow fragments — is MATLAB
% passthrough: it joins setupStmts verbatim, runs in the parse-time eval, and
% is replayed in the generated files (old-toolbox behavior). One hard error:
% a command-form statement whose first word is a near-miss of a declaration
% keyword is almost certainly a typo'd declaration; fail early with a
% suggestion instead of obscurely at eval.
for j = 1:numel(assignIdx)
    st = stmts{assignIdx(j)};
    lhs = assignLHS(st);
    if isempty(lhs)
        checkKeywordTypo(st, DECL_KINDS);
        decls.setupStmts{end+1} = st; %#ok<AGROW>
        continue;
    end
    rhs = strtrim(assignRHS(st));
    if ismember(lhs, decls.interpNames)
        decls.interpUpdate{end+1} = struct('name', lhs, 'expr', rhs); %#ok<AGROW>
    elseif ismember(lhs, decls.tensorNames)
        % var_tensor assignment: capture as opaque text (codegen emits it with
        % ndgrid broadcasting). Plain eval would dimension-mismatch.
        decls.tensorAssign{end+1} = struct('name', lhs, 'expr', rhs); %#ok<AGROW>
    else
        decls.setupStmts{end+1} = st; %#ok<AGROW>
        if ismember(lhs, decls.stateNames)
            decls.gridText{end+1} = struct('name', lhs, 'expr', rhs); %#ok<AGROW>
        end
    end
end

% ---- Pass C: positional sectioning ---------------------------------------
% Group setup statements into ordered {kind, body} sections by the preceding
% declaration. Body holds exactly the statements that go to setupStmts (plain
% value assignments), so the concatenation of bodies in source order equals
% the flat setup text emitSetup replays. Interp-updates and tensor-assigns are
% carried separately (NOT in any body), exactly as setupStmts excludes them.
setupSet = containers.Map('KeyType','char','ValueType','logical');
for k = 1:numel(decls.setupStmts); setupSet(decls.setupStmts{k}) = true; end

sections = {};
curKind = 'global';
curBody = {};
headerOpen = false;
for i = 1:numel(stmts)
    st = stmts{i};
    kw = regexp(st, '^([A-Za-z_]\w*)', 'tokens', 'once');
    key = ''; if ~isempty(kw); key = kw{1}; end
    if ismember(key, DECL_KINDS)
        if headerOpen && strcmp(key, curKind)
            continue;   % contiguous same-kind declaration line — extend header
        end
        % flush the current section (skip an empty leading 'global' preamble).
        % A kind may re-open later (old-toolbox layouts interleave declaration
        % blocks); sections are positional and kinds may repeat in ir.setup.
        if ~(strcmp(curKind, 'global') && isempty(curBody))
            sections{end+1} = struct('kind', curKind, 'body', joinBody(curBody)); %#ok<AGROW>
        end
        curBody = {};
        curKind = key; headerOpen = true;
    else
        headerOpen = false;
        expected = expectedSection(st, key, decls);
        if isKey(setupSet, st)
            curBody{end+1} = st; %#ok<AGROW>
        end
        if ~isempty(expected) && ~strcmp(expected, curKind)
            lhs = assignLHS(st);
            warning('gdsge:parser:setupBlockMismatch', ...
                ['"%s" is defined in the "%s" block but is associated with "%s". ' ...
                 'Move this statement so it appears after the "%s" declaration ' ...
                 '(and before the next declaration) to silence this warning.'], ...
                lhs, curKind, expected, expected);
        end
    end
end
if ~(strcmp(curKind, 'global') && isempty(curBody))
    sections{end+1} = struct('kind', curKind, 'body', joinBody(curBody));
end
decls.sections = sections;
end

% Re-terminate setup statements with ';' (splitStatements strips it), matching
% the eval'd setupBody convention so emitSetup's replay does not echo.
function b = joinBody(stmts)
if isempty(stmts); b = ''; return; end
b = [strjoin(stmts, sprintf(';\n')), ';'];
end

% ===== locals =============================================================
function r = rest(st, key)
r = strtrim(st(numel(key)+1:end));
end

function names = splitNames(s)
names = regexp(strtrim(s), '\s+', 'split');
names = names(~cellfun(@isempty, names));
end

function items = parseSizedNames(s)
toks = splitNames(s);
items = cell(1, numel(toks));
for i = 1:numel(toks)
    m = regexp(toks{i}, '^(\w+)(\[(\d+)\])?$', 'tokens', 'once');
    if isempty(m)
        error('gdsge:parser:unknownDeclaration', 'Bad variable token "%s"', toks{i});
    end
    len = 1;
    if numel(m) >= 2 && ~isempty(m{2})
        dig = regexp(m{2}, '\d+', 'match', 'once');
        if ~isempty(dig); len = str2double(dig); end
    end
    items{i} = struct('name', m{1}, 'length', len);
end
end

function b = parseBound(s)
toks = splitNames(s);
if numel(toks) < 3
    error('gdsge:parser:badBound', ...
        'inbound needs "<name> <lower> <upper>", got "%s"', strtrim(s));
end
b = struct('name', toks{1}, 'lower', toks{2}, 'upper', toks{3});
adp = regexp(s, 'adaptive\(\s*([^)]*?)\s*\)', 'tokens', 'once');
if ~isempty(adp)
    v = str2double(adp{1});
    if isnan(v)
        % a setup-variable name (e.g. ADAPTIVE_FACTOR); parseDeclarations
        % resolves it against the eval'd setup workspace
        b.adaptiveFactor = adp{1};
    else
        b.adaptiveFactor = v;
    end
end
end

function o = parseInitial(s)
[nm, rem] = strtok(strtrim(s));
o = struct('name', nm, 'expr', strtrim(rem));
end

function lhs = assignLHS(st)
m = regexp(firstLine(st), '^\s*([A-Za-z_]\w*)\s*(\([^)]*\)|\[[^\]]*\])?\s*=', 'tokens', 'once');
if isempty(m); lhs = ''; else; lhs = m{1}; end
end

function rhs = assignRHS(st)
idx = strfind(st, '=');
rhs = st(idx(1)+1:end);
end

function l = firstLine(st)
parts = regexp(st, '\n', 'split');
l = parts{1};
end

function kind = expectedSection(st, key, decls)
% The section a setup/associate statement *belongs* to, by its LHS. '' = no
% expectation (options, helper intermediates, declaration lines, inbound,
% interp initial/updates).
kind = '';
% 'initial' and interp-name assignments are extracted by keyword/LHS name
% regardless of position (never enter a section body), and the canonical
% layout forces them after model_init -> positionally in var_aux_init.
if ismember(key, {'initial','inbound','inbound_init'}); return; end
lhs = assignLHS(st);
if isempty(lhs); return; end
if any(strcmp(lhs, {'shock_trans','shock_num'}));  kind = 'var_shock'; return; end
if ismember(lhs, decls.shockNames);  kind = 'var_shock';  return; end
if ismember(lhs, decls.stateNames);  kind = 'var_state';  return; end
if ismember(lhs, decls.interpNames); return; end
if ismember(lhs, decls.tensorNames); kind = 'var_tensor'; return; end
if ismember(lhs, decls.paramNames);  kind = 'parameters'; return; end
end

function checkKeywordTypo(st, DECL_KINDS)
% Guard against typo'd declaration keywords. Fires only on command-form
% statements — "<word> <word> ..." with no '=' on the first line, the shape a
% declaration takes — whose first word is within edit distance 2 of a keyword.
% Function-call syntax ("rouwen(...)"), assignments, control-flow headers
% (contain '='), and single-word statements ("end") never match.
l = firstLine(st);
if any(l == '='); return; end
m = regexp(l, '^\s*([A-Za-z_]\w*)\s+[A-Za-z_]', 'tokens', 'once');
if isempty(m); return; end
word = m{1};
KEYWORDS = [DECL_KINDS, {'inbound', 'inbound_init', 'initial'}];
for k = 1:numel(KEYWORDS)
    if levenshtein(word, KEYWORDS{k}) <= 2
        error('gdsge:parser:probableTypo', ...
            ['Unknown keyword "%s" in "%s" — did you mean "%s"? (To pass a ' ...
             'command-form MATLAB statement through to the generated code, ' ...
             'use function syntax: %s(...).)'], ...
            word, strtrim(l), KEYWORDS{k}, word);
    end
end
end

function d = levenshtein(a, b)
% Plain Levenshtein distance. The caller only tests <= 2, so bail out early
% when the length gap alone exceeds that.
la = numel(a); lb = numel(b);
if abs(la - lb) > 2; d = 3; return; end
D = zeros(la + 1, lb + 1);
D(:, 1) = (0:la)';
D(1, :) = 0:lb;
for i = 1:la
    for j = 1:lb
        D(i+1, j+1) = min([D(i, j+1) + 1, D(i+1, j) + 1, D(i, j) + (a(i) ~= b(j))]);
    end
end
d = D(la + 1, lb + 1);
end
