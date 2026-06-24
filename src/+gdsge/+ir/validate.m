function report = validate(ir)
% VALIDATE  Check an IR struct against the schema descriptor.
%   Returns struct with .pass (logical) and .errors (cellstr of located messages).
%   Checks structural shape, kinds, enums, tagged-union dispatch, ref resolution,
%   slot layout, and option invariants. Model-soundness checks that need parse
%   context (square system, bounds completeness) live in the parser, not here.
s = gdsge.ir.schema();
nodes = s.nodes;
stmts = s.stmts;
errors = {};
pools = buildPools(ir);

vField(ir, s.root, 'ir');
checkSlots();
checkOptionInvariants();
checkVersion();

report = struct('pass', isempty(errors), 'errors', {errors});

    function checkSlots()
        checkCat('variables.policy', getCat('policy'));
        checkCat('variables.aux',    getCat('aux'));
    end

    function items = getCat(cat)
        items = {};
        if isfield(ir,'variables') && isfield(ir.variables, cat) && iscell(ir.variables.(cat))
            items = ir.variables.(cat);
        end
    end

    function checkCat(path, items)
        if isempty(items); return; end
        starts = zeros(numel(items),1); stops = zeros(numel(items),1); good = true;
        for i = 1:numel(items)
            it = items{i};
            if ~(isstruct(it) && all(isfield(it, {'slot','length','name'})))
                err(sprintf('%s{%d}', path, i), 'missing name/length/slot'); good = false; continue;
            end
            sl = it.slot;
            if numel(sl) ~= 2
                err(sprintf('%s{%d}.slot', path, i), 'slot must be [start stop]'); good = false; continue;
            end
            starts(i) = sl(1); stops(i) = sl(2);
            if sl(1) > sl(2)
                err(sprintf('%s{%d}.slot', path, i), 'start > stop'); good = false; continue;
            end
            if sl(2) ~= sl(1) + it.length - 1
                err(sprintf('%s{%d}', path, i), ...
                    sprintf('length %g inconsistent with slot span %g', it.length, sl(2)-sl(1)+1));
                good = false;
            end
        end
        if ~good; return; end
        [starts, ord] = sort(starts); stops = stops(ord);
        if starts(1) ~= 1
            err(path, sprintf('slots must start at 1 (got %g)', starts(1)));
        end
        for i = 2:numel(starts)
            if starts(i) ~= stops(i-1) + 1
                err(path, sprintf('slot gap/overlap near index %g', starts(i)));
            end
        end
    end

    function checkOptionInvariants()
        if ~isfield(ir, 'options'); return; end
        o = ir.options;
        if isfield(o, 'interpOrder') && ~ismember(o.interpOrder, [2 4])
            err('options.interpOrder', 'must be 2 or 4');
        end
        if isfield(o, 'interpMethod') && strcmp(o.interpMethod, 'asg')
            if isfield(ir,'variables') && isfield(ir.variables,'tensor') ...
                    && iscell(ir.variables.tensor) && ~isempty(ir.variables.tensor)
                err('variables.tensor', 'must be empty when interpMethod = asg');
            end
            if ~isfield(o, 'asgMaxLevel')
                err('options.asgMaxLevel', 'required when interpMethod = asg');
            end
        end
        if isfield(o, 'jacobianBackend') && strcmp(o.jacobianBackend, 'sympy')
            if isfield(o,'interpMethod') && ~ismember(o.interpMethod, {'spline', 'asg'})
                err('options.jacobianBackend', ...
                    'sympy backend supports interpMethod = spline or asg');
            end
            % Phase 9a: MATLAB-side var_tensor (bounds/initial) never enters the
            % C++ model body the sympy backend differentiates, so it is allowed.
            % A tensor used in the body is rejected earlier by analyzeModel.
            if isfield(ir,'hooks') && isfield(ir.hooks,'cxx') ...
                    && ischar(ir.hooks.cxx) && ~isempty(ir.hooks.cxx)
                err('options.jacobianBackend', ...
                    'sympy backend cannot differentiate cxx blocks');
            end
        end
    end

    function checkVersion()
        if ~isfield(ir,'irVersion') || ~ischar(ir.irVersion); return; end
        want = majorOf(s.irVersion); got = majorOf(ir.irVersion);
        if isnan(got)
            err('irVersion', sprintf('unparseable version "%s"', ir.irVersion));
        elseif got ~= want
            err('irVersion', sprintf('incompatible major version %d (this toolbox supports %d)', got, want));
        end
    end

    function vField(val, spec, path)
        switch spec.kind
            case 'scalar'
                if ~(isnumeric(val) && isscalar(val))
                    err(path, 'expected numeric scalar');
                end
            case 'matrix'
                if ~(isnumeric(val) && ismatrix(val))
                    err(path, 'expected numeric matrix');
                end
            case 'text'
                if ~(ischar(val) || isstring(val))
                    err(path, 'expected text');
                end
            case 'enum'
                if isstring(val) && isscalar(val); val = char(val); end
                if ~ischar(val)
                    err(path, 'expected text (enum)');
                elseif ~ismember(val, spec.values)
                    err(path, sprintf('expected one of {%s}, got "%s"', ...
                        strjoin(spec.values, ','), val));
                end
            case 'ref'
                if ~ischar(val)
                    err(path, 'ref must be text');
                elseif ~isfield(pools, spec.pool)
                    err(path, sprintf('schema bug: unknown pool "%s"', spec.pool));
                elseif ~ismember(val, pools.(spec.pool))
                    err(path, sprintf('unresolved ref "%s" (pool %s)', val, spec.pool));
                end
            case 'list'
                if ~iscell(val)
                    err(path, 'expected list (cell array)');
                else
                    for i = 1:numel(val)
                        vField(val{i}, spec.item, sprintf('%s{%d}', path, i));
                    end
                end
            case 'map'
                if ~isstruct(val) || ~isscalar(val)
                    err(path, 'expected map (scalar struct)');
                else
                    fn = fieldnames(val);
                    for i = 1:numel(fn)
                        vField(val.(fn{i}), spec.value, sprintf('%s.%s', path, fn{i}));
                    end
                end
            case 'tagged'
                vTagged(val, spec, path);
            case 'taggedList'
                if ~iscell(val)
                    err(path, 'expected list (cell array)');
                else
                    for i = 1:numel(val)
                        vTagged(val{i}, spec, sprintf('%s{%d}', path, i));
                    end
                end
            case 'struct'
                vStruct(val, spec, path);
            otherwise
                err(path, sprintf('schema bug: unknown kind "%s"', spec.kind));
        end
    end

    function vStruct(val, spec, path)
        if ~isstruct(val)
            err(path, 'expected struct'); return;
        end
        fn = fieldnames(spec.fields);
        for i = 1:numel(fn)
            nm = fn{i}; fs = spec.fields.(nm); p = sprintf('%s.%s', path, nm);
            if ~isfield(val, nm)
                if fs.required; err(p, 'required field missing'); end
            else
                vField(val.(nm), fs, p);
            end
        end
    end

    function vTagged(val, spec, path)
        switch spec.registry
            case 'nodes'; registry = nodes;
            case 'stmts'; registry = stmts;
            case 'eqkinds'; registry = s.eqkinds;
            otherwise
                err(path, sprintf('schema bug: unknown registry "%s"', spec.registry));
                return;
        end
        if ~isstruct(val) || ~isfield(val, spec.tag)
            err(path, sprintf('expected tagged record with "%s"', spec.tag)); return;
        end
        tag = val.(spec.tag);
        if ~ischar(tag) || ~isfield(registry, tag)
            err(path, sprintf('unknown %s "%s"', spec.tag, toStr(tag))); return;
        end
        vStruct(val, registry.(tag), path);
    end

    function err(path, msg)
        errors{end+1} = sprintf('%s: %s', path, msg); %#ok<AGROW>
    end
end

function s = toStr(x)
if ischar(x); s = x; else; s = class(x); end
end

function pools = buildPools(ir)
pools.transitions = {};
if isfield(ir,'shocks') && isfield(ir.shocks,'transitions') && isstruct(ir.shocks.transitions)
    pools.transitions = fieldnames(ir.shocks.transitions);
end
% CaoKS2016 surface: square shock_num x shock_num parameters are legal
% reduction-pipe targets (the C++ pop defines <name>(idx) for array params).
% Guarded against malformed IR: count must be a real scalar, each param a
% struct with a numeric value (this runs before the shape checks complete).
if isfield(ir,'params') && iscell(ir.params) && isfield(ir,'shocks') ...
        && isfield(ir.shocks,'count') && isnumeric(ir.shocks.count) ...
        && isscalar(ir.shocks.count)
    n = ir.shocks.count;
    for iP = 1:numel(ir.params)
        p = ir.params{iP};
        if isstruct(p) && isfield(p,'name') && isfield(p,'value') ...
                && isnumeric(p.value) && numel(p.value) == n*n
            pools.transitions{end+1} = p.name; %#ok<AGROW>
        end
    end
end
pools.states = {};
if isfield(ir,'states') && isfield(ir.states,'names') && iscell(ir.states.names)
    nm = reshape(ir.states.names, 1, []);
    pools.states = nm(cellfun(@ischar, nm));
end
pools.policyAux = {};
if isfield(ir,'variables')
    pools.policyAux = [listNames(getfieldor(ir.variables,'policy')), ...
                       listNames(getfieldor(ir.variables,'aux'))];
end
pools.interpNames = listNames(getfieldor(ir, 'interp'));
end

function names = listNames(lst)
names = {};
if iscell(lst)
    for i = 1:numel(lst)
        if isstruct(lst{i}) && isfield(lst{i}, 'name'); names{end+1} = lst{i}.name; end %#ok<AGROW>
    end
end
end

function v = getfieldor(s, f)
if isfield(s, f); v = s.(f); else; v = {}; end
end

function m = majorOf(v)
t = regexp(v, '^(\d+)\.', 'tokens', 'once');
if isempty(t); m = NaN; else; m = str2double(t{1}); end
end
