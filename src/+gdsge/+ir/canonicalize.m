function out = canonicalize(ir)
% CANONICALIZE  Normalize any IR (in-memory or freshly jsondecode'd) to a
%   canonical MATLAB form: object-arrays as row cell arrays of scalar structs,
%   non-finite scalars restored from tag strings, text as char. Idempotent.
s = gdsge.ir.schema();
nodes = s.nodes; stmts = s.stmts;
out = cField(ir, s.root);

    function v = cField(val, spec)
        switch spec.kind
            case 'scalar'
                v = canonScalar(val);
            case 'matrix'
                v = double(val);
            case {'text','enum','ref'}
                v = char(val);
            case 'list'
                c = ensureCell(val); v = cell(1, numel(c));
                for i = 1:numel(c); v{i} = cField(c{i}, spec.item); end
            case 'map'
                v = struct();
                if isstruct(val)
                    fn = fieldnames(val);
                    for i = 1:numel(fn); v.(fn{i}) = cField(val.(fn{i}), spec.value); end
                end
            case 'tagged'
                v = cTagged(val, spec);
            case 'taggedList'
                c = ensureCell(val); v = cell(1, numel(c));
                for i = 1:numel(c); v{i} = cTagged(c{i}, spec); end
            case 'struct'
                v = cStruct(val, spec.fields);
            otherwise
                error('gdsge:ir:canonicalize:kind', 'unknown kind %s', spec.kind);
        end
    end

    function v = cStruct(val, fields)
        v = struct();
        fn = fieldnames(fields);
        for i = 1:numel(fn)
            nm = fn{i};
            if isfield(val, nm)
                v.(nm) = cField(val.(nm), fields.(nm));
            end
        end
    end

    function v = cTagged(val, spec)
        switch spec.registry
            case 'stmts';   registry = stmts;
            case 'eqkinds'; registry = s.eqkinds;
            otherwise;      registry = nodes;
        end
        if ~isfield(val, spec.tag)
            error('gdsge:ir:canonicalize:noTag', 'missing %s on a %s record', spec.tag, spec.registry);
        end
        tag = char(val.(spec.tag));
        if ~isfield(registry, tag)
            error('gdsge:ir:canonicalize:unknownKind', 'unknown %s "%s" in %s', spec.tag, tag, spec.registry);
        end
        v = cStruct(val, registry.(tag).fields);
    end
end

function v = canonScalar(val)
if ischar(val) || isstring(val)
    switch char(val)
        case 'Inf';  v = Inf;
        case '-Inf'; v = -Inf;
        case 'NaN';  v = NaN;
        otherwise;   v = str2double(char(val));
    end
else
    v = double(val);
end
end

function c = ensureCell(val)
if iscell(val)
    c = reshape(val, 1, []);
elseif isstruct(val)
    c = reshape(num2cell(val), 1, []);
elseif isempty(val)
    c = {};
else
    c = {val};
end
end
