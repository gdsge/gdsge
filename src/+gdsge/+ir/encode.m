function txt = encode(ir)
% ENCODE  IR struct -> pretty JSON text. Matrices are emitted as arrays-of-rows
%   (orientation-preserving); non-finite scalars are tagged strings.
s = gdsge.ir.schema();
nodes = s.nodes; stmts = s.stmts; eqkinds = s.eqkinds;
try
    canon = gdsge.ir.canonicalize(ir);
catch ME
    if strcmp(ME.identifier, 'gdsge:ir:canonicalize:unknownKind')
        error('gdsge:ir:encode:unknownKind', '%s', ME.message);
    end
    rethrow(ME);
end
jr = jField(canon, s.root);
txt = jsonencode(jr, 'PrettyPrint', true);

    function v = jField(val, spec)
        switch spec.kind
            case 'scalar'
                v = jScalar(val);
            case 'matrix'
                v = wrapMatrix(val);
            case {'text','enum','ref'}
                v = char(val);
            case 'list'
                v = cell(1, numel(val));
                for i = 1:numel(val); v{i} = jField(val{i}, spec.item); end
            case 'map'
                v = struct(); fn = fieldnames(val);
                for i = 1:numel(fn); v.(fn{i}) = jField(val.(fn{i}), spec.value); end
            case 'tagged'
                v = jTagged(val, spec);
            case 'taggedList'
                v = cell(1, numel(val));
                for i = 1:numel(val); v{i} = jTagged(val{i}, spec); end
            case 'struct'
                v = jStruct(val, spec.fields);
            otherwise
                error('gdsge:ir:encode:kind', 'unknown kind %s', spec.kind);
        end
    end

    function v = jStruct(val, fields)
        v = struct(); fn = fieldnames(fields);
        for i = 1:numel(fn)
            nm = fn{i};
            if isfield(val, nm); v.(nm) = jField(val.(nm), fields.(nm)); end
        end
    end

    function v = jTagged(val, spec)
        switch spec.registry
            case 'stmts';   registry = stmts;
            case 'eqkinds'; registry = eqkinds;
            otherwise;      registry = nodes;
        end
        if ~isfield(val, spec.tag)
            error('gdsge:ir:encode:noTag', 'missing %s on a %s record', spec.tag, spec.registry);
        end
        tag = char(val.(spec.tag));
        if ~isfield(registry, tag)
            error('gdsge:ir:encode:unknownKind', 'unknown %s "%s" in %s', spec.tag, tag, spec.registry);
        end
        v = jStruct(val, registry.(tag).fields);
    end
end

function v = jScalar(x)
if ~isfinite(x)
    if isnan(x); v = 'NaN';
    elseif x > 0; v = 'Inf';
    else; v = '-Inf';
    end
else
    v = x;
end
end

function out = wrapMatrix(M)
% Emit as a 1xR cell of row vectors so jsonencode produces [[...],[...]]
% and jsondecode reconstructs the exact 2-D shape.
% Precondition: M has r>0 rows, OR M is 0x0 empty. A 0xN (N>0) matrix is not
% round-trip safe (it decodes as 0x0); no schema matrix field is ever 0-row.
[r, ~] = size(M);
out = cell(1, r);
for i = 1:r; out{i} = M(i, :); end
end
