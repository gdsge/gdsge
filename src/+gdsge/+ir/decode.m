function ir = decode(txt)
% DECODE  JSON text -> canonical IR struct. Rejects an incompatible MAJOR version.
ir = gdsge.ir.canonicalize(jsondecode(txt));
s = gdsge.ir.schema();
if isfield(ir,'irVersion') && ischar(ir.irVersion)
    want = sscanf(s.irVersion, '%d', 1);
    got  = sscanf(ir.irVersion, '%d', 1);
    if isempty(got) || got ~= want
        error('gdsge:ir:incompatibleVersion', ...
            'IR irVersion "%s" is incompatible (this toolbox supports major %d).', ...
            ir.irVersion, want);
    end
end
end
