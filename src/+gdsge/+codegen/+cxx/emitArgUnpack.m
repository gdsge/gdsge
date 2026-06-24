function txt = emitArgUnpack(ir)
% EMITARGUNPACK  ARG_CODE: bind policy unknowns to the solver vector GDSGE_x
%   by IR slot. Scalars are adouble references; array policies (e.g. w1n[8])
%   are pointers with 1-based _GRID accessor macros (old conventions).
w = gdsge.codegen.codeWriter();
for i = 1:numel(ir.variables.policy)
    p = ir.variables.policy{i};
    if p.length == 1
        w.add('adouble& %s = GDSGE_x[%d];', p.name, p.slot(1) - 1);
    else
        w.add('adouble* %s_GDSGE_GRID = GDSGE_x+(%d);', p.name, p.slot(1) - 1);
        w.add('#define %s_GRID(idx) %s_GDSGE_GRID[int(idx)-1]', p.name, p.name);
        w.add('#define %s(idx) %s_GDSGE_GRID[int(idx)-1]', p.name, p.name);
    end
end
txt = w.str();
end
