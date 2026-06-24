function txt = emitAux(ir)
% EMITAUX  HEADER_AUX_ASSIGN_CODE: write aux values out of the model lambda
%   (adept value(); the double lambda's value() is the vendored identity).
w = gdsge.codegen.codeWriter();
for i = 1:numel(ir.variables.aux)
    a = ir.variables.aux{i};
    if a.length == 1
        w.add('GDSGE_aux[%d]=value(%s);', a.slot(1) - 1, a.name);
    else
        w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++)', a.length);
        w.add('{');
        w.add('GDSGE_aux[%d+GDSGE_iter]=value(%s_GDSGE_GRID[GDSGE_iter-1]);', ...
            a.slot(1) - 2, a.name);
        w.add('}');
    end
end
txt = w.str();
end
