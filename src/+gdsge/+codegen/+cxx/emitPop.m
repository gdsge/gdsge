function txt = emitPop(ir)
% EMITPOP  POP_CODE: sequential POPN/POPNARRAY unpacking of one GDSGE_DATA
%   column into stack locals, derived from gdsge.codegen.dataLayout — the
%   same descriptor that drives the MATLAB-side packer (emitDataPack), so
%   the two sides cannot drift.
L = gdsge.codegen.dataLayout(ir);
w = gdsge.codegen.codeWriter();
for i = 1:numel(L.entries)
    en = L.entries(i);
    switch en.role
        case 'shockCount'
            w.add('POPN(GDSGE_NUM_SHOCKS_D);');
            w.add('int GDSGE_NUM_SHOCKS = (int) GDSGE_NUM_SHOCKS_D;');
        case 'param'
            if en.length == 1
                w.add('POPN(%s);', en.name);
            else
                w.add('POPNARRAY(%s_GDSGE_GRID,%d);', en.name, en.length);
                w.add('#define %s(idx) %s_GDSGE_GRID[int(idx)-1]', en.name, en.name);
            end
        case 'trans'
            if strcmp(en.name, 'shock_trans')
                w.add('POPNARRAY(GDSGE_shock_trans,%d);', en.length);
                w.add('#define shock_trans(idx) GDSGE_shock_trans[int(idx)-1]');
            else
                w.add('POPNARRAY(%s_GDSGE_GRID,%d);', en.name, en.length);
                w.add('#define %s(idx) %s_GDSGE_GRID[int(idx)-1]', en.name, en.name);
            end
        case 'shockVals'
            w.add('POPNARRAY(%s_GDSGE_GRID,%d);', en.name, en.length);
            w.add('#define %s_GRID(idx) %s_GDSGE_GRID[int(idx)-1]', en.name, en.name);
        case 'shockIdx'
            w.add('POPN(GDSGE_shockIdx_d);');
            w.add('int shock = (int)GDSGE_shockIdx_d;');
            for j = 1:numel(ir.shocks.names)
                w.add('double %s = %s_GRID(shock);', ir.shocks.names{j}, ir.shocks.names{j});
            end
        case 'state'
            w.add('POPN(%s);', en.name);
    end
end
txt = w.str();
end
