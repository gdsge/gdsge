function txt = emitDeclare(ir)
% EMITDECLARE  DECLARE_CODE: adouble temporaries for the model lambda —
%   scalar temps as adouble x(0); per-shock future arrays stack-allocated
%   (heap vector<adouble> fallback when MAXDIM exceeds MAX_STACK_DIM), with
%   1-based accessor macros. Mirrors the old generator's declare block.
V = gdsge.codegen.cxx.modelVars(ir);
n = ir.shocks.count;
w = gdsge.codegen.codeWriter();
if ~isempty(V.scalarTemps)
    w.add('adouble %s;', strjoin(strcat(V.scalarTemps, '(0)'), ', '));
end
if ~isempty(V.futureVars)
    w.add('#if MAXDIM>MAX_STACK_DIM');
    w.add('vector<adouble> %s;', ...
        strjoin(strcat(V.futureVars, sprintf('_GDSGE_GRID(%d)', n)), ', '));
    w.add('#else');
    w.add('adouble %s;', ...
        strjoin(strcat(V.futureVars, sprintf('_GDSGE_GRID[%d]', n)), ', '));
    w.add('#endif');
end
for i = 1:numel(V.futureVars)
    f = V.futureVars{i};
    w.add('#define %s_GRID(idx) %s_GDSGE_GRID[int(idx)-1]', f, f);
    w.add('#define %s(idx) %s_GDSGE_GRID[int(idx)-1]', f, f);
end
txt = w.str();
end
