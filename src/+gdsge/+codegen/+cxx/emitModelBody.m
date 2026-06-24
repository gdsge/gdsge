function txt = emitModelBody(ir)
% EMITMODELBODY  MODEL_CODE: IR model statements -> C++ in the adouble model
%   lambda. Names are emitted NEUTRAL (GDSGE_INTERP_VEC, GDSGE_INTERP_RSLT);
%   emitModel suffixes them _adouble/_double per precision (old-generator
%   trick). Reductions lower to explicit shock loops with stack accumulators
%   (parent spec section 10).
import gdsge.codegen.cxx.emitExpr
V = gdsge.codegen.cxx.modelVars(ir);
n = ir.shocks.count;
numInterp = numel(ir.variables.interp);
transNames = fieldnames(ir.shocks.transitions);
for iP = 1:numel(ir.params)
    if numel(ir.params{iP}.value) == ir.shocks.count^2
        transNames{end+1} = ir.params{iP}.name; %#ok<AGROW>
    end
end
w = gdsge.codegen.codeWriter();
mbody = gdsge.codegen.cxx.modelBody(ir);
for i = 1:numel(mbody.statements)
    st = mbody.statements{i};
    switch st.type
        case 'interpCall'
            if ~st.primed
                error('gdsge:codegen:unsupported', ...
                    'unprimed GDSGE_INTERP_VEC: Phase 7b');
            end
            args = cellfun(@(a) emitExpr(a, V.ctxLoop), st.args, 'UniformOutput', false);
            if strcmp(st.interpRef, 'GDSGE_INTERP_VEC')
                w.add('memset(INTERP_VEC_FLAG,true,sizeof(bool)*%d);', numInterp);
                w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++)', n);
                w.add('{');
                w.add('GDSGE_INTERP_VEC(GDSGE_iter,%s);', strjoin(args, ','));
                for t = 1:numel(st.targets)
                    w.add('%s_GDSGE_GRID[GDSGE_iter-1]=GDSGE_INTERP_RSLT[%d];', ...
                        st.targets{t}, t - 1);
                end
                w.add('}');
            else
                % Named scalar interp call (Barro: log_u1n' = log_u1future'(omega1n')):
                % per-shock loop calling the named evaluator lambda. Emitted
                % NEUTRAL; emitModel suffixes _adouble/_double (interp names
                % are on its suffix list). analyzeModel guarantees interpRef
                % is a declared var_interp.
                if numel(st.targets) ~= 1
                    error('gdsge:codegen:badInterpCall', ...
                        'named interp call %s'' must have exactly 1 target (got %d)', ...
                        st.interpRef, numel(st.targets));
                end
                w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++)', n);
                w.add('{');
                w.add('%s_GDSGE_GRID[GDSGE_iter-1]=%s(GDSGE_iter,%s);', ...
                    st.targets{1}, st.interpRef, strjoin(args, ','));
                w.add('}');
            end
        case 'reduction'
            if ~ismember(st.transRef, transNames)
                error('gdsge:codegen:badTransRef', 'unknown transRef %s', st.transRef);
            end
            body = emitExpr(st.body, V.ctxLoop);
            switch st.kind
                case 'EXPECT', w.add('%s = 0;', st.target);
                case 'MIN',    w.add('%s = 1e20;', st.target);
                case 'MAX',    w.add('%s = -1e20;', st.target);
                case 'PROD',   w.add('%s = 1;', st.target);
                otherwise
                    error('gdsge:codegen:unsupported', 'reduction kind %s', st.kind);
            end
            w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++)', n);
            w.add('{');
            switch st.kind
                case 'EXPECT'
                    w.add('%s = %s+%s((shock)+%d*(GDSGE_iter-1))*(%s);', ...
                        st.target, st.target, st.transRef, n, body);
                case 'MIN'
                    w.add('%s = MIN(%s,%s);', st.target, st.target, body);
                case 'MAX'
                    w.add('%s = MAX(%s,%s);', st.target, st.target, body);
                case 'PROD'
                    w.add('%s = %s*(%s);', st.target, st.target, body);
            end
            w.add('}');
        case 'assign'
            if st.primed
                w.add('for(int GDSGE_iter=1; GDSGE_iter<=%d; GDSGE_iter++)', n);
                w.add('{');
                w.add('%s_GDSGE_GRID[GDSGE_iter-1]=%s;', st.target, ...
                    emitExpr(st.expr, V.ctxLoop));
                w.add('}');
            else
                w.add('%s=%s;', st.target, emitExpr(st.expr, V.ctx));
            end
        otherwise
            error('gdsge:codegen:unsupported', ...
                'model statement type %s: Phase 7b', st.type);
    end
end
txt = w.str();
end
