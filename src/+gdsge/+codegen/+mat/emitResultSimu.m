function frag = emitResultSimu(ir)
% EMITRESULTSIMU  SimuRslt fragments: prealloc, initial conditions, the
%   GDSGE_OPTIONS.init overwrite, and per-period assignment + state transition.
w1 = gdsge.codegen.codeWriter();   % prealloc (states then varSimu)
for i = 1:numel(ir.states.names)
    w1.add('SimuRslt.%s=zeros(num_samples,num_periods);', ir.states.names{i});
end
for i = 1:numel(ir.simulate.varSimu)
    w1.add('SimuRslt.%s=zeros(num_samples,num_periods);', ir.simulate.varSimu{i});
end

w2 = gdsge.codegen.codeWriter();   % initial conditions
for i = 1:numel(ir.simulate.initial)
    it = ir.simulate.initial{i};
    w2.add('SimuRslt.%s(:,1)=%s;', it.var, it.value);
end

w3 = gdsge.codegen.codeWriter();   % GDSGE_OPTIONS.init overwrite
for i = 1:numel(ir.simulate.initial)
    n = ir.simulate.initial{i}.var;
    w3.add('SimuRslt.%s(:,1:size(GDSGE_OPTIONS.init.%s,2))=GDSGE_OPTIONS.init.%s;', n, n, n);
end

w4 = gdsge.codegen.codeWriter();   % per-period assigns + transitions
for i = 1:numel(ir.simulate.varSimu)
    n = ir.simulate.varSimu{i};
    w4.add('SimuRslt.%s(:,GDSGE_t)=%s;', n, n);
end
for i = 1:numel(ir.simulate.transitions)
    t = ir.simulate.transitions{i};
    if t.primed
        % Primed RHS: policy is shock×state array → index next-period shock
        % (parity: old generator emits GDSGE_SHOCK_VAR_LINEAR_INDEX, line ~330)
        w4.add('SimuRslt.%s(:,GDSGE_t+1) = %s(GDSGE_SHOCK_VAR_LINEAR_INDEX);', t.state, t.expr);
    else
        % Unprimed RHS: already a per-sample scalar/vector → broadcast with (:)
        % Known simplification vs old generator: the old code had a THIRD branch
        % (gdsge_parser.m:599-605) — if the unprimed RHS already contained '('
        % it was emitted as-is (no (:) appended). No Phase-7a model uses an
        % unprimed parenthesized transition RHS, so the binary branch suffices.
        % Revisit if a future model hits that case.
        w4.add('SimuRslt.%s(:,GDSGE_t+1) = %s(:);', t.state, t.expr);
    end
end
frag = struct('prealloc', w1.str(), 'init', w2.str(), ...
    'initOverwrite', w3.str(), 'assign', w4.str());
end
