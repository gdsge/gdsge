function dec = resolveBackend(ir, pyAvailableFn)
% RESOLVEBACKEND  Decide the C++ Jacobian backend for this codegen run.
%   dec = struct('backend',b,'mode',m,'reason',r) with b in {'autodiff','sympy'},
%   m in {'explicit','env','auto'}. Precedence:
%     1) in-gmod UseAutoDiff (ir.options.jacobianBackend) — authoritative
%     2) GDSGE_BACKEND env (adept|autodiff|sympy|auto)
%     3) auto-detect: Python available => sympy (with codegen fallback), else autodiff
%   Auto mode guards on smoothness: a model that calls max/min/abs, or uses a
%   relational operator (<,>,<=,>=,==,~=), gets adept even with Python
%   present — SymPy Jacobians of non-smooth functions are Heaviside-laden and
%   blow up symbolically (a max-heavy CCP model produced a 1.4 GB C++ file).
%   Explicit and env pins are honored as given.
%   pyAvailableFn (default @gdsge.codegen.sympy.ensurePyenv) is injected for tests.
if nargin < 2; pyAvailableFn = @gdsge.codegen.sympy.ensurePyenv; end

% (1) explicit in-gmod choice
if isfield(ir, 'options') && isfield(ir.options, 'jacobianBackend')
    b = ir.options.jacobianBackend;
    dec = struct('backend', b, 'mode', 'explicit', ...
        'reason', sprintf('UseAutoDiff=%d', double(strcmp(b, 'autodiff'))));
    return;
end

% (2) GDSGE_BACKEND env override
env = lower(strtrim(getenv('GDSGE_BACKEND')));
switch env
    case {'adept', 'autodiff'}
        dec = struct('backend', 'autodiff', 'mode', 'env', 'reason', 'GDSGE_BACKEND=adept');
        return;
    case 'sympy'
        dec = struct('backend', 'sympy', 'mode', 'env', 'reason', 'GDSGE_BACKEND=sympy');
        return;
    case {'', 'auto'}
        % fall through to auto-detection
    otherwise
        error('gdsge:codegen:badBackendEnv', ...
            'GDSGE_BACKEND=%s is not recognized (use adept | sympy | auto).', env);
end

% (3) auto-detection
hit = nonSmoothIn(ir);
if ~isempty(hit)
    dec = struct('backend', 'autodiff', 'mode', 'auto', 'reason', ...
        sprintf('model %s (non-smooth; SymPy Jacobian impractical)', hit));
elseif pyAvailableFn()
    dec = struct('backend', 'sympy', 'mode', 'auto', 'reason', 'Python detected');
else
    dec = struct('backend', 'autodiff', 'mode', 'auto', ...
        'reason', 'SymPy Python unavailable (run gdsge_setup_sympy)');
end
end

function hit = nonSmoothIn(ir)
% Human-readable fragment ('calls max' / 'uses >') for the first non-smooth
% construct anywhere in the model or model_init expressions ('' if none):
% a max/min/abs call node, or a binop node with a relational operator.
% Deep-walks structs/cells so every statement form (assign expr, reduction
% body, interpCall args, plain/conditional equations) is covered without
% enumerating them.
hit = '';
if isfield(ir, 'model');     hit = scan(ir.model);     end
if isempty(hit) && isfield(ir, 'modelInit'); hit = scan(ir.modelInit); end
end

function hit = scan(x)
hit = '';
relOps = {'==', '~=', '<', '>', '<=', '>='};
if iscell(x)
    for i = 1:numel(x)
        hit = scan(x{i});
        if ~isempty(hit); return; end
    end
elseif isstruct(x)
    for j = 1:numel(x)
        s = x(j);
        if isfield(s, 'kind') && strcmp(s.kind, 'call') ...
                && any(strcmp(s.fn, {'max','min','abs'}))
            hit = sprintf('calls %s', s.fn); return
        end
        if isfield(s, 'kind') && strcmp(s.kind, 'binop') ...
                && any(strcmp(s.op, relOps))
            hit = sprintf('uses %s', s.op); return
        end
        names = fieldnames(s);
        for i = 1:numel(names)
            hit = scan(s.(names{i}));
            if ~isempty(hit); return; end
        end
    end
end
end
