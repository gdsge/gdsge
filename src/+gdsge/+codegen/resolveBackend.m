function dec = resolveBackend(ir, pyAvailableFn)
% RESOLVEBACKEND  Decide the C++ Jacobian backend for this codegen run.
%   dec = struct('backend',b,'mode',m,'reason',r) with b in {'autodiff','sympy'},
%   m in {'explicit','env','auto'}. Precedence:
%     1) in-gmod UseAutoDiff (ir.options.jacobianBackend) — authoritative
%     2) GDSGE_BACKEND env (adept|autodiff|sympy|auto)
%     3) auto-detect: Python available => sympy (with codegen fallback), else autodiff
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
if pyAvailableFn()
    dec = struct('backend', 'sympy', 'mode', 'auto', 'reason', 'Python detected');
else
    dec = struct('backend', 'autodiff', 'mode', 'auto', 'reason', 'Python not detected');
end
end
