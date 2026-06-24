function msg = backendMessage(dec)
% BACKENDMESSAGE  One-line user-facing description of a resolveBackend decision.
switch dec.backend
    case 'sympy';    name = 'SymPy analytic Jacobian';
    case 'autodiff'; name = 'adept autodiff';
    otherwise; error('gdsge:codegen:badBackend', 'unknown backend %s', dec.backend);
end
switch dec.mode
    case {'explicit', 'env'}
        msg = sprintf('Backend: %s (%s).', name, dec.reason);
    case 'auto'
        if strcmp(dec.backend, 'sympy')
            msg = sprintf(['Backend: %s — using %s (auto; set UseAutoDiff=1 to ' ...
                'force adept).'], dec.reason, name);
        else
            msg = sprintf('Backend: %s — using %s (auto).', dec.reason, name);
        end
    otherwise
        error('gdsge:codegen:badBackendMode', 'unknown mode %s', dec.mode);
end
end
