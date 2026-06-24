function [files, finalBackend] = generateCxxWithFallback(ir, outDir, dec, genFn)
% GENERATECXXWITHFALLBACK  Generate the C++ for the resolved backend, falling back to
%   adept autodiff ONLY when an auto-mode sympy run throws during code generation.
%   Explicit/env-pinned backends are honored (the error propagates). genFn is injected
%   for tests (default @gdsge.codegen.generateCxx).
if nargin < 4; genFn = @gdsge.codegen.generateCxx; end
finalBackend = dec.backend;
try
    files = genFn(ir, outDir, struct('backend', dec.backend));
catch err
    if strcmp(dec.mode, 'auto') && strcmp(dec.backend, 'sympy')
        lines = strsplit(err.message, newline);
        fprintf(['Backend: SymPy codegen failed (%s) — falling back to adept ' ...
            'autodiff.\n'], lines{1});
        finalBackend = 'autodiff';
        files = genFn(ir, outDir, struct('backend', 'autodiff'));
    else
        rethrow(err);
    end
end
end
