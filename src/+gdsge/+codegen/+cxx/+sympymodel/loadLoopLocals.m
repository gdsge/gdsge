function loadLoopLocals(w, freeSymbols, shockNames, perShockBases)
% LOADLOOPLOCALS  Declare the per-shock C++ locals a shock-loop body references.
%   perShockBases = names that have a per-shock <name>_GDSGE_GRID array (interp
%   targets, primed-assign targets, AND array-policy unknowns). For each free
%   symbol the Python body uses:
%     <base>__prime, base a shock         -> double <sym> = <base>_GRID(GDSGE_iter);
%     <base>__prime, base a per-shock arr -> double <sym> = <base>_GDSGE_GRID[GDSGE_iter-1];
%     <name>, name a per-shock array      -> double <name> = <name>_GDSGE_GRID[GDSGE_iter-1];
%   Plain names that are scalar-scope (current unknowns, scalar assigns,
%   reduction results) are already in C++ scope and need no load.
for i = 1:numel(freeSymbols)
    s = freeSymbols{i};
    if endsWith(s, '__prime')
        base = s(1:end-7);
        if ismember(base, shockNames)
            w.add('double %s = %s_GRID(GDSGE_iter);', s, base);
        elseif ismember(base, perShockBases)
            w.add('double %s = %s_GDSGE_GRID[GDSGE_iter-1];', s, base);
        else
            error('gdsge:codegen:sympyUnknownPrimed', ...
                'primed symbol %s: base %s is neither a shock nor a per-shock array', s, base);
        end
    elseif ismember(s, perShockBases)
        w.add('double %s = %s_GDSGE_GRID[GDSGE_iter-1];', s, s);
    end
end
end
