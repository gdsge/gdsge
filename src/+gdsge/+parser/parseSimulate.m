function sim = parseSimulate(blockText)
% PARSESIMULATE  The simulate; ... end; block -> IR simulate section. All exprs
%   are carried as opaque text (no AST). The future-marker quote is stripped
%   from transition state names and right-hand sides.
sim = struct('numPeriods', [], 'numSamples', [], ...
    'initial', {{}}, 'varSimu', {{}}, 'transitions', {{}});
stmts = gdsge.parser.splitStatements(blockText);
for i = 1:numel(stmts)
    st = stmts{i};
    if ~isempty(regexp(st, '^num_periods\s*=', 'once'))
        sim.numPeriods = str2double(strtrim(afterEq(st)));
    elseif ~isempty(regexp(st, '^num_samples\s*=', 'once'))
        sim.numSamples = str2double(strtrim(afterEq(st)));
    elseif ~isempty(regexp(st, '^initial\s+', 'once'))
        m = regexp(st, '^initial\s+(\w+)\s+(.+)$', 'tokens', 'once');
        sim.initial{end+1} = struct('var', m{1}, 'value', strtrim(m{2})); %#ok<AGROW>
    elseif ~isempty(regexp(st, '^var_simu\s+', 'once'))
        toks = regexp(strtrim(st(numel('var_simu')+1:end)), '\s+', 'split');
        sim.varSimu = [sim.varSimu, toks(~cellfun(@isempty, toks))];
    else
        m = regexp(st, '^([A-Za-z_]\w*)''?\s*=\s*(.+)$', 'tokens', 'once');
        if ~isempty(m)
            rhs = strtrim(m{2});
            sim.transitions{end+1} = struct('state', m{1}, ...
                'expr', strtrim(regexprep(rhs, '''', '')), ...
                'primed', double(contains(rhs, ''''))); %#ok<AGROW>
        end
    end
end
% Old-toolbox default: a simulate block that omits num_samples runs 2 samples
% (Cao2011EZ relies on this). numPeriods has no default — it must be declared.
if isempty(sim.numSamples) || isnan(sim.numSamples)
    sim.numSamples = 2;
end
end

function r = afterEq(st)
idx = strfind(st, '=');
r = st(idx(1)+1:end);
end
