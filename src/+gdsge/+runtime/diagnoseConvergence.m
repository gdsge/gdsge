function msg = diagnoseConvergence(sol, lb, ub, f, eqVal, cfg, level)
% DIAGNOSECONVERGENCE  Human-readable diagnostics for a stuck cartesian solve.
%   Pure formatter — returns the report text; the caller decides fprintf/error.
%   level : 'summary' (tier-1 headlines) | 'full' (all signals + point table)
%   cfg   : .tolSol .probSize .solNames (cell, numel == size(sol,1));
%           optional .stateNames/.stateGrids (point table appended iff present).
%   Signals: NaN count; bounds pinned among unconverged (per named component);
%   near-bound among converged (full only); dominant equation residual vs the
%   other equations; worst-K grid points (full only — iter path).
tolSol = cfg.tolSol;
[nSol, nProb] = size(sol);

if isfield(cfg,'solNames') && numel(cfg.solNames) == nSol
    names = cfg.solNames;
else
    names = arrayfun(@(r) sprintf('x%d', r), 1:nSol, 'UniformOutput', false);
end

unc  = (f > tolSol) | isnan(f);          % 1 × nProb
conv = ~unc;
nUnc = nnz(unc);
span = max(1, ub - lb);                  % nSol × nProb; guards tiny ranges

L = {};
L{end+1} = sprintf('%d of %d problems unconverged (TolSol=%g).', nUnc, nProb, tolSol);

% --- NaN ---
nNan = nnz(isnan(f));
if nNan > 0
    L{end+1} = sprintf('  NaN residuals: %d (suspect log/sqrt of a negative or divide-by-zero in the model body).', nNan);
end

% --- bounds pinned among unconverged (per component) ---
atLB = abs(sol - lb) <= 1e-6 .* span;
atUB = abs(ub - sol) <= 1e-6 .* span;
pinLB = sum(atLB(:, unc), 2);
pinUB = sum(atUB(:, unc), 2);
[~, ord] = sort(pinLB + pinUB, 'descend');
shown = false;
for k = 1:min(5, nSol)
    r = ord(k);
    if pinLB(r) == 0 && pinUB(r) == 0; break; end
    if ~shown
        L{end+1} = '  Bounds (unconverged points pinned within 1e-6 of span):';
        shown = true;
    end
    p = {};
    if pinUB(r) > 0; p{end+1} = sprintf('%d at UB', pinUB(r)); end %#ok<AGROW>
    if pinLB(r) > 0; p{end+1} = sprintf('%d at LB', pinLB(r)); end %#ok<AGROW>
    L{end+1} = sprintf('    %s: %s (of %d unconverged) -- bound may be too tight', ...
        names{r}, strjoin(p, ', '), nUnc);
end

% --- dominant equation residual (vs the other equations) ---
if ~isempty(eqVal) && nUnc > 0
    rowMax = max(abs(eqVal(:, unc)), [], 2);   % nEq × 1 (max ignores NaN)
    rowMax(isnan(rowMax)) = 0;
    [sorted, eord] = sort(rowMax, 'descend');
    if numel(sorted) > 1; medOthers = median(sorted(2:end)); else; medOthers = sorted(1); end
    if medOthers == 0; medOthers = eps; end
    L{end+1} = sprintf('  Largest equation residuals over unconverged points (vs median of the others = %.3g):', medOthers);
    for k = 1:min(3, numel(sorted))
        e = eord(k);
        tag = '';
        if k == 1 && sorted(1) > 5 * medOthers
            tag = '  <-- dominant; likely a bug or needs normalization';
        end
        L{end+1} = sprintf('    eq#%d: %.3g (%.1fx others)%s', e, rowMax(e), rowMax(e) / medOthers, tag);
    end
end

if strcmp(level, 'summary')
    msg = strjoin(L, newline);
    return;
end

% --- full-only: near-bound among converged (per component) ---
nearBand = 0.01;
nearLB = sum(((sol - lb) <= nearBand .* span) & conv, 2);
nearUB = sum(((ub - sol) <= nearBand .* span) & conv, 2);
[~, nord] = sort(nearLB + nearUB, 'descend');
shown = false;
for k = 1:min(5, nSol)
    r = nord(k);
    if nearLB(r) == 0 && nearUB(r) == 0; break; end
    if ~shown
        L{end+1} = '  Bounds (converged points within 1% of span):';
        shown = true;
    end
    p = {};
    if nearUB(r) > 0; p{end+1} = sprintf('%d near UB', nearUB(r)); end %#ok<AGROW>
    if nearLB(r) > 0; p{end+1} = sprintf('%d near LB', nearLB(r)); end %#ok<AGROW>
    L{end+1} = sprintf('    %s: %s', names{r}, strjoin(p, ', '));
end

% --- full-only: worst-K grid points (iter path supplies stateGrids) ---
if isfield(cfg,'stateGrids') && ~isempty(cfg.stateGrids)
    L{end+1} = gdsge.runtime.reportUnconverged(unc, f, cfg.probSize, ...
        cfg.stateNames, cfg.stateGrids, 5);
end

msg = strjoin(L, newline);
end
