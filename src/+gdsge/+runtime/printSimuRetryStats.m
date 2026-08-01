function stats = printSimuRetryStats(stats)
% PRINTSIMURETRYSTATS  One aggregate resolve line at a SimuPrintFreq checkpoint.
%   The simulate counterpart of the "resolve:" line printIterProgress emits for
%   the iteration loop. Models where bulk retrying is normal (TechBubblesV4nblogKRAv3
%   needs ~5300 restarts on its hardest periods) would otherwise report either
%   nothing at all or one diagnostics block per period.
%
%   Prints only when some period in the window retried, and ALWAYS returns a
%   fresh accumulator, so counts cover exactly the window since the last
%   checkpoint. Feed it the stats built by retryStatsInit/Add/EndIter.
if isempty(stats)
    stats = gdsge.runtime.retryStatsInit();
    return;
end
if stats.itersRetried > 0
    fprintf('  resolve: %d/%d periods retried, %d rounds total, worst period %d (%d rounds)\n', ...
        stats.itersRetried, stats.itersSeen, stats.totalRounds, ...
        stats.worstIter, stats.worstRounds);
end
stats = gdsge.runtime.retryStatsInit();
end
