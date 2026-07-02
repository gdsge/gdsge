function stats = retryStatsEndIter(stats, iter)
% RETRYSTATSENDITER  Close the current major iteration's retry bucket.
%   Call exactly once per major iteration, after all solves and before the
%   printIterProgress checkpoint, so itersSeen counts iterations (not solve
%   calls) and the N/M checkpoint denominator is right.
if isempty(stats); stats = gdsge.runtime.retryStatsInit(); end
stats.itersSeen = stats.itersSeen + 1;
if stats.curRounds > 0
    stats.itersRetried = stats.itersRetried + 1;
    stats.totalRounds = stats.totalRounds + stats.curRounds;
    if stats.curRounds > stats.worstRounds
        stats.worstRounds = stats.curRounds;
        stats.worstIter = iter;
    end
end
stats.curRounds = 0;
end
