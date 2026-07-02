function stats = retryStatsAdd(stats, minorIters)
% RETRYSTATSADD  Add one solve call's retry rounds to the current bucket.
%   Pass [] to start a fresh accumulator. Call once per solveProblems /
%   solveProblemsAsg call (the ASG refinement loop calls several times per
%   major iteration); retryStatsEndIter closes the bucket.
if isempty(stats); stats = gdsge.runtime.retryStatsInit(); end
stats.curRounds = stats.curRounds + minorIters;
end
