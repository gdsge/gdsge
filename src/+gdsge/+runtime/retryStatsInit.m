function stats = retryStatsInit()
% RETRYSTATSINIT  Zeroed retry-stats accumulator.
%   Tracks resolve-retry activity per major iteration between PrintFreq
%   checkpoints; see retryStatsAdd / retryStatsEndIter / printIterProgress.
stats = struct('itersSeen', 0, 'itersRetried', 0, 'totalRounds', 0, ...
    'worstIter', 0, 'worstRounds', 0, 'curRounds', 0);
end
