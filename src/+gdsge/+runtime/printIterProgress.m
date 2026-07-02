function [printed, retryStats] = printIterProgress(iter, metric, maxF, nUnconverged, elapsedSec, printFreq, noPrint, stopFlag, retryStats)
% PRINTITERPROGRESS  One structured progress line, PrintFreq/NoPrint gated.
%   Returns true if it printed (caller resets its elapsed timer on true).
%   Optional retryStats (see retryStatsInit/Add/EndIter): when the checkpoint
%   prints and any covered iteration retried, a second indented "resolve:"
%   line summarizes the window; the returned stats are reset whenever the
%   checkpoint prints, so counts cover exactly the window since the last
%   printed checkpoint.
if nargin < 9; retryStats = []; end
printed = false;
if noPrint; return; end
if mod(iter, printFreq) == 0 || stopFlag
    fprintf('Iter:%d, Metric:%g, maxF:%g, unconverged:%d, elapsed:%.1fs\n', ...
        iter, metric, maxF, nUnconverged, elapsedSec);
    if ~isempty(retryStats)
        if retryStats.itersRetried > 0
            fprintf('  resolve: %d/%d iters retried, %d rounds total, worst iter %d (%d rounds)\n', ...
                retryStats.itersRetried, retryStats.itersSeen, retryStats.totalRounds, ...
                retryStats.worstIter, retryStats.worstRounds);
        end
        retryStats = gdsge.runtime.retryStatsInit();
    end
    printed = true;
end
end
