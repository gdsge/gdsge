function printed = printIterProgress(iter, metric, maxF, nUnconverged, elapsedSec, printFreq, noPrint, stopFlag)
% PRINTITERPROGRESS  One structured progress line, PrintFreq/NoPrint gated.
%   Returns true if it printed (caller resets its elapsed timer on true).
printed = false;
if noPrint; return; end
if mod(iter, printFreq) == 0 || stopFlag
    fprintf('Iter:%d, Metric:%g, maxF:%g, unconverged:%d, elapsed:%.1fs\n', ...
        iter, metric, maxF, nUnconverged, elapsedSec);
    printed = true;
end
end
