function printResolveProgress(mode, label, majorIter, minorIter, step, ...
        nUnconverged, worstF, freq, verboseRetry, hitCap)
% PRINTRESOLVEPROGRESS  Throttled resolve-round line + end-of-loop summary.
%   Every line carries the major iteration as a `[iter N]` prefix so the resolve
%   chatter stays anchored to its outer iteration and you can see majors advance.
%   mode = 'round'   : per-round line, printed only when a `freq` window
%                      boundary is crossed this round (step = the increment to
%                      minorIter this round; 1 for the +1 paths, the actual
%                      trials performed for the in-MEX randomize path). Prints at
%                      most once per window even when one step spans several.
%   mode = 'summary' : end-of-loop line, printed ONLY when hitCap is true
%                      ("stopped at MaxMinorIter=N: M still unconverged");
%                      converged solves are summarized at the PrintFreq
%                      checkpoint instead (printIterProgress + retryStats*).
%   label  : '' for cartesian, 'asg ' for ASG (prefixes "resolve ...").
%   Master gate: prints nothing when verboseRetry is false.
if ~verboseRetry
    return;
end
switch mode
    case 'round'
        if floor(minorIter / freq) ~= floor((minorIter - step) / freq)
            fprintf('  [iter %d] %sresolve round %d: %d unconverged, worst residual %g\n', ...
                majorIter, label, minorIter, nUnconverged, worstF);
        end
    case 'summary'
        % Converged retries print nothing here: they are summarized at the
        % PrintFreq checkpoint (printIterProgress + retryStats*). Only the
        % MaxMinorIter-exhaustion line remains a per-solve print.
        if ~hitCap
            return;
        end
        fprintf(['  [iter %d] %sresolve stopped at MaxMinorIter=%d: ' ...
            '%d points still unconverged, worst residual %g\n'], ...
            majorIter, label, minorIter, nUnconverged, worstF);
    otherwise
        error('gdsge:runtime:printResolveProgress:badMode', ...
            'mode must be ''round'' or ''summary'', got ''%s''.', mode);
end
end
