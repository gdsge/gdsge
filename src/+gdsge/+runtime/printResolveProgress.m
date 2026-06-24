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
%   mode = 'summary' : one end-of-loop line. hitCap selects the wording:
%                      false -> "converged after N rounds"
%                      true  -> "stopped at MaxMinorIter=N: M still unconverged".
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
        if minorIter <= 0
            return;
        end
        if hitCap
            fprintf(['  [iter %d] %sresolve stopped at MaxMinorIter=%d: ' ...
                '%d points still unconverged, worst residual %g\n'], ...
                majorIter, label, minorIter, nUnconverged, worstF);
        else
            fprintf('  [iter %d] %sresolve converged after %d rounds\n', ...
                majorIter, label, minorIter);
        end
    otherwise
        error('gdsge:runtime:printResolveProgress:badMode', ...
            'mode must be ''round'' or ''summary'', got ''%s''.', mode);
end
end
