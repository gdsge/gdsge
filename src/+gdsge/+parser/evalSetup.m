function GDSGE_WS = evalSetup(GDSGE_SCRIPT)
% EVALSETUP  Eval a setup script in this isolated function workspace and capture
%   every resulting variable into a struct. On error, raise
%   gdsge:parser:setupEvalFailed echoing the numbered script lines. Internal
%   locals use a GDSGE_ prefix (a reserved prefix in gmod) and are excluded from
%   the captured workspace.
try
    eval(GDSGE_SCRIPT);
catch GDSGE_ME
    GDSGE_LINES = regexp(GDSGE_SCRIPT, '\n', 'split');
    GDSGE_MSG = sprintf('Error evaluating setup code:\n');
    for GDSGE_I = 1:numel(GDSGE_LINES)
        GDSGE_MSG = [GDSGE_MSG, sprintf('%3d\t%s\n', GDSGE_I, GDSGE_LINES{GDSGE_I})]; %#ok<AGROW>
    end
    GDSGE_MSG = [GDSGE_MSG, sprintf('\n%s', GDSGE_ME.message)];
    error('gdsge:parser:setupEvalFailed', '%s', GDSGE_MSG);
end
GDSGE_VARS = who;
GDSGE_WS = struct();
for GDSGE_I = 1:numel(GDSGE_VARS)
    GDSGE_NAME = GDSGE_VARS{GDSGE_I};
    if strncmp(GDSGE_NAME, 'GDSGE_', 6); continue; end
    GDSGE_WS.(GDSGE_NAME) = eval(GDSGE_NAME);
end
end
