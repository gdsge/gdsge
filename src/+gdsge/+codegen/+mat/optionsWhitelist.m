function lit = optionsWhitelist(ir, extra)
% OPTIONSWHITELIST  The GDSGE_OPTIONS_VALID cell literal shared by all
%   generated entry points: frozen option surface + model names + gmod
%   setup-assigned names (old semantics: GDSGE_OPTIONS may override anything
%   the declaration-region code defined, e.g. Bianchi's bMin/bMax) + names
%   read by the MATLAB pre/post_iter and pre/post_simulate_iter hooks.
base = {'TolEq','TolSol','TolFun','PrintFreq','ResolvePrintFreq','NoPrint','SaveFreq','NoSave', ...
    'SimuPrintFreq','SimuSaveFreq','MaxIter','MaxMinorIter','num_samples', ...
    'num_periods','SolMaxIter','UseBroyden','FiniteDiffDelta','GDSGE_USE_BROYDEN', ...
    'GDSGE_DEBUG_EVAL_ONLY','INTERP_ORDER','EXTRAP_ORDER','OutputInterpOrder', ...
    'IterSaveAll','SkipModelInit','UseAdaptiveBound','UseAdaptiveBoundInSol', ...
    'EnforceSimuStateInbound','REUSE_WARMUP_SOL','INTERP_WARMUP_SOL', ...
    'CONSTRUCT_OUTPUT','NumThreads','WarmUp','UseMexResolve','DiagnoseMinorIter', ...
    'UseMexRandomize','MexRandomizeBatch','MexRandomSeed'};
modelNames = {};
for i = 1:numel(ir.params); modelNames{end+1} = ir.params{i}.name; end %#ok<AGROW>
modelNames = [modelNames, ir.shocks.names, fieldnames(ir.shocks.transitions)', ir.states.names];
setupNames = {};
if isfield(ir, 'setupNames'); setupNames = reshape(ir.setupNames, 1, []); end
% identifiers read by the gmod's MATLAB hooks: an option may feed one even when
% no setup line defines it (gdsge.com huggett1997: options.T / r_t / w_t read
% in pre_iter) — the old v2struct(GDSGE_OPTIONS) accepted them
hookNames = {};
if isfield(ir, 'hooks')
    for f = {'preIter', 'postIter', 'preSimulateIter', 'postSimulateIter'}
        if isfield(ir.hooks, f{1}) && ~isempty(ir.hooks.(f{1}))
            hookNames = [hookNames, regexp(ir.hooks.(f{1}), '[A-Za-z]\w*', 'match')]; %#ok<AGROW>
        end
    end
end
hookNames = hookNames(~startsWith(hookNames, 'GDSGE_') & ~cellfun(@iskeyword, hookNames) ...
    & ~strcmp(hookNames, 'OUTPUT_CONSTRUCT_CODE'));   % placeholder, not a variable
wl = unique([base, modelNames, setupNames, hookNames, reshape(extra, 1, [])], 'stable');
lit = ['{''' strjoin(wl, ''',''') '''}'];
end
