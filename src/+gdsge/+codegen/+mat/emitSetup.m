function txt = emitSetup(ir, target)
% EMITSETUP  Defaults + gmod setup replay block.
%   target: 'iter' or 'simulate' (simulate takes num_periods/num_samples from
%   ir.simulate; iter keeps the old generator's 1/1000 defaults).
%   The gmod declaration-region code is replayed verbatim from ir.setup
%   (old-generator setParamsCode semantics), preserving original order and
%   any intermediate expressions (e.g. Re_n bounds, Ngrid/a1_lb grid text).
w = gdsge.codegen.codeWriter();
opt = ir.options;

w.add('%% ---- toolbox options (frozen public surface; defaults as in the old generator)');
w.add('TolEq = %s;', mat2str(getOpt(opt, 'tolEq', 1e-6), 17));
w.add('TolSol = 1e-8;');
w.add('TolFun = 1e-8;');
w.add('PrintFreq = %s;', mat2str(getOpt(opt, 'printFreq', 10), 17));
w.add('ResolvePrintFreq = %s;', mat2str(getOpt(opt, 'resolvePrintFreq', 100), 17));
w.add('NoPrint = 0;');
w.add('SaveFreq = %s;', mat2str(getOpt(opt, 'saveFreq', 10), 17));
w.add('NoSave = 0;');
w.add('SimuPrintFreq = %s;', mat2str(getOpt(opt, 'simuPrintFreq', 1000), 17));
w.add('SimuSaveFreq = %s;', mat2str(getOpt(opt, 'simuSaveFreq', inf), 17));
w.add('MaxIter = inf;');
w.add('MaxMinorIter = 200;');        % bounded retries; was inf (could loop forever)
w.add('DiagnoseMinorIter = 10;');    % print a diagnostic heartbeat every N retries while stuck
if strcmp(target, 'simulate')
    w.add('num_samples = %s;', mat2str(ir.simulate.numSamples, 17));
    w.add('num_periods = %s;', mat2str(ir.simulate.numPeriods, 17));
else
    w.add('num_samples = 1;');
    w.add('num_periods = 1000;');
end
w.add('SolMaxIter = 200;');
w.add('UseBroyden = 0;');
w.add('FiniteDiffDelta = 1e-6;');
w.add('GDSGE_USE_BROYDEN = 1;');
w.add('GDSGE_DEBUG_EVAL_ONLY = 0;');
w.add('INTERP_ORDER = %s;', mat2str(getOpt(opt, 'interpOrder', 4), 17));
w.add('EXTRAP_ORDER = %s;', mat2str(getOpt(opt, 'extrapOrder', 2), 17));
w.add('OutputInterpOrder = 2;');
w.add('IterSaveAll = 0;');
w.add('SkipModelInit = 0;');
w.add('GDSGE_EMPTY = [];');
w.add('UseAdaptiveBound = 1;');
w.add('UseAdaptiveBoundInSol = 0;');
w.add('UseMexResolve = 1;');
w.add('UseMexRandomize = 1;');     % restart loop runs in the MEX (C++ RNG); 0 = MATLAB path
w.add('MexRandomizeBatch = 100;'); % minor-iterations per MEX call before returning to MATLAB
w.add('MexRandomSeed = 0;');       % fixed seed; reproducible out of the box
w.add('EnforceSimuStateInbound = 1;');
w.add('REUSE_WARMUP_SOL = 1;');
w.add('INTERP_WARMUP_SOL = 1;');
w.add('CONSTRUCT_OUTPUT = 1;');
if strcmp(getOpt(opt, 'interpMethod', 'spline'), 'asg')
    % ASG defaults from old default_mod.nmod (overridable via the gmod
    % setup replay below, exactly as the old generator). GDSGE_ASG_FIX_GRID
    % is seeded here so gmod setup sections that reference it stay valid;
    % the fixed-grid branch itself is not emitted (spec 7b §2 — no corpus
    % driver exercises it). Spline models keep byte-identical output.
    w.add('AsgMinLevel = 4;');
    w.add('AsgMaxLevel = 10;');
    w.add('AsgThreshold = 1e-2;');
    w.add('AsgOutputMaxLevel = 10;');
    w.add('AsgOutputThreshold = 1e-2;');
    w.add('GDSGE_ASG_FIX_GRID = 0;');
end
w.add('MEX_TASK_INIT = 0;');
w.add('MEX_TASK_INF_HORIZON = 1;');
if getOpt(opt, 'numThreads', 0) == 0
    w.add('NumThreads = feature(''numcores'');');
else
    w.add('NumThreads = %s;', mat2str(opt.numThreads, 17));
end
w.blank();

w.add('%% ---- gmod setup (declaration-region code, grouped by declaration block)');
setupPresent = isfield(ir, 'setup') && iscell(ir.setup) && ~isempty(ir.setup) ...
    && ~isempty(strtrim(gdsge.ir.setupText(ir)));
if ~setupPresent
    % Silently skip only for genuinely setup-free IRs (no params, shocks, or states).
    % If the IR declares params/shocks/states but has no setup text the generated file
    % would silently lack all its definitions — a hard-to-diagnose runtime failure.
    hasDecls = ~isempty(ir.params) || ...
               ~isempty(ir.shocks.names) || ...
               ~isempty(ir.states.names);
    if hasDecls
        error('gdsge:codegen:missingSetup', ...
            ['emitSetup: ir.setup is absent or empty but the IR declares params, ' ...
             'shocks, or states. The generated file would lack all variable ' ...
             'definitions. Populate ir.setup with the gmod declaration-region sections.']);
    end
else
    for i = 1:numel(ir.setup)
        sec = ir.setup{i};
        if ~isstruct(sec) || ~isfield(sec, 'body') || isempty(strtrim(sec.body))
            continue;   % declaration-only block (e.g. inbound-only var_policy)
        end
        w.add('%% -- %s', sec.kind);
        w.addRaw(sec.body);
    end
end
txt = w.str();
end

function v = getOpt(opt, name, dflt)
if isfield(opt, name); v = opt.(name); else; v = dflt; end
end
