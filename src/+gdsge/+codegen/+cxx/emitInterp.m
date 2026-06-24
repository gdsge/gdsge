function frag = emitInterp(ir)
% EMITINTERP  CXX interp sections from the IR interp list (spline + asg).
%   frag.getCode   -> INTERP_GET_CODE  (per-task setup)
%   frag.threadCode-> INTERP_GET_THREAD_CODE (per-thread scratch + evaluators)
%   Spline: GDSGE_SPLINE_VEC ingest + per-interp lambdas (both precisions).
%   ASG:    class-handle evaluator; see emitInterpAsg below.
import gdsge.codegen.cxx.fillTemplate
import gdsge.codegen.cxx.readTemplate

if strcmp(ir.options.interpMethod, 'asg')
    frag = emitInterpAsg(ir);
    return;
end

numInterp = numel(ir.interp);
if numInterp == 0
    frag = struct('getCode', '', 'threadCode', '');
    return;
end
args1 = ir.interp{1}.args;            % all HL1996 interps share {w1}
xdim  = numel(args1);
assert(all(cellfun(@(it) numel(it.args) == xdim, ir.interp)), ...
    'gdsge:codegen:mixedXdim', ...
    'all interp vars must share the same argument dimension (Phase 7b)');
order = ir.options.interpOrder;

% --- construct block (once per task) ---
getCode = fillTemplate(readTemplate('interp_spline_construct.tpl.cpp'), ...
    {'VAR_NUM',    num2str(xdim); ...
     'NUM_INTERP', num2str(numInterp)});

% --- per-interp get block ---
getTpl = readTemplate('interp_spline_get.tpl.cpp');
for j = 1:numInterp
    it   = ir.interp{j};
    a    = it.args;
    adoubleArgs = strjoin(cellfun(@(v) ['adouble ' v], a, 'UniformOutput', false), ', ');
    doubleArgs  = strjoin(cellfun(@(v) ['double '  v], a, 'UniformOutput', false), ', ');
    varNames    = strjoin(a, ', ');
    % compound tokens like GDSGE_PP_INTERP_NAME are invisible to fillTemplate's
    % word boundary (underscore = word char); substitute them literally first.
    filled = strrep(getTpl,  'GDSGE_PP_INTERP_NAME',  ['GDSGE_PP_'  it.name]);
    filled = strrep(filled,  'GDSGE_CPP_INTERP_NAME', ['GDSGE_CPP_' it.name]);
    filled = strrep(filled,  'INTERP_NAME_adouble',   [it.name '_adouble']);
    filled = strrep(filled,  'INTERP_NAME_double',    [it.name '_double']);
    filled = strrep(filled,  'eval_INTERP_ORDER',     ['eval_' num2str(order)]);
    % Now fill the remaining standalone word-boundary holes.
    % INTERP_ORDER still appears standalone in MatlabPp<VAR_NUM,INTERP_ORDER>
    % (fillTemplate matches it there); eval_INTERP_ORDER was already handled.
    filled = fillTemplate(filled, { ...
        'ADOUBLE_VAR_NAME',  adoubleArgs; ...
        'DOUBLE_VAR_NAME',   doubleArgs; ...
        'VAR_NAME',          varNames; ...
        'GDSGE_INTERP_XDIM', num2str(numel(a)); ...
        'INTERP_ORDER',      num2str(order); ...
        'VAR_NUM',           num2str(numel(a))});
    getCode = [getCode newline filled]; %#ok<AGROW>
end

% --- prepare-space block (once per thread) ---
adoubleArgs1 = strjoin(cellfun(@(v) ['adouble ' v], args1, 'UniformOutput', false), ', ');
doubleArgs1  = strjoin(cellfun(@(v) ['double '  v], args1, 'UniformOutput', false), ', ');
varNames1    = strjoin(args1, ', ');

threadCode = fillTemplate(readTemplate('interp_spline_prepare_space.tpl.cpp'), { ...
    'GDSGE_NUM_INTERP',  num2str(numInterp); ...
    'GDSGE_INTERP_XDIM', num2str(xdim); ...
    'ADOUBLE_VAR_NAME',  adoubleArgs1; ...
    'DOUBLE_VAR_NAME',   doubleArgs1; ...
    'VAR_NAME',          varNames1});

frag = struct('getCode', getCode, 'threadCode', threadCode);
end

function frag = emitInterpAsg(ir)
% ASG variant: one shared evaluator over the class handle. getCode (task
% scope) converts GDSGE_ASG_HANDLE; threadCode declares per-thread scratch
% (DIM/NVEC substituted in place; ASG_MAX_LEVEL stays a compile-time define)
% plus the GDSGE_INTERP_VEC lambdas and one named lambda per interp var.
% The evaluator sites are always the full state vector — parseDeclarations
% universally sets interp.args = stateNames regardless of interpMethod,
% so partial-state interp calls cannot be constructed from the front end.
import gdsge.codegen.cxx.readTemplate
states = ir.states.names;
numInterp = numel(ir.interp);
% No interpolants (e.g. the init view): emit nothing — the OLD cpp's
% task_init has an empty interp section (no GDSGE_ASG_HANDLE input, no
% scratch arrays), and the spline branch above guards the same way. Without
% this, ASG_MAX_NVEC fills as 0 and the scratch arrays get a zero size.
if numInterp == 0
    frag = struct('getCode', '', 'threadCode', '');
    return;
end
getCode = readTemplate('interp_asg_construct.tpl.cpp');
threadCode = readTemplate('interp_asg_prepare_space.tpl.cpp');
% Substitute DIM/NVEC in place BEFORE appending the per-interp get-blocks, so
% those blocks keep nothing to substitute but names (old-parser parity).
% ASG_MAX_LEVEL is a distinct literal and survives as a compile-time define.
threadCode = strrep(threadCode, 'ASG_MAX_DIM',  num2str(numel(states)));
threadCode = strrep(threadCode, 'ASG_MAX_NVEC', num2str(numInterp));
getTpl = readTemplate('interp_asg_get.tpl.cpp');
for j = 1:numInterp
    filled = strrep(getTpl, 'INTERP_NAME', ir.interp{j}.name);
    filled = strrep(filled, 'INTERP_IDX', num2str(j - 1));
    threadCode = [threadCode newline filled]; %#ok<AGROW>
end
% State names fill the evaluator signatures/sites. Order matters: ADOUBLE
% before DOUBLE before VAR (they nest as substrings).
adoubleArgs = strjoin(cellfun(@(v) ['adouble ' v], states, 'UniformOutput', false), ',');
doubleArgs  = strjoin(cellfun(@(v) ['double '  v], states, 'UniformOutput', false), ',');
threadCode = strrep(threadCode, 'ADOUBLE_VAR_NAME', adoubleArgs);
threadCode = strrep(threadCode, 'DOUBLE_VAR_NAME',  doubleArgs);
threadCode = strrep(threadCode, 'VAR_NAME', strjoin(states, ','));
frag = struct('getCode', getCode, 'threadCode', threadCode);
end
