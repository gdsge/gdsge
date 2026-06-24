function txt = emitSimulateMex(ir)
% EMITSIMULATEMEX  simulate_<model>_mex.cpp from templates/simulate_mex.tpl.cpp.
%   Whole SIMU_INTERP period loop in C++ (stacked uniform-order interp). Only
%   emitted for spline (non-ASG) SIMU_INTERP models whose simulate block is
%   MEX-expressible (var_simu = output vars; transitions = single output-var
%   reference, primed or unprimed) — see gdsge.codegen.cxx.isSimuMexExpressible.
%   Row offsets (each output var's 0-based stacked row) are baked in from the IR
%   slot layout, matching emitResultIter's output_var_index ordering.
import gdsge.codegen.cxx.fillTemplate
import gdsge.codegen.cxx.readTemplate
nStates = numel(ir.states.names);
% The interpolation order is dispatched at runtime in the template from the pp
% (output_interp is built with OutputInterpOrder, default 2), so it is NOT baked.

% 0-based stacked row of each output var (sum of preceding output lengths)
rowOf = containers.Map('KeyType','char','ValueType','double');
row = 0;
for i = 1:numel(ir.variables.output)
    n = ir.variables.output{i};
    rowOf(n) = row;
    row = row + outputLen(ir, n);
end

% state name -> 0-based slot index (SimuRslt field order: states then var_simu)
stateSlot = containers.Map(ir.states.names, num2cell(0:nStates-1));

% RECORD_AND_TRANSITION_CODE ------------------------------------------------
L = {};
for k = 1:numel(ir.simulate.varSimu)
    n = ir.simulate.varSimu{k};
    L{end+1} = sprintf(['varsimu[%d][(size_t)(t-1)*num_samples+i] = ' ...
        'ev.nosearch_eval(xLeft, cell, base + %d);'], k-1, rowOf(n)); %#ok<AGROW>
end
for k = 1:numel(ir.simulate.transitions)
    tr = ir.simulate.transitions{k};
    d  = stateSlot(tr.state);
    n  = strtrim(tr.expr);                        % single output-var name (validated)
    % Next-state column is 0-based index t (= MATLAB column t+1), written for
    % every t in [t0,t1]; state fields are preallocated with num_periods+1 cols.
    if tr.primed
        L{end+1} = sprintf(['state[%d][(size_t)t*num_samples+i] = ' ...
            'ev.nosearch_eval(xLeft, cell, base + %d + (nextShock-1));'], d, rowOf(n)); %#ok<AGROW>
    else
        L{end+1} = sprintf(['state[%d][(size_t)t*num_samples+i] = ' ...
            'ev.nosearch_eval(xLeft, cell, base + %d);'], d, rowOf(n)); %#ok<AGROW>
    end
end
recordCode = strjoin(L, sprintf('\n            '));

% STATE_FIELD_SETUP / VARSIMU_FIELD_SETUP -----------------------------------
sf = {sprintf('double* state[%d];', nStates)};
for j = 1:nStates
    sf{end+1} = sprintf('state[%d] = mxGetPr(mxGetField(simu, 0, "%s"));', j-1, ir.states.names{j}); %#ok<AGROW>
end
stateSetup = strjoin(sf, sprintf('\n    '));

nVarSimu = numel(ir.simulate.varSimu);
vf = {sprintf('double* varsimu[%d];', max(nVarSimu,1))};
for k = 1:nVarSimu
    vf{end+1} = sprintf('varsimu[%d] = mxGetPr(mxGetField(simu, 0, "%s"));', k-1, ir.simulate.varSimu{k}); %#ok<AGROW>
end
varsimuSetup = strjoin(vf, sprintf('\n    '));

% PRINT_FIELDS_CODE — sample-1 value of each field (states then var_simu) ----
allNames = [ir.states.names(:); ir.simulate.varSimu(:)];
pf = {};
for j = 1:numel(allNames)
    pf{end+1} = sprintf(['mexPrintf("%%8.4g", ' ...
        'mxGetPr(mxGetField(simu,0,"%s"))[(size_t)(t-1)*num_samples]);'], allNames{j}); %#ok<AGROW>
end
printCode = strjoin([pf, {'mexPrintf("\n");'}], sprintf('\n            '));

txt = readTemplate('simulate_mex.tpl.cpp');
txt = strrep(txt, 'MODEL_NAME', ir.modelName);
txt = fillTemplate(txt, { ...
    'GDSGE_NSTATES', num2str(nStates); ...
    'STATE_FIELD_SETUP',   stateSetup; ...
    'VARSIMU_FIELD_SETUP', varsimuSetup; ...
    'RECORD_AND_TRANSITION_CODE', recordCode; ...
    'PRINT_FIELDS_CODE', printCode});
end

function L = outputLen(ir, name)
% output var length from the IR (1 for scalars/aux; shock_num for primed/future
% policy vars). Mirrors emitResultIter's lenOf.
L = 1;
for i = 1:numel(ir.variables.policy)
    if strcmp(ir.variables.policy{i}.name, name); L = ir.variables.policy{i}.length; return; end
end
for i = 1:numel(ir.variables.aux)
    if strcmp(ir.variables.aux{i}.name, name); L = ir.variables.aux{i}.length; return; end
end
end
