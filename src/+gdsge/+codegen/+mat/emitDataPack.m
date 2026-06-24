function frag = emitDataPack(ir)
% EMITDATAPACK  The GDSGE_DATA packing expressions — THE layout contract with
%   the MEX, derived from gdsge.codegen.dataLayout (shared with the C++ POP
%   sequence). Must match the old generator so generated .m files can drive
%   either MEX.
%   .asgPack    per-problem DATA pack for ASG (evalArrayIdx;evalGrids replace tensor rows)
L = gdsge.codegen.dataLayout(ir);
inner = {}; perProb = {}; simuRows = {};
for i = 1:numel(L.entries)
    en = L.entries(i);
    switch en.role
        case 'shockCount'
            inner{end+1} = 'shock_num'; %#ok<AGROW>
        case {'param','trans','shockVals'}
            inner{end+1} = [en.name '(:)']; %#ok<AGROW>
        case 'shockIdx'
            perProb{end+1} = 'GDSGE_TENSOR_shockIdx(:)'''; %#ok<AGROW>
            simuRows{end+1} = 'shock(:)'''; %#ok<AGROW>
        case 'state'
            perProb{end+1} = ['GDSGE_TENSOR_' en.name '(:)''']; %#ok<AGROW>
            simuRows{end+1} = [en.name '(:)''']; %#ok<AGROW>
    end
end
innerStr = strjoin(inner, ';');
frag = struct();
frag.iterPack = sprintf('GDSGE_DATA(:) = [repmat([%s],1,GDSGE_NPROB);%s];', ...
    innerStr, strjoin(perProb, ';'));
frag.simuData0 = sprintf('GDSGE_data0 = repmat([%s],1,GDSGE_NPROB);', innerStr);
frag.simuPack = sprintf('GDSGE_DATA = [GDSGE_data0;%s];', strjoin(simuRows, ';'));
frag.asgPack = sprintf('GDSGE_DATA(:) = [repmat([%s],1,GDSGE_NPROB);GDSGE_evalArrayIdx;GDSGE_evalGrids];', ...
    innerStr);
frag.maxData = L.nData;
end
