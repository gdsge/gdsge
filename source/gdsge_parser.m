% Copyright (C) 2020 GDSGE Dan Cao, Wenlan Luo and Guangyu Nie
%
% gdsge.cln2020@gmail.com 
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
% 
% The license allows you to use, compose and distribute the GDSGE compiler 
% or generated codes freely. However, it is requested that the companion paper 
% be cited:
% 
% **Cao, Dan and Luo, Wenlan and Nie, Guangyu, Global DSGE Models (April 1, 2020). 
% Available at SSRN: https://ssrn.com/abstract=3569013 
% or http://dx.doi.org/10.2139/ssrn.3569013**

function [model,matlabCode,cppCode,simulateCode,compileCode,codeSegment] = gdsge_parser(modelName)
% This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE

try
    code = fileread([modelName '.gmod']);
catch ME
    code = fileread([modelName '.nmod']);
end

parser_folder = mfilename('fullpath');
gdsge_folder = fileparts(parser_folder);
template_folder = [gdsge_folder '/code_template'];

LINE_BREAK = newline;
BEFORE_A_WORD = '((?<=^)|(?<=\W))';
AFTER_A_WORD = '((?=$)|(?=\W))';

defaultMod = fileread([gdsge_folder '/code_template/default_mod.nmod']);
code = [defaultMod,LINE_BREAK,code];

% Process GNDSGE to GDSGE
code = regexprep(code,'GNDSGE','GDSGE');
code = process_deprecate(code);

%% Process macro: cinclude
cincludeList = regexp(code,'(?<=(\n|^)[\s]*)cinclude\(');
includeString = cell(1,length(cincludeList));
cxxIncludeString = '';
for i=1:length(cincludeList)
    % Find the closing bracket
    locStart = cincludeList(i);
    locEnd = locStart;
    while code(locEnd)~=')'
        locEnd = locEnd+1;
    end
    if code(locEnd+1)==';'
        locEnd = locEnd+1;
    end
    includeString{i} = code(locStart:locEnd);
    % Get the file inside the bracket
    includeName = regexp(includeString{i},'(?<='').*(?='')','match');
    % Construct include strings in cxx
    if ~isempty(includeName)
        cxxIncludeString = [cxxIncludeString, ...
            '#include "', includeName{1}, '"', LINE_BREAK, ...
            ];
    end
end
for i=1:length(cincludeList)
    code = strrep(code,includeString{i},'');
end

%% Process macro: cinclude<>
cincludeList = regexp(code,'(?<=(\n|^)[\s]*)cinclude \<');
includeString = cell(1,length(cincludeList));
for i=1:length(cincludeList)
    % Find the closing bracket
    locStart = cincludeList(i);
    locEnd = locStart;
    while code(locEnd)~='>'
        locEnd = locEnd+1;
    end
    if code(locEnd+1)==';'
        locEnd = locEnd+1;
    end
    includeString{i} = code(locStart:locEnd);
    % Get the file inside the bracket
    includeName = regexp(includeString{i},'(?<=<).*(?=>)','match');
    % Construct include strings in cxx
    if ~isempty(includeName)
        cxxIncludeString = [cxxIncludeString, ...
            '#include <', includeName{1}, '>', LINE_BREAK, ...
            ];
    end
end
for i=1:length(cincludeList)
    code = strrep(code,includeString{i},'');
end

%% Process macro: include
allIncludeStrings = regexp(code,'(?<=(\n|^)[\s]*)include\(');
includeString = cell(1,length(allIncludeStrings));
% Insert include strings
for i=1:length(allIncludeStrings)
    % Find the closing bracket
    locStart = allIncludeStrings(i);
    locEnd = locStart;
    while code(locEnd)~=')'
        locEnd = locEnd+1;
    end
    if code(locEnd+1)==';'
        locEnd = locEnd+1;
    end
    includeString{i} = code(locStart:locEnd);
    % Get the file inside the bracket
    includeName = regexp(includeString{i},'(?<='').*(?='')','match');
    if ~isempty(includeName)
        includeFileContent{i} = fileread(includeName{1});
    end

end
for i=1:length(allIncludeStrings)
    code = strrep(code,includeString{i},includeFileContent{i});
end

%% Process macro: define
allDefineStrings = regexp(code,'#define');
macroName = cell(1,length(allDefineStrings));
macroValue = cell(1,length(allDefineStrings));
defineStringList = cell(1,length(allDefineStrings));
for i=1:length(allDefineStrings)
    % Find the closing bracket
    locStart = allDefineStrings(i);
    locEnd = locStart;
    while code(locEnd)~=char(13) && code(locEnd)~=char(10)
        locEnd = locEnd+1;
    end
    locEnd = locEnd-1;
    defineString = code(locStart:locEnd);
    defineStringList{i} = defineString;
    if strcmp(defineString(end),';')
        defineString = defineString(1:end-1);
    end
    % Split space
    defineStringSplit = strsplit(defineString,' ');
    if length(defineStringSplit)~=3
        error('#define should be followed by macro name, then macro value');
    end
    macroName{i} = defineStringSplit{2};
    macroValue{i} = defineStringSplit{3};
end
% Replace
for i=1:length(allDefineStrings)
    code = regexprep(code,defineStringList{i},'');
    code = regexprep(code,[BEFORE_A_WORD,macroName{i},AFTER_A_WORD],macroValue{i});
end

%% Process macro: #foreach ... #end block
code = rec_extract_loop_seg(code);

%% Process macro: #for ... #end block
allForStrings = regexp(code,'#for');
macroForString = cell(1,length(allForStrings));
macroExpands = cell(1,length(allForStrings));
for i=1:length(allForStrings)
    % Find the closing bracket
    locStart = allForStrings(i);
    locEnd = locStart+3;
    while strcmp(code(locEnd-3:locEnd),'#end')~=1
        locEnd = locEnd+1;
        if strcmp(code(locEnd-3:locEnd),'#for')==1
            error('Nested #for not supported yet');
        end
        
        if locEnd>=length(code)
            error('No closing #end found for #for block');
        end
    end

    forString = code(locStart:locEnd);
    macroForString{i} = forString;
    loopLines = strsplit(forString,{char(10),char(13)});
    % Extract the iterator
    loopHeader = loopLines{1};
    loopHeaderSplit = strsplit(loopHeader,' ');
    loopHeaderSplit = strsplit(loopHeaderSplit{2},'=');
    iteratorName = loopHeaderSplit{1};
    numList = eval(loopHeaderSplit{2:end});
    % for body
    forBody = strjoin(loopLines(2:end-1),LINE_BREAK);
    expandedForBody = '';
    for ii=1:length(numList)
        currentForBody = regexprep(forBody,['#',iteratorName],int2str(numList(ii)));
        currentForBody = strrep(currentForBody,['#(',iteratorName,'+1)'],int2str(numList(ii)+1));
        expandedForBody = [expandedForBody,currentForBody,LINE_BREAK];
    end
    macroExpands{i} = expandedForBody;
end
% Replace
for i=1:length(allForStrings)
    code = strrep(code, macroForString{i}, macroExpands{i});
end

%% Process #if ... #endif macro block
code = process_if_macro(code);

%% Extract blocks
% [equationCode, code] = extract_equation_segment(code);
[preModelCode,code] = extract_segment(code,'pre_model;');
[startLoopCode,code] = extract_segment(code,'start_loop;');
[finishLoopCode,code] = extract_segment(code,'finish_loop;');
[preJacCode,code] = extract_segment(code,'pre_jac_code;');
[postJacCode,code] = extract_segment(code,'post_jac_code;');
[modelCodes,modelConditionCodes,equationCodes,code] = extract_model_segment(code,'model',1);
[modelInitCodes,modelInitConditionCodes,equationInitCodes,code] = extract_model_segment(code,'model_init',0);
[simulateCode,code] = extract_simulate_segment(code);
[postInitCode,code] = extract_segment(code,'post_init;');
[preIterCode,code] = extract_segment(code,'pre_iter;');
[postIterCode,code] = extract_segment(code,'post_iter;');
[preInfCode,code] = extract_segment(code,'pre_inf;');
[postSolCode,code] = extract_segment(code,'post_sol;');
[preCallMexCode,code] = extract_segment(code,'pre_call_mex;');
[postCallMexCode,code] = extract_segment(code,'post_call_mex;');
[preMinorCode,code] = extract_segment(code,'pre_minor;');
[preUpdateCode,code] = extract_segment(code,'pre_update;');
[simuPreIterCode,code] = extract_segment(code,'pre_simulate_iter;');
[simuPostIterCode,code] = extract_segment(code,'post_simulate_iter;');

%% Extract var declaration
lines = strsplit(code,{char(10),char(13)});

parameters_name = {};
var_state_name = {};
var_shock_name = {};
var_tensor_name = {};

var_policy_name = {};
var_policy_length = [];
var_policy_loc = [];

var_interp_name = {};
var_interp_arg_name = {};
var_interp_arg_pos = {};
var_interp_noconverge = {};

var_aux_name = {};
var_aux_length = [];
var_aux_loc = [];

var_policy_init_name = {};
var_policy_init_length = [];
var_policy_init_loc = [];

var_aux_init_name = {};
var_aux_init_length = [];
var_aux_init_loc = [];

var_others_name = {};

var_output_name = {};

policy_loc = 1;
aux_loc = 1;
policy_init_loc = 1;
aux_init_loc = 1;

%% Process variables
% Rewrite code
code2 = '';
for j=1:length(lines)
    line = lines{j};
    % Remove comments
    [~,loc_start,loc_end] = regexp(line,'%.*$','tokenExtents','all');
    line(loc_start:loc_end)=[];
    
    % Split by simi-colon
    segs = strsplit(line,';');
    
    for i=1:length(segs)
        seg = segs{i};
        % remove trailing space
        seg = strtrim(seg);
        % Split by space
        words = strsplit(seg,' ');
        
        switch words{1}
            case 'parameters'
                for k=2:length(words)
                    parameters_name{end+1} = words{k};
                end
            case 'var_others'
                for k=2:length(words)
                    var_others_name{end+1} = words{k};
                end
            case 'var_output'
                for k=2:length(words)
                    var_output_name{end+1} = words{k};
                end
            case 'var_state'
                for k=2:length(words)
                    var_state_name{end+1} = words{k};
                end
            case 'var_shock'
                for k=2:length(words)
                    var_shock_name{end+1} = words{k};
                end
            case 'var_tensor'
                for k=2:length(words)
                    var_tensor_name{end+1} = words{k};
                end
            case 'var_policy'
                for k=2:length(words)
                    % Search length
                    locLeftBracket = regexp(words{k},'[');
                    if (locLeftBracket)
                        for locRightBracket=locLeftBracket+1:length(words{k})
                            if words{k}(locRightBracket)==']'
                                break;
                            end
                        end
                        var_length = eval(words{k}(locLeftBracket+1:locRightBracket-1));
                        var_policy_length(end+1) = var_length;
                        var_policy_name{end+1} = words{k}(1:locLeftBracket-1);
                        var_policy_loc(end+1) = policy_loc;
                        policy_loc = policy_loc+var_length;
                    else
                        var_policy_length(end+1) = 0;
                        var_policy_name{end+1} = words{k};
                        var_policy_loc(end+1) = policy_loc;
                        policy_loc = policy_loc+1;
                    end
                end
            case 'var_policy_init'
                for k=2:length(words)
                    % Search length
                    locLeftBracket = regexp(words{k},'[');
                    if (locLeftBracket)
                        for locRightBracket=locLeftBracket+1:length(words{k})
                            if words{k}(locRightBracket)==']'
                                break;
                            end
                        end
                        var_length = eval(words{k}(locLeftBracket+1:locRightBracket-1));
                        var_policy_init_length(end+1) = var_length;
                        var_policy_init_name{end+1} = words{k}(1:locLeftBracket-1);
                        var_policy_init_loc(end+1) = policy_init_loc;
                        policy_init_loc = policy_init_loc+var_length;
                    else
                        var_policy_init_length(end+1) = 0;
                        var_policy_init_name{end+1} = words{k};
                        var_policy_init_loc(end+1) = policy_init_loc;
                        policy_init_loc = policy_init_loc+1;
                    end
                end
            case 'var_interp'
                for k=2:length(words)
                    var_interp_arg_name{end+1} = {};
                    % Search interp state
                    locLeftBracket = regexp(words{k},'(');
                    if (locLeftBracket)
                        for locRightBracket=locLeftBracket+1:length(words{k})
                            if words{k}(locRightBracket)==')'
                                break;
                            end
                        end
                        % split the argument list
                        argWords = words{k}(locLeftBracket+1:locRightBracket-1);
                        var_interp_name{end+1} = words{k}(1:locLeftBracket-1);
                        % split with |
                        argWords = strsplit(argWords,'|');
                        for i_arg=1:length(argWords)
                            argWord = argWords{i_arg};
                            switch argWord
                                case 'noconverge'
                                    var_interp_noconverge{end+1} = var_interp_name{end};
                                otherwise
                                    var_interp_arg_name{end} = strsplit(argWord,',');
                            end
                        end
                    else
                        var_interp_name{end+1} = words{k};
                        var_interp_arg_name{end} = {};
                    end
                end
            case 'var_aux'
                for k=2:length(words)
                    % Search length
                    locLeftBracket = regexp(words{k},'[');
                    if (locLeftBracket)
                        for locRightBracket=locLeftBracket+1:length(words{k})
                            if words{k}(locRightBracket)==']'
                                break;
                            end
                        end
                        var_length = eval(words{k}(locLeftBracket+1:locRightBracket-1));
                        var_aux_length(end+1) = var_length;
                        var_aux_name{end+1} = words{k}(1:locLeftBracket-1);
                        var_aux_loc(end+1) = aux_loc;
                        aux_loc = aux_loc+var_length;
                    else
                        var_aux_length(end+1) = 0;
                        var_aux_name{end+1} = words{k};
                        var_aux_loc(end+1) = aux_loc;
                        aux_loc = aux_loc+1;
                    end
                end
            case 'var_aux_init'
                for k=2:length(words)
                    % Search length
                    locLeftBracket = regexp(words{k},'[');
                    if (locLeftBracket)
                        for locRightBracket=locLeftBracket+1:length(words{k})
                            if words{k}(locRightBracket)==']'
                                break;
                            end
                        end
                        var_length = eval(words{k}(locLeftBracket+1:locRightBracket-1));
                        var_aux_init_length(end+1) = var_length;
                        var_aux_init_name{end+1} = words{k}(1:locLeftBracket-1);
                        var_aux_init_loc(end+1) = aux_init_loc;
                        aux_init_loc = aux_init_loc+var_length;
                    else
                        var_aux_init_length(end+1) = 0;
                        var_aux_init_name{end+1} = words{k};
                        var_aux_init_loc(end+1) = aux_init_loc;
                        aux_init_loc = aux_init_loc+1;
                    end
                end
            case ''
            otherwise
                code2 = [code2 seg ';' LINE_BREAK];
        end
    end
end

% process automatic INTERP
for i=1:length(var_interp_name)
    current_var_interp_name = var_interp_name{i};
    if length(current_var_interp_name)>7 && strcmp(current_var_interp_name(end-6:end),'_INTERP')==1
        current_var_interp_name_short = current_var_interp_name(1:end-7);
        % Plug in var_aux if not yet specified
        if ~any(ismember(current_var_interp_name_short, var_aux_name)) && ~any(ismember(current_var_interp_name_short, var_policy_name))
            var_aux_length(end+1) = 0;
            var_aux_name{end+1} = current_var_interp_name_short;
            var_aux_loc(end+1) = aux_loc;
            aux_loc = aux_loc+1;
        end

        if ~any(ismember(current_var_interp_name_short, var_aux_init_name)) && ~any(ismember(current_var_interp_name_short, var_policy_init_name))
            var_aux_init_length(end+1) = 0;
            var_aux_init_name{end+1} = current_var_interp_name_short;
            var_aux_init_loc(end+1) = aux_init_loc;
            aux_init_loc = aux_init_loc+1;
        end
    end

    if length(current_var_interp_name)>9 && strcmp(current_var_interp_name(end-8:end),'_INTERP_E')==1
        current_var_interp_name_short = current_var_interp_name(1:end-9);
        % Plug in var_aux if not yet specified
        if ~any(ismember(current_var_interp_name_short, var_aux_name)) && ~any(ismember(current_var_interp_name_short, var_policy_name))
            var_aux_length(end+1) = 0;
            var_aux_name{end+1} = current_var_interp_name_short;
            var_aux_loc(end+1) = aux_loc;
            aux_loc = aux_loc+1;
        end

        if ~any(ismember(current_var_interp_name_short, var_aux_init_name)) && ~any(ismember(current_var_interp_name_short, var_policy_init_name))
            var_aux_init_length(end+1) = 0;
            var_aux_init_name{end+1} = current_var_interp_name_short;
            var_aux_init_loc(end+1) = aux_init_loc;
            aux_init_loc = aux_init_loc+1;
        end
    end
end

%% Information about variables
num_parameters = length(parameters_name);
num_state = length(var_state_name);
num_shock_var = length(var_shock_name);
num_tensor = length(var_tensor_name);
num_policy = length(var_policy_name);
num_policy_total = sum(max(var_policy_length,1));
num_interp = length(var_interp_name);
num_aux = length(var_aux_name);
num_aux_total = sum(max(var_aux_length,1));

% Add dummy if no state or no shock
if num_state==0 && num_shock_var==0
    % static problem
    is_static_problem = true;
else
    is_static_problem = false;
end

if num_state==0
    var_state_name = {'GDSGE_dummy_state'};
    num_state = 1;
end
if num_shock_var==0
    var_shock_name = {'GDSGE_dummy_shock'};
    num_shock_var = 1;
end

num_policy_init = length(var_policy_init_name);
num_policy_init_total = sum(max(var_policy_init_length,1));
num_aux_init = length(var_aux_init_name);
num_aux_init_total = sum(max(var_aux_init_length,1));

parameters_vec_name = strcat(parameters_name,{'(:)'});

var_pre_tensor_name = [var_state_name var_shock_name];

var_all_name = [parameters_name var_state_name var_shock_name var_tensor_name var_policy_name var_interp_name var_aux_name];
assert_valid_var_name(var_all_name);

var_init_all_name = [parameters_name var_state_name var_shock_name var_tensor_name var_policy_init_name var_aux_init_name];
assert_valid_var_name(var_init_all_name);

% var_policy is always in var_output
% var_output_name = unique([var_policy_name,var_output_name],'stable');
num_output = length(var_output_name);
var_output_length = zeros(1,num_output);
var_output_loc = zeros(1,num_output);
output_loc = 1;

if ~isempty(simulateCode)
    lines = strsplit(simulateCode,{char(10),char(13)});
    var_simulate_name = {};
    simulateCode2 = '';
    simulateCode1 = '';
    simuInitCode = '';
    simuInitOverwriteCode = '';
    shockInitialized = 0;
    stateInitialized = zeros(1,num_state);
    stateAssigned = zeros(1,num_state);
    for j=1:length(lines)
        % Remove comments
        line = lines{j};
        [~,loc_start,loc_end] = regexp(line,'%.*$','tokenExtents','all');
        line(loc_start:loc_end)=[];
        
        % Split by simi-colon
        segs = strsplit(line,';');
        for i=1:length(segs)
            seg = segs{i};
            seg = strtrim(seg);
            words = strsplit(seg,' ');
            switch words{1}
                case 'var_simu'
                    for k=2:length(words)
                        var_simulate_name{end+1} = words{k};
                    end
                    break;
                case 'initial'
                    simuInitCode = [simuInitCode ...
                        'SimuRslt.' words{2} '(:,1)=' words{3} ';' LINE_BREAK];
                    simuInitOverwriteCode = [simuInitOverwriteCode ...
                        'SimuRslt.' words{2} '(:,1:size(GDSGE_OPTIONS.init.' words{2} ',2))=GDSGE_OPTIONS.init.' words{2} ';' LINE_BREAK];
                    if strcmp(words{2},'shock')
                        shockInitialized = 1;
                        break;
                    end
                    whichMember = find(ismember(var_state_name,words{2}));
                    stateInitialized(whichMember) = 1;
                case ''
                case 'num_periods'
                    simulateCode1 = [simulateCode1 seg ';' LINE_BREAK];
                case 'num_samples'
                    simulateCode1 = [simulateCode1 seg ';' LINE_BREAK];
                otherwise
                    for k=1:num_state
                        if regexp(words{1},['(?<=(\W|^))' var_state_name{k} ''''])
                            stateAssigned(k) = 1;
                            
                            loopCode = regexprep(seg,['(?<=(\W|^))' var_state_name{k} ''''], ['SimuRslt.' var_state_name{k} '(:,GDSGE_t+1)']);
                            
                            % Split by =
                            loopCodeSplits = strsplit(loopCode,'=');
                            if length(loopCodeSplits)~=2
                                error('assignment error in simulation');
                            end
                            assignObject = loopCodeSplits{2};
                            
                            if regexp(assignObject,'''')
                                loopCode = regexprep(loopCode, '''', '(GDSGE_SHOCK_VAR_LINEAR_INDEX)');
                            elseif isempty(regexp(assignObject,'(','once'))
                                loopCode = [loopCode '(:)'];
                            else
                                loopCode = [loopCode];
                            end
                            
                            fullCode = [loopCode ';' LINE_BREAK];
                            simulateCode2 = [simulateCode2 fullCode LINE_BREAK];
                        end
                    end
            end
        end
    end
    
    assert_valid_var_name(var_simulate_name);
    num_simulate = length(var_simulate_name);
    
    % Write to var_output
    simulateAddToOutput = {''};
    for i=1:num_simulate
        var_simulate_name_i = var_simulate_name{i};
        simulate_in_state = find(strcmp(var_simulate_name_i, var_state_name));
        simulate_in_output = find(strcmp(var_simulate_name_i, var_output_name));
        if ~any(simulate_in_output) && ~any(simulate_in_state)
            var_output_name = [var_output_name, var_simulate_name_i];
            if num_output == 0
                var_output_loc = [var_output_loc, 1];
            else
                var_output_loc = [var_output_loc, var_output_loc(end)+1];
            end
            var_output_length = var_output_length + 1;
            num_output = num_output + 1;
            simulateAddToOutput = [simulateAddToOutput, var_simulate_name_i];
        end
    end
    % Print added
    if length(simulateAddToOutput)>1
        fprintf('The following var_simulate are added to var_output: %s\n', strjoin(simulateAddToOutput, ' '));
    end
end

% Write output var index
outputAddToAux = {''};
for i=1:num_output
    var_output_name_i = var_output_name{i};
    output_in_state = find(strcmp(var_output_name_i, var_state_name));
    if output_in_state
        error('state %s should not be in output var', var_output_name_i);
    end
    output_in_policy = find(strcmp(var_output_name_i, var_policy_name));
    if any(output_in_policy)
        var_output_length(i) = var_policy_length(output_in_policy);
        var_output_loc(i) = output_loc;
        output_loc = output_loc + max(var_output_length(i),1);
        continue;
    end
    output_in_aux = find(strcmp(var_output_name_i, var_aux_name));
    if any(output_in_aux)
        var_output_length(i) = var_aux_length(output_in_aux);
        var_output_loc(i) = output_loc;
        output_loc = output_loc + max(var_output_length(i),1);
        continue;
    end
    % Else, add output in aux
    
    var_aux_name = [var_aux_name, var_output_name_i];
    if num_aux==0
        var_aux_loc = [var_aux_loc, 1];
    else
        var_aux_loc = [var_aux_loc, var_aux_loc(end)+max(var_aux_length(end),1)];
    end
    var_aux_length = [var_aux_length, 0];
    num_aux = num_aux+1;
    num_aux_total = num_aux_total+1;
    outputAddToAux = [outputAddToAux, var_output_name_i];
    var_output_length(i) = 0;
    var_output_loc(i) = output_loc;
    output_loc = output_loc + 1;
    % error('output var %s not found in policy or aux', var_output_name_i);
end
% Print added
if length(outputAddToAux)>1
    fprintf('The following var_output are added to var_aux: %s\n', strjoin(outputAddToAux, ' '));
end
num_output_total = sum(max(var_output_length,1));

%% Construct tensors
tensorPrefix = 'GDSGE_TENSOR_';
tensorConstructCode = '';
tensor_state_name = var_state_name;
% Stack state variable in the order of definition
for j=1:length(var_state_name)
    tensor_state_name{j} = [tensorPrefix var_state_name{j} '_MODEL_ID'];
end
tensorConstructCode = [tensorConstructCode '[' tensorPrefix 'shockIdx_MODEL_ID,' ...
    my_strjoin(tensor_state_name,',') ']=ndgrid(1:shock_num,' ...
    my_strjoin(var_state_name,',') ');' LINE_BREAK];
% Stacking shock to the beginning
for j=1:length(var_shock_name)
    tensorConstructCode = [tensorConstructCode ...
        tensorPrefix var_shock_name{j} '=ndgrid(' ...
        var_shock_name{j} ',' my_strjoin(var_state_name,',') ');' LINE_BREAK];
end
tensorConstructCode = [tensorConstructCode ...
    'GDSGE_NPROB_MODEL_ID=numel(GDSGE_TENSOR_shockIdx_MODEL_ID);' LINE_BREAK];

%% policy init assign
policyInitAssignCode = '';
for j=1:num_policy_init
    policyInitAssignCode = [policyInitAssignCode ...
        var_policy_init_name{j} '=GDSGE_SOL_MODEL_ID(' num2str(var_policy_init_loc(j)) ':' num2str(var_policy_init_loc(j)+max(var_policy_init_length(j),1)-1) ...
        ',:);' LINE_BREAK];
end

%% policy assign
policyAssignCode = '';
for j=1:num_policy
    policyAssignCode = [policyAssignCode ...
        var_policy_name{j} '=GDSGE_SOL_MODEL_ID(' num2str(var_policy_loc(j)) ':' num2str(var_policy_loc(j)+max(var_policy_length(j),1)-1) ...
        ',:);' LINE_BREAK];
end

%% interp assign
interpAssignCode = '';
for j=1:num_interp
    interpAssignCode = [interpAssignCode ...
        var_interp_name{j} '=GDSGE_NEW_' var_interp_name{j} ';' LINE_BREAK];
end

%% output index assign code
outputIndexAssignCode = ['output_var_index=struct();',LINE_BREAK];
for j=1:num_output
    outputIndexAssignCode = [outputIndexAssignCode ...
        'output_var_index.' var_output_name{j} '=' num2str(var_output_loc(j)) ':' num2str(var_output_loc(j)+max(var_output_length(j),1)-1) ...
        ';' LINE_BREAK];
end

%% inbound and initial
lines = strsplit(code2,{LINE_BREAK,char(13)});
code3 = '';
tensorAssignCode = '';
interpNewAssignCode = '';
interpAsgAssignCode = '';
inboundCode = '';
inboundInitCode = '';
initialCode = '';

initialCode = [initialCode ...
    'GDSGE_X0_MODEL_ID = (GDSGE_LB_MODEL_ID+GDSGE_UB_MODEL_ID)/2;' LINE_BREAK ...
    ];

inboundInitCode = [inboundInitCode ...
    'GDSGE_LB_MODEL_ID = zeros(' num2str(num_policy_init_total) ',GDSGE_NPROB_MODEL_ID);' LINE_BREAK ...
    'GDSGE_UB_MODEL_ID = 1e3*ones(' num2str(num_policy_init_total) ',GDSGE_NPROB_MODEL_ID);' LINE_BREAK ...
    ];

inboundCode = [inboundCode ...
    'GDSGE_LB_MODEL_ID = zeros(' num2str(num_policy_total) ',GDSGE_NPROB_MODEL_ID);' LINE_BREAK ...
    'GDSGE_UB_MODEL_ID = 1e3*ones(' num2str(num_policy_total) ',GDSGE_NPROB_MODEL_ID);' LINE_BREAK ...
    ];
warmupBoundAdaptiveCode = '';
inboundAdaptiveCode = '';
inboundAdaptiveTightCode = '';

policyInitInitializeCode = '';
policyInitInitializeCode = [policyInitInitializeCode ...
    'GDSGE_X0_MODEL_ID = (GDSGE_LB_MODEL_ID+GDSGE_UB_MODEL_ID)/2;' LINE_BREAK ...
    ];

policyInitializeCode = '';
%{
policyInitializeCode = [policyInitializeCode ...
    'GDSGE_X0_MODEL_ID = (GDSGE_LB_MODEL_ID+GDSGE_UB_MODEL_ID)/2;' LINE_BREAK ...
    ];
%}
interpInitializeCode = '';
for j=1:num_interp
    % Check whether add default initialization and assign
    if length(var_interp_name{j})>7 && strcmp(var_interp_name{j}(end-6:end),'_INTERP')==1
        % Plug in default initialization and assignment
        interpInitializeCode = [interpInitializeCode ...
            var_interp_name{j} '=' var_interp_name{j}(1:end-7) ';' LINE_BREAK];
        interpNewAssignCode = [interpNewAssignCode ...
            'GDSGE_NEW_' var_interp_name{j} '=' var_interp_name{j}(1:end-7) ';' LINE_BREAK];
        continue;
    else
        interpInitializeCode = [interpInitializeCode ...
            var_interp_name{j} '=zeros(GDSGE_SIZE);' LINE_BREAK];
        continue;
    end
    
    if length(var_interp_name{j})>9 && strcmp(var_interp_name{j}(end-8:end),'_INTERP_E')==1
        % Plug in default initialization and assignment
        interpInitializeCode = [interpInitializeCode ...
            var_interp_name{j} '=shock_trans*reshape(' var_interp_name{j}(1:end-9) ',shock_num,[]);' LINE_BREAK];
        interpNewAssignCode = [interpNewAssignCode ...
            'GDSGE_NEW_' var_interp_name{j} '=shock_trans*reshape(' var_interp_name{j}(1:end-9) ',shock_num,[]);' LINE_BREAK];
        continue;
    else
        interpInitializeCode = [interpInitializeCode ...
            var_interp_name{j} '=zeros(GDSGE_SIZE);' LINE_BREAK];
        continue;
    end
end
asgInitializeCode = '';
tensor_tensor_name = {};
for j=1:length(var_tensor_name)
    tensor_tensor_name{j} = [var_tensor_name{j}];
end

for j=1:length(lines)
    try
    % Split by simi-colon
    line = lines{j};
    segs = strsplit(line,';');
    
    for i=1:length(segs)
        seg = segs{i};
        % Split by space
        words = strsplit(seg,{' ','='});
        % Assign tensor var
        if find(ismember(var_tensor_name,words{1}))
            % Replace tensor name
            for k=1:length(var_pre_tensor_name)
                if ~ismember(var_pre_tensor_name{k}, var_tensor_name)
                    seg = regexprep(seg,['((?<=^)|(?<=\W))' var_pre_tensor_name{k} '((?=$)|(?=\W))'],[tensorPrefix var_pre_tensor_name{k}],'all');
                end
            end
            
            tensorAssignCode = [tensorAssignCode ...
                seg ';' LINE_BREAK];
            continue;
        end
        
        % Assign interp var
        if find(ismember(var_interp_name,words{1}))
            % Replace tensor name
            for k=1:length(var_pre_tensor_name)
                if ~ismember(var_pre_tensor_name{k}, var_tensor_name)
                    seg = regexprep(seg,['((?<=^)|(?<=\W))' var_pre_tensor_name{k} '((?=$)|(?=\W))'],[tensorPrefix var_pre_tensor_name{k}],'all');
                end
            end
            
            % Replace policy name
            for k=1:length(var_policy_name)
                seg = regexprep(seg,['((?<=^)|(?<=\W))' var_policy_name{k} '((?=$)|(?=\W))'],[var_policy_name{k}],'all');
            end
            
            interpAsgAssignCode = [interpAsgAssignCode ...
                seg ';' LINE_BREAK];
            
            segSplitBySpace = strsplit(seg,'=');
            
            new_seg = [regexprep(segSplitBySpace{1}, words{1}, ['GDSGE_NEW_' words{1}]), '=', segSplitBySpace{2:end}];
            interpNewAssignCode = [interpNewAssignCode ...
                new_seg ';' LINE_BREAK];
            continue;
        end
        
        switch words{1}
            case 'inbound_init'
                if length(words)~=4
                    error('inbound_init needs exact 4 arguments');
                end
                var_name = words{2};
                var_name_idx = find(ismember(var_policy_init_name,var_name));
                % Substitute tensor
                words{3} = replace_tensor_name(words{3},var_pre_tensor_name,tensorPrefix);
                words{3} = replace_tensor_name(words{3},var_tensor_name,'');
                words{4} = replace_tensor_name(words{4},var_pre_tensor_name,tensorPrefix);
                words{4} = replace_tensor_name(words{4},var_tensor_name,'');
                % Account for policy length
                for loc=1:max(1,var_policy_init_length(var_name_idx))
                    inboundInitCode = [inboundInitCode ...
                        'GDSGE_LB_MODEL_ID(' num2str(var_policy_init_loc(var_name_idx)+loc-1) ',:)=' words{3} ';' LINE_BREAK...
                        ];
                    inboundInitCode = [inboundInitCode ...
                        'GDSGE_UB_MODEL_ID(' num2str(var_policy_init_loc(var_name_idx)+loc-1) ',:)=' words{4} ';' LINE_BREAK...
                        ];
                end
            case 'inbound'
                if ~(length(words)>=4)
                    error('inbound needs more than 4 arguments');
                end
                var_name = words{2};
                var_name_idx = find(ismember(var_policy_name,var_name));
                % Substitute tensor
                words{3} = replace_tensor_name(words{3},var_pre_tensor_name,tensorPrefix);
                words{3} = replace_tensor_name(words{3},var_tensor_name,'');
                words{4} = replace_tensor_name(words{4},var_pre_tensor_name,tensorPrefix);
                words{4} = replace_tensor_name(words{4},var_tensor_name,'');
                % Account for policy length
                inboundCode = [inboundCode ...
                        'GDSGE_LB_MODEL_ID(' num2str(var_policy_loc(var_name_idx)) ':' num2str(var_policy_loc(var_name_idx)+max(var_policy_length(var_name_idx),1)-1) ',:)=' words{3} ';' LINE_BREAK...
                        ];
                inboundCode = [inboundCode ...
                        'GDSGE_UB_MODEL_ID(' num2str(var_policy_loc(var_name_idx)) ':' num2str(var_policy_loc(var_name_idx)+max(var_policy_length(var_name_idx),1)-1) ',:)=' words{4} ';' LINE_BREAK...
                        ];
                if (length(words)>4)
                    % Have some additional options
                    for iwords=5:length(words)
                        found = regexp(words{iwords},'adaptive');
                        if found
                            multiplier = regexp(words{iwords},'(?<=adaptive\().*(?=\))','match');
                        else
                            error('Lower and upper bounds following inbound should not contain spaces');
                        end
                        if length(multiplier)>0
                            multiplier = multiplier{1};
                            % Adaptive upper and lower bound
                            for loc=1:max(1,var_policy_length(var_name_idx))
                                inboundAdaptiveCode = [inboundAdaptiveCode ...
                                    'GDSGE_LB_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)=' ...
                                    'min(GDSGE_LB_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:),' ...
                                    'GDSGE_SOL_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)*' multiplier ');' LINE_BREAK...
                                    ];
                                
                                inboundAdaptiveCode = [inboundAdaptiveCode ...
                                    'GDSGE_UB_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)=' ...
                                    'max(GDSGE_UB_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:),' ...
                                    'GDSGE_SOL_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)*' multiplier ');' LINE_BREAK...
                                    ];
                                
                                inboundAdaptiveTightCode = [inboundAdaptiveTightCode ...
                                    'GDSGE_LB_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)=' ...
                                    'GDSGE_SOL_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)/' multiplier ';' LINE_BREAK...
                                    ];
                                
                                inboundAdaptiveTightCode = [inboundAdaptiveTightCode ...
                                    'GDSGE_UB_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)=' ...
                                    'GDSGE_SOL_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)*' multiplier ';' LINE_BREAK...
                                    ];
                                
                                warmupBoundAdaptiveCode = [warmupBoundAdaptiveCode ...
                                    'GDSGE_LB_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)=' ...
                                    'GDSGE_LB_NEW_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:);' LINE_BREAK...
                                    ];
                                warmupBoundAdaptiveCode = [warmupBoundAdaptiveCode ...
                                    'GDSGE_UB_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)=' ...
                                    'GDSGE_UB_NEW_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:);' LINE_BREAK...
                                    ];
                                    
                                %{
                                inboundAdaptiveCode = [inboundAdaptiveCode ...
                                    'GDSGE_UB_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)=' ...
                                    'GDSGE_SOL_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)*' multiplier ';' LINE_BREAK...
                                    ];
                                %}
                            end
                        else
                            error('Missing argument in adaptive');
                        end
                    end
                end
            case 'initial'
                if ~(length(words)>=3)
                    error('initial needs more than 3 arguments');
                end
                var_name = words{2};
                var_name_idx = find(ismember(var_policy_name,var_name));
                words{3} = strcat(words{3:end});
                words = words(1:3);
                if var_name_idx
                    % Substitute tensor
                    words{3} = replace_tensor_name(words{3},var_pre_tensor_name,tensorPrefix);
                    words{3} = replace_tensor_name(words{3},var_tensor_name,'');
                    for loc=1:max(1,var_policy_length(var_name_idx))
                        policyInitializeCode = [policyInitializeCode ...
                            'GDSGE_X0_MODEL_ID(' num2str(var_policy_loc(var_name_idx)+loc-1) ',:)=' words{3:end} ';' LINE_BREAK...
                            ];
                    end
                else
                    var_name_idx = find(ismember(var_interp_name,var_name));
                    if var_name_idx
                        % Substitute tensor
                        words{3} = replace_tensor_name(words{3},var_pre_tensor_name,tensorPrefix);
                        words{3} = replace_tensor_name(words{3},var_tensor_name,'');
                        interpInitializeCode = [interpInitializeCode ...
                            words{2} '(:)=' words{3:end} ';' LINE_BREAK...
                            ];
                    else
                        error('undefined variable in initial: %s', words{2});
                    end
                end
            case 'shock_num'
                if ~(length(words)==2)
                    error('Line starting with shock_num should specify shock_num with a constant literal.');
                end
                shock_num = str2num(words{2});
                code3 = [code3 seg ';' LINE_BREAK];
            case ''
            otherwise
                code3 = [code3 seg ';' LINE_BREAK];
        end
    end
    catch ME
        fprintf(2,'Error found outside model block, at Line\n%s\n', line);
        rethrow(ME);
    end
end

%% Evaluate pre Code
[parameters_size, pre_params,params] = get_parameters(code3,parameters_name);
v2struct(pre_params);
shock_num = pre_params.shock_num;

%% Check validity of parameters
% Check valid of preCode
if USE_SPLINE+USE_ASG+USE_PCHIP~=1
    error('specify one of spline, pchip or asg');
end

if USE_ASG && num_tensor>0
    error('cannot use var_tensor when using adpative grids');
end

if USE_PCHIP && num_state>1
    error('pchip can only used with single state');
end

if INTERP_ORDER~=4 && INTERP_ORDER~=2
    error('INTERP_ORDER should either 4 (cubic) or 2 (linear)');
end

%% Compute metric between interpolation object
metricCode = ['GDSGE_Metric = max([' char(13)];
for j=1:num_interp
    if ~any(strcmp(var_interp_name{j},var_interp_noconverge))
        metricCode = [metricCode ...
            'max(abs(' ...
            'GDSGE_NEW_' var_interp_name{j} '(:)-' var_interp_name{j} '(:)))' char(13)];
    end
end
metricCode = [metricCode ']);' LINE_BREAK];

%% aux var code
auxAssignCode = '';
headerAuxAssignCode = '';
for j=1:num_aux
    if var_aux_length(j)==0
        auxAssignCode = [auxAssignCode ...
            var_aux_name{j} '=GDSGE_AUX_MODEL_ID(' num2str(var_aux_loc(j)) ',:);' LINE_BREAK];
        
        headerAuxAssignCode = [headerAuxAssignCode ...
            'GDSGE_aux[' num2str(var_aux_loc(j)-1) ']=value(' var_aux_name{j} ');' LINE_BREAK];
    else
        auxAssignCode = [auxAssignCode ...
            var_aux_name{j} '=GDSGE_AUX_MODEL_ID(' num2str(var_aux_loc(j)) ':' num2str(var_aux_loc(j)+max(var_aux_length(j),1)-1) ...
            ',:);' LINE_BREAK];
        
        loopCode = ['for(int GDSGE_iter=1; GDSGE_iter<=' num2str(var_aux_length(j)) '; GDSGE_iter++)' LINE_BREAK '{' LINE_BREAK];
        loopCode = [loopCode 'GDSGE_aux[' num2str(var_aux_loc(j)-2) '+GDSGE_iter]=value(' var_aux_name{j} '_GDSGE_GRID[GDSGE_iter-1]);' LINE_BREAK];
        loopCode = [loopCode '}' LINE_BREAK];
        headerAuxAssignCode = [headerAuxAssignCode loopCode];
    end
end

%% aux init code
auxInitAssignCode = '';
headerAuxInitAssignCode = '';
for j=1:num_aux_init
    if var_aux_init_length(j)==0
        auxInitAssignCode = [auxInitAssignCode ...
            var_aux_init_name{j} '=GDSGE_AUX_MODEL_ID(' num2str(var_aux_init_loc(j)) ',:);' LINE_BREAK];
        
        headerAuxInitAssignCode = [headerAuxInitAssignCode ...
            'GDSGE_aux[' num2str(var_aux_init_loc(j)-1) ']=value(' var_aux_init_name{j} ');' LINE_BREAK];
    else
        auxInitAssignCode = [auxInitAssignCode ...
            var_aux_init_name{j} '=GDSGE_AUX_MODEL_ID(' num2str(var_aux_init_loc(j)) ':' num2str(var_aux_init_loc(j)+max(var_aux_init_length(j),1)-1) ...
            ',:);' LINE_BREAK];
        
        loopCode = ['for(int GDSGE_iter=1; GDSGE_iter<=' num2str(shock_num) '; GDSGE_iter++)' LINE_BREAK '{' LINE_BREAK];
        loopCode = [loopCode 'GDSGE_aux[' num2str(var_aux_init_loc(j)-2) '+GDSGE_iter]=value(' var_aux_init_name{j} '_GDSGE_GRID[GDSGE_iter-1]);' LINE_BREAK];
        loopCode = [loopCode '}' LINE_BREAK];
        headerAuxInitAssignCode = [headerAuxInitAssignCode loopCode];
    end
end

%% Reshape code
reshapeCode = '';
for j=1:num_policy
    if var_policy_length(j)==0
        reshapeCode = [reshapeCode ...
            var_policy_name{j} '=reshape(' var_policy_name{j} ',GDSGE_SIZE);' LINE_BREAK];
    else
        reshapeCode = [reshapeCode ...
            var_policy_name{j} '=reshape(' var_policy_name{j} ',[' num2str(var_policy_length(j)) ' GDSGE_SIZE]);' LINE_BREAK];
    end
end
for j=1:num_interp
    reshapeCode = [reshapeCode ...
        var_interp_name{j} '=reshape(' var_interp_name{j} ',GDSGE_SIZE);' LINE_BREAK];
end
for j=1:num_aux
    if var_aux_length(j)==0
        reshapeCode = [reshapeCode ...
            var_aux_name{j} '=reshape(' var_aux_name{j} ',GDSGE_SIZE);' LINE_BREAK];
    else
        reshapeCode = [reshapeCode ...
            var_aux_name{j} '=reshape(' var_aux_name{j} ',[' num2str(var_aux_length(j)) ' GDSGE_SIZE]);' LINE_BREAK];
    end
end

% For init
reshapeInitCode = '';
for j=1:num_policy_init
    if var_policy_init_length(j)==0
        reshapeInitCode = [reshapeInitCode ...
            var_policy_init_name{j} '=reshape(' var_policy_init_name{j} ',GDSGE_SIZE);' LINE_BREAK];
    else
        reshapeInitCode = [reshapeInitCode ...
            var_policy_init_name{j} '=reshape(' var_policy_init_name{j} ',[' num2str(var_policy_init_length(j)) ' GDSGE_SIZE]);' LINE_BREAK];
    end
end
for j=1:num_aux_init
    if var_aux_init_length(j)==0
        reshapeInitCode = [reshapeInitCode ...
            var_aux_init_name{j} '=reshape(' var_aux_init_name{j} ',GDSGE_SIZE);' LINE_BREAK];
    else
        reshapeInitCode = [reshapeInitCode ...
            var_aux_init_name{j} '=reshape(' var_aux_init_name{j} ',[' num2str(var_aux_init_length(j)) ' GDSGE_SIZE]);' LINE_BREAK];
    end
end

%% SubReshape code
subReshapeCode = '';
for j=1:num_policy
    subReshapeCode = [subReshapeCode ...
        var_policy_name{j} '=reshape(' var_policy_name{j} ',shock_num,[]);' LINE_BREAK];
end
for j=1:num_aux
    subReshapeCode = [subReshapeCode ...
        var_aux_name{j} '=reshape(' var_aux_name{j} ',shock_num,[]);' LINE_BREAK];
end

%% Assert something
assertCode = '';
%
make_assertion = @(condition) ['assert(' condition ')'] ;
% assert length of shock var
assertCode = [assertCode make_assertion('exist(''shock_num'',''var'')==1') ';' LINE_BREAK];
for j=1:num_shock_var
    assertCode = [assertCode ...
        make_assertion(['length(' var_shock_name{j} ')==' num2str(shock_num)]) ';' LINE_BREAK ...
        ];
end
assertCode = [assertCode make_assertion(['size(shock_trans,1)==' num2str(shock_num)]) ';' LINE_BREAK];
assertCode = [assertCode make_assertion(['size(shock_trans,2)==' num2str(shock_num)]) ';' LINE_BREAK];

%% Macro for interpolation
% Process interp arg
% var_interp_arg_pos = cell(1,num_interp);
interpGetCode = '';
interpGetThreadCode = '';
interpConstructCode = '';
if num_interp>0
for j=1:num_interp
    if isempty(var_interp_arg_name{j})
        var_interp_arg_name{j} = var_state_name;
    else
        if USE_ASG
            error('interp var with position arguments not implemented with ASG yet');
        end
    end
end

% Template
if ~USE_ASG
    interpGetTemplate = fileread([template_folder '/interp_spline_get_template.cpp']);
else
    interpGetTemplate = fileread([template_folder '/interp_asg_get_template.cpp']);
end
if USE_ASG
    interpGetCode = fileread([template_folder '/interp_asg_construct_template.cpp']);
    interpGetThreadCode = fileread([template_folder '/interp_asg_prepare_space_template.cpp']);
    % Do an inplace substitution
    interpGetThreadCode = regexprep(interpGetThreadCode,'ASG_MAX_DIM',num2str(num_state));
    interpGetThreadCode = regexprep(interpGetThreadCode,'ASG_MAX_NVEC',num2str(num_interp));
elseif USE_SPLINE || USE_PCHIP
    interpGetCode = fileread([template_folder '/interp_spline_construct_template.cpp']);
    interpGetThreadCode = fileread([template_folder '/interp_spline_prepare_space_template.cpp']);
    % Do an inplace substitution
    first_interp_arg_length = length(var_interp_arg_name{1});
    
    interpGetCode = regexprep(interpGetCode,'VAR_NUM',num2str(first_interp_arg_length));
    interpGetCode = regexprep(interpGetCode,'NUM_INTERP',num2str(num_interp));
    
    interpGetThreadCode = regexprep(interpGetThreadCode,'GDSGE_INTERP_XDIM',num2str(first_interp_arg_length));
    interpGetThreadCode = regexprep(interpGetThreadCode,'GDSGE_NUM_INTERP',num2str(num_interp));
    
    % Ref interpolation objects
    refInterpVarList = my_strjoin(strcat('&GDSGE_CPP_',var_interp_name),',');
    interpGetThreadCode = regexprep(interpGetThreadCode,'REF_INTERP_VAR_LIST',refInterpVarList);
    
    % Call the search function
    interpGetThreadCode = regexprep(interpGetThreadCode,'FIRST_INTERP_CPP_OBJECT',['GDSGE_CPP_' var_interp_name{1}]);
    
    % call each interp var evaluation
    interpEvalTemplate = [...
        'if (INTERP_VEC_FLAG[interpIdx++]) {' LINE_BREAK ...
        '  if (!hasSearched) { GDSGE_CPP_INTERP_NAME.search(xSite, GDSGE_INTERP_CELL, GDSGE_XSITE_TO_LEFT); hasSearched=true; }' LINE_BREAK ...
        '  GDSGE_INTERP_RSLT[rsltIdx++]=GDSGE_CPP_INTERP_NAME.nosearch_eval_INTERP_ORDER(GDSGE_XSITE_TO_LEFT, GDSGE_INTERP_CELL, shockIdx-1);' LINE_BREAK ...
        '}' LINE_BREAK ...
        ];
    interpEvalCodeAll = '';
    % Do vector evaluation only for first several full state
    % TODO(wenlan): this is left to be better implemented
    first_interp_arg_length = length(var_interp_arg_name{1});
    for i=1:num_interp
        if length(var_interp_arg_name{i})~=first_interp_arg_length
            break;
        end
        interpEvalCode = ['' regexprep(interpEvalTemplate,'INTERP_NAME',var_interp_name{i})];
        interpEvalCodeAll = [interpEvalCodeAll,interpEvalCode,LINE_BREAK];
    end
    
    interpGetThreadCode = regexprep(interpGetThreadCode,'LOOP_OVER_ALL_INTERP_EVALUATION_CODE',interpEvalCodeAll);
end

adouble_var_state_name = strcat({'adouble '}, var_state_name);
double_var_state_name = strcat({'double '}, var_state_name);

for j=1:num_interp
    fillInTemplate = interpGetTemplate;
    fillInTemplate = regexprep(fillInTemplate,'INTERP_IDX',num2str(j-1));
    fillInTemplate = regexprep(fillInTemplate,'INTERP_NAME',var_interp_name{j});
    if USE_ASG
        interpGetThreadCode = [interpGetThreadCode fillInTemplate LINE_BREAK];
    else
        fillInTemplate = regexprep(fillInTemplate,'ADOUBLE_VAR_NAME',my_strjoin(strcat({'adouble '}, var_interp_arg_name{j}),','));
        fillInTemplate = regexprep(fillInTemplate,'DOUBLE_VAR_NAME',my_strjoin(strcat({'double '}, var_interp_arg_name{j}),','));
        fillInTemplate = regexprep(fillInTemplate,'VAR_NAME',my_strjoin(var_interp_arg_name{j},','));
        fillInTemplate = regexprep(fillInTemplate,'VAR_NUM',num2str(length(var_interp_arg_name{j})));
        fillInTemplate = regexprep(fillInTemplate,'NUM_INTERP',num2str(num_interp));
        fillInTemplate = regexprep(fillInTemplate,'INTERP_ORDER',num2str(INTERP_ORDER));
        
        interpGetCode = [interpGetCode fillInTemplate LINE_BREAK];
    end
end
if USE_ASG
    interpGetThreadCode = regexprep(interpGetThreadCode,'ADOUBLE_VAR_NAME',my_strjoin(adouble_var_state_name,','));
    interpGetThreadCode = regexprep(interpGetThreadCode,'DOUBLE_VAR_NAME',my_strjoin(double_var_state_name,','));
    interpGetThreadCode = regexprep(interpGetThreadCode,'VAR_NAME',my_strjoin(var_state_name,','));
    interpGetThreadCode = regexprep(interpGetThreadCode,'VAR_NUM',num2str(num_state));
else
    adouble_var_interp_arg_name_first = strcat({'adouble '}, var_interp_arg_name{1});
    double_var_interp_arg_name_first = strcat({'double '}, var_interp_arg_name{1});
    
    interpGetThreadCode = regexprep(interpGetThreadCode,'ADOUBLE_VAR_NAME',my_strjoin(adouble_var_interp_arg_name_first,','));
    interpGetThreadCode = regexprep(interpGetThreadCode,'DOUBLE_VAR_NAME',my_strjoin(double_var_interp_arg_name_first,','));
    interpGetThreadCode = regexprep(interpGetThreadCode,'VAR_NAME',my_strjoin(var_interp_arg_name{1},','));
    interpGetThreadCode = regexprep(interpGetThreadCode,'VAR_NUM',num2str(num_state));
    interpGetThreadCode = regexprep(interpGetThreadCode,'INTERP_ORDER',num2str(INTERP_ORDER));
end

% Template for construction
if USE_SPLINE==1
    if INTERP_ORDER==4 && EXTRAP_ORDER==2
        interpConstructTemplate = [
            'GDSGE_PP_INTERP_NAME=struct(''form'',''pp'',''breaks'',{{STATE_COMMA}},''Values'',reshape(INTERP_NAME,[],GDSGE_SIZE_STATE_MODEL_ID{:}),''coefs'',[],''order'',GDSGE_INTERP_ORDER_MODEL_ID,''Method'',[],''ExtrapolationOrder'',GDSGE_EXTRAP_ORDER,''thread'',NumThreads,''orient'',''curvefit'');' LINE_BREAK ...
            'GDSGE_PP_INTERP_NAME=myppual(myppual(GDSGE_PP_INTERP_NAME));' LINE_BREAK ...
            ];
    elseif (INTERP_ORDER==4 && EXTRAP_ORDER==4) || INTERP_ORDER==2
        interpConstructTemplate = [
            'GDSGE_PP_INTERP_NAME=struct(''form'',''MKL'',''breaks'',{{STATE_COMMA}},''Values'',reshape(INTERP_NAME,[],GDSGE_SIZE_STATE_MODEL_ID{:}),''coefs'',[],''order'',GDSGE_INTERP_ORDER_MODEL_ID,''Method'',[],''ExtrapolationOrder'',GDSGE_EXTRAP_ORDER,''thread'',NumThreads,''orient'',''curvefit'');' LINE_BREAK ...
            'GDSGE_PP_INTERP_NAME=myppual(GDSGE_PP_INTERP_NAME);' LINE_BREAK ...
            ];
    end
elseif USE_PCHIP==1
    interpConstructTemplate = [
        'GDSGE_PP_INTERP_NAME=pchip(STATE_COMMA,reshape(INTERP_NAME,[],GDSGE_SIZE_STATE_MODEL_ID{:}));' LINE_BREAK ...
        'GDSGE_PP_INTERP_NAME=myppual(GDSGE_PP_INTERP_NAME);' LINE_BREAK ...
        ];
else
    interpConstructTemplate = '';
end
pp_interp_name = strcat({'GDSGE_PP_'}, var_interp_name);
interpConstructCode = '';
for j=1:num_interp
    fillInTemplate = interpConstructTemplate;
    fillInTemplate = regexprep(fillInTemplate,'INTERP_NAME',var_interp_name{j});
    fillInTemplate = regexprep(fillInTemplate,'STATE_COMMA',my_strjoin(var_interp_arg_name{j},','));
    interpConstructCode = [interpConstructCode fillInTemplate LINE_BREAK];
end
else
    pp_interp_name = {'GDSGE_EMPTY'};
end

%% Pop code
% pop numShocks
popCode = '';
popCode = [popCode 'POPN(GDSGE_NUM_SHOCKS_D);' LINE_BREAK];
popCode = [popCode 'int GDSGE_NUM_SHOCKS = (int) GDSGE_NUM_SHOCKS_D;' LINE_BREAK];
% pop parameters
for j=1:num_parameters
    if parameters_size(j)==1
        popCode = [popCode ...
            'POPN(' parameters_name{j} ');' LINE_BREAK
            ];
    else
        popCode = [popCode ...
            'POPNARRAY(' parameters_name{j} '_GDSGE_GRID' ',' num2str(parameters_size(j)) ');' LINE_BREAK ...
            '#define ' parameters_name{j} '(idx) ' parameters_name{j} '_GDSGE_GRID' '[int(idx)-1]' LINE_BREAK
            ];
    end
end
% pop transition matrix
popCode = [popCode 'POPNARRAY(GDSGE_shock_trans,',num2str(shock_num),'*',num2str(shock_num),');' LINE_BREAK ...
    '#define ' 'shock_trans(idx) GDSGE_shock_trans[int(idx)-1]' LINE_BREAK];
for j=1:num_shock_var
    popCode = [popCode 'POPNARRAY(' var_shock_name{j} '_GDSGE_GRID,' num2str(shock_num) ');' LINE_BREAK ...
        '#define ' var_shock_name{j} '_GRID(idx) ' var_shock_name{j} '_GDSGE_GRID[int(idx)-1]' LINE_BREAK];
end

% pop shock code
popCode = [popCode 'POPN(GDSGE_shockIdx_d);' LINE_BREAK];
popCode = [popCode 'int shock = (int)GDSGE_shockIdx_d;' LINE_BREAK];
for j=1:num_shock_var
    popCode = [popCode 'double ' var_shock_name{j} ' = ' var_shock_name{j} '_GRID(shock);' LINE_BREAK];
end

% pop state variable
for j=1:num_state
    popCode = [popCode ...
        'POPN(' var_state_name{j} ');' LINE_BREAK
        ];
end

% pop tensor variable
for j=1:num_tensor
    popCode = [popCode ...
        'POPN(' var_tensor_name{j} ');' LINE_BREAK
        ];
end

%% Process model code
gen_cxx_model_code_inline = @(modelCode,equationCode) gen_cxx_model_code(modelCode,equationCode,...
    num_interp,var_interp_name,num_policy,var_policy_name,num_policy_total,var_policy_length,var_policy_loc,...
    var_aux_name, var_aux_length, num_shock_var,var_shock_name,shock_num,LINE_BREAK);
[modelAllCodes,callFMinCodes,num_equations] = process_model_code(modelCodes,equationCodes,modelConditionCodes,gen_cxx_model_code_inline,...
    headerAuxAssignCode,template_folder,USE_FINITE_DIFF,LINE_BREAK);
% patch it right way


gen_cxx_model_code_no_var_aux_inline = @(modelCode,equationCode) gen_cxx_model_code(modelCode,equationCode,...
    num_interp,var_interp_name,num_policy,var_policy_name,num_policy_total,var_policy_length,var_policy_loc,...
    {''}, 0, num_shock_var,var_shock_name,shock_num,LINE_BREAK);
preJacCode = gen_cxx_model_code_no_var_aux_inline(preJacCode,'');
postJacCode = gen_cxx_model_code_no_var_aux_inline(postJacCode,'');
modelAllCodes = my_regexprep(modelAllCodes, 'PRE_JAC_CODE', preJacCode);
modelAllCodes = my_regexprep(modelAllCodes, 'POST_JAC_CODE', postJacCode);
[preModelCode,~,~,preModelDeclareCode] = gen_cxx_model_code_no_var_aux_inline(preModelCode,'');
preModelCode = regexprep(preModelCode,'((?<=^)|(?<=\W)|(?<=_))adouble(?=\W)','double');
% Manually process preModelCode
preModelDeclareCode = regexprep(preModelDeclareCode,'((?<=^)|(?<=\W)|(?<=_))adouble(?=\W)','double');
var_interp_name_expanded = [var_interp_name,'GDSGE_INTERP_VEC','GDSGE_INTERP_RSLT'];
for ii = 1:length(var_interp_name_expanded)
    preModelCode = regexprep(preModelCode,['((?<=^)|(?<=\W)',var_interp_name_expanded{ii},'(?=\W)'],[var_interp_name_expanded{ii},'_double']);
end
preModelAllCode = [preModelDeclareCode,LINE_BREAK,preModelCode];

startLoopCode = gen_cxx_model_code_no_var_aux_inline(startLoopCode,'');
finishLoopCode = gen_cxx_model_code_no_var_aux_inline(finishLoopCode,'');

if num_policy_init>0
    gen_cxx_model_init_code_inline = @(modelCode,equationCode) gen_cxx_model_code(modelCode,equationCode,...
        0,{},num_policy_init,var_policy_init_name,num_policy_init_total,var_policy_init_length,var_policy_init_loc,...
        var_aux_init_name, var_aux_init_length, num_shock_var,var_shock_name,shock_num,LINE_BREAK);
    [modelAllInitCodes,callFMinInitCodes,num_equations_init] = process_model_code(modelInitCodes,equationInitCodes,modelInitConditionCodes,gen_cxx_model_init_code_inline,...
        headerAuxInitAssignCode,template_folder,USE_FINITE_DIFF,LINE_BREAK);
    
    modelAllInitCodes = my_regexprep(modelAllInitCodes, 'PRE_JAC_CODE', '');
    modelAllInitCodes = my_regexprep(modelAllInitCodes, 'POST_JAC_CODE', '');
else
    modelAllInitCodes='';
    callFMinInitCodes='';
    num_equations_init=0;
end

%% Simulate code
num_data = 1 + sum(parameters_size) + shock_num*shock_num + shock_num*num_shock_var + 1 + num_state + num_tensor;

%% modify output code
output_name_reshape = cell(1,num_output);
for j=1:num_output
    output_name_reshape{j} = ['reshape(',var_output_name{j},',',num2str(max(var_output_length(j),1)),',[])'];
end

%% Read params preset code
paramsPresetCode = fileread([template_folder '/params_template.m']);

%% Attach all code - matlab iter
matlabCode = fileread([template_folder '/iter_template.m']);

% Commonly shared code for solve the problem
if USE_ASG==1
    solveProblemCode = fileread([template_folder '/iter_solve_problem_asg_template.m']);
else
    solveProblemCode = fileread([template_folder '/iter_solve_problem_template.m']);
end

if USE_ASG==1
    solveInitProblemCode = fileread([template_folder '/iter_solve_init_problem_asg_template.m']);
else
    solveInitProblemCode = fileread([template_folder '/iter_solve_problem_template.m']);
end

% Read last period iter template
if USE_ASG==1
    iterInitCode = fileread([template_folder '/iter_init_asg_template.m']);
    iterInitCode = my_regexprep(iterInitCode,'SOLVE_PROBLEM_CODE',solveInitProblemCode);
elseif num_policy_init>0
    iterInitCode = fileread([template_folder '/iter_init_template.m']);
    iterInitCode = my_regexprep(iterInitCode,'SOLVE_PROBLEM_CODE',solveInitProblemCode);
else
    iterInitCode = '';
end

% Read the inifite horizon iter template
outputConstructCode='';
if USE_ASG==1
    iterInfHorizonCode = fileread([template_folder '/iter_inf_horizon_asg_template.m']);
    iterInfHorizonCode = my_regexprep(iterInfHorizonCode,'PRE_ITER_CODE',preIterCode);
    iterInfHorizonCode = my_regexprep(iterInfHorizonCode,'POST_ITER_CODE',postIterCode);
    if num_output>0
        outputConstructCode = fileread([template_folder '/iter_inf_output_construct_asg_template.m']);
    end
    iterInfHorizonCode = my_regexprep(iterInfHorizonCode,'OUTPUT_CONSTRUCT_CODE',outputConstructCode);
    asgProposeGridsAndSolveCode = fileread([template_folder '/asg_propose_grids_and_solve_template.m']);
    iterInfHorizonCode = my_regexprep(iterInfHorizonCode,'ASG_PROPOSE_GRIDS_AND_SOLVE',asgProposeGridsAndSolveCode);
    iterInfHorizonCode = my_regexprep(iterInfHorizonCode,'SOLVE_PROBLEM_CODE',solveProblemCode);
else
    iterInfHorizonCode = fileread([template_folder '/iter_inf_horizon_template.m']);
    iterInfHorizonCode = my_regexprep(iterInfHorizonCode,'SOLVE_PROBLEM_CODE',solveProblemCode);
    iterInfHorizonCode = my_regexprep(iterInfHorizonCode,'PRE_ITER_CODE',preIterCode);
    iterInfHorizonCode = my_regexprep(iterInfHorizonCode,'POST_ITER_CODE',postIterCode);
    if num_output>0
        outputConstructCode = fileread([template_folder '/iter_inf_output_construct_template.m']);
    end
    iterInfHorizonCode = my_regexprep(iterInfHorizonCode,'OUTPUT_CONSTRUCT_CODE',outputConstructCode);
end

% Fill in code segments
iterInitCode = my_regexprep(iterInitCode,'POLICY_INIT_INBOUND_CODE',inboundInitCode);
iterInitCode = my_regexprep(iterInitCode,'POLICY_INIT_INITIALIZE_CODE',policyInitInitializeCode);
iterInitCode = my_regexprep(iterInitCode,'POLICY_INIT_ASSIGN_CODE',policyInitAssignCode);
iterInitCode = my_regexprep(iterInitCode,'AUX_INIT_ASSIGN_CODE',auxInitAssignCode);
iterInitCode = my_regexprep(iterInitCode,'POST_INIT_CODE',postInitCode);
iterInitCode = my_regexprep(iterInitCode,'GDSGE_NUM_AUX_INIT',num2str(num_aux_init_total));
iterInitCode = my_regexprep(iterInitCode,'GDSGE_MAXDIM_INIT',num2str(num_policy_init_total));
iterInitCode = my_regexprep(iterInitCode,'GDSGE_MAXDATA',num2str(num_data));

iterInitCode = my_regexprep(iterInitCode,'TENSOR_CONSTRUCT_CODE',tensorConstructCode);
iterInitCode = my_regexprep(iterInitCode,'TENSOR_ASSIGN_CODE',tensorAssignCode);

iterInitCode = my_regexprep(iterInitCode,'MODEL_NAME',modelName);
iterInitCode = my_regexprep(iterInitCode,'SHOCK_GRID_SEMI_COLON',my_strjoin([var_shock_name ' '],'(:);'));
iterInitCode = my_regexprep(iterInitCode,'TENSOR_SEMI_COLON',my_strjoin([tensor_tensor_name ' '],'(:)'';'));
iterInitCode = my_regexprep(iterInitCode,'PARAMS_SEMI_COLON',my_strjoin(parameters_vec_name,';'));
iterInitCode = my_regexprep(iterInitCode,'STATE_SEMI_COLON',my_strjoin([tensor_state_name ' '],'(:)'';'));

if UseModelId
    iterInitCode = regexprep(iterInitCode,'MODEL_ID',modelName);
else
    iterInitCode = regexprep(iterInitCode,'_MODEL_ID','');
end

% Prepare space
prepareSpaceCode = fileread([template_folder '/iter_inf_prepare_space_template.m']);
prepareSpaceCode = my_regexprep(prepareSpaceCode,'TENSOR_CONSTRUCT_CODE',tensorConstructCode);
prepareSpaceCode = my_regexprep(prepareSpaceCode,'TENSOR_ASSIGN_CODE',tensorAssignCode);
prepareSpaceCode = my_regexprep(prepareSpaceCode,'POLICY_INBOUND_CODE',inboundCode);
prepareSpaceCode = my_regexprep(prepareSpaceCode,'POLICY_INITIALIZE_CODE',policyInitializeCode);

prepareSpaceCode = my_regexprep(prepareSpaceCode,'GDSGE_NUM_AUX',num2str(num_aux_total));
prepareSpaceCode = my_regexprep(prepareSpaceCode,'GDSGE_NUM_INTERP',num2str(num_interp));
prepareSpaceCode = my_regexprep(prepareSpaceCode,'GDSGE_MAXDATA',num2str(num_data));
prepareSpaceCode = my_regexprep(prepareSpaceCode,'GDSGE_MAXDIM',num2str(num_policy_total));

if UseModelId
    prepareSpaceCode = regexprep(prepareSpaceCode,'MODEL_ID',modelName);
else
    prepareSpaceCode = regexprep(prepareSpaceCode,'_MODEL_ID','');
end

% Construct spline
constructSplineCode = fileread([template_folder '/iter_construct_spline_template.m']);
constructSplineCode = my_regexprep(constructSplineCode,'CONSTRUCT_INTERP_CODE',interpConstructCode);

if UseModelId
    constructSplineCode = regexprep(constructSplineCode,'MODEL_ID',modelName);
else
    constructSplineCode = regexprep(constructSplineCode,'_MODEL_ID','');
end

constructSplineCode = my_regexprep(constructSplineCode,'GDSGE_INTERP_ALL',my_strjoin(strcat('GDSGE_PP_',var_interp_name),','));

% Solve and assign
solveAndAssignCode = fileread([template_folder '/iter_inf_solve_and_assign_template.m']);
solveAndAssignCode = my_regexprep(solveAndAssignCode,'SOLVE_PROBLEM_CODE',solveProblemCode);
solveAndAssignCode = my_regexprep(solveAndAssignCode,'POLICY_ASSIGN_CODE',policyAssignCode);
solveAndAssignCode = my_regexprep(solveAndAssignCode,'AUX_ASSIGN_CODE',auxAssignCode);

solveAndAssignCode = my_regexprep(solveAndAssignCode,'MODEL_NAME',modelName);
solveAndAssignCode = my_regexprep(solveAndAssignCode,'SHOCK_GRID_SEMI_COLON',my_strjoin([var_shock_name ' '],'(:);'));
solveAndAssignCode = my_regexprep(solveAndAssignCode,'TENSOR_SEMI_COLON',my_strjoin([tensor_tensor_name ' '],'(:)'';'));
solveAndAssignCode = my_regexprep(solveAndAssignCode,'PARAMS_SEMI_COLON',my_strjoin(parameters_vec_name,';'));
solveAndAssignCode = my_regexprep(solveAndAssignCode,'STATE_SEMI_COLON',my_strjoin([tensor_state_name ' '],'(:)'';'));

if UseModelId
    solveAndAssignCode = regexprep(solveAndAssignCode,'MODEL_ID',modelName);
else
    solveAndAssignCode = regexprep(solveAndAssignCode,'_MODEL_ID','');
end

% Plug into the matlab code
matlabCode = strrep(matlabCode,'POST_INIT_CODE',postInitCode);
matlabCode = strrep(matlabCode,'PRE_INF_CODE',preInfCode);
matlabCode = strrep(matlabCode,'ITER_INIT',iterInitCode);
matlabCode = strrep(matlabCode,'ITER_INF_HORIZON',iterInfHorizonCode);
matlabCode = strrep(matlabCode,'PREPARE_SPACE_CODE',prepareSpaceCode);
matlabCode = strrep(matlabCode,'CONSTRUCT_SPLINE_CODE',constructSplineCode);
matlabCode = strrep(matlabCode,'SOLVE_AND_ASSIGN_CODE',solveAndAssignCode);

% Fill in template for the infinite horizon problem
setParamsCode = [paramsPresetCode LINE_BREAK code3];
matlabCode = my_regexprep(matlabCode,'PRE_CODE',setParamsCode);
matlabCode = my_regexprep(matlabCode,'ASSERT_CODE',assertCode);
matlabCode = my_regexprep(matlabCode,'PRE_MINOR_CODE',preMinorCode);
matlabCode = my_regexprep(matlabCode,'POST_SOL_CODE',postSolCode);
matlabCode = my_regexprep(matlabCode,'PRE_CALL_MEX',preCallMexCode);
matlabCode = my_regexprep(matlabCode,'POST_CALL_MEX',postCallMexCode);

matlabCode = my_regexprep(matlabCode,'TENSOR_CONSTRUCT_CODE',tensorConstructCode);
matlabCode = my_regexprep(matlabCode,'TENSOR_ASSIGN_CODE',tensorAssignCode);

matlabCode = my_regexprep(matlabCode,'POLICY_INBOUND_CODE',inboundCode);
matlabCode = my_regexprep(matlabCode,'POLICY_INITIALIZE_CODE',policyInitializeCode);
matlabCode = my_regexprep(matlabCode,'INTERP_INITIALIZE_CODE',interpInitializeCode);
matlabCode = my_regexprep(matlabCode,'CONSTRUCT_INTERP_CODE',interpConstructCode);

matlabCode = my_regexprep(matlabCode,'POLICY_ASSIGN_CODE',policyAssignCode);
matlabCode = my_regexprep(matlabCode,'INTERP_ASSIGN_CODE',interpAssignCode);
matlabCode = my_regexprep(matlabCode,'INTERP_NEW_ASSIGN_CODE',interpNewAssignCode);
matlabCode = my_regexprep(matlabCode,'INTERP_ASG_ASSIGN_CODE',interpAsgAssignCode);
matlabCode = my_regexprep(matlabCode,'METRIC_CODE',metricCode);
matlabCode = my_regexprep(matlabCode,'RESHAPE_CODE',reshapeCode);
matlabCode = my_regexprep(matlabCode,'RESHAPE_SUB_CODE',subReshapeCode);
matlabCode = my_regexprep(matlabCode,'AUX_ASSIGN_CODE',auxAssignCode);
matlabCode = my_regexprep(matlabCode,'INBOUND_ADAPTIVE_CODE',inboundAdaptiveCode);
matlabCode = my_regexprep(matlabCode,'INBOUND_ADAPTIVE_TIGHT_CODE',inboundAdaptiveTightCode);
matlabCode = my_regexprep(matlabCode,'WARMUP_BOUND_ADAPTIVE_CODE',warmupBoundAdaptiveCode);
matlabCode = my_regexprep(matlabCode,'PRE_UPDATE_CODE',preUpdateCode);
matlabCode = my_regexprep(matlabCode,'OUTPUT_INDEX_ASSIGN_CODE',outputIndexAssignCode);

matlabCode = my_regexprep(matlabCode,'MODEL_NAME',modelName);
matlabCode = my_regexprep(matlabCode,'SHOCK_GRID_SEMI_COLON',my_strjoin([var_shock_name ' '],'(:);'));
matlabCode = my_regexprep(matlabCode,'TENSOR_SEMI_COLON',my_strjoin([tensor_tensor_name ' '],'(:)'';'));
matlabCode = my_regexprep(matlabCode,'PARAMS_SEMI_COLON',my_strjoin(parameters_vec_name,';'));
matlabCode = my_regexprep(matlabCode,'STATE_SEMI_COLON',my_strjoin([tensor_state_name ' '],'(:)'';'));

matlabCode = my_regexprep(matlabCode,'INTERP_VAR_SEMI_COLON',my_strjoin([var_interp_name ' '],';'));
matlabCode = my_regexprep(matlabCode,'STATE_COMMA',my_strjoin(var_state_name,','));

matlabCode = my_regexprep(matlabCode,'GDSGE_NUM_OUTPUT',num2str(num_output_total));
matlabCode = my_regexprep(matlabCode,'OUTPUT_VAR_SEMI_COLON',my_strjoin([var_output_name ' '],';'));
matlabCode = my_regexprep(matlabCode,'OUTPUT_VAR_RESHAPE_COMMA',my_strjoin(output_name_reshape,','));

matlabCode = my_regexprep(matlabCode,'RSLT_PARAMS',my_strjoin([parameters_name 'GDSGE_EMPTY'],','));
matlabCode = my_regexprep(matlabCode,'RSLT_VAR_OTHERS',my_strjoin([var_others_name 'GDSGE_EMPTY'],','));
matlabCode = my_regexprep(matlabCode,'RSLT_SHOCK',my_strjoin(var_shock_name,','));
matlabCode = my_regexprep(matlabCode,'RSLT_STATE',my_strjoin(var_state_name,','));
matlabCode = my_regexprep(matlabCode,'RSLT_POLICY',my_strjoin(var_policy_name,','));

temp = my_strjoin(var_interp_name,',');
if temp
matlabCode = my_regexprep(matlabCode,'RSLT_INTERP',temp);
else
matlabCode = my_regexprep(matlabCode,'RSLT_INTERP','GDSGE_EMPTY');
end

matlabCode = my_regexprep(matlabCode,'RSLT_AUX',my_strjoin([var_aux_name 'GDSGE_EMPTY'],','));
matlabCode = my_regexprep(matlabCode,'RSLT_TENSOR',my_strjoin([var_tensor_name 'GDSGE_EMPTY'],','));
matlabCode = my_regexprep(matlabCode,'RSLT_PP',my_strjoin(pp_interp_name,','));
matlabCode = my_regexprep(matlabCode,'GDSGE_NUM_AUX',num2str(num_aux_total));
matlabCode = my_regexprep(matlabCode,'GDSGE_NUM_INTERP',num2str(num_interp));
matlabCode = my_regexprep(matlabCode,'GDSGE_NUM_SOL',num2str(num_policy_total));
matlabCode = my_regexprep(matlabCode,'GDSGE_MAXDATA',num2str(num_data));
matlabCode = my_regexprep(matlabCode,'GDSGE_MAXDIM',num2str(num_policy_total));

if UseModelId
    matlabCode = regexprep(matlabCode,'MODEL_ID',modelName);
else
    matlabCode = regexprep(matlabCode,'_MODEL_ID','');
end

%% Simulation code
% Read simulateCode
if ~isempty(simulateCode)
    if ~shockInitialized
        error('in simulate: shock is not initalized');
    end
    
    for j=1:num_state
        if ~stateInitialized(j)
            error('in simualte: state var %s is not initalized', var_state_name{j});
        end
        if ~stateAssigned(j)
            error('in simulate: transition for state var %s is not specified', var_state_name{j});
        end
    end
    
    % Construct the space
    simuConstructCode = [];
    for j=1:num_state
        simuConstructCode = [simuConstructCode ...
            'SimuRslt.' var_state_name{j} '=zeros(num_samples,num_periods);' LINE_BREAK];
    end
    for j=1:num_simulate
        simuConstructCode = [simuConstructCode ...
            'SimuRslt.' var_simulate_name{j} '=zeros(num_samples,num_periods);' LINE_BREAK];
    end
    
    % Simu inbound code
    lines = strsplit(inboundCode,{LINE_BREAK,char(13)});
    simuOuterInboundCode = '';
    simuInnerInboundCode = '';
    for j=1:length(lines)
        line = lines{j};
        if regexp(line,'GDSGE_TENSOR_')
            line = regexprep(line,'GDSGE_TENSOR_','');
            simuInnerInboundCode = [simuInnerInboundCode ...
                line LINE_BREAK];
        else
            simuOuterInboundCode = [simuOuterInboundCode ...
                line LINE_BREAK];
        end
    end
    
    lines = strsplit(inboundAdaptiveCode,{LINE_BREAK,char(13)});
    simuInboundAdaptiveCode = '';
    for j=1:length(lines)
        line = lines{j};
        if ~strcmp(line,'')
            words = strsplit(line,'=');
            simuInboundAdaptiveCode = [simuInboundAdaptiveCode ...
                words{1} '=max(IterRslt.GDSGE_PROB.' words{1} ');' LINE_BREAK
                ];
        end
    end
    
    lines = strsplit(initialCode,{LINE_BREAK,char(13)});
    simuOuterInitialCode = '';
    simuInnerInitialCode = '';
    for j=1:length(lines)
        line = lines{j};
        if regexp(line,'GDSGE_TENSOR_')
            line = regexprep(line,'GDSGE_TENSOR_','');
            simuInnerInitialCode = [simuInnerInitialCode ...
                line LINE_BREAK];
        else
            simuOuterInitialCode = [simuOuterInitialCode ...
                line LINE_BREAK];
        end
    end
    
    %% Simu state inbound code
    simuStateInboundCode = '';
    for j=1:num_state
        simuStateInboundCode = [simuStateInboundCode ...
            'SimuRslt.' var_state_name{j} '(:,GDSGE_t)' ...
            '=max(min(IterRslt.var_state.' var_state_name{j} '),SimuRslt.' var_state_name{j} '(:,GDSGE_t));' LINE_BREAK ...
            'SimuRslt.' var_state_name{j} '(:,GDSGE_t)' ...
            '=min(max(IterRslt.var_state.' var_state_name{j} '),SimuRslt.' var_state_name{j} '(:,GDSGE_t));' LINE_BREAK];
    end
    
    %% Simu Pre Assign
    simuPreAssignCode = '';
    for j=1:num_state
        simuPreAssignCode = [simuPreAssignCode ...
            var_state_name{j} '=SimuRslt.' var_state_name{j} '(:,GDSGE_t);' LINE_BREAK];
    end
    %{
    for j=1:num_shock_var
        simuPreAssignCode = [simuPreAssignCode ...
            var_shock_name{j} '=IterRslt.var_shock.' var_shock_name{j} '(shock)'';' LINE_BREAK];
    end
    %}
    
    %% Simu Tensor Code
    simuTensorCode = regexprep(tensorAssignCode,'GDSGE_TENSOR_','');
    
    %% Simu construction and evaluation
    simulateInterpConstructionCode = '';
    simulateInterpEvalCode = '';
    if USE_ASG
        simulateInterpConstructionCode = 'GDSGE_ASG_INTERP = asg.construct_from_struct(IterRslt.asg_output_struct);';
        simulateInterpEvalCode = 'GDSGE_INTERP_RESULTS = GDSGE_ASG_INTERP.eval_vec(SimuRslt.shock(:,GDSGE_t)'',[SIMU_RSLT_STATE_SEMI_COLON]);';
    elseif USE_SPLINE
        simulateInterpConstructionCode = 'GDSGE_PP = IterRslt.output_interp;';
        if shock_num>1
            simulateInterpEvalCode = 'GDSGE_INTERP_RESULTS = myppual_mex(int32(NumThreads),GDSGE_PP.breaks,GDSGE_PP.coefs,int32(GDSGE_PP.pieces),int32(GDSGE_PP.order),int32(GDSGE_PP.dim),''not-a-knot'',[SimuRslt.shock(:,GDSGE_t)'';SIMU_RSLT_STATE_SEMI_COLON],[],[],[]);';
        else
            simulateInterpEvalCode = 'GDSGE_INTERP_RESULTS = myppual_mex(int32(NumThreads),GDSGE_PP.breaks,GDSGE_PP.coefs,int32(GDSGE_PP.pieces),int32(GDSGE_PP.order),int32(GDSGE_PP.dim),''not-a-knot'',[SIMU_RSLT_STATE_SEMI_COLON],[],[],[]);';
        end
    end
    
    %% Simu Post Assign
    simuPostAssignCode = '';
    for j=1:num_simulate
        simuPostAssignCode = [simuPostAssignCode ...
            'SimuRslt.' var_simulate_name{j} '(:,GDSGE_t)=' var_simulate_name{j} ';' LINE_BREAK];
    end
    simuPostAssignCode = [simuPostAssignCode ...
        simulateCode2 LINE_BREAK];
    
    %% state variable in simulate result
    simuRslt_state_name = cell(1,num_state);
    for j=1:num_state
        simuRslt_state_name{j} = ['SimuRslt.' var_state_name{j} '(:,GDSGE_t)'''];
    end
    
    %% Assign simulate results to output variables
    assignInterpToOutputCode = '';
    for j=1:num_output
        output_name_j = var_output_name{j};
        codeLine = [output_name_j '=GDSGE_INTERP_RESULTS(output_var_index.' output_name_j ',:);'];
        assignInterpToOutputCode = [assignInterpToOutputCode, codeLine, LINE_BREAK];
    end
    
    if SIMU_INTERP
        simulateCode = fileread([template_folder '/simulate_interp_template.m']);
        simulateCode = regexprep(simulateCode,'PARAMS_PRESET_CODE',paramsPresetCode);
        simulateCode = regexprep(simulateCode,'MODEL_NAME',modelName);
        simulateCode = regexprep(simulateCode,'SIMU_STATE_INBOUND_CODE',simuStateInboundCode);
        simulateCode = regexprep(simulateCode,'SIMU_PRE_CODE',[code3 simulateCode1]);
        simulateCode = regexprep(simulateCode,'SIMU_CONSTRUCT_CODE',simuConstructCode);
        simulateCode = regexprep(simulateCode,'SIMULATE_INIT_CODE',simuInitCode);
        simulateCode = regexprep(simulateCode,'SIMULATE_INIT_OVERWRITE_CODE',simuInitOverwriteCode);
        simulateCode = regexprep(simulateCode,'SIMULATE_INTERP_CONSTRUCT_CODE',simulateInterpConstructionCode);
        simulateCode = regexprep(simulateCode,'SIMULATE_INTERP_EVAL_CODE',simulateInterpEvalCode);
        simulateCode = regexprep(simulateCode,'SIMU_POST_ASSIGN_CODE',simuPostAssignCode);
        simulateCode = regexprep(simulateCode,'SIMU_PRE_ITER_CODE',simuPreIterCode);
        simulateCode = regexprep(simulateCode,'SIMU_POST_ITER_CODE',simuPostIterCode);
        simulateCode = regexprep(simulateCode,'SIMU_RSLT_STATE_SEMI_COLON',my_strjoin(simuRslt_state_name,';'));
        simulateCode = regexprep(simulateCode,'ASSIGN_INTERP_TO_OUTPUT_CODE',assignInterpToOutputCode);
    elseif SIMU_RESOLVE
        if USE_ASG
            simulateCode = fileread([template_folder '/simulate_resolve_asg_template.m']);
        else
            simulateCode = fileread([template_folder '/simulate_resolve_template.m']);
        end
        simulateCode = regexprep(simulateCode,'SIMU_RSLT_STATE_SEMI_COLON',my_strjoin(simuRslt_state_name,';'));
        simulateCode = regexprep(simulateCode,'PARAMS_PRESET_CODE',paramsPresetCode);
        simulateCode = regexprep(simulateCode,'MODEL_NAME',modelName);
        simulateCode = regexprep(simulateCode,'SIMU_PRE_CODE',[code3 simulateCode1]);
        simulateCode = regexprep(simulateCode,'SIMU_CONSTRUCT_CODE',simuConstructCode);
        simulateCode = regexprep(simulateCode,'POLICY_INBOUND_CODE',simuOuterInboundCode);
        simulateCode = regexprep(simulateCode,'POLICY_INIT_CODE',simuOuterInitialCode);
        simulateCode = regexprep(simulateCode,'SIMULATE_INBOUND_ADAPTIVE_CODE',simuInboundAdaptiveCode);
        simulateCode = regexprep(simulateCode,'SIMULATE_INIT_CODE',simuInitCode);
        simulateCode = regexprep(simulateCode,'SIMULATE_INIT_OVERWRITE_CODE',simuInitOverwriteCode);
        simulateCode = regexprep(simulateCode,'PARAMS_SEMI_COLON',my_strjoin(parameters_vec_name,';'));
        simulateCode = regexprep(simulateCode,'SHOCK_SEMI_COLON',my_strjoin([var_shock_name ' '],'(:);'));
        simulateCode = regexprep(simulateCode,'SIMU_PRE_ASSIGN_CODE',simuPreAssignCode);
        simulateCode = regexprep(simulateCode,'SIMU_TENSOR',simuTensorCode);
        simulateCode = regexprep(simulateCode,'SIMU_INBOUND_CODE',simuInnerInboundCode);
        simulateCode = regexprep(simulateCode,'SIMU_INIT_CODE',simuInnerInitialCode);
        simulateCode = regexprep(simulateCode,'POLICY_ASSIGN_CODE',policyAssignCode);
        simulateCode = regexprep(simulateCode,'AUX_ASSIGN_CODE',auxAssignCode);
        simulateCode = regexprep(simulateCode,'INBOUND_ADAPTIVE_CODE',inboundAdaptiveCode);
        simulateCode = regexprep(simulateCode,'INBOUND_ADAPTIVE_TIGHT_CODE',inboundAdaptiveTightCode);
        simulateCode = regexprep(simulateCode,'SIMU_POST_ASSIGN_CODE',simuPostAssignCode);
        simulateCode = regexprep(simulateCode,'SIMU_PRE_ITER_CODE',simuPreIterCode);
        simulateCode = regexprep(simulateCode,'SIMU_POST_ITER_CODE',simuPostIterCode);
        
        simulateCode = regexprep(simulateCode,'STATE_COMMA',my_strjoin(var_state_name,','));
        simulateCode = regexprep(simulateCode,'STATE_SEMI_COLON',my_strjoin([var_state_name ' '],'(:)'';'));
        simulateCode = regexprep(simulateCode,'TENSOR_SEMI_COLON',my_strjoin([var_tensor_name ' '],'(:)'';'));
        simulateCode = regexprep(simulateCode,'GDSGE_NUM_AUX',num2str(num_aux));
        simulateCode = regexprep(simulateCode,'GDSGE_MAXDIM',num2str(num_policy_total));
        simulateCode = regexprep(simulateCode,'GDSGE_MAXDATA',num2str(num_data));
        
        if UseModelId
            simulateCode = regexprep(simulateCode,'MODEL_ID',modelName);
        else
            simulateCode = regexprep(simulateCode,'_MODEL_ID','');
        end
    end
else
    simulateCode = '';
end

%% Attach all code together - header file
% Start with template
cppCode = fileread([template_folder '/mex_template.cpp']);
cppCode = regexprep(cppCode,'GDSGE_OTHER_INCLUDE',cxxIncludeString);
% init
taskInitCode = fileread([template_folder '/task_template.cpp']);
taskInitCode = my_regexprep(taskInitCode,'PRE_MODEL_CODE','');
taskInitCode = my_regexprep(taskInitCode,'START_LOOP_CODE','');
taskInitCode = my_regexprep(taskInitCode,'FINISH_LOOP_CODE','');
taskInitCode = my_regexprep(taskInitCode,'MODEL_CODE',modelAllInitCodes);
taskInitCode = my_regexprep(taskInitCode,'CALL_FMIN_CODE;',callFMinInitCodes);
taskInitCode = my_regexprep(taskInitCode,'NUM_EQUATIONS',num2str(num_equations_init));
taskInitCode = my_regexprep(taskInitCode,'NUM_AUX',num2str(num_aux_init_total));
taskInitCode = my_regexprep(taskInitCode,'INTERP_GET_CODE','');
taskInitCode = my_regexprep(taskInitCode,'INTERP_GET_THREAD_CODE','');
taskInitCode = my_regexprep(taskInitCode,'TASK_NAME','task_init');
cppCode = regexprep(cppCode,'TASK_INIT_CODE',taskInitCode);

% infinite horizon
% init
taskInfHorizonCode = fileread([template_folder '/task_template.cpp']);
taskInfHorizonCode = my_regexprep(taskInfHorizonCode,'PRE_MODEL_CODE',preModelAllCode);
taskInfHorizonCode = my_regexprep(taskInfHorizonCode,'START_LOOP_CODE',startLoopCode);
taskInfHorizonCode = my_regexprep(taskInfHorizonCode,'FINISH_LOOP_CODE',finishLoopCode);
taskInfHorizonCode = my_regexprep(taskInfHorizonCode,'MODEL_CODE',modelAllCodes);
taskInfHorizonCode = my_regexprep(taskInfHorizonCode,'CALL_FMIN_CODE;',callFMinCodes);
taskInfHorizonCode = my_regexprep(taskInfHorizonCode,'NUM_EQUATIONS',num2str(num_equations));
taskInfHorizonCode = my_regexprep(taskInfHorizonCode,'NUM_AUX',num2str(num_aux_total));
taskInfHorizonCode = my_regexprep(taskInfHorizonCode,'INTERP_GET_CODE',interpGetCode);
taskInfHorizonCode = my_regexprep(taskInfHorizonCode,'INTERP_GET_THREAD_CODE',interpGetThreadCode);
taskInfHorizonCode = my_regexprep(taskInfHorizonCode,'TASK_NAME','task_inf_horizon');
cppCode = regexprep(cppCode,'TASK_INF_HORIZON_CODE',taskInfHorizonCode);

cppCode = regexprep(cppCode,'POP_CODE',popCode);
cppCode = regexprep(cppCode,'GDSGE_INTERP_XDIM',num2str(num_state));


%{
% print header to file
fileID = fopen(['mex_' modelName '.cpp'],'w');
fprintf(fileID,'%s',headerCode);
fclose(fileID);
%}

%% Compilation code
% Get the folder of gdsge
include_folder = [gdsge_folder '/include'];
% Calculate some extra definitions
extraDef = '';
% First is if has model_init
if num_policy_init>0
    extraDef = [extraDef  ' -DHAS_INIT'];
end
if USE_ASG
    % Get the max size
    [ASG_MAX_DIM,ASG_MAX_NVEC,ASG_MAX_LEVEL] = asg.get_mex_constants();
    if num_state>ASG_MAX_DIM
        error('num of states excceeds MAX_DIM of ASG.');
    end
    if AsgMaxLevel>ASG_MAX_LEVEL
        error('asg level excceeds MAX_LEVEL.');
    end
    
    extraDef = [extraDef ' -DUSE_ASG' ' -DASG_MAX_DIM=' num2str(ASG_MAX_DIM) ' -DASG_MAX_NVEC=' num2str(ASG_MAX_NVEC) ' -DASG_MAX_LEVEL=' num2str(ASG_MAX_LEVEL) ];
end
if REMOVE_NULL_STATEMENTS
    extraDef = [extraDef ' -DADEPT_REMOVE_NULL_STATEMENTS'];
end
if USE_SPARSE_JACOBIAN
    extraDef = [extraDef ' -DUSE_SPARSE_JACOBIAN'];
end
if USE_FINITE_DIFF
    extraDef = [extraDef ' -DUSE_FINITE_DIFF'];
end

compileCode = fileread([template_folder '/compile_template.m']);
compileCode = strrep(compileCode,'EXTRA_DEF',extraDef);
compileCode = strrep(compileCode,'GDSGE_FOLDER',gdsge_folder);
compileCode = strrep(compileCode,'MODEL_NAME',modelName);
compileCode = strrep(compileCode,'INCLUDE_FOLDER',include_folder);
compileCode = strrep(compileCode,'GDSGE_MAXDIM',num2str(max(num_policy_total,num_policy_init_total)+4));
compileCode = strrep(compileCode,'GDSGE_INTERP_ORDER',num2str(INTERP_ORDER));

model = v2struct(parameters_name,var_state_name,var_shock_name,var_tensor_name,var_interp_name,var_aux_name,shock_num,shock_trans,params);
codeSegment = v2struct(setParamsCode,iterInitCode,prepareSpaceCode,constructSplineCode,solveAndAssignCode);
end

function codeNew = my_regexprep(codeOld,keyword,filledCode)
allButUnderscore = '[^a-zA-Z0-9]';

BEFORE_A_WORD = ['((?<=^)|(?<=' allButUnderscore '))'];
AFTER_A_WORD = ['((?=$)|(?=' allButUnderscore '))'];

codeNew = regexprep(codeOld,[BEFORE_A_WORD keyword AFTER_A_WORD],filledCode);
end

function codeNew = my_regexp(codeOld,keyword)
allButUnderscore = '[^a-zA-Z0-9]';

BEFORE_A_WORD = ['((?<=^)|(?<=' allButUnderscore '))'];
AFTER_A_WORD = ['((?=$)|(?=' allButUnderscore '))'];

codeNew = regexp(codeOld,[BEFORE_A_WORD keyword AFTER_A_WORD],'all');
end

function [modelAllCodes,callFMinCodes,num_equations] = process_model_code(modelCodes,equationCodes,modelConditionCodes,gen_cxx_model_code_inline,...
    headerAuxAssignCode,template_folder,USE_FINITE_DIFF,LINE_BREAK)
modelTemplate = fileread([template_folder '/model_template.cpp']);
modelAllCodes = '';
modelEvalTemplate = fileread([template_folder '/model_eval_template.cpp']);
modelEvalCode = modelEvalTemplate;

callFMinTemplate = fileread([template_folder '/call_fmin_template.cpp']);
callFMinCodes = '';
for j=1:length(modelCodes)
    modelAllCode = modelTemplate;
    [modelCode,equationCode,argCode,declareCode,num_equations,var_interp_name] = gen_cxx_model_code_inline(modelCodes{j},equationCodes{j});
    var_interp_name = [var_interp_name,'GDSGE_INTERP_VEC','GDSGE_INTERP_RSLT'];
    for i=1:length(var_interp_name)
        modelCode = regexprep(modelCode,['((?<=^)|(?<=\W)|(?<=_))', var_interp_name{i} '(?=\W)'],[var_interp_name{i},'_adouble']);
    end
    
    modelAllCode = regexprep(modelAllCode,'ARG_CODE',argCode);
    modelAllCode = regexprep(modelAllCode,'DECLARE_CODE',declareCode);
    modelAllCode = regexprep(modelAllCode,'MODEL_CODE',modelCode);
    modelAllCode = regexprep(modelAllCode,'EQUATION_CODE',equationCode);
    modelAllCode = regexprep(modelAllCode,'HEADER_AUX_ASSIGN_CODE',headerAuxAssignCode);
    modelAllCodeDouble = regexprep(modelAllCode,'((?<=^)|(?<=\W)|(?<=_))adouble(?=\W)','double');
    if USE_FINITE_DIFF
        modelAllCode = modelTemplate;
        modelAllCode = regexprep(modelAllCode,'ARG_CODE','');
        modelAllCode = regexprep(modelAllCode,'DECLARE_CODE','');
        modelAllCode = regexprep(modelAllCode,'MODEL_CODE','');
        modelAllCode = regexprep(modelAllCode,'EQUATION_CODE','');
        modelAllCode = regexprep(modelAllCode,'HEADER_AUX_ASSIGN_CODE','');
        modelAllCode = [modelAllCode,LINE_BREAK,modelAllCodeDouble];
    else
        modelAllCode = [modelAllCode,LINE_BREAK,modelAllCodeDouble];
    end
    modelAllCode = [modelAllCode,LINE_BREAK,modelEvalCode];
    modelAllCode = regexprep(modelAllCode,'NUM_EQUATIONS',num2str(num_equations));
    modelAllCode = regexprep(modelAllCode,'MODEL_NUMBER',num2str(j));
    
    modelAllCodes = [modelAllCodes LINE_BREAK modelAllCode];
    
    callFMinCode = callFMinTemplate;
    callFMinCode = regexprep(callFMinCode,'MODEL_NUMBER',num2str(j));
    callFMinCode = regexprep(callFMinCode,'NUM_EQUATIONS',num2str(num_equations));
    modelCondition = modelConditionCodes{j};
    if strcmp(modelCondition,'')
        modelCondition = '1';
    end
    callFMinCode = regexprep(callFMinCode,'MODEL_CONDITION',modelCondition);
    callFMinCodes = [callFMinCodes LINE_BREAK callFMinCode];
end
end

function word = replace_tensor_name(word,var_pre_tensor_name,tensorPrefix)
% Replace tensor name with TENSOR_ prefix
for k=1:length(var_pre_tensor_name)
    word = regexprep(word,['((?<=^)|(?<=\W))' var_pre_tensor_name{k} '((?=$)|(?=\W))'],[tensorPrefix var_pre_tensor_name{k} '(:)'''],'all');
end
end

function assert_valid_var_name(var_all_name)
% Check duplicates and reserve words in the var name
[unique_var_name,~,var_count_bin] = unique(var_all_name,'stable');
var_all_name_count = histc(var_count_bin,1:numel(unique_var_name));
duplicate_var = find(var_all_name_count>1);
if (duplicate_var)
    error('GDSGE:duplicateVar','Duplicate var found: %s.',unique_var_name{duplicate_var});
end
check_reserve_word(unique_var_name);
end

function gridCode = replace_grid_variable(gridCode,num_interp,var_interp_name,num_policy,var_policy_name,num_shock_var,var_shock_name,f_var_name,fVarPrefix)
% Replace cFuture'(Xp) to cFuture(GDSGE_iter,Xp)
for k=1:num_interp
    gridCode = regexprep(gridCode,['((?<=^)|(?<=\W))' var_interp_name{k} '''('], [var_interp_name{k} '(GDSGE_iter,'],'all');
end

%{
% Replace r' to r_GDSGE_GRID[GDSGE_iter]
% The reason for keeping r_GDSGE_GRID[GDSGE_iter] notation is because matlab parser
% cannot recognize r(GDSGE_iter) as a left hand side variable
gridCode = regexprep(gridCode,'''',[fVarPrefix '[GDSGE_iter]']);
%}
% Replace r' to r_GDSGE_GRID_GDSGE_iter
% The commented out block is obsoleted by sym2str
% so we have to do a repladcement after calling sym2str
gridCode = regexprep(gridCode,'''',[fVarPrefix '_GDSGE_iter']);

% Replace Xp_GDSGE_GRID[GDSGE_iter] to Xp_GRID(GDSGE_iter)
for k=1:num_policy
    gridCode = regexprep(gridCode,['((?<=^)|(?<=\W))' var_policy_name{k} fVarPrefix '\[GDSGE_iter\]'],[var_policy_name{k} '_GRID(GDSGE_iter)']);
end

% Replace Xp(xx) to Xp_GRID(xx)
for k=1:num_policy
    gridCode = regexprep(gridCode,['((?<=^)|(?<=\W))' var_policy_name{k} '(?=\()'],[var_policy_name{k} '_GRID']);
end

% Replace A(xx) to A_GRID(xx)
for k=1:num_shock_var
    gridCode = regexprep(gridCode,['((?<=^)|(?<=\W))' var_shock_name{k} '(?=\()'],[var_shock_name{k} '_GRID']);
end

% Replace qNext(xx) to qNext_GRID(xx)
for k=1:length(f_var_name)
    gridCode = regexprep(gridCode,['((?<=^)|(?<=\W))' f_var_name{k} '(?=\()'],[f_var_name{k} '_GRID']);
end

end

function [reductionCode,reductionCounter,seg] = extract_reduction_operator(seg,operator,var_interp_name,num_interp,var_policy_name,num_policy,num_shock_var,var_shock_name,f_var_name,fVarPrefix,reductionCounter,shock_num,LINE_BREAK)
reductionCode = '';
parseComment = '';
% Search for GDSGE_EXPECT
loc_start = regexp(seg,[operator '{'],'all');
% Search for end
for loc_end=loc_start:length(seg)
    if seg(loc_end) == '}'
        break;
    end
end
if loc_end>length(seg)
    error(['Unmatched {} in ' operator]);
end

while (loc_start)
    loopCode = seg(loc_start:loc_end);
    [~,exp_start,exp_end] = regexp(loopCode,['(?<=' operator '{).*(?=})'],'tokenExtents');
    loopCode = loopCode(exp_start:exp_end);
    loopCodeArg = strsplit(loopCode,'|');
    loopCode = loopCodeArg{1};
    
    loopCode = replace_grid_variable(loopCode,num_interp,var_interp_name,num_policy,var_policy_name,num_shock_var,var_shock_name,f_var_name,fVarPrefix);
    
    if length(loopCodeArg)>1
        transName = loopCodeArg{2};
    else
        transName = 'shock_trans';
    end
    
    seg = [seg(1:loc_start-1) 'GDSGE_REDUCTION_VALUE' num2str(reductionCounter) seg(loc_end+1:end)];
    
    switch operator
        case 'GDSGE_EXPECT'
            loopCode = ['GDSGE_REDUCTION_VALUE' num2str(reductionCounter) '=' 'GDSGE_REDUCTION_VALUE' num2str(reductionCounter) '+' transName '((shock)+' num2str(shock_num) '*(GDSGE_iter-1)) * (' loopCode ');'];
            parsed_code = my_ccode(my_sym(loopCode));
            parsed_code = strrep(parsed_code, '_GDSGE_iter', '[GDSGE_iter-1]');
            reductionCode = [reductionCode ...
                'adouble GDSGE_REDUCTION_VALUE' num2str(reductionCounter) '=0;' LINE_BREAK ...
                'for(int GDSGE_iter=1; GDSGE_iter<=' num2str(shock_num) '; GDSGE_iter++)' ...
                LINE_BREAK '{' LINE_BREAK ...
                parsed_code parseComment LINE_BREAK ...
                '}'
                ];
        case 'GDSGE_MIN'
            loopCode = ['GDSGE_REDUCTION_VALUE' num2str(reductionCounter) '=' 'MIN(GDSGE_REDUCTION_VALUE' num2str(reductionCounter) ',' loopCode ');'];
            parsed_code = my_ccode(my_sym(loopCode));
            parsed_code = strrep(parsed_code, '_GDSGE_iter', '[GDSGE_iter-1]');
            reductionCode = [reductionCode ...
                'adouble GDSGE_REDUCTION_VALUE' num2str(reductionCounter) '=1e20;' LINE_BREAK ...
                'for(int GDSGE_iter=1; GDSGE_iter<=' num2str(shock_num) '; GDSGE_iter++)' ...
                LINE_BREAK '{' LINE_BREAK ...
                parsed_code parseComment LINE_BREAK ...
                '}'
                ];
        case 'GDSGE_MAX'
            loopCode = ['GDSGE_REDUCTION_VALUE' num2str(reductionCounter) '=' 'MAX(GDSGE_REDUCTION_VALUE' num2str(reductionCounter) ',' loopCode ');'];
            parsed_code = my_ccode(my_sym(loopCode));
            parsed_code = strrep(parsed_code, '_GDSGE_iter', '[GDSGE_iter-1]');
            reductionCode = [reductionCode ...
                'adouble GDSGE_REDUCTION_VALUE' num2str(reductionCounter) '=-1e20;' LINE_BREAK ...
                'for(int GDSGE_iter=1; GDSGE_iter<=' num2str(shock_num) '; GDSGE_iter++)' ...
                LINE_BREAK '{' LINE_BREAK ...
                parsed_code parseComment LINE_BREAK ...
                '}'
                ];
        case 'GDSGE_PROD'
            loopCode = ['GDSGE_REDUCTION_VALUE' num2str(reductionCounter) '=' 'GDSGE_REDUCTION_VALUE' num2str(reductionCounter) '*' loopCode ';'];
            parsed_code = my_ccode(my_sym(loopCode));
            parsed_code = strrep(parsed_code, '_GDSGE_iter', '[GDSGE_iter-1]');
            reductionCode = [reductionCode ...
                'adouble GDSGE_REDUCTION_VALUE' num2str(reductionCounter) '=1;' LINE_BREAK ...
                'for(int GDSGE_iter=1; GDSGE_iter<=' num2str(shock_num) '; GDSGE_iter++)' ...
                LINE_BREAK '{' LINE_BREAK ...
                parsed_code parseComment LINE_BREAK ...
                '}'
                ];
    end
    
    reductionCounter = reductionCounter+1;
    loc_start = regexp(seg,operator,'all');
    % Search for end
    for loc_end=loc_start:length(seg)
        if seg(loc_end) == '}'
            break;
        end
    end
    if loc_end>length(seg)
        error(['Unmatched {} in ' operator]);
    end
end
end

function [modelCodes, modelConditions, equationCodes, remainingCode] = extract_model_segment(code,modelTypeName,assertNotEmpty)
searchString = [modelTypeName '(\(.*?\))?;.*?equations;.*?end;.*?end;'];
[blockStart,blockEnd] = regexp(code,searchString);
if assertNotEmpty
    if isempty(blockStart)
        error('GDSGE:noModel','No model segment found');
    end
end

modelCodes = {};
modelConditions = {};
equationCodes = {};

while(blockStart)
    % Deal with each model segment
    seg = code(blockStart(1):blockEnd(1));
    code(blockStart(1):blockEnd(1))=[];
    
    % Extract model conditions
    [findStart,findEnd] = regexp(seg,'model(\(.*?\))?;');
    conditionSeg = seg(findStart:findEnd);
    modelCondition = regexp(conditionSeg,'(?<=model\().*(?=\);)','match');
    if isempty(modelCondition)
        modelCondition = {''};
    end
    modelConditions = [modelConditions,modelCondition];
    
    % Remove model header
    seg(findStart:findEnd) = [];
    % Remove model ender
    seg(end-3:end) = [];
    
    % Extract equation
    [findStart,findEnd] = regexp(seg,'equations;.*end;');
    equationSeg = seg(findStart:findEnd);
    seg(findStart:findEnd) = [];
    
    equationCode = regexp(equationSeg,'(?<=equations;).*(?=end;)','match');
    equationCodes = [equationCodes,equationCode];
    
    modelCodes = [modelCodes seg];
    
    [blockStart,blockEnd] = regexp(code,'model\(.*?\);.*?equations;.*?end;.*?end;');
end
remainingCode = code;
end

function [modelCode,equationCode,argCode,declareCode,num_equations,var_interp_name] = gen_cxx_model_code(modelCode,equationCode,num_interp,var_interp_name,num_policy,var_policy_name,num_total_policy,var_policy_length,var_policy_loc,...
    var_aux_name,var_aux_length,num_shock_var,var_shock_name,shock_num,LINE_BREAK)
%% Process model code
% Split lines and seg
lines = strsplit(modelCode,{char(10),char(13)});
c_var_name = {};
f_var_name = {};
vec_var_name = {};
vec_var_length = [];
fVarPrefix = '_GDSGE_GRID';
modelCode = '';
warning('off','symbolic:generate:FunctionNotVerifiedToBeValid');
reductionCounter = 0;
j=0;
while j<length(lines)
    j=j+1;
    line = lines{j};
    line0 = line;
    parseComment = [' //parsed from gmod Line: ' line0];
    %{
    try
    %}
        % Remove comments
        [~,loc_start,loc_end] = regexp(line,'%.*$','tokenExtents','all');
        line(loc_start:loc_end)=[];
        
        line = strtrim(line);
        words = strsplit(line,{' '});
        
        if strcmp(words{1},'cxx;')
            % Begin a cxx block
            % Find the next end or cxxend;
            line = '';
            words = strsplit(line,{' '});
            firstWord = words{1};
            while ~strcmp(firstWord,'endcxx;')
                modelCode = [modelCode, line, LINE_BREAK];
                j = j+1;
                line = lines{j};
                line = strtrim(line);
                words = strsplit(line,{' '});
                firstWord = words{1};
            end
            continue;
        end
        
        if strcmp(words{1},'vector')
            % Declare a vector
            if length(words)<2
                error('no vector declaration');
            end
            for i_word = 2:length(words)
                currentWord = words{i_word};
                % Split name and size information
                pat = '(?<varname>\w+)\[(?<varsize>.+)\]';
                patfind = regexp(currentWord,pat,'names');
                if isempty(patfind)
                    error('vector declararion format error, use varname[varlength]');
                end
                vec_var_name{end+1} = patfind(1).varname;
                vec_var_length{end+1} = patfind(1).varsize;
            end
        end
        
        if strcmp(words{1},'for')
            % Begins a for block
            words = strsplit(line,{' ','='});
            % Get counter of the for block
            counter_var_name = words{2};
            range_name = words{3};
            
            % Split range
            range_name_words = strsplit(range_name,':');
            if length(range_name_words)==2
                % This is a step 1 range
                step = '1';
                rangeStart = range_name_words{1};
                rangeEnd = range_name_words{2};
            elseif length(range_name_words)==3
                % This is a step n range
                step = range_name_words{2};
                rangeStart = range_name_words{1};
                rangeEnd = range_name_words{3};
            else
                error('Range not specified correctly in for block');
            end
            
            % Write the for block header
            blockCode = ['for (int ' counter_var_name '=' rangeStart ';' counter_var_name '<=' rangeEnd ';' ...
                counter_var_name '=' counter_var_name '+' step ')' LINE_BREAK ...
                '{' LINE_BREAK];
            modelCode = [modelCode ...
                blockCode];
            continue;
        end
        
        if strcmp(words{1},'if')
            % This is a if block
            % Get conditions of the if block
            condition_idx = 1;
            conditions = {};
            eq_assignments = {};
            
            conditions{condition_idx} = strcat(words{2:end});
            
            conditions{condition_idx} = replace_grid_variable(conditions{condition_idx},num_interp,var_interp_name,num_policy,var_policy_name,num_shock_var,var_shock_name,f_var_name,fVarPrefix);
            
            % Find else and else if block
            found_else = false;
            in_if_block = true;
            
            blockCode = ['if (' conditions{condition_idx} ')' LINE_BREAK '{' LINE_BREAK];
            modelCode = [modelCode ...
                blockCode];
            continue;
        end
        
        if strcmp(words{1},'else')
            % else block found
            found_else = true;
            
            condition_idx = condition_idx+1;
            conditions{condition_idx} = '1';
            
            blockCode = ['}' LINE_BREAK 'else if (' conditions{condition_idx} ')' LINE_BREAK '{' LINE_BREAK];
            modelCode = [modelCode ...
                blockCode];
            continue;
        end
        
        if strcmp(words{1},'elseif')
            % elseif block found
            condition_idx = condition_idx+1;
            conditions{condition_idx} = strcat(words{2:end});
            
            conditions{condition_idx} = replace_grid_variable(conditions{condition_idx},num_interp,var_interp_name,num_policy,var_policy_name,num_shock_var,var_shock_name,f_var_name,fVarPrefix);
            
            blockCode = ['}' LINE_BREAK 'else if (' conditions{condition_idx} ')' LINE_BREAK '{' LINE_BREAK];
            modelCode = [modelCode ...
                blockCode];
            continue;
        end
        
        if strcmp(words{1},'end')
            blockCode = ['}' LINE_BREAK];
            modelCode = [modelCode ...
                blockCode];
            
            continue;
        end
        
        % If it's an assignment
        if contains(line,'=')
            % Split by space
            words = strsplit(line,{' ','='});
            
            % deal with GDSGE_VEC operator
            if contains(line,'GDSGE_INTERP_VEC')
                % split to lhs and rhs
                splitsByEqual = strsplit(line,{'='});
                if ~(length(splitsByEqual)==2)
                    error('GDSGE_INTERP_VEC needs to return values');
                end
                lhs = strtrim(splitsByEqual{1});
                rhs = strtrim(splitsByEqual{2});
                %
                % Each entry on lhs is a future variable
                lhs = regexp(lhs,'(?<=[).*(?=])','match');
                if isempty(lhs)
                    error('GDSGE_INTERP_VEC lhs error');
                end
                lhsSplitsByComma = strsplit(lhs{1},{','});
                if isempty(lhsSplitsByComma)
                    warning('GDSGE_INTERP_VEC no assignments');
                end
                
                % Find optional arguments in a bracket
                leftBracketFound = false;
                rightBracketFound = false;
                for i_char = 1:length(rhs)
                    if rhs(i_char)=='['
                        % Left bracket found
                        leftBracketPosition = i_char;
                        leftBracketFound = true;
                        break;
                    end
                end
                if leftBracketFound
                    for i_char=leftBracketPosition:length(rhs)
                        if rhs(i_char)==']'
                            %
                            rightBracketPosition = i_char;
                            rightBracketFound = true;
                            break;
                        end
                    end
                    if ~rightBracketFound
                        error('Matched bracket not found in GDSGE_INTERP_VEC');
                    end
                    
                    % Generate index based on contents in the bracket
                    bracketContents = rhs(leftBracketPosition:rightBracketPosition);
                    bracketRange = eval(bracketContents);
                    % Push bracket range into evaluation index
                    % Reset to zero
                    modelCode = [modelCode 'memset(INTERP_VEC_FLAG,0,sizeof(bool)*' num2str(num_interp) ');' LINE_BREAK];
                    for i=1:length(bracketRange)
                        modelCode = [modelCode 'INTERP_VEC_FLAG[' num2str(bracketRange(i)-1) ']=true;' LINE_BREAK];
                    end
                    
                    % Remove the bracket options
                    rhs(leftBracketPosition:rightBracketPosition) = [];
                else % if leftBracket Found
                    % left bracket not found
                    modelCode = [modelCode 'memset(INTERP_VEC_FLAG,true,sizeof(bool)*' num2str(num_interp) ');' LINE_BREAK];
                end
                
                % if this evaluates future
                if contains(rhs,'GDSGE_INTERP_VEC''')
                    loopCodes = '';
                    
                    % The evaluation loop
                    % Prepare a GDSGE_INTERP_VEC as a temp interp
                    var_interp_name_extended = [var_interp_name,'GDSGE_INTERP_VEC'];
                    num_interp_extended = num_interp+1;
                    loopCode = replace_grid_variable(rhs,num_interp_extended,var_interp_name_extended,num_policy,var_policy_name,num_shock_var,var_shock_name,f_var_name,fVarPrefix);
                    loopCode = my_ccode(my_sym(loopCode));
                    % PATCHED for str2sym in R2024a
                    loopCode = strrep(loopCode, '_GDSGE_iter', '[GDSGE_iter-1]');
                    % Remove the equal bracket
                    loopCode = loopCode(8:end);
                    loopCode = [loopCode, parseComment];
                    loopCodes = [loopCodes,loopCode,LINE_BREAK];
                    
                    for i_lhs = 1:length(lhsSplitsByComma)
                        word = lhsSplitsByComma{i_lhs};
                        if ~contains(word,'''')
                            error('GDSGE_INTERP_VEC lhs should agree with rhs prime');
                        end
                        word = regexprep(word,'''','');
                        f_var_name{end+1} = word;
                        % loopCode = [word '_GDSGE_GRID[GDSGE_iter]=GDSGE_INTERP_RSLT[' num2str(i_lhs) ']'];
                        loopCode = [word '_GDSGE_GRID_GDSGE_iter=GDSGE_INTERP_RSLT_GDSGE_bracket'];
                        loopCode = my_ccode(my_sym(loopCode));
                        loopCode = strrep(loopCode, '_GDSGE_iter', '[GDSGE_iter-1]');
                        loopCode = strrep(loopCode, '_GDSGE_bracket', ['[' num2str(i_lhs-1) ']']);
                        loopCode = [loopCode, parseComment];
                        loopCodes = [loopCodes,LINE_BREAK,...
                            loopCode];
                    end
                    
                    % Replace future to a loop
                    forLoopCode = [
                        'for(int GDSGE_iter=1; GDSGE_iter<=' ...
                        num2str(shock_num) '; GDSGE_iter++)' ...
                        LINE_BREAK '{' LINE_BREAK...
                        loopCodes LINE_BREAK...
                        '}' LINE_BREAK ...
                        ];
                    modelCode = [modelCode ...
                        forLoopCode];
                else
                    % Check the first word is valid
                    rhsWords = strsplit(rhs,{',','(',' '});
                    numberCharList = '1234567890';
                    if ~(strcmp(rhsWords{2},'shock') || all(ismember(rhsWords{2}, numberCharList)))
                        error('GDSGE_INTERP_VEC without prime must take first argument as "shock" or an integer');
                    end
                    loopCodes = '';
                    loopCode = my_ccode(my_sym(rhs));
                    % Remove the equal bracket
                    loopCode = loopCode(8:end);
                    loopCodes = [loopCodes,loopCode,LINE_BREAK];
                    for i_lhs = 1:length(lhsSplitsByComma)
                        word = lhsSplitsByComma{i_lhs};
                        if contains(word,'''')
                            error('GDSGE_INTERP_VEC lhs does not have prime but rhs has prime');
                        end
                        word = regexprep(word,'''','');
                        c_var_name{end+1} = word;
                        % loopCode = [word '=GDSGE_INTERP_RSLT[' num2str(i_lhs) ']'];
                        loopCode = [word '=GDSGE_INTERP_RSLT_GDSGE_bracket'];
                        loopCode = my_ccode(my_sym(loopCode));
                        loopCode = strrep(loopCode, '_GDSGE_bracket', ['[' num2str(i_lhs-1) ']']);
                        loopCode = [loopCode, parseComment];
                        loopCodes = [loopCodes,LINE_BREAK,...
                            loopCode];
                    end
                    modelCode = [modelCode ...
                        loopCodes];
                end
            elseif contains(words{1},'''')
                % if it's a future variable
                % Replace ' to GDSGE_GRID
                loopCode = line;
                f_var_name{end+1} = regexprep(words{1},'''','');
                
                loopCode = replace_grid_variable(loopCode,num_interp,var_interp_name,num_policy,var_policy_name,num_shock_var,var_shock_name,f_var_name,fVarPrefix);
                loopCode = my_ccode(my_sym(loopCode));
                % PATCHED for str2sym in R2024a
                loopCode = strrep(loopCode, '_GDSGE_iter', '[GDSGE_iter-1]');
                loopCode = [loopCode, parseComment];
                
                % Replace future to a loop
                forLoopCode = [
                    'for(int GDSGE_iter=1; GDSGE_iter<=' ...
                    num2str(shock_num) '; GDSGE_iter++)' ...
                    LINE_BREAK '{' LINE_BREAK...
                    loopCode LINE_BREAK...
                    '}' LINE_BREAK ...
                    ];
                modelCode = [modelCode ...
                    forLoopCode];
            else
                if ~contains(words{1},'(')
                    c_var_name{end+1} = words{1};
                end
                
                [reductionCode1,reductionCounter,line] = extract_reduction_operator(line,'GDSGE_EXPECT',var_interp_name,num_interp,var_policy_name,num_policy,num_shock_var,var_shock_name,f_var_name,fVarPrefix,reductionCounter,shock_num,LINE_BREAK);
                [reductionCode2,reductionCounter,line] = extract_reduction_operator(line,'GDSGE_MIN',var_interp_name,num_interp,var_policy_name,num_policy,num_shock_var,var_shock_name,f_var_name,fVarPrefix,reductionCounter,shock_num,LINE_BREAK);
                [reductionCode3,reductionCounter,line] = extract_reduction_operator(line,'GDSGE_MAX',var_interp_name,num_interp,var_policy_name,num_policy,num_shock_var,var_shock_name,f_var_name,fVarPrefix,reductionCounter,shock_num,LINE_BREAK);
                [reductionCode4,reductionCounter,line] = extract_reduction_operator(line,'GDSGE_PROD',var_interp_name,num_interp,var_policy_name,num_policy,num_shock_var,var_shock_name,f_var_name,fVarPrefix,reductionCounter,shock_num,LINE_BREAK);
                
                
                line = replace_grid_variable(line,num_interp,var_interp_name,num_policy,var_policy_name,num_shock_var,var_shock_name,f_var_name,fVarPrefix);
                
                if ~contains(words{1},'(')
                    parsed_code = my_ccode(my_sym(line));
                    % PATCHED for str2sym in R2024a
                    parsed_code = strrep(parsed_code, '_GDSGE_iter', '[GDSGE_iter-1]');
                    modelCode = [modelCode ...
                        reductionCode1 reductionCode2 reductionCode3 reductionCode4 LINE_BREAK ...
                        parsed_code parseComment ...
                        LINE_BREAK];
                else
                    % only process rhs
                    assign_words = strsplit(line,'=');
                    assign_lhs = assign_words{1};
                    assign_rhs = assign_words{2};
                    
                    assign_rhs_c_code = my_ccode(my_sym(assign_rhs));
                    assign_rhs_c_code_words = strsplit(assign_rhs_c_code,'=');
                    
                    modelCode = [modelCode ...
                        reductionCode1 reductionCode2 reductionCode3 reductionCode4 LINE_BREAK ...
                        assign_lhs '=' assign_rhs_c_code_words{2} parseComment LINE_BREAK];
                end
            end % Assign to future variable
        end % Assignment
        %{
    catch ME
        fprintf(2,'Error found in model block, at Line\n%s\n', line0);
        if strcmp(ME.identifier,'symbolic:kernel:RecursiveDefinitionMaxdepth')
            fprintf(2, 'Identifiers T, I, i are preserved by the MATLAB symbolic toolbox and cannot be used as variables in the model block.\n');
        end
        rethrow(ME);
    end
        %}
end

%% Declar variables
declareCode = ['adouble ' ...
    my_strjoin(strcat(unique(c_var_name),'(0)'),', ') ...
    ';' LINE_BREAK
    ];
declareCodeStack = '';
declareCodeStack = [declareCodeStack 'adouble ' ...
    my_strjoin(strcat(unique(f_var_name),[fVarPrefix '[' num2str(shock_num) ']']), ', ') ...
    ';' LINE_BREAK
    ];
declareCodeStack = [declareCodeStack 'adouble ' ...
    my_strjoin(strcat(vec_var_name,[fVarPrefix '['],vec_var_length,']'), ', ') ...
    ';' LINE_BREAK
    ];
declareCodeHeap = '';
declareCodeHeap = [declareCodeHeap 'vector<adouble> ' ...
    my_strjoin(strcat(unique(f_var_name),[fVarPrefix '(' num2str(shock_num) ')']), ', ') ...
    ';' LINE_BREAK
    ];
declareCodeHeap = [declareCodeHeap 'vector<adouble> ' ...
    my_strjoin(strcat(vec_var_name,[fVarPrefix '('],vec_var_length,')'), ', ') ...
    ';' LINE_BREAK
    ];
declareCode = [declareCode '#if MAXDIM>MAX_STACK_DIM' LINE_BREAK ...
    declareCodeHeap ...
    '#else' LINE_BREAK ...
    declareCodeStack ...
    '#endif' LINE_BREAK
    ];
for j=1:length(f_var_name)
    declareCode = [declareCode ...
        '#define ' f_var_name{j} '_GRID(idx) ' f_var_name{j} fVarPrefix '[int(idx)-1]' LINE_BREAK ...
        '#define ' f_var_name{j} '(idx) ' f_var_name{j} fVarPrefix '[int(idx)-1]' LINE_BREAK];
end
for j=1:length(vec_var_name)
    declareCode = [declareCode ...
        '#define ' vec_var_name{j} '(idx) ' vec_var_name{j} fVarPrefix '[int(idx)-1]' LINE_BREAK];
end

%% Declare var_aux
for j=1:length(var_aux_name)
    var_aux_current = var_aux_name{j};
    if ~ismember(var_aux_current,c_var_name) && ~ismember(var_aux_current,f_var_name)
        if var_aux_length(j)>1
            declareCode = [declareCode 'adouble ' ...
                var_aux_current fVarPrefix '[' num2str(var_aux_length(j)) '];' LINE_BREAK];
            declareCode = [declareCode ...
                '#define ' var_aux_current '_GRID(idx) ' var_aux_current fVarPrefix '[int(idx)-1]' LINE_BREAK ...
                '#define ' var_aux_current '(idx) ' var_aux_current fVarPrefix '[int(idx)-1]' LINE_BREAK];
        else
%             declareCode = [declareCode 'adouble ' ...
%                 var_aux_current ';' LINE_BREAK];
        end
    end
end

%% Argument code
argCode = '';
for j=1:num_policy
    if var_policy_length(j)
        argCode = [argCode 'adouble* ' var_policy_name{j} '_GDSGE_GRID = GDSGE_x+(' num2str(var_policy_loc(j)-1) ');' LINE_BREAK ...
            '#define ' var_policy_name{j} '_GRID(idx) ' var_policy_name{j} '_GDSGE_GRID[int(idx)-1]' LINE_BREAK ...
            '#define ' var_policy_name{j} '(idx) ' var_policy_name{j} '_GDSGE_GRID[int(idx)-1]' LINE_BREAK];
    else
        argCode = [argCode 'adouble& ' var_policy_name{j} ' = GDSGE_x[' num2str(var_policy_loc(j)-1) '];' LINE_BREAK];
    end
end

%% Equation code
if ~isempty(equationCode)
    lines = strsplit(equationCode,{char(10),char(13)});
    equationCode = '';
    num_equations = 0;
    equationPrefix = 'GDSGE_EQ';
    
    process_equation_code_inline = @(seg,num_equations) process_equation_line(seg,num_equations,...
        shock_num,num_interp,var_interp_name,num_policy,var_policy_name,num_shock_var,var_shock_name,f_var_name,equationPrefix,fVarPrefix,LINE_BREAK);
    
    line_loc = 0;
    in_if_block = false;
    in_for_block = false;
    while line_loc < length(lines)
        line_loc = line_loc+1;
        
        line = lines{line_loc};
        % Remove comments
        [~,loc_start,loc_end] = regexp(line,'%.*$','tokenExtents','all');
        line(loc_start:loc_end)=[];
        line = strtrim(line);
        if strcmp(line,'')
            continue;
        end
        
        % Process if block
        words = strsplit(line,' ');
        if strcmp(words{1},'for')
            % Begins a for block
            words = strsplit(line,{' ','='});
            % Get counter of the for block
            counter_var_name = words{2};
            range_name = words{3};
            
            % Split range
            range_name_words = strsplit(range_name,':');
            if length(range_name_words)==2
                % This is a step 1 range
                step = '1';
                rangeStart = range_name_words{1};
                rangeEnd = range_name_words{2};
            elseif length(range_name_words)==3
                % This is a step n range
                step = range_name_words{2};
                rangeStart = range_name_words{1};
                rangeEnd = range_name_words{3};
            else
                error('Range not specified correctly in for block');
            end
            
            % Write the for block header
            blockCode = ['for (int ' counter_var_name '=' rangeStart ';' counter_var_name '<=' rangeEnd ';' ...
                counter_var_name '=' counter_var_name '+' step ')' LINE_BREAK ...
                '{' LINE_BREAK];
            forBlockLoopLength = length(eval(range_name));
            eq_assignments = {};
            in_for_block = true;
            continue;
        end
        
        
        if strcmp(words{1},'if')
            % This is a if block
            % Get conditions of the if block
            condition_idx = 1;
            conditions = {};
            eq_assignments = {};
            
            conditions{condition_idx} = strcat(words{2:end});
            eq_assignments{condition_idx} = {};
            
            % Find else and else if block
            found_else = false;
            in_if_block = true;
            continue;
        end
        
        if strcmp(words{1},'else')
            % else block found
            found_else = true;
            
            condition_idx = condition_idx+1;
            conditions{condition_idx} = '1';
            eq_assignments{condition_idx} = {};
            continue;
        end
        
        if strcmp(words{1},'elseif')
            % elseif block found
            
            condition_idx = condition_idx+1;
            conditions{condition_idx} = strcat(words{2:end});
            eq_assignments{condition_idx} = {};
            continue;
        end
        
        if strcmp(words{1},'end') && in_if_block
            in_if_block = false;
            if ~found_else
                error('else not found in if block in equation definitions');
            end
            
            % Process codes if the if condition
            total_condition_branches = condition_idx;
            num_assignments_in_first_branch = length(eq_assignments{1});
            blockCode = ['if (' conditions{1} ')' LINE_BREAK '{' LINE_BREAK];
            num_equations_before_if = num_equations;
            for i_assignment = 1:num_assignments_in_first_branch
                line = eq_assignments{1}{i_assignment};
                [line,num_equations] = process_equation_code_inline(line,num_equations);
                blockCode = [blockCode ...
                    line,LINE_BREAK];
            end
            blockCode = [blockCode '}' LINE_BREAK];
            equationCode = [equationCode blockCode];
            
            for i_branch = 2:total_condition_branches
                num_assignments = length(eq_assignments{i_branch});
                if num_assignments ~= num_assignments_in_first_branch
                    error('num of assignments not consistent in if...else block in equation definitions');
                end
                
                num_equations = num_equations_before_if;
                blockCode = ['else if (' conditions{i_branch} ')' LINE_BREAK '{' LINE_BREAK];
                for i_assignment = 1:num_assignments
                    line = eq_assignments{i_branch}{i_assignment};
                    [line,num_equations] = process_equation_code_inline(line,num_equations);
                    blockCode = [blockCode ...
                        line,LINE_BREAK];
                end
                blockCode = [blockCode '}' LINE_BREAK];
                equationCode = [equationCode blockCode];
            end
            
            continue;
        end
        
        if strcmp(words{1},'end') && in_for_block
            in_for_block = false;
            num_assignments = length(eq_assignments);
            for i_assignment = 1:num_assignments
                line = eq_assignments{i_assignment};
                [line,~] = process_equation_code_inline(line,num_equations);
                % Replace counters
                pat = ['(?<=',equationPrefix,'\[)[0-9]+(?=\])'];
                replaceStr = [num2str(num_equations), '+', num2str(num_assignments), '*(', counter_var_name,'-1)+', num2str(i_assignment-eval(rangeStart))];
                line = regexprep(line, pat, replaceStr);
                blockCode = [blockCode ...
                    line,LINE_BREAK];
            end
            blockCode = [blockCode '}' LINE_BREAK];
            equationCode = [equationCode blockCode];
            num_equations = num_equations + num_assignments*forBlockLoopLength;
            
            continue;
        end
        
        if in_if_block
            eq_assignments{condition_idx} = [eq_assignments{condition_idx},line];
        elseif in_for_block
            eq_assignments = [eq_assignments,line];
        else
            % Not in a if block
            [line,num_equations] = process_equation_code_inline(line,num_equations);
            %
            equationCode = [equationCode ...
                line LINE_BREAK];
        end
    end
    if num_equations~=num_total_policy
        errMsg = sprintf('Number of equations (%d) is not equal to number of endogenous variables (%d)',num_equations,num_total_policy);
        error(errMsg);
    end
end
% equationCode = ['adouble ' equationPrefix '[' num2str(num_equations) '];' char(13) equationCode];
end

function [line,num_equations] = process_equation_line(line,num_equations,...
    shock_num,num_interp,var_interp_name,num_policy,var_policy_name,num_shock_var,var_shock_name,f_var_name,equationPrefix,fVarPrefix,LINE_BREAK)
% Replace the prime variables
line0 = line;
parseComment = [' //parsed from gmod Line: ' line0];
try
    if contains(line,'''')
        % This contains future variables
        line = regexprep(line,'''','(GDSGE_iter)');
        line = replace_grid_variable(line,num_interp,var_interp_name,num_policy,var_policy_name,num_shock_var,var_shock_name,f_var_name,fVarPrefix);
        % line = [equationPrefix '[' num2str(num_equations) '+GDSGE_iter]=' line];
        % PATCHED for str2sym in R2024a; first call str2sym then a replacement
        line = [equationPrefix '=' line];
        line = my_ccode(my_sym(line));
        line = strrep(line, '_GDSGE_iter', '[GDSGE_iter-1]');
        line = my_regexprep(line, equationPrefix, [equationPrefix '[' num2str(num_equations-1) '+GDSGE_iter]']);
        line = get_for_loop_code(shock_num,line,LINE_BREAK);
        num_equations = num_equations+shock_num;
    else
        line = replace_grid_variable(line,num_interp,var_interp_name,num_policy,var_policy_name,num_shock_var,var_shock_name,f_var_name,fVarPrefix);
        % line = [equationPrefix '[' num2str(num_equations+1) ']=' line];
        % PATCHED for str2sym in R2024a; first call str2sym then a replacement
        line = [equationPrefix '=' line];
        line = my_ccode(my_sym(line));
        line = strrep(line, '_GDSGE_iter', '[GDSGE_iter-1]');
        line = my_regexprep(line, equationPrefix, [equationPrefix '[' num2str(num_equations+1-1) ']']);
        line = [line, parseComment];
        num_equations = num_equations+1;
    end
catch ME
    fprintf(2,'Error found in model block, at Line\n%s\n', line0);
    rethrow(ME);
end
end

function forLoopCode = get_for_loop_code(shock_num,loopCodes,LINE_BREAK)
forLoopCode = [
    'for(int GDSGE_iter=1; GDSGE_iter<=' ...
    num2str(shock_num) '; GDSGE_iter++)' ...
    LINE_BREAK '{' LINE_BREAK...
    loopCodes LINE_BREAK...
    '}' LINE_BREAK ...
    ];
end

function [blockCode,remainingCode] = extract_segment(code,blockName)
% Extract block segment
blockStart = my_regexp(code,blockName);
if length(blockStart)>1
    error('GDSGE:multipleModel','Multiple segment not allowed');
end
% Search for end from equationStart
for j=blockStart:length(code)
    if strcmp(code(j:j+3),'end;')
        break;
    end
end
blockCode = code(blockStart:j+3);
code(blockStart:j+3) = [];
blockCode = regexprep(blockCode,blockName,'');
blockCode = regexprep(blockCode,'end;','');

remainingCode = code;
end

function [simulateCode, remainingCode] = extract_simulate_segment(code)
% Extract simulate segment
simulateStart = regexp(code,'simulate;','all');
if length(simulateStart)>1
    error('GDSGE:multipleModel','Multiple simulate segment not allowed');
end
% Search for end from modelStart
for j=simulateStart:length(code)
    if strcmp(code(j:j+3),'end;')
        break;
    end
end
simulateCode = code(simulateStart:j+3);
code(simulateStart:j+3) = [];
simulateCode = regexprep(simulateCode,'simulate;','');
simulateCode = regexprep(simulateCode,'end;','');

remainingCode = code;
end

function check_reserve_word(words)
reserveWordsList = {'parameters','var_state','shock_num','var_shock','var_tensor','var_policy','inbound','initial','var_interp','model','end','equations','f','grad','x','i','data','otherData','j','MAX',...
    'shock_trans','shock_num','shock','diff','cxx','E','pi','I'};
for j=1:length(words)
    if find(ismember(reserveWordsList,words{j}))
        error('GDSGE:nameConflict','reserved name not allowed: %s',words{j});
    end
end
end

function str = my_strjoin(str_array,delim)
if length(str_array)>=1
    str = strjoin(str_array,delim);
else
    str = '';
end
end

function [parameters_size,pre_params,params] = get_parameters(preCode,parameters_name)
% Print preCode to a file
if ~isempty(preCode)
    fileID = fopen(['gdsge_precode.m'],'w');
    fprintf(fileID,'%s',preCode);
    fclose(fileID);
end
try
    gdsge_precode;
catch ME
    %{
    errorMsg = strsplit(ME.message,{char(10),char(13)});
    errorLine = regexp(ME.message,'(?<=File: temp_precode.m Line: )(\w*)','match');
    errorCol = regexp(ME.message,'(?<=File: temp_precode.m Line: \w* Column: )(\w*)','match');
    %}
    preCodeLine = strsplit(preCode,{char(10),char(13)});
    fprintf('\n');
    for i=1:length(preCodeLine)
        fprintf(1,'%d\t%s\n', i,preCodeLine{i});
    end
    %{
    error(['Error found outside the model block. Line: ', errorLine{1}, ' Column: ', errorCol{1}]);
    %}
    rethrow(ME)
end

parameters_size = [];
for j=1:length(parameters_name)
    % Check if parameters exist
    if ~exist(parameters_name{j},'var')
        warning(['Parameter ', parameters_name{j}, ' not defined']);
        % Default parameter
        evalString = [parameters_name{j},'=nan;'];
        eval(evalString);
    end
    evalString = ['parameters_size(j) = numel(' parameters_name{j} ');'];
    eval(evalString);
    %{
    try
        eval(evalString);
    catch ME
        error(['Parameter ', parameters_name{j}, ' not defined']);
    end
    %}
end

pre_params = v2struct(INTERP_ORDER,USE_FINITE_DIFF,EXTRAP_ORDER,USE_SPLINE,USE_ASG,USE_SPARSE_JACOBIAN,USE_PCHIP,AsgMaxLevel,SIMU_INTERP,SIMU_RESOLVE,UseModelId,shock_num,shock_trans,REMOVE_NULL_STATEMENTS);
evalString = ['params = v2struct(', my_strjoin(parameters_name,','), ');'];
eval(evalString);
end

function code = process_deprecate_single(code,oldWord,newWord)
startIdx = regexp(code,['((?<=^)|(?<=\W))' oldWord], 'ONCE');
if isempty(startIdx)
    return
end

warning('%s is going to be deprecated, use %s instead\n',oldWord, newWord);
code = regexprep(code,['((?<=^)|(?<=\W))' oldWord],newWord);
end

function code = process_deprecate(code)
deprecatedList = {
    'InterpOrder', 'INTERP_ORDER'
    'SimuInterp',   'SIMU_INTERP'
    'SimuResolve',  'SIMU_RESOLVE'
    'ExtrapOrder',  'EXTRAP_ORDER'
    };

for i=1:size(deprecatedList,1)
    code = process_deprecate_single(code,deprecatedList{i,1},deprecatedList{i,2});
end
end

function y = my_ccode(x)
y = ccode(x);
y = strrep(y,'3.141592653589793','pi');

%{
temp = regexp(y,'GDSGE_iter\*(?<digit>\d+\.0)\+shock\-\k<digit>','tokenExtents');
for i=1:length(temp)
% Prepare string replace
digit = y(temp{i}(1):temp{i}(2));
str_to_be_replaced = ['GDSGE_iter*',digit,'+shock-',digit];
digit_trimed = int2str(floor(str2double(digit)));
str_replacing = ['GDSGE_iter*',digit_trimed ,'+shock-',digit_trimed];
y = strrep(y,str_to_be_replaced,str_replacing);
end
%}

end

function allCode = rec_extract_loop_seg(code)
allButUnderscore = '[^a-zA-Z0-9]';

BEFORE_A_WORD = ['((?<=^)|(?<=' allButUnderscore '))'];
AFTER_A_WORD = ['((?=$)|(?=' allButUnderscore '))'];
code = [newline,code];
regex = '#foreach (?<id>\w*)\s*in\s*([^\n]*)(?:\n)(.*?)#endfor \k<id>';
subcodeCells = regexp(code,regex,'tokens');
allCode = [];
if ~isempty(subcodeCells)
    subcodeStarts = regexp(code,regex,'start');
    subcodeEnds = regexp(code,regex,'end');
    for i=1:length(subcodeCells)
        subcodeOriginal = subcodeCells{i}{3};
        id = subcodeCells{i}{1};
        valueList = strsplit(subcodeCells{i}{2},' ');
        % Expand macro
        subcode = [];
        for j=1:length(valueList)
            subcode = [subcode,newline,regexprep(subcodeOriginal,[BEFORE_A_WORD,'#',id,AFTER_A_WORD],strtrim(valueList{j}))];
        end
        % Attach head
        if i==1
            allCode = [allCode,code(1:subcodeStarts(i)-1),newline];
        else
            allCode = [allCode,code(subcodeEnds(i-1)+1:subcodeStarts(i)-1),newline];
        end
        % Attach body
        allCode = [allCode,rec_extract_loop_seg(subcode),newline];
        % Attach tail
        if i==length(subcodeCells)
            allCode = [allCode,code(subcodeEnds(i)+1:end)];
        end
    end
else
    allCode = code;
end
end

function code = process_if_macro(code)
code0 = code;
regex = '#if ([^\n]*)(.*?)#endif';
conds = regexp(code,regex,'tokens');
condStarts = regexp(code,regex,'start');
condEnds = regexp(code,regex,'end');
for i=1:length(conds)
    condeval = eval(conds{i}{1});
    if ~condeval
        code = strrep(code,code0(condStarts(i):condEnds(i)),'');
    else
        code = strrep(code,code0(condStarts(i):condEnds(i)),conds{i}{2});
    end
end

end

function expr = my_sym(str)
expr = str2sym(str);
end
