function [model,iterCode,cppCache,cppCode,codeSegment] = gdsge_codegen(modelName,options)
% This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE

fprintf('Parsing gmod file: ');
[model,iterCode,cppCode,simulateCode,compileCode,codeSegment] = gdsge_parser(modelName);

fileID = fopen(['iter_' modelName '.m'],'w');
fprintf(fileID,'%s',iterCode);
fclose(fileID);

if ~isempty(simulateCode)
    fileID = fopen(['simulate_' modelName '.m'],'w');
    fprintf(fileID,'%s',simulateCode);
    fclose(fileID);
end

cppCacheName = ['mex_' modelName '.cache'];
cppCache = '';
if exist(cppCacheName, 'file') ~= 0
    cppCache = fileread(cppCacheName);
end

fileID = fopen(['mex_' modelName '.cpp'],'w');
fprintf(fileID,'%s',cppCode);
fclose(fileID);
cppCode = fileread(['mex_' modelName '.cpp']);

fileID = fopen(['compile_' modelName '.m'],'w');
fprintf(fileID,'%s',compileCode);
fclose(fileID);
fprintf('ok\n');

if ~strcmp(cppCache, cppCode)
    fprintf('Compile mex file: \n');
    eval(['compile_' modelName]);
    fileID = fopen(cppCacheName,'w');
    fprintf(fileID,'%s',cppCode);
    fclose(fileID);
    fprintf('ok\n');
else
    fprintf('Header file not changed, skip compiling. \n');
end

if nargin>1 && ~isempty(options)
    if isfield(options,'GenCodeSegment') && options.GenCodeSegment==1
        fileID = fopen(['init_iter_' modelName '.m'],'w');
        fprintf(fileID,'%s',codeSegment.iterInitCode);
        fclose(fileID);
        
        fileID = fopen(['set_params_' modelName '.m'],'w');
        fprintf(fileID,'%s',codeSegment.setParamsCode);
        fclose(fileID);
        
        fileID = fopen(['prepare_space_' modelName '.m'],'w');
        fprintf(fileID,'%s',codeSegment.prepareSpaceCode);
        fclose(fileID);
        
        fileID = fopen(['construct_spline_' modelName '.m'],'w');
        fprintf(fileID,'%s',codeSegment.constructSplineCode);
        fclose(fileID);
        
        fileID = fopen(['solve_and_assign_' modelName '.m'],'w');
        fprintf(fileID,'%s',codeSegment.solveAndAssignCode);
        fclose(fileID);
    end
end
end