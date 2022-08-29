function eq = gdsge(modelName)
% This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE

fprintf('GDSGE: A Toolbox for Solving Global DSGE Models:\n');

% Code gen
[model,iterCode,cppCache,cppCode] = gdsge_codegen(modelName);

iterCacheName = ['iter_' modelName '.cache'];
iterCache = '';
if exist(iterCacheName,'file')~=0
    iterCache = fileread(iterCacheName);
end

iterRsltName = ['IterRslt_' modelName '.mat'];
if strcmp(iterCache,iterCode) && strcmp(cppCache,cppCode) && exist(iterRsltName,'file')~=0
    fprintf('Iter file not changed. IterRslt found. Load results: \n');
    load(iterRsltName);
    fprintf('ok\n');
else
    fprintf('Solve policies: \n');
    eval(['IterRslt = iter_' modelName ';']);
    save(iterRsltName,'IterRslt');
    fileID = fopen(iterCacheName,'w');
    fprintf(fileID,'%s',iterCode);
    fclose(fileID);
    fprintf('ok\n');
end

fprintf('Simulate: \n');
eval(['SmltRslt = simulate_' modelName '(IterRslt);']);
fprintf('ok\n');

eq = v2struct(model,IterRslt,SmltRslt);
end