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