% This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE
function [IterRslt,IterFlag] = iter_MODEL_NAME(GDSGE_OPTIONS)
%% Add path
if ispc
    BLAS_FILE = 'essential_blas.dll';
    PATH_DELIMITER = ';';
    
    GDSGE_TOOLBOX_ROOT = fileparts(which(BLAS_FILE));
    if ~any(strcmp(strsplit(getenv('PATH'),PATH_DELIMITER),GDSGE_TOOLBOX_ROOT))
        setenv('PATH',[getenv('PATH'),PATH_DELIMITER,GDSGE_TOOLBOX_ROOT]);
    end
    
    clear BLAS_FILE PATH_DELIMITER GDSGE_TOOLBOX_ROOT
elseif ismac
    if exist('./essential_blas.dylib','file') == 0
        copyfile(which('essential_blas.dylib'),'./');
    end
end

%% Iter code starts here
PRE_CODE

if nargin>=1
    if isfield(GDSGE_OPTIONS,'WarmUp')
        if isfield(GDSGE_OPTIONS.WarmUp,'var_tensor')
            v2struct(GDSGE_OPTIONS.WarmUp.var_tensor);
        end
    end
    v2struct(GDSGE_OPTIONS)
end
  
ASSERT_CODE

%% Solve the last period problem
ITER_INIT

%% Solve the infinite horizon problem
ITER_INF_HORIZON    

%% Return the success flag
IterFlag = 0;
end
