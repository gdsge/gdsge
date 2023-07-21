% This file is a part of GDSGE. License is Apache License, Version 2.0: http://github.com/gdsge/gdsge/LICENSE
function compile_MODEL_NAME(GDSGE_OPTIONS)
clear mex;
current_folder = fileparts(mfilename('fullpath'));
gdsge_folder = 'GDSGE_FOLDER';

cd(gdsge_folder);

if nargin>=1
    v2struct(GDSGE_OPTIONS)
end

try
flag0 = '';
flag1 = '';
flag2 = '';
flag3 = '';
flag1 = [flag1 ' -DMAXDIM=GDSGE_MAXDIM -DINTERP_ORDER=GDSGE_INTERP_ORDER -DEIGEN_DONT_PARALLELIZE ' 'EXTRA_DEF'];

cppFileName = [' ' current_folder '/mex_MODEL_NAME.cpp'];

ME = [];

if ispc
    compiler_name = mex.getCompilerConfigurations('C++').Name;
    if contains(compiler_name, 'MinGW64')
        compiler_name = 'MinGW64';
    elseif contains(compiler_name, 'Microsoft Visual C++')
        compiler_name = 'MSVC';
    else
        compiler_name = 'Intel';
    end
elseif isunix
    % only support g++ under linux/macOS
    compiler_name = 'g++';
end

if ispc
    mexCommand = 'mex';
    switch compiler_name
        case 'MinGW64'
            flag2 = [' CFLAGS="$CFLAGS -w -fopenmp -fpermissive -DADEPT_THREAD_LOCAL=__thread"'];
            flag3 = [sprintf(' LDFLAGS="$LDFLAGS -w -fopenmp"')];
            copyfile("mingw/asg_mex.mexw64", current_folder);
        case 'MSVC'
            flag2 = [' -DUSE_MKL OPTIMFLAGS="/O2 /DNDEBUG" COMPFLAGS="$COMPFLAGS /wd4267 /wd4068 /wd4091 /diagnostics:caret /openmp /Z7"'];
            flag3 = [' LINKOPTIMFLAGS="/DEBUG:FULL"'];
        otherwise
            % default to intel
            flag2 = [' -DUSE_MKL OPTIMFLAGS="/O2 /DNDEBUG" COMPFLAGS="$COMPFLAGS /wd4267 /wd4068 /wd4091 /Qansi-alias /openmp"'];
    end

    if exist('DEBUG','var')~=0
        switch compiler_name
            case 'MSVC'
                flag2 = [' -v OPTIMFLAGS="/Og /Z7 /Zo /debug:full " COMPFLAGS="$COMPFLAGS "'];
                flag3 = [' LINKOPTIMFLAGS="/DEBUG:FULL"'];
            otherwise
                flag2 = [' -v OPTIMFLAGS="/Og /Z7 /Zo /debug:full " COMPFLAGS="$COMPFLAGS "'];
                flag3 = [' LINKOPTIMFLAGS="/DEBUG:FULL"'];
        end
    end
    
    switch compiler_name
        case 'MinGW64'
            link_to_lib = '';
        case 'MSVC'
            link_to_lib = sprintf(' "%s/essential_blas.lib"',gdsge_folder);
        otherwise
            error('compiler not supported');
    end

    compileString = [mexCommand ' ' flag0 flag1 flag2 flag3 cppFileName link_to_lib ' -outdir "' current_folder '"' ' -I"INCLUDE_FOLDER"'];

elseif isunix && ~ismac
    mexCommand = 'mex CXX=icc CXXOPTIMFLAGS= LDOPTIMFLAGS=';
    
    flag2 = [' CXXFLAGS="$CXXFLAGS -std=c++11 -openmp -no-inline-max-size -no-inline-max-total-size -ansi-alias -w"'];
    
    flag3 = [' LDFLAGS="$LDFLAGS -static-intel -liomp5"'];
    
    compileString = [mexCommand ' ' flag0 flag1 flag2 flag3 cppFileName LIBF2C_SRC SNOPT_SRC ' -outdir "' current_folder '"' ' -I"INCLUDE_FOLDER"'];
elseif ismac
    mexCommand = 'mex CC=/usr/local/bin/gcc-8 CXX=/usr/local/bin/g++-8';
    
    flag2 = [' CFLAGS="$CFLAGS -w -fopenmp -fpermissive -fexceptions -DUSE_MKL -DADEPT_THREAD_LOCAL=__thread"'];
    flag3 = [sprintf(' LDFLAGS="$LDFLAGS essential_blas.dylib -w -fopenmp"')];
    
    compileString = [mexCommand ' ' flag0 flag1 flag2 flag3 cppFileName ' -outdir "' current_folder '"' ' -I"INCLUDE_FOLDER"'];
end
    eval(compileString);
catch ME
end
cd(current_folder);

if ~isempty(ME)
    rethrow(ME);
end
end