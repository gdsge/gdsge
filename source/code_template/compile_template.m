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
MISC = 'sig_die.o cabs.o ctype.o exit_.o ';
POW = 'pow_di.o pow_dd.o ';
CX = '';
DCX = 'z_div.o ';
REAL = 'r_mod.o ';
DBL = 'd_sign.o ';
INT = 'i_len.o i_indx.o ';
HALF = '';
CMP = 'l_ge.o l_le.o l_lt.o ';
EFL = '';
CHAR = 'f77_aloc.o s_cat.o s_cmp.o s_copy.o ';
I77 =	['backspac.o close.o dfe.o dolio.o due.o endfile.o err.o '...
    'fmt.o fmtlib.o ftell_.o iio.o ilnw.o inquire.o lread.o lwrite.o '...
    'open.o rdfmt.o rewind.o rsfe.o rsli.o rsne.o sfe.o sue.o '...
    'typesize.o uio.o util.o wref.o wrtfmt.o wsfe.o wsle.o wsne.o xwsne.o '];
QINT = '';
TIME = 'etime_.o';
OFILES = [MISC POW CX DCX REAL DBL INT HALF CMP EFL CHAR I77 TIME];
OFILES = strsplit(OFILES, ' ');
LIBF2C_SRC = strrep(OFILES, '.o', '.c');
LIBF2C_SRC = strcat([' dep/snopt/libf2c/'], LIBF2C_SRC);
LIBF2C_SRC = strjoin(LIBF2C_SRC);

SNOPT_SRC = dir('dep/snopt/snopt/*.c');
SNOPT_SRC = strcat([' dep/snopt/snopt/'], {SNOPT_SRC.name});
SNOPT_SRC = strjoin(SNOPT_SRC);

flag0 = '';
flag1 = '';
flag2 = '';
flag3 = '';
flag1 = [flag1 ' -DMAXDIM=GDSGE_MAXDIM -DINTERP_ORDER=GDSGE_INTERP_ORDER -DEIGEN_DONT_PARALLELIZE ' 'EXTRA_DEF'];

cppFileName = [' ' current_folder '/mex_MODEL_NAME.cpp'];

ME = [];

if ispc
    mexCommand = 'mex';
    flag0 = [flag0 ' -DNO_TRUNCATE'];
    % err.c
    flag0 = [flag0 ' -DNO_ISATTY'];
    % open.c
    flag0 = [flag0 ' -DMSDOS'];
    
    compilerWin = mex.getCompilerConfigurations('C++');
    isIntel = strcmp(compilerWin.Manufacturer,'Intel');
    if isIntel>0
        %{
        flag2 = [' OPTIMFLAGS="/O2 /DNDEBUG /QxCORE-AVX2" COMPFLAGS="$COMPFLAGS /wd4267 /wd4068 /wd4091 /Qansi-alias /openmp /QxCORE-AVX2"'];
        %}
        flag2 = [' OPTIMFLAGS="/O2 /DNDEBUG" COMPFLAGS="$COMPFLAGS /wd4267 /wd4068 /wd4091 /Qansi-alias /openmp"'];
    else
        flag2 = [' OPTIMFLAGS="/O2 /DNDEBUG" COMPFLAGS="$COMPFLAGS /wd4267 /wd4068 /wd4091 /diagnostics:caret /openmp /Z7"'];
        flag3 = [' LINKOPTIMFLAGS="/DEBUG:FULL"'];
        %{
        flag2 = [' OPTIMFLAGS="/O2 /DNDEBUG /arch:AVX2" COMPFLAGS="$COMPFLAGS /wd4267 /wd4068 /wd4091 /diagnostics:caret /openmp /arch:AVX2"'];
        %}
    end
    
    if exist('DEBUG','var')~=0
        flag2 = [' -v OPTIMFLAGS="/Og /Z7 /Zo /debug:full " COMPFLAGS="$COMPFLAGS "'];
        flag3 = [' LINKOPTIMFLAGS="/DEBUG:FULL"'];
    end
    
    compileString = [mexCommand ' ' flag0 flag1 flag2 flag3 cppFileName ' ./lib/essential_blas.lib' ' -outdir "' current_folder '"' ' -I"INCLUDE_FOLDER"'];

elseif isunix && ~ismac
    mexCommand = 'mex CXX=icc CXXOPTIMFLAGS= LDOPTIMFLAGS=';
    
    flag2 = [' CXXFLAGS="$CXXFLAGS -std=c++11 -openmp -no-inline-max-size -no-inline-max-total-size -ansi-alias -w"'];
    
    flag3 = [' LDFLAGS="$LDFLAGS -static-intel -liomp5"'];
    
    compileString = [mexCommand ' ' flag0 flag1 flag2 flag3 cppFileName LIBF2C_SRC SNOPT_SRC ' -outdir "' current_folder '"' ' -I"INCLUDE_FOLDER"'];
elseif ismac
    mexCommand = 'mex CXX=icc CXXOPTIMFLAGS="-O1 -DNDEBUG -QxCORE-AVX2" LDOPTIMFLAGS=';
    
    flag2 = [' CXXFLAGS="$CXXFLAGS -O1 -DADEPT_THREAD_LOCAL=__thread -DNDEBUG -QxCORE-AVX2 -std=c++14 -qopenmp -qopenmp-link=static -ansi-alias -w"'];
    
    flag3 = [' LDFLAGS="$LDFLAGS essential_blas.dylib -liomp5 -lpthread -qopenmp-link=static"'];
    
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