clear;
addpath ../source
testList = {
    'HeatonLucas1996'
    'Bianchi2011_asg'
    'Mendoza2010'
    'CaoKS2016'
    'GLSW2020'
    'Barro_et_al_2017'
    };
for i=1:length(testList)
    fprintf('===========================================\n');
    testFolder = testList{i};
    fprintf(['Testing ',testFolder,'\n']);
    fprintf('===========================================\n');
    cd(testFolder);
    delete gdsge_*;
    delete mex_*;
    delete iter_*.m;
    delete simulate_*.m;
    run('test.m');
    fprintf('===========================================\n');
    fprintf(['Testing ',testFolder,' done!\n']);
    fprintf('===========================================\n');
    drawnow;
    cd ..
end