TolEq = 1e-6;
TolSol = 1e-8;
TolFun = 1e-8;
PrintFreq = 10;
NoPrint = 0;
SaveFreq = 10;
NoSave = 0;
SimuPrintFreq = 1000;
SimuSaveFreq = inf;
MaxIter = inf;
MaxMinorIter = inf;
num_samples = 1;
num_periods = 1000;
SolMaxIter = 200;

% task constants
MEX_TASK_INIT = 0;
MEX_TASK_INF_HORIZON = 1;

% Solver
UseBroyden = 0;
FiniteDiffDelta = 1e-6;

% DEBUG flag
GDSGE_DEBUG_EVAL_ONLY = 0;
GDSGE_USE_BROYDEN = 1;
GDSGE_USE_BROYDEN_NOW = 0;