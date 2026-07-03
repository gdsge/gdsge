function [ok, why] = ensurePyenv(port)
% ENSUREPYENV  Point MATLAB's pyenv at the uv-managed venv interpreter and make
%   the gdsge_sympy package importable. Idempotent; callable from matlab -batch.
%   ok is true iff `import sympy` and `import gdsge_sympy` succeed; why is ''
%   on success, else a short human-readable reason.
%   The SymPy backend is the ONLY consumer; the default autodiff path never
%   calls this, so a missing/unsynced venv only affects sympy requests.
%
%   Python runs OUT-OF-PROCESS. The bridge only exchanges JSON strings (see
%   callSympy), so the IPC boundary is transparent and results are identical.
%   This isolates MATLAB from a native crash (access violation in MATLAB's
%   error-handling path) observed when the Python interpreter is loaded
%   in-process alongside the OpenMP solver MEX over a long run.
%
%   Other code (startup.m, another toolbox) may have loaded a foreign Python
%   into the session before us; pyenv Version/ExecutionMode are only settable
%   while NotLoaded/Terminated. Strategy when Loaded: probe by importing (any
%   SymPy-capable interpreter works); on failure terminate-and-repoint if
%   OutOfProcess, warn (restart MATLAB) if InProcess.
%
%   port (tests only): struct of effect functions {venvPython,getEnv,setEnv,
%   terminate,tryImports} so unit tests drive every branch without touching
%   the process-global pyenv state.
here = fileparts(mfilename('fullpath'));                 % src/+gdsge/+codegen/+sympy
repoRoot = fileparts(fileparts(fileparts(fileparts(here))));
pyextDir = fullfile(repoRoot, 'pyext');
if nargin < 1
    port = struct( ...
        'venvPython', @() gdsge.codegen.sympy.venvPython(pyextDir), ...
        'getEnv',     @() pyenv, ...
        'setEnv',     @(p) pyenv('Version', p, 'ExecutionMode', 'OutOfProcess'), ...
        'terminate',  @() terminate(pyenv), ...
        'tryImports', @() tryImports(pyextDir));
end

ok = false;
venvPy = port.venvPython();
if isempty(venvPy)
    why = 'the SymPy environment is not set up (run gdsge_setup_sympy)';
    return;
end

try
    pe = port.getEnv();
    if pe.Status ~= "Loaded"                 % NotLoaded or Terminated: settable
        port.setEnv(venvPy);
        ok = port.tryImports();
        why = importsWhy(ok);
        return;
    end

    % Python already loaded (possibly a foreign interpreter loaded by other
    % code). Probe by importing rather than comparing Executable paths — on
    % macOS the venv python is a symlink and pyenv may report a resolved path.
    if port.tryImports()
        ok = true; why = '';
        return;
    end

    if pe.ExecutionMode == "OutOfProcess"
        % Wrong interpreter, but replaceable: terminate and re-point.
        port.terminate();
        port.setEnv(venvPy);
        ok = port.tryImports();
        why = importsWhy(ok);
    else
        % InProcess: cannot be unloaded without restarting MATLAB.
        why = sprintf( ...
            ['a different Python (%s) was already loaded in-process by other ' ...
             'code in this MATLAB session and does not provide SymPy. Restart ' ...
             'MATLAB and re-run; GDSGE will then load its own Python ' ...
             'out-of-process.'], pe.Executable);
        warning('gdsge:sympy:inProcessMismatch', '%s', why);
    end
catch ME
    ok = false;
    why = ME.message;
end
end

function w = importsWhy(ok)
if ok
    w = '';
else
    w = 'the venv interpreter failed to import sympy/gdsge_sympy (re-run gdsge_setup_sympy)';
end
end

function ok = tryImports(pyextDir)
% Make pyext/ importable (the package lives at pyext/gdsge_sympy) and probe.
ok = false;
try
    if count(py.sys.path, pyextDir) == 0
        insert(py.sys.path, int32(0), pyextDir);
    end
    py.importlib.import_module('sympy');
    py.importlib.import_module('gdsge_sympy');
    ok = true;
catch
end
end
