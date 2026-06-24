function ok = ensurePyenv()
% ENSUREPYENV  Point MATLAB's pyenv at the uv-managed venv interpreter and make
%   the gdsge_sympy package importable. Idempotent; callable from matlab -batch.
%   Returns true iff `import gdsge_sympy` and `import sympy` both succeed.
%   The SymPy backend is the ONLY consumer; the default autodiff path never
%   calls this, so a missing/unsynced venv only affects sympy requests.
ok = false;
try
    here = fileparts(mfilename('fullpath'));                 % src/+gdsge/+codegen/+sympy
    repoRoot = fileparts(fileparts(fileparts(fileparts(here))));
    pyextDir = fullfile(repoRoot, 'pyext');
    venvPy = gdsge.codegen.sympy.venvPython(pyextDir);   % pyext/.venv or local fallback
    if isempty(venvPy)
        return;     % env not synced: caller raises the informative guard error
    end
    pe = pyenv;
    if pe.Status ~= "Loaded"
        % Run Python OUT-OF-PROCESS. The SymPy bridge only exchanges JSON strings
        % (see callSympy), so the IPC boundary is transparent and results are
        % identical. This isolates MATLAB from a native crash (access violation in
        % MATLAB's error-handling path) observed when the Python interpreter is
        % loaded in-process alongside the OpenMP solver MEX over a long run.
        % ExecutionMode is only settable while Python is NotLoaded.
        pyenv('Version', venvPy, 'ExecutionMode', 'OutOfProcess');
    end
    % Make pyext/ importable (the package lives at pyext/gdsge_sympy).
    if count(py.sys.path, pyextDir) == 0
        insert(py.sys.path, int32(0), pyextDir);
    end
    py.importlib.import_module('sympy');
    py.importlib.import_module('gdsge_sympy');
    ok = true;
catch
    ok = false;
end
end
