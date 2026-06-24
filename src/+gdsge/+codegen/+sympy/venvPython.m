function p = venvPython(pyextDir)
% VENVPYTHON  Absolute path to the SymPy venv interpreter, or '' if the env is
%   not set up. Checks pyext/.venv (the normal location) first, then the local
%   fallback dir used on filesystems that cannot host a venv (see localVenvDir).
cands = { fullfile(pyextDir, '.venv'), gdsge.codegen.sympy.localVenvDir() };
for i = 1:numel(cands)
    if ispc
        p = fullfile(cands{i}, 'Scripts', 'python.exe');
    else
        p = fullfile(cands{i}, 'bin', 'python');
    end
    if exist(p, 'file') == 2; return; end
end
p = '';
end
