function d = localVenvDir()
% LOCALVENVDIR  Fallback location for the SymPy venv when it cannot live next to
%   pyext/. A virtual environment always needs symlinks (the interpreter, plus a
%   lib64 symlink on Linux), so on filesystems without symlink support — notably
%   MATLAB Online's "/MATLAB Drive" — the venv is placed here instead. Derived
%   from the user's home so gdsge_setup_sympy (creates it) and
%   gdsge.codegen.sympy.venvPython (finds it) agree on the path.
if ispc; base = getenv('USERPROFILE'); else; base = getenv('HOME'); end
if isempty(base); base = tempdir; end
d = fullfile(base, '.gdsge', 'sympy-venv');
end
