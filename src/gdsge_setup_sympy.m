function gdsge_setup_sympy()
% GDSGE_SETUP_SYMPY  One-time setup for the optional SymPy analytic-Jacobian
%   backend. Installs the `uv` package manager if it is missing, then creates
%   the pyext/.venv environment (a managed Python 3.12 + SymPy). Works on
%   Windows, macOS and Linux. The default adept-autodiff path needs none of this.
%
%   Usage (run once):
%       addpath('src');
%       gdsge_setup_sympy
%
%   Afterwards the SymPy backend is auto-selected when you run a model; force it
%   per-model with `UseAutoDiff=0;` (or `UseAutoDiff=1;` to force adept).

here     = fileparts(mfilename('fullpath'));   % src
repoRoot = fileparts(here);
pyextDir = fullfile(repoRoot, 'pyext');
if exist(pyextDir, 'dir') ~= 7
    error('gdsge:setup:noPyext', 'pyext/ not found next to src/ (looked in %s).', pyextDir);
end

% 1) Find uv; install it if absent.
uv = locateUv();
if isempty(uv)
    fprintf('[gdsge] uv not found — installing it ...\n');
    installUv();
    uv = locateUv();
    if isempty(uv)
        error('gdsge:setup:uvMissing', ...
            ['Could not install uv automatically. Install it from ' ...
             'https://docs.astral.sh/uv/ , open a NEW MATLAB session, then re-run ' ...
             'gdsge_setup_sympy.']);
    end
end
fprintf('[gdsge] using uv: %s\n', uv);

% 2) Create the environment. A virtual environment always needs symlinks (the
% interpreter, plus a lib64 symlink on Linux — even `venv --copies` makes it),
% which fail on filesystems without symlink support (MATLAB Online's
% "/MATLAB Drive": os error 38). So sync into pyext/.venv first; if symlinks are
% unsupported there, relocate the venv to a symlink-capable local directory via
% UV_PROJECT_ENVIRONMENT (gdsge.codegen.sympy.ensurePyenv looks in both places).
setenv('UV_LINK_MODE', 'copy');
fprintf('[gdsge] creating the SymPy environment (uv sync) ...\n');
[status, out] = system(sprintf('"%s" sync --project "%s"', uv, pyextDir));
fprintf('%s\n', strtrim(out));
if status ~= 0
    if symlinkUnsupported(out)
        localVenv = gdsge.codegen.sympy.localVenvDir();
        fprintf(['[gdsge] this filesystem does not support symlinks (cloud/network drive) — '...
                 'placing the environment at\n        %s\n'], localVenv);
        if exist(localVenv, 'dir') == 7; rmdir(localVenv, 's'); end
        broken = fullfile(pyextDir, '.venv');               % clear the half-made venv
        if exist(broken, 'dir') == 7; rmdir(broken, 's'); end
        setenv('UV_PROJECT_ENVIRONMENT', localVenv);
        [status, out] = system(sprintf('"%s" sync --project "%s"', uv, pyextDir));
        fprintf('%s\n', strtrim(out));
        if status ~= 0
            error('gdsge:setup:uvSyncFailed', 'uv sync (relocated) failed (exit %d):\n%s', status, out);
        end
    else
        error('gdsge:setup:uvSyncFailed', 'uv sync failed (exit %d):\n%s', status, out);
    end
end

% 3) Confirm the interpreter exists (same resolver ensurePyenv uses).
venvPy = gdsge.codegen.sympy.venvPython(pyextDir);
if isempty(venvPy)
    error('gdsge:setup:noVenv', 'uv sync ran but no interpreter was found.');
end
fprintf(['[gdsge] SymPy backend ready. It is auto-selected when you run a model;\n' ...
         '        force it per-model with UseAutoDiff=0; (UseAutoDiff=1; forces adept).\n']);
end

% ---------------------------------------------------------------------------
function uv = locateUv()
% Return a usable uv path, or '' if none. uv does not refresh the current
% process PATH after install, so also probe its standard install locations.
uv = '';
exe = 'uv'; if ispc; exe = 'uv.exe'; end
[st, ~] = system([exe ' --version']);
if st == 0; uv = exe; return; end          % on PATH

if ispc; home = getenv('USERPROFILE'); else; home = getenv('HOME'); end
cands = { fullfile(home, '.local', 'bin', exe), ...
          fullfile(home, '.cargo', 'bin', exe) };
for i = 1:numel(cands)
    if exist(cands{i}, 'file') == 2; uv = cands{i}; return; end
end
end

% ---------------------------------------------------------------------------
function installUv()
if ispc
    cmd = 'powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"';
else
    cmd = 'curl -LsSf https://astral.sh/uv/install.sh | sh';
end
if system(cmd) ~= 0
    warning('gdsge:setup:installReturnedNonzero', ...
        'The uv install command returned a non-zero exit code; continuing to probe for uv.');
end
end

% ---------------------------------------------------------------------------
function tf = symlinkUnsupported(out)
% True when uv's failure output indicates the target filesystem cannot create
% symlinks (e.g. MATLAB Online's /MATLAB Drive returns ENOSYS / os error 38).
tf = contains(out, 'symlink') || contains(out, 'os error 38') || ...
     contains(lower(out), 'not implemented');
end
