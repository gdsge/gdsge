# pyext — SymPy analytic-Jacobian backend (alternative path only)

This uv-managed Python environment exists **only** for the alternative symbolic backend
(spec §9.3, Phase 8): it turns the IR's expression ASTs into a stack-allocated analytic
Jacobian via SymPy + `cse()` + `ccode()`. The **default** GDSGE path (MATLAB parser + adept
autodiff MEX) needs no Python at all.

## Usage

```
uv sync --project pyext           # create the env (downloads a managed Python 3.12)
uv run --project pyext pytest pyext/tests
```

`pyext/.venv/`, `.python-version`, and `__pycache__/` are git-ignored; `pyproject.toml` and
`uv.lock` are tracked.
