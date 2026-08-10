# Global var_aux are constructed outside the MEX and feed back only through var_interp

A global var_aux (e.g. `MPK = GDSGE_DER(Y, K);`) is a function of the whole solved
policy, so it cannot be computed state-by-state inside the model MEX. We decided:
(1) globals are constructed in the generated `iter_<model>.m` every iteration —
after the solution tensors are unpacked, before the var_interp update expressions
run — using the existing C++ spline-construction and analytic-gradient kernels;
(2) they are **excluded from the `GDSGE_AUX` slot layout** (the MEX ABI never sees
them), tracked in the IR as `variables.auxGlobal` but merged into `IterRslt.var_aux`
so users see one uniform surface; (3) they enter next-iteration equations **only**
when the user explicitly routes them through a var_interp update
(`MPK_interp = MPK;`) — there is no implicit feedback.

## Considered Options

- Teaching the MEX about globals (a second solve pass or in-MEX spline eval):
  rejected — it would grow the task-MEX ABI and duplicate spline construction that
  already lives in MATLAB-callable kernels.
- A dedicated "derivative interp" object (implicit feedback): rejected — the
  explicit `var_interp` update reuses existing convergence/metric semantics and
  keeps the user in control of what iterates.

## Consequences

- Globals must not be assigned in the model block (parse error), and ordinary
  var_aux gain the previously-missing "must be assigned in the model block" check.
- Simulate can never re-solve a global; it interpolates them via
  `IterRslt.output_interp` in both simulate modes (with a codegen warning under
  SIMU_RESOLVE).
- ASG grids and `var_aux_init` are unsupported (namespaced errors) until the ASG
  MEX exposes gradient evaluation.
