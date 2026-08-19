# GDSGE Toolbox

A MATLAB toolbox that solves global DSGE economic models: a `.gmod` model file is
parsed into a JSON IR, from which MATLAB and C++ (MEX) solver code is generated.

## Language

### Model structure

**Model region**:
One `model;` or `model(<guard>);` block. A model may declare several regions; at each
grid point the toolbox runs the first region whose guard holds. Regions partition the
state space — they do not compose.
_Avoid_: conditional model, model branch

**Conditional statement**:
An `if/elseif/else/end` block among the assignments of a model body, selecting which
assignments run at a grid point. Every case is part of the same region and fills the
same equation slots.
_Avoid_: model-body conditional, if block

**Conditional equation**:
An `if/elseif/else/end` block inside `equations;…end;`, selecting which residual
expressions fill a slot range. Every case must fill the same number of slots.
_Avoid_: equation branch, conditional residual

**Guard**:
The boolean expression that selects a case of a model region, a conditional statement,
or a conditional equation. A guard may reference only var_state, var_shock, and
parameters — never a var_policy or a value computed in the model body.
_Avoid_: condition, branch condition, predicate

### Variables

**Ordinary var_aux**:
An auxiliary variable declared in `var_aux` and assigned point-by-point inside the
model block; the model MEX computes it state-by-state alongside the policies.
_Avoid_: local var_aux, pointwise var_aux

**Global var_aux**:
An auxiliary variable declared in `var_aux` but *defined* by an operator (e.g.
`GDSGE_DER`) in the var_aux declaration section, constructed from the whole solved
function after each iteration rather than state-by-state. Users see it as a plain
var_aux; the toolbox tracks it separately (IR: `variables.auxGlobal`). It must NOT
be assigned in the model block.
_Avoid_: derived var_aux, post-iteration var_aux

**GDSGE_DER**:
The global-var_aux operator `GDSGE_DER(v, s)`: the analytical first derivative,
with respect to var_state `s`, of the tensor spline fitted to variable `v` (a
var_policy, an ordinary var_aux, or an earlier global var_aux), evaluated at every
grid point, per shock. The spline uses the model's var_interp specification.

**Definition section**:
The region of the declaration block following a `var_aux` line where global
var_aux definitions (`name = GDSGE_DER(...)`) must appear. Placement elsewhere is
an error (stricter than var_interp updates, which tolerate any position).

### Simulation

**Resolve simulate (SIMU_RESOLVE)**:
The simulate mode that re-solves the model's nonlinear system each period, so
policies and ordinary var_aux are exact at the simulated point. Global var_aux
cannot be re-solved and are interpolated even in this mode (with a warning).

**Interp simulate (SIMU_INTERP)**:
The simulate mode that never re-solves: every simulated variable is read off the
output interpolant. The two modes are mutually exclusive per model run.

### Feedback

**Feedback via var_interp**:
The only way a global var_aux enters the next iteration's equations: the user
declares a var_interp and assigns the global var_aux in its update
(e.g. `MPK_interp = MPK;`). Global var_aux values are therefore computed every
iteration, before the var_interp update expressions run.
