# Global var_aux and the GDSGE_DER operator

Approved design from the 2026-08-09 grilling session (task20260809_extend_global_aux.md).
See also `docs/adr/0001-global-var-aux-outside-the-mex.md` and
`docs/adr/0002-always-generate-both-simulate-variants.md`, and the glossary in
`CONTEXT.md`.

## Feature

A new class of auxiliary variable — **global var_aux** — declared as a plain
`var_aux` but *defined* by an operator over the whole solved function instead of
state-by-state in the model block. First operator:

```
var_aux Y MPK_analytical MPK;
MPK = GDSGE_DER(Y, K);
```

`GDSGE_DER(v, s)` = the analytic partial derivative w.r.t. var_state `s` of the
full tensor spline fitted to `v` over all continuous states (per shock), evaluated
at every grid point. Result has the same tensor shape as any var_aux.

## Semantics (all settled)

- **First argument `v`**: any var_policy, any ordinary var_aux, or any
  *earlier-defined* global var_aux (chaining ⇒ second derivatives allowed;
  accuracy caveat documented, not policed).
- **Second argument `s`**: any var_state.
- **Spline spec**: shares the model's var_interp configuration
  (`INTERP_ORDER`/`EXTRAP_ORDER`, default cubic not-a-knot, linear extrap).
  Extrapolation is observably irrelevant in iter (evaluation at grid points).
- **Cadence**: computed **every iteration**, in generated `iter_<model>.m`, after
  `unpack`/`reshapeShock` and **before** the var_interp update expressions — so
  `MPK_interp = MPK;` updates work. Same treatment at init/WarmUp construction
  sites wherever update expressions run.
- **Feedback**: only via explicit var_interp updates (ADR-0001). No implicit use
  in the model block.
- **Kernel route**: MATLAB-side call into the existing C++ machinery —
  `interp_construct_mex` for construction; a new MATLAB-facing gradient-eval
  kernel over `InterpEval::search_eval_with_grad_vec` /
  `MatlabInterpEval::search_eval_with_grad_vec_at_array` (runtime XDIM dispatch
  like `interp_eval_mex.cpp`), since `interp_eval_mex` exposes values only.

## Syntax & errors

- Definition statements live in the **var_aux definition section** — the
  declaration-region section following a `var_aux` line (Pass B/C name-routing in
  `parseVarDecls.m`, like var_interp updates / var_tensor assignments, but with
  **enforced placement**: a global definition outside a var_aux section is an
  error, via the `expectedSection` machinery).
- New namespaced errors (register in `docs/deferred-features.md` where deferred):
  - global var_aux assigned in the model block;
  - ordinary var_aux **never** assigned in the model block (new check, replacing
    today's opaque C++ "undeclared identifier" failure);
  - `GDSGE_DER` on `var_aux_init` / in `model_init` context;
  - global var_aux in an ASG model (`interpMethod ~= 'spline'`) — until
    `asg_mex` exposes `eval_vec_with_grad`;
  - unknown names in `GDSGE_DER(v, s)`.
- Users need not know the distinction; the toolbox tracks it.

## IR & result surfaces

- New optional IR list `variables.auxGlobal`: `{name, op: "der", of, wrt}` in
  definition order. **Excluded** from `variables.aux` slot layout (`GDSGE_AUX`,
  `NUM_AUX`, `emitAux`/`emitDeclare` untouched by globals).
- Minor IR version bump; `schema.m` + `validate.m` (name resolution: `of` ∈
  policy ∪ aux ∪ earlier auxGlobal; `wrt` ∈ states) + `ir-schema.md` regen +
  model JSON snapshot regen (known fan-out checklist).
- `IterRslt.var_aux.<name>` carries globals alongside ordinary aux (frozen shape
  preserved; globals are additive fields). Globals participate in
  `var_output`/`output_interp` stacking like any aux.

## Simulate

- Globals cannot be re-solved; both modes interpolate them via
  `IterRslt.output_interp` (`OutputInterpOrder`, runtime-dispatched), so both
  modes agree on simulated values.
- `SIMU_RESOLVE` + global in `var_simu` ⇒ **codegen-stage warning** that these
  values are interpolated; the resolve-variant simulate file embeds the
  interp-eval for exactly those aux.
- `SIMU_INTERP` non-expressible fallback (slow MATLAB loop) now also **warns at
  codegen** (was silent).
- `CONSTRUCT_OUTPUT=0` + globals in `var_simu` ⇒ clear namespaced error at
  simulate entry naming `CONSTRUCT_OUTPUT` (never silently override the user).
- **Always-generate restructuring (ADR-0002)**: `isSimuMexExpressible` drops its
  `options.simuInterp == 1` conjunct; `generateCxx` always emits
  `simulate_<model>_mex.cpp` for expressible spline models; single
  `simulate_<model>.m` remains, rewritten per toggle; compile cache becomes
  toggle-invariant with zero cache-code changes.
- Latent pre-existing bug (state name in `var_simu` never promoted to
  `var_output` → grid/vector dim mismatch in the fallback loop): **logged in
  `docs/deferred-features.md`, not fixed here.**

## rbc example & tests (`tests/rbc/`, standard 5-part layout)

Adapted from `base_package/gdsge/tests/rbc/rbc_simple.gmod`:

- Add ordinary aux `Y = z*K^alpha;` and `MPK_analytical = alpha*z*K^(alpha-1);`
  (exact, solution-independent); global `MPK = GDSGE_DER(Y, K);`; keep `w`; add
  `MPK` to `var_simu` (exercises the resolve-mode interpolation warning path).
- **Goldens from the new pipeline** (old toolbox cannot parse `GDSGE_DER`;
  precedent: `CaoKS2016_simu_interp`). Correctness gate = self-validating
  `MPK ≈ MPK_analytical` at **all** grid points (boundaries included), one
  empirically calibrated tolerance documented in the test.
- Second, feedback-variant gmod in `tests/rbc/` (declares `MPK_interp` with
  update `MPK_interp = MPK;`) run end-to-end — loop ordering is untestable at
  unit level.
- Plus parser unit tests (placement enforcement, all new errors), emit tests,
  kernel gradient test, and snapshot/artifact updates across the corpus for the
  always-generate change. TDD throughout; end-to-end gates run sequentially.
