# Notes for agents — facts about the old codebase

Concrete findings from exploring `base_package/` so you don't re-derive them. Paths are
relative to the repo root. `base_package/` is git-ignored but present in the working dir.

## Old toolbox map (`base_package/gdsge/source/`)

- `gdsge.m` — top orchestrator: `gdsge_codegen` → solve (`iter_<model>`) → `simulate_<model>`,
  with caching of generated code and `IterRslt`.
- `gdsge_codegen.m` — writes generated `iter_*.m`, `simulate_*.m`, `mex_*.cpp`,
  `compile_*.m`; triggers MEX compilation when the C++ changed.
- `gdsge_parser.m` — **the 3027-line monolith** to replace. Stages, in order: macro
  expansion (`#define/#for/#foreach/#if/#mat/#strcat_comma/include/cinclude`) → block
  extraction (`pre_model`, `model`, `model_init`, `simulate`, `pre_iter`, `post_iter`,
  `pre_jac_code`, `post_jac_code`, …) → variable-declaration parsing (with flat slot
  indexing) → parameter `eval` → interpolation/tensor setup → model-body translation to C++
  (`gen_cxx_model_code`, `extract_reduction_operator`) → template substitution. It leans on
  fixed regex + `str2sym`/`ccode` and word-boundary `my_regexprep`.
- `source/code_template/` — the `.m` and `.cpp`/`.h` templates filled by placeholder
  substitution. Key ones: `iter_inf_horizon_template.m` (main loop; uses `v2struct`),
  `mex_template.cpp`, `task_template.cpp` (OpenMP per-grid-point; `POPN`/`POPNARRAY` stack
  locals; `Stack` adept), `model_template.cpp` (`adouble` model fn), `call_fmin_template.cpp`
  (`CoDoSol::solve`).
- `source/include/` — vendored C++: `adept*` (autodiff), `Eigen/`, `codosol.h`
  (constrained dogleg solver), `asg*.h` (adaptive sparse grids), `cubic_spline_notaknot.h`.

## Runtime shape to preserve

- Each (state × shock) grid point is an independent nonlinear solve, run in parallel under
  OpenMP. Inputs are "popped" from a flat `GDSGE_DATA` array into stack locals via
  `POPN(var)` / `POPNARRAY(var,len)`. `CoDoSol::solve` drives the solve; the model function
  returns residuals and (today) an adept-recorded Jacobian. Stack allocation throughout
  (e.g. `double GDSGE_EQ[NUM_EQUATIONS]` when small).
- `v2struct` is used to pack/unpack many workspace variables in the generated `.m`. The
  refactor drops it (copies + silent overwrite) in favor of explicit named locals and direct
  struct dot-assignment.

Model expressions support relational operators on the adept path (`emitExpr` two-tier
precedence: `==,~=` below `<,>,<=,>=` below arithmetic, matching C++; `~=` → `!=`).
Logical operators (`&`, `|`, `&&`, `||`) parse into binop nodes but the emitter rejects
them (`unsupported binop`). The SymPy backend does not support relationals — auto mode
guards on this (resolveBackend), but an explicit sympy pin fails mid-codegen with a raw
Python error (known, deferred).

## Reference for the analytic-Jacobian backend (`base_package/hans/`)

- `hans_parser.py` — SymPy-based: defines custom `Function`s with `fdiff`, computes
  gradients/Hessians via `diff`, runs `cse(...)` to extract repeated subexpressions, and
  emits stack-allocated `double helper_i;` via `ccode()`. See `template/vfi_template_sympy.h`
  for the C++ wrapper shape.
- It uses **Lark** grammars (`hans/lark/*.lark`) to parse `.hmod`. Note: it treats `EXPECT`
  as the **identity** and aggregates expectations in MATLAB. GDSGE instead needs **true
  loop fusion** of reductions with the symbolic gradient — see spec §10.
- `py_package/{sympy,lark,mpmath}` are vendored deps in hans — ignore their internals; the
  new package uses uv-managed SymPy under `pyext/`.

## The gmod language is a MATLAB superset

Lines like `k = exp(linspace(log(kMin+1e-3),log(kMax+1e-3),kPts)) - 1e-3;` and
`shock_trans = shock_trans ./ repmat(sum(shock_trans,2),[1,shock_num]);` are real MATLAB,
`eval`'d during parsing to obtain grids/params. This is why the parser must run in MATLAB.

The declaration region accepts *arbitrary* MATLAB, not just `name = value`
assignments: multi-output calls (`[a,b] = ndgrid(...)`), bare calls (`rng(0)`),
struct-field LHS, and `for`/`if` blocks all pass through verbatim
(`parseVarDecls` Pass B), and declaration kinds may re-open
(`parameters ... var_shock ... parameters ...`, Pass C — `ir.setup` kinds can
repeat). Passthrough code runs twice: once in the parse-time eval and again in
the generated iter/simulate files after the default flags. Two guardrails
replace the old toolbox's blind copying: a command-form statement whose first
word is edit-distance ≤ 2 from a declaration keyword raises
`gdsge:parser:probableTypo`, and broken MATLAB still fails as
`gdsge:parser:setupEvalFailed` with numbered lines. Limitation: state grids and
interp updates are only captured from single `name = rhs` assignments.

## Test models (`base_package/gdsge/tests/`)

Folders: `HeatonLucas1996`, `Bianchi2011_asg`, `Mendoza2010`, `CaoKS2016`, `GLSW2020`,
`Barro_et_al_2017`. Each has `<model>.gmod` + a `test.m`. The old driver is
`base_package/gdsge/tests/runtests.m`. (Huggett1997 was removed by the user.)

### HeatonLucas1996 — the first vertical slice

`HL1996.gmod`: 1 state `w1` (wealth share, 201 pts on `[-0.05,1.05]`), 8 shocks
(`g,d,eta1`) with a full 8×8 `shock_trans`. 12 policy vars incl. array `w1n[8]` (one future
wealth share per shock). Uses: the `'` future marker, `GDSGE_INTERP_VEC'(w1n')` (primed
vector interp), `GDSGE_EXPECT{...}` reductions over the 8 shocks, a **future-indexed
equation** `w1_consis'` (expands to 8 equations), `adaptive(1.5)` bounds, KKT
complementarity equations. Square 19×19 system. Converges ≈ iteration 209. `test.m` runs
`rng(0823)` before `simulate_HL1996`. No ASG, no `var_tensor`, no `model_init`, no
`pre/post_iter`, no macros, no `cxx` block — a clean "core" model.

## Full gmod feature inventory (backward-compat surface)

Blocks: `parameters`, `var_state`, `var_shock`, `var_policy[/_init]`, `var_aux[/_init]`,
`var_interp`, `var_tensor`, `var_output`, `var_others`, `model[/_init]`, `equations`,
`simulate`, `pre_model`, `pre_iter`, `post_iter`, `pre_jac_code`, `post_jac_code`.
Declarations: `inbound[/_init]` (+ `adaptive(factor)`), `initial`, `shock_num`,
`shock_trans`. Operators: `GDSGE_EXPECT`, `GDSGE_MIN`, `GDSGE_MAX`, `GDSGE_PROD`,
`GDSGE_INTERP_VEC` (and primed `'`). Markers: trailing `'` = future/next-state; `[N]` =
shock-indexed array var. Built-ins: `GDSGE_Iter`, `TASK`, `OUTPUT_CONSTRUCT_CODE`. Macros:
`#define`, `#for`, `#foreach`, `#if`, `#mat{}`, `#strcat_comma`, `include`, `cinclude`.
Options: `USE_SPLINE/USE_ASG/USE_PCHIP`, `INTERP_ORDER`, `ExtrapOrder`,
`SIMU_RESOLVE/SIMU_INTERP`, `TolEq`, `SaveFreq`, `PrintFreq`, `NumThreads`, `AsgMaxLevel`,
`AsgThreshold`, and more.

## Driving MATLAB headless

`matlab -batch "<expr>"` runs and exits; the exit code reflects uncaught errors or an
explicit `exit(n)`. Use it to isolate old vs new (one source `addpath`ed per process, with a
controlled `cd`).

## Verifying ASG (and other MEX) changes — compile + run SEQUENTIALLY

Each end-to-end gate does its own `gdsge_codegen` → **compiles a model MEX (~30 s)** before it
runs, and MATLAB/OpenMP must run **one process at a time** (each saturates all cores; concurrent
runs contend and distort timing). So a batch of N ASG gates is inherently **minutes** of mostly
sequential compilation — that is **not a hang**.

- **Run gates one at a time (or one `runtests({...})` call) with progress visible.** Filtering the
  output down to only the final `PASS=`/summary hides the per-test `Running…/Done…` and the `Iter:`
  lines, which makes a normal (compiling) run look stuck. Show `Running t`/`Done t` and periodic
  `Iter:` milestones so "compiling" is distinguishable from "looping".
- **Compile-slow vs a real solver hang look different.** Compilation: no `Iter:` output at all,
  high CPU. A genuine non-convergence (e.g. a value-iteration cycle): `Iter:` keeps climbing but the
  **`Metric` does not decrease** (and `maxF`≈0, i.e. each minor solve converges but the fixed point
  doesn't). To probe the latter, run the model's `iter` directly with a **capped `MaxIter`** and
  watch the `Metric` trend instead of letting `MaxIter=inf` spin forever.
- The `asg_mex` kernel build is cache-gated (`ensureAsgMex`, keyed on `asg_mex.cpp` + `asg.h`); the
  **model MEX cache keys only on the model `.cpp`**, so an `asg.h`/`asg_adouble.h` edit does NOT
  auto-rebuild a cached model MEX — gates rebuild fresh in temp dirs, but a manual reuse can go
  stale. See the [[asg-kernel-perf]] memory.

## CoDoSol Broyden mode — dead plumbing, now fixable; profiled net-negative on HL1996

- **The Broyden flag is disconnected.** `call_fmin.tpl.cpp` / `call_fmin_sympy.tpl.cpp` pass
  the static `UseBroyden` (default 0) to `CoDoSol::solve` as `flagBroyden`. The per-iteration
  `GDSGE_USE_BROYDEN_NOW = (GDSGE_Iter>1)*GDSGE_USE_BROYDEN` (`emitIter.m`, with
  `GDSGE_USE_BROYDEN` default 1) is plumbed all the way into the MEX via `GET_INT` but **never
  used** — in both the old and new toolbox. So "Broyden after the first major iteration" has
  never actually run. To enable it, pass `(UseBroyden || GDSGE_USE_BROYDEN_NOW)` to `solve`.
- **`A_plus_x_times_yp` (the rank-1 Broyden Jacobian update) had a loop-variable bug**
  (`for (int j=0; i<n; i++)` — used `i` not `j`); fixed 2026-06-24. Latent because the default
  corpus never exercises Broyden (`UseBroyden=0`).
- **Profiling (HL1996, bug fixed, flag wired): Broyden is *slower*, both backends.** Full
  Newton vs Broyden-after-iter-1: −13%/−15% at the default `broydenthresh=1e-6`, −30%/−44% with
  an aggressive threshold (Broyden every Newton step). All converge identically to golden (209
  iters). Reason: HL1996's per-point solves are well warm-started (few Newton steps), so the
  rank-1 update overhead exceeds the Jacobian-eval savings. **Caveat:** HL1996-specific —
  Broyden's payoff scales with Jacobian-eval cost and Newton-steps-per-solve, so a large,
  cold-started model could still benefit. See [[eigen-openmp-malloc-and-solve-dim]].
