# Phase 7b Design — ASG Widening (CaoKS2016, Bianchi2011_asg)

- **Date:** 2026-06-13
- **Status:** Approved (design); pending implementation
- **Parent spec:** `2026-06-11-refactor-gdsge-design.md` (§15, Phase 7); split decided in
  `2026-06-12-phase7a-spline-widening-design.md` §2
- **Owner:** Wenlan Luo

---

## 1. Context

Phase 7a closed the spline path: all four spline-based corpus models run green end-to-end
through the public API against fresh old-toolbox goldens. The two remaining corpus models
(CaoKS2016, Bianchi2011_asg) set `USE_ASG=1; USE_SPLINE=0;` and use the **adaptive sparse
grid** architecture, which differs structurally from the spline path:

- An adaptive grid-refinement loop runs *inside* each fixed-point iteration: propose grids
  level by level (`get_eval_grids`), solve the model at the proposed points via the MEX,
  push results back (`push_eval_results`), refine until `AsgMaxLevel`/`AsgThreshold`.
- A public MATLAB class `asg` (over `asg_mex`) implements the interpolant. **Both test
  drivers call it directly** (`asg.construct_from_struct`, `get_grids_info`, `eval`,
  `eval_vec`) — it is frozen public surface, not an internal detail.
- `IterRslt` carries new fields: `asg_interp_struct` (var_interp interpolant),
  `sol_asg_interp_struct` (policy solutions, used to warm-start simulate's re-solves), and
  `asg_output_struct` (var_output interpolant) + `output_var_index`.
- The model MEX evaluates the previous iteration's interpolant inside equations through a
  **class handle** (uint64 → `AsgInterpArray*`) with an adept-aware evaluator
  (`asg_adouble.h`), instead of the spline path's `GDSGE_PP_*` caller-workspace structs.
- Bianchi's driver additionally exercises **runtime shock overrides**
  (`options.shock_trans/yT/yN` replace gmod placeholder values), `MaxIter` staging,
  state-bound widening (`options.b = [-1.1, 0]` plus `bMin`/`bMax` setup-variable
  overrides), `SkipModelInit=1`, and ASG `WarmUp`.

Old-toolbox reference points (all under `base_package/gdsge/source/`): templates
`iter_inf_horizon_asg_template.m`, `iter_init_asg_template.m`,
`iter_solve_problem_asg_template.m`, `iter_solve_init_problem_asg_template.m`,
`asg_propose_grids_and_solve_template.m`, `iter_inf_output_construct_asg_template.m`,
`simulate_resolve_asg_template.m`; C++ fragments
`interp_asg_{construct,get,prepare_space}_template.cpp`; runtime `asg.m`, `asg_mex.cpp`,
`include/asg.h`, `include/asg_adouble.h`, `include/asg_matlab.h`, `class_handle.hpp`;
ASG logic in `gdsge_parser.m` (option validation ~line 1005, template selection ~1412,
simulate selection ~1738, compile defines ~1847).

Groundwork already in place in the new package: the IR schema defines
`interpMethod ∈ {spline, asg, pchip}` plus `asgMaxLevel`/`asgThreshold`; the validator
enforces asg ⊥ `var_tensor`; `generateCxx` raises a deliberate "Phase 7b" error for
non-spline methods; `model_init`, `var_policy_init`/`var_aux_init`/`inbound_init`,
`var_output`, `SkipModelInit`, and state-grid overrides landed in 7a.

## 2. Scope decisions

Two cuts, both approved by the owner (2026-06-13):

- **`var_tensor` moves to 7c.** No corpus model uses it (CaoKS2016 has the block commented
  out), and the old parser hard-errors when it is combined with ASG
  (`gdsge_parser.m`: "cannot use var_tensor when using adpative grids"). 7c is the
  golden-less-surface phase with synthetic fixtures — exactly its charter. `PROGRESS.md`
  is updated accordingly.
- **ASG + `SIMU_INTERP` moves to 7c.** The old toolbox supports simulate-by-evaluating
  `asg_output_struct`, but neither ASG driver uses it (both use the default
  `SIMU_RESOLVE`). 7b raises an informative error if a gmod combines `USE_ASG` with
  `SIMU_INTERP`.

Also out of scope: `GDSGE_ASG_FIX_GRID` and other ASG option branches no driver exercises
(each gets an informative error or a faithfully replicated default — decided per-option in
the plan); PCHIP; macros/hooks/`cxx` blocks (7c); SymPy backend (Phase 8).

## 3. Goal & success criteria

Phase 7b is done when, for each model, a Slow-tagged gate drives the public API exactly as
the model's `test.m` does and matches a freshly captured golden:

1. `gdsge_codegen('<model>')` honors the §3.1 artifact contract from the 7a gates
   (`iter_*.m`, `simulate_*.m`, `mex_*.cpp`, `<model>.gdsge.json` schema-valid and
   round-tripping, cache-gated MEX recompilation) — now with `interpMethod = 'asg'`.
2. **CaoKS2016**: `iter_CaoKS2016(options)` (PrintFreq/SaveFreq overrides) converges at
   `TolEq=1e-4` with iteration count and metric matching the golden; the seeded reduced
   simulate matches (shock path bit-exact, fields within tolerance); the driver's
   post-processing works against the new results —
   `asg.construct_from_struct(IterRslt.asg_interp_struct)`, `get_grids_info`, `eval_vec`.
3. **Bianchi2011_asg**: the full driver sequence reproduces the golden — stage-1 iter with
   runtime shock overrides (`options.shock_trans/yT/yN` from `shock_process.mat`) and
   `MaxIter=50`; stage-2 re-solve to convergence with `options.WarmUp`, widened bounds
   (`options.b=[-1.1,0]`, `bMin`/`bMax`), `SkipModelInit=1`; seeded reduced simulate;
   driver post-processing on `asg_output_struct` (`get_grids_info`, per-shock `eval`).
4. All HL1996 and 7a gates stay green throughout.

## 4. Approach — 7a's hybrid, adapted

Considered: (a) model-by-model verticals — would design the shared ASG machinery
mid-stream; (b) layer-by-layer across both models — end-to-end feedback arrives last.
**Chosen (c), the 7a recipe:**

1. **Capture both goldens** in one old-toolbox session (§9).
2. **Vendor the public ASG runtime** (§5) and settle the small front-end/IR delta (§6).
3. **Vertical per model, CaoKS2016 first** — its vanilla driver validates the core ASG
   iteration loop in isolation; Bianchi then layers the options/WarmUp/override surface on
   a proven loop.

## 5. Vendored public ASG runtime

- `asg.m` lands **flat** in `src/` beside `gdsge_codegen.m` — it is frozen public surface
  (both drivers construct and query it), matching the parent spec's policy that deliberate
  backward-compat shims sit at top level.
- `asg_mex.cpp` and the headers it needs (`asg.h`, `asg_adouble.h`, `asg_matlab.h`,
  `class_handle.hpp`) join the vendored kernel/include tree **as-is** (parent spec:
  numerical kernels are moved, not rewritten).
- **Auto-compile, cache-gated:** when the IR says `interpMethod='asg'`,
  `gdsge.codegen.codegen` ensures an up-to-date `asg_mex` binary exists before model
  compilation — gated by a source+header hash cache (same mechanism as
  `mex_<model>.cache`). This is needed anyway because codegen reads
  `asg.get_mex_constants()` to obtain the `-DASG_MAX_DIM/NVEC/LEVEL` compile defines.
  (Deliberate improvement over the old toolbox's ship-prebuilt-binaries + manual
  `compile_asg.m`; behavior change is invisible to callers.)

## 6. Front-end / IR delta (small)

- **Options:** `resolveOptions` + schema gain `AsgMinLevel` (default 4),
  `AsgOutputMaxLevel` (10), `AsgOutputThreshold` (1e-2) alongside the existing
  `asgMaxLevel`/`asgThreshold`; the USE_ASG ⊕ USE_SPLINE ⊕ USE_PCHIP exclusivity and the
  asg ⊥ var_tensor invariant are already enforced. ASG + `SIMU_INTERP` becomes a parser
  error (§2).
- **Confirm-through-gates (no new machinery expected):** Bianchi's degenerate init system
  (`equations; 0; end;` — a bare numeric-literal equation balancing the `dummy` policy)
  and CaoKS's named-transition reduction pipe `GDSGE_EXPECT{ … | shock_trans2 }` where the
  transition is a declared *parameter*. Front-end gates assert both.
- **Setup-assigned names:** the IR setup section records the names assigned by the gmod's
  setup code (the parser already evals it), so the runtime whitelist (§8) can accept
  overrides like Bianchi's `bMin`/`bMax`.

## 7. Codegen

### 7.1 MATLAB emitters (`gdsge.codegen.mat`, selected on `interpMethod='asg'`)

- **ASG iteration loop** mirroring `iter_inf_horizon_asg_template.m` +
  `asg_propose_grids_and_solve_template.m`: each fixed-point iteration builds the next
  interpolant by level-by-level refinement — `get_eval_grids(threshold)` (threshold 0
  below `AsgMinLevel`, else `AsgThreshold`) → pack states/shocks for the proposed points →
  MEX solve → `push_eval_results` — then metric via `asg.compute_inf_metric` against the
  previous interpolant; `sol_asg_interp_struct` accumulates policy solutions.
- **ASG init segment** (`SkipModelInit`-gated) reusing the `initView` emitter pattern from
  7a, mirroring `iter_init_asg_template.m`: solve the init problem on ASG-proposed grids
  to seed `var_interp`.
- **Output construction** mirroring `iter_inf_output_construct_asg_template.m`: refine a
  fresh ASG over `var_output` up to `AsgOutputMaxLevel`/`AsgOutputThreshold`.
- **ASG `SIMU_RESOLVE` simulate** mirroring `simulate_resolve_asg_template.m`: reconstruct
  `GDSGE_ASG_INTERP` and `GDSGE_SOL_ASG_INTERP` from the structs; per period, `eval_vec`
  the solution interpolant at current states to warm-start, then MEX re-solve.
- Generated files stay **thin** over new `gdsge.runtime` ASG helpers (refinement-loop
  driver, warm-up reconstruction, metric) — Phase 4 philosophy: explicit locals, no
  v2struct, structured progress/unconverged reporting reused.

### 7.2 C++ emitters / templates

- New `templates/cxx/` fragments mirroring `interp_asg_{construct,get,prepare_space}`:
  the MEX receives the previous interpolant as a **class handle** and evaluates it inside
  equations through `AsgInterpArrayAdoubleEvaluator` (adept-aware), with stack buffers
  sized by `ASG_MAX_*` defines.
- The caller-workspace contract widening (interp handle + ASG buffers replacing
  `GDSGE_PP_*` splines) flows through the shared `gdsge.codegen.dataLayout` descriptor so
  the MATLAB packer and the C++ `POP*` macros stay generated from one source.
- Compile step gains `-DUSE_ASG -DASG_MAX_DIM=<num_state> -DASG_MAX_NVEC=<num_interp>
  -DASG_MAX_LEVEL=<from asg.get_mex_constants()>` (old `gdsge_parser.m` ~line 1847
  convention).

## 8. Runtime options widening

`gdsge.runtime.unpackOptions`' whitelist grows, preserving informative errors for true
unknowns and the **old ordering** (overrides applied after the setup-code replay):

- `shock_trans` and declared `var_shock` names (Bianchi replaces gmod placeholders at
  runtime);
- setup-assigned variable names from the IR setup section (§6) — covers `bMin`/`bMax`;
- `MaxIter` (finite stage-1 → `inf` stage-2 staging);
- ASG `WarmUp`: reconstruct from `WarmUp.asg_interp_struct` /
  `WarmUp.sol_asg_interp_struct` via `asg.construct_from_struct`, resume `WarmUp.Iter` —
  exact old-template semantics, including interaction with `SkipModelInit=1` and
  state-bound overrides (for ASG only state min/max matter; the interpolant carries its
  own grids).

## 9. Result-struct contract (frozen)

`IterRslt` gains `asg_interp_struct`, `sol_asg_interp_struct`, `asg_output_struct` (when
`var_output` is present), and `output_var_index`, alongside the existing
`Metric`/`MetricVec`/`Iter`/`shock_num`/`shock_trans`/`var_shock`/`var_state`/`params`
fields. The `asg` struct layout is whatever `asg.convert_to_struct` emits (vendored
unmodified — `numDim`, `numVec`, `numArray`, `inputGrids`, `stateMin/Max/Range`,
`currentLevel`, `gridsCurrentLevel`, `gridsNextLevel`, `unscaledGrids`, `surplus`,
`levels`). `SimuRslt` shapes are unchanged from the spline path.

## 10. Golden capture protocol

One session; each model in its **own `matlab -batch` process** with controlled cwd and
only `base_package/gdsge/source` on the path. Committed capture scripts record exact
settings and log wall-clock timings.

- **CaoKS2016:** full-convergence iter as the driver runs it (`PrintFreq=10`,
  `SaveFreq=100`, gmod `TolEq=1e-4`, `AsgMaxLevel=10`, `AsgThreshold=1e-4`); simulate
  reduced to a seeded 6 samples × 1,000 periods (7a sizes; gmod default is 100×10,000 —
  far too large for a committed golden given per-period re-solves).
- **Bianchi2011_asg:** the **full two-stage driver sequence captured as-is** — stage-1
  iter (`MaxIter=50`, shock overrides from `shock_process.mat`), stage-2 WarmUp re-solve
  to convergence with widened bounds and `SkipModelInit=1` — because the staging itself is
  surface under test. `shock_process.mat` is copied into the new test fixture. Simulate:
  seeded 6×1,000.
- Goldens land under `tests/<Model>/golden/` with a golden-integrity test each. ASG
  goldens should be far smaller than Mendoza's 99.8 MB outlier; sizes are checked before
  commit. If a full-convergence run proves prohibitively long, pause and revisit the
  protocol for that model before building gates on it.

## 11. Testing

- Harness unchanged: headless `matlab.unittest` via `matlab -batch "cd('tests');
  run_tests"`; `junit.xml` + exit code authoritative; one-source-per-process path policy.
- Per model: golden-integrity test; front-end gate (validator pass, JSON round-trip,
  structural assertions incl. §6 confirmations; IR JSON committed as snapshot once the
  end-to-end gate is green); Slow end-to-end gate replicating the driver (§3).
- Unit level: vendored `asg` class smoke (construct → push → eval → struct round-trip on
  a known function); asg_mex cache-gated auto-compile (compiles once, skips when fresh);
  option-resolution and whitelist-widening tests (accept shock/setup overrides, error on
  unknowns, ASG+SIMU_INTERP parser error); MATLAB/C++ emitter snapshots.
- Existing HL1996 + 7a gates stay green throughout.
- Tolerance discipline as established: seeded shock paths bit-exact; fields within the
  comparison utility's bounds; iteration count and metric matched to the golden.

## 12. Risks

- **Handle-based MEX interop** (class_handle + adept-aware evaluator) is the deepest new
  C++ contract. Mitigation: headers vendored unmodified; CaoKS golden gate pins behavior
  before Bianchi's option surface is added.
- **Unknown ASG wall-clock** (Bianchi: `AsgMaxLevel=14`, 16 shocks; CaoKS simulate
  re-solves 1,000 periods × 6 samples). Capture scripts log timings; gates Slow-tagged;
  protocol revisited if prohibitive (§10).
- **Setup-variable override surface** could over-widen the whitelist. Constrained to
  names actually assigned in the gmod's setup code; everything else still errors
  informatively.

## 13. Out of scope

`var_tensor` and ASG+`SIMU_INTERP` (7c, per §2); `GDSGE_ASG_FIX_GRID` and unexercised ASG
option branches (informative error or replicated default, per-option in the plan); PCHIP;
macro engine, `gdsge.m`, `GenCodeSegment`, legacy 5-output signature, hook/`cxx` blocks
(7c); SymPy backend (Phase 8); performance work beyond parity (Phase 9).
