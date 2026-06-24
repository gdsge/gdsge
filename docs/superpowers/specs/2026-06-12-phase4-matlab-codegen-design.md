# Phase 4 — MATLAB Codegen Design

- **Date:** 2026-06-12
- **Status:** Approved (design); pending implementation
- **Parent spec:** `2026-06-11-refactor-gdsge-design.md` (§9.1, §15 Phase 4)
- **Scope:** IR → `iter_<model>.m` / `simulate_<model>.m` for HL1996, plus the
  `gdsge.runtime` library and the vendored MATLAB kernels they depend on.

---

## 1. Goal

Consume the validated HL1996 IR (Phase 3's output) and generate the two runtime
MATLAB files, rewritten per the parent spec: **no `v2struct`**, explicit error
checking, informative reporting, frozen public surface. The phase gate is
**functional**: the new generated files must drive the **old compiled MEX**
(`mex_HL1996.mexw64`) to convergence and reproduce the committed goldens within
tolerance.

Decisions made during brainstorming (with the user):

1. **Test gate = run against the old MEX**, not snapshot-only. This forces the
   `GDSGE_DATA` layout and the MEX caller-workspace contract to be honored now,
   making Phase 6 a formality.
2. **Keep all public behavior.** Everything reachable via documented options
   (WarmUp, MaxIter, IterSaveAll, adaptive bounds, the resolve cascade) is
   preserved behaviorally; only modernized in form.
3. **Thin generated file + runtime library.** Model-agnostic logic moves into
   unit-testable `gdsge.runtime.*` functions; the generated file keeps only
   model-specific code.
4. **Pure emitters, no .m templates.** Section-by-section emitter functions
   build the file text; no placeholder substitution for the MATLAB output.

## 2. Architecture

```
IR (validated struct)
  └─ gdsge.codegen.generateMatlab(ir, outDir)
       ├─ gdsge.codegen.mat.emitIter(ir)      ──▶ iter_<model>.m
       └─ gdsge.codegen.mat.emitSimulate(ir)  ──▶ simulate_<model>.m

generated files at runtime
  ├─ call gdsge.runtime.* (hand-written, unit-tested library)
  └─ call flat kernels in src/kernels/ (vendored as-is: myppual, …)
```

### 2.1 Components

```
src/+gdsge/+codegen/
  generateMatlab.m       % driver: IR -> two file strings, writes them to outDir
  +mat/                  % MATLAB-emitting backend ("mat" avoids shadowing matlab.*)
    emitIter.m           %   assembles iter file from the sections below
    emitSimulate.m       %   assembles simulate file
    emitSetup.m          %   defaults + params/grids/shocks/shock_trans literals
    emitBounds.m         %   GDSGE_LB/UB blocks + adaptive-bound update fragments
    emitInterpInitial.m  %   interp-var initialization (state refs -> GDSGE_TENSOR_*)
    emitDataPack.m       %   the GDSGE_DATA packing expression (layout contract)
    emitSolUnpack.m      %   c1=GDSGE_SOL(1:1,:); ... from IR slot layout
    emitResultIter.m     %   explicit IterRslt assembly (no v2struct)
    emitResultSimu.m     %   SimuRslt preallocation + per-period assignments
    rewriteNames.m       %   word-boundary identifier rewriting in opaque exprs
    codeWriter.m         %   tiny indented-line buffer
src/+gdsge/+runtime/
  unpackOptions.m        % whitelisted GDSGE_OPTIONS unpacking; error on unknown
  ensurePath.m           % essential_blas.dll PATH setup (replaces header dance)
  solveProblems.m        % mex-call wrapper + full resolve cascade
  applyWarmUp.m          % WarmUp SOL/LB/UB interpolation transfer
  constructSplines.m     % GDSGE_PP construction via myppual + SPLINE_VEC assembly
  computeMetric.m        % max-abs diff across interp vars
  printIterProgress.m    % structured progress line, PrintFreq/NoPrint gated
  reportUnconverged.m    % failed points: count, (shock,state) coords, residuals
src/kernels/             % vendored flat MATLAB/MEX runtime kernels, moved as-is:
                         %   myppual.m, myppual_mex.<mexext>, essential_blas.dll,
                         %   convert_to_interp_eval_array.m, gen_discrete_markov_rn.m,
                         %   get_scalar.m (+ dependency closure found during impl.)
```

The kernels keep flat names because generated code and `IterRslt.pp` consumers
reference them by those names. The "new source" for path purposes is now
`src/` **plus** `src/kernels/`; the harness addpaths both (still exactly one
toolbox source per process — the golden rule is about never mixing old and new).

### 2.2 The MEX caller-workspace contract (made explicit)

The old MEX reads, via `mexGetVariable("caller", ...)`, from the workspace of
the function that invokes it:

- Scalars: `MEX_TASK_NAME`, `MEX_TASK_INIT`, `MEX_TASK_INF_HORIZON`, `TolFun`,
  `TolSol`, `SolMaxIter`, `NumThreads`, `GDSGE_DEBUG_EVAL_ONLY`, `UseBroyden`,
  `FiniteDiffDelta`, `GDSGE_USE_BROYDEN_NOW`
- Struct: `GDSGE_SPLINE_VEC` (spline path; `GET_STRUCT` in
  `interp_spline_construct_template.cpp`)

With the thin-runtime design, **`gdsge.runtime.solveProblems` is the direct
caller**: it receives a config struct, defines these names as locals in its own
workspace, then invokes the model MEX through a function handle. This contract
is documented here and pinned by a test (fake mex handle that asserts the
variables exist via `evalin('caller', ...)`).

Positional MEX interface (unchanged):
`[SOL,F,AUX,EQVAL,OPT_INFO] = mex_<model>(SOL,LB,UB,DATA,SKIP,F,AUX,EQVAL)`.

## 3. The generated `iter_<model>.m`

Skeleton (each block = emitted model-specific code or a runtime call):

1. **Header**: function signature `[IterRslt,IterFlag] = iter_<model>(GDSGE_OPTIONS)`;
   `gdsge.runtime.ensurePath()`.
2. **Defaults block**: the frozen public option surface (`TolEq`, `TolSol`,
   `TolFun`, `PrintFreq`, `NoPrint`, `SaveFreq`, `NoSave`, `MaxIter`,
   `MaxMinorIter`, `SolMaxIter`, `UseBroyden`, `FiniteDiffDelta`,
   `UseAdaptiveBound`, `UseAdaptiveBoundInSol`, `IterSaveAll`, `SkipModelInit`,
   `NumThreads`, …) with identical names and default values as the old
   generated file. MEX task constants (`MEX_TASK_INIT = 0`,
   `MEX_TASK_INF_HORIZON = 1`) included.
3. **Setup from the IR**: params, shock values, and the normalized
   `shock_trans` are emitted as numeric literals with 17-significant-digit
   round-trip formatting (`mat2str(v,17)`). State grids are carried in the IR
   as **opaque MATLAB text** (the two-expression-worlds rule:
   `states.grids` is `fMap(fText())` in the schema) and are emitted verbatim
   (`w1 = linspace(-0.05,1.05,201);`) *after* the param literals, so grid text
   may reference params. The original gmod setup text is otherwise not
   replayed.
4. **Options unpacking**: `gdsge.runtime.unpackOptions(GDSGE_OPTIONS, whitelist)`
   replaces `v2struct(GDSGE_OPTIONS)`. Whitelist = runtime options + model
   params + state grids + simulate controls (preserves the old param-override
   workflow, e.g. `options.beta = 0.96`). Unknown fields → clear error listing
   valid names. Declaration asserts (shock lengths, `shock_trans` shape) kept.
5. **Tensor construct**: `ndgrid` over `(1:shock_num, states…)` producing
   `GDSGE_TENSOR_*`; `GDSGE_NPROB`, `GDSGE_SIZE`, `GDSGE_SIZE_STATE`.
6. **Bounds**: `GDSGE_LB`/`GDSGE_UB` filled per IR slot ranges from the opaque
   bound expressions (emitted verbatim; state/shock refs rewritten to
   `GDSGE_TENSOR_*` by `rewriteNames` where the old generator did so).
7. **Solution-space prep**: `GDSGE_EQVAL/F/SOL/X0/AUX/SKIP/DATA` allocations,
   `rand`-based X0 (behavior parity, RNG unseeded as before).
8. **Interp initialization** (skipped under WarmUp, as before): interp vars
   initialized from opaque initial expressions with `GDSGE_TENSOR_*`
   rewriting; spline construction via `gdsge.runtime.constructSplines`.
9. **WarmUp path**: explicit assigns of `WarmUp.var_interp` fields (names known
   from IR — no v2struct), Iter restore, `GDSGE_SOL` size check, and the
   SOL/LB/UB interpolation transfer via `gdsge.runtime.applyWarmUp` (with the
   adaptive-bound slot indices passed in). Behaviorally identical to old.
10. **Main loop**:
    - `GDSGE_DATA(:) = [repmat([shock_num; params…; shock_trans(:); shocks…],1,NPROB); shockIdx(:)'; state tensors…]`
      — **exact old layout** (emitDataPack).
    - `gdsge.runtime.solveProblems(...)` — full cascade: initial solve,
      nearest-neighbor sweeps along each dimension (both directions), randomize
      restarts, `UseAdaptiveBoundInSol` hook, `MaxMinorIter` cap, Broyden flag
      handling. Returns SOL/F/AUX/EQVAL + diagnostics.
    - slot unpack (emitSolUnpack), aux unpack, `reshape(…, shock_num, [])`.
    - `GDSGE_NEW_<interp>` assigns; metric via `gdsge.runtime.computeMetric`
      (try/catch → NaN preserved); interp update; adaptive-bound widening
      (`min`/`max` against `factor*SOL`, old semantics).
    - spline reconstruct; `stopFlag`; progress print
      (`gdsge.runtime.printIterProgress`); periodic/final save with reshape to
      `GDSGE_SIZE` (`[len GDSGE_SIZE]` for array vars), output construction
      (`output_interp`, `output_var_index`), and explicit `IterRslt` assembly.
11. **Result struct** (identical field set): `shock_num`, `shock_trans`,
    `params`, `var_shock`, `var_policy`, `var_aux`, `var_tensor`, `pp`,
    `Metric0`, `Metric`, `Iter`, `var_state`, `var_interp`, `GDSGE_PROB`
    (`GDSGE_LB/UB/SOL/F/SIZE`), `var_others`, `NeedResolved`, plus
    `output_interp`/`output_var_index` under `CONSTRUCT_OUTPUT`. All built by
    explicit dot-assignment. `get_scalar` still applied to policy/aux structs.
    Save names `IterRslt_<model>_<iter>.mat`; `IterSaveAll=1` keeps
    save-everything semantics (the variable set inside the new file differs
    from old — acceptable: its purpose is a debug dump; `WarmUp` consumes only
    `IterRslt`).

HL1996 has no `model_init`, so the init-problem block is empty this phase; the
emitter takes the IR's modelInit presence as a real `if` for Phase 7.

## 4. The generated `simulate_<model>.m`

`SIMU_RESOLVE` variant only (HL1996's configuration; `SIMU_INTERP` is
parse-time-selected by gmod options no Phase-4 model uses — Phase 7):

1. Header, `ensurePath`, defaults, setup literals (as in iter), then overrides
   from `IterRslt` (`shock_trans`, params, `var_shock`, `var_state` — explicit
   assigns, no v2struct) and `unpackOptions(GDSGE_OPTIONS, …)`.
2. Solution-interp construction from `IterRslt.GDSGE_PROB.GDSGE_SOL` for
   initial guesses (both `ispc`/`ismac` branches preserved as-is for parity).
3. `SimuRslt` preallocation from IR `simulate.varSimu` + `shock`; initial
   conditions from IR `simulate.initial`; `GDSGE_OPTIONS.init` overwrite; shock
   path via `gen_discrete_markov_rn`.
4. Per-period loop: map current states from `SimuRslt`, rebuild `GDSGE_DATA`
   (same layout, states from sample vector), interp warm start, resolve via
   `gdsge.runtime.solveProblems` (same cascade semantics as the old inline
   randomize-retry loop), slot unpack, assignment of simulated vars, state
   transition through `GDSGE_SHOCK_VAR_LINEAR_INDEX` (old semantics for
   shock-indexed policies like `w1n`), periodic print/save.

## 5. Error handling & reporting

- **`unpackOptions`**: unknown option field → error
  `gdsge:options:unknownField` listing valid names. No more silent workspace
  pollution.
- **`solveProblems`**: returns diagnostics (`NeedResolved`, worst residual,
  minor-iter count). Default `MaxMinorIter = inf` keeps old retry-forever
  behavior, but retry-progress printing (unconverged points remaining per
  cascade round) is config-gated: **on** in `iter_*` (so a stuck solve is never
  a silent spin), **off** in `simulate_*` (the old simulate cascade is silent,
  and per-period retry lines over 10,000 periods would be noise). A finite
  `MaxMinorIter` exhausting triggers `reportUnconverged` + `warning`, then
  continues (old data flow) — in both files.
- **`reportUnconverged`**: maps linear problem indices to (shock index, state
  coordinates/values), prints count and the worst offenders with residuals.
- **`printIterProgress`**: single line, extended:
  `Iter:%d, Metric:%g, maxF:%g, unconverged:%d, elapsed:%gs` — gated by
  `PrintFreq`/`NoPrint` as before.
- Metric `try/catch → NaN` preserved.

## 6. IR amendment (small)

`options.numThreads` gains a **dynamic sentinel**: `0` means "resolve at
runtime via `feature('numcores')`". Parser (`resolveOptions`) emits `0` when
the gmod does not set `NumThreads`; codegen emits `feature('numcores')` for the
sentinel and a numeric literal otherwise. Schema doc + validator + HL1996
fixture/golden regen accordingly. (Without this, the IR pins the parse
machine's core count — wrong on other machines and a behavior change.)

## 7. Testing (TDD order)

1. **Runtime unit tests (no MEX)**: cascade semantics via a fake in-MATLAB
   "mex" function handle (nearest-neighbor sweep order, skip flags, randomize
   restarts, termination, caller-workspace variables present);
   `unpackOptions` whitelist/error; `computeMetric`; `applyWarmUp`;
   `constructSplines` against `myppual` outputs.
2. **Emitter unit tests**: each `emit*` vs expected text for HL1996 IR
   sections; numeric-fidelity test (`eval` of `emitSetup` output `isequal`s
   the IR values — proves the `mat2str(v,17)` round-trip).
3. **Snapshot test**: full generated `iter_HL1996.m` / `simulate_HL1996.m`
   equal committed golden text; `checkcode`/`mtree` report no syntax errors.
4. **Functional phase gate**: `capture_HL1996.m` extended to also keep the old
   compiled `mex_HL1996.mexw64` under `tests/HeatonLucas1996/oldmex/`
   (untracked build artifact; the gate test **fails with "run capture_HL1996
   first"** if absent — no silent skip). The test runs in an isolated
   process/work dir with `src/` + `src/kernels/` + the old MEX: `iter_HL1996`
   to convergence, then `rng(0823); simulate_HL1996`; both compared to the
   committed goldens via `compareNumericClose` (policy/aux/interp arrays,
   `Metric < TolEq`, `Iter` within a window — same tolerances as
   `tGoldenHL1996`).
5. **IR amendment tests**: `resolveOptions` sentinel, validator, gendoc
   no-drift, fixture regen.

## 8. Out of scope (explicit)

- C++ MEX codegen (`mex_*.cpp`, `compile_*.m`) — Phase 5.
- ASG / PCHIP / `var_tensor` / `model_init` / `SIMU_INTERP` / macro-using
  models — Phase 7 (emitters branch on the IR with real `if`s; no template
  forking).
- The flat `gdsge_codegen.m` / `gdsge.m` shims — Phase 6 (end-to-end wiring).
- Any change to numerical kernels — `myppual` & co. are moved verbatim.

## 9. Risks

- **Old-MEX coupling**: the functional gate depends on a local build artifact
  (`mex_HL1996.mexw64`). Mitigated by the capture script regenerating it and
  the gate failing loudly when missing.
- **Behavior parity of the extracted cascade**: subtle ordering differences in
  the nearest-neighbor sweeps could change iteration paths. Mitigated by
  porting the cascade line-by-line, unit tests pinning sweep order, and the
  tolerance-based golden gate (the fixed point is unique; iteration count gets
  a window).
- **Kernel dependency closure**: `myppual` may pull more files than listed
  (e.g. spline helpers). Discovered at implementation by running the
  functional gate; the design accepts adding files to `src/kernels/`.
