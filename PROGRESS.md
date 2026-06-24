# PROGRESS — GDSGE toolbox refactor

Living tracker for the refactor. Update this file as work completes. The authoritative
design is `docs/superpowers/specs/2026-06-11-refactor-gdsge-design.md`.

**Summary:** Rewrite the GDSGE generation layer around an explicit JSON **IR**: a modular
MATLAB parser emits the IR, and interchangeable backends generate the runtime MATLAB and the
C++ MEX solver (default = adept autodiff, no Python; alternative = SymPy stack-allocated
analytic Jacobian). Drop `v2struct`, add error checking and informative output, stay
backward compatible. Build vertical-slice first on **HeatonLucas1996**, autodiff backend.

**Status legend:** ☐ not started · ◐ in progress · ☑ done

---

## Phases

- ☑ **Phase 0 — Infrastructure** (done 2026-06-11)
  - ☑ git init, `.gitignore`, design spec committed
  - ☑ Package skeleton (`src/+gdsge/{+parser,+ir,+codegen,+runtime}`, `templates/`)
  - ☑ `CLAUDE.md`
  - ☑ `PROGRESS.md` (this file)
  - ☑ `docs/notes-for-agents.md`, `docs/ir-schema.md` stub
  - ☑ MATLAB test harness (`tests/run_tests.m`, `tests/run.ps1`, smoke test) — verified pass & fail paths
  - ☑ Golden comparison utility + unit tests (`gdsgetest.compareNumericClose`, 6 tests)
  - ☑ uv Python env + SymPy smoke test (Python 3.12.11, sympy 1.14.0)
  - ☑ Old-toolbox build verification + HL1996 golden capture (**risk gate passed**: MSVC 2022, Iter=209, Metric=9.58e-07)
  - ☑ HL1996 golden-integrity test
- ☑ **Phase 1 — IR schema + MATLAB scaffolding** (done 2026-06-12)
  - ☑ Declarative schema descriptor (`gdsge.ir.schema`) — single source of truth
  - ☑ AST node API (`gdsge.ir.node.*`) + accessors
  - ☑ Descriptor-driven validator (shape, refs, slots, option invariants)
  - ☑ Type-aware round-trip (`canonicalize`/`encode`/`decode`/`isequalIR`); the 4 quirks defeated
  - ☑ Doc generator (`gdsge.ir.gendoc` -> `docs/ir-schema.md`, no-drift test)
  - ☑ HL1996 reference IR + JSON golden (Phase 3's forward target)
- ☑ **Phase 2 — Parser front-end** (macros-stub + lexer + block-split + declaration parser → partial IR) (done 2026-06-12)
  - ☑ `params` IR section added; HL1996 options corrected to real toolbox defaults
  - ☑ `gdsge.parser`: preprocess, splitBlocks, splitStatements, eval sandbox
    (`defaultSetupCode`/`evalSetup`), parseVarDecls, parseDeclarations (+ slot layout),
    parseSimulate, resolveOptions, assemblePartialIR, parseFrontEnd
  - ☑ HL1996 partial-IR golden (section-by-section vs the reference IR) green
  - ☑ Full macro engine (`#define/#for/...`) deferred to Phase 7 (HL1996 uses none)
- ☑ **Phase 3 — Model-expression parser** (done 2026-06-12)
  - ☑ `tokenize` (positions, primes, reduction braces) + `parseExpr` (full MATLAB
    precedence: `^` left-assoc, unary between power and multiplicative; all 7 node kinds)
  - ☑ `parseModel`: assign / primed assign / `GDSGE_INTERP_VEC'` / all 4 reductions
    (+ `| trans` pipe syntax); inline reductions hoisted as `GDSGE_<KIND>_<n>`
  - ☑ `analyzeModel`: name resolution, square system, interp arity, transRef
  - ☑ Phase gate green: `parseFrontEnd(HL1996.gmod)` `isequalIR` the reference IR
  - ☑ Fixture temp renamed `__ep_expect` → `GDSGE_EXPECT_1` (+ JSON golden regen)
- ☑ **Phase 4 — MATLAB codegen** (done 2026-06-12)
  - ☑ Vendored flat kernels (`src/kernels/`: myppual, essential_blas, …)
  - ☑ `gdsge.runtime`: unpackOptions (whitelist + error), ensurePath, solveProblems
    (resolve cascade, owns the MEX caller-workspace contract incl. `GDSGE_PP_*`),
    applyWarmUp, constructSplines, computeMetric, printIterProgress, reportUnconverged
  - ☑ `gdsge.codegen.mat` pure emitters + `gdsge.codegen.generateMatlab` driver
    (no v2struct; explicit IterRslt/SimuRslt assembly; frozen result shapes)
  - ☑ IR amendment: `numThreads` 0 = dynamic sentinel → `feature('numcores')`
  - ☑ Snapshot goldens + **functional gate green**: new `iter_HL1996.m` +
    OLD `mex_HL1996` converge at Iter=209, Metric=9.58e-07 (matches golden
    capture exactly); simulate reproduces the golden shock path bit-exactly
- ☑ **Phase 5 — C++ codegen (autodiff)** (done 2026-06-12)
  - ☑ Vendored `include/` (adept headers) + `essential_blas.lib`; MSVC build proven
  - ☑ `templates/cxx/` cleaned templates + `fillTemplate`/`readTemplate` engine
  - ☑ `emitExpr`: direct AST→C++ printer, no Symbolic Toolbox dependency
  - ☑ Shared `gdsge.codegen.dataLayout` contract: MATLAB packer + C++ `POP*` macros from one descriptor
  - ☑ Section emitters (`emitTask`, `emitCompile`, `emitArgDeclare`, `emitPop`, `emitInterp`, `emitModelBody`, …) + `generateCxx` driver + snapshot goldens
  - ☑ MEX equivalence gate (old-vs-new on identical inputs, eval-only + solve): residuals match to 1e-8
  - ☑ Functional gate green: NEW generated MATLAB + NEW compiled MEX → Iter=209, Metric=9.57709e-07; simulate reproduces golden shock path bit-exactly and all field moments within 5e-3
- ☑ **Phase 6 — Vertical slice green** (done 2026-06-12)
  - ☑ `src/gdsge_codegen.m` flat shim + `gdsge.codegen.codegen` driver (parse →
    `<model>.gdsge.json` → generateMatlab/Cxx into cwd → cache-gated MEX compile)
  - ☑ Old cache semantics replicated (`mex_<model>.cache`; cache written only after a
    successful compile); options whitelist errors; `GenCodeSegment` deferred error
  - ☑ Gate green: `tEndToEndHL1996` drives the public API like the old `test.m`,
    asserts artifacts + IR JSON round-trip + second-run compile skip + goldens;
    replaces `tFunctionalCxxHL1996`
  - ☑ Deferred to Phase 7: `gdsge.m` orchestrator, `GenCodeSegment`, legacy 5-output
    signature
- ☑ **Phase 7a — Spline-path widening** (done 2026-06-13)
  - ☑ Goldens captured (old toolbox): safe_assets (Iter=1271/Metric 9.92e-7),
    mendoza2010 coarse (Iter=117) → fine-grid 80×80 WarmUp re-solve (Iter=161/Metric
    9.57e-7), GLSW_interp (Iter=903/Metric 9.71e-9); reduced seeded simulates
    6×1000 (all three) + GLSW init-override batch 6×100
  - ☑ IR: `setup` section (gmod declaration code replayed verbatim — old setParamsCode
    parity; fixes Barro `Re_n` and GLSW `Ngrid`/`a1_lb`/`a1_ub` undefined in bounds);
    transition `primed` flag (fixes unprimed-transition mis-indexing in Mendoza/GLSW);
    init-bound item with text name
  - ☑ Parser: model_init block + `var_policy_init`/`var_aux_init`/`inbound_init`
    declarations; `var_simu` → `var_output` → `var_aux` promotion (old parity, prints
    same informational lines); `SimuPrintFreq`/`SimuSaveFreq` flow into IR options;
    named scalar interp calls (surface HL1996 never exercised); `adaptive()` factors
    naming setup variables
  - ☑ Codegen: `task_init` via `initView` emitter reuse; `-DHAS_INIT` + `max()` MAXDIM
    in compile step; iter init segment (`SkipModelInit`-gated, runs under WarmUp —
    old `iter_template.m` parity); `SIMU_INTERP` simulate variant (`emitSimulateInterp`:
    `GDSGE_PP = IterRslt.output_interp` + per-period `myppual_mex` eval; console-print
    deliberate simplification — prints Periods+elapsed, not all SimuRslt fields — SimuRslt
    contents unaffected); `emitExpr` future-var indexed calls
  - ☑ Gates green: `tEndToEndSafeAssets`, `tEndToEndMendoza2010`, `tEndToEndGLSW`
    (all passed first run after front-loaded emitter work) + all HL1996 gates stay green;
    all on branch `phase7a-spline-widening`
- ☑ **Phase 7b — ASG widening** (done 2026-06-13)
  - ☑ Goldens captured (old toolbox): CaoKS2016 (Iter=281/Metric 9.92e-5, iter 14.4s);
    Bianchi2011_asg full two-stage driver (stage-1 MaxIter=50/Metric 2.35e-3, stage-2
    Iter=69/Metric 9.12e-7); seeded simulates 6×1000; all goldens ≤0.6 MB
  - ☑ Public `asg` class + `asg_mex` vendored as-is (`src/asg.m`, `src/kernels/`);
    `gdsge.codegen.ensureAsgMex` cache-gated auto-compile (improvement over the old
    ship-prebuilt + manual compile_asg.m); committed mexw64 (myppual precedent)
  - ☑ Front-end/IR: Asg option set (AsgMinLevel/AsgOutputMaxLevel/AsgOutputThreshold);
    ASG+SIMU_INTERP parser error (7c); `setupNames` section (gmod setup-assigned names →
    runtime override whitelist — Bianchi's bMin/bMax); square shock_num×shock_num
    parameters as reduction-pipe targets (CaoKS `| shock_trans2`); shared
    `optionsWhitelist` extracted from emitIter/emitSimulate
  - ☑ C++ backend: interp_asg_* templates vendored; class-handle interop
    (GDSGE_ASG_HANDLE → AsgInterpArrayAdoubleEvaluator); -DUSE_ASG + ASG_MAX_* defines
    from asg.get_mex_constants(); init task emits empty interp section (old parity)
  - ☑ Runtime: `solveProblemsAsg` (old cascade: solve → warm-from-NEW → restore-from-OLD
    → randomize; init = solve+randomize); solveProblems carries GDSGE_ASG_HANDLE
  - ☑ MATLAB codegen: `emitIterAsg` (level-by-level refinement loop, SkipModelInit-gated
    ASG init segment, output construct, frozen ASG IterRslt field set incl. MetricVec/
    asg_*_structs), `emitSimulateAsg` (per-period sol-interp warm start + re-solve);
    ASG defaults block gated on interpMethod (spline output stays byte-identical)
  - ☑ Gates green: `tEndToEndCaoKS2016` (Iter=281, Metric matches golden exactly) and
    `tEndToEndBianchi2011` (two-stage driver: stage-1 Iter=50 bit-match, stage-2 Iter=69/
    Metric 9.12e-7 bit-match; shock overrides/WarmUp/SkipModelInit/MaxIter staging) —
    Bianchi passed first run; IR JSON snapshots committed for both models
- ☑ **Phase 7c — Legacy orchestration & ASG-interp simulate** (`gdsge.m` orchestrator,
  legacy 5-output `gdsge_codegen` signature, `GenCodeSegment` → informative error,
  ASG+`SIMU_INTERP` simulate variant) — design:
  `docs/superpowers/specs/2026-06-13-phase7c-legacy-orchestration-design.md`. The smaller
  legacy items: first three reuse existing corpus goldens; the fourth needs one small new
  capture (CaoKS2016 + `SIMU_INTERP`). The original Phase 7c grab-bag was split into 7c/7d/7e
  on 2026-06-13 (owner-approved), following the 7a/7b focused-sub-phase precedent — no corpus
  model exercises any of these items; synthetic/differential fixtures.
- ☑ **Phase 7d — Macro engine** (done 2026-06-13) — `#define/#for/#foreach/#if/#mat{}/#strcat_comma/include/cinclude` as a single `gdsge.parser.expandMacros` stage; parity + two cleanups (multi-token #define, nested #for); `gdsge:parser:macro*` error taxonomy; cinclude via optional IR `hooks.cxxIncludes` -> `GDSGE_OTHER_INCLUDE`; verified by expandMacros unit tests + cinclude C++-injection gate + IR-equivalence gate (macro-ized HL1996 isequalIR plain HL1996) — zero new golden capture.
- ☑ **Phase 7e — `var_tensor` deferred error** (done 2026-06-13) — re-scoped (owner-approved):
  instead of building the full variable kind, a non-commented `var_tensor` now raises a clear
  parse-time error (`gdsge:parser:varTensorUnsupported`, gate in `parseDeclarations` before
  `evalSetup`), with a defense-in-depth codegen invariant
  (`gdsge.codegen.assertSupportedIR` → `gdsge:codegen:varTensorUnsupported`, called by both
  generators). Mirrors Phase 7c's `GenCodeSegment` deferral; replaces today's silent
  broken-codegen path. No corpus gmod declares a non-commented tensor (CaoKS2016's is a
  stripped comment), so nothing regresses; zero golden capture. IR schema unchanged
  (`variables.tensor` stays the always-empty placeholder). Verified by `tVarTensorUnsupported`
  (parser gate + no-over-fire + both-generator invariant). Design:
  `docs/superpowers/specs/2026-06-13-phase7e-var-tensor-error-design.md`.
- ◐ **Phase 7e-full — `var_tensor`** — **MATLAB-side support landed in Phase 9a** (ndgrid tensors
  feeding inbound bounds + initial interp; IR `variables.tensor` enriched to `{name,expr}`).
  Still deferred: **C++-body tensors** (the `POPN` path, needed only if a tensor enters the model
  residual — none in the corpus; `analyzeModel` rejects it with `varTensorInBodyUnsupported`) and
  the `IterRslt.var_tensor` result field (no corpus golden carries it). ⊥ ASG. 7e spec §7 sketch.
- ☑ **Phase 8 — SymPy analytic-Jacobian backend** (done 2026-06-13) — Approach C:
  MATLAB owns the fused shock-loop structure + a gradient registry (sparse slot→C++
  rows, forward-mode chain rule, per-shock "templated" slots); SymPy (pyenv bridge,
  JSON in/out) differentiates each body and returns value+partials as shared-CSE C++.
  Chain rule closes through the spline `search_eval_with_grad_vec_at_array` (double).
  Selected by `UseAutoDiff=0` → `options.jacobianBackend='sympy'` (emitted only when
  non-default; zero golden regen). Reductions EXPECT/MIN/MAX/PROD; single-state +
  multi-state `GDSGE_INTERP_VEC'`; named scalar interp; primed assigns/equations;
  constant-indexed per-shock access (`Re_n(2)`). Gated three-way (sympy=adept=FD
  Jacobian, HL1996) + end-to-end golden on ALL FOUR spline models
  (HL1996/safe_assets/Mendoza/GLSW); synthetic MIN/MAX/PROD FD checks.
  cxx/asg/pchip under sympy → clear errors. ASG + pchip → Phase 8b (deferred).
- ◐ **Phase 8b — SymPy backend for ASG (done 2026-06-15); pchip still deferred** — the ASG interp
  derivative for the chain rule already exists in pure double (`asg_adouble.h::eval_vec_adept`); the
  SymPy path reuses it (new `eval_vec_with_grad`) so CaoKS2016 / Bianchi2011_asg now run under
  `UseAutoDiff=0`, matching their autodiff goldens. **pchip** remains deferred (no C++ backend at all).
- ☑ **Phase 9a — Cao2011EZ corpus model** (done 2026-06-14) — Epstein–Zin two-agent asset
  pricing, green end-to-end on **both** backends (autodiff + SymPy) vs a captured golden.
  Four new features: multiple inline reductions per statement (EZ Euler eqs carry two `EXPECT`
  each — generalized `parseModel`'s single-reduction hoist); embedded named-interp call hoisting
  (`qp' = qFuture'(Wp')+d'` → `GDSGE_INTERP_n'`); `qp_GRID(n)` per-shock indexing (parser `_GRID`
  normalization — both backends already handle `qp(const)`); MATLAB-side `var_tensor` (Phase
  7e-full subset). Both backends converge at Iter 275 / Metric 9.56e-7 matching the golden.
- ☑ **Phase 9b — CaoNie2016 corpus model** (done 2026-06-14) — two-agent production economy
  with a collateral constraint, green end-to-end on **both** backends (autodiff + SymPy) vs a
  regenerated old-toolbox golden (Iter 1077 / Metric 9.90e-7, bit-matched by both backends).
  **Two new control-structure constructs** (probe-confirmed — every expression form already
  parsed): conditional model regions `model(X>0)`/`model(X==0)` and `if/else/end` conditional
  equations on the shock value; plus `var_others` wiring. **IR evolved to 1.1.0**: `model.regions`
  list + tagged `plain|conditional` equations; a `regionView` adapter + `modelBody` helper keep the
  per-body emitters region-agnostic so only the drivers loop regions; all 8 model IR snapshots
  regenerated once (migration differential proves shape-only). Shared `lowerCondition` lowers the
  guard to C++; both backends emit guarded region bodies + guarded equation slots (SymPy = per-branch
  value+Jacobian, no new symbolic math). **Latent bug found & fixed:** the model templates wrote
  `ARG_CODE;`/`DECLARE_CODE;`, dropping a stray `;` onto the trailing `#define name(idx)` macro —
  harmless until CaoNie2016 used `Xp(1)`/`qp(1)` *mid-expression* (`consis=.../Xp(1)-1.0` silently
  lost the `-1.0`, so nothing converged). Full suite 447/447 green.
- ☑ **Phase 10 — Polish** (done 2026-06-14) — docs (README + architecture + user guides), IR
  version freeze (policy + changelog + `tIrVersionFreeze`), perf check vs old (harness + report),
  cleanup (Mendoza golden 99.8MB→9.5MB, deferred-error map). Full suite 450/450 green. The only
  remaining open items are the deferred sub-phases **8b** (SymPy ASG+pchip) and **7e-full
  remainder** (C++-body `var_tensor`), plus a logged perf follow-up (safe_assets iter overhead).

Each phase gets its own spec→plan→implement cycle (see `docs/superpowers/plans/`).

---

## Decisions log

- **Test harness:** native `matlab.unittest`, run headless via `matlab -batch` (no Python in
  the correctness loop). MATLAB is R2025b.
- **First vertical slice:** **HeatonLucas1996** (Huggett1997 was removed). Autodiff (adept)
  backend first as the backward-compat anchor; SymPy backend follows in Phase 8.
- **IR format:** versioned JSON (`jsonencode`/`jsondecode`).
- **MATLAB path policy:** never persistent; one source per process; see `CLAUDE.md`.
- **Backward compat:** public flat entry points preserved as thin shims; only internal
  modules namespaced under `+gdsge`.
- **IR schema:** one declarative descriptor (`gdsge.ir.schema`) drives validation,
  round-trip, and doc-gen. Round-trip quirks (orientation, array-kind, empty/scalar,
  non-finite) handled type-aware in `canonicalize`/`toJsonReady`. Two expression worlds:
  model body → AST; grids/bounds/initial/transitions → opaque MATLAB text. Reductions are
  flat statements (nested ones hoisted by the parser). Validator checks IR-computable
  properties; model semantics (square system, bounds-completeness) stay parser-side.

## Risks / open

- **Old toolbox must build here** to capture goldens (needs a MATLAB-configured C++
  compiler). Phase 0 verifies this first; if it fails, pause for a user decision.
- **SymPy reduction-to-loop fusion** (spec §10) is the most novel piece; deferred to Phase 8
  and cross-checked against autodiff goldens.
- **Reproducibility** is tolerance-based, not bit-exact.

## Changelog

- 2026-06-16: **Backend auto-selection + informative compile output.** The C++ Jacobian
  backend is now chosen at codegen time by `gdsge.codegen.resolveBackend` with precedence
  in-gmod `UseAutoDiff` > `GDSGE_BACKEND` env (`adept|sympy|auto`) > auto-detect (Python
  present ⇒ SymPy, else adept). `UseAutoDiff` became tri-state in `resolveOptions` (unset ⇒
  no IR field/auto; =1 ⇒ `'autodiff'`; =0 ⇒ `'sympy'`), so default-path IR stays
  byte-identical (no corpus gmod sets it). In **auto mode only**, a SymPy *codegen* failure
  falls back to adept (`generateCxxWithFallback`); explicit/env pins error instead.
  `generateCxx` gained an `opts.backend` override; the kernel-ensure branch was extracted to
  `gdsge.codegen.ensureKernels`. The driver prints a one-line `Backend: …` message every run
  and numbered `[i/5]` phase banners (parsing → MATLAB → kernels → C++ → solver MEX), and
  `compile_<model>.m` now labels the simulate-MEX step. The test harness pins
  `GDSGE_BACKEND=adept` to keep the correctness loop Python-free; sympy gates opt in via
  `UseAutoDiff=0`. Two pre-existing tests updated to the new semantics
  (`tJacobianBackendOption` autodiff-on-flag-on; `tEnsureSplineConstructMex` call-chain via
  `ensureKernels`). Full suite **538/538** green. Spec/plan:
  `docs/superpowers/{specs,plans}/2026-06-16-backend-auto-selection*`. Branch
  `feature/backend-auto-selection`.
- 2026-06-16: **In-MEX cartesian SIMU_INTERP simulation + myppual retired.** Cartesian
  `SIMU_INTERP=1` now runs the whole period loop in a generated `simulate_<model>_mex`
  (hans-style: search-once per site, OpenMP over samples, in-place writes; chunked by
  `SimuSaveFreq`). `IterRslt.output_interp` and the resolve warm-start interp moved to a
  **stacked uniform-order** layout (shock = stacked vector index, cubic over states),
  built by `interp_construct_mex` and evaluated by a new generic `interp_eval_mex` (both
  on the double-only `include/interp_eval_double.h`). `myppual.m`, `myppual_mex`,
  `convert_to_interp_eval_array.m`, `get_scalar.m`, and the `useFusedConstruct` toggle are
  **deleted**; non-expressible simulate blocks fall back to a per-period `interp_eval_mex`
  loop. Whole-loop MEX is **2.2–2.5× faster** than the fallback (see `docs/perf-report.md`).
  GLSW (unprimed) + Mendoza (primed) e2e green; full suite 495 tests, 0 failures. Branch
  `feat/simu-interp-mex`. Spec/plan: `docs/superpowers/{specs,plans}/2026-06-15-cartesian-simu-interp-mex*`.
- 2026-06-15: **SymPy analytic-Jacobian backend for ASG (Phase 8b, ASG slice).** Lifted
  `gdsge:codegen:sympyInterpUnsupported` for `interpMethod='asg'` (in **both** the IR validator and
  `generateCxx` — the end-to-end path hits the validator first), so ASG models run on the SymPy backend
  (`UseAutoDiff=0`) with the same support as the spline path. **Key fact:** the ASG value+gradient is
  already computed in pure double (`include/asg_adouble.h::eval_vec_adept`; `eval_vec_adouble` just scales
  `/= stateRange` and injects into adept), so the chain rule needed only a double entry point — new
  `eval_vec_with_grad` (mirrors `eval_vec_adouble` minus the adept injection) — plus an ASG branch in
  `emitInterpSympy` emitting the shared `GDSGE_INTERP_VEC_double_grad` lambda with the native `[v*DIM+d]`
  → spline `[v + NVEC*d]` gradient transpose. Everything downstream (the `+sympymodel` registry/chain-rule,
  `emitTask` region loop, `emitCompile` ASG defines) was already in place. **Latent SymPy bug fixed en
  route** (surfaced by CaoKS2016, the first ASG model on SymPy): the plain-assign emitter in
  `sympymodel/generate.m` always wrote `double <name>`, so a reassigned scalar
  (`rhs1 = lambda1; ... rhs1 = rhs1 + EXPECT{...}`) redeclared `rhs1` → C2374; now it tracks declared
  targets and reassigns bare. Additive (no spline SymPy model reassigns a scalar). **Verified:** CaoKS2016
  + Bianchi2011_asg end-to-end on SymPy match their autodiff goldens (Iter 281 / two-stage), and a
  **sympy↔adept↔FD Jacobian cross-check** (`tSympyJacobianCaoKS2016`) agrees per equation/unknown at one
  interior off-node point (sympy==adept 1e-6, sympy==FD 1e-4; the point is kept off ASG grid nodes because
  the piecewise-linear interpolant's gradient jumps there). New: `tEmitInterpSympyAsg`,
  `tGenerateCxxAsgSympy`, `tEndToEndCaoKS2016Sympy`, `tEndToEndBianchi2011Sympy`, `tSympyJacobianCaoKS2016`
  (+ `gdsgetest.buildAsgMexInputs`). **Descoped:** the planned **CaoNie2016_asg** capstone (ASG × Phase-9b
  conditional regions) was dropped — that model does not converge cleanly under the reference toolbox, so it
  is not a reliable golden; the ASG × conditional-region combination is consequently unexercised by the
  corpus (no other model uses it). **pchip** under the C++ backend remains a separate deferral.
  Spec `docs/superpowers/specs/2026-06-15-phase8b-sympy-asg-design.md`, plan
  `docs/superpowers/plans/2026-06-15-phase8b-sympy-asg.md`. Branch `phase8b-sympy-asg`.

- 2026-06-15: **In-MEX randomized restart (`UseMexRandomize`, default on).** The randomized-restart
  trial loop — and the outer-iter-1 initial guess — moved out of MATLAB into the C++ MEX for the
  **cartesian iter path** (both backends; shared `task.tpl.cpp`/`mex.tpl.cpp`). A **counter-based
  per-grid-point RNG** (`gdsge_splitmix64` in `mex.tpl.cpp`) makes every draw a pure function of
  `(MexRandomSeed, salt=GDSGE_Iter, point, trial, component)` — **thread-independent** (identical at
  any `NumThreads`), **batch-resumable**, and reproducible from the seed alone. With the feature on,
  the iter initial guess is the deterministic **bound midpoint** (no MATLAB `rand`), so MATLAB's `rng`
  no longer affects iter results (it still drives `simulate` shock draws). The MEX runs up to
  `MexRandomizeBatch` (default **100**) full minor-iterations internally, returning to MATLAB only at
  batch boundaries for diagnostics, the `MaxMinorIter` cap, and the converge/exit decision. New public
  options: `UseMexRandomize=1`, `MexRandomizeBatch=100`, `MexRandomSeed=0` (whitelisted; `emitSetup`
  defaults). `solveProblems` gained the batched-call driver branch (`useMexRandomize`); `=0` keeps the
  legacy MATLAB restart path byte-for-byte; `UseAdaptiveBoundInSol=1` forces that fallback. **Scope:**
  iter cartesian (spline/pchip) only — **ASG and `simulate`/`init` keep the MATLAB restart path**.
  **Critical perf fix found during benchmarking:** the first cut ran the full-grid in-MEX
  `neighbor_fixpoint` on **every** trial, launching thousands of no-op OpenMP grid sweeps when a
  restart made no progress — making safe_assets ~**3× SLOWER** than legacy despite *identical solve
  counts and fewer MEX calls* (measured: 3605 vs 3582 solves, 4 vs 6 calls). Gating the sweep on a
  change in the unconverged count (mirrors `solveProblems.m`'s `numNeedResolvedAfter` gate; skipping a
  no-op sweep is result-preserving → bit-exact) flipped it to ~**3× FASTER**: safe_assets **2.5–3.0s
  vs legacy 5.6–9.7s**, with far less seed-variance (`tests/perf/inmex_randomize_bench.m`). **Golden
  re-baseline (owner-approved):** safe_assets (multi-root) converges to the **same equilibrium** under
  the new default (Iter=1271/Metric 9.92e-7; value functions ~1e-8, all policies ~1e-6, 9/10 aux ~1e-8)
  **except auxiliary `x1` at 2/1002 grid points**, where the seed-0 RNG selects a different valid root
  (residual <1e-8); `tests/Barro_et_al_2017/golden/IterRslt.mat` re-captured from the new pipeline
  (reproducer `tests/golden/rebaseline_safe_assets.m`); SimuRslt golden unchanged. **Verification:**
  new `tInMexRandomizeDeterminism` (bit-identical across `NumThreads`, seed-driven), `tInMexRandomizeAB_HL1996`
  + `tInMexRandomizeSafeAssets` (both modes reach the equilibrium within the project's 1e-4 iter
  tolerance), `tInMexRandomizeDriver` (fast fake-MEX driver), plus codegen/snapshot tests. Adjusted:
  `tFunctionalHL1996` pins `UseMexRandomize=0` (it drives the OLD reference MEX, which has no in-MEX
  randomize loop); `tInMexResolve{HL1996,SafeAssets}` pin `UseMexRandomize=0` in both arms (clean
  neighbor-sweep A/B). HL1996 iter+cpp snapshots regenerated. Full suite **491/491** green. Spec
  `docs/superpowers/specs/2026-06-15-in-mex-randomize-design.md`, plan
  `docs/superpowers/plans/2026-06-15-in-mex-randomize.md`. Branch `inmex-randomize`.

- 2026-06-15: **Bounded retries + tiered convergence diagnostics (cartesian path).** The inner
  resolve loop in `gdsge.runtime.solveProblems` no longer retries random restarts forever:
  `MaxMinorIter` default `inf → 200`, plus a new `DiagnoseMinorIter` (default 10) **heartbeat
  interval** and a new pure formatter `gdsge.runtime.diagnoseConvergence` that reports, per named
  policy component, NaN count; bounds **pinned** among unconverged points (hints bounds too tight);
  **near-bound** among converged points; the **dominant equation residual** vs the other equations
  (hints a bug/normalization need); and the worst-K grid points (iter only). While stuck, a summary
  prints every `DiagnoseMinorIter` retries; at the cap the **iter** path prints the full block and
  throws `gdsge:runtime:solveNotConverged` (`cfg.errorOnNonconvergence=true`), while **simulate**
  prints the same block but warns and continues (now bounded). All new behavior is gated on
  `isfield(cfg,'diagnoseAt')`, so converging models are byte-identical (functional goldens unchanged;
  only the two HL1996 text snapshots regenerated). New helper `gdsge.codegen.mat.solComponentNames`
  maps SOL rows to named components (e.g. `w1n(3)`), threaded into emitIter/emitSimulate `cfg`.
  **Cap revised 20 → 200 mid-implementation (owner decision):** the corpus surfaced that
  **safe_assets** (multi-root grid points) and the **Mendoza2010 sympy** gate legitimately need a
  heavy-tailed, RNG-dependent number of restarts — the old default was `inf`, and the
  golden-reproducing seed (`rng('default')`) needs **>1000** restarts on safe_assets' first outer
  iteration. So no fixed hard-error cap is safe for them; the 5 restart-heavy gates
  (`tEndToEndSafeAssets`, `tEndToEndSafeAssetsSympy`, `tFusedConstructSafeAssets`,
  `tInMexResolveSafeAssets`, `tEndToEndMendoza2010Sympy`) pass `MaxMinorIter=inf` +
  `DiagnoseMinorIter=inf` (and `rng('default')` where unpinned), restoring the unbounded regime their
  goldens were captured under. Real users of such models raise `MaxMinorIter` themselves — the
  diagnostics tell them to. The diagnostics also revealed/relied on: `NeedResolved` stays an
  `IterRslt` field (only the dead post-return warning block was dropped from `emitIter`). Full suite
  **482/482** green. ASG path + outer-loop `MaxIter` stalls remain out of scope. Spec
  `docs/superpowers/specs/2026-06-14-convergence-diagnostics-design.md`, plan
  `docs/superpowers/plans/2026-06-14-convergence-diagnostics.md`. Branch `convergence-diagnostics`.

- 2026-06-14: **Build `GDSGE_CFG` once before the iter loop (code tidy; perf negligible).** Following
  the `GDSGE_DATA` hoist, `emitIter` now builds the solver-config struct once before the `while` and
  updates only the three per-iteration fields inside (`splineVec`, `ppCell`, `useBroydenNow`); the other
  ~14 fields are set from options and are loop-invariant (`solveProblems` treats `cfg` read-only, so it
  persists). **Bit-exact** (full suite **471/471** green; only the HL1996 iter snapshot regenerated).
  **Honest perf note:** this was expected to help the cheap-solve/high-iteration marshal but the measured
  win is **negligible** (~0.004 s on GLSW, ~0.008 s on HL1996 — within noise): MATLAB builds the 17-field
  struct in ~4 µs, so it was never the bottleneck. Kept as a clean-code change. The remaining iter-loop
  marshal (~0.12 s / ~0.13 ms per iter on GLSW) is interpreter overhead spread across ~40 small
  per-iteration statements (`GDSGE_SOL`/`GDSGE_AUX` row-slices, reshapes, interp update) with no single
  dominant line; chasing it further is diminishing returns on an already-fast model, so the marshal hunt
  stops here. The two real wins remain the **`GDSGE_DATA` hoist** (~18% on HL1996) and the **fused spline
  construction** (1D rebuild ~halved). Branch `iter-cfg-hoist`.

- 2026-06-14: **Hoist the loop-invariant `GDSGE_DATA` pack out of the iter loop (perf).** A line-level
  re-profile of the iter loop (`scratch/profile_fused.m`/`bench_data.m`) showed the per-iteration
  `GDSGE_DATA(:) = [repmat([params],1,NPROB); GDSGE_TENSOR_grid]` rebuild is the **dominant marshal cost
  on models with a large param block × NPROB** (HL1996: **1.24 ms/iter × 209 ≈ 0.26 s**, ~83% of its
  marshal bucket), while negligible on small ones (GLSW 0.011 s). But `GDSGE_DATA` is **loop-invariant**:
  the param block + the fixed solve grid (`GDSGE_TENSOR_*`) never change across iterations, and the data
  that *does* change (the policy/value functions) reaches the MEX via `GDSGE_SPLINE_VEC`/`GDSGE_PP_*`, not
  `GDSGE_DATA`. `emitDataPack` already builds the param block once for the **simulate** path
  (`GDSGE_data0`); `emitIter` now does the same for **iter** — emit `pack.iterPack` **once before the
  `while`** instead of inside it, gated on no `pre_iter` hook (`hoistData`; no corpus gmod has one, and a
  `pre_iter` that could mutate a param keeps the in-loop rebuild). **Bit-exact** (identical constant data,
  built once → identical MEX input → identical `IterRslt`); only the HL1996 iter **text snapshot** was
  regenerated (`regen_snapshots.m` — the `GDSGE_DATA(:)` line moved above the loop), no functional golden
  changed. Shared by **both** backends (the MATLAB iter file is backend-agnostic); ASG is untouched (its
  `emitIterAsg` grid genuinely changes per level). **Verified:** full suite **471/471** green. **Perf:**
  HL1996 iter-loop marshal bucket **~25% → ~8.6%** (≈0.31 s → 0.087 s), HL1996 iter wall ~1.20 s → ~0.99 s
  (~18%); GLSW barely moves (its `GDSGE_DATA` was already tiny — its remaining marshal is the
  output-stack cat/reshapes, a separate lever). Branch `iter-data-hoist`.

- 2026-06-14: **Fused in-MEX spline construction (post-Phase-10 perf enhancement).** Replaced the
  per-outer-iteration cartesian spline/linear rebuild — per-`var_interp` `myppual` dispatch +
  `convert_to_interp_eval_array` (a coefficient **concat copy + `permute`** into eval order) — with
  **one** dedicated construct-MEX call (`src/kernels/interp_construct_mex.cpp`). It takes all interp
  `Values` via a **struct** (MATLAB shares the buffers, no concat), runs the **identical**
  construction math (the `mkl_start`/`mkl_start_with_extrap` driver copied verbatim from the owner's
  `myppual_mex.cpp`, over `mkl_dummy_interp.h`'s `construct_cubic_spline_notaknot` /
  `construct_linear_interp` — *not* real MKL), and writes the eval-order `GDSGE_SPLINE_VEC.coefs`
  **directly** (no permute) plus the per-variable natural-order pp structs — in a single call.
  **Bit-exact by construction:** the permute is a pure reordering of the same doubles, so the result
  is byte-identical to the legacy path. Scope = the cartesian **spline/linear** path only (the
  construction-only kernel is minimal/self-contained; the prebuilt **`myppual_mex` is retained,
  unchanged, for the simulate / `output_interp` eval path**); ASG, pchip, and the C++ solver MEX
  are untouched (same `GDSGE_SPLINE_VEC` + per-`GDSGE_PP_<name>` contracts → no solver/template
  change). `gdsge.runtime.constructSplines` makes the single call (signature unchanged → generated
  `iter_<model>.m` unchanged); the legacy `myppual`+`convert` path stays reachable behind
  `gdsge.runtime.useFusedConstruct(false)` for the A/B gates. New `ensureSplineConstructMex`
  cache-gated compile (mirrors `ensureAsgMex`; flags from the owner's `compile_myppual.m` + `/O2`,
  MSVC `/fp:precise` to preserve bit-exactness), invoked at `generateCxx` entry for spline models;
  committed prebuilt `.mexw64`. **Owner refinement:** `extrap_order==4` (cubic) attaches **no** extra
  knots — a cubic spline extrapolates as a cubic naturally — so it reuses the plain cubic
  construction (the old `myppual.m` errors on it); covered by a dedicated `extrap=4 == plain cubic`
  test. **Verified:** full suite **471/471** green (prior 454 + 17 new). New tests: `tInterpConstructMex`
  (driver — 1D/2D/3D, cubic/linear, extrap=2, extrap=4, frozen pp field set: `splineVec.coefs` and
  per-var pp **bit-identical** to `myppual`/`convert_to_interp_eval_array`); `tConstructSplinesEquivalence`
  (fused == legacy); A/B end-to-end `tFusedConstruct{HL1996,SafeAssets}` (`iter` bit-identical to the
  legacy constructor — `var_policy`/`var_aux`/`var_interp` + `Iter` + `Metric` + `IterRslt.pp` frozen
  shape — rng-pinned); `tEnsureSplineConstructMex`. **Perf** (`tests/perf/fused_construct_bench.m`,
  rng-pinned, identical Iter/Metric): **GLSW ~6%** faster (0.62s vs 0.66s), **safe_assets ~9%** faster
  (3.07s vs 3.37s) — the savings is the eliminated permute/concat/dispatch marshalling (the underlying
  `mkl_start` compute is unchanged). Spec
  `docs/superpowers/specs/2026-06-14-fused-spline-construct-design.md`, plan
  `docs/superpowers/plans/2026-06-14-fused-spline-construct.md`. All on branch `fused-spline-construct`.
  **Follow-up — cache-friendly eval-order pack.** An iter-composition re-profile (`scratch/profile_fused.m`,
  fused vs legacy on the three Phase-10 spline cases) showed the spline-rebuild bucket roughly **halves**
  on 1D models (GLSW 30.6%→16.2%, HL1996 10.5%→4.9%) but initially **regressed ~4% on the 2D model
  (mendoza)**: the original eval-order writer did `numVec` full **scattered** sweeps of the coef array
  (`evalCoefs[iv+numVec*baseMap[natIdx]]`) plus a per-element modulo/divide map build, which loses to
  MATLAB's cache-blocked `permute` on the larger 2D coef array. Replaced it with a **single-pass odometer**:
  walk the natural buffer sequentially (cache-friendly per-var reads), accumulate the eval base
  incrementally (no modulo/divide), and write each `numVec` value as a contiguous burst — one scattered
  write-line per `natIdx` instead of `numVec` sweeps. mendoza's fused rebuild dropped **0.561s→0.248s**
  (2.3×), flipping it to a win; now **all three cases** are faster than legacy on both the rebuild bucket
  and wall-clock (HL1996 1.37→1.20s, GLSW 0.65→0.56s, mendoza 3.65→3.46s). Bit-exact preserved
  (`tInterpConstructMex`/`tConstructSplinesEquivalence` 12/12, A/B gates 2/2).

- 2026-06-14: **In-MEX nearby-warmup resolve (post-Phase-10 perf enhancement).** Ported the
  hans/dpopt `solver_use_resolve` idea: the nearby-neighbor warm-start sweep (resolve cascade
  step 2) now runs **inside the C++ MEX** (`templates/cxx/task.tpl.cpp`) instead of via
  MATLAB↔MEX round-trips, collapsing the whole sweep into **one** MEX call per outer iteration;
  the random restart stays in MATLAB. Default on (`UseMexResolve=1`, an undocumented runtime
  option / internal switch; the MATLAB sweep is kept behind it). Scope: the cartesian-grid
  (spline/pchip) path, covering **both** backends via the shared task template; ASG and the init
  task pass no strides and are untouched. **Design: deterministic snapshot semantics** (freeze
  the converged set per direction-pass, then fused copy-neighbor+solve) — *not* libdpopt's racy
  live-read — so it is **bit-exact** with the MATLAB sweep. `solveProblems.m` keeps its control
  flow byte-for-byte: only the sweep *body* is replaced by one strides-enabled MEX call, guarded
  by the persistent `numNeedResolvedAfter` entry-gate; the per-point C++ solve became a reusable
  lambda and `call_fmin` dropped its `sol0→sol` memcpy (warm start now pre-placed by the caller).
  **Verified:** full suite **454/454** green on both backends (all goldens match, incl.
  safe_assets multi-root and the ASG models); a new in-process **A/B differential gate**
  (`tInMexResolve{SafeAssets,HL1996}`) asserts `iter(UseMexResolve=1)` vs `(=0)` are
  **bit-identical** (`maxabsdiff=0` on every `var_policy`/`var_aux`/`var_interp` field, `Metric`
  diff 0) **once `rng` is pinned before each run** — the A/B comparison surfaced the same
  RNG-fragility the end-to-end gates pin for: the iter loop draws `rand()` for the initial guess
  and every restart, so unsynchronized runs diverge to ~1e-12 (a test artifact, not a resolve
  difference). Plus a fast driver unit test (`tInMexResolveDriver`, nested fake MEX) locking the
  call structure. **Perf:** safe_assets benchmark (`tests/perf/inmex_resolve_bench.m`, rng-pinned,
  identical Iter=1271/Metric) ~**16% faster** (3.00s vs 3.59s) — modest because solve compute
  dominates marshalling on this small model, but it partially addresses the logged Phase-10
  safe_assets follow-up (the sweep round-trips are eliminated; the multi-root random restarts
  remain in MATLAB by design). Spec `docs/superpowers/specs/2026-06-14-in-mex-nearby-resolve-design.md`,
  plan `docs/superpowers/plans/2026-06-14-in-mex-nearby-resolve.md`. All on branch
  `inmex-nearby-resolve`. Deferred/possible-future: an in-MEX scheme for ASG (interpolation-based
  warm start, no cartesian neighbors).
- 2026-06-14: **Phase 10 complete (the refactor's finishing phase).** Four polish workstreams on
  branch `phase10-polish`, full suite **450/450** green (the prior 447 + 3 new `tIrVersionFreeze`).
  **(1) Cleanup** — the Mendoza2010 golden `IterRslt.mat` trimmed **99.8 MB → 9.5 MB** by dropping the
  three uncompared bulk fields (`output_interp` ~74 MB, `pp` ~28 MB, `GDSGE_PROB` ~16 MB) and keeping
  only the 8 fields the gates read (`tests/golden/trim_mendoza_golden.m`, committed reproducer); the
  TODO scan found only one marker, inside the vendored `myppual` kernel (no action);
  `docs/deferred-features.md` consolidates every honest-error ID (grep-verified — the IDs differed
  from the spec draft: the live ones are `varTensorInBodyUnsupported`/`varTensorAsgUnsupported`, not
  the removed `varTensorUnsupported`). **(2) IR version freeze** — `docs/ir-versioning.md` (major/minor
  policy + 1.0.0→1.1.0 changelog); `decode` already enforced the major gate, so `tIrVersionFreeze`
  locks it (accepts any 1.x, rejects 2.0.0, tripwire on the schema major). No schema change → no
  golden regen. **(3) Perf check vs old** — a reusable one-source-per-process harness (`tests/perf/`:
  `perf_worker.m` + `run_perf.m`/`.ps1`, honoring the MATLAB-path golden rule, file-existence as the
  success signal since MEX teardown can crash a child *after* it writes its result) and a committed
  `docs/perf-report.md` over a representative subset (HL1996/safe_assets/GLSW/CaoKS2016 × {old,
  autodiff, sympy}; CaoKS2016 ASG = autodiff only). **3 of 4 at/above parity** (CaoKS2016 ~6% faster
  on new autodiff; SymPy faster than autodiff on the small spline models); all iter counts/metrics
  match. The harness flagged `safe_assets` at ~3× (18.7s vs 6.3s), but **follow-up profiling showed
  that is RNG-fragile, not a fixed regression**: safe_assets' multi-root restart depends on the RNG
  state at iter start, and the same build profiles at 6.8s ≈ old when the RNG is perturbed — the real
  driver is MEX *solve-call count* from the resolve cascade (≈12.5 solves/iter vs ~1 on clean models),
  not MATLAB overhead. Profiler attribution (corrected in `docs/perf-report.md`): the dominant MATLAB
  overhead is the **generated `iter_<model>.m` marshalling loop** (~13–26%, backend-agnostic) plus the
  **per-iteration spline rebuild** (~30% on cheap/high-iteration GLSW); `solveProblems` is minor (~5%)
  except where a resolve cascade makes it loop extra real solves; `computeMetric`/`printIterProgress`
  are negligible. SymPy's analytic Jacobian is the *fastest* MEX solve (so MATLAB overhead is a larger
  fraction of its runtime). Follow-ups (non-blocking): pin `rng` before iter for safe_assets and
  re-measure; if optimizing, target the iter-loop marshalling + spline rebuild, not the backend.
  **(4) Docs** —
  root `README.md`, `docs/architecture.md` (pipeline, module map, how-to-add-a-model, golden rule),
  `docs/user-guide.md` (gmod authoring, options, SymPy backend, reading results). Spec
  `docs/superpowers/specs/2026-06-14-phase10-polish-design.md`, plan
  `docs/superpowers/plans/2026-06-14-phase10-polish.md`. **The GDSGE refactor is functionally
  complete**; remaining open items are the deferred sub-phases 8b (SymPy for ASG+pchip) and 7e-full
  remainder (C++-body `var_tensor`), plus the logged safe_assets perf follow-up.
- 2026-06-14: **Phase 9b complete.** CaoNie2016 (two-agent production economy with a collateral
  constraint) added to the backward-compat corpus, green end-to-end through the public API on **both**
  C++ Jacobian backends vs a regenerated old-toolbox golden (Iter 1077 / Metric 9.90e-7, bit-matched
  by autodiff *and* SymPy). **Two new control-structure constructs**, both confirmed first by a
  model-body probe (`scratch/probe_caonie2016.m`) that narrowed the inventory — every expression form
  already parsed (array policy `Xp[3]`, `qp(1)`/`Xp(1)` const-index, multi-output primed interp, `Eqp`
  reuse all covered by HL1996/9a): (1) **conditional model regions** `model(X>0)`/`model(X==0)` —
  `splitBlocks` learns the `model(<cond>)` opener and emits `modelRegions`; `parseFrontEnd` iterates
  them; each region is an independent square system selected per grid point by its state guard; (2)
  **`if/else/end` conditional equations** on the shock value — `parseModel`'s equations parser groups
  control lines into a `conditional` entry (cases of plain equations); both backends emit a runtime
  branch filling the slot. Plus **`var_others`** → `IterRslt` (MATLAB codegen). **IR evolved to
  1.1.0** (owner-approved clean schema evolution over additive-optional): `model` becomes a `regions`
  list, equations a tagged `plain|conditional` list (new `eqkinds` registry threaded through
  schema/validate/canonicalize/encode/decode/gendoc). A `gdsge.ir.regionView` adapter + a
  `gdsge.codegen.cxx.modelBody` helper keep the per-body emitters (modelVars/emitModelBody/
  emitEquations/sympy generate) region-agnostic, so only the **drivers** (emitModel, sympy emitTask,
  analyzeModel) loop regions — the migration touched builders/gates but not the per-body emitter
  internals. All 8 existing model IR snapshots regenerated once; a migration differential
  (`tRoundtrip/existingModelsAreSingleUnconditionalRegion`) proves the change is shape-only (regions[0],
  condition '', plain tag). New shared `gdsge.codegen.cxx.lowerCondition` lowers the opaque guard to a
  C++ boolean (states/shocks are C++ locals, params `#define`d, `~=`→`!=`); both backends emit guarded
  region bodies + guarded equation slots. SymPy = per-region/per-branch guarded value+Jacobian (the
  existing `diff_body` machinery reused; **no new symbolic capability**). **Golden:** the vendored old
  `gdsge_parser` regenerates CaoNie2016 cleanly (≤1 literal `EXPECT`/line, so the 9a two-EXPECT blocker
  is absent), so the golden was captured from a fresh old-toolbox build (Iter 1077, ~10s iter, rng 0823)
  — strongest parity; reduced seeded simulate 6×1000. **Latent codegen bug found & fixed:** the model
  templates wrote `ARG_CODE;`/`DECLARE_CODE;`, dropping a stray `;` onto the *last* substituted line —
  for an array-policy / future-var section that line is a `#define name(idx) name_GDSGE_GRID[int(idx)-1]`
  macro, so its body gained a `;`. Harmless until used MID-expression: HL1996 never uses `w1n(idx)`
  mid-expression and Cao2011EZ only at expression-end, but CaoNie2016's `consis=.../Xp(1)-1.0` silently
  dropped the `-1.0`, so the solver thrashed (2103/2103 unconverged). Fixed in both autodiff + SymPy
  templates; HL1996 cxx snapshot regenerated (4 macro lines, EOL-normalized comparison). Full suite
  **447/447 green**. Spec `docs/superpowers/specs/2026-06-14-phase9b-caonie2016-design.md`, plan
  `docs/superpowers/plans/2026-06-14-phase9b-caonie2016.md`. All on branch `phase9b-caonie2016`.
  CaoNie2016 is the last corpus model → next is **Phase 10 (polish)**.
- 2026-06-14: **Phase 9a complete.** Cao2011EZ (Epstein–Zin two-agent asset pricing) added to the
  backward-compat corpus, green end-to-end through the public API on **both** C++ Jacobian backends
  (adept autodiff + SymPy) vs a freshly captured golden. **Four new features**, all confirmed first
  by a model-body probe (`scratch/probe_cao2011ez.m`): (1) **multiple inline reductions per
  statement** — the EZ Euler equations carry two `GNDSGE_EXPECT` each; generalized `parseModel`'s
  single-reduction hoist to hoist every reduction to `GDSGE_<KIND>_<n>` (whole-RHS single-reduction
  fast path kept, so existing IR snapshots stay byte-identical); (2) **embedded named-interp call
  hoisting** — `qp' = qFuture'(Wp')+d'` hoists the interp call to `GDSGE_INTERP_<n>'` (whole-RHS
  named call unchanged); (3) **`qp_GRID(n)` per-shock indexing** — a one-line tokenizer `_GRID`
  normalization (`qp_GRID(1)`→`qp(1)`); **zero backend code** — autodiff's future-var call branch
  and SymPy's `liftIndexed` already handle `qp(const)`; (4) **MATLAB-side `var_tensor`** (the
  Phase 7e-full subset) — `variables.tensor` enriched to `{name,expr}`; ndgrid tensors emitted in
  iter, feeding `inbound` bounds + `initial` interp; per-period stripped recompute in `SIMU_RESOLVE`
  simulate; tensors never enter the C++ body so both backends are unaffected (no `IterRslt.var_tensor`
  — the golden carries none). Both backends converge at **Iter 275 / Metric 9.55841e-07**, matching
  the golden (var_policy/aux/interp within 1e-4). **Unplanned fixes found en route:** `parseSimulate`
  now defaults `num_samples=2` (Cao2011EZ omits it — old-toolbox parity); two stale deferral guards
  removed (`generateCxx` var_tensor rejection; IR validator sympy+tensor rejection). **Golden
  capture note:** the vendored old `gdsge_parser` *cannot* regenerate Cao2011EZ (it rejects 2
  `EXPECT`/line — the limit we lifted), so the golden was captured by running the committed pre-built
  2016 mex. The simulate gate compares period-1 `q`/`p` only (the 2016 golden's shock-draw RNG
  predates `gen_discrete_markov_rn`) and strips the `GDSGE_EMPTY` placeholder (golden predates it).
  Full suite **431/431 green**. Spec `docs/superpowers/specs/2026-06-13-phase9a-cao2011ez-var-tensor-design.md`,
  plan `docs/superpowers/plans/2026-06-13-phase9a-cao2011ez.md`. All on branch `phase9a-cao2011ez`.
  Next: Phase 9b (CaoNie2016) or Phase 10 (polish).
- 2026-06-13: **Phase 8 complete.** SymPy analytic-Jacobian C++ backend (alternative to
  adept autodiff), selected by `UseAutoDiff=0` → `options.jacobianBackend='sympy'` (optional
  IR field, emitted only when non-default, so the default path stays byte-identical — zero
  golden regen). **Approach C:** MATLAB owns the fused shock-loop structure and a gradient
  registry (name→sparse gradRow of `{slot,expr,templated}`, forward-mode chain rule; templated
  slots are per-shock array-unknown elements `base+(iter-1)`); SymPy (in-process via `pyenv`
  at the uv venv, one JSON string in/out) differentiates each body and returns value+partials
  as shared-CSE C++ (`helper_i`). The §10 reduction-fusion chain rule closes through the
  spline kernel's `search_eval_with_grad_vec_at_array` (pure double — no kernel work).
  Handles: EXPECT/MIN/MAX/PROD reductions (fused value+grad loops; MIN/MAX extremal-j
  subgradient; PROD product rule); single- and multi-state `GDSGE_INTERP_VEC'`; named scalar
  interp calls (evaluated via the shared `GDSGE_SPLINE_VEC`); primed assigns/equations
  (per-shock value+grad arrays); constant-indexed per-shock access (`Re_n(2)` lifted to a
  scalar with a shock-specialized row). New: `gdsge.codegen.sympy.{ensurePyenv,callSympy}`,
  the `pyext/gdsge_sympy/` package (`ast_to_sympy`, `diff_body`, `ccode_helpers`),
  `gdsge.codegen.cxx.+sympymodel/*` (registry + emitters), `emitInterpSympy`,
  `model_sympy.tpl.cpp` + `call_fmin_sympy.tpl.cpp`, and an additive 6th-output analytic
  Jacobian (`GDSGE_DEBUG_EVAL_ONLY==2`, shared `task.tpl.cpp`) for cross-checking. Guards:
  sympy + cxx/asg/pchip → clear errors; Python-unavailable → `uv sync` hint; all sympy tests
  `assumeTrue(sympyAvailable)` so the default Python-free loop stays green. **Verification:**
  HL1996 three-way Jacobian cross-check (sympy == adept to 1e-6 == finite-diff to 1e-4 at the
  golden solution, all grid points); end-to-end golden match on ALL FOUR spline models —
  HL1996, safe_assets (MIN + named interp + `Re_n(2)`), Mendoza2010 (2-state interp + WarmUp),
  GLSW (named interp + model_init + SIMU_INTERP); synthetic MIN/MAX/PROD FD checks; Python
  unit tests. The cross-check caught a real bug: a reduction's collapsed row double-counted a
  slot when a templated (interp-chain) entry aliased a fixed (collapsed-MIN) entry onto the
  same `omega1n` slot (MIN-shift pattern) — fixed by registering the unique touched-slot set
  (`dacc` already accumulates both); safe_assets went from 14-min thrashing to Iter=1271 in
  ~11s. ASG + pchip under sympy deferred to Phase 8b. All on branch
  `phase8-sympy-jacobian-backend`. Next: Phase 8b (ASG) or Phase 9 (polish).
- 2026-06-13: **Phase 7e complete (re-scoped).** `var_tensor` deferred to a future phase;
  declaring one now fails fast instead of silently generating broken code. Parser gate in
  `parseDeclarations` (after `parseVarDecls`, before `evalSetup`) raises
  `gdsge:parser:varTensorUnsupported`, naming the offending variables and pointing to the
  inline-into-equations workaround — purely structural, so it fires ahead of any setup-eval
  error. Defense-in-depth: a shared `gdsge.codegen.assertSupportedIR(ir)` invariant, called at
  the top of `generateMatlab` and `generateCxx`, raises `gdsge:codegen:varTensorUnsupported`
  if a hand-authored IR carries a tensor (`dataLayout`/`emitDataPack`/`emitPop` only handle
  states). Mirrors the Phase 7c `GenCodeSegment` architecture-honest deferral. No shipped gmod
  declares a non-commented `var_tensor` (CaoKS2016's `% var_tensor ...` is stripped by
  `preprocess` before parsing), so the corpus stays green; IR schema unchanged
  (`variables.tensor` remains the always-empty forward-compat placeholder — no golden
  regeneration). Zero old-toolbox capture. New suite `tVarTensorUnsupported` (5 tests: parser
  gate, message-names-vars, commented-no-over-fire, both-generator invariant). Full suite green
  (390/390). All on branch `phase7e-var-tensor-error`. Next: Phase 8 (SymPy analytic-Jacobian
  backend), with full `var_tensor` support tracked as a deferred sub-phase (7e spec §7).
- 2026-06-13: **Phase 7d complete.** Macro engine replaces the preprocess.m
  guard: new `gdsge.parser.expandMacros` stage runs the ordered passes
  (cinclude/include/#define/#strcat_comma/#mat/#foreach/#for/#if) on raw text
  before line-cleaning. Fidelity = parity where the old engine worked, plus two
  sanctioned cleanups (multi-token #define values; nested #for via balanced
  #for/#end matching) and a `gdsge:parser:macro*` error taxonomy
  (macroDefineMalformed/macroForUnterminated/macroForeachUnterminated/
  macroIfUnterminated/macroEvalFailed/macroIncludeNotFound). #foreach keeps the
  old boundary rule (#id substituted only at non-alphanumeric boundaries, so
  K#i stays literal); #for substitutes with no leading boundary (x#i->x1) but a
  trailing boundary so nested prefix-named iterators (i/i2) don't collide;
  #mat is char-aware (scalar/string ok, array errors). cinclude flows through a
  new optional IR `hooks.cxxIncludes` field into the existing
  GDSGE_OTHER_INCLUDE placeholder (no template change; absent-when-empty so
  macro-free model IR snapshots are byte-identical — no golden regeneration).
  Verified by direct expandMacros unit tests (tExpandMacros), a cinclude
  C++-injection structural gate (tCincludeInjection, incl. an end-to-end
  parseFrontEnd->generateCxx case reusing HL1996), and an IR-equivalence
  differential (tMacroEquivIR: a macro-ized HL1996 using #define/#mat/#foreach
  parses isequalIR the plain HL1996 IR) — ZERO new golden capture, since macros
  reduce to an already-covered model. All HL1996/7a/7b/7c gates stay green; full
  suite green (385/385). Macro error taxonomy includes macroForeachUnterminated
  (parity with #for/#if unterminated handling). Two deliberate spec refinements:
  the GNDSGE/deprecate
  rewrite stays in preprocess (behavior-equivalent); hooks.cxxIncludes is
  optional/present-only-when-nonempty. All on branch phase7d-macro-engine.
  Next: Phase 7e (var_tensor).
- 2026-06-13: **Phase 7c complete.** Legacy orchestration surface + ASG-interp
  simulate. `gdsge.codegen.codegen` gained a 2nd output (generated code strings +
  pre-overwrite cache); the flat `gdsge_codegen` re-exposes the legacy 5 outputs
  `[model,iterCode,cppCache,cppCode,codeSegment]` (single-output callers
  unaffected); `codeSegment` carries the real generated strings (no old template
  fragments). New flat `gdsge.m` orchestrator: codegen -> cache-gated solve
  (IterRslt_<model>.mat + iter_<model>.cache) -> simulate -> eq{model,IterRslt,
  SmltRslt}, no v2struct; HL1996 gate proves the golden match and the second-run
  skip-solve. `options.GenCodeSegment` now raises an architecture-honest error
  (gdsge:codegen:genCodeSegmentUnsupported) — the thin-file architecture has no
  segment decomposition. ASG+SIMU_INTERP: removed the 7b parser guard; extended
  emitSimulateInterp with an ASG branch (asg.construct_from_struct(
  asg_output_struct) + per-period eval_vec, no re-solve) and routed it via
  generateMatlab; validated differentially against a fresh old-toolbox golden
  (CaoKS2016 with SIMU_INTERP flipped: own fixture/capture/integrity/front-end/
  end-to-end gates; iter Iter=281/Metric 9.92e-5, simulate shock path bit-exact,
  asg_output_struct within 1e-3). Full suite green (339/339). All on branch
  phase7c-legacy-orchestration. Next: Phase 7d (macro engine).
- 2026-06-13: **Phase 7b complete.** Both ASG models green end-to-end through the public
  API vs fresh old-toolbox goldens: CaoKS2016 (Iter=281/Metric 9.92e-5, iter 14.4s) and
  Bianchi2011_asg two-stage driver (stage-1 MaxIter=50/Metric 2.35e-3, stage-2 Iter=69/
  Metric 9.12e-7); seeded simulate shock paths bit-exact, IterRslt fields within tolerance.
  Both gates' iteration counts and metrics match their goldens bit-exactly.
  Architecture additions: the public `asg` class and `asg_mex` are vendored as-is under
  `src/`; `gdsge.codegen.ensureAsgMex` provides cache-gated auto-compile (an improvement
  over the old ship-prebuilt + manual `compile_asg.m` workflow; committed mexw64 on the
  myppual precedent). Front-end/IR widening: `setupNames` section captures gmod
  setup-assigned names as the runtime override whitelist (enabling Bianchi's `bMin`/`bMax`
  params to flow through); square shock_num×shock_num parameters accepted as
  reduction-pipe targets (CaoKS `| shock_trans2`); the Asg option set
  (AsgMinLevel/AsgOutputMaxLevel/AsgOutputThreshold) added to IR; ASG+SIMU_INTERP
  combination raises a parser error (deferred to 7c); shared `optionsWhitelist` extracted
  from emitIter/emitSimulate. C++ backend: interp_asg_* templates vendored; class-handle
  interop (GDSGE_ASG_HANDLE → AsgInterpArrayAdoubleEvaluator); -DUSE_ASG + ASG_MAX_*
  defines sourced from `asg.get_mex_constants()`; init task emits empty interp section
  (old parity). Runtime: `solveProblemsAsg` implements the old cascade (solve →
  warm-from-NEW → restore-from-OLD → randomize; init path = solve+randomize);
  `solveProblems` carries GDSGE_ASG_HANDLE through to the MEX workspace. MATLAB codegen:
  `emitIterAsg` (level-by-level refinement loop, SkipModelInit-gated ASG init segment,
  frozen ASG IterRslt field set including MetricVec/asg_*_structs) and `emitSimulateAsg`
  (per-period sol-interp warm start + re-solve); ASG defaults block emitted only for ASG
  models so spline output stays byte-identical. Deliberate deviations recorded: ASG
  console print via `printIterProgress` (7a precedent); unconditional minorIter cap fixes
  an old unbounded-loop quirk (identical at inf default); GDSGE_ASG_INTERP_OUTPUT rename.
  IR JSON snapshots committed for both models. Full suite: 330/330 tests green. All on
  branch `phase7b-asg-widening`. Next: Phase 7c (golden-less legacy surface).
- 2026-06-13: **Phase 7a complete.** Three spline models green end-to-end through the
  public API vs fresh old-toolbox goldens: safe_assets (Iter=1271/Metric 9.92e-7),
  mendoza2010 coarse→fine 80×80 (Iter=161/Metric 9.57e-7), GLSW_interp (Iter=903/Metric
  9.71e-9); reduced seeded simulates 6×1000 + GLSW init-override 6×100. Front-end
  widening: `model_init` block + `var_policy_init`/`var_aux_init`/`inbound_init`
  declarations; `var_simu`→`var_output`→`var_aux` promotion (old-parser parity);
  `SimuPrintFreq`/`SimuSaveFreq` into IR options; named scalar interp calls (surface
  HL1996 never exercised); `adaptive()` factors naming setup variables. IR additions:
  `setup` replay section (fixes Barro `Re_n` and GLSW grid intermediates dropping from
  generated files — two latent bugs found via the new models) and transition `primed`
  flag (fixes unprimed-transition RHS mis-indexing in Mendoza/GLSW). Codegen: `task_init`
  via `initView` emitter reuse; `-DHAS_INIT` + `max()` MAXDIM; `SkipModelInit`-gated iter
  init segment; `SIMU_INTERP` simulate variant; `emitExpr` future-var indexed calls.
  Console-print in the SIMU_INTERP variant deliberately simplified (prints Periods+elapsed,
  not all SimuRslt fields — SimuRslt contents unaffected). All three gates passed first run
  after the front-loaded emitter work. Full suite: 294/294 tests pass (74 suites). All on
  branch `phase7a-spline-widening`. Next: Phase 7b (ASG widening).
  Note: the Mendoza golden `tests/Mendoza2010/golden/IterRslt.mat` is ~99.8 MB — fine for
  this local-only repo, but needs Git LFS or trimming before any GitHub remote is added.
- 2026-06-12: **Phase 6 complete.** Public wiring: flat `gdsge_codegen` shim +
  `gdsge.codegen.codegen` driver (gmod from cwd → IR → `<model>.gdsge.json` via
  `gdsge.ir.encode` → existing generators → `mex_<model>.cache`-gated compile,
  cache written only after success). Options whitelist with informative errors;
  `GenCodeSegment` raises a deferred-to-Phase-7 error. New Slow gate
  `tEndToEndHL1996` proves the architecture through the public surface (artifacts,
  IR round-trip, second-run skip, golden match) and replaces `tFunctionalCxxHL1996`.
  Shared `writeText` extracted (was triplicated). All on branch
  `phase6-vertical-slice`. Next: Phase 7 (widen for backward-compat).
- 2026-06-12: **Phase 5 complete.** C++ codegen (autodiff): vendored `include/`
  (adept) + `essential_blas.lib`; `templates/cxx/` + `fillTemplate` engine;
  `emitExpr` direct AST→C++ printer (no Symbolic Toolbox); shared
  `gdsge.codegen.dataLayout` contract (MATLAB packer + C++ `POP*` from one
  descriptor); section emitters + `generateCxx` driver + snapshot goldens; MEX
  equivalence gate (old-vs-new, eval-only + solve, residuals to 1e-8); functional
  gate green: NEW MATLAB + NEW MEX → Iter=209, Metric=9.57709e-07; simulate
  reproduces golden shock path bit-exactly, all field moments within 5e-3. All on
  branch `phase5-cxx-codegen`. Next: Phase 6 (vertical slice wiring).
- 2026-06-12: **Phase 4 complete.** MATLAB codegen: thin generated files
  (`gdsge.codegen.mat` section emitters + `generateMatlab` driver) over a
  unit-tested `gdsge.runtime` library (resolve cascade with the MEX
  caller-workspace contract — 11 scalars + `GDSGE_SPLINE_VEC` + per-interp
  `GDSGE_PP_*` structs — options whitelist with errors on unknown fields,
  WarmUp transfer, structured progress + unconverged reporting). Flat kernels
  vendored to `src/kernels/`. `numThreads` 0-sentinel added to the IR.
  Functional gate: new `iter_HL1996.m`/`simulate_HL1996.m` drive the OLD
  compiled MEX to Iter=209/Metric=9.58e-07 (identical to the golden capture)
  and match the goldens (shock path bit-exact). Found & noted: vendored
  myppual mis-evaluates boundary points on 2-point grids (tests use ≥3).
  All on branch `phase4-matlab-codegen`. Next: Phase 5 (C++ codegen).
- 2026-06-12: **Phase 3 complete.** Model-expression parser: `tokenize` →
  `parseExpr` (precedence climbing, MATLAB semantics) → `parseModel` (statement
  classification, inline-reduction hoisting to `GDSGE_<KIND>_<n>`) →
  `analyzeModel` (names, square system, interp arity, transRef). `parseFrontEnd`
  now emits the complete HL1996 IR, equal to the hand-authored reference
  (`tFullIRHL1996` gate + JSON round-trip). Fixture hoisted-temp renamed to
  `GDSGE_EXPECT_1`. All on branch `phase3-expression-parser`. Next: Phase 4
  (MATLAB codegen).
- 2026-06-12: **Phase 2 complete.** `gdsge.parser` front-end turns `HL1996.gmod` into a
  schema-valid partial IR (everything but the model body), proven section-by-section against
  the corrected reference IR (golden matched on first run). Added the `params` IR section;
  corrected the HL1996 options to the real toolbox defaults (`default_mod.nmod` +
  `params_template.m`). Stages: preprocess, splitBlocks, splitStatements, eval sandbox,
  parseVarDecls, parseDeclarations (+ slot layout), parseSimulate, resolveOptions,
  assemblePartialIR, parseFrontEnd. Full `#…` macro engine deferred to Phase 7. All on branch
  `phase2-parser-frontend`. Next: Phase 3 (model-expression parser).
- 2026-06-12: **Phase 1 complete.** `gdsge.ir` package: descriptor-driven schema,
  validator, type-aware JSON round-trip, doc generator, AST node API, and a
  hand-authored HL1996 reference IR (Phase 3's golden target). All on branch
  `phase1-ir-schema`. Next: Phase 2 (parser front-end).
- 2026-06-11: Phase 0 started. Spec + Phase 0 plan written and committed; package skeleton +
  CLAUDE.md + this tracker added.
- 2026-06-11: **Phase 0 complete.** Headless MATLAB unittest harness (9 tests green), uv +
  SymPy env, tolerance comparison utility, and HL1996 goldens captured from the old toolbox
  (MSVC 2022 build verified; converged at Iter 209). All work on branch
  `phase0-infrastructure`. Next: Phase 1 (IR schema).
