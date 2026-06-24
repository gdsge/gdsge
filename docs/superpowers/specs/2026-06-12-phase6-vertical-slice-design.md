# Phase 6 — Vertical Slice Green: Design Spec

- **Date:** 2026-06-12
- **Status:** Approved (design); pending implementation
- **Parent spec:** `2026-06-11-refactor-gdsge-design.md` (§15, Phase 6)

---

## 1. Goal

HL1996 runs end-to-end exactly the way the old toolbox's `test.m` drives it:

```matlab
gdsge_codegen('HL1996');
IterRslt = iter_HL1996(options);
SimuRslt = simulate_HL1996(IterRslt);
```

and matches the committed goldens. Phases 0–5 proved every pipeline stage with
hand-wired internal calls (`gdsge.parser.parseFrontEnd` → `gdsge.codegen.generate*`
→ manual compile) inside tests. Phase 6 adds the **public wiring** — the
backward-compat entry point, the IR-on-disk artifact, and the compile cache — and
proves the architecture through the public surface.

## 2. Scope decisions (settled with the user)

1. **Entry points:** `gdsge_codegen` only. The `gdsge.m` orchestrator (codegen →
   iterate → simulate with `IterRslt` caching) is deferred to Phase 7 — HL1996's
   `test.m` does not use it.
2. **IR on disk:** `gdsge_codegen` writes `<model>.gdsge.json` on every run. The IR
   is the pipeline's contract; materializing it makes runs debuggable and future
   non-MATLAB backends possible.
3. **Gate consolidation:** the new public-API gate **replaces**
   `tFunctionalCxxHL1996.m` (it strictly subsumes it: same artifacts and golden
   assertions, plus shim/cache/JSON wiring). `tFunctionalHL1996.m` (new MATLAB +
   old MEX) stays as a bisection aid; the MEX-equivalence test stays.
4. **Driver structure (Approach A):** thin flat shim + namespaced driver; backends
   consume the **in-memory** IR; the gate test asserts the written JSON decodes
   `isequalIR`-equal to it. Contract verification without putting decode on the
   hot path. (Rejected: JSON-as-data-path — encode/decode every run for no user
   value now; fat flat function — breaks the "internal modules are namespaced,
   shims are thin" rule.)

## 3. Components

Two files of product code:

### 3.1 `src/gdsge_codegen.m` — flat backward-compat shim

~3 lines: forwards `(modelName, options)` to `gdsge.codegen.codegen` and returns
its output. Sits next to `+gdsge` so the single `addpath(src)` exposes both (the
test harness already does this).

### 3.2 `gdsge.codegen.codegen(modelName, options)` — the driver

All file I/O happens in the **caller's cwd** (old behavior: generated files land
where the user runs). Steps:

1. Resolve `<modelName>.gmod` in cwd; clear error if absent.
2. `ir = gdsge.parser.parseFrontEnd(fileread(gmodFile), modelName)`.
3. Write `<modelName>.gdsge.json` via the existing `gdsge.ir.encode`.
4. `gdsge.codegen.generateMatlab(ir, pwd)` and `gdsge.codegen.generateCxx(ir, pwd)`
   — writes `iter_<model>.m`, `simulate_<model>.m`, `mex_<model>.cpp`,
   `compile_<model>.m`. Existing Phase 4/5 code, unchanged.
5. **Compile cache**, replicating old semantics:
   - Compare the freshly generated cpp text against `mex_<model>.cache` (if any).
   - Different (or no cache) → `gdsge.runtime.ensurePath()`, run
     `compile_<model>`, then write the cache **only after a successful compile**.
   - Identical → print the skip message; no recompile.
6. Return the IR struct as the single output (plays the old `model` output's
   role; see §6 for the deferred legacy outputs).

**Progress output** stays in the spirit of the old driver but structured:
`Parsing gmod file: ok`, `Compile mex file:` / `C++ source unchanged, skip
compiling.` — consistent with the refactor's informative-output goal.

## 4. Error handling

- Missing gmod → `error('gdsge:codegen:gmodNotFound', ...)` naming the file and
  the cwd searched.
- The options struct is whitelist-validated (`unpackOptions` style): unknown
  fields error informatively. `GenCodeSegment` specifically errors with
  "deferred to Phase 7" rather than being silently ignored.
- A failed MEX compile propagates the error and does **not** write the cache, so
  the next run retries (matches old behavior).
- Parser/codegen errors pass through untouched — they already carry located
  messages.

## 5. Testing

### 5.1 New gate — `tests/HeatonLucas1996/codegen/tEndToEndHL1996.m` (tag `Slow`)

Replaces `tFunctionalCxxHL1996.m`. In a temp work dir with `HL1996.gmod` copied
in and cwd set there (harness path policy unchanged: `src/` + `src/kernels`
only):

1. `ir = gdsge_codegen('HL1996')` — through the **flat shim**, as `test.m` does.
2. Assert all artifacts exist: `iter_HL1996.m`, `simulate_HL1996.m`,
   `mex_HL1996.cpp`, `compile_HL1996.m`, `HL1996.gdsge.json`,
   `mex_HL1996.cache`, `mex_HL1996.<mexext>`.
3. Assert `gdsge.ir.decode(fileread('HL1996.gdsge.json'))` is
   `gdsge.ir.isequalIR` to the returned IR.
4. Call `gdsge_codegen('HL1996')` again; assert the skip path (MEX binary
   mtime unchanged).
5. Iterate + simulate with the same golden assertions as the Phase 5 gate:
   `Metric < 1e-6`, `Iter ∈ (100, 400)`, `shock_trans`/`var_state` tight
   (1e-12), `var_policy`/`var_aux`/`var_interp` within 1e-4, bit-exact shock
   path under `rng(0823)`, early simulated paths (first 100 periods) within
   1e-2, long-run mean/std within 5e-3.

### 5.2 Fast unit tests for the driver

Cheap logic, no compile required: missing-gmod error, unknown-option error,
`GenCodeSegment` deferral error, cache skip/invalidate decision (structure the
driver so the cache comparison is testable without invoking `mex`).

### 5.3 Done means

Full suite green via `matlab -batch "cd('tests'); run_tests"` (junit.xml + exit
code authoritative), including the new gate. `PROGRESS.md` Phase 6 checked off
with a changelog entry.

## 6. Explicitly deferred

- `gdsge.m` orchestrator shim → Phase 7.
- `GenCodeSegment` code-segment outputs → Phase 7.
- The old 5-output `gdsge_codegen` signature (`iterCode`, `cppCache`, `cppCode`,
  `codeSegment`) — no test model consumes the extra outputs; revisit in Phase 7
  if a model does.
- Everything already deferred by the parent spec (macros, ASG, `var_tensor`,
  `model_init`, pre/post hooks → Phase 7; SymPy backend → Phase 8).
