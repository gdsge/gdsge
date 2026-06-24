# Phase 7e Design — `var_tensor` → Informative Deferred Error

- **Date:** 2026-06-13
- **Status:** Approved (design); pending implementation
- **Parent spec:** `2026-06-11-refactor-gdsge-design.md` (§13 backward-compat surface,
  §15 Phase 7); the 7c split (`2026-06-13-phase7c-legacy-orchestration-design.md` §1)
  scheduled `var_tensor` as sub-phase 7e
- **Owner:** Wenlan Luo

---

## 1. Context

Phases 7a/7b brought all six corpus models green; 7c/7d cleared the golden-less legacy
surface (orchestration, ASG-interp simulate, macro engine). `var_tensor` was the last item
the 7c split scheduled — originally as a full new variable kind (decl → IR → MATLAB+C++
codegen → runtime, spline-path only, ⊥ ASG).

**Re-scope (owner-approved, 2026-06-13):** rather than build the full feature now, Phase 7e
makes `var_tensor` a **clean, identified parse-time error** and **defers** the real
implementation to a future phase. This mirrors the architecture-honest deferral Phase 7c
applied to `GenCodeSegment` (`gdsge:codegen:genCodeSegmentUnsupported`): the old toolbox
supported the feature, but the refactored pipeline consciously raises a clear error instead
of shipping a half-built path.

**What `var_tensor` is** (old toolbox, `base_package/gdsge/source/gdsge_parser.m`): a variable
precomputed in MATLAB at every (shock × state-grid) point via `ndgrid` broadcasting — a pure
function of states/shocks, independent of the policy unknowns (e.g. cash-on-hand
`y = z·f(k) + (1-δ)k`) — then fed per-problem into the C++ solve as a popped `double`
(`POPN(<tensor>)`; `num_data += num_tensor`). It is an *input* (unlike `var_aux`, an output),
computed once rather than recomputed inside every per-point solve. The old toolbox already
forbids `USE_ASG && num_tensor>0` (`gdsge_parser.m:1009-1010`).

**Why an error now is safe and correct.** No shipped gmod declares a non-commented
`var_tensor` — only CaoKS2016 carries a commented example
(`% var_tensor r w budget1 budget2;`), which `preprocess.m` strips at the first `%` before
declaration parsing. So nothing regresses; the corpus stays green.

**Current behavior is a silent failure.** `parseVarDecls.m:36-37` collects `var_tensor`
names into `decls.tensorNames`; `parseDeclarations.m:55` assembles them into the IR
`variables.tensor` list — but no backend consumes that list
(`dataLayout`/`emitDataPack`/`emitPop` only handle states), so a real `var_tensor` produces
**broken generated code with no diagnostic**. Phase 7e replaces that silent failure with a
fail-fast error.

## 2. Scope decisions

**In scope:**

1. A parse-time **gate**: a non-empty `var_tensor` declaration raises
   `gdsge:parser:varTensorUnsupported` with an informative message.
2. A codegen-side **internal invariant** assert (defense-in-depth) so a hand-authored IR
   JSON carrying a tensor cannot silently slip past into ignored codegen.
3. Tests + `PROGRESS.md` bookkeeping (re-scope note + a deferred entry for full support).

**Deferred to a future phase (recorded in `PROGRESS.md`):** the full `var_tensor` feature —
IR enrichment (names → names + lengths + opaque-text assignment expressions), the MATLAB
tensor construct/assign/pack, the C++ `POPN` sequence, per-period recompute in `SIMU_RESOLVE`
simulate, the frozen `IterRslt.var_tensor` result shape, and the `⊥ ASG` guard. Spec §7
sketches this so the future phase has a starting point.

**Out of scope:** SymPy backend (Phase 8); performance work beyond parity (Phase 9); any
change to the `var_tensor` `.gmod` surface syntax.

## 3. Goal & success criteria

Phase 7e is done when:

1. `parseFrontEnd` on a gmod with a non-commented `var_tensor` declaration raises
   `gdsge:parser:varTensorUnsupported`, and the message names the offending variable(s) and
   points to the inline-into-equations workaround.
2. A gmod with only a **commented** `var_tensor` (or none) parses cleanly — the gate does not
   over-fire.
3. The IR-consuming generators (`generateMatlab`/`generateCxx`) assert `ir.variables.tensor`
   is empty (internal invariant), raising an identified error if violated.
4. All existing HL1996 / 7a / 7b / 7c / 7d gates stay green. Zero new old-toolbox golden
   capture.

## 4. Approach

Established TDD discipline, smallest surface first:

1. **Failing test → gate.** Add a synthetic fixture gmod declaring `var_tensor`; assert the
   error identifier + that the message lists the names. Then add the guard in
   `parseDeclarations.m`.
2. **No-over-fire test.** Assert a commented-only `var_tensor` parses (reusing/representing
   the CaoKS2016 case) so the gate's narrowness is intentional and regression-guarded.
3. **Internal invariant.** Add the one-line assert in the codegen driver + a direct unit test
   feeding a synthetic IR with a non-empty `variables.tensor`.

Own branch `phase7e-var-tensor-error`, merged when green.

## 5. The change (detail)

### 5.1 Parser gate (primary, user-facing)

In `src/+gdsge/+parser/parseDeclarations.m`, immediately after `parseVarDecls` collects the
declarations and **before** `evalSetup` (i.e. right after line 7, before the IR `variables`
struct is assembled), raise when the tensor list is non-empty. Gating pre-eval rejects
`var_tensor` purely structurally — no dependence on the eval'd workspace — so the specific
diagnostic fires ahead of any generic setup error and a tensor-only fixture needs no valid
surrounding setup:

```matlab
if ~isempty(decls.tensorNames)
    error('gdsge:parser:varTensorUnsupported', ...
        ['var_tensor is not yet supported by the refactored GDSGE pipeline ' ...
         '(deferred to a future phase). Found tensor variable(s): %s. ' ...
         'Remove the var_tensor declaration and inline the precomputed ' ...
         'state/shock quantity into the model equations.'], ...
        strjoin(decls.tensorNames, ', '));
end
```

`parseVarDecls.m` stays a **pure collector** — it continues to gather `var_tensor` names into
`decls.tensorNames` (no feature-policy decision embedded in the tokenizing stage). Placing the
gate in `parseDeclarations` follows the per-stage-owner precedent set by
`resolveOptions.m`'s `gdsge:parser:simuModeConflict` check: the stage that owns the assembled
data owns its semantic validation.

After the gate, `decls.tensorNames` is provably empty wherever the IR is built, so the IR
`variables.tensor` list is always `{}`.

### 5.2 Codegen internal invariant (defense-in-depth)

`codegen.m` always re-parses the gmod itself (`parseFrontEnd`, line 22), so the §5.1 gate
already fires there. The invariant therefore belongs in the **IR-consuming generators** —
`generateMatlab` and `generateCxx` — which are the entry points reachable directly with a
hand-authored or hand-edited IR (bypassing the parser). Factor the guard into a small shared
helper (e.g. `gdsge.codegen.assertSupportedIR(ir)`) called at the top of both, to avoid
duplication:

```matlab
function assertSupportedIR(ir)
if ~isempty(ir.variables.tensor)
    error('gdsge:codegen:varTensorUnsupported', ...
        ['IR carries var_tensor variable(s) (%s), which the refactored ' ...
         'pipeline does not yet generate code for.'], ...
        strjoin(ir.variables.tensor, ', '));
end
end
```

This closes the path where an IR struct/JSON authored or edited by hand could otherwise reach
a generator and be silently ignored (`dataLayout`/`emitDataPack`/`emitPop` only handle
states) — consistent with the parent spec's no-silent-fallthrough principle (§11).

### 5.3 IR schema — unchanged

`variables.tensor` keeps its current spec `fList(fText())` (a name list). It remains the
forward-compat placeholder for the deferred feature; the §5.1 gate guarantees it is always
empty in any IR the parser emits. **No schema change, no doc-gen drift, no golden
regeneration.**

## 6. Testing

Harness unchanged: headless `matlab.unittest` via `matlab -batch "cd('tests'); run_tests"`;
`junit.xml` + exit code authoritative; one-source-per-process path policy.

- **`tVarTensorUnsupported` (new suite):**
  - *Positive:* a minimal synthetic fixture gmod with `var_tensor` (plus the minimum valid
    surrounding declarations to reach `parseDeclarations`) → `parseFrontEnd` throws
    `gdsge:parser:varTensorUnsupported`; assert the identifier and that the message contains
    the declared names.
  - *No over-fire:* a gmod identical except the `var_tensor` line is commented (`%`) parses
    without error (mirrors the live CaoKS2016 case), so the gate's narrowness is explicit.
  - *Codegen invariant:* build a small valid IR (e.g. clone the HL1996 reference IR), inject a
    non-empty `variables.tensor`, and assert `generateMatlab` and `generateCxx` (via the
    shared `assertSupportedIR` helper) throw `gdsge:codegen:varTensorUnsupported`.
- **Regression:** all existing gates (HL1996, safe_assets, Mendoza2010, GLSW, CaoKS2016,
  Bianchi2011, macro/orchestration suites) stay green — no corpus model declares a real
  tensor, and the commented CaoKS2016 example is stripped before parsing.

No old-toolbox golden capture (consistent with the 7c/7d golden-minimization precedent — the
error path produces no numerical result to compare).

## 7. Deferred: full `var_tensor` (future-phase sketch)

Recorded here so the future phase starts from the analysis already done, not from scratch.
Full support, spline-path only (⊥ ASG), would add:

- **IR:** enrich `variables.tensor` from a name list to items carrying `name`, `length`
  (scalar per the old `POPN`-one-per-tensor convention), and the **opaque-MATLAB-text
  assignment expression** (the RHS, a pure function of states/shocks — same opaque-text world
  as grids/bounds/setup, *not* the model-body AST).
- **dataLayout / emitDataPack / emitPop:** append per-problem tensor rows after the state
  rows; `<tensor>(:)'` in the MATLAB packer, `POPN(<tensor>)` in C++; `nData += num_tensor`.
- **emitIter (MATLAB):** after the existing `GDSGE_TENSOR_<state>`/`shockIdx` ndgrid
  construction, emit the tensor assignments with state/shock names rewritten to
  `GDSGE_TENSOR_<name>` (the broadcast forms), producing `[shock_num × state-grid]` arrays.
- **Simulate (`SIMU_RESOLVE`):** recompute tensors per period from the simulated
  states/shocks (the old `simuTensorCode = regexprep(tensorAssignCode,'GDSGE_TENSOR_','')`
  transform — no ndgrid, scalar/vector ops) and pack their rows. (`SIMU_INTERP` needs no
  tensor handling — it evaluates the output interpolant without re-solving.)
- **Result:** assemble `IterRslt.var_tensor` (frozen shape, old `v2struct(RSLT_TENSOR)`).
- **Guard:** keep an explicit `USE_ASG + var_tensor` rejection (subsumed today by the blanket
  error).
- **Fixture:** a purpose-built minimal spline model where the tensor is load-bearing (its
  value enters a residual), differential-tested against a freshly captured old-toolbox golden.

## 8. Risks

- **Over-fire on commented examples.** Mitigation: `preprocess.m` strips `%` comments before
  `parseVarDecls`; the no-over-fire test pins this. Low risk.
- **Bypassing the parser.** A hand-authored IR JSON could carry a tensor. Mitigation: the §5.2
  codegen invariant catches it. Low risk.
- **Future-phase drift.** The deferred sketch (§7) could go stale against the codebase.
  Mitigation: it is explicitly a starting point, re-validated when the future phase opens.

## 9. Result-struct & API contract

Unchanged. No model in the corpus produces an `IterRslt.var_tensor`; the field is introduced
only when full support lands (§7). Public entry points and all existing result shapes are
untouched by this phase.
