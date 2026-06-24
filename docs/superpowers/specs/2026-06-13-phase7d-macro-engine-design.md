# Phase 7d Design — Macro engine

- **Date:** 2026-06-13
- **Status:** Approved (design); pending implementation
- **Parent spec:** `2026-06-11-refactor-gdsge-design.md` (§13, §15 Phase 7); sibling of the
  7c split (`2026-06-13-phase7c-legacy-orchestration-design.md`)
- **Owner:** Wenlan Luo

---

## 1. Context

Phases 7a (spline widening), 7b (ASG widening), and 7c (legacy orchestration & ASG-interp
simulate) brought all six corpus models green end-to-end through the public API and closed
most of the backward-compat surface. The original Phase 7c grab-bag was split (owner-approved,
2026-06-13) into 7c / 7d / 7e. This spec is **7d — the macro engine**.

The new front-end currently has only a **macro guard**: `src/+gdsge/+parser/preprocess.m`
raises `gdsge:parser:macroUnsupported` on any line beginning with `#`. 7d replaces that guard
with a real engine covering the full backward-compat macro surface.

**No corpus model uses any macro.** A grep of every `base_package/gdsge/tests/**/*.gmod` for
`#define / #for / #foreach / #if / #mat / #strcat / include( / cinclude` is empty. So 7d's
fixtures are **synthetic**, and (per §6) verification is by **IR-equivalence to already-covered
models** rather than fresh golden captures.

**Old-toolbox reference** (`base_package/gdsge/source/gdsge_parser.m`):
- Macro passes, in order, lines 45–220: `GNDSGE→GDSGE` + `process_deprecate` (45–47) →
  `cinclude("…")` (49–75) → `cinclude <…>` (77–102) → `include("…")` (104–128) → `#define`
  (130–160) → `#strcat_comma{}` (162–165) → `#mat{}` (167–171) → `#foreach…#endfor`
  (173–174) → `#for…#end` (176–217) → `#if…#endif` (219–220).
- Helpers (3025-line file tail): `process_deprecate` (2913), `rec_extract_loop_seg`
  (`#foreach`, recursive, 2944), `process_if_macro` (2983), `process_strcat_comma` (2999),
  `process_mat_str` (3014).
- `cinclude` collects C++ `#include` lines into `cxxIncludeString`, injected into the
  generated `.cpp` via the `GDSGE_OTHER_INCLUDE` placeholder (line 1787).

**Groundwork already in the new package:**
- `src/+gdsge/+parser/preprocess.m` — the macro guard to remove (lines 31–38); also does
  `GNDSGE→GDSGE` (line 42), which moves into the new engine.
- `src/+gdsge/+parser/parseFrontEnd.m` — the pipeline entry; `expandMacros` slots in as its
  first step; `hooksFromBlocks` assembles the `hooks` IR section (gains a `cxxIncludes` field).
- `src/+gdsge/+parser/evalSetup.m` / `defaultSetupCode.m` — the MATLAB eval sandbox reused by
  the eval-bearing macros.
- `src/+gdsge/+codegen/generateCxx.m:41` — the `GDSGE_OTHER_INCLUDE` placeholder, already
  present in `templates/cxx/mex.tpl.cpp:44`, currently filled `''`. `cinclude` wires into it
  with no template change.

## 2. Scope decisions

**In scope — the full macro surface:**
`include`, `cinclude("…")`, `cinclude <…>`, `#define`, `#strcat_comma{}`, `#mat{}`,
`#foreach…#endfor`, `#for…#end`, `#if…#endif`, and the deprecated-spelling normalization
(`GNDSGE→GDSGE`, `InterpOrder/SimuInterp/SimuResolve/ExtrapOrder`).

**Fidelity stance — "compat + minor cleanup" (owner-approved):** byte-identical expansion for
every input the old engine handled correctly (the backward-compat contract); informative
`gdsge:parser:macro*` errors where the old engine silently mangled or accepted-then-broke;
plus exactly two lifted limitations (§5). No other syntax extensions.

**Out of scope:** `var_tensor` (7e); SymPy backend & reduction-loop fusion (Phase 8);
`pre/post_jac_code` and raw `cxx` blocks beyond `cinclude` (already collected as `hooks`);
`#else`, `#if` nesting, `#for` arithmetic beyond `#(id+1)` (parity, not extended);
performance work (Phase 9).

## 3. Goal & success criteria

Phase 7d is done when:

1. `expandMacros(rawText, baseDir)` returns `{text, cxxIncludes}`; `parseFrontEnd` calls it
   first and the macro guard is gone from `preprocess`.
2. Each macro produces the agreed expansion (§4–§5), verified by direct unit tests on
   `expandMacros` (exact expanded text + `cxxIncludes`).
3. A macro-ized `HL1996_macro.gmod` fixture that generates the same declarations as
   `HL1996.gmod` parses to an IR `isequalIR` the existing HL1996 reference IR (§6).
4. A `cinclude` fixture injects its `#include` line(s) into the generated `mex_<model>.cpp`
   via `GDSGE_OTHER_INCLUDE`.
5. The error taxonomy (§7) is exercised: malformed macros raise identified, informative
   errors.
6. All HL1996 / 7a / 7b / 7c gates stay green. **No new old-toolbox golden capture.**

## 4. Architecture — `gdsge.parser.expandMacros`

New stage `expandMacros(rawText, baseDir) → struct('text',…, 'cxxIncludes', {cellstr})`,
called as the **first** step of `parseFrontEnd`, before line-cleaning:

```
parseFrontEnd(gmodText, modelName, baseDir):
  m     = gdsge.parser.expandMacros(gmodText, baseDir);   % NEW
  clean = gdsge.parser.preprocess(m.text);                % macro-guard removed
  sb    = gdsge.parser.splitBlocks(clean);
  …                                                        % unchanged
  ir.hooks.cxxIncludes = m.cxxIncludes;                   % carried into the IR
```

Internally, ordered passes as local sub-functions mirror the old engine's order **exactly**.
This order is load-bearing: macros expand on raw, comment-bearing text **before** `preprocess`
strips comments and joins continuations.

| # | Pass | Old ref | eval? |
|---|------|---------|-------|
| 0 | deprecated-spelling normalize (`GNDSGE→GDSGE`, `InterpOrder/SimuInterp/SimuResolve/ExtrapOrder`) | 45–47 | no |
| 1 | `cinclude("f")` / `cinclude <f>` — collect, strip | 49–102 | no |
| 2 | `include("f")` — `fileread` rel. `baseDir`, splice inline | 104–128 | no |
| 3 | `#define NAME VALUE` — word-boundary substitution | 130–160 | no |
| 4 | `#strcat_comma{args, range}` | 162–165, 2999 | yes (range) |
| 5 | `#mat{expr}` — eval, interpolate as string | 167–171, 3014 | yes |
| 6 | `#foreach id in v1 v2 … \n body #endfor id` — recursive | 173–174, 2944 | no |
| 7 | `#for id=range \n body #end` — recursive (cleanup, §5) | 176–217 | yes (range) |
| 8 | `#if cond \n body #endif` — eval cond, keep/drop | 219–220, 2983 | yes |

Pass 0 moves the `GNDSGE→GDSGE` rewrite out of `preprocess` so all later passes (and the
eval-bearing macros) see canonical text; `preprocess` keeps only comment-strip +
continuation-join.

## 5. The two agreed cleanups (everything else strict parity)

1. **Multi-token `#define` value.** Old split the `#define` line on spaces and errored unless
   it produced *exactly* 3 tokens, so a value containing spaces silently failed. New grammar:
   `#define NAME <rest-of-line>` — `NAME` must be a valid MATLAB identifier; the value is the
   trimmed, comment-stripped remainder of the line. Word-boundary substitution
   (`((?<=^)|(?<=\W))NAME((?=$)|(?=\W))`) is unchanged. `#define` with no `NAME` →
   `macroDefineMalformed`.

2. **Nested `#for`.** Old scanned for the first `#end` and made nested `#for` a hard error
   ("Nested #for not supported yet"). New: balanced `#for`/`#end` matching + recursive
   expansion (the shape `#foreach` already has via `rec_extract_loop_seg`). `#id` →
   `int2str(value)` and `#(id+1)` → `int2str(value+1)` substitution stays at parity (no new
   arithmetic forms).

**Strict parity (not changed):** `#if` (no `#else`, no nesting — non-greedy to first
`#endif`), `#strcat_comma` reproduces the old *implementation's* output (not its stale inline
comment), `#mat`/`#foreach` semantics, and the eval scope (literals/built-ins + `#define`d
values only — macros expand before parameter/grid eval, exactly as old).

## 6. The `cinclude` coupling (the one piece leaving the parser)

- `include("f")` is pure front-end: read `f` relative to `baseDir`, splice its text in place.
- `cinclude("f")` / `cinclude <f>`: strip the directive from the gmod; collect the names into
  a `cxxIncludes` cellstr as `#include "f"` / `#include <f>` lines.
- `cxxIncludes` is carried on a **new IR `hooks.cxxIncludes`** field (cellstr; `{}` default).
  `generateCxx` joins it into the existing `GDSGE_OTHER_INCLUDE` placeholder
  (`generateCxx.m:41`) — **no template change** (`mex.tpl.cpp:44` already has the slot).
- `baseDir`: `parseFrontEnd` gains an optional trailing `baseDir` arg (default `pwd`, matching
  the old cwd-relative `fileread`), threaded from `gdsge.codegen.codegen`, which knows the
  gmod's directory. Existing 2-arg `parseFrontEnd(text, name)` callers are unaffected.

## 7. Testing — synthetic + IR-equivalence, zero new capture

The macro engine is pure gmod-text→gmod-text (plus `cxxIncludes`). That lets 7d verify by
**equivalence to already-covered models** instead of fresh old-toolbox captures. Harness
unchanged: headless `matlab.unittest` via `matlab -batch "cd('tests'); run_tests"`;
`junit.xml` + exit code authoritative; one-source-per-process path policy.

- **Unit gate — `expandMacros` directly.** Synthetic snippets → assert exact expanded `text`
  and `cxxIncludes`, covering every macro, both cleanups (§5), and the error taxonomy (§7).
  Hand-authored expectations encode the agreed semantics.
- **IR-equivalence gate — the differential.** A `HL1996_macro.gmod` fixture uses
  `#define`/`#for`/`#foreach`/`#mat` to generate the **same declarations** as `HL1996.gmod`;
  assert `parseFrontEnd(HL1996_macro)` `isequalIR` the existing HL1996 reference IR
  (`tests/HeatonLucas1996/golden/`). Plain HL1996 is already differential-verified end-to-end
  against the old toolbox, so the macro fixture is covered **transitively** — no solve, no
  capture.
- **`cinclude` structural gate.** A fixture with `cinclude` of a trivial no-op header asserts
  the generated `mex_<model>.cpp` contains the injected `#include` line(s) — the one thing not
  covered transitively (plain HL1996 has no includes). Structural only; no compile/solve
  required.
- **Regression.** All HL1996 / 7a / 7b / 7c gates stay green.

This makes 7d the first widening sub-phase needing **zero new golden capture**, because macros
reduce to models already covered. (If belt-and-suspenders is later wanted, one fresh
old-toolbox capture of a macro fixture could be added — but it is redundant with the
IR-equivalence argument and is therefore omitted.)

## 8. Error taxonomy (`gdsge:parser:macro*`)

Where the old engine silently mangled input, 7d raises informative, identified errors; where
it already errored, 7d matches with a clearer message:

| Identifier | Trigger |
|---|---|
| `gdsge:parser:macroDefineMalformed` | `#define` with no valid `NAME` |
| `gdsge:parser:macroForUnterminated` | `#for` with no balanced `#end` |
| `gdsge:parser:macroForeachUnterminated` | `#foreach id` with no matching `#endfor id` |
| `gdsge:parser:macroIfUnterminated` | `#if` with no `#endif` |
| `gdsge:parser:macroEvalFailed` | `#mat` / `#if` / range eval throws (wraps the cause) |
| `gdsge:parser:macroIncludeNotFound` | `include` / `cinclude` file unreadable |

## 9. Approach & sequencing

Established TDD + differential discipline, low-risk first, on a dedicated branch
`phase7d-macro-engine`, merged when green:

1. **Scaffold `expandMacros`** with pass 0 + a pass-through, wire it into `parseFrontEnd`
   (optional `baseDir`), remove the `preprocess` macro guard and `GNDSGE` rewrite. Existing
   gates must stay green (no-op for macro-free models).
2. **Text-only passes** (`#define`, `#strcat_comma`, `#mat`, `#foreach`, `#for`, `#if`,
   `include`) with unit tests per pass, including both cleanups and the error taxonomy.
3. **`cinclude`** → `hooks.cxxIncludes` → `generateCxx`; structural gate.
4. **IR-equivalence gate**: author `HL1996_macro.gmod`; assert `isequalIR` the reference IR.
5. Full-suite green; commit the IR snapshot for the macro fixture; update `PROGRESS.md`.

## 10. Result-struct & API contract

- `parseFrontEnd(text, name)` keeps working; the new optional `baseDir` defaults to `pwd`.
- The IR gains `hooks.cxxIncludes` (cellstr, `{}` default); all other IR sections and the
  `IterRslt`/`SimuRslt` shapes are unchanged. JSON round-trip covers the new field via the
  existing schema-descriptor machinery (descriptor entry added).
- `generateCxx` output is byte-identical for macro-free / include-free models (empty
  `cxxIncludes` → empty `GDSGE_OTHER_INCLUDE`, as today).

## 11. Risks

- **Macro/comment/continuation ordering.** Old expands macros on raw text, then strips
  comments; new mirrors this (expand first, `preprocess` after). The multi-token `#define`
  cleanup strips a trailing `%…` comment from the captured value so comments can't leak in.
  Mitigation: unit tests pin the ordering on fixtures with comments inside/around macros.
- **Recursive `#for` matching.** Balanced `#for`/`#end` scanning replaces the old naive
  first-`#end` scan. Mitigation: explicit nested-`#for` unit tests plus a single-level test
  proving parity with the old expansion.
- **IR-equivalence fixture drift.** `HL1996_macro.gmod` must expand to *exactly* HL1996's
  declarations or the `isequalIR` assertion fails loudly — which is the desired tripwire, not
  a risk to correctness. Mitigation: build it incrementally against the reference IR.
- **`baseDir` resolution.** `include`/`cinclude` resolve relative to `baseDir` (= the gmod's
  directory under `codegen`, `pwd` otherwise), matching the old cwd-relative `fileread`.
  Mitigation: the include fixtures live beside their gmod; the structural gate runs under a
  controlled cwd.
