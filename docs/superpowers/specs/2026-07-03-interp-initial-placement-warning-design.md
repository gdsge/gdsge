# Exempt interp `initial`/update statements from placement warnings — design

**Date:** 2026-07-03
**Status:** Approved.
**Area:** `src/+gdsge/+parser/parseVarDecls.m` (positional-sections warning pass).
**Amends:** `2026-06-16-decl-region-positional-sections-design.md` (expected-section table).

## Problem

The canonical GDSGE layout declares `var_interp` early, then `var_policy_init` /
`inbound_init` / `var_aux_init`, then the `model_init; ... end;` block, and only *after*
it the `initial X_interp EXPR;` and `X_interp = EXPR;` statements. That placement is
semantically forced: the `initial` expressions reference `var_aux_init` names that
`model_init` computes, so they cannot sit directly under the `var_interp` declaration.

But the positional-sections pass never sees the `model_init` block — `splitBlocks`
excises it before `parseVarDecls` runs. So every post-`model_init` interp statement lands
positionally in the `var_aux_init` section, and `expectedSection()` says it belongs to
`var_interp`, producing one `gdsge:parser:setupBlockMismatch` warning per statement
(56 of them on a real user model, `debug_gdsge_re`). The warning fires on the
documented-correct layout — a false positive.

## Key fact

Interp `initial` and interp-name assignments are extracted by keyword / LHS name in
Pass A/B **regardless of position** (`interpInitial` / `interpUpdate`) and never enter a
section body. Their placement is functionally inert. The design already exempts
`inbound` / `inbound_init` from placement warnings for exactly this reason ("low signal;
structurally bound regardless of section"); interp initial/updates have the same
character.

## Change

In `expectedSection()` (`parseVarDecls.m`), the two returns that currently claim
`var_interp`:

- `initial X expr` keyword statements, and
- assignments whose LHS is a declared interp name,

become "no expectation" (`''`), joining `inbound` / `inbound_init` in the not-warned
category. Everything else — the sectioning pass, Pass A/B extraction, eval order, IR
sections — is untouched, so no source goldens or IR snapshots change.

## Testing

One new test in `tests/parser/tParseVarDecls.m`: a decl text with the canonical shape

```
var_interp Ev;
var_aux_init E;
initial Ev E;
Ev = E;
```

parses with **no** `setupBlockMismatch` warning while still populating
`interpInitial` / `interpUpdate`. Existing warning tests (`shock_trans` misplacement,
actionable-message) stay as-is and must keep passing.

## Docs

Update the expected-section table in
`docs/superpowers/specs/2026-06-16-decl-region-positional-sections-design.md`: drop the
`var_interp` row (declared interp LHS, `initial X`) and fold it into the not-warned note
alongside `inbound` / `inbound_init`.

## Verification

- `matlab -batch "cd('tests'); run_tests"` — full suite passes.
- Re-run `gdsge_codegen` on
  `D:\macro_local\debug_gdsge_re\large13calibrationv5utilitiesCCPsREEtweak.gmod` and
  confirm zero `setupBlockMismatch` warnings.

## Risk

Minimal. Warning-only behavior change; no eval, IR, or codegen output differs. The lost
signal (an interp assignment placed under, say, `var_state`) was cosmetic — the statement
is still extracted and used correctly wherever it appears.
