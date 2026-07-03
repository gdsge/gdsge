# Design: MATLAB passthrough in the gmod declaration region

**Date:** 2026-07-03
**Status:** Approved
**Motivating model:** `D:\macro_local\exchange_rate_dynamics_debug\incomp.gmod` (runs on the
old toolbox, rejected by the new parser).

## 1. Problem and goal

The old parser copied any unrecognized declaration-region line verbatim into the generated
MATLAB (`gdsge_parser.m`: the `otherwise → code2/code3` arms), so the preamble was
effectively arbitrary MATLAB. The new parser narrowed acceptance to single-identifier
assignments: `assignLHS` in `src/+gdsge/+parser/parseVarDecls.m` requires
`name = rhs` (identifier optionally followed by one `(...)`/`[...]` index), and Pass B
raises `gdsge:parser:unknownDeclaration` for anything else. That rejects legal old-toolbox
gmod code:

- multi-output assignments — `[uip_shock_trans, uip_shock_grid] = rouwen(...)`,
  `[z,z_hat,uip_shock,omega] = ndgrid(...)`, `[~,~,~,is_crisis_state] = ndgrid(...)`
- bare expression / function-call statements — `rng(0);`
- struct-field LHS — `opts.x = 1;`
- control flow — `for` / `if` / `end` in the preamble

Separately, Pass C raises `gdsge:parser:reopenedDeclBlock` when a declaration kind is
re-opened after other statements. The old toolbox allowed declarations anywhere;
`incomp.gmod` re-opens `parameters` three times (source lines 20, 37, 91) and `var_policy`
twice (131, 148).

**Goal:** any preamble MATLAB the old toolbox ran, the new toolbox runs — with better
diagnostics than blind copying — and `incomp.gmod` solves end-to-end with results matching
the old toolbox's saved `IterRslt_incomp.mat`.

**Non-goals:** no new gmod syntax (no explicit `matlab; ... end;` verbatim block — can be
added later, purely additive); no change to model-block parsing; no IR schema change.

## 2. Approach decision

Two options were considered for widening acceptance:

1. **Permissive residual with guardrails (chosen).** Anything that is not a recognized
   declaration and not a routed single assignment passes through as setup text. Safety
   comes from a near-miss keyword guard plus the existing parse-time eval, not from
   enumerating statement forms.
2. **Enumerated statement forms (rejected).** Regex-recognize multi-output assignments and
   bare calls specifically. Crisper rejection messages, but brittle: every unenumerated
   legal form (struct-field LHS, control flow, …) is another old-toolbox model the new
   parser rejects. Backward compatibility is a hard constraint of the refactor, so
   enumeration loses.

The passthrough plumbing already exists and is untouched: setup statements are grouped
into ordered `{kind, body}` sections (`ir.setup`), eval'd once at parse time
(`parseDeclarations.m` → `evalSetup.m`) to obtain numeric param/shock/grid values, and
replayed verbatim by `emitSetup.m` into generated `iter_*.m`/`simulate_*.m` after the
default flags. This design only changes the acceptance filter.

## 3. Parser changes (all in `src/+gdsge/+parser/parseVarDecls.m`)

### 3.1 Pass B — permissive residual

When `assignLHS` returns empty (the statement is not a `name = rhs` assignment):

1. **Near-miss keyword guard.** If the statement is command-form — an identifier followed
   by whitespace and more identifiers, with no `=` — and the first word is within
   Damerau/Levenshtein edit distance 2 of a declaration keyword (`DECL_KINDS` plus
   `inbound`, `inbound_init`, `initial`), raise a new error
   `gdsge:parser:probableTypo` naming the offending word and the suggested keyword.
   This catches `paramters beta gamma;` at parse time; the old toolbox copied it verbatim
   and failed later at eval with a confusing message.
2. **Otherwise: passthrough.** Append the statement to `setupStmts`, exactly as accepted
   single assignments are today. No error.

Single-assignment routing is unchanged: LHS ∈ interp names → `interpUpdate`; LHS ∈ tensor
names → `tensorAssign`; else `setupStmts`, with state-name LHS also recorded as
`gridText`.

Control flow works without special handling: `splitStatements` yields `for i=1:3`,
`x(i)=i`, `end` as consecutive statements in source order; `x(i)=i` routes as a normal
assignment and the loop header/footer pass through; `joinBody` re-terminates each with
`;`, which is valid MATLAB (`for i=1:3;`), so both the parse-time eval and the generated
replay reassemble the loop correctly. This relies on section bodies preserving source
order, which Pass C guarantees.

**Documented limitation (unchanged behavior):** a state grid or interp update assigned via
multi-output LHS (e.g. `[K, K_hat] = deal(...)`) is not captured for those special roles —
only plain `name = rhs` assignments are. The old corpus does not use this form.

### 3.2 Pass C — allow re-opened declaration blocks

Remove the `seenKinds` tracking and the `gdsge:parser:reopenedDeclBlock` error. A
declaration keyword that re-opens a previously seen kind flushes the current section and
starts a new section of that kind. `ir.setup` becomes a positional list in which kinds may
repeat.

No downstream change is needed: the IR schema already types `setup` as a list of
`{kind, body}` with no uniqueness constraint (no `irVersion` bump), and
`emitSetup.m` already replays section bodies positionally. The
`gdsge:parser:setupBlockMismatch` **warning** (statement positioned under a different
block than its LHS suggests) stays as-is — a non-fatal nudge toward canonical layout.
`incomp.gmod` will emit a few of these (e.g. `b_last = [b_min,b_max]` appears before its
`var_state` declaration); that is acceptable.

### 3.3 Validation backstop (already exists — no change)

The concatenated preamble is eval'd once at parse time; genuinely broken MATLAB fails
there as `gdsge:parser:setupEvalFailed`, and `evalSetup.m` already echoes the numbered
script lines plus the underlying MATLAB error.

## 4. Semantics of passthrough code (documented contract, not a change)

Preamble passthrough code runs **twice**, as in the old toolbox:

1. at parse time, inside `evalSetup`, to produce the numeric values stored in the IR
   (`params{i}.value`, `shocks.count`, `shocks.transitions`, …);
2. at runtime, replayed inside generated `iter_*` / `simulate_*` after the default option
   flags — so gmod assignments override defaults.

Consequences: preamble code must be side-effect-safe to run twice, and helper functions it
calls (`rouwen`, `bivariateNormalProbMidpoints`, …) must be on the MATLAB path both when
`gdsge_codegen` runs and when the generated files run. Same contract as the old toolbox.

## 5. Testing

TDD; synthetic fixtures only — the motivating research model is not committed.

**Parser unit tests** (new small gmod fixtures under the existing parser-test layout):

- multi-output assignment, including `~` placeholders, reaches `ir.setup` and evals;
- bare function-call statement (`rng(0);`) passes through;
- struct-field LHS (`opts.x = 1;`) passes through;
- `for` loop spanning several preamble lines survives split → section → replay → eval;
- re-opened `parameters` / `var_policy` blocks parse, with `ir.setup` sections in source
  order and duplicate kinds present;
- `paramters beta;` → `gdsge:parser:probableTypo` naming `parameters`;
- a legitimate command-form call far from any keyword (e.g. `format long`) passes through
  (no false positive);
- broken MATLAB in the preamble → existing `gdsge:parser:setupEvalFailed` with numbered
  lines.

**Codegen check:** generated setup replay preserves source order across re-opened sections
(targeted assertion or snapshot).

**Regression:** full suite via `matlab -batch "cd('tests'); run_tests"`; existing model IR
snapshots and codegen goldens must be unchanged (no existing corpus model re-opens blocks
or uses non-assignment preamble statements).

**End-to-end acceptance (manual, outside the repo):** replicate `main_incomp.m` on the new
toolbox — `gdsge_codegen('incomp')`, warm-start from the old `IterRslt_incomp.mat`
(`options.WarmUp.asg_interp_struct`, `SkipModelInit=1`, `AsgMaxLevel=3`,
`TolEq=1e-5`). The saved `.mat` cannot be regenerated (it took a long, heavily fine-tuned
old-toolbox run), so it is a warm-up input, not a reproduction target: success is that the
first iteration(s) from the warm start solve cleanly and stay close to the warm-up
solution (small `Metric`, policies near the interpolated warm-up values) — the same check
`main_incomp.m` embodies. No full from-scratch convergence run.

## 6. Known downstream risks (in scope; fixed as separate small steps if they surface)

Per the agreed scope, success is `incomp.gmod` running end-to-end; the parser widening is
the known fix, and further gaps get their own small fix-with-test steps:

- `incomp.gmod` line 260 uses `GNDSGE_INTERP_VEC` (the `GNDSGE_` spelling); confirm the
  new model-block parser accepts that alias of `GDSGE_INTERP_VEC`.
- Warm-starting the new `iter_incomp` from an **old-toolbox** `asg_interp_struct` is the
  first real cross-toolbox warm-start; the struct shape is frozen by the compat contract
  but this path is untested.
- Scale: 54 shock states × 5 continuous states with ASG — perf or convergence surprises
  possible; pin `rng` before iter if the solve shows restart sensitivity.
- `bivariateNormalProbMidpoints` is not in the model folder; locate it on the user's
  MATLAB path at implementation time (required for the parse-time eval, as it was for the
  old toolbox).
