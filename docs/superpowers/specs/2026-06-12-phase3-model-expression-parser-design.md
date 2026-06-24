# Phase 3 — Model-Expression Parser — Design Spec

- **Date:** 2026-06-12
- **Status:** Approved
- **Parent spec:** `2026-06-11-refactor-gdsge-design.md` (§7 stages 4–5, §15 Phase 3)

---

## 1. Goal

Fill the one hole the Phase-2 front-end left: `ir.model`. Phase 3 parses the
`model;…equations;…end;…end;` body into typed statements and equations, runs semantic
analysis, and produces the **complete IR for HL1996**.

**The phase gate:** `gdsge.parser.parseFrontEnd` on `tests/HeatonLucas1996/HL1996.gmod`
must satisfy `gdsge.ir.isequalIR` against the hand-authored reference
`tests/HeatonLucas1996/ir/buildHL1996IR.m`, and round-trip against the regenerated JSON
golden `HL1996.gdsge.json`.

## 2. Scope

**In scope (approved choice: "full expressions, HL1996 statements"):**

- The **complete expression grammar** — all 7 schema node kinds (`num`, `name`, `primed`,
  `unop`, `binop`, `call`, `index`), full MATLAB operator precedence, parentheses.
- All **4 reduction kinds** (`EXPECT`, `MIN`, `MAX`, `PROD`) and the
  `GDSGE_*{body | transName}` custom-transition syntax (default `shock_trans`).
- The statement surface HL1996 needs: plain/primed assignment, the primed
  `GDSGE_INTERP_VEC'` multi-target call, reduction assignment, inline-reduction hoisting,
  the `equations;…end;` residual list with primed (future-indexed) equations.
- Semantic analysis: name resolution, square-system check, interp-call arity,
  `transRef` validity.

**Out of scope (Phase 7):** `model_init` bodies, conditional `model(cond);` blocks,
named-interp calls (`ps_future'(x)`), `OUTPUT_CONSTRUCT_CODE`/`TASK` expansion, macros,
absolute file-line error positions, string-aware tokenizing.

## 3. Module layout

```
sb.blocks.model ──▶ parseModel ──▶ {statements, equations} ──▶ analyzeModel(ir) ──▶ full IR
```

New files in `src/+gdsge/+parser/`, one stage per file:

| File | Purpose |
|---|---|
| `tokenize.m` | model-body text → token array (type, text, line/col) |
| `parseExpr.m` | token stream → expression AST (the 7 node kinds) |
| `parseModel.m` | model-block text → ordered statement list + equation list |
| `analyzeModel.m` | semantic checks over the assembled full IR |

Small edits: `parseFrontEnd.m` calls `parseModel` + `analyzeModel`;
`assemblePartialIR.m` accepts the model section and drops its "partial" caveat. The
Phase-1 schema, node constructors, and validator are **unchanged**.

## 4. Tokenizer (`tokenize.m`)

Single character-driven pass over the model-block text (comments already stripped,
continuations already joined by `preprocess`). Token types:

- `name` (`[A-Za-z]\w*`), `number` (decimal with optional fraction/exponent: `1`, `0.95`,
  `1e-3`, `1.5e+2`).
- Operators: `+ - * / ^ .* ./ .^ == ~= <= >= < > && || & | ~`.
- Punctuation: `( ) [ ] { } , ; =`. The reduction braces `{ }` and the pipe `|` are
  ordinary tokens; the *parser* disambiguates pipe-as-or vs pipe-as-transition-separator
  by context (inside a reduction brace, a depth-0 `|` separates the body from the
  transition name).
- `prime` (`'`). In gmod model bodies `'` is always the future marker — never transpose,
  never a string quote — so no ambiguity handling is needed.

Each token records line and column **relative to the model block**. Unknown characters
raise `gdsge:parser:badToken` with that position.

## 5. Expression grammar (`parseExpr.m`)

Recursive descent with precedence climbing. MATLAB semantics:

- Binary precedence, low → high: `||` < `&&` < `|` < `&` < comparisons
  (`== ~= < > <= >=`) < additive (`+ -`) < multiplicative (`* / .* ./`) < power (`^ .^`).
- **All binary operators are left-associative, including `^`** (MATLAB: `2^3^2 = 64`).
- Unary `- + ~` bind between power and multiplicative (MATLAB: `-2^2 = -4` because power
  binds tighter than unary, but `-a*b` is `(-a)*b` because unary binds tighter than `*`).
  Power additionally accepts a unary expression as its immediate right operand (`2^-3`).
- Postfix `'` on a name → `primed` node. `name(args)` → `call` node. Call vs index is
  syntactically identical in MATLAB; the parser emits `call` uniformly (`index` stays in
  the schema for future semantic reclassification; HL1996 uses neither).
- Elementwise operators are **normalized to scalar forms** in the AST (`.*`→`*`, `./`→`/`,
  `.^`→`^`): model bodies are scalar per grid point, and the IR stays backend-friendly.

Errors: `gdsge:parser:badExpr` with the offending token and its position.

## 6. Statement layer (`parseModel.m`)

1. Extract the nested `equations;…end;` region first (`splitBlocks` deliberately leaves
   it inside the model body), then split both regions with the existing
   `splitStatements`.
2. Classify each model statement by leading tokens:
   - `[a',b',…] = GDSGE_INTERP_VEC'(args)` → `interpCall{targets, primed, args,
     interpRef}`. The `'` on the function name sets `primed`; every target must carry a
     matching prime.
   - `target = GDSGE_<KIND>{body | trans}` with the RHS exactly one reduction →
     `reduction{kind, target, body, transRef}`; `transRef` defaults to `shock_trans`.
   - `target['] = expr` → `assign{target, primed, expr}`.
3. **Inline-reduction hoisting:** a reduction nested inside a larger expression becomes
   its own `reduction` statement inserted immediately **before** the host statement; the
   call site is replaced by a generated name **`GDSGE_<KIND>_<n>`** (`GDSGE_EXPECT_1`, …)
   with one global counter across kinds. The `GDSGE_` prefix is the toolbox's reserved
   namespace (collision-free, legal in MATLAB and C++).
4. Equation statements parse as residual expressions with a `primed` flag (a trailing
   `'` marks a future-indexed equation, e.g. `w1_consis'`).
5. Statement order is **source order** (hoisted reductions directly before their hosts).

Parser errors (matching old-toolbox behavior where it had any): more than one reduction
in a statement, nested reductions, unterminated `{`, malformed interp-call target lists.

**Fixture change:** the reference IR's hoisted temp `__ep_expect` is renamed
`GDSGE_EXPECT_1` (`__`-prefixed identifiers are reserved in C++), in
`buildHL1996IR.m` and the regenerated `HL1996.gdsge.json`.

## 7. Semantic analysis (`analyzeModel.m`)

Runs on the assembled full IR — parser-side, per the Phase-1 split (validator = shape;
parser = model semantics):

1. **Name resolution.** Every `name`/`primed` leaf resolves to: param, state, shock,
   policy, aux, interp-call output, reduction target, or a model-local defined by an
   earlier `assign`. Built-ins are whitelisted: `GDSGE_Iter`, `TASK`, and `shock` (the
   current shock index, available to model code in the old toolbox).
   Unknown names → `gdsge:parser:unknownName` naming the statement. Primed
   names must be shock names, interp-call outputs, or policy variables.
2. **Square system.** Unknown count = sum of policy slot widths (HL1996: 19). Equation
   count = plain equations + `shock_num` per primed equation (HL1996: 11 + 8 = 19).
   Mismatch → `gdsge:parser:notSquare` reporting both counts.
3. **Interp-call arity.** `interpCall` target count must equal the number of declared
   `var_interp` objects, in declaration order (HL1996: 4).
4. **Reduction `transRef`** must name a declared transition matrix.

Bounds-completeness stays in `assemblePartialIR` (Phase 2).

## 8. Testing (TDD order)

1. `tests/parser/tTokenize.m` — token kinds, positions, number forms, primes, `{}`/`|`,
   bad-character error.
2. `tests/parser/tParseExpr.m` — precedence/associativity table (`-2^2`, `2^3^2`,
   `a/b/c`, `2^-3`), primes, calls, parens, elementwise normalization; node-equality
   against hand-built `gdsge.ir.node.*` trees, including the exact `es1`, `budget_1`,
   `w1_consis` expressions from the fixture.
3. `tests/parser/tParseModel.m` — classification of every HL1996 statement form;
   hoisting position + generated name; reduction pipe syntax; error cases.
4. `tests/parser/tAnalyzeModel.m` — mutated-IR negatives: unknown name, non-square
   system, wrong interp arity, bad transRef.
5. `tests/HeatonLucas1996/parser/tFullIRHL1996.m` — **the phase gate** (§1).

All errors use `gdsge:parser:*` IDs with block-relative positions, consistent with
Phase 2. Run via `pwsh -File tests/run.ps1`.

## 9. Risks

- **MATLAB operator-precedence quirks** (`^` left-associativity, unary-minus binding) are
  easy to get subtly wrong — covered by a dedicated precedence test table checked against
  real MATLAB evaluations.
- **Fixture edit** (`__ep_expect` → `GDSGE_EXPECT_1`) touches the committed golden; it is
  a pure rename, and the JSON golden is regenerated by the existing round-trip tooling in
  the same commit.
