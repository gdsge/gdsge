# Declaration-region positional section model — design

**Date:** 2026-06-16
**Status:** Approved (design); implementation plan pending.
**Area:** `src/+gdsge/+parser/` (front end) + `src/+gdsge/+codegen/+mat/emitSetup.m` + IR schema.

## Problem

The declaration region of a `.gmod` (everything outside the named `model`/`simulate`/
hook blocks) is currently classified by **inspecting each statement**:

- `parseVarDecls.m` Pass A classifies any line whose leading word is a declare keyword
  (`parameters`, `var_state`, …).
- Pass B classifies every remaining assignment by its **LHS name** — LHS in `interpNames`
  → interp update; in `tensorNames` → tensor assign; in `stateNames` → grid text; otherwise
  it is appended to one monolithic `setupStmts` blob.
- That blob becomes `ir.setup` (a single text string), eval'd all at once and replayed
  **verbatim in one place** by `emitSetup.m`.

We want a **positional** model instead: the statements between two declaration statements
belong to the **preceding** declaration's block, and each block feeds a distinct coding
block (e.g. `parameters` → params section, `var_shock` → shock section defining
`shock_num`/`shock_trans`/shock values, `var_state` → grid section). Position becomes the
authoritative classifier.

## Decisions (locked)

1. **Position is authoritative**, but backward compatibility is preserved by *warning*
   (not erroring) when a known statement lands in a block different from the one its
   declaration is in. The statement is **not** physically moved.
2. **IR representation:** replace the single `ir.setup` text blob with a structured,
   ordered list of sections. (IR schema version bumps; golden IR JSONs regenerate.)
3. **Re-opening a declaration kind non-contiguously is a hard error.** No current corpus
   model does this, and the contiguity guarantee is what lets block-grouping never reorder
   statements.

## Core concept — sections

The declaration region becomes an **ordered list of sections**. A section is:

- a **header** — one or more *consecutive* declaration lines of the **same kind**, plus
- a **body** — every statement after the header up to the next header of a *different* kind.

**Declaration kinds (block openers):** `parameters`, `var_shock`, `var_state`,
`var_policy`, `var_aux`, `var_interp`, `var_tensor`, `var_output`, `var_others`,
`var_policy_init`, `var_aux_init`.

**Body statements** (everything else): plain assignments (parameter values, shock values,
`shock_num`/`shock_trans`, state grids, interp updates, tensor assigns, toolbox options,
helper intermediates) and the keyword-led associates `inbound` / `inbound_init` /
`initial`.

Statements appearing **before the first declaration** form an optional **`global`
(preamble)** section, emitted first.

## Parsing algorithm

Replaces the Pass-A/Pass-B structure in `parseVarDecls.m`. Single forward pass over
`splitStatements` output, tracking `currentSection`, `headerOpen` (can the current header
still extend?), and a `seenKinds` set:

- **Declaration line of kind K:**
  - if `headerOpen && K == currentSection.kind` → extend the current header.
  - else if `K ∈ seenKinds` → **hard error** `gdsge:parser:reopenedDeclBlock`
    ("declaration kind 'K' re-opened; make all 'K' declarations contiguous").
  - else → close the current section, open a new one of kind K, `seenKinds += K`,
    `headerOpen = true`.
- **Body statement:** append to `currentSection.body`; set `headerOpen = false`.

Comments and blank lines are already stripped by `preprocess.m` (comment cut at first `%`,
continuations joined) before statements are split, so they never appear as statements and
never break header contiguity.

This guarantees each kind maps to exactly one **contiguous** source region, so grouping by
section never reorders statements.

## Semantic validation — warnings, never reroute

After sectioning, each body statement gets an **expected** section computed from its LHS:

| Statement LHS / form                         | Expected section |
|----------------------------------------------|------------------|
| a declared parameter name                    | `parameters`     |
| a declared shock var name, `shock_num`, `shock_trans` | `var_shock` |
| a declared state name (grid assignment)      | `var_state`      |
| a declared `var_interp` name, or `initial X` | `var_interp`     |
| a declared `var_tensor` name                 | `var_tensor`     |
| unknown LHS (options, helper intermediates)  | *(none)*         |

- `shock_trans` / `shock_num` always count as shock-associated even if the name was also
  declared elsewhere (handles CaoNie's `var_others shock_trans;`).
- Unknown-LHS statements (e.g. `USE_ASG`, `TolEq`, CaoKS's `EZTrans` intermediates) have
  no expectation and never warn.
- `inbound` / `inbound_init` **placement is not warned** — low signal; they are already
  structurally bound to policy/aux bounds regardless of section.

If a statement's actual section ≠ its expected section → emit
`gdsge:parser:setupBlockMismatch`. **The warning must be actionable**: it names the
variable, the block it is currently in, the block it belongs to, the source line, and a
concrete fix instruction telling the author to move the statement under the right
declaration. For example:

```
gdsge:parser:setupBlockMismatch: 'shock_trans' (line 41) is defined in the 'var_state'
block but is associated with 'var_shock'. Move this assignment so it appears after the
'var_shock' declaration (and before the next declaration) to silence this warning.
```

**The statement is not moved** — it stays in its source position for both eval and code
emission.

## Eval — identical to today

The eval text is the source-ordered concatenation of all section bodies, minus
interp-updates and tensor-assigns (still carried as opaque text, exactly as now). This
reproduces today's `setupBody` byte-for-byte, so **numeric parity is guaranteed** for every
existing model. Grid / interp / tensor text extraction is unchanged, now additionally
cross-checked for placement.

## IR shape — replace `ir.setup`

`ir.setup` (single text blob) → an **ordered list of sections**, each:

```
{ kind:  'parameters' | 'var_shock' | ... | 'global',
  names: { declared names in this section's header },   % {} for global / empty headers
  body:  '<verbatim, eval-able body text in source order>' }
```

- The source-ordered concatenation of every `body` (minus interp/tensor text) equals the
  old `ir.setup` eval text.
- `ir.setupNames` (workspace-derived; consumed by `optionsWhitelist`) is unchanged.

**Ripple:** `schema.m`, `ir-schema.md`, `ir-versioning.md` (version bump),
`assemblePartialIR.m`, and **every committed `*.gdsge.json` golden** (regenerate). The C++
backend does **not** read `ir.setup` (verified — only `emitSetup.m` does), so it is
unaffected.

## Codegen — `emitSetup.m` emits labeled coding blocks

Rewrite `emitSetup.m` to walk the section list and emit each body under a `%% ---- <kind>`
comment, in source order. Net generated text = today's setup replay + the section-header
comments; numeric behavior is unchanged. Regenerate `iter_HL1996_golden.txt` /
`simulate_HL1996_golden.txt`; update `tEmitSetup.m`.

## Edge cases (all observed in the real corpus)

- **Preamble** statements before the first declaration → `global` section, emitted first.
- **Empty declaration** line (`var_tensor;` in CaoNie) → opens a section with an empty
  header name set.
- **`shock_num` / `shock_trans` before `var_shock`** (CaoNie, CaoKS, Cao2011EZ) → land in
  the preceding section → warn, still eval correctly.
- **Parameter values / `shock_trans2` assigned in a later section** (Cao2011EZ `ql,qh,pl,ph`
  in shock section; CaoKS `shock_trans2 = shock_trans;`) → warn, still eval correctly.
- **`var_others shock_trans`** (CaoNie) → allowed; `shock_trans` keeps shock precedence.
- **Interspersed toolbox options** (CaoKS `NumThreads`, `USE_ASG`, `TolEq`, …) → unknown
  LHS, no warning; eval'd and read by `resolveOptions` as today.
- **Comments / blank lines between declarations** → stripped by `preprocess`; contiguity
  safe.
- **`inbound` expressions referencing setup vars** (safe_assets `inbound Rf Re_n(2)
  Re_n(1)`) → unchanged.
- **Non-canonical kind order** (Bianchi declares `var_state` before `var_shock`) → allowed;
  each kind still opens exactly once.
- **Contiguous multi-line same-kind declarations** (Mendoza's two `parameters` lines) →
  one header.
- **Re-opened kind** (none in corpus) → hard error.

## Testing

- **Unit tests for the sectioning pass:** canonical order; reversed kind order
  (state before shock); re-open hard error; each misplacement warning (captured via
  `lastwarn`); empty declaration; preamble.
- Update `tests/parser/tParseVarDecls.m` and `tests/codegen/tEmitSetup.m`.
- Regenerate golden source snapshots (`iter_HL1996_golden.txt`,
  `simulate_HL1996_golden.txt`) and all `*.gdsge.json` IR goldens.
- Run all end-to-end gates (HL1996, Mendoza, GLSW, Bianchi, CaoKS[/simu_interp], CaoNie,
  Cao2011EZ, safe_assets) to confirm: numeric parity preserved **and** misplacement
  warnings fire on Cao2011EZ / CaoKS / CaoNie while those models still solve.

## Risk

Low. Eval order is preserved exactly (source-ordered body concatenation), so numeric
results cannot change. The visible surface is the IR schema version bump plus regenerated
source/IR goldens, and the new parser warnings. The only hard rejection added
(non-contiguous re-open) matches no existing corpus model.
