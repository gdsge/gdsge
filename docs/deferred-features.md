# Deferred features & honest-error map

The refactor fails fast — with a clear, namespaced error — on `.gmod` constructs and option
combinations it does not yet generate code for, rather than emitting silently-broken output. Each
is a deliberate deferral with its own future spec→gate cycle (see `PROGRESS.md`). This table is the
single source of truth for "what raises, where, and why".

| Error identifier | Raised in | Trigger | Workaround / tracked in |
|---|---|---|---|
| `gdsge:codegen:genCodeSegmentUnsupported` | `+codegen/codegen.m` | `options.GenCodeSegment` set | The thin-file architecture has no segment decomposition. Phase 7c. |
| `gdsge:parser:varTensorInBodyUnsupported` | `+parser/analyzeModel.m` | a `var_tensor` enters the model residual (the C++ `POPN` body path) | Inline the tensor expression into the equations. Phase 7e-full. |
| `gdsge:codegen:varTensorAsgUnsupported` | `+codegen/assertSupportedIR.m` | `var_tensor` together with the ASG interpolation backend | Phase 7e-full. |
| `gdsge:parser:conditionalModelInitUnsupported` | `+parser/splitBlocks.m` | a conditional `model_init(<cond>)` region | Use an unconditional `model_init`. Future. |
| `gdsge:codegen:sympyInterpUnsupported` | `+codegen/generateCxx.m` (and `+ir/validate.m`) | SymPy backend (`UseAutoDiff=0`) with the **pchip** interpolant | Use the autodiff backend (ASG is now supported); pchip is a separate gap. |
| `gdsge:codegen:unsupported` | `+codegen/generateCxx.m` | pchip interpolant, or `cxx` hook blocks, under the C++ backend | Phase 7c / 8b. |
| `gdsge:ir:incompatibleVersion` | `+ir/decode.m` | an IR JSON whose **major** version differs from `gdsge.ir.schema().irVersion` | Regenerate the IR with the matching toolbox major. See `docs/ir-versioning.md`. |

Internal "should never happen" invariants (`gdsge:codegen:sympyUnsupportedStmt`,
`gdsge:codegen:unsupported` inside the model-body emitter) guard against malformed IR reaching a
backend; they are not user-facing deferrals and indicate a bug if hit.

## Open deferred sub-phases

- **Phase 8b** — SymPy analytic-Jacobian backend: **ASG landed 2026-06-15** (CaoKS2016 /
  Bianchi2011_asg run under `UseAutoDiff=0`); **pchip** remains deferred (it has no C++ backend at
  all — `gdsge:codegen:unsupported` — a separate gap from ASG).
- **Phase 7e-full (remainder)** — C++-body `var_tensor` (the `POPN` path) and the
  `IterRslt.var_tensor` result field. The MATLAB-side `var_tensor` subset (ndgrid tensors feeding
  inbound bounds / initial interp) landed in Phase 9a.

No corpus model exercises either, so neither blocks backward compatibility.
