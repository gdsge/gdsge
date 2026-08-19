# Conditional statements are guard-restricted and adept-only

The old toolbox let a model body carry `if/elseif/else/end`, lowering it to literal C++
after running `replace_grid_variable` on the condition — so the branch could test a
var_policy and flip mid-Newton-iteration. Restoring the feature, we narrowed it on two
axes: a guard may reference **only var_state, var_shock, and parameters** (the same rule
`lowerCondition` already enforces for conditional equations), and the **SymPy Jacobian
backend rejects conditional statements** outright, with `resolveBackend` steering
auto-mode to adept for such models.

## Consequences

- Branch selection is constant across a solve, so the residual stays smooth *within* a
  branch and the Jacobian adept tapes is the Jacobian of the branch actually taken. A
  policy-dependent guard would have produced a silently wrong Jacobian at the switching
  boundary — the failure mode we chose to make impossible rather than document.
- No model in the corpus loses anything: model-body `if` was never used in an old
  `.gmod` (only conditional *equations* were, in CaoNie2016), so the narrowing costs no
  backward compatibility. The parity claim we are honouring is the feature's existence,
  not its exact permissiveness.
- SymPy is a backend the old toolbox never had, so rejecting it here breaks no parity.
  Supporting it means materialising each branch target's value *and* all its partial
  derivatives into outer-scope variables at every branch exit — `sympymodel/generate`
  emits helper temporaries whose C++ scope dies at the closing brace while the
  derivative registry keeps referencing them. Deferred until a real model needs it.
- `hoistCommonSubtrees` treats a conditional statement as a barrier: nothing is hoisted
  into or out of a branch, so a subtree guarded against (say) a negative logarithm is
  never lifted to the unguarded path. The memo survives the join minus anything the
  branches assign, so CSE continues for the rest of the body.
- Rejected: allowing policy-dependent guards behind an opt-in flag — it makes the wrong
  Jacobian a supported configuration; and capping nesting at one level to match
  conditional equations — the recursive IR makes nesting free, so the cap would have
  been extra code buying nothing.
