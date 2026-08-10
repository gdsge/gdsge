# Both simulate variants are always generated; the mode toggle only rewrites the entry .m

`SIMU_RESOLVE`/`SIMU_INTERP` used to gate whether `simulate_<model>_mex.cpp` was
emitted at all, so toggling the mode changed the concatenated C++ text and forced a
full recompile of both MEX binaries. We decided to emit the simulate MEX source
whenever the model is expressible, **regardless of mode** (expressibility no longer
consults `options.simuInterp`), keeping exactly one public `simulate_<model>.m`
whose content is the selected variant. The C++ text — and therefore the verbatim
string-compare compile cache — becomes invariant to the toggle: flipping the option
re-runs codegen but only rewrites the entry `.m`; no recompilation occurs once
cached.

## Consequences

- SIMU_RESOLVE-only expressible models now compile `simulate_<model>_mex` once,
  even if it is never called — the accepted price for free mode switching.
- A future reader will find a simulate MEX binary next to models that never use
  it; that is deliberate, not a leftover.
- Rejected: three-file layout (`_resolve.m`/`_interp.m` + thin entry) — no
  precedent for generated helper `.m` files; and runtime mode selection via the
  options whitelist — a public-API change the entry-rewrite approach avoids.
