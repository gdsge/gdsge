# IR versioning policy

The JSON IR (`<model>.gdsge.json`) is the stable contract between the MATLAB front-end and the code
generators: backends never re-parse `.gmod`. `gdsge.ir.schema().irVersion` carries a semantic
version, currently **2.0.0**, and `gdsge.ir.decode` enforces major-version compatibility
(`gdsge:ir:incompatibleVersion`).

## Policy

- **Major (`N.x.x`)** — incremented only on a breaking shape change. `decode` rejects an IR whose
  major differs from the toolbox's. Generated artifacts and committed goldens are tied to a major.
- **Minor (`1.N.x`)** — additive / backward-compatible enrichment. A newer toolbox reads an older
  minor; the schema gains shapes without removing or repurposing existing ones.
- **Patch (`1.1.N`)** — editorial, no semantic change.
- **"Frozen"** means: the major-1 schema is stable. Corpus work may add minors; it never silently
  breaks major 1. Bumping to major 2 is a deliberate act that must also update
  `tests/ir/tIrVersionFreeze.m` (the tripwire).

## IR changelog

- **1.0.0** — initial IR (Phases 1–8). Sections: `params`, `setup`, `setupNames`, `variables`,
  `model` (single unconditional body), `equations`, `simulate`, `options`, `hooks`.
- **1.1.0** — Phase 9b. `model` becomes a `regions` list (`{condition, statements, ...}`);
  `equations` become a tagged `plain | conditional` list (the `eqkinds` registry). Existing
  single-region models migrate as `regions[0]` / condition `''` / `plain` tag — proven shape-only
  by `tRoundtrip/existingModelsAreSingleUnconditionalRegion`.
- **2.0.0** — declaration-region positional sections. `setup` changes from a single
  verbatim text blob to an ordered list of `{kind, body}` sections (one per declaration
  block, in source order). The concatenation of section bodies equals the old `setup`
  string, so eval semantics are unchanged. Breaking shape change → major bump; all
  committed `*.gdsge.json` goldens regenerate and `tIrVersionFreeze` flips to major 2.

See `docs/ir-schema.md` for the full field-by-field schema (auto-generated; no-drift tested).
