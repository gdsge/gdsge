# Phase 10 — Polish Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Finish the GDSGE refactor: trim the oversized golden, freeze the IR contract, measure performance vs the old toolbox, and write the documentation — leaving a clean, self-explanatory package.

**Architecture:** Four independent workstreams committed in order on branch `phase10-polish` (already created; the design spec is already committed there): Cleanup → IR version freeze → Perf harness+report → Docs. Docs go last so they describe the final tree and cite the perf numbers. No `.gmod`, IR-shape, or backend changes — this phase documents and tidies what Phases 0–9b built.

**Tech Stack:** MATLAB R2025b (`matlab.unittest`, `matlab -batch`), JSON IR (`gdsge.ir`), markdown docs. Old toolbox under `base_package/gdsge/source` (git-ignored, present). No new Python.

**Spec:** `docs/superpowers/specs/2026-06-14-phase10-polish-design.md`

**Pre-flight (run once before Task 1):**
```bash
git rev-parse --abbrev-ref HEAD     # expect: phase10-polish
git log --oneline -1                # expect the Phase 10 design-spec commit
```
If not on `phase10-polish`, `git checkout phase10-polish` first.

---

## Workstream 1 — Cleanup

### Task 1: Trim the 97 MB Mendoza2010 golden

The on-disk `IterRslt.mat` is 97 MB (130 MB in memory). Measured field sizes: `output_interp`
73.7 MB, `pp` 27.7 MB, `GDSGE_PROB` 16.4 MB — **none read by any test**. The only fields any test
reads from this golden are `{Iter, Metric, shock_num, shock_trans, var_state, var_policy, var_aux,
var_interp}` (confirmed from `tEndToEndMendoza2010.m`, `tEndToEndMendoza2010Sympy.m`,
`tGoldenMendoza2010.m`). Keep exactly those; drop the rest; re-save `-v7` (compressed).

**Files:**
- Create: `tests/golden/trim_mendoza_golden.m` (committed one-shot reproducer)
- Modify (in place, regenerated): `tests/Mendoza2010/golden/IterRslt.mat`

- [ ] **Step 1: Baseline — confirm the two Mendoza gates are green before touching the golden**

Run:
```
matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); addpath(fullfile(fileparts(pwd),'src','kernels')); r=runtests('Mendoza2010/tGoldenMendoza2010.m'); disp(table(r)); exit(double(any([r.Failed])))"
```
Expected: PASS (the integrity test loads the current big golden). Record that it is green.

- [ ] **Step 2: Write the trim reproducer**

Create `tests/golden/trim_mendoza_golden.m`:
```matlab
function trim_mendoza_golden()
% Trim tests/Mendoza2010/golden/IterRslt.mat to only the fields the Mendoza
% gates read, re-saving -v7. Drops the uncompared bulk (output_interp ~74MB,
% pp ~28MB, GDSGE_PROB ~16MB). One-shot; committed for reproducibility.
here = fileparts(mfilename('fullpath'));            % tests/golden
goldenFile = fullfile(fileparts(here), 'Mendoza2010', 'golden', 'IterRslt.mat');

g = load(goldenFile);
full = g.IterRslt;
keep = {'Iter','Metric','shock_num','shock_trans', ...
        'var_state','var_policy','var_aux','var_interp'};
IterRslt = struct();
for i = 1:numel(keep)
    assert(isfield(full, keep{i}), 'golden missing field %s', keep{i});
    IterRslt.(keep{i}) = full.(keep{i});
end

before = dir(goldenFile);
save(goldenFile, 'IterRslt', '-v7');
after = dir(goldenFile);
fprintf('Mendoza golden: %.1f MB -> %.1f MB (%d fields kept)\n', ...
    before.bytes/1e6, after.bytes/1e6, numel(keep));
end
```

- [ ] **Step 3: Run the trim**

Run:
```
matlab -batch "cd('tests/golden'); trim_mendoza_golden"
```
Expected: prints `Mendoza golden: 97.x MB -> <~9> MB (8 fields kept)`. (Floor is set by the
compared `var_policy`+`var_aux` arrays, ~11 MB in memory → single-digit MB on disk.)

- [ ] **Step 4: Re-run BOTH Mendoza gates (integrity + end-to-end, autodiff + sympy) against the trimmed golden**

Run:
```
matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); addpath(fullfile(fileparts(pwd),'src','kernels')); r=runtests({'Mendoza2010/tGoldenMendoza2010.m','Mendoza2010/codegen/tEndToEndMendoza2010.m','Mendoza2010/codegen/tEndToEndMendoza2010Sympy.m'}); disp(table(r)); exit(double(any([r.Failed])))"
```
Expected: all PASS — the comparisons read only the kept fields, so values are byte-identical to
before; only the dropped bulk is gone. (Slow: the end-to-end gates re-solve Mendoza; minutes.)

Note: `tGoldenMendoza2010.m` (the integrity test) needs **no code change** — every field it reads
(`shock_num`, `var_policy.cTildeNext`, `var_state.cTilde/k`, `Metric`, `var_aux.gdp`) is in the
kept set. The spec's "update the integrity test" turns out to be a no-op; the plan only re-runs it.

- [ ] **Step 5: Commit**

```bash
git add tests/golden/trim_mendoza_golden.m tests/Mendoza2010/golden/IterRslt.mat
git commit -m "cleanup(phase10): trim Mendoza golden 97MB->~9MB (drop uncompared output_interp/pp/GDSGE_PROB)"
```

---

### Task 2: TODO / dead-code pass

**Files:** none changed (audit + record). The scan target is first-party code only.

- [ ] **Step 1: Scan for markers**

Run:
```
grep -rnE "TODO|FIXME|XXX|HACK" src/ templates/ pyext/gdsge_sympy/
```
Expected: exactly one hit — `src/kernels/myppual.m:63` ("TODO List for the version 3.0: C++
implementation ..."), inside the **vendored** myppual kernel docstring.

- [ ] **Step 2: Triage**

Disposition: the single hit is in vendored third-party code (rewriting kernels is a non-goal per
the parent spec §3). No action. There are **no actionable first-party TODO/FIXME markers**. (No
commit for this task on its own; recorded in the Task 3 deferred-features doc and the final
PROGRESS update.)

---

### Task 3: Consolidate the deferred-error taxonomy

Collect every architecture-honest deferral the code raises into one reference table. The exact IDs
(verified by grep) are below — use these, not the stale ones from the spec draft.

**Files:**
- Create: `docs/deferred-features.md`

- [ ] **Step 1: Verify the deferral sites**

Run:
```
grep -rnE "error\('gdsge:[a-zA-Z:]*(Unsupported|genCodeSegment|incompatible)" src/
```
Expected sites (IDs to put in the table):
- `gdsge:codegen:genCodeSegmentUnsupported` — `src/+gdsge/+codegen/codegen.m`
- `gdsge:parser:varTensorInBodyUnsupported` — `src/+gdsge/+parser/analyzeModel.m`
- `gdsge:codegen:varTensorAsgUnsupported` — `src/+gdsge/+codegen/assertSupportedIR.m`
- `gdsge:parser:conditionalModelInitUnsupported` — `src/+gdsge/+parser/splitBlocks.m`
- `gdsge:codegen:sympyInterpUnsupported` — `src/+gdsge/+codegen/generateCxx.m`
- `gdsge:codegen:unsupported` (pchip interp; cxx hook blocks under the cxx backend) — `src/+gdsge/+codegen/generateCxx.m`
- `gdsge:ir:incompatibleVersion` — `src/+gdsge/+ir/decode.m`

- [ ] **Step 2: Write `docs/deferred-features.md`**

Create `docs/deferred-features.md`:
```markdown
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
| `gdsge:codegen:sympyInterpUnsupported` | `+codegen/generateCxx.m` | SymPy backend (`UseAutoDiff=0`) with the ASG or pchip interpolant | Use the autodiff backend, or wait for Phase 8b. |
| `gdsge:codegen:unsupported` | `+codegen/generateCxx.m` | pchip interpolant, or `cxx` hook blocks, under the C++ backend | Phase 7c / 8b. |
| `gdsge:ir:incompatibleVersion` | `+ir/decode.m` | an IR JSON whose **major** version differs from `gdsge.ir.schema().irVersion` | Regenerate the IR with the matching toolbox major. See `docs/ir-versioning.md`. |

## Open deferred sub-phases

- **Phase 8b** — SymPy analytic-Jacobian backend for the ASG and pchip interpolants (unblocks
  CaoKS2016 / Bianchi2011_asg under `UseAutoDiff=0`).
- **Phase 7e-full (remainder)** — C++-body `var_tensor` (the `POPN` path) and the
  `IterRslt.var_tensor` result field. The MATLAB-side `var_tensor` subset (ndgrid tensors feeding
  inbound bounds / initial interp) landed in Phase 9a.

No corpus model exercises either, so neither blocks backward compatibility.
```
(When writing, confirm each row's identifier against the Step-1 grep output; fix any drift.)

- [ ] **Step 3: Commit**

```bash
git add docs/deferred-features.md
git commit -m "docs(phase10): deferred-features.md — single map of honest-error IDs + open sub-phases"
```

---

## Workstream 2 — IR version freeze

### Task 4: IR versioning policy doc + changelog

**Files:**
- Create: `docs/ir-versioning.md`

- [ ] **Step 1: Confirm the current version and the major-gate behavior**

Run:
```
matlab -batch "addpath('src'); disp(gdsge.ir.schema().irVersion)"
```
Expected: `1.1.0`. (Confirms the doc states the right current version.)

- [ ] **Step 2: Write `docs/ir-versioning.md`**

Create `docs/ir-versioning.md`:
```markdown
# IR versioning policy

The JSON IR (`<model>.gdsge.json`) is the stable contract between the MATLAB front-end and the code
generators: backends never re-parse `.gmod`. `gdsge.ir.schema().irVersion` carries a semantic
version, currently **1.1.0**, and `gdsge.ir.decode` enforces major-version compatibility
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

See `docs/ir-schema.md` for the full field-by-field schema (auto-generated; no-drift tested).
```

- [ ] **Step 3: Commit**

```bash
git add docs/ir-versioning.md
git commit -m "docs(phase10): ir-versioning.md — major/minor policy + 1.0.0->1.1.0 changelog"
```

---

### Task 5: IR version-freeze test (regression lock)

Lock the major-version gate so a future accidental break is caught. Uses string substitution on
encoded JSON (robust regardless of `encode`/`canonicalize` internals).

**Files:**
- Create: `tests/ir/tIrVersionFreeze.m`

- [ ] **Step 1: Write the test**

Create `tests/ir/tIrVersionFreeze.m`:
```matlab
classdef tIrVersionFreeze < matlab.unittest.TestCase
    % PHASE 10: freeze the IR major-version contract. decode accepts the
    % current major (any 1.x), rejects a different major (2.0.0). The schema
    % version starting '1.' is a tripwire: a deliberate major bump must edit
    % this test.
    methods (Test)
        function schemaIsMajor1(tc)
            v = gdsge.ir.schema().irVersion;
            tc.verifyTrue(startsWith(v, '1.'), ...
                sprintf('irVersion %s left major 1 — update the freeze test deliberately', v));
        end
        function decodeAcceptsSameMajorNewerMinor(tc)
            txt = gdsge.ir.encode(gdsgefix.minimalIR());
            txt = regexprep(txt, '"irVersion"\s*:\s*"[^"]*"', '"irVersion": "1.9.0"');
            ir = gdsge.ir.decode(txt);              % must not error
            tc.verifyTrue(isstruct(ir));
        end
        function decodeRejectsDifferentMajor(tc)
            txt = gdsge.ir.encode(gdsgefix.minimalIR());
            txt = regexprep(txt, '"irVersion"\s*:\s*"[^"]*"', '"irVersion": "2.0.0"');
            tc.verifyError(@() gdsge.ir.decode(txt), 'gdsge:ir:incompatibleVersion');
        end
    end
end
```

- [ ] **Step 2: Run it**

Run:
```
matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); addpath(fullfile(fileparts(pwd),'src','kernels')); r=runtests('ir/tIrVersionFreeze.m'); disp(table(r)); exit(double(any([r.Failed])))"
```
Expected: 3 PASS. (All three behaviors already exist in `decode.m`; this is a regression lock, so
green-on-first-run is correct.)

- [ ] **Step 3: Sanity — prove the reject test can fail**

Temporarily change `'2.0.0'` to `'1.0.0'` in `decodeRejectsDifferentMajor`, re-run that one test,
confirm it now FAILS (no error raised for a same-major version), then revert. This proves the test
has teeth.
Run:
```
matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); addpath(fullfile(fileparts(pwd),'src','kernels')); r=runtests('ir/tIrVersionFreeze/decodeRejectsDifferentMajor'); exit(double(~any([r.Failed])))"
```
Expected with the temporary edit: FAIL. Then revert the edit and re-run Step 2 (3 PASS).

- [ ] **Step 4: Commit**

```bash
git add tests/ir/tIrVersionFreeze.m
git commit -m "test(phase10): tIrVersionFreeze locks the IR major-version contract"
```

---

## Workstream 3 — Performance harness + report

### Task 6: Perf worker (one model × one backend, isolated process)

A parameterized version of the golden-capture scripts: addpath ONE source, codegen+compile, time
iter (median of reps) and simulate, write a JSON result. Invoked once per `(model, mode)` in its
own `matlab -batch` child (the MATLAB-path golden rule — old and new share flat names).

**Files:**
- Create: `tests/perf/perf_worker.m`

- [ ] **Step 1: Write the worker**

Create `tests/perf/perf_worker.m`:
```matlab
function perf_worker(model, mode, repoRoot, gmodPath, simToken, outFile)
% PERF_WORKER  Time one (model, backend) end-to-end and write a JSON result.
%   mode      : 'old' | 'autodiff' | 'sympy'
%   repoRoot  : absolute repo root (forward slashes ok)
%   gmodPath  : absolute path to the model's .gmod
%   simToken  : 'noarg' | '6x1000'  (simulate call shape)
%   outFile   : absolute path for the JSON result
% One source on the path per process. SymPy is selected by prepending
% 'UseAutoDiff=0;' to the gmod. SaveFreq=inf for all (no disk dumps -> less
% timing noise; convergence is unaffected).
REPS = 3;

oldPath = path; restorePath = onCleanup(@() path(oldPath)); %#ok<NASGU>
switch mode
    case 'old'
        addpath(fullfile(repoRoot, 'base_package', 'gdsge', 'source'));
    case {'autodiff', 'sympy'}
        addpath(fullfile(repoRoot, 'src'));
        addpath(fullfile(repoRoot, 'src', 'kernels'));
    otherwise
        error('perf_worker:badMode', 'unknown mode %s', mode);
end

work = tempname; mkdir(work);
oldCd = pwd; restoreCd = onCleanup(@() cd(oldCd)); %#ok<NASGU>

gtxt = fileread(gmodPath);
if strcmp(mode, 'sympy')
    gtxt = ['UseAutoDiff=0;' newline gtxt];
end
[~, gname, gext] = fileparts(gmodPath);
fid = fopen(fullfile(work, [gname gext]), 'w'); fwrite(fid, gtxt); fclose(fid);
cd(work);

iterFn = str2func(['iter_' model]);
simFn  = str2func(['simulate_' model]);
iterOpts = struct('SaveFreq', inf);

tc = tic; gdsge_codegen(model); tCodegen = toc(tc);

iterTimes = zeros(1, REPS); IterRslt = [];
for r = 1:REPS
    tt = tic; IterRslt = iterFn(iterOpts); iterTimes(r) = toc(tt);
end

rng(0823); ts = tic;
switch simToken
    case 'noarg';  simFn(IterRslt);
    case '6x1000'; simFn(IterRslt, struct('num_samples', 6, 'num_periods', 1000));
    otherwise; error('perf_worker:badSim', 'unknown simToken %s', simToken);
end
tSim = toc(ts);

res = struct('model', model, 'mode', mode, ...
    'tCodegen', tCodegen, 'iterTimes', iterTimes, ...
    'iterMedian', median(iterTimes), 'tSim', tSim, ...
    'Iter', double(IterRslt.Iter), 'Metric', double(IterRslt.Metric));
fid = fopen(outFile, 'w'); fwrite(fid, jsonencode(res)); fclose(fid);
fprintf('[perf_worker] %s/%s: iter median %.2fs (Iter=%d Metric=%g)\n', ...
    model, mode, res.iterMedian, res.Iter, res.Metric);
end
```

- [ ] **Step 2: Smoke-test the worker on the fastest model/backend**

Run (HL1996 / autodiff, the quickest):
```
matlab -batch "cd('tests/perf'); perf_worker('HL1996','autodiff',strrep(fileparts(fileparts(pwd)),'\','/'),strrep(fullfile(fileparts(fileparts(pwd)),'tests','HeatonLucas1996','HL1996.gmod'),'\','/'),'noarg',strrep(fullfile(tempdir,'hl_ad.json'),'\','/'))"
```
Expected: prints `[perf_worker] HL1996/autodiff: iter median ~5-18s (Iter=209 Metric=~9.6e-07)`
and writes the JSON. Confirms codegen→iter→simulate→JSON works end-to-end.

- [ ] **Step 3: Commit**

```bash
git add tests/perf/perf_worker.m
git commit -m "feat(phase10): perf_worker — isolated one-source timing of codegen/iter/simulate"
```

---

### Task 7: Perf driver + report generator

Loops the subset × backends, spawns one child per `(model, mode)`, aggregates the JSON results into
a pivoted `docs/perf-report.md`. Plus a thin `.ps1` wrapper mirroring `tests/run.ps1` (the harness
must run via `matlab -batch`, since pwsh is absent on this machine).

**Files:**
- Create: `tests/perf/run_perf.m`
- Create: `tests/perf/run_perf.ps1`

- [ ] **Step 1: Write the driver**

Create `tests/perf/run_perf.m`:
```matlab
function run_perf()
% RUN_PERF  Benchmark new (autodiff, sympy) vs old toolbox on a representative
%   subset, one source per child process; write docs/perf-report.md.
%   Usage: matlab -batch "cd('tests/perf'); run_perf"  (or run_perf.ps1)
% SLOW: ~11 child MATLAB processes, each codegen+compile+iter. Minutes.
here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(here));
matlabExe = fullfile(matlabroot, 'bin', 'matlab.exe');
outDir = tempname; mkdir(outDir);
fwd = @(p) strrep(p, '\', '/');

% {model, gmodRelPath, simToken, runSympy}
cfg = {
  'HL1996',      'tests/HeatonLucas1996/HL1996.gmod',       'noarg',  true
  'safe_assets', 'tests/Barro_et_al_2017/safe_assets.gmod', '6x1000', true
  'GLSW_interp', 'tests/GLSW2020/GLSW_interp.gmod',         '6x1000', true
  'CaoKS2016',   'tests/CaoKS2016/CaoKS2016.gmod',          '6x1000', false   % ASG: autodiff only
};

rows = {};
for i = 1:size(cfg, 1)
    model = cfg{i,1};
    gmodPath = fullfile(repoRoot, cfg{i,2});
    simToken = cfg{i,3};
    modes = {'old', 'autodiff'};
    if cfg{i,4}; modes{end+1} = 'sympy'; end %#ok<AGROW>
    for m = 1:numel(modes)
        mode = modes{m};
        outFile = fullfile(outDir, sprintf('%s_%s.json', model, mode));
        cmd = sprintf(['"%s" -batch "cd(''%s''); perf_worker(''%s'',''%s'',''%s'',''%s'',''%s'',''%s'')"'], ...
            matlabExe, fwd(here), model, mode, fwd(repoRoot), fwd(gmodPath), simToken, fwd(outFile));
        fprintf('[run_perf] %s / %s ...\n', model, mode);
        status = system(cmd);
        if status ~= 0 || exist(outFile, 'file') ~= 2
            warning('run_perf:childFailed', '%s/%s failed (status %d) — skipping', model, mode, status);
            continue;
        end
        rows{end+1} = jsondecode(fileread(outFile)); %#ok<AGROW>
    end
end

writeReport(fullfile(repoRoot, 'docs', 'perf-report.md'), rows);
fprintf('[run_perf] wrote docs/perf-report.md (%d results)\n', numel(rows));
end

function writeReport(outPath, rows)
% Pivot rows into one line per model: old / autodiff / sympy iter medians.
byKey = containers.Map('KeyType', 'char', 'ValueType', 'any');
models = {};
for i = 1:numel(rows)
    r = rows{i};
    byKey(sprintf('%s|%s', r.model, r.mode)) = r;
    if ~any(strcmp(models, r.model)); models{end+1} = r.model; end %#ok<AGROW>
end
g = @(model, mode) localGet(byKey, model, mode);

L = {};
L{end+1} = '# Phase 10 — Performance report (new vs old toolbox)';
L{end+1} = '';
L{end+1} = sprintf('Point-in-time measurement on `%s`, %d compute threads. Iter time is the median of 3 reps; SaveFreq=inf (no disk dumps). Regenerate with `tests/perf/run_perf.ps1` or `matlab -batch "cd(''tests/perf''); run_perf"`.', ...
    computer, maxNumCompThreads);
L{end+1} = '';
L{end+1} = '| Model | Iter count (old / AD / SymPy) | old iter s | AD iter s | SymPy iter s | AD/old | codegen+compile s (AD) |';
L{end+1} = '|---|---|---|---|---|---|---|';
for i = 1:numel(models)
    md = models{i};
    o = g(md,'old'); a = g(md,'autodiff'); s = g(md,'sympy');
    iters = sprintf('%s / %s / %s', numOrDash(o,'Iter'), numOrDash(a,'Iter'), numOrDash(s,'Iter'));
    ratio = '—';
    if ~isempty(o) && ~isempty(a) && o.iterMedian > 0
        ratio = sprintf('%.2f', a.iterMedian / o.iterMedian);
    end
    L{end+1} = sprintf('| %s | %s | %s | %s | %s | %s | %s |', md, iters, ...
        secOrDash(o), secOrDash(a), secOrDash(s), ratio, secCodegen(a)); %#ok<AGROW>
end
L{end+1} = '';
L{end+1} = '**Reading:** the autodiff backend emits the same adept/CoDoSol/kernel structure as the old toolbox, so AD/old ≈ 1 is the no-regression bar. The SymPy column quantifies the analytic-Jacobian backend vs autodiff. Iter counts must match across columns (same convergence).';
L{end+1} = '';

fid = fopen(outPath, 'w'); fwrite(fid, strjoin(L, newline)); fclose(fid);
end

function r = localGet(byKey, model, mode)
k = sprintf('%s|%s', model, mode);
if isKey(byKey, k); r = byKey(k); else; r = []; end
end
function s = numOrDash(r, f); if isempty(r); s = '—'; else; s = sprintf('%d', r.(f)); end; end
function s = secOrDash(r);  if isempty(r); s = '—'; else; s = sprintf('%.2f', r.iterMedian); end; end
function s = secCodegen(r); if isempty(r); s = '—'; else; s = sprintf('%.1f', r.tCodegen); end; end
```

- [ ] **Step 2: Write the `.ps1` wrapper**

Create `tests/perf/run_perf.ps1` (mirrors `tests/run.ps1`):
```powershell
# Run the performance harness headless. Spawns one MATLAB child per (model,backend).
$ErrorActionPreference = 'Stop'
$matlab = 'C:\Program Files\MATLAB\R2025b\bin\matlab.exe'
if (-not (Test-Path $matlab)) { throw "MATLAB not found at $matlab" }
$perfDir = $PSScriptRoot
& $matlab -batch "cd('$perfDir'); run_perf"
exit $LASTEXITCODE
```

- [ ] **Step 3: Commit (harness only; report generated in Task 8)**

```bash
git add tests/perf/run_perf.m tests/perf/run_perf.ps1
git commit -m "feat(phase10): run_perf driver + ps1 wrapper (subset x {old,AD,SymPy}, pivoted report)"
```

---

### Task 8: Run the harness, commit the report

**Files:**
- Create: `docs/perf-report.md` (generated)

- [ ] **Step 1: Run the full harness**

Run:
```
matlab -batch "cd('tests/perf'); run_perf"
```
Expected: per-child `[run_perf] <model> / <mode> ...` lines, each child printing its iter median +
`Iter`/`Metric`; finally `wrote docs/perf-report.md (N results)`. Cross-check the printed `Iter`
values against the known goldens (HL1996≈209, safe_assets≈1271, GLSW≈903, CaoKS2016≈281) — and
confirm old/AD/SymPy agree per model. Slow (minutes); if a child fails, the driver warns and
continues, and the report shows `—` for that cell.

- [ ] **Step 2: Eyeball the report**

Run:
```
cat docs/perf-report.md
```
Expected: a table with AD/old ≈ 1 (no regression) and a populated SymPy column for the three spline
models (`—` for CaoKS2016 SymPy, which is ASG). Iter counts match across columns.

- [ ] **Step 3: Commit**

```bash
git add docs/perf-report.md
git commit -m "docs(phase10): perf-report — AD parity vs old + SymPy-vs-AD on the subset"
```

---

## Workstream 4 — Docs

### Task 9: Root README

**Files:**
- Create: `README.md`

- [ ] **Step 1: Write `README.md`**

Create `README.md` with these sections (real content, drawn from the verified facts below):

1. **Title + one-paragraph what-this-is** — "A ground-up refactor of the GDSGE MATLAB toolbox
   (solves global DSGE models) around an explicit JSON IR. Drop-in backward compatible: every
   existing `.gmod` runs unchanged and `IterRslt`/`SimuRslt` shapes are frozen. Two C++ Jacobian
   backends: adept autodiff (default, no Python) and a SymPy analytic Jacobian."
2. **Quickstart** — from a directory containing `<model>.gmod`, with `src/` on the path:
   ```matlab
   addpath('src'); addpath('src/kernels');
   gdsge_codegen('HL1996');          % parse -> IR -> generate iter_/simulate_/mex_ + compile
   IterRslt = iter_HL1996;           % solve
   SimuRslt = simulate_HL1996(IterRslt);
   ```
   Or the one-call orchestrator: `eq = gdsge('HL1996');` (returns `model`/`IterRslt`/`SmltRslt`).
   Point at `tests/HeatonLucas1996/HL1996.gmod` as a worked example.
3. **The two backends** — default autodiff needs only MATLAB + a C++ MEX compiler; the SymPy
   analytic-Jacobian backend is selected by `UseAutoDiff=0;` in the gmod and needs the `uv` Python
   env under `pyext/`. Link `docs/user-guide.md` for when to use which. Note ASG/pchip are
   autodiff-only today (link `docs/deferred-features.md`).
4. **Running the tests** — `matlab -batch "cd('tests'); run_tests"` (exit 0 = all pass; report at
   `tests/results/junit.xml`). Note: `tests/run.ps1` needs PowerShell 7 (`pwsh`), absent on this
   machine — use `matlab -batch` directly. SymPy tests: `uv run --project pyext pytest pyext/tests`.
5. **Repo layout** — table: `src/+gdsge/{+parser,+ir,+codegen,+runtime}` (the pipeline),
   `src/{gdsge.m,gdsge_codegen.m,asg.m,kernels/}` (flat public shims + vendored kernels),
   `templates/cxx/` (C++ templates), `pyext/` (SymPy backend), `tests/` (per-model gates + suites),
   `docs/`, `base_package/` (old toolbox, git-ignored reference only).
6. **Architecture in one line** — `model.gmod → MATLAB parser → JSON IR → { MATLAB codegen, C++
   codegen (adept | SymPy) }`; the IR is the contract.
7. **Links** — `docs/superpowers/specs/2026-06-11-refactor-gdsge-design.md` (design),
   `docs/architecture.md`, `docs/user-guide.md`, `docs/ir-schema.md`, `docs/ir-versioning.md`,
   `docs/perf-report.md`, `docs/deferred-features.md`, `PROGRESS.md`.

- [ ] **Step 2: Sanity-check internal links resolve**

Run:
```
grep -oE "\]\(([^)]+\.md[^)]*)\)" README.md | sed -E 's/.*\(([^)]+)\)/\1/' | while read p; do f="${p%%#*}"; if [ -f "$f" ]; then echo "OK   $f"; else echo "MISS $f"; fi; done
```
Expected: every listed `.md` prints `OK` (all targets exist by the time docs are written last).

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "docs(phase10): root README — what/quickstart/backends/tests/layout/links"
```

---

### Task 10: Architecture / developer guide

**Files:**
- Create: `docs/architecture.md`

- [ ] **Step 1: Write `docs/architecture.md`**

Create `docs/architecture.md` consolidating the pipeline knowledge (real content from these facts):

1. **The pipeline** — `gmod → gdsge.parser → JSON IR → { gdsge.codegen.mat, gdsge.codegen.generateCxx
   (adept | sympy) }`. One diagram/ordered list. State that the IR is the contract — backends never
   re-parse gmod.
2. **Module map** — `gdsge.parser.*` (preprocess, expandMacros, splitBlocks, splitStatements,
   tokenize, parseExpr, parseModel, analyzeModel, parseFrontEnd), `gdsge.ir.*` (schema → validate,
   canonicalize, encode/decode/roundtrip, gendoc), `gdsge.codegen.*` (codegen driver,
   generateMatlab, generateCxx, dataLayout, the `+cxx` and `+cxx/+sympymodel` emitters),
   `gdsge.runtime.*` (unpackOptions, solveProblems / solveProblemsAsg, constructSplines,
   computeMetric, progress/reporting).
3. **How each backend consumes the IR** — the shared `gdsge.codegen.dataLayout` contract (MATLAB
   packer + C++ `POPN`/`POPNARRAY` macros from one descriptor); the autodiff `emitExpr` AST→C++
   printer; the SymPy registry + `pyenv` bridge (`gdsge.codegen.sympy.{ensurePyenv,callSympy}` →
   `pyext/gdsge_sympy`), with the §10 reduction-fusion chain rule closing through the spline kernel.
   Note the region-agnostic per-body emitters (Phase 9b): only the drivers loop `model.regions`.
4. **How to add a model** — the model-body probe-first discipline (Phases 9a/9b): write a probe that
   enumerates every expression/operator form the gmod uses, check each against what the parser/
   backends already cover, write the gap list, THEN spec → plan → implement, both backends to
   parity, vs a captured golden. Reference `scratch/probe_*.m` as examples.
5. **The MATLAB-path golden rule** — old and new share flat public names; each run addpaths exactly
   ONE source in its own `matlab -batch` process with a controlled cwd, restoring the path on
   cleanup. Why (shadowing). Cite the golden-capture scripts and `tests/perf/perf_worker.m` as the
   canonical pattern.
6. **Testing model** — native `matlab.unittest`, headless via `matlab -batch`; per-model gates
   (`tGolden*`, `tFrontEnd*`, `tEndToEnd*[Sympy]`); goldens captured from the old toolbox; tolerance
   (not bit-exact) comparison via `gdsgetest.compareNumericClose`.
7. **Pointers** — `docs/ir-schema.md`, `docs/ir-versioning.md`, `docs/deferred-features.md`,
   `docs/notes-for-agents.md`.

- [ ] **Step 2: Verify the module names referenced actually exist**

Run:
```
ls src/+gdsge/+parser/ src/+gdsge/+ir/ src/+gdsge/+codegen/ src/+gdsge/+runtime/
```
Cross-check that every module named in the guide appears (fix any name that drifted).

- [ ] **Step 3: Commit**

```bash
git add docs/architecture.md
git commit -m "docs(phase10): architecture guide — pipeline, module map, add-a-model, golden rule"
```

---

### Task 11: User guide

**Files:**
- Create: `docs/user-guide.md`

- [ ] **Step 1: Write `docs/user-guide.md`**

Create `docs/user-guide.md` (self-contained; complements gdsge.com). Source the inventory from
`docs/notes-for-agents.md` §"Full gmod feature inventory" and the parent spec §13:

1. **Authoring a model** — the gmod blocks (`parameters`, `var_state`, `var_shock`,
   `var_policy[/_init]`, `var_aux[/_init]`, `var_interp`, `var_tensor`, `var_output`, `var_others`,
   `model[/_init]`, `equations`, `simulate`, `pre_model`, `pre_iter`, `post_iter`,
   `pre_jac_code`, `post_jac_code`); declarations (`inbound[/_init]` + `adaptive(factor)`,
   `initial`, `shock_num`, `shock_trans`); operators (`GDSGE_EXPECT/MIN/MAX/PROD`,
   `GDSGE_INTERP_VEC` and the primed `'` form); markers (trailing `'` = next-state/future, `[N]` =
   shock-indexed array var); built-ins (`GDSGE_Iter`, `TASK`, `OUTPUT_CONSTRUCT_CODE`); the Phase 9b
   conditional `model(<cond>)` regions + `if/else/end` conditional equations; macros
   (`#define`, `#for`, `#foreach`, `#if`, `#mat{}`, `#strcat_comma`, `include`, `cinclude`).
   Note the gmod is a MATLAB superset (grid/param lines are real MATLAB, `eval`'d during parsing).
2. **Options reference** — interpolation (`USE_SPLINE`/`USE_ASG`/`USE_PCHIP`, `INTERP_ORDER`,
   `ExtrapOrder`), simulate (`SIMU_RESOLVE`/`SIMU_INTERP`, `num_samples`, `num_periods`, `SimuPrintFreq`,
   `SimuSaveFreq`), solver/iteration (`TolEq`, `MaxIter`, `SaveFreq`, `PrintFreq`, `NumThreads`),
   ASG (`AsgMinLevel`, `AsgMaxLevel`, `AsgThreshold`, `AsgOutputMaxLevel`, `AsgOutputThreshold`),
   and `UseAutoDiff` (1 = adept autodiff default; 0 = SymPy analytic Jacobian — spline interpolants
   only today). One row each: name, meaning, default.
3. **Selecting the SymPy backend** — add `UseAutoDiff=0;` to the gmod; requires the `uv` Python env
   (`uv sync --project pyext`); spline-only (ASG/pchip → clear error, see deferred-features). Same
   results within tolerance as autodiff (proven on 6 models).
4. **Reading results** — `IterRslt` (`Iter`, `Metric`, `var_state`, `var_policy`, `var_aux`,
   `var_interp`, `shock_trans`, …) and `SimuRslt` (`shock`, per-variable paths) field shapes — frozen
   and backward compatible.
5. **Deferred features** — short pointer to `docs/deferred-features.md`.

- [ ] **Step 2: Sanity-check option names against the codebase**

Run:
```
grep -roiE "AsgMinLevel|AsgMaxLevel|AsgThreshold|SIMU_INTERP|SIMU_RESOLVE|UseAutoDiff|TolEq|PrintFreq|SaveFreq|ExtrapOrder|INTERP_ORDER" src/ | sort -u | head -40
```
Cross-check that every option named in the guide appears in the runtime/parser/codegen (fix drift).

- [ ] **Step 3: Commit**

```bash
git add docs/user-guide.md
git commit -m "docs(phase10): user guide — gmod authoring, options, SymPy backend, reading results"
```

---

## Workstream 5 — Finalize

### Task 12: Update PROGRESS, run the full suite, merge

**Files:**
- Modify: `PROGRESS.md`

- [ ] **Step 1: Run the FULL suite (the real gate for the phase)**

Run:
```
matlab -batch "cd('tests'); run_tests"
```
Expected: exit 0; `tests/results/junit.xml` shows all tests pass (≥ the prior 447, plus the 3 new
`tIrVersionFreeze` tests → ~450). The trimmed Mendoza golden and new freeze test are exercised here.
If any fail, STOP and fix before proceeding (do not mark the phase done with a red suite).

- [ ] **Step 2: Update `PROGRESS.md`**

Mark Phase 10 ☑ (was ☐ at line ~189), and add a Changelog entry at the top of the Changelog section
dated 2026-06-14 summarizing: Mendoza golden trimmed 97MB→~9MB; `docs/deferred-features.md`,
`docs/ir-versioning.md` + `tIrVersionFreeze`; perf harness `tests/perf/` + `docs/perf-report.md`
(AD≈old, SymPy quantified); `README.md` + `docs/architecture.md` + `docs/user-guide.md`. State the
only remaining open items are the deferred sub-phases **8b** (SymPy ASG+pchip) and **7e-full
remainder** (C++-body `var_tensor`). Note the full suite is green.

- [ ] **Step 3: Commit**

```bash
git add PROGRESS.md
git commit -m "docs(phase10): PROGRESS — Phase 10 complete (docs, IR freeze, perf, cleanup)"
```

- [ ] **Step 4: Finish the branch**

Use the `superpowers:finishing-a-development-branch` skill to choose merge/PR/cleanup for
`phase10-polish` (no GitHub remote exists, so the likely path is a local merge to `main`). Do NOT
merge before Step 1 is green.

---

## Self-review notes (author)

- **Spec coverage:** WS1 §4 → Tasks 1–3; WS2 §5 → Tasks 4–5; WS3 §6 → Tasks 6–8; WS4 §7 → Tasks
  9–11; finalize/PROGRESS → Task 12. The spec's "scratch left untouched" and "no perf gate in the
  suite" are honored (no task touches scratch; Task 8 commits a report, adds no timing assertion).
- **Concrete IDs:** the deferred-error identifiers in Task 3 are the grep-verified current IDs (the
  spec draft's `varTensorUnsupported`/`codegen:varTensorUnsupported` rows were stale and are
  replaced by `varTensorInBodyUnsupported`/`varTensorAsgUnsupported`).
- **Perf subset:** HL1996 / safe_assets / GLSW_interp / CaoKS2016, with simulate shapes and Iter
  expectations matching the golden-capture scripts; CaoKS2016 SymPy intentionally absent (ASG).
- **Naming consistency:** worker `perf_worker(model,mode,repoRoot,gmodPath,simToken,outFile)` is
  called with exactly those 6 args by `run_perf.m`; modes `old|autodiff|sympy`; simTokens
  `noarg|6x1000` match between driver config and worker switch.
