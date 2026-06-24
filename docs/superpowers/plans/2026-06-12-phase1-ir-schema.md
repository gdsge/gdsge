# Phase 1 — IR Schema + MATLAB Scaffolding — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the `gdsge.ir` package — a single declarative schema descriptor that drives validation, MATLAB↔JSON round-trip, and doc generation — plus the AST node API and a hand-authored HL1996 reference IR that round-trips, validates, and serves as Phase 3's golden target.

**Architecture:** One descriptor (`gdsge.ir.schema`) is the source of truth. A small set of *field kinds* (`scalar`, `matrix`, `text`, `enum`, `ref`, `list`, `map`, `tagged`, `taggedList`, `struct`) is the vocabulary a generic walker understands. `validate`, `canonicalize`, `encode` (via `toJsonReady`), and `gendoc` each walk the descriptor. AST nodes and model statements are *tagged unions* (dispatch on a `kind`/`type` field against the `nodes`/`stmts` registries). Round-trip quirks (orientation, array-kind, empty/scalar, non-finite) are defeated in `canonicalize`/`toJsonReady`, one rule per kind.

**Tech Stack:** MATLAB R2025b (`jsonencode`/`jsondecode`, `matlab.unittest`), the Phase-0 headless harness (`tests/run_tests.m`, `tests/run.ps1`).

**Spec:** `docs/superpowers/specs/2026-06-12-phase1-ir-schema-design.md`.

---

## File structure

Created in `src/+gdsge/+ir/`:

- `schema.m` — THE descriptor: top-level sections, `nodes` registry, `stmts` registry, `irVersion`. Field-spec helpers are local functions.
- `validate.m` — `validate(ir) -> {pass, errors}`; descriptor-driven shape/kind/enum checks + node/stmt dispatch + ref pools + slot + option-invariant checks.
- `canonicalize.m` — `canonicalize(ir) -> ir`; normalizes any IR (in-memory or freshly decoded) to canonical MATLAB form. Idempotent.
- `decode.m` — `decode(txt) -> ir` = `canonicalize(jsondecode(txt))`.
- `encode.m` — `encode(ir) -> txt`; `toJsonReady` (local) wraps matrices and tags non-finite, then `jsonencode(...,'PrettyPrint',true)`.
- `roundtrip.m` — `roundtrip(ir) -> decode(encode(ir))`.
- `isequalIR.m` — `isequalIR(a,b) -> logical`; `deepEqual` (local) over canonicalized inputs.
- `gendoc.m` — `gendoc()` returns markdown; `gendoc(path)` also writes it.
- `+node/{num,name,primed,unop,binop,call,index,kindOf,children}.m` — AST constructors + accessors.
- `Contents.m` — package help.

Modified:

- `tests/run_tests.m` — add `src/` to the path so `gdsge.ir.*` resolves.

Tests (plain folders so `TestSuite.fromFolder` finds the classes):

- `tests/ir/tNodeConstructors.m`, `tests/ir/tSchema.m`, `tests/ir/tValidateShape.m`, `tests/ir/tValidateSemantics.m`, `tests/ir/tRoundtrip.m`, `tests/ir/tGendoc.m`
- `tests/HeatonLucas1996/ir/buildHL1996IR.m` (fixture source), `tests/HeatonLucas1996/ir/HL1996.gdsge.json` (generated golden), `tests/HeatonLucas1996/ir/tIrHL1996.m`

Generated artifacts committed: `docs/ir-schema.md` (regenerated from `schema.m`), `tests/HeatonLucas1996/ir/HL1996.gdsge.json`.

**Descriptor vocabulary note (refines spec §5.1):** the design listed `scalar/matrix/enum/text/ref/list/nodeList/struct/optional`. Implementation uses two generalizations: `optional` becomes a `required:false` flag on any spec; `node`/`nodeList` become the generic `tagged`/`taggedList` kinds (so the *same* machinery serves AST nodes *and* model statements), plus `map` for dynamic-key structs (`shocks.transitions`, `states.grids`). These are mechanical refinements, not scope changes.

---

## Task 1: AST node API + test-harness path

**Files:**
- Modify: `tests/run_tests.m`
- Create: `src/+gdsge/+ir/+node/num.m`, `name.m`, `primed.m`, `unop.m`, `binop.m`, `call.m`, `index.m`, `kindOf.m`, `children.m`
- Test: `tests/ir/tNodeConstructors.m`

- [ ] **Step 1: Write the failing test**

`tests/ir/tNodeConstructors.m`:
```matlab
classdef tNodeConstructors < matlab.unittest.TestCase
    methods (Test)
        function leafNodes(tc)
            import gdsge.ir.node.*
            tc.verifyEqual(num(5),    struct('kind','num','value',5));
            tc.verifyEqual(name('c1'),struct('kind','name','id','c1'));
            tc.verifyEqual(primed('g'),struct('kind','primed','id','g'));
        end
        function unaryAndBinary(tc)
            import gdsge.ir.node.*
            u = unop('-', num(1));
            tc.verifyEqual(u.kind, 'unop');
            tc.verifyEqual(u.op, '-');
            tc.verifyEqual(u.arg, num(1));
            b = binop('+', name('a'), num(2));
            tc.verifyEqual(b.kind, 'binop');
            tc.verifyEqual(b.op, '+');
            tc.verifyEqual(b.lhs, name('a'));
            tc.verifyEqual(b.rhs, num(2));
        end
        function callAndIndexHoldCellArgs(tc)
            import gdsge.ir.node.*
            c = call('exp', {num(1)});
            tc.verifyEqual(c.kind, 'call');
            tc.verifyEqual(c.fn, 'exp');
            tc.verifyTrue(iscell(c.args));
            tc.verifyEqual(c.args{1}, num(1));
            ix = index(name('w1n'), {num(3)});
            tc.verifyEqual(ix.kind, 'index');
            tc.verifyEqual(ix.base, name('w1n'));
            tc.verifyTrue(iscell(ix.args));
        end
        function accessors(tc)
            import gdsge.ir.node.*
            b = binop('*', name('a'), name('b'));
            tc.verifyEqual(gdsge.ir.node.kindOf(b), 'binop');
            kids = gdsge.ir.node.children(b);
            tc.verifyEqual(numel(kids), 2);
            tc.verifyEqual(kids{1}, name('a'));
            tc.verifyEqual(gdsge.ir.node.children(num(7)), {});
            c = call('max', {num(1), num(2)});
            tc.verifyEqual(numel(gdsge.ir.node.children(c)), 2);
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `pwsh -File tests/run.ps1`
Expected: FAIL — `gdsge.ir.node.num` is undefined (and `src/` not yet on path).

- [ ] **Step 3: Put `src/` on the test path**

In `tests/run_tests.m`, after line 11 (`addpath(thisDir); ...`), add the repo `src/` directory:
```matlab
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);   % so the +gdsgetest package (utility fns) resolves inside tests
addpath(fullfile(fileparts(thisDir), 'src'));   % new-pipeline source: resolves gdsge.ir.*
```
(Path policy: the *one* toolbox source a new-pipeline test process adds is `src/`. Golden capture from the old toolbox runs in its own separate process — untouched.)

- [ ] **Step 4: Implement the node constructors**

`src/+gdsge/+ir/+node/num.m`:
```matlab
function n = num(value)
% NUM  Numeric-literal AST node.
n.kind = 'num';
n.value = value;
end
```

`src/+gdsge/+ir/+node/name.m`:
```matlab
function n = name(id)
% NAME  Current-state identifier AST node.
n.kind = 'name';
n.id = id;
end
```

`src/+gdsge/+ir/+node/primed.m`:
```matlab
function n = primed(id)
% PRIMED  Future/next-state (trailing-quote) identifier AST node.
n.kind = 'primed';
n.id = id;
end
```

`src/+gdsge/+ir/+node/unop.m`:
```matlab
function n = unop(op, arg)
% UNOP  Unary-operator AST node.
n.kind = 'unop';
n.op = op;
n.arg = arg;
end
```

`src/+gdsge/+ir/+node/binop.m`:
```matlab
function n = binop(op, lhs, rhs)
% BINOP  Binary-operator AST node.
n.kind = 'binop';
n.op = op;
n.lhs = lhs;
n.rhs = rhs;
end
```

`src/+gdsge/+ir/+node/call.m`:
```matlab
function n = call(fn, args)
% CALL  Function-call AST node. ARGS is a cell array of nodes.
n.kind = 'call';
n.fn = fn;
n.args = args;
end
```

`src/+gdsge/+ir/+node/index.m`:
```matlab
function n = index(base, args)
% INDEX  Indexing AST node. ARGS is a cell array of nodes.
n.kind = 'index';
n.base = base;
n.args = args;
end
```

`src/+gdsge/+ir/+node/kindOf.m`:
```matlab
function k = kindOf(node)
% KINDOF  Return the kind tag of an AST node.
k = node.kind;
end
```

`src/+gdsge/+ir/+node/children.m`:
```matlab
function kids = children(node)
% CHILDREN  Return a cell array of an AST node's child nodes (empty for leaves).
switch node.kind
    case {'num','name','primed'}
        kids = {};
    case 'unop'
        kids = {node.arg};
    case 'binop'
        kids = {node.lhs, node.rhs};
    case 'call'
        kids = node.args;
    case 'index'
        kids = [{node.base}, reshape(node.args, 1, [])];
    otherwise
        error('gdsge:ir:node:unknownKind', 'Unknown node kind "%s".', node.kind);
end
end
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `pwsh -File tests/run.ps1`
Expected: PASS — `tNodeConstructors` (4 tests) green; Phase-0 tests still green.

- [ ] **Step 6: Commit**

```bash
git add tests/run_tests.m src/+gdsge/+ir/+node tests/ir/tNodeConstructors.m
git commit -m "feat(ir): AST node constructors + accessors; wire src onto test path"
```

---

## Task 2: Schema descriptor

**Files:**
- Create: `src/+gdsge/+ir/schema.m`
- Test: `tests/ir/tSchema.m`

- [ ] **Step 1: Write the failing test**

`tests/ir/tSchema.m`:
```matlab
classdef tSchema < matlab.unittest.TestCase
    methods (Test)
        function versionIsSemver(tc)
            s = gdsge.ir.schema();
            tc.verifyTrue(ischar(s.irVersion));
            tc.verifyNotEmpty(regexp(s.irVersion, '^\d+\.\d+\.\d+$', 'once'));
        end
        function rootIsStructWithAllSections(tc)
            s = gdsge.ir.schema();
            tc.verifyEqual(s.root.kind, 'struct');
            expected = {'irVersion','modelName','options','shocks','states', ...
                        'variables','bounds','interp','model','modelInit', ...
                        'simulate','hooks'};
            tc.verifyEqual(sort(fieldnames(s.root.fields))', sort(expected));
        end
        function nodeRegistryHasSevenKinds(tc)
            s = gdsge.ir.schema();
            tc.verifyEqual(sort(fieldnames(s.nodes))', ...
                sort({'num','name','primed','unop','binop','call','index'}));
        end
        function stmtRegistryHasThreeTypes(tc)
            s = gdsge.ir.schema();
            tc.verifyEqual(sort(fieldnames(s.stmts))', ...
                sort({'assign','interpCall','reduction'}));
        end
        function everyFieldHasAKnownKind(tc)
            s = gdsge.ir.schema();
            known = {'scalar','matrix','text','enum','ref','list','map', ...
                     'tagged','taggedList','struct'};
            walk(tc, s.root, known);
            cellfun(@(k) walk(tc, s.nodes.(k), known), fieldnames(s.nodes));
            cellfun(@(k) walk(tc, s.stmts.(k), known), fieldnames(s.stmts));
        end
    end
end

function walk(tc, spec, known)
    tc.verifyTrue(ismember(spec.kind, known), ...
        sprintf('unknown kind "%s"', spec.kind));
    switch spec.kind
        case 'struct'
            fn = fieldnames(spec.fields);
            for i = 1:numel(fn); walk(tc, spec.fields.(fn{i}), known); end
        case 'list'
            walk(tc, spec.item, known);
        case 'map'
            walk(tc, spec.value, known);
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('ir','tSchema.m')); disp(table(r))"`
Expected: FAIL — `gdsge.ir.schema` undefined.

- [ ] **Step 3: Implement the descriptor**

`src/+gdsge/+ir/schema.m`:
```matlab
function s = schema()
% SCHEMA  The single source of truth for the GDSGE IR.
%   Returns a descriptor struct with fields:
%     .irVersion  semver string
%     .root       a 'struct' field-spec whose .fields are the top-level sections
%     .nodes      registry of AST node specs, keyed by node kind
%     .stmts      registry of model-statement specs, keyed by statement type
%   Field kinds: scalar matrix text enum ref list map tagged taggedList struct.

s.irVersion = '1.0.0';

% ---- reusable item specs -------------------------------------------------
varItem = fStruct(structOf('name', fText(), 'length', fScalar(), 'slot', fMatrix()));

boundItem = fStruct(structOf( ...
    'name', fRef('policyAux'), ...
    'lower', fText(), ...
    'upper', fText(), ...
    'adaptiveFactor', opt(fScalar())));

interpItem = fStruct(structOf( ...
    'name', fText(), ...
    'args', fList(fRef('states')), ...
    'initialExpr', fText(), ...
    'updateExpr', fText()));

eqItem = fStruct(structOf('expr', fNode(), 'primed', fScalar()));

initItem  = fStruct(structOf('var', fText(), 'value', fText()));
transItem = fStruct(structOf('state', fRef('states'), 'expr', fText()));

% ---- top-level sections --------------------------------------------------
options = fStruct(structOf( ...
    'interpMethod', fEnum({'spline','asg','pchip'}), ...
    'interpOrder',  opt(fScalar()), ...
    'extrapOrder',  opt(fScalar()), ...
    'asgMaxLevel',  opt(fScalar()), ...
    'asgThreshold', opt(fScalar()), ...
    'tolEq',        opt(fScalar()), ...
    'numThreads',   opt(fScalar()), ...
    'simuResolve',  opt(fScalar()), ...
    'simuInterp',   opt(fScalar()), ...
    'printFreq',    opt(fScalar()), ...
    'saveFreq',     opt(fScalar()), ...
    'simuPrintFreq',opt(fScalar()), ...
    'simuSaveFreq', opt(fScalar())));

shocks = fStruct(structOf( ...
    'names', fList(fText()), ...
    'count', fScalar(), ...
    'values', fMap(fMatrix()), ...
    'transitions', fMap(fMatrix())));

states = fStruct(structOf( ...
    'names', fList(fText()), ...
    'grids', fMap(fText())));

variables = fStruct(structOf( ...
    'policy', fList(varItem), ...
    'aux',    fList(varItem), ...
    'interp', fList(fText()), ...
    'tensor', fList(fText()), ...
    'output', fList(fRef('policyAux')), ...
    'others', fList(fText())));

model = fStruct(structOf( ...
    'statements', fStmtList(), ...
    'equations',  fList(eqItem)));

modelInit = fStruct(structOf( ...
    'variables', fStruct(structOf('policyInit', fList(varItem), 'auxInit', fList(varItem))), ...
    'bounds',    fList(boundItem), ...
    'statements', fStmtList(), ...
    'equations',  fList(eqItem)));

simulate = fStruct(structOf( ...
    'numPeriods', fScalar(), ...
    'numSamples', fScalar(), ...
    'initial',    fList(initItem), ...
    'varSimu',    fList(fRef('policyAux')), ...
    'transitions', fList(transItem)));

hooks = fStruct(structOf( ...
    'preModel',    opt(fText()), ...
    'preIter',     opt(fText()), ...
    'postIter',    opt(fText()), ...
    'preJacCode',  opt(fText()), ...
    'postJacCode', opt(fText()), ...
    'cxx',         opt(fText())));

s.root = fStruct(structOf( ...
    'irVersion', fText(), ...
    'modelName', fText(), ...
    'options',   options, ...
    'shocks',    shocks, ...
    'states',    states, ...
    'variables', variables, ...
    'bounds',    fList(boundItem), ...
    'interp',    fList(interpItem), ...
    'model',     model, ...
    'modelInit', opt(modelInit), ...
    'simulate',  simulate, ...
    'hooks',     hooks));

% ---- AST node registry ---------------------------------------------------
s.nodes = structOf( ...
    'num',    fStruct(structOf('kind', fEnum({'num'}),    'value', fScalar())), ...
    'name',   fStruct(structOf('kind', fEnum({'name'}),   'id',    fText())), ...
    'primed', fStruct(structOf('kind', fEnum({'primed'}), 'id',    fText())), ...
    'unop',   fStruct(structOf('kind', fEnum({'unop'}),   'op', fText(), 'arg', fNode())), ...
    'binop',  fStruct(structOf('kind', fEnum({'binop'}),  'op', fText(), 'lhs', fNode(), 'rhs', fNode())), ...
    'call',   fStruct(structOf('kind', fEnum({'call'}),   'fn', fText(), 'args', fNodeList())), ...
    'index',  fStruct(structOf('kind', fEnum({'index'}),  'base', fNode(), 'args', fNodeList())));

% ---- model statement registry --------------------------------------------
s.stmts = structOf( ...
    'assign', fStruct(structOf( ...
        'type', fEnum({'assign'}), 'target', fText(), 'primed', fScalar(), 'expr', fNode())), ...
    'interpCall', fStruct(structOf( ...
        'type', fEnum({'interpCall'}), 'targets', fList(fText()), 'primed', fScalar(), ...
        'args', fNodeList(), 'interpRef', fText())), ...
    'reduction', fStruct(structOf( ...
        'type', fEnum({'reduction'}), 'kind', fEnum({'EXPECT','MIN','MAX','PROD'}), ...
        'target', fText(), 'body', fNode(), 'transRef', fRef('transitions'))));
end

% ===== field-spec constructors (local) ====================================
function s = fScalar(),  s.kind = 'scalar';  s.required = true; end
function s = fMatrix(),  s.kind = 'matrix';  s.required = true; end
function s = fText(),    s.kind = 'text';    s.required = true; end
function s = fEnum(vals),s.kind = 'enum';    s.required = true; s.values = vals; end
function s = fRef(pool), s.kind = 'ref';     s.required = true; s.pool = pool; end
function s = fNode(),    s.kind = 'tagged';     s.required = true; s.registry = 'nodes'; s.tag = 'kind'; end
function s = fNodeList(),s.kind = 'taggedList'; s.required = true; s.registry = 'nodes'; s.tag = 'kind'; end
function s = fStmtList(),s.kind = 'taggedList'; s.required = true; s.registry = 'stmts'; s.tag = 'type'; end
function s = fList(item),s.kind = 'list';    s.required = true; s.item = item; end
function s = fMap(val),  s.kind = 'map';     s.required = true; s.value = val; end
function s = fStruct(flds), s.kind = 'struct'; s.required = true; s.fields = flds; end
function s = opt(spec),  spec.required = false; s = spec; end

% structOf('a',specA,'b',specB,...) builds a struct mapping names to specs,
% preserving insertion order (used for deterministic field order in JSON/doc).
function st = structOf(varargin)
st = struct();
for i = 1:2:numel(varargin)
    st.(varargin{i}) = varargin{i+1};
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('ir','tSchema.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 5 tests.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+ir/schema.m tests/ir/tSchema.m
git commit -m "feat(ir): declarative schema descriptor (source of truth)"
```

---

## Task 3: Validator — structural shape, kinds, tagged dispatch

**Files:**
- Create: `src/+gdsge/+ir/validate.m`
- Test: `tests/ir/tValidateShape.m`

- [ ] **Step 1: Write the failing test**

`tests/ir/tValidateShape.m`:
```matlab
classdef tValidateShape < matlab.unittest.TestCase
    methods (Test)
        function minimalValidIrPasses(tc)
            ir = gdsgefix.minimalIR();
            r = gdsge.ir.validate(ir);
            tc.verifyTrue(r.pass, strjoin(r.errors, ' | '));
        end
        function missingRequiredFieldFails(tc)
            ir = gdsgefix.minimalIR();
            ir = rmfield(ir, 'modelName');
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
            tc.verifyTrue(any(contains(r.errors, 'modelName')));
        end
        function wrongKindForScalarFails(tc)
            ir = gdsgefix.minimalIR();
            ir.shocks.count = 'eight';   % should be scalar
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
            tc.verifyTrue(any(contains(r.errors, 'shocks.count')));
        end
        function badEnumFails(tc)
            ir = gdsgefix.minimalIR();
            ir.options.interpMethod = 'banana';
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
            tc.verifyTrue(any(contains(r.errors, 'interpMethod')));
        end
        function unknownNodeKindFails(tc)
            ir = gdsgefix.minimalIR();
            ir.model.equations{1}.expr = struct('kind','bogus');
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
            tc.verifyTrue(any(contains(r.errors, 'bogus')));
        end
        function malformedBinopFails(tc)
            ir = gdsgefix.minimalIR();
            bad = struct('kind','binop','op','+','lhs',gdsge.ir.node.num(1)); % no rhs
            ir.model.equations{1}.expr = bad;
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
            tc.verifyTrue(any(contains(r.errors, 'rhs')));
        end
    end
end
```

This test depends on a tiny shared fixture builder `gdsgefix.minimalIR`. Create it now in a package folder under `tests/` so every IR test can reuse it.

`tests/+gdsgefix/minimalIR.m`:
```matlab
function ir = minimalIR()
% A smallest valid IR: 1 shock, 1 state, 1 policy, 1 equation. Used by IR unit tests.
import gdsge.ir.node.*
ir.irVersion = '1.0.0';
ir.modelName = 'mini';
ir.options = struct('interpMethod','spline','interpOrder',4);
ir.shocks = struct('names', {{'z'}}, 'count', 1, ...
    'values', struct('z', 1.0), ...
    'transitions', struct('shock_trans', 1.0));
ir.states = struct('names', {{'k'}}, 'grids', struct('k','linspace(0,1,3)'));
ir.variables = struct( ...
    'policy', {{ struct('name','c','length',1,'slot',[1 1]) }}, ...
    'aux',    {{}}, ...
    'interp', {{}}, ...
    'tensor', {{}}, ...
    'output', {{'c'}}, ...
    'others', {{}});
ir.bounds = { struct('name','c','lower','0','upper','1') };
ir.interp = {};
ir.model = struct( ...
    'statements', {{}}, ...
    'equations',  {{ struct('expr', binop('-', name('c'), num(1)), 'primed', 0) }});
ir.simulate = struct('numPeriods',10,'numSamples',1, ...
    'initial', {{ struct('var','k','value','0.5'), struct('var','shock','value','1') }}, ...
    'varSimu', {{'c'}}, ...
    'transitions', {{ struct('state','k','expr','c') }});
ir.hooks = struct('preModel','','preIter','','postIter','', ...
                  'preJacCode','','postJacCode','','cxx','');
end
```
(`adaptiveFactor` is optional in the schema; non-adaptive bounds simply **omit** it. The validator rejects a *present* `[]` for a scalar field — `isscalar([])` is false — so always omit rather than set empty.)

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('ir','tValidateShape.m')); disp(table(r))"`
Expected: FAIL — `gdsge.ir.validate` undefined.

- [ ] **Step 3: Implement the structural validator**

`src/+gdsge/+ir/validate.m`:
```matlab
function report = validate(ir)
% VALIDATE  Check an IR struct against the schema descriptor.
%   Returns struct with .pass (logical) and .errors (cellstr of located messages).
%   Phase 1: structural shape, kinds, enums, tagged-union dispatch (this task);
%   refs, slots, and option invariants are added in Task 4.
s = gdsge.ir.schema();
nodes = s.nodes;
stmts = s.stmts;
errors = {};

vField(ir, s.root, 'ir');

report = struct('pass', isempty(errors), 'errors', {errors});

    function vField(val, spec, path)
        switch spec.kind
            case 'scalar'
                if ~(isnumeric(val) && isscalar(val))
                    err(path, 'expected numeric scalar');
                end
            case 'matrix'
                if ~(isnumeric(val) && ismatrix(val))
                    err(path, 'expected numeric matrix');
                end
            case 'text'
                if ~(ischar(val) || isstring(val))
                    err(path, 'expected text');
                end
            case 'enum'
                if ~ischar(val) || ~ismember(val, spec.values)
                    err(path, sprintf('expected one of {%s}', strjoin(spec.values, ',')));
                end
            case 'ref'
                % Resolution checked in Task 4; here only require text.
                if ~ischar(val)
                    err(path, 'ref must be text');
                end
            case 'list'
                if ~iscell(val)
                    err(path, 'expected list (cell array)');
                else
                    for i = 1:numel(val)
                        vField(val{i}, spec.item, sprintf('%s{%d}', path, i));
                    end
                end
            case 'map'
                if ~isstruct(val)
                    err(path, 'expected map (struct)');
                else
                    fn = fieldnames(val);
                    for i = 1:numel(fn)
                        vField(val.(fn{i}), spec.value, sprintf('%s.%s', path, fn{i}));
                    end
                end
            case 'tagged'
                vTagged(val, spec, path);
            case 'taggedList'
                if ~iscell(val)
                    err(path, 'expected list (cell array)');
                else
                    for i = 1:numel(val)
                        vTagged(val{i}, spec, sprintf('%s{%d}', path, i));
                    end
                end
            case 'struct'
                vStruct(val, spec, path);
            otherwise
                err(path, sprintf('schema bug: unknown kind "%s"', spec.kind));
        end
    end

    function vStruct(val, spec, path)
        if ~isstruct(val)
            err(path, 'expected struct'); return;
        end
        fn = fieldnames(spec.fields);
        for i = 1:numel(fn)
            nm = fn{i}; fs = spec.fields.(nm); p = sprintf('%s.%s', path, nm);
            if ~isfield(val, nm)
                if fs.required; err(p, 'required field missing'); end
            else
                vField(val.(nm), fs, p);
            end
        end
    end

    function vTagged(val, spec, path)
        registry = nodes; if strcmp(spec.registry, 'stmts'); registry = stmts; end
        if ~isstruct(val) || ~isfield(val, spec.tag)
            err(path, sprintf('expected tagged record with "%s"', spec.tag)); return;
        end
        tag = val.(spec.tag);
        if ~ischar(tag) || ~isfield(registry, tag)
            err(path, sprintf('unknown %s "%s"', spec.tag, toStr(tag))); return;
        end
        vStruct(val, registry.(tag), path);
    end

    function err(path, msg)
        errors{end+1} = sprintf('%s: %s', path, msg); %#ok<AGROW>
    end
end

function s = toStr(x)
if ischar(x); s = x; else; s = class(x); end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('ir','tValidateShape.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 6 tests.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+ir/validate.m tests/ir/tValidateShape.m tests/+gdsgefix/minimalIR.m
git commit -m "feat(ir): structural validator (shape, kinds, tagged dispatch) + shared minimalIR fixture"
```

---

## Task 4: Validator — refs, slots, option invariants

**Files:**
- Modify: `src/+gdsge/+ir/validate.m`
- Test: `tests/ir/tValidateSemantics.m`

- [ ] **Step 1: Write the failing test**

`tests/ir/tValidateSemantics.m`:
```matlab
classdef tValidateSemantics < matlab.unittest.TestCase
    methods (Test)
        function danglingTransRefFails(tc)
            ir = gdsgefix.minimalIR();
            ir.model.statements = { struct('type','reduction','kind','EXPECT', ...
                'target','e','body',gdsge.ir.node.primed('z'),'transRef','nope') };
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
            tc.verifyTrue(any(contains(r.errors, 'transRef')) || any(contains(r.errors,'nope')));
        end
        function resolvedTransRefPasses(tc)
            ir = gdsgefix.minimalIR();
            ir.model.statements = { struct('type','reduction','kind','EXPECT', ...
                'target','e','body',gdsge.ir.node.primed('z'),'transRef','shock_trans') };
            r = gdsge.ir.validate(ir);
            tc.verifyTrue(r.pass, strjoin(r.errors,' | '));
        end
        function unknownOutputRefFails(tc)
            ir = gdsgefix.minimalIR();
            ir.variables.output = {'c','ghost'};
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
            tc.verifyTrue(any(contains(r.errors, 'ghost')));
        end
        function overlappingSlotsFail(tc)
            ir = gdsgefix.minimalIR();
            ir.variables.policy = { ...
                struct('name','c','length',1,'slot',[1 1]), ...
                struct('name','d','length',1,'slot',[1 1]) };  % overlap at 1
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
            tc.verifyTrue(any(contains(r.errors, 'variables.policy')));
        end
        function lengthSlotMismatchFails(tc)
            ir = gdsgefix.minimalIR();
            ir.variables.policy = { struct('name','c','length',4,'slot',[1 2]) };
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
        end
        function badInterpOrderFails(tc)
            ir = gdsgefix.minimalIR();
            ir.options.interpOrder = 3;
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
            tc.verifyTrue(any(contains(r.errors, 'interpOrder')));
        end
        function asgWithTensorFails(tc)
            ir = gdsgefix.minimalIR();
            ir.options.interpMethod = 'asg';
            ir.options.asgMaxLevel = 10;
            ir.variables.tensor = {'r'};
            r = gdsge.ir.validate(ir);
            tc.verifyFalse(r.pass);
            tc.verifyTrue(any(contains(r.errors, 'tensor')));
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('ir','tValidateSemantics.m')); disp(table(r))"`
Expected: FAIL — refs not resolved, slots/options not checked yet.

- [ ] **Step 3: Extend the validator**

In `src/+gdsge/+ir/validate.m`, build ref pools before the walk and add the semantic checks after it. Replace the body between the `errors = {};` line and the `report = ...` line with:

```matlab
errors = {};
pools = buildPools(ir);

vField(ir, s.root, 'ir');
checkSlots();
checkOptionInvariants();

report = struct('pass', isempty(errors), 'errors', {errors});
```

Change the `'ref'` case inside `vField` to resolve against the pool:
```matlab
            case 'ref'
                if ~ischar(val)
                    err(path, 'ref must be text');
                elseif ~ismember(val, pools.(spec.pool))
                    err(path, sprintf('unresolved ref "%s" (pool %s)', val, spec.pool));
                end
```

Add these nested functions (place them alongside `vField`/`vStruct`, before the final `err` definition):
```matlab
    function checkSlots()
        checkCat('variables.policy', getCat('policy'));
        checkCat('variables.aux',    getCat('aux'));
    end

    function items = getCat(cat)
        items = {};
        if isfield(ir,'variables') && isfield(ir.variables, cat) && iscell(ir.variables.(cat))
            items = ir.variables.(cat);
        end
    end

    function checkCat(path, items)
        if isempty(items); return; end
        starts = zeros(numel(items),1); stops = zeros(numel(items),1); good = true;
        for i = 1:numel(items)
            it = items{i};
            if ~(isstruct(it) && all(isfield(it, {'slot','length','name'})))
                err(sprintf('%s{%d}', path, i), 'missing name/length/slot'); good = false; continue;
            end
            sl = it.slot;
            if numel(sl) ~= 2
                err(sprintf('%s{%d}.slot', path, i), 'slot must be [start stop]'); good = false; continue;
            end
            starts(i) = sl(1); stops(i) = sl(2);
            if sl(1) > sl(2)
                err(sprintf('%s{%d}.slot', path, i), 'start > stop'); good = false;
            end
            if sl(2) ~= sl(1) + it.length - 1
                err(sprintf('%s{%d}', path, i), ...
                    sprintf('length %g inconsistent with slot span %g', it.length, sl(2)-sl(1)+1));
                good = false;
            end
        end
        if ~good; return; end
        [starts, ord] = sort(starts); stops = stops(ord);
        if starts(1) ~= 1
            err(path, sprintf('slots must start at 1 (got %g)', starts(1)));
        end
        for i = 2:numel(starts)
            if starts(i) ~= stops(i-1) + 1
                err(path, sprintf('slot gap/overlap near index %g', starts(i)));
            end
        end
    end

    function checkOptionInvariants()
        if ~isfield(ir, 'options'); return; end
        o = ir.options;
        if isfield(o, 'interpOrder') && ~ismember(o.interpOrder, [2 4])
            err('options.interpOrder', 'must be 2 or 4');
        end
        if isfield(o, 'interpMethod') && strcmp(o.interpMethod, 'asg')
            if isfield(ir,'variables') && isfield(ir.variables,'tensor') ...
                    && iscell(ir.variables.tensor) && ~isempty(ir.variables.tensor)
                err('variables.tensor', 'must be empty when interpMethod = asg');
            end
            if ~isfield(o, 'asgMaxLevel')
                err('options.asgMaxLevel', 'required when interpMethod = asg');
            end
        end
    end
```

Add the pool builder as a plain local function at the end of the file (after the `toStr` function):
```matlab
function pools = buildPools(ir)
pools.transitions = {};
if isfield(ir,'shocks') && isfield(ir.shocks,'transitions') && isstruct(ir.shocks.transitions)
    pools.transitions = fieldnames(ir.shocks.transitions);
end
pools.states = {};
if isfield(ir,'states') && isfield(ir.states,'names')
    pools.states = cellstr(reshape(ir.states.names, 1, []));
end
pools.policyAux = {};
if isfield(ir,'variables')
    pools.policyAux = [listNames(getfieldor(ir.variables,'policy')), ...
                       listNames(getfieldor(ir.variables,'aux'))];
end
end

function names = listNames(lst)
names = {};
if iscell(lst)
    for i = 1:numel(lst)
        if isstruct(lst{i}) && isfield(lst{i}, 'name'); names{end+1} = lst{i}.name; end %#ok<AGROW>
    end
end
end

function v = getfieldor(s, f)
if isfield(s, f); v = s.(f); else; v = {}; end
end
```

- [ ] **Step 4: Run both validator tests to verify they pass**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests('ir'); assert(all(~[r.Failed]))"`
Expected: PASS — `tValidateShape` (6) and `tValidateSemantics` (7), plus `tNodeConstructors`/`tSchema`.

- [ ] **Step 5: Commit**

```bash
git add src/+gdsge/+ir/validate.m tests/ir/tValidateSemantics.m
git commit -m "feat(ir): validator refs, slot layout, and option-invariant checks"
```

---

## Task 5: Round-trip — canonicalize, decode, encode, isequalIR

**Files:**
- Create: `src/+gdsge/+ir/canonicalize.m`, `decode.m`, `encode.m`, `roundtrip.m`, `isequalIR.m`
- Test: `tests/ir/tRoundtrip.m`

- [ ] **Step 1: Write the failing test**

`tests/ir/tRoundtrip.m`:
```matlab
classdef tRoundtrip < matlab.unittest.TestCase
    methods (Test)
        function rowVectorKeepsOrientation(tc)
            ir = gdsgefix.minimalIR();
            ir.shocks.values.z = [0.99 1.05 0.99 1.05];   % 1x4 row
            ir2 = gdsge.ir.roundtrip(ir);
            tc.verifyEqual(size(ir2.shocks.values.z), [1 4]);
            tc.verifyEqual(ir2.shocks.values.z, [0.99 1.05 0.99 1.05], 'AbsTol', 1e-12);
        end
        function matrixSurvives(tc)
            ir = gdsgefix.minimalIR();
            M = magic(4) / 34;
            ir.shocks.transitions.shock_trans = M;
            ir2 = gdsge.ir.roundtrip(ir);
            tc.verifyEqual(size(ir2.shocks.transitions.shock_trans), [4 4]);
            tc.verifyEqual(ir2.shocks.transitions.shock_trans, M, 'AbsTol', 1e-12);
        end
        function infSurvives(tc)
            ir = gdsgefix.minimalIR();
            ir.options.saveFreq = Inf;
            ir2 = gdsge.ir.roundtrip(ir);
            tc.verifyEqual(ir2.options.saveFreq, Inf);
        end
        function nanSurvives(tc)
            ir = gdsgefix.minimalIR();
            ir.options.tolEq = NaN;
            ir2 = gdsge.ir.roundtrip(ir);
            tc.verifyTrue(isnan(ir2.options.tolEq));
        end
        function heterogeneousNodeListBecomesCell(tc)
            ir = gdsgefix.minimalIR();
            ir.model.statements = { struct('type','interpCall', ...
                'targets', {{'an'}}, 'primed', 1, ...
                'args', {{ gdsge.ir.node.primed('kn') }}, 'interpRef','GDSGE_INTERP_VEC') };
            ir2 = gdsge.ir.roundtrip(ir);
            tc.verifyTrue(iscell(ir2.model.statements));
            tc.verifyEqual(ir2.model.statements{1}.type, 'interpCall');
            tc.verifyTrue(iscell(ir2.model.statements{1}.args));
        end
        function emptyListRoundtrips(tc)
            ir = gdsgefix.minimalIR();
            ir2 = gdsge.ir.roundtrip(ir);
            tc.verifyTrue(iscell(ir2.variables.aux));
            tc.verifyEmpty(ir2.variables.aux);
        end
        function fullRoundtripIsEqual(tc)
            ir = gdsgefix.minimalIR();
            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));
        end
        function decodedIrStillValidates(tc)
            ir = gdsgefix.minimalIR();
            ir2 = gdsge.ir.roundtrip(ir);
            r = gdsge.ir.validate(ir2);
            tc.verifyTrue(r.pass, strjoin(r.errors,' | '));
        end
        function differentIrsAreNotEqual(tc)
            a = gdsgefix.minimalIR();
            b = gdsgefix.minimalIR(); b.modelName = 'other';
            tc.verifyFalse(gdsge.ir.isequalIR(a, b));
        end
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('ir','tRoundtrip.m')); disp(table(r))"`
Expected: FAIL — `gdsge.ir.canonicalize`/`encode`/`decode`/`roundtrip`/`isequalIR` undefined.

- [ ] **Step 3: Implement `canonicalize`**

`src/+gdsge/+ir/canonicalize.m`:
```matlab
function out = canonicalize(ir)
% CANONICALIZE  Normalize any IR (in-memory or freshly jsondecode'd) to a
%   canonical MATLAB form: object-arrays as row cell arrays of scalar structs,
%   non-finite scalars restored from tag strings, text as char. Idempotent.
s = gdsge.ir.schema();
nodes = s.nodes; stmts = s.stmts;
out = cField(ir, s.root);

    function v = cField(val, spec)
        switch spec.kind
            case 'scalar'
                v = canonScalar(val);
            case 'matrix'
                v = double(val);
            case {'text','enum','ref'}
                v = char(val);
            case 'list'
                c = ensureCell(val); v = cell(1, numel(c));
                for i = 1:numel(c); v{i} = cField(c{i}, spec.item); end
            case 'map'
                v = struct();
                if isstruct(val)
                    fn = fieldnames(val);
                    for i = 1:numel(fn); v.(fn{i}) = cField(val.(fn{i}), spec.value); end
                end
            case 'tagged'
                v = cTagged(val, spec);
            case 'taggedList'
                c = ensureCell(val); v = cell(1, numel(c));
                for i = 1:numel(c); v{i} = cTagged(c{i}, spec); end
            case 'struct'
                v = cStruct(val, spec.fields);
            otherwise
                error('gdsge:ir:canonicalize:kind', 'unknown kind %s', spec.kind);
        end
    end

    function v = cStruct(val, fields)
        v = struct();
        fn = fieldnames(fields);
        for i = 1:numel(fn)
            nm = fn{i};
            if isfield(val, nm)
                v.(nm) = cField(val.(nm), fields.(nm));
            end
        end
    end

    function v = cTagged(val, spec)
        registry = nodes; if strcmp(spec.registry, 'stmts'); registry = stmts; end
        tag = val.(spec.tag);
        v = cStruct(val, registry.(char(tag)).fields);
    end
end

function v = canonScalar(val)
if ischar(val) || isstring(val)
    switch char(val)
        case 'Inf';  v = Inf;
        case '-Inf'; v = -Inf;
        case 'NaN';  v = NaN;
        otherwise;   v = str2double(char(val));
    end
else
    v = double(val);
end
end

function c = ensureCell(val)
if iscell(val)
    c = reshape(val, 1, []);
elseif isstruct(val)
    c = reshape(num2cell(val), 1, []);
elseif isempty(val)
    c = {};
else
    c = {val};
end
end
```

- [ ] **Step 4: Implement `decode`, `encode`, `roundtrip`, `isequalIR`**

`src/+gdsge/+ir/decode.m`:
```matlab
function ir = decode(txt)
% DECODE  JSON text -> canonical IR struct.
ir = gdsge.ir.canonicalize(jsondecode(txt));
end
```

`src/+gdsge/+ir/encode.m`:
```matlab
function txt = encode(ir)
% ENCODE  IR struct -> pretty JSON text. Matrices are emitted as arrays-of-rows
%   (orientation-preserving); non-finite scalars are tagged strings.
s = gdsge.ir.schema();
nodes = s.nodes; stmts = s.stmts;
jr = jField(gdsge.ir.canonicalize(ir), s.root);
txt = jsonencode(jr, 'PrettyPrint', true);

    function v = jField(val, spec)
        switch spec.kind
            case 'scalar'
                v = jScalar(val);
            case 'matrix'
                v = wrapMatrix(val);
            case {'text','enum','ref'}
                v = char(val);
            case 'list'
                v = cell(1, numel(val));
                for i = 1:numel(val); v{i} = jField(val{i}, spec.item); end
            case 'map'
                v = struct(); fn = fieldnames(val);
                for i = 1:numel(fn); v.(fn{i}) = jField(val.(fn{i}), spec.value); end
            case 'tagged'
                v = jTagged(val, spec);
            case 'taggedList'
                v = cell(1, numel(val));
                for i = 1:numel(val); v{i} = jTagged(val{i}, spec); end
            case 'struct'
                v = jStruct(val, spec.fields);
            otherwise
                error('gdsge:ir:encode:kind', 'unknown kind %s', spec.kind);
        end
    end

    function v = jStruct(val, fields)
        v = struct(); fn = fieldnames(fields);
        for i = 1:numel(fn)
            nm = fn{i};
            if isfield(val, nm); v.(nm) = jField(val.(nm), fields.(nm)); end
        end
    end

    function v = jTagged(val, spec)
        registry = nodes; if strcmp(spec.registry, 'stmts'); registry = stmts; end
        v = jStruct(val, registry.(char(val.(spec.tag))).fields);
    end
end

function v = jScalar(x)
if ~isfinite(x)
    if isnan(x); v = 'NaN';
    elseif x > 0; v = 'Inf';
    else; v = '-Inf';
    end
else
    v = x;
end
end

function out = wrapMatrix(M)
% Emit as a 1xR cell of row vectors so jsonencode produces [[...],[...]]
% and jsondecode reconstructs the exact 2-D shape.
[r, ~] = size(M);
out = cell(1, r);
for i = 1:r; out{i} = M(i, :); end
end
```

`src/+gdsge/+ir/roundtrip.m`:
```matlab
function ir2 = roundtrip(ir)
% ROUNDTRIP  decode(encode(ir)); the workhorse of the serialization tests.
ir2 = gdsge.ir.decode(gdsge.ir.encode(ir));
end
```

`src/+gdsge/+ir/isequalIR.m`:
```matlab
function tf = isequalIR(a, b)
% ISEQUALIR  Structural IR equality through the canonical lens
%   (orientation- and array-kind-tolerant; numeric leaves compared with a
%   small tolerance to absorb JSON precision).
tf = deepEqual(gdsge.ir.canonicalize(a), gdsge.ir.canonicalize(b));
end

function tf = deepEqual(a, b)
if isstruct(a) && isstruct(b)
    fa = sort(fieldnames(a)); fb = sort(fieldnames(b));
    if ~isequal(fa, fb); tf = false; return; end
    tf = true;
    for i = 1:numel(fa)
        if ~deepEqual(a.(fa{i}), b.(fa{i})); tf = false; return; end
    end
elseif iscell(a) && iscell(b)
    if numel(a) ~= numel(b); tf = false; return; end
    tf = true;
    for i = 1:numel(a)
        if ~deepEqual(a{i}, b{i}); tf = false; return; end
    end
elseif (ischar(a) || isstring(a)) && (ischar(b) || isstring(b))
    tf = strcmp(char(a), char(b));
elseif isnumeric(a) && isnumeric(b)
    if ~isequal(size(a), size(b)); tf = false; return; end
    if isempty(a); tf = true; return; end
    tf = isequaln(a, b) || all(abs(a(:) - b(:)) <= 1e-12 + 1e-12 .* abs(b(:)));
else
    tf = isequaln(a, b);
end
end
```

- [ ] **Step 5: Run the round-trip tests to verify they pass**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('ir','tRoundtrip.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 9 tests.

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+ir/canonicalize.m src/+gdsge/+ir/decode.m src/+gdsge/+ir/encode.m src/+gdsge/+ir/roundtrip.m src/+gdsge/+ir/isequalIR.m tests/ir/tRoundtrip.m
git commit -m "feat(ir): descriptor-driven canonicalize/encode/decode round-trip + isequalIR"
```

---

## Task 6: Doc generator

**Files:**
- Create: `src/+gdsge/+ir/gendoc.m`
- Modify (regenerate): `docs/ir-schema.md`
- Test: `tests/ir/tGendoc.m`

- [ ] **Step 1: Write the failing test**

`tests/ir/tGendoc.m`:
```matlab
classdef tGendoc < matlab.unittest.TestCase
    methods (Test)
        function coversEverySection(tc)
            md = gdsge.ir.gendoc();
            for sec = {'options','shocks','states','variables','bounds', ...
                       'interp','model','simulate','hooks'}
                tc.verifyTrue(contains(md, sec{1}), ...
                    sprintf('doc missing section %s', sec{1}));
            end
            tc.verifyTrue(contains(md, 'binop'));  % node registry rendered
            tc.verifyTrue(contains(md, 'reduction'));  % stmt registry rendered
        end
        function committedDocIsUpToDate(tc)
            md = gdsge.ir.gendoc();
            here = fileparts(mfilename('fullpath'));        % tests/ir
            docPath = fullfile(fileparts(fileparts(here)), 'docs', 'ir-schema.md');
            onDisk = fileread(docPath);
            tc.verifyEqual(normalizeNewlines(md), normalizeNewlines(onDisk), ...
                'docs/ir-schema.md is stale — regenerate with gdsge.ir.gendoc(path).');
        end
    end
end

function s = normalizeNewlines(s)
s = strrep(s, sprintf('\r\n'), sprintf('\n'));
s = strtrim(s);
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('ir','tGendoc.m')); disp(table(r))"`
Expected: FAIL — `gdsge.ir.gendoc` undefined.

- [ ] **Step 3: Implement `gendoc`**

`src/+gdsge/+ir/gendoc.m`:
```matlab
function md = gendoc(targetPath)
% GENDOC  Render the IR schema descriptor as Markdown.
%   md = gdsge.ir.gendoc()        returns the markdown string.
%   gdsge.ir.gendoc(targetPath)   also writes it to TARGETPATH (utf-8, LF).
s = gdsge.ir.schema();
L = {};
L{end+1} = '# IR Schema (generated — do not edit by hand)';
L{end+1} = '';
L{end+1} = sprintf('Generated from `src/+gdsge/+ir/schema.m` by `gdsge.ir.gendoc`. irVersion: `%s`.', s.irVersion);
L{end+1} = '';
L{end+1} = '## Document sections';
L{end+1} = '';
L = renderStruct(L, s.root, 0);
L{end+1} = '';
L{end+1} = '## AST nodes';
L{end+1} = '';
L = renderRegistry(L, s.nodes);
L{end+1} = '';
L{end+1} = '## Model statements';
L{end+1} = '';
L = renderRegistry(L, s.stmts);
md = strjoin(L, sprintf('\n'));
if nargin >= 1 && ~isempty(targetPath)
    fid = fopen(targetPath, 'w', 'n', 'UTF-8');
    if fid < 0; error('gdsge:ir:gendoc:open', 'cannot write %s', targetPath); end
    fwrite(fid, md);
    fclose(fid);
end
end

function L = renderRegistry(L, reg)
fn = fieldnames(reg);
for i = 1:numel(fn)
    L{end+1} = sprintf('### `%s`', fn{i}); %#ok<AGROW>
    L{end+1} = ''; %#ok<AGROW>
    L = renderStruct(L, reg.(fn{i}), 0);
    L{end+1} = ''; %#ok<AGROW>
end
end

function L = renderStruct(L, spec, depth)
fn = fieldnames(spec.fields);
for i = 1:numel(fn)
    nm = fn{i}; fs = spec.fields.(nm);
    L{end+1} = line(depth, nm, fs); %#ok<AGROW>
    L = renderChildren(L, fs, depth + 1);
end
end

function L = renderChildren(L, fs, depth)
switch fs.kind
    case 'struct'
        L = renderStruct(L, fs, depth);
    case 'list'
        if strcmp(fs.item.kind, 'struct')
            L{end+1} = pad(depth, '- *(each item)*'); %#ok<AGROW>
            L = renderStruct(L, fs.item, depth + 1);
        end
    case 'map'
        if strcmp(fs.value.kind, 'struct')
            L{end+1} = pad(depth, '- *(each value)*'); %#ok<AGROW>
            L = renderStruct(L, fs.value, depth + 1);
        end
end
end

function s = line(depth, nm, fs)
extra = '';
switch fs.kind
    case 'enum'; extra = sprintf(' {%s}', strjoin(fs.values, ', '));
    case 'ref';  extra = sprintf(' -> pool:%s', fs.pool);
    case 'list'; extra = sprintf(' of %s', fs.item.kind);
    case 'map';  extra = sprintf(' of %s', fs.value.kind);
    case {'tagged','taggedList'}; extra = sprintf(' (%s)', fs.registry);
end
req = 'required'; if isfield(fs,'required') && ~fs.required; req = 'optional'; end
s = pad(depth, sprintf('- `%s` : %s%s (%s)', nm, fs.kind, extra, req));
end

function s = pad(depth, txt)
s = [repmat('  ', 1, depth), txt];
end
```

- [ ] **Step 4: Generate the committed doc**

From the repo root:
```bash
matlab -batch "addpath('src'); gdsge.ir.gendoc(fullfile('docs','ir-schema.md'));"
```
Expected: `docs/ir-schema.md` is overwritten with the generated reference (replacing the Phase-0 stub).

- [ ] **Step 5: Run the gendoc test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('ir','tGendoc.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 2 tests (the on-disk doc now matches `gendoc()`).

- [ ] **Step 6: Commit**

```bash
git add src/+gdsge/+ir/gendoc.m docs/ir-schema.md tests/ir/tGendoc.m
git commit -m "feat(ir): schema doc generator; regenerate docs/ir-schema.md (no-drift test)"
```

---

## Task 7: HL1996 reference IR (end-to-end contract)

**Files:**
- Create: `tests/HeatonLucas1996/ir/buildHL1996IR.m`
- Create (generated): `tests/HeatonLucas1996/ir/HL1996.gdsge.json`
- Test: `tests/HeatonLucas1996/ir/tIrHL1996.m`

- [ ] **Step 1: Write the failing test**

`tests/HeatonLucas1996/ir/tIrHL1996.m`:
```matlab
classdef tIrHL1996 < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addThisFolder(tc)
            here = fileparts(mfilename('fullpath'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(here));
        end
    end
    methods (Test)
        function referenceIrValidates(tc)
            ir = buildHL1996IR();
            r = gdsge.ir.validate(ir);
            tc.verifyTrue(r.pass, strjoin(r.errors, ' | '));
        end
        function squareSystemNineteen(tc)
            % Sanity: 19 unknowns and (after prime expansion) 19 equations.
            ir = buildHL1996IR();
            nUnknown = 0;
            for i = 1:numel(ir.variables.policy); nUnknown = nUnknown + ir.variables.policy{i}.length; end
            tc.verifyEqual(nUnknown, 19);
            nEq = 0;
            for i = 1:numel(ir.model.equations)
                e = ir.model.equations{i};
                if e.primed; nEq = nEq + ir.shocks.count; else; nEq = nEq + 1; end
            end
            tc.verifyEqual(nEq, 19);
        end
        function roundtripIsEqual(tc)
            ir = buildHL1996IR();
            tc.verifyTrue(gdsge.ir.isequalIR(ir, gdsge.ir.roundtrip(ir)));
        end
        function encodingMatchesGolden(tc)
            ir = buildHL1996IR();
            here = fileparts(mfilename('fullpath'));
            onDisk = fileread(fullfile(here, 'HL1996.gdsge.json'));
            tc.verifyEqual(norm_(gdsge.ir.encode(ir)), norm_(onDisk), ...
                'HL1996.gdsge.json is stale — regenerate it.');
        end
    end
end

function s = norm_(s)
s = strtrim(strrep(s, sprintf('\r\n'), sprintf('\n')));
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('HeatonLucas1996','ir','tIrHL1996.m')); disp(table(r))"`
Expected: FAIL — `buildHL1996IR` undefined.

- [ ] **Step 3: Implement the reference-IR builder**

`tests/HeatonLucas1996/ir/buildHL1996IR.m`:
```matlab
function ir = buildHL1996IR()
% BUILDHL1996IR  Hand-authored reference IR for HeatonLucas1996.
%   This is the canonical Phase-1 fixture AND the forward target the Phase-3
%   parser must reproduce (parser output must gdsge.ir.isequalIR this struct).
%   Option values below are the resolved defaults for this model; reconcile the
%   exact toolbox defaults against the parser in Phase 3.
N  = @gdsge.ir.node.num;
NM = @gdsge.ir.node.name;
P  = @gdsge.ir.node.primed;
B  = @gdsge.ir.node.binop;
U  = @gdsge.ir.node.unop;

ir.irVersion = '1.0.0';
ir.modelName = 'HL1996';

ir.options = struct('interpMethod','spline','interpOrder',4,'extrapOrder',2, ...
    'tolEq',1e-6,'numThreads',8,'simuResolve',0,'simuInterp',1, ...
    'printFreq',100,'saveFreq',Inf);

% ----- shocks -------------------------------------------------------------
g    = [.9904 1.0470 .9904 1.0470 .9904 1.0470 .9904 1.0470];
d    = [.1402 .1437 .1561 .1599 .1402 .1437 .1561 .1599];
eta1 = [.3772 .3772 .3772 .3772 .6228 .6228 .6228 .6228];
rawTrans = [ ...
  0.3932 0.2245 0.0793 0.0453 0.1365 0.0779 0.0275 0.0157
  0.3044 0.3470 0.0425 0.0484 0.1057 0.1205 0.0147 0.0168
  0.0484 0.0425 0.3470 0.3044 0.0168 0.0147 0.1205 0.1057
  0.0453 0.0793 0.2245 0.3932 0.0157 0.0275 0.0779 0.1365
  0.1365 0.0779 0.0275 0.0157 0.3932 0.2245 0.0793 0.0453
  0.1057 0.1205 0.0147 0.0168 0.3044 0.3470 0.0425 0.0484
  0.0168 0.0147 0.1205 0.1057 0.0484 0.0425 0.3470 0.3044
  0.0157 0.0275 0.0779 0.1365 0.0453 0.0793 0.2245 0.3932 ];
shock_trans = rawTrans ./ sum(rawTrans, 2);   % normalized, as the gmod does
ir.shocks = struct('names', {{'g','d','eta1'}}, 'count', 8, ...
    'values', struct('g', g, 'd', d, 'eta1', eta1), ...
    'transitions', struct('shock_trans', shock_trans));

% ----- states -------------------------------------------------------------
ir.states = struct('names', {{'w1'}}, ...
    'grids', struct('w1', 'linspace(-0.05,1.05,201)'));

% ----- variables + slot layout -------------------------------------------
pol = { ...
    var('c1',1,[1 1]),   var('c2',1,[2 2]),   var('s1p',1,[3 3]), ...
    var('nb1p',1,[4 4]), var('nb2p',1,[5 5]), var('ms1',1,[6 6]), ...
    var('ms2',1,[7 7]),  var('mb1',1,[8 8]),  var('mb2',1,[9 9]), ...
    var('ps',1,[10 10]), var('pb',1,[11 11]), var('w1n',8,[12 19]) };
ir.variables = struct( ...
    'policy', {pol}, ...
    'aux',    {{ var('equity_premium',1,[1 1]) }}, ...
    'interp', {{ 'ps_future','pb_future','c1_future','c2_future' }}, ...
    'tensor', {{}}, ...
    'output', {{ 'c1','c2','ps','pb','equity_premium','w1n' }}, ...
    'others', {{}});

% ----- bounds -------------------------------------------------------------
ir.bounds = { ...
    bnd('c1','0','1'),   bnd('c2','0','1'),   bnd('s1p','0.0','1.0'), ...
    bnd('nb1p','0.0','1.0'), bnd('nb2p','0.0','1.0'), ...
    bnd('ms1','0','1'),  bnd('ms2','0','1'),  bnd('mb1','0','1'), bnd('mb2','0','1'), ...
    bndA('ps','0','3',1.5), bndA('pb','0','3',1.5), bnd('w1n','-0.5','1.5') };

% ----- interpolants -------------------------------------------------------
ir.interp = { ...
    interpObj('ps_future','0.0','ps'), ...
    interpObj('pb_future','0.0','pb'), ...
    interpObj('c1_future','w1.*d+eta1','c1'), ...
    interpObj('c2_future','(1-w1).*d+1-eta1','c2') };

% ----- model statements ---------------------------------------------------
gPow1mGamma = B('^', P('g'), B('-', N(1), NM('gamma')));   % g'^(1-gamma)
gPowMGamma  = B('^', P('g'), U('-', NM('gamma')));         % g'^(-gamma)
psnPlusD    = B('+', P('psn'), P('d'));                    % (psn'+d')
% es1 body: g'^(1-gamma)*(c1n'/c1)^(-gamma)*(psn'+d')/ps
es1Body = B('/', B('*', B('*', gPow1mGamma, ...
            B('^', B('/', P('c1n'), NM('c1')), U('-', NM('gamma')))), psnPlusD), NM('ps'));
es2Body = B('/', B('*', B('*', gPow1mGamma, ...
            B('^', B('/', P('c2n'), NM('c2')), U('-', NM('gamma')))), psnPlusD), NM('ps'));
eb1Body = B('/', B('*', gPowMGamma, ...
            B('^', B('/', P('c1n'), NM('c1')), U('-', NM('gamma')))), NM('pb'));
eb2Body = B('/', B('*', gPowMGamma, ...
            B('^', B('/', P('c2n'), NM('c2')), U('-', NM('gamma')))), NM('pb'));
epBody  = B('*', B('/', psnPlusD, NM('ps')), P('g'));      % (psn'+d')/ps*g'

% budget_1 = w1*(ps+d)+eta1 - c1 - ps*s1p - pb*b1p
budget1 = B('-', B('-', B('-', B('+', B('*', NM('w1'), B('+', NM('ps'), NM('d'))), NM('eta1')), ...
            NM('c1')), B('*', NM('ps'), NM('s1p'))), B('*', NM('pb'), NM('b1p')));
% budget_2 = (1-w1)*(ps+d)+(1-eta1) - c2 - ps*s2p - pb*b2p
budget2 = B('-', B('-', B('-', B('+', B('*', B('-', N(1), NM('w1')), B('+', NM('ps'), NM('d'))), ...
            B('-', N(1), NM('eta1'))), NM('c2')), B('*', NM('ps'), NM('s2p'))), B('*', NM('pb'), NM('b2p')));
% w1_consis' = (s1p*(psn'+d') + b1p/g')/(psn'+d') - w1n'
w1consis = B('-', B('/', B('+', B('*', NM('s1p'), psnPlusD), B('/', NM('b1p'), P('g'))), psnPlusD), P('w1n'));

ir.model.statements = { ...
    stInterp({'psn','pbn','c1n','c2n'}, { P('w1n') }), ...
    stReduce('EXPECT','es1', es1Body), ...
    stReduce('EXPECT','es2', es2Body), ...
    stReduce('EXPECT','eb1', eb1Body), ...
    stReduce('EXPECT','eb2', eb2Body), ...
    stAssign('b1p', 0, B('+', NM('nb1p'), NM('Kb'))), ...
    stAssign('b2p', 0, B('+', NM('nb2p'), NM('Kb'))), ...
    stAssign('s2p', 0, B('-', N(1), NM('s1p'))), ...
    stAssign('budget_1', 0, budget1), ...
    stAssign('budget_2', 0, budget2), ...
    stAssign('w1_consis', 1, w1consis), ...
    stReduce('EXPECT','__ep_expect', epBody), ...
    stAssign('equity_premium', 0, B('-', NM('__ep_expect'), B('/', N(1), NM('pb')))) };

% ----- equations (residuals) ---------------------------------------------
ir.model.equations = { ...
    eqn(B('+', B('+', U('-', N(1)), B('*', NM('beta'), NM('es1'))), NM('ms1')), 0), ...
    eqn(B('+', B('+', U('-', N(1)), B('*', NM('beta'), NM('es2'))), NM('ms2')), 0), ...
    eqn(B('+', B('+', U('-', N(1)), B('*', NM('beta'), NM('eb1'))), NM('mb1')), 0), ...
    eqn(B('+', B('+', U('-', N(1)), B('*', NM('beta'), NM('eb2'))), NM('mb2')), 0), ...
    eqn(B('*', NM('ms1'), NM('s1p')), 0), ...
    eqn(B('*', NM('ms2'), NM('s2p')), 0), ...
    eqn(B('*', NM('mb1'), NM('nb1p')), 0), ...
    eqn(B('*', NM('mb2'), NM('nb2p')), 0), ...
    eqn(B('+', NM('b1p'), NM('b2p')), 0), ...
    eqn(NM('budget_1'), 0), ...
    eqn(NM('budget_2'), 0), ...
    eqn(NM('w1_consis'), 1) };

% ----- simulate -----------------------------------------------------------
ir.simulate = struct('numPeriods', 10000, 'numSamples', 24, ...
    'initial', {{ struct('var','w1','value','0.5'), struct('var','shock','value','1') }}, ...
    'varSimu', {{ 'c1','c2','ps','pb','equity_premium' }}, ...
    'transitions', {{ struct('state','w1','expr','w1n') }});

% ----- hooks --------------------------------------------------------------
ir.hooks = struct('preModel','','preIter','','postIter','', ...
                  'preJacCode','','postJacCode','','cxx','');
end

% ===== small builders =====================================================
function v = var(name, len, slot)
v = struct('name', name, 'length', len, 'slot', slot);
end
function b = bnd(name, lo, hi)
b = struct('name', name, 'lower', lo, 'upper', hi);   % non-adaptive: omit optional field
end
function b = bndA(name, lo, hi, f)
b = struct('name', name, 'lower', lo, 'upper', hi, 'adaptiveFactor', f);
end
function o = interpObj(name, initExpr, updExpr)
o = struct('name', name, 'args', {{'w1'}}, 'initialExpr', initExpr, 'updateExpr', updExpr);
end
function s = stAssign(target, primed, expr)
s = struct('type','assign','target',target,'primed',primed,'expr',expr);
end
function s = stReduce(kind, target, body)
s = struct('type','reduction','kind',kind,'target',target,'body',body,'transRef','shock_trans');
end
function s = stInterp(targets, args)
s = struct('type','interpCall','targets',{targets},'primed',1,'args',{args}, ...
           'interpRef','GDSGE_INTERP_VEC');
end
function e = eqn(expr, primed)
e = struct('expr', expr, 'primed', primed);
end
```

- [ ] **Step 4: Verify the builder produces a valid, square IR (before generating the golden)**

Run: `matlab -batch "addpath('src'); addpath('tests/HeatonLucas1996/ir'); ir=buildHL1996IR(); r=gdsge.ir.validate(ir); assert(r.pass, strjoin(r.errors,' | ')); disp('valid');"`
Expected: prints `valid`. If validation fails, fix the builder (slots, refs, node shapes) before proceeding — do **not** generate a golden from an invalid IR.

- [ ] **Step 5: Generate the committed JSON golden**

Run: `matlab -batch "addpath('src'); addpath('tests/HeatonLucas1996/ir'); ir=buildHL1996IR(); txt=gdsge.ir.encode(ir); fid=fopen('tests/HeatonLucas1996/ir/HL1996.gdsge.json','w','n','UTF-8'); fwrite(fid,txt); fclose(fid); disp('wrote golden');"`
Expected: `tests/HeatonLucas1996/ir/HL1996.gdsge.json` is written.

- [ ] **Step 6: Run the HL1996 reference-IR test to verify it passes**

Run: `matlab -batch "cd('tests'); addpath(pwd); addpath(fullfile(fileparts(pwd),'src')); r=runtests(fullfile('HeatonLucas1996','ir','tIrHL1996.m')); assert(all(~[r.Failed]))"`
Expected: PASS — 4 tests (validates, square 19×19, round-trips, golden matches).

- [ ] **Step 7: Commit**

```bash
git add tests/HeatonLucas1996/ir/buildHL1996IR.m tests/HeatonLucas1996/ir/HL1996.gdsge.json tests/HeatonLucas1996/ir/tIrHL1996.m
git commit -m "test(ir): hand-authored HL1996 reference IR + JSON golden (Phase 3 target)"
```

---

## Task 8: Full suite green + PROGRESS update

**Files:**
- Modify: `PROGRESS.md`

- [ ] **Step 1: Run the entire suite headless**

Run: `pwsh -File tests/run.ps1`
Expected: exit 0. All green: Phase-0 (`tSmoke`, `tCompareNumericClose`, `tGoldenHL1996`) **and** Phase-1 (`tNodeConstructors`, `tSchema`, `tValidateShape`, `tValidateSemantics`, `tRoundtrip`, `tGendoc`, `tIrHL1996`). `tests/results/junit.xml` written.

- [ ] **Step 2: Update PROGRESS.md**

In `PROGRESS.md`, change the Phase 1 line from `◐` to `☑` with a completion date, and flip Phase 2 to `NEXT`:
```markdown
- ☑ **Phase 1 — IR schema + MATLAB scaffolding** (done 2026-06-12)
  - ☑ Declarative schema descriptor (`gdsge.ir.schema`) — single source of truth
  - ☑ AST node API (`gdsge.ir.node.*`) + accessors
  - ☑ Descriptor-driven validator (shape, refs, slots, option invariants)
  - ☑ Type-aware round-trip (`canonicalize`/`encode`/`decode`/`isequalIR`); the 4 quirks defeated
  - ☑ Doc generator (`gdsge.ir.gendoc` → `docs/ir-schema.md`, no-drift test)
  - ☑ HL1996 reference IR + JSON golden (Phase 3's forward target)
- ◐ **Phase 2 — Parser front-end** (macros + lexer + block-split + declaration parser → partial IR) — NEXT
```

Add to the Changelog section:
```markdown
- 2026-06-12: **Phase 1 complete.** `gdsge.ir` package: descriptor-driven schema,
  validator, type-aware JSON round-trip, doc generator, AST node API, and a
  hand-authored HL1996 reference IR (Phase 3's golden target). All on branch
  `phase1-ir-schema`. Next: Phase 2 (parser front-end).
```

Add to the Decisions log section:
```markdown
- **IR schema:** one declarative descriptor (`gdsge.ir.schema`) drives validation,
  round-trip, and doc-gen. Round-trip quirks (orientation, array-kind, empty/scalar,
  non-finite) handled type-aware in `canonicalize`/`toJsonReady`. Two expression worlds:
  model body → AST; grids/bounds/initial/transitions → opaque MATLAB text. Reductions are
  flat statements (nested ones hoisted by the parser). Validator checks IR-computable
  properties; model semantics (square system, bounds-completeness) stay parser-side.
```

- [ ] **Step 3: Commit**

```bash
git add PROGRESS.md
git commit -m "docs: mark Phase 1 complete in PROGRESS.md"
```

---

## Phase 1 done-when

- `pwsh -File tests/run.ps1` exits 0 with all Phase-0 and Phase-1 tests green; `tests/results/junit.xml` written.
- `gdsge.ir.validate` passes on the HL1996 reference IR; `gdsge.ir.isequalIR(ir, roundtrip(ir))` holds; `HL1996.gdsge.json` golden matches `encode`.
- `docs/ir-schema.md` is regenerated from `schema.m` and the `tGendoc` no-drift test passes.
- `PROGRESS.md` marks Phase 1 done and Phase 2 next; IR decisions logged.

## Self-review notes

- **Spec coverage:** module surface §3 → Tasks 1–6; IR document shape §4 → schema (Task 2) + HL1996 reference (Task 7); two-expression-worlds + reduction-hoisting §4 → encoded in schema (`text` vs `tagged` fields) and the HL1996 builder (`__ep_expect` hoist); round-trip quirks §5.2–5.5 → Task 5 tests (orientation, matrix, Inf, NaN, heterogeneous list, empty); validator boundary §6 → Tasks 3–4 (square-system/bounds-completeness deliberately *absent*; the HL1996 test asserts squareness as a fixture sanity check, not via the `+ir` validator); versioning §7 → schema `irVersion` + `tSchema` semver check; testing §8 → Tasks 3–7; done-when §9 → Task 8.
- **Vocabulary deviation flagged:** `optional`→`required` flag, `node/nodeList`→`tagged/taggedList`, added `map`. Documented in the file-structure note; functionally covers the spec's vocabulary.
- **Type consistency:** `validate` returns `{pass, errors}` (used identically in every test). Node constructors `gdsge.ir.node.{num,name,primed,unop,binop,call,index}` and accessors `kindOf`/`children` match across Tasks 1 and 7. Field kinds emitted by `schema` (Task 2) are exactly those dispatched in `validate` (Tasks 3–4), `canonicalize`/`encode` (Task 5), and `gendoc` (Task 6). Ref pools `transitions`/`states`/`policyAux` defined in Task 4's `buildPools` match the `fRef(...)` pools declared in Task 2's schema.
- **No placeholders:** every code step contains complete, runnable content; every run step has an exact command and expected result.
