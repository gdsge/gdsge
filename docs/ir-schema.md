# IR Schema (generated — do not edit by hand)

Generated from `src/+gdsge/+ir/schema.m` by `gdsge.ir.gendoc`. irVersion: `2.0.0`.

## Document sections

- `irVersion` : text (required)
- `modelName` : text (required)
- `params` : list of struct (required)
  - *(each item)*
    - `name` : text (required)
    - `value` : matrix (required)
- `setup` : list of struct (optional)
  - *(each item)*
    - `kind` : text (required)
    - `body` : text (required)
- `setupNames` : list of text (optional)
- `options` : struct (required)
  - `interpMethod` : enum {spline, asg, pchip} (required)
  - `jacobianBackend` : enum {autodiff, sympy} (optional)
  - `interpOrder` : scalar (optional)
  - `extrapOrder` : scalar (optional)
  - `asgMaxLevel` : scalar (optional)
  - `asgThreshold` : scalar (optional)
  - `asgMinLevel` : scalar (optional)
  - `asgOutputMaxLevel` : scalar (optional)
  - `asgOutputThreshold` : scalar (optional)
  - `tolEq` : scalar (optional)
  - `numThreads` : scalar (optional)
  - `simuResolve` : scalar (optional)
  - `simuInterp` : scalar (optional)
  - `printFreq` : scalar (optional)
  - `resolvePrintFreq` : scalar (optional)
  - `saveFreq` : scalar (optional)
  - `simuPrintFreq` : scalar (optional)
  - `simuSaveFreq` : scalar (optional)
- `shocks` : struct (required)
  - `names` : list of text (required)
  - `count` : scalar (required)
  - `values` : map of matrix (required)
  - `transitions` : map of matrix (required)
- `states` : struct (required)
  - `names` : list of text (required)
  - `grids` : map of text (required)
- `variables` : struct (required)
  - `policy` : list of struct (required)
    - *(each item)*
      - `name` : text (required)
      - `length` : scalar (required)
      - `slot` : matrix (required)
  - `aux` : list of struct (required)
    - *(each item)*
      - `name` : text (required)
      - `length` : scalar (required)
      - `slot` : matrix (required)
  - `interp` : list of ref -> pool:interpNames (required)
  - `tensor` : list of struct (required)
    - *(each item)*
      - `name` : text (required)
      - `expr` : text (required)
  - `output` : list of ref -> pool:policyAux (required)
  - `others` : list of text (required)
- `bounds` : list of struct (required)
  - *(each item)*
    - `name` : ref -> pool:policyAux (required)
    - `lower` : text (required)
    - `upper` : text (required)
    - `adaptiveFactor` : scalar (optional)
- `interp` : list of struct (required)
  - *(each item)*
    - `name` : text (required)
    - `args` : list of ref -> pool:states (required)
    - `initialExpr` : text (required)
    - `updateExpr` : text (required)
- `model` : struct (required)
  - `regions` : list of struct (required)
    - *(each item)*
      - `condition` : text (required)
      - `statements` : taggedList (stmts) (required)
      - `equations` : taggedList (eqkinds) (required)
- `modelInit` : struct (optional)
  - `variables` : struct (required)
    - `policyInit` : list of struct (required)
      - *(each item)*
        - `name` : text (required)
        - `length` : scalar (required)
        - `slot` : matrix (required)
    - `auxInit` : list of struct (required)
      - *(each item)*
        - `name` : text (required)
        - `length` : scalar (required)
        - `slot` : matrix (required)
  - `bounds` : list of struct (required)
    - *(each item)*
      - `name` : text (required)
      - `lower` : text (required)
      - `upper` : text (required)
      - `adaptiveFactor` : scalar (optional)
  - `statements` : taggedList (stmts) (required)
  - `equations` : taggedList (eqkinds) (required)
- `simulate` : struct (required)
  - `numPeriods` : scalar (required)
  - `numSamples` : scalar (required)
  - `initial` : list of struct (required)
    - *(each item)*
      - `var` : text (required)
      - `value` : text (required)
  - `varSimu` : list of ref -> pool:policyAux (required)
  - `transitions` : list of struct (required)
    - *(each item)*
      - `state` : ref -> pool:states (required)
      - `expr` : text (required)
      - `primed` : scalar (required)
- `hooks` : struct (required)
  - `preModel` : text (optional)
  - `preIter` : text (optional)
  - `postIter` : text (optional)
  - `preJacCode` : text (optional)
  - `postJacCode` : text (optional)
  - `cxx` : text (optional)
  - `cxxIncludes` : list of text (optional)

## AST nodes

### `num`

- `kind` : enum {num} (required)
- `value` : scalar (required)

### `name`

- `kind` : enum {name} (required)
- `id` : text (required)

### `primed`

- `kind` : enum {primed} (required)
- `id` : text (required)

### `unop`

- `kind` : enum {unop} (required)
- `op` : text (required)
- `arg` : tagged (nodes) (required)

### `binop`

- `kind` : enum {binop} (required)
- `op` : text (required)
- `lhs` : tagged (nodes) (required)
- `rhs` : tagged (nodes) (required)

### `call`

- `kind` : enum {call} (required)
- `fn` : text (required)
- `args` : taggedList (nodes) (required)

### `index`

- `kind` : enum {index} (required)
- `base` : tagged (nodes) (required)
- `args` : taggedList (nodes) (required)


## Model statements

### `assign`

- `type` : enum {assign} (required)
- `target` : text (required)
- `primed` : scalar (required)
- `expr` : tagged (nodes) (required)

### `interpCall`

- `type` : enum {interpCall} (required)
- `targets` : list of text (required)
- `primed` : scalar (required)
- `args` : taggedList (nodes) (required)
- `interpRef` : text (required)

### `reduction`

- `type` : enum {reduction} (required)
- `kind` : enum {EXPECT, MIN, MAX, PROD} (required)
- `target` : text (required)
- `body` : tagged (nodes) (required)
- `transRef` : ref -> pool:transitions (required)


## Equation kinds

### `plain`

- `kind` : enum {plain} (required)
- `expr` : tagged (nodes) (required)
- `primed` : scalar (required)

### `conditional`

- `kind` : enum {conditional} (required)
- `cases` : list of struct (required)
  - *(each item)*
    - `cond` : text (required)
    - `equations` : taggedList (eqkinds) (required)
