% +PARSER  GDSGE gmod front-end: text -> validated IR.
%   parseFrontEnd   - top-level: gmod text + model name -> validated full IR
%   preprocess      - strip comments, join continuations, rewrite deprecated keywords
%   splitBlocks     - separate declaration lines from named blocks
%   splitStatements - bracket-aware split into logical statements (+ start lines)
%   defaultSetupCode- parser-owned default flag values
%   evalSetup       - eval a setup script in an isolated workspace; capture vars
%   parseVarDecls   - structural parse of the declaration region (no eval)
%   parseDeclarations - eval + assemble params/shocks/states/variables/bounds/interp
%   parseSimulate   - the simulate block -> IR simulate section
%   resolveOptions  - flag workspace -> curated IR options
%   tokenize        - model-body text -> token stream (name/number/op/punct/prime)
%   parseExpr       - precedence-climbing expression parser -> AST nodes
%   parseModel      - model block -> typed statements + equations (reduction hoisting)
%   analyzeModel    - semantic checks: names, square system, interp arity, transRefs
%   assemblePartialIR - combine sections (incl. model) into a validated IR
