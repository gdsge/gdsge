% GDSGE  Refactored Global DSGE toolbox (package namespace for internal modules).
%
% Subpackages:
%   +parser   - gmod -> IR front end (preprocess, lex, declarations, model expressions)
%   +ir       - IR struct definition, JSON encode/decode, validators
%   +codegen  - IR -> MATLAB and IR -> C++ generators
%   +runtime  - helpers used by generated MATLAB (error reporting, printing)
%
% Public flat entry points (gdsge.m, gdsge_codegen.m) are thin shims added later;
% they are kept at top level for backward compatibility.
