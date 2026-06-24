function [model, iterCode, cppCache, cppCode, codeSegment] = gdsge_codegen(modelName, options)
% GDSGE_CODEGEN  Backward-compatible flat entry point (frozen public API).
%   Thin shim over gdsge.codegen.codegen: parses <modelName>.gmod from the
%   current directory, writes the generated files there, and compiles the MEX
%   when its C++ changed. Returns the model IR struct. The legacy positional
%   outputs are re-exposed for the gdsge.m orchestrator:
%     iterCode    - text of the generated iter_<model>.m
%     cppCache    - mex_<model>.cache content before this run ('' if none)
%     cppCode     - text of the generated mex_<model>.cpp
%     codeSegment - struct of the real generated strings (iterCode/simulateCode/
%                   cppCode/compileCode); NOT the old template fragments, which
%                   have no counterpart in the thin-file architecture.
if nargin < 2; options = []; end
if nargout <= 1
    model = gdsge.codegen.codegen(modelName, options);
else
    [model, gen] = gdsge.codegen.codegen(modelName, options);
    iterCode    = gen.iterCode;
    cppCache    = gen.cppCache;
    cppCode     = gen.cppCode;
    codeSegment = gen.codeSegment;
end
end
