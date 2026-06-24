% +IR  GDSGE intermediate-representation contract.
%
% Schema & docs:
%   schema       - the single declarative descriptor (source of truth)
%   gendoc       - render docs/ir-schema.md from the schema
%
% Validation:
%   validate     - check an IR struct against the schema (returns a report)
%
% Serialization (MATLAB <-> JSON, descriptor-driven, round-trip safe):
%   encode       - IR struct  -> pretty JSON text
%   decode       - JSON text   -> canonical IR struct (rejects incompatible major version)
%   canonicalize - normalize any IR to canonical MATLAB form (idempotent)
%   roundtrip    - decode(encode(ir))
%   isequalIR    - structural IR equality (orientation/array-kind tolerant)
%
% AST node constructors live in the +node subpackage.
