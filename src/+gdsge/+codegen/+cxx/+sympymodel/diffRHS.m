function d = diffRHS(node, freeNames, helperPrefix)
% DIFFRHS  Ask Python (diff_body) for the value + partials of an expression
%   w.r.t. every currently-registered symbol name. Returns a normalized struct:
%     .value       char  (C++ expr for the value)
%     .helpers     struct array with .lhs/.rhs (CSE temporaries; possibly empty)
%     .partialsMap containers.Map symbolName -> C++ partial expr (nonzero only)
%     .freeSymbols cellstr of all symbol names in the body (incl. constants)
%   helperPrefix makes CSE temporary names unique per call so function-scope
%   statements don't collide on helper_0. Symbols NOT in freeNames (shock
%   realizations, params) are constants and never appear in partials.
if nargin < 3; helperPrefix = 'helper_'; end
req = struct('body', node, 'diffVars', {reshape(freeNames, 1, [])}, ...
    'helperPrefix', helperPrefix);
raw = gdsge.codegen.sympy.callSympy('diff_body_json', req);
d.value = raw.value;
d.helpers = normalizeHelpers(raw);
d.partialsMap = partialsToMap(raw);
d.freeSymbols = normalizeStrList(raw, 'freeSymbols');
end

function L = normalizeStrList(raw, fld)
% jsondecode: ["a"] -> 1x1 cell? no -> char row; ["a","b"] -> 2x1 cell of char;
% [] -> []. Normalize to a cellstr row.
L = {};
if ~isfield(raw, fld) || isempty(raw.(fld)); return; end
v = raw.(fld);
if ischar(v); L = {v}; elseif iscell(v); L = reshape(v, 1, []); end
end

function H = normalizeHelpers(raw)
% jsondecode: [] (empty list) -> [] ; [{...}] -> 1x1 struct ; [{...},{...}] ->
% Nx1 struct array. Normalize the empty case to an empty struct array.
if ~isfield(raw, 'helpers') || isempty(raw.helpers)
    H = struct('lhs', {}, 'rhs', {});
else
    H = raw.helpers;
end
end

function m = partialsToMap(raw)
m = containers.Map('KeyType', 'char', 'ValueType', 'char');
if isfield(raw, 'partials') && isstruct(raw.partials)
    fn = fieldnames(raw.partials);
    for i = 1:numel(fn)
        m(fn{i}) = raw.partials.(fn{i});
    end
end
end
