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
    % Tolerance absorbs JSON serialization jitter (~1 ULP). Safe against false
    % equality because no IR scalar lies in (0, 1e-12); option values are >= 1e-6.
    tf = isequaln(a, b) || all(abs(a(:) - b(:)) <= 1e-12 + 1e-12 .* abs(b(:)));
else
    tf = isequaln(a, b);
end
end
