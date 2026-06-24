function e = mkEntry(slot, expr, templated)
% MKENTRY  One gradient-row entry. slot is a 0-based base slot (numeric).
%   templated=false: the actual GDSGE_x slot is `slot` (fixed). templated=true:
%   the actual slot at shock iteration GDSGE_iter is `slot + (GDSGE_iter-1)`
%   (a per-shock array-unknown element); such entries are only valid inside a
%   shock loop. expr is the C++ partial value (may reference per-shock locals).
if nargin < 3; templated = false; end
e = struct('slot', slot, 'expr', expr, 'templated', templated);
end
