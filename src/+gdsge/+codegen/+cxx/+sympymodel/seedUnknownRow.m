function row = seedUnknownRow(slotStart0, len, k)
% SEEDUNKNOWNROW  Gradient row of the k-th (0-based) element of an unknown
%   occupying [slotStart0 .. slotStart0+len-1] in GDSGE_x: unit at its own slot.
assert(k >= 0 && k < len);
row = gdsge.codegen.cxx.sympymodel.mkEntry(slotStart0 + k, '1', false);
end
