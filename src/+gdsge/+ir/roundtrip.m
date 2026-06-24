function ir2 = roundtrip(ir)
% ROUNDTRIP  decode(encode(ir)); the workhorse of the serialization tests.
ir2 = gdsge.ir.decode(gdsge.ir.encode(ir));
end
