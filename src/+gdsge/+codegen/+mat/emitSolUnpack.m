function frag = emitSolUnpack(ir)
% EMITSOLUNPACK  Slot unpacking + reshape fragments from the IR slot layout.
%   .unpack       <name>=GDSGE_SOL(lo:hi,:);  (policy, decl order)
%   .unpackAux    <name>=GDSGE_AUX(lo:hi,:);
%   .reshapeShock <name>=reshape(<name>,shock_num,[]);  (policy then aux)
%   .reshapeSize  save-time reshape to GDSGE_SIZE / [len GDSGE_SIZE]
%                 (policy, then interp vars, then aux — old reference order)
wU = gdsge.codegen.codeWriter();
wA = gdsge.codegen.codeWriter();
wS = gdsge.codegen.codeWriter();
wZ = gdsge.codegen.codeWriter();
for i = 1:numel(ir.variables.policy)
    v = ir.variables.policy{i};
    wU.add('%s=GDSGE_SOL(%d:%d,:);', v.name, v.slot(1), v.slot(2));
    wS.add('%s=reshape(%s,shock_num,[]);', v.name, v.name);
    if v.length == 1
        wZ.add('%s=reshape(%s,GDSGE_SIZE);', v.name, v.name);
    else
        wZ.add('%s=reshape(%s,[%d GDSGE_SIZE]);', v.name, v.name, v.length);
    end
end
for i = 1:numel(ir.variables.interp)
    n = ir.variables.interp{i};
    wZ.add('%s=reshape(%s,GDSGE_SIZE);', n, n);
end
for i = 1:numel(ir.variables.aux)
    v = ir.variables.aux{i};
    wA.add('%s=GDSGE_AUX(%d:%d,:);', v.name, v.slot(1), v.slot(2));
    wS.add('%s=reshape(%s,shock_num,[]);', v.name, v.name);
    if v.length == 1
        wZ.add('%s=reshape(%s,GDSGE_SIZE);', v.name, v.name);
    else
        wZ.add('%s=reshape(%s,[%d GDSGE_SIZE]);', v.name, v.name, v.length);
    end
end
frag = struct('unpack', wU.str(), 'unpackAux', wA.str(), ...
    'reshapeShock', wS.str(), 'reshapeSize', wZ.str());
end
