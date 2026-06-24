function n = binop(op, lhs, rhs)
% BINOP  Binary-operator AST node.
n.kind = 'binop';
n.op = op;
n.lhs = lhs;
n.rhs = rhs;
end
