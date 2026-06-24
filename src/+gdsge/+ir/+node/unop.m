function n = unop(op, arg)
% UNOP  Unary-operator AST node.
n.kind = 'unop';
n.op = op;
n.arg = arg;
end
