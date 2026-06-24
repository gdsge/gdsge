function n = call(fn, args)
% CALL  Function-call AST node. ARGS is a cell array of nodes.
n.kind = 'call';
n.fn = fn;
n.args = args;
end
