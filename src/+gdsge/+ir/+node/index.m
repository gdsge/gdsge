function n = index(base, args)
% INDEX  Indexing AST node. ARGS is a cell array of nodes.
n.kind = 'index';
n.base = base;
n.args = args;
end
