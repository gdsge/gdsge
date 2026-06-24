function kids = children(node)
% CHILDREN  Return a cell array of an AST node's child nodes (empty for leaves).
switch node.kind
    case {'num','name','primed'}
        kids = {};
    case 'unop'
        kids = {node.arg};
    case 'binop'
        kids = {node.lhs, node.rhs};
    case 'call'
        kids = reshape(node.args, 1, []);
    case 'index'
        kids = [{node.base}, reshape(node.args, 1, [])];
    otherwise
        error('gdsge:ir:node:unknownKind', 'Unknown node kind "%s".', node.kind);
end
end
