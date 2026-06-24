function s = emitExpr(node, ctx)
% EMITEXPR  IR expression AST -> C++ text (adept-compatible; no Symbolic
%   Toolbox). ctx: shockNames, futureVars, inLoop. See tEmitExpr for the
%   printing rules; conventions follow the old generator (_GRID accessors,
%   _GDSGE_GRID arrays, pow()).
%
%   Node field API (from +gdsge/+ir/+node/):
%     num   : .kind='num',   .value
%     name  : .kind='name',  .id
%     primed: .kind='primed',.id
%     binop : .kind='binop', .op, .lhs, .rhs
%     unop  : .kind='unop',  .op, .arg
%     call  : .kind='call',  .fn,  .args (cell array of nodes)
s = emit(node, 0);

    function s = emit(n, parentPrec)
        switch n.kind
            case 'num'
                s = numLit(n.value);

            case 'name'
                if ctx.inLoop && ismember(n.id, ctx.futureVars)
                    s = [n.id '_GDSGE_GRID[GDSGE_iter-1]'];
                else
                    s = n.id;
                end

            case 'primed'
                if ~ctx.inLoop
                    error('gdsge:codegen:primedOutsideLoop', ...
                        'primed reference %s'' outside a shock loop', n.id);
                end
                if ismember(n.id, ctx.shockNames)
                    s = [n.id '_GRID(GDSGE_iter)'];
                else
                    s = [n.id '_GDSGE_GRID[GDSGE_iter-1]'];
                end

            case 'binop'
                if strcmp(n.op, '^')
                    % Power: always emit as pow(lhs,rhs); children at prec 0
                    s = ['pow(' emit(n.lhs, 0) ',' emit(n.rhs, 0) ')'];
                else
                    p = binopPrec(n.op);
                    % Right child of left-associative - and / needs higher
                    % threshold to force parens on equal-precedence rhs
                    rp = p;
                    if any(strcmp(n.op, {'-', '/'}))
                        rp = p + 1;
                    end
                    s = [emit(n.lhs, p) n.op emit(n.rhs, rp)];
                    if p < parentPrec
                        s = ['(' s ')'];
                    end
                end

            case 'unop'
                if ~strcmp(n.op, '-')
                    error('gdsge:codegen:badNode', 'unsupported unop %s', n.op);
                end
                % Operand printed at prec 99 so any compound rhs gets parens;
                % the whole -x is wrapped when it appears under a parent op.
                s = ['-' emit(n.arg, 99)];
                if parentPrec > 0
                    s = ['(' s ')'];
                end

            case 'call'
                % call node stores function name in .fn (not .id)
                if ismember(n.fn, ctx.futureVars)
                    % Future-var indexed call — general over the index
                    % expression; the constant case (Barro: Re_n(2)) is the
                    % typical use. Emits a direct POPNARRAY grid access using
                    % the 1-based accessor convention; valid in or out of the
                    % shock loop because the index is explicit.
                    %
                    % This branch deliberately shadows emitDeclare's
                    %   #define <name>(idx) <name>_GDSGE_GRID[int(idx)-1]
                    % macro: the macro accepts any arity at the C preprocessor
                    % level, whereas this branch enforces the single-argument
                    % guard that the macro lacks, catching multi-arg mistakes
                    % at IR emit time.
                    if numel(n.args) ~= 1
                        error('gdsge:codegen:badNode', ...
                            'future var %s indexed with %d args (want 1)', ...
                            n.fn, numel(n.args));
                    end
                    s = [n.fn '_GDSGE_GRID[int(' emit(n.args{1}, 0) ')-1]'];
                else
                    argStrs = cellfun(@(a) emit(a, 0), n.args, 'UniformOutput', false);
                    s = [n.fn '(' strjoin(argStrs, ',') ')'];
                end

            otherwise
                error('gdsge:codegen:badNode', 'unknown node kind: %s', n.kind);
        end
    end
end

% ---------------------------------------------------------------------------
function p = binopPrec(op)
switch op
    case {'+', '-'}, p = 1;
    case {'*', '/'}, p = 2;
    otherwise
        error('gdsge:codegen:badNode', 'unsupported binop: %s', op);
end
end

% ---------------------------------------------------------------------------
function s = numLit(v)
if ~isfinite(v)
    error('gdsge:codegen:badNode', 'non-finite numeric literal %g', v);
end
if v == floor(v)
    s = sprintf('%.1f', v);   % %.1f never uses exponent; always yields N.0
else
    s = sprintf('%.17g', v);
end
end
