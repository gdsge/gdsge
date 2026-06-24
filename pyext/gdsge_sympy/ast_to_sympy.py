"""IR expression node (dict) -> SymPy expression. Inverse of emitExpr.m.

Each `name`/`primed` node becomes a SymPy Symbol; the caller supplies a `syms`
dict that maps the canonical symbol name to a sympy.Symbol so the same name
reuses one symbol. Primed names get a distinct '__prime' suffix so a current
unknown `g` and its shock realization `g'` never collide.
"""
import sympy

# gmod built-ins -> sympy. Anything not here raises (caught upstream).
_BUILTINS = {
    "exp": sympy.exp, "log": sympy.log, "sqrt": sympy.sqrt,
    "abs": sympy.Abs, "sin": sympy.sin, "cos": sympy.cos, "tan": sympy.tan,
    "max": sympy.Max, "min": sympy.Min,
}


def symbol_name(node):
    if node["kind"] == "name":
        return node["id"]
    if node["kind"] == "primed":
        return node["id"] + "__prime"
    raise ValueError("symbol_name expects a name/primed node")


def _sym(node, syms):
    nm = symbol_name(node)
    if nm not in syms:
        syms[nm] = sympy.Symbol(nm, real=True)
    return syms[nm]


def to_sympy(node, syms):
    k = node["kind"]
    if k == "num":
        v = node["value"]
        return sympy.Integer(int(v)) if v == int(v) else sympy.Float(v)
    if k in ("name", "primed"):
        return _sym(node, syms)
    if k == "unop":
        if node["op"] != "-":
            raise ValueError(f"unsupported unop {node['op']}")
        return -to_sympy(node["arg"], syms)
    if k == "binop":
        lhs = to_sympy(node["lhs"], syms)
        rhs = to_sympy(node["rhs"], syms)
        op = node["op"]
        if op == "+": return lhs + rhs
        if op == "-": return lhs - rhs
        if op == "*": return lhs * rhs
        if op == "/": return lhs / rhs
        if op == "^": return lhs ** rhs
        raise ValueError(f"unsupported binop {op}")
    if k == "call":
        fn = node["fn"]
        if fn == "pow":
            a, b = (to_sympy(x, syms) for x in node["args"])
            return a ** b
        if fn in _BUILTINS:
            return _BUILTINS[fn](*[to_sympy(x, syms) for x in node["args"]])
        raise ValueError(f"unsupported call {fn}")
    raise ValueError(f"unknown node kind {k}")
