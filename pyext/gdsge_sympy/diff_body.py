"""Differentiate one body expression w.r.t. a set of free symbols, CSE the value
and all partials together, and emit C++ (helpers + value + partials). Spec §5.3.
"""
import json
import sympy
from sympy.simplify.cse_main import cse
from sympy.utilities.iterables import numbered_symbols
from .ast_to_sympy import to_sympy
from .ccode_helpers import ccode


def diff_body(body_node, diff_var_names, helper_prefix="helper_"):
    syms = {}
    value_expr = to_sympy(body_node, syms)

    # Differentiation targets: only those that actually appear in the body.
    targets = [(nm, syms[nm]) for nm in diff_var_names if nm in syms]

    partial_exprs = []
    partial_keys = []
    for nm, s in targets:
        d = sympy.diff(value_expr, s)
        if d != 0:
            partial_keys.append(nm)
            partial_exprs.append(d)

    # One CSE over value + all nonzero partials so they share helpers. The
    # prefix makes helper names unique per call so function-scope statements
    # (scalar assigns, equations) don't collide on helper_0.
    all_exprs = [value_expr] + partial_exprs
    replacements, reduced = cse(all_exprs, symbols=numbered_symbols(helper_prefix))

    # helpers as a list of {lhs,rhs} objects (not [lhs,rhs] pairs): decodes to a
    # clean MATLAB struct array via jsondecode, avoiding cell/char-matrix ambiguity.
    helpers = [{"lhs": str(lhs), "rhs": ccode(rhs)} for lhs, rhs in replacements]
    value_cpp = ccode(reduced[0])
    partials = {nm: ccode(reduced[1 + i]) for i, nm in enumerate(partial_keys)}
    # All symbols appearing in the body (incl. constants like shock realizations
    # that are NOT differentiation targets). The MATLAB reduction emitter uses
    # this to know which per-shock locals it must load inside the loop.
    free_symbols = sorted(str(s) for s in value_expr.free_symbols)
    return {"helpers": helpers, "value": value_cpp, "partials": partials,
            "freeSymbols": free_symbols}


def diff_body_json(payload: str) -> str:
    req = json.loads(payload)
    return json.dumps(diff_body(req["body"], req["diffVars"],
                                req.get("helperPrefix", "helper_")))
