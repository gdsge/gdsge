import sympy
from gdsge_sympy.ast_to_sympy import to_sympy, symbol_name


def num(v):  return {"kind": "num", "value": v}
def name(i): return {"kind": "name", "id": i}
def primed(i): return {"kind": "primed", "id": i}
def binop(op, a, b): return {"kind": "binop", "op": op, "lhs": a, "rhs": b}
def unop(op, a): return {"kind": "unop", "op": op, "arg": a}
def call(fn, *a): return {"kind": "call", "fn": fn, "args": list(a)}


def test_basic_arithmetic():
    syms = {}
    e = to_sympy(binop("+", name("a"), binop("*", num(2), name("b"))), syms)
    a, b = sympy.symbols("a b", real=True)
    assert sympy.simplify(e - (a + 2 * b)) == 0


def test_power_maps_to_pow():
    syms = {}
    e = to_sympy(binop("^", name("x"), num(3)), syms)
    x = sympy.Symbol("x", real=True)
    assert sympy.simplify(e - x ** 3) == 0


def test_unary_minus():
    syms = {}
    e = to_sympy(unop("-", name("x")), syms)
    assert sympy.simplify(e + sympy.Symbol("x", real=True)) == 0


def test_builtin_call_exp_log():
    syms = {}
    e = to_sympy(call("exp", call("log", name("x"))), syms)
    # exp(log(x)) == x for the positive economic quantities we differentiate;
    # check numerically to avoid branch-cut assumptions.
    assert abs(float(e.subs(sympy.Symbol("x", real=True), 2.0)) - 2.0) < 1e-12


def test_primed_distinct_symbol():
    syms = {}
    e = to_sympy(binop("+", name("g"), primed("g")), syms)
    # name g and primed g are DISTINCT symbols (g vs g__prime)
    assert symbol_name(name("g")) != symbol_name(primed("g"))
    assert len(e.free_symbols) == 2
