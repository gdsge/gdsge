import json
import math
import re
from gdsge_sympy import diff_body_json


def binop(op, a, b): return {"kind": "binop", "op": op, "lhs": a, "rhs": b}
def name(i): return {"kind": "name", "id": i}
def num(v): return {"kind": "num", "value": v}
def call(fn, *a): return {"kind": "call", "fn": fn, "args": list(a)}


def _fix(s):
    """Rewrite emitted C++ expr into something Python's eval can run."""
    return re.sub(r"\bpow\(", "__pow(", s)


def _scope(res, env):
    scope = dict(env)
    scope["__pow"] = lambda a, b: a ** b
    scope["exp"] = math.exp
    scope["log"] = math.log
    for h in res["helpers"]:
        scope[h["lhs"]] = eval(_fix(h["rhs"]), {"__builtins__": {}}, scope)
    return scope


def _ev(expr, scope):
    return eval(_fix(expr), {"__builtins__": {}}, scope)


def test_value_and_partial_simple():
    # body = a*b ; d/da = b ; d/db = a
    body = binop("*", name("a"), name("b"))
    res = json.loads(diff_body_json(json.dumps(
        {"body": body, "diffVars": ["a", "b"]})))
    scope = _scope(res, {"a": 3.0, "b": 5.0})
    assert abs(_ev(res["value"], scope) - 15.0) < 1e-12
    assert abs(_ev(res["partials"]["a"], scope) - 5.0) < 1e-12
    assert abs(_ev(res["partials"]["b"], scope) - 3.0) < 1e-12


def test_zero_partial_omitted():
    # body = a + 7 ; d/db = 0 -> omitted
    body = binop("+", name("a"), num(7))
    res = json.loads(diff_body_json(json.dumps(
        {"body": body, "diffVars": ["a", "b"]})))
    assert "b" not in res["partials"]
    assert "a" in res["partials"]


def test_free_symbols_includes_constants():
    # body = c * a^2 ; differentiate only w.r.t. a; c is a constant but must
    # still appear in freeSymbols so the caller can load it.
    body = binop("*", name("c"), binop("^", name("a"), num(2)))
    res = json.loads(diff_body_json(json.dumps(
        {"body": body, "diffVars": ["a"]})))
    assert set(res["freeSymbols"]) == {"a", "c"}
    assert "c" not in res["partials"]   # not a diff target
    assert "a" in res["partials"]


def test_power_partial():
    # body = a^3 ; d/da = 3*a^2
    body = binop("^", name("a"), num(3))
    res = json.loads(diff_body_json(json.dumps(
        {"body": body, "diffVars": ["a"]})))
    scope = _scope(res, {"a": 2.0})
    assert abs(_ev(res["partials"]["a"], scope) - 12.0) < 1e-9
