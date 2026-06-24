"""SymPy analytic-Jacobian code generation for the GDSGE alt backend.

Entry points take a single JSON string and return a single JSON string, so the
MATLAB<->Python boundary is one char in / one char out (see callSympy.m).
"""
import json


def ping(payload: str) -> str:
    """Round-trip sanity check used by the bridge test."""
    obj = json.loads(payload)
    return json.dumps({"ok": True, "echo": obj})


from .diff_body import diff_body_json  # noqa: E402,F401
