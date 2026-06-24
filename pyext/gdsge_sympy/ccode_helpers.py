"""ccode emission tuned for the GDSGE C++ templates (hans recipe, generalized).
pow() is emitted for ** (matches emitExpr.m). Symbol names pass through verbatim
so they bind to the C++ locals MATLAB declares for them."""
from sympy.printing.c import C99CodePrinter


class _GdsgePrinter(C99CodePrinter):
    def _print_Pow(self, expr):
        if expr.exp == -1:
            return f"1.0/({self._print(expr.base)})"
        return f"pow({self._print(expr.base)},{self._print(expr.exp)})"


_PRINTER = _GdsgePrinter()


def ccode(expr):
    return _PRINTER.doprint(expr)
