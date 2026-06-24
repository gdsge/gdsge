def test_sympy_diff_and_ccode():
    import sympy

    x = sympy.symbols("x", real=True)
    assert sympy.diff(x**2, x) == 2 * x
    # ccode is the codegen primitive the alt backend depends on
    assert "x" in sympy.ccode(x * x)
