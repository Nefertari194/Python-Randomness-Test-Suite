from __future__ import print_function

try:
    from scipy.special import gamma, gammainc, gammaincc
except Exception:
    try:
        import mpmath

        def gamma(a):
            return float(mpmath.gamma(a))

        def gammainc(a, x):
            return float(mpmath.gammainc(a, 0, x) / mpmath.gamma(a))

        def gammaincc(a, x):
            return float(mpmath.gammainc(a, x, mpmath.inf) / mpmath.gamma(a))
    except Exception:
        from math import gamma
        from main import gammainc, gammaincc

