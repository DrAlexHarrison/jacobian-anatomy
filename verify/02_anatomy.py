#!/usr/bin/env python3
"""Geometric anatomy of the counterexample: eliminants, monodromy, image.

Verifies:
  1. x-eliminant: the cubic  D*x^3 + (4-3bc)*x - 2c  with
     D = 27a^2c^2 - 18abc + 16a + b^3c - b^2  (depressed: no x^2 term).
  2. y-eliminant: 2y^3 - 3by^2 + 18ay + (27a^2c - 18ab + b^3), leading
     coefficient CONSTANT (y never escapes: y is integral over C[a,b,c]).
  3. disc_x factors as const * (27ac^2 - 9bc + 8)^2 * D  => monodromy S3.
  4. Fiber stratification 3/1/0: generic 3; on V(D)\\Gamma exactly 1;
     on Gamma(t) = (4/(27t^2), 4/(3t), t) exactly 0  => image = C^3 minus Gamma.
"""
import itertools
import sympy as sp
from sympy.polys.orderings import monomial_key

x, y, z, a, b, c = sp.symbols('x y z a b c')
R = sp.Rational
F = [(1 + x*y)**3 * z + y**2 * (1 + x*y) * (4 + 3*x*y),
     y + 3*x*(1 + x*y)**2 * z + 3*x*y**2 * (4 + 3*x*y),
     2*x - 3*x**2*y - x**3*z]
D = 27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2

# F3 = c is linear in z: substitute and clear denominators.
zsub = (2*x - 3*x**2*y - c) / x**3
G1 = sp.expand(sp.together(F[0].subs(z, zsub) - a) * x**3)
G2 = sp.expand(sp.together(F[1].subs(z, zsub) - b) * x**2)

# --- 1. x-eliminant -------------------------------------------------------
res_x = sp.Poly(G1, y).resultant(sp.Poly(G2, y))
res_x = res_x.as_expr() if isinstance(res_x, sp.Poly) else res_x
Px = sp.factor_list(sp.expand(res_x))
cubic = None
for f, m in Px[1]:
    if sp.degree(f, x) == 3:
        cubic = sp.expand(f)
assert cubic is not None
expected = sp.expand(D*x**3 + (4 - 3*b*c)*x - 2*c)
ratio = sp.simplify(cubic / expected)
assert ratio.is_constant() if hasattr(ratio, 'is_constant') else ratio.is_number, cubic
assert sp.simplify(cubic - ratio*expected) == 0
print("PASS 1 : x-eliminant = D*x^3 + (4-3bc)*x - 2c (depressed cubic)")

# --- 2. y-eliminant -------------------------------------------------------
res_y = sp.Poly(G1, x).resultant(sp.Poly(G2, x))
res_y = res_y.as_expr() if isinstance(res_y, sp.Poly) else res_y
Py = sp.factor_list(sp.expand(res_y))
ycub = None
for f, m in Py[1]:
    if sp.degree(f, y) == 3:
        ycub = sp.expand(f)
assert ycub is not None
lc = sp.Poly(ycub, y).LC()
assert lc.is_number, lc                      # constant leading coefficient
exp_y = sp.expand(2*y**3 - 3*b*y**2 + 18*a*y + 27*a**2*c - 18*a*b + b**3)
r2 = sp.simplify(ycub / exp_y)
assert r2.is_number and sp.simplify(ycub - r2*exp_y) == 0
print("PASS 2 : y-eliminant = 2y^3 - 3by^2 + 18ay + (27a^2c - 18ab + b^3); LC constant")

# --- 3. discriminant factorization => S3 ----------------------------------
disc = sp.factor_list(sp.discriminant(sp.Poly(expected, x)))
fs = {sp.expand(f): m for f, m in disc[1] if f.free_symbols}
even = sp.expand(27*a*c**2 - 9*b*c + 8)
assert fs.get(even) == 2 or fs.get(sp.expand(-even)) == 2, fs
assert fs.get(sp.expand(D)) == 1 or fs.get(sp.expand(-D)) == 1, fs
print("PASS 3 : disc_x = const*(27ac^2-9bc+8)^2 * D  (D to ODD power => S3 monodromy)")

# --- 4. fiber stratification 3/1/0 ---------------------------------------
key = monomial_key('grevlex')
def fiber_dim(t):
    eqs = [sp.expand(F[i] - t[i]) for i in range(3)]
    Gb = sp.groebner(eqs, x, y, z, order='grevlex')
    if list(Gb.exprs) == [sp.S.One]:
        return 0
    lm = [max(sp.Poly(g, x, y, z).monoms(), key=key) for g in Gb.exprs]
    bound = [None]*3
    for e in lm:
        for i in range(3):
            if e[i] > 0 and all(e[k] == 0 for k in range(3) if k != i):
                bound[i] = e[i] if bound[i] is None else min(bound[i], e[i])
    assert None not in bound
    return sum(1 for m in itertools.product(*(range(bb) for bb in bound))
               if not any(all(m[i] >= e[i] for i in range(3)) for e in lm))

assert fiber_dim((R(3, 7), R(-2, 5), R(1, 3))) == 3          # generic
assert fiber_dim((R(-16, 27), 0, 1)) == 1                    # on V(D), off Gamma
assert D.subs({a: R(-16, 27), b: 0, c: 1}) == 0
for t0 in (1, 2, -3):                                        # on Gamma: empty
    gpt = (R(4, 27)/t0**2, R(4, 3)/t0, R(t0))
    assert sp.simplify(D.subs(dict(zip((a, b, c), gpt)))) == 0
    assert fiber_dim(gpt) == 0
sol = sp.solve([sp.expand(F[i] - t) for i, t in
                zip(range(3), (R(-16, 27), 0, 1))], [x, y, z], dict=True)
assert sol == [{x: R(1, 2), y: R(-8, 3), z: 16}]
print("PASS 4 : fibers 3 (generic) / 1 (V(D)\\Gamma, witness (1/2,-8/3,16)) / 0 (Gamma)")

print("ALL ANATOMY CERTIFICATES PASS")
