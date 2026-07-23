#!/usr/bin/env python3
"""
Spot-check of the jacobianfun uniform-in-n construction at n = 7
(d = 6), a value beyond the previously verified n = 3, 4, 5. Builds F_7 from their
section-3 recipe with seed p_6, then verifies SYMBOLICALLY:
  1. alpha, beta are polynomials in (v,t) (the lift conditions close);
  2. F_7 = (alpha/x^2, beta/x, x*gamma) is a polynomial map;
  3. det J F_7 == 1 identically (b = c = 1);
  4. the inverse equation R(w) = w*P - Q holds as a polynomial identity;
  5. at a random rational target the degree-7 inverse polynomial is
     squarefree (7 distinct roots) and gamma != 0 at each root
     (resultant check) => generic fiber degree 7.
"""
import sympy as sp

w_, v, t, x, y, z = sp.symbols("w v t x y z")
d = 6
c = sp.Integer(1)
b = sp.Integer(1)

# seed p_d, their section 5
p = sp.expand(2*w_ - 3*w_**2 + w_*(1 - w_)*(w_**(d - 2) - sp.Rational(6, d*(d + 1))))
assert p.subs(w_, 0) == 0
assert p.subs(w_, 1) == -1
assert sp.integrate(p, (w_, 0, 1)) == 0
assert sp.degree(p, w_) == d and sp.LC(p, w_) == -1
kappa = sp.diff(p, w_).subs(w_, 1) / c
assert kappa == -5 + sp.Rational(6, d*(d + 1)) and kappa != -2
a = -(1 + kappa) / (2 + kappa)

# q' = w p'/c, q(0) = 0
q = sp.integrate(w_*sp.diff(p, w_)/c, (w_, 0, w_))
assert q.subs(w_, 0) == 0 and sp.diff(q, w_).subs(w_, 0) == 0

u = 1 + v
gam = 1 + a*v + b*t
W = u*gam

# 1. polynomiality of the lift pieces
beta = sp.expand(c + sp.cancel(p.subs(w_, W)/gam))
alpha = sp.expand(u + sp.cancel(q.subs(w_, W)/gam**2))
assert not sp.denom(sp.together(beta)).free_symbols, "beta not polynomial"
assert not sp.denom(sp.together(alpha)).free_symbols, "alpha not polynomial"
print("1. alpha, beta polynomial in (v,t): PASS")

# 2. the 3D lift is polynomial
av = alpha.subs({v: x*y, t: x**2*z})
bv = beta.subs({v: x*y, t: x**2*z})
gv = gam.subs({v: x*y, t: x**2*z})
F1 = sp.cancel(sp.expand(av)/x**2)
F2 = sp.cancel(sp.expand(bv)/x)
F3 = sp.expand(x*gv)
assert all(not sp.denom(sp.together(sp.expand(Fi))).free_symbols for Fi in (F1, F2, F3)), \
    "F_7 not polynomial"
print("2. F_7 = (alpha/x^2, beta/x, x gamma) polynomial: PASS  degrees:",
      [sp.total_degree(sp.expand(Fi)) for Fi in (F1, F2, F3)])

# 3. Keller
J = sp.Matrix([[sp.diff(Fi, s) for s in (x, y, z)] for Fi in (F1, F2, F3)])
detJ = sp.expand(J.det())
assert detJ == b*c, f"det J F_7 = {detJ}"
print("3. det J F_7 == 1 identically: PASS")

# 4. inverse equation R(w) = w P - c Q as identity in (w, gamma)
g_ = sp.Symbol("g_")
P = c*g_ + p
Q = w_*g_ + q
R = sp.integrate(p, (w_, 0, w_))
assert sp.expand(w_*P - c*Q - R) == 0
print("4. inverse identity  int_0^w p = w*P - c*Q: PASS")

# 5. random rational target: 7 simple roots, gamma != 0 at each
from random import Random
rng = Random(20260720)
ok_targets = 0
for _ in range(3):
    P0, Q0 = sp.Rational(rng.randint(-9, 9), rng.randint(1, 5)), \
             sp.Rational(rng.randint(-9, 9), rng.randint(1, 5))
    inv = sp.expand(R - w_*P0 + c*Q0)
    assert sp.degree(inv, w_) == 7
    disc = sp.discriminant(inv, w_)
    gam_at_root = sp.resultant(inv, P0 - p, w_)   # 0 iff some root has gamma=0
    if disc != 0 and gam_at_root != 0:
        ok_targets += 1
assert ok_targets == 3
print("5. three random targets: 7 distinct roots, gamma != 0 at all: PASS")
print()
print("n=7 spot-check: all checks pass; uniform-in-n recipe verified beyond n=3,4,5")
