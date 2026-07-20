#!/usr/bin/env python3
"""Basic certificates for Alpöge's Jacobian-conjecture counterexample.

Verifies, with exact arithmetic and a sympy-independent cross-check:
  1. det Jacobian(F) == -2 identically (sympy symbolic AND hand-rolled
     Fraction-based polynomial engine AND numeric finite differences).
  2. The three announced points are distinct and share the image (-1/4,0,0).
  3. That fiber is COMPLETE (quotient dimension 3) and generic fibers also
     have exactly 3 points: F is generically 3-to-1.
  4. The induced Weyl-algebra endomorphism is well-defined: A = J^{-1} has
     polynomial entries and the derivations D_j = sum_k A[k,j] d_k commute.
"""
import itertools
from fractions import Fraction as Q

import sympy as sp
from sympy.polys.orderings import monomial_key

x, y, z = sp.symbols('x y z')
V = (x, y, z)
F = [(1 + x*y)**3 * z + y**2 * (1 + x*y) * (4 + 3*x*y),
     y + 3*x*(1 + x*y)**2 * z + 3*x*y**2 * (4 + 3*x*y),
     2*x - 3*x**2*y - x**3*z]

# --- 1a. sympy symbolic determinant ---------------------------------------
J = sp.Matrix(F).jacobian(V)
assert sp.expand(J.det()) == -2
print("PASS 1a: sympy det J == -2")

# --- 1b. independent Fraction-based engine (no sympy) ---------------------
def pmul(a, b):
    r = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(i + j for i, j in zip(ka, kb))
            r[k] = r.get(k, Q(0)) + va * vb
    return {k: v for k, v in r.items() if v != 0}

def padd(*ps):
    r = {}
    for p in ps:
        for k, v in p.items():
            r[k] = r.get(k, Q(0)) + v
    return {k: v for k, v in r.items() if v != 0}

def pscale(a, c):
    return {k: v * Q(c) for k, v in a.items()}

def pdiff(a, i):
    r = {}
    for k, v in a.items():
        if k[i] > 0:
            kk = list(k); kk[i] -= 1
            r[tuple(kk)] = r.get(tuple(kk), Q(0)) + v * k[i]
    return {k: v for k, v in r.items() if v != 0}

X = {(1, 0, 0): Q(1)}; Y = {(0, 1, 0): Q(1)}; Z = {(0, 0, 1): Q(1)}
ONE = {(0, 0, 0): Q(1)}
U = padd(ONE, pmul(X, Y)); U2 = pmul(U, U); U3 = pmul(U2, U)
G = padd(pscale(ONE, 4), pscale(pmul(X, Y), 3))
FF = [padd(pmul(U3, Z), pmul(pmul(pmul(Y, Y), U), G)),
      padd(Y, pscale(pmul(pmul(X, U2), Z), 3), pscale(pmul(pmul(X, pmul(Y, Y)), G), 3)),
      padd(pscale(X, 2), pscale(pmul(pmul(X, X), Y), -3), pscale(pmul(pmul(pmul(X, X), X), Z), -1))]
Jm = [[pdiff(f, i) for i in range(3)] for f in FF]
det = padd(
    pmul(Jm[0][0], padd(pmul(Jm[1][1], Jm[2][2]), pscale(pmul(Jm[1][2], Jm[2][1]), -1))),
    pscale(pmul(Jm[0][1], padd(pmul(Jm[1][0], Jm[2][2]), pscale(pmul(Jm[1][2], Jm[2][0]), -1))), -1),
    pmul(Jm[0][2], padd(pmul(Jm[1][0], Jm[2][1]), pscale(pmul(Jm[1][1], Jm[2][0]), -1))))
assert det == {(0, 0, 0): Q(-2)}
print("PASS 1b: independent Fraction engine det J == -2")

# --- 1c. numeric finite differences at random complex points --------------
import random
def Fnum(p):
    a, b, c = p
    u = 1 + a*b
    return [(u**3)*c + b*b*u*(4 + 3*a*b),
            b + 3*a*(u**2)*c + 3*a*b*b*(4 + 3*a*b),
            2*a - 3*a*a*b - a**3*c]
random.seed(7)
for _ in range(3):
    p = [complex(random.uniform(-2, 2), random.uniform(-2, 2)) for _ in range(3)]
    h = 1e-6
    Jn = [[0]*3 for _ in range(3)]
    for j in range(3):
        pp, pm = list(p), list(p)
        pp[j] += h; pm[j] -= h
        fp, fm = Fnum(pp), Fnum(pm)
        for i in range(3):
            Jn[i][j] = (fp[i] - fm[i]) / (2*h)
    d = (Jn[0][0]*(Jn[1][1]*Jn[2][2] - Jn[1][2]*Jn[2][1])
         - Jn[0][1]*(Jn[1][0]*Jn[2][2] - Jn[1][2]*Jn[2][0])
         + Jn[0][2]*(Jn[1][0]*Jn[2][1] - Jn[1][1]*Jn[2][0]))
    assert abs(d + 2) < 1e-3, d
print("PASS 1c: numeric det == -2 at random complex points")

# --- 2. the collision -----------------------------------------------------
R = sp.Rational
pts = [(0, 0, R(-1, 4)), (1, R(-3, 2), R(13, 2)), (-1, R(3, 2), R(13, 2))]
assert len(set(pts)) == 3
for p in pts:
    img = tuple(f.subs(dict(zip(V, p))) for f in F)
    assert img == (R(-1, 4), 0, 0), (p, img)
print("PASS 2 : three distinct points -> (-1/4, 0, 0)")

# --- 3. fiber completeness and generic degree -----------------------------
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
    assert None not in bound, "fiber not zero-dimensional"
    return sum(1 for m in itertools.product(*(range(b) for b in bound))
               if not any(all(m[i] >= e[i] for i in range(3)) for e in lm))

assert fiber_dim((R(-1, 4), 0, 0)) == 3          # announced fiber is COMPLETE
for t in [(R(3, 7), R(-2, 5), R(1, 3)), (R(11, 4), R(5, 9), R(-7, 2)),
          (R(-8, 3), R(13, 11), R(1, 6))]:
    assert fiber_dim(t) == 3
print("PASS 3 : special fiber complete (3 pts); generic fibers = 3 pts")

# --- 4. Weyl-algebra endomorphism integrability ---------------------------
A = sp.simplify(J.inv())
assert all(e.is_polynomial(x, y, z) for e in A)
for i in range(3):
    for j in range(i + 1, 3):
        for l in range(3):
            expr = sum(A[k, i]*sp.diff(A[l, j], V[k]) - A[k, j]*sp.diff(A[l, i], V[k])
                       for k in range(3))
            assert sp.expand(expr) == 0
print("PASS 4 : J^{-1} polynomial; derivations D_j commute (Weyl endomorphism exists)")

print("ALL BASIC CERTIFICATES PASS")
