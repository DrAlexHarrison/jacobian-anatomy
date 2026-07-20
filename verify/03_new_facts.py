#!/usr/bin/env python3
"""Equivariance, the Punt covering cubic, explicit escape curve, computed S3.

Verifies:
  1. C*-equivariance: F(tx, y/t, z/t^2) = (F1/t^2, F2/t, t*F3).
     Source weights (1,-1,-2), target weights (-2,-1,1).
  2. Weight bookkeeping: D is tau-homogeneous of weight -2; 27ac^2-9bc+8 is
     tau-INVARIANT; Gamma is a single C*-orbit: tau_s(Gamma(t)) = Gamma(st);
     the triple-collision fibers are the orbit of the announced one:
     F(lambda_s(p_i)) = tau_s(-1/4,0,0) = (-1/(4s^2),0,0).
  3. Punt covering: r = x/(1+xy) satisfies 2*F1*r^3 - F2*r^2 + 2r - F3 = 0
     identically; the covering cubic 2a*r^3 - b*r^2 + 2r - c is equivariant
     of weight 1. Its r-discriminant is computed and factored.
  4. Explicit non-properness witness: F(t, -1/t, 5/t^2) = (0, 2/t, 0), an
     escaping arc over the origin, which lies on V(D).
  5. S3 is COMPUTED, not inferred: at (a,b,c)=(3/7,-2/5,1/3) the x-cubic is
     irreducible over Q with non-square discriminant, so its Galois group is
     S3; the specialized group embeds in the generic group, hence generic S3.
"""
import sympy as sp

x, y, z, a, b, c, t, s, r = sp.symbols('x y z a b c t s r')
R = sp.Rational
F = [(1 + x*y)**3 * z + y**2 * (1 + x*y) * (4 + 3*x*y),
     y + 3*x*(1 + x*y)**2 * z + 3*x*y**2 * (4 + 3*x*y),
     2*x - 3*x**2*y - x**3*z]
D = 27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2

# --- 1. equivariance ------------------------------------------------------
sub = {x: t*x, y: y/t, z: z/t**2}
scaled = [sp.simplify(f.subs(sub)) for f in F]
assert sp.simplify(scaled[0] - F[0]/t**2) == 0
assert sp.simplify(scaled[1] - F[1]/t) == 0
assert sp.simplify(scaled[2] - t*F[2]) == 0
print("PASS 1 : F(tx, y/t, z/t^2) = (F1/t^2, F2/t, t*F3)")

# --- 2. weight bookkeeping ------------------------------------------------
tau = {a: a/s**2, b: b/s, c: s*c}
assert sp.simplify(D.subs(tau) - D/s**2) == 0
even = 27*a*c**2 - 9*b*c + 8
assert sp.simplify(even.subs(tau) - even) == 0
gam = lambda u: (R(4, 27)/u**2, R(4, 3)/u, u)
g_t = gam(t)
tau_g = (g_t[0]/s**2, g_t[1]/s, s*g_t[2])
g_st = gam(s*t)
assert all(sp.simplify(u - v) == 0 for u, v in zip(tau_g, g_st))
pts = [(0, 0, R(-1, 4)), (1, R(-3, 2), R(13, 2)), (-1, R(3, 2), R(13, 2))]
for p in pts:
    lam_p = {x: s*p[0], y: p[1]/s, z: p[2]/s**2}
    img = [sp.simplify(f.subs(x, lam_p[x]).subs(y, lam_p[y]).subs(z, lam_p[z]))
           for f in F]
    assert sp.simplify(img[0] + 1/(4*s**2)) == 0 and img[1] == 0 and img[2] == 0
print("PASS 2 : D wt -2; 27ac^2-9bc+8 invariant; Gamma one orbit; collision fibers = orbit of a-axis")

# --- 3. Punt covering cubic ----------------------------------------------
rr = x/(1 + x*y)
ident = 2*F[0]*rr**3 - F[1]*rr**2 + 2*rr - F[2]
assert sp.simplify(sp.together(ident)) == 0
cub = 2*a*r**3 - b*r**2 + 2*r - c
tau_r = cub.subs(tau).subs(r, s*r)
assert sp.simplify(tau_r - s*cub) == 0            # weight-1 equivariant
disc_r = sp.factor(sp.discriminant(sp.Poly(cub, r)))
print("PASS 3 : 2*F1*r^3 - F2*r^2 + 2r - F3 == 0 with r = x/(1+xy); covering cubic wt 1")
print("        disc_r(2a r^3 - b r^2 + 2r - c) =", disc_r)

# --- 4. explicit escape arc ----------------------------------------------
arc = {x: t, y: -1/t, z: 5/t**2}
img = [sp.simplify(f.subs(arc)) for f in F]
assert img[0] == 0 and sp.simplify(img[1] - 2/t) == 0 and img[2] == 0
assert D.subs({a: 0, b: 0, c: 0}) == 0            # origin lies on V(D)
print("PASS 4 : F(t, -1/t, 5/t^2) = (0, 2/t, 0) -> origin; origin on V(D)")

# --- 5. computed S3 -------------------------------------------------------
spec = {a: R(3, 7), b: R(-2, 5), c: R(1, 3)}
cub_x = sp.expand((D*x**3 + (4 - 3*b*c)*x - 2*c).subs(spec))
fl = sp.factor_list(cub_x)
assert len(fl[1]) == 1 and fl[1][0][1] == 1 and sp.degree(fl[1][0][0], x) == 3
d0 = sp.discriminant(sp.Poly(cub_x, x))
assert d0 != 0 and not sp.sqrt(d0).is_rational
print("PASS 5 : specialized cubic irreducible/Q, disc non-square => Galois S3;")
print("        specialization embeds in generic group => generic monodromy = S3")

print("ALL NEW-FACT CERTIFICATES PASS")
