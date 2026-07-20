#!/usr/bin/env python3
"""Three presentations of the covering cubic, and rigidity of the map.

Verifies:
  1. slowform333's presentation: v = y + 1/x satisfies
         c v^3 - 2 v^2 + b v - 2a = 0
     on fibers, and v = 1/r exactly: this cubic is the reciprocal-root
     (coefficient-reversed) form of Punt's C(r) = 2a r^3 - b r^2 + 2r - c.
     The x-eliminant D x^3 + (4-3bc) x - 2c is the pushforward of C under
     the (non-Moebius) section transform x = r/(1 - r y(r)).
  2. Siphon's verbatim family (@NinzaK40245): for all r != 0,
         F(0, 0, -1/(4r^2)) = F(r, -3/(2r), 13/(2r^2))
                            = F(-r, 3/(2r), 13/(2r^2)) = (-1/(4r^2), 0, 0)
     -- exactly the lambda_r-orbit of the announced fiber.
  3. RIGIDITY: in the ansatz q(w) = t + s w (w = xy),
         F1 = u^3 z + y^2 u q,  F2 = y + 3x u^2 z + 3x y^2 q,  F3 unchanged,
     requiring det J to be constant forces (t, s) = (4, 3) uniquely, where
     det = -2.  The counterexample is an isolated point of its natural family.
"""
import sympy as sp

x, y, z, a, b, c, r, v, t, s = sp.symbols('x y z a b c r v t s')
R = sp.Rational
F1 = (1+x*y)**3*z + y**2*(1+x*y)*(4+3*x*y)
F2 = y + 3*x*(1+x*y)**2*z + 3*x*y**2*(4+3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
Ccub = 2*a*r**3 - b*r**2 + 2*r - c

# --- 1. reciprocal presentation ------------------------------------------
vv = y + 1/x
assert sp.simplify(sp.together(F3*vv**3 - 2*vv**2 + F2*vv - 2*F1)) == 0
assert sp.simplify(vv - (1+x*y)/x) == 0                      # v = u/x = 1/r
rev = sp.expand(-(Ccub.subs(r, 1/v)) * v**3)
assert sp.expand(rev - (c*v**3 - 2*v**2 + b*v - 2*a)) == 0
print("PASS 1 : v = y + 1/x = 1/r; slowform333's cubic is the reversed C(r)")

# --- 2. Siphon's family = the lambda-orbit --------------------------------
tgt = (-1/(4*r**2), 0, 0)
for p in [(0, 0, -1/(4*r**2)), (r, -3/(2*r), 13/(2*r**2)), (-r, 3/(2*r), 13/(2*r**2))]:
    img = [sp.simplify(f.subs([(x, p[0]), (y, p[1]), (z, p[2])])) for f in (F1, F2, F3)]
    assert sp.simplify(img[0] - tgt[0]) == 0 and img[1] == 0 and img[2] == 0
print("PASS 2 : Siphon family verbatim == lambda_r-orbit of the announced fiber")

# --- 3. rigidity ----------------------------------------------------------
q = t + s*x*y
G1 = (1+x*y)**3*z + y**2*(1+x*y)*q
G2 = y + 3*x*(1+x*y)**2*z + 3*x*y**2*q
G3 = 2*x - 3*x**2*y - x**3*z
Jd = sp.expand(sp.Matrix([G1, G2, G3]).jacobian([x, y, z]).det())
poly = sp.Poly(Jd, x, y, z)
noncon = [cf for mono, cf in zip(poly.monoms(), poly.coeffs()) if mono != (0, 0, 0)]
assert noncon                                   # generic (t,s): det NOT constant
sols = sp.solve([sp.Eq(cf, 0) for cf in noncon], [t, s], dict=True)
assert sols == [{s: 3, t: 4}], sols             # unique solution
assert sp.expand(Jd.subs({t: 4, s: 3})) == -2
print("PASS 3 : constant det in the q = t + s*xy ansatz forces (t,s) = (4,3); det = -2")

print("ALL PRESENTATION/RIGIDITY CERTIFICATES PASS")
