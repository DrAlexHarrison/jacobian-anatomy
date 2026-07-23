#!/usr/bin/env /usr/bin/python3
"""From-scratch verification of the sheet-parametrization claims
(addendum to equivariance_checks.py; same independence rules -- no reuse of
verify/03_new_facts.py or verify/04_parametrization.py).

Claims under test:
  P0. Covering identity: with r = x/(1+xy),  2F1 r^3 - F2 r^2 + 2r - F3 == 0
      identically; equivalently C_q(r) = 2a r^3 - b r^2 + 2r - c vanishes at
      the r-coordinate of every fiber point over q = (a,b,c).
  P1. WALL IDENTITY: with y(r) = -(b r^2 + 3c - 6r)/(2 r^2),
      1 - r*y(r) = C'(r)/2 on the cubic. We verify the stronger EXACT form
      with explicit cofactor: (1 - r*y(r)) - C'(r)/2 == -3 C(r) / (2r).
  P2. SECTION: x(r) = r/(1 - r*y(r)), z(r) = (2x - 3x^2 y - c)/x^3 solved from
      F3 = c; then F(x(r),y(r),z(r)) == (a,b,c) identically mod C(r) in
      Q(a,b,c)[r].
  P3. r = x/(1+xy) and x = r/(1-ry) are mutually inverse substitutions.
  P4. disc_r(C) == -4 D  (D the branch quartic).
  P5. Over Gamma(t), C is a perfect cube (8/(27 t^2)) (r - 3t/2)^3; and the
      triple-root locus of the family is exactly Gamma (coefficient matching).
  P6. LIFTING LAW "#fiber = #simple roots of C": tested at a
      battery of exact rational points, including the degenerate strata a = 0
      (degree drop; sheet at r = infinity = the surface {1+xy=0}) and c = 0
      (r = 0 root = the x = 0 stratum). We compare Groebner fiber counts vs
      simple-root counts in BOTH the affine and the projective (homogenized)
      reading, to pin down the exactly-correct statement.
  P7. x = 0 stratum: F(0,y,z) = (z+4y^2, y, 0), an isomorphism onto {c=0}.
  P8. Infinity sheet: F(x, -1/x, z) = (0, 2/x, 5x - x^3 z), covering {a=0,b!=0}.
"""
import itertools
import sys

import sympy as sp
from sympy import Rational as R

x, y, z, t, r, rho = sp.symbols('x y z t r rho')
a, b, c = sp.symbols('a b c')

F1 = (1 + x*y)**3 * z + y**2 * (1 + x*y) * (4 + 3*x*y)
F2 = y + 3*x*(1 + x*y)**2 * z + 3*x*y**2 * (4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z

C = 2*a*r**3 - b*r**2 + 2*r - c
Cp = sp.diff(C, r)
D = 27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2

failures = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"{tag} {name}" + (f"  [{detail}]" if detail else ""))
    if not ok:
        failures.append(name)


# ---------------------------------------------------------------- P0
rr = x / (1 + x*y)
cov = 2*F1*rr**3 - F2*rr**2 + 2*rr - F3
num0, _den0 = sp.fraction(sp.together(cov))
check("P0 covering identity 2F1 r^3 - F2 r^2 + 2r - F3 == 0 at r=x/(1+xy)",
      sp.expand(num0) == 0)

# ---------------------------------------------------------------- P1
yr = -(b*r**2 + 3*c - 6*r) / (2*r**2)
wall_lhs = 1 - r*yr
delta = sp.together(wall_lhs - Cp/2 - (-3*C/(2*r)))
check("P1 WALL IDENTITY exact form: (1 - r y(r)) - C'/2 == -3C/(2r)",
      sp.expand(sp.fraction(delta)[0]) == 0)
# hence on C(r)=0 (r != 0): 1 - r y(r) = C'(r)/2, so the section's x = r/(1-ry)
# blows up exactly at multiple roots (C = C' = 0). Also verify 1 - r*y(r)
# simplifies to the claimed wall numerator:
check("P1 1 - r y(r) == (b r^2 - 4r + 3c)/(2r)",
      sp.cancel(wall_lhs - (b*r**2 - 4*r + 3*c)/(2*r)) == 0)

# ---------------------------------------------------------------- P2
xr = sp.cancel(r / (1 - r*yr))
check("P2 x(r) simplifies to 2r^2/(b r^2 - 4r + 3c)",
      sp.cancel(xr - 2*r**2/(b*r**2 - 4*r + 3*c)) == 0)
zr = sp.cancel((2*xr - 3*xr**2*yr - c) / xr**3)
sec = {x: xr, y: yr, z: zr}

dom = 'QQ(a,b,c)'
Cpoly = sp.Poly(C, r, domain=dom)

def residue_mod_C(expr):
    """numerator of (expr) as a poly in r, reduced mod C over QQ(a,b,c)."""
    num, _den = sp.fraction(sp.together(sp.cancel(expr)))
    p = sp.Poly(sp.expand(num), r, domain=dom)
    _q, rem = sp.div(p, Cpoly)
    return sp.expand(num), rem

n3, _ = sp.fraction(sp.together(sp.cancel(F3.subs(sec) - c)))
check("P2 F3(section) == c EXACTLY (z solved from F3)", sp.expand(n3) == 0)

num1, rem1 = residue_mod_C(F1.subs(sec) - a)
check("P2 F1(section) - a: nonzero pre-reduction, == 0 mod C(r)",
      num1 != 0 and rem1.is_zero, f"rem = {rem1.as_expr()}")
num2, rem2 = residue_mod_C(F2.subs(sec) - b)
check("P2 F2(section) - b == 0 mod C(r)", rem2.is_zero,
      "in fact EXACTLY zero (stronger than claimed)" if num2 == 0
      else f"rem = {rem2.as_expr()}")
# structural note: y(r), z(r) make F2 = b and F3 = c hold EXACTLY; C(r) = 0 is
# precisely the remaining condition F1 = a. Verify that reading directly:
num1_only = sp.fraction(sp.together(sp.cancel(F1.subs(sec) - a)))[0]
check("P2 structure: F1(section) - a is a unit multiple of C(r) alone",
      sp.expand(num1_only) != 0 and rem1.is_zero and num2 == 0 and sp.expand(n3) == 0)

# ---------------------------------------------------------------- P3
r_of_x = x / (1 + x*y)                       # fiber-side coordinate
x_of_r = r / (1 - r*y)                       # section-side inverse (same y)
comp1 = sp.cancel(x_of_r.subs(r, r_of_x) - x)
comp2 = sp.cancel(r_of_x.subs(x, x_of_r) - r)
check("P3 x(r(x)) == x and r(x(r)) == r (mutually inverse)",
      comp1 == 0 and comp2 == 0)

# ---------------------------------------------------------------- P4
disc = sp.discriminant(C, r)
check("P4 disc_r(C) == -4 D exactly", sp.expand(disc + 4*D) == 0,
      f"disc = {sp.expand(disc)}")

# ---------------------------------------------------------------- P5
Gam = lambda u: {a: R(4, 27)/u**2, b: R(4, 3)/u, c: u}
CGam = C.subs(Gam(t))
cube = R(8, 27)/t**2 * (r - 3*t/2)**3
check("P5 C over Gamma(t) == (8/(27 t^2)) (r - 3t/2)^3 (perfect cube)",
      sp.expand(sp.together(CGam - cube)) == 0 or
      sp.cancel(CGam - cube) == 0)
# converse: C = 2a (r - rho)^3 forces (a,b,c) on Gamma. Match coefficients:
# expand 2a(r-rho)^3 = 2a r^3 - 6a rho r^2 + 6a rho^2 r - 2a rho^3
sols = sp.solve([sp.Eq(b, 6*a*rho), sp.Eq(2, 6*a*rho**2), sp.Eq(c, 2*a*rho**3)],
                [a, b, c], dict=True)
ok5 = (len(sols) == 1
       and sp.cancel(sols[0][a] - R(4, 27)/( (2*rho/3) )**2 * 1) is not None)
# substitute t = c = 2 a rho^3 -> with a = 1/(3 rho^2): t = 2 rho/3
if len(sols) == 1:
    s0 = sols[0]
    tval = sp.Rational(2, 3)*rho
    ok5 = (sp.cancel(s0[a] - R(4, 27)/tval**2) == 0
           and sp.cancel(s0[b] - R(4, 3)/tval) == 0
           and sp.cancel(s0[c] - tval) == 0)
else:
    ok5 = False
check("P5 triple-root locus == Gamma (coefficient matching, t = 2 rho/3)",
      ok5, f"solve -> {sols}")

# ---------------------------------------------------------------- P6
def fiber_count(pt):
    """exact #points of F^{-1}(pt) over C: Groebner quotient dimension;
    reduced because det J == -2 everywhere (etale)."""
    polys = [sp.expand(F1 - pt[0]), sp.expand(F2 - pt[1]), sp.expand(F3 - pt[2])]
    G = sp.groebner(polys, x, y, z, order='grevlex')
    if list(G.exprs) == [1]:
        return 0
    lms = [tuple(g.monoms(order='grevlex')[0]) for g in G.polys]
    bounds = []
    for v in range(3):
        pures = [m[v] for m in lms if all(m[u] == 0 for u in range(3) if u != v)]
        if not pures:
            return None  # not zero-dimensional: would itself be a finding
        bounds.append(min(pures))
    cnt = 0
    for e in itertools.product(*(range(bd) for bd in bounds)):
        if not any(all(e[i] >= m[i] for i in range(3)) for m in lms):
            cnt += 1
    return cnt

def root_counts(pt):
    """(affine simple roots of C_q, projective simple roots incl. r=infinity)."""
    a0, b0, c0 = pt
    Cq = sp.Poly(2*a0*r**3 - b0*r**2 + 2*r - c0, r)
    rts = sp.roots(Cq)
    assert sum(rts.values()) == Cq.degree(), "roots() missed a root"
    n_aff = sum(1 for m in rts.values() if m == 1)
    inf_mult = 3 - Cq.degree()
    return n_aff, n_aff + (1 if inf_mult == 1 else 0)

battery = [
    ("generic (2,1,1)",            (2, 1, 1)),
    ("generic (1,1,1)",            (1, 1, 1)),
    ("announced fiber (-1/4,0,0)", (R(-1, 4), 0, 0)),
    ("a-axis (1,0,0)",             (1, 0, 0)),
    ("V(D) (-16/27,0,1)",          (R(-16, 27), 0, 1)),
    ("V(D), a=0 (0,0,1)",          (0, 0, 1)),
    ("origin (0,0,0)",             (0, 0, 0)),
    ("Gamma(1) (4/27,4/3,1)",      (R(4, 27), R(4, 3), 1)),
    ("Gamma(2) (1/27,2/3,2)",      (R(1, 27), R(2, 3), 2)),
    ("a=0,b!=0 (0,1,0)",           (0, 1, 0)),
    ("a=0,b!=0 (0,2,3)",           (0, 2, 3)),
    ("a=0,bc=1 (0,1,1)",           (0, 1, 1)),
    ("c=0 generic (5,7,0)",        (5, 7, 0)),
    ("c=0 on V(D) (1,4,0)",        (1, 4, 0)),
]
print()
print("P6 lifting-law battery: fiber count vs simple roots (affine | projective)")
affine_law_holds = True
projective_law_holds = True
for label, pt in battery:
    fc = fiber_count(pt)
    n_aff, n_proj = root_counts(pt)
    Dval = D.subs({a: pt[0], b: pt[1], c: pt[2]})
    ok_aff = (fc == n_aff)
    ok_proj = (fc == n_proj)
    affine_law_holds &= ok_aff
    projective_law_holds &= ok_proj
    print(f"  {label:28s} fiber={fc}  simple_aff={n_aff}  simple_proj={n_proj}"
          f"  D={Dval}  {'' if ok_proj else '<-- PROJECTIVE MISMATCH'}")
check("P6 lifting law, PROJECTIVE reading (simple roots in P^1), all 14 points",
      projective_law_holds)
print(f"NOTE P6: affine reading holds at all points: {affine_law_holds} "
      "(expected False -- fails exactly on a=0, b!=0 where the r=infinity "
      "sheet {1+xy=0} contributes)")

# explicit witnesses for the two boundary sheets over (0,1,0):
w_inf = {x: 2, y: R(-1, 2), z: R(5, 4)}     # 1+xy = 0 sheet
w_sec = {x: -2, y: 1, z: 2}                  # section at r = 2
w_zero = {x: 0, y: 1, z: -4}                 # x = 0 sheet (r = 0)
for name, w in [("infinity-sheet (2,-1/2,5/4)", w_inf),
                ("section r=2 (-2,1,2)", w_sec),
                ("x=0 sheet (0,1,-4)", w_zero)]:
    img = tuple(sp.simplify(Fi.subs(w)) for Fi in (F1, F2, F3))
    check(f"P6 witness {name} maps to (0,1,0)", img == (0, 1, 0), f"img={img}")

# ---------------------------------------------------------------- P7
img0 = tuple(sp.expand(Fi.subs(x, 0)) for Fi in (F1, F2, F3))
check("P7 F(0,y,z) == (z + 4y^2, y, 0)",
      img0 == (sp.expand(z + 4*y**2), y, 0), f"img={img0}")
# iso onto {c=0}: (a,b,0) <- (0, b, a - 4b^2), unique x=0 preimage
back = tuple(sp.expand(Fi.subs({x: 0, y: b, z: a - 4*b**2})) for Fi in (F1, F2, F3))
check("P7 (0, b, a-4b^2) maps to (a, b, 0) (explicit inverse on the stratum)",
      back == (sp.expand(a), b, 0), f"img={back}")

# ---------------------------------------------------------------- P8
img_inf = tuple(sp.cancel(Fi.subs(y, -1/x)) for Fi in (F1, F2, F3))
check("P8 F(x, -1/x, z) == (0, 2/x, 5x - x^3 z)",
      (sp.cancel(img_inf[0]) == 0
       and sp.cancel(img_inf[1] - 2/x) == 0
       and sp.expand(img_inf[2] - (5*x - x**3*z)) == 0),
      f"img={img_inf}")

# ---------------------------------------------------------------- P9
# UNIQUENESS: on any fiber point (x,y,z) with r = x/(1+xy) finite and nonzero,
# y is FORCED to equal the section value y(r) at b = F2, c = F3. Together with
# x = r/(1-ry) (P3) and z solved from F3 (P2), the fiber point IS section(r):
# at most one fiber point per affine nonzero root. This closes the lifting law
# into a complete proof (see report for the assembled argument).
y_forced = -(F2*rr**2 + 3*F3 - 6*rr) / (2*rr**2)
check("P9 UNIQUENESS identity: y == y(r)|_{b=F2,c=F3,r=x/(1+xy)} identically",
      sp.cancel(sp.together(y_forced - y)) == 0)

# ---------------------------------------------------------------- summary
print()
if failures:
    print(f"RESULT: {len(failures)} FAILURE(S): {failures}")
    sys.exit(1)
print("RESULT: ALL CHECKS PASS")
