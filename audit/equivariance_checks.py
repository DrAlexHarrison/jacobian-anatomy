#!/usr/bin/env /usr/bin/python3
"""From-scratch verification of the C*-equivariance claim-set for
Alpoge's Jacobian counterexample F: C^3 -> C^3.

Written independently, from the claim statements alone. Deliberately does
NOT import, read, or reuse verify/03_new_facts.py (the artifact under audit).
All checks are exact symbolic computations over Q; no floating point anywhere.

Claims under test (from the orchestrator relay):
  C1. F(t x, y/t, z/t^2) = (F1/t^2, F2/t, t F3)  identically in t.
  C2. The source weight vector (1,-1,-2) is FORCED (unique up to scaling) by
      the requirement that each component of F be quasi-homogeneous; induced
      target weights are (-2,-1,1).
  C3. Gamma(t) = (4/(27 t^2), 4/(3 t), t) is a single orbit of
      tau_t: (a,b,c) -> (a/t^2, b/t, t c), and Gamma is Zariski-CLOSED in C^3
      (no extra closure points from t -> 0 or t -> oo).
  C4. The strata are tau-invariant: the quartic D = 27a^2c^2 - 18abc + 16a
      + b^3c - b^2 (quasi-invariant), the surface S = 27ac^2 - 9bc + 8
      (genuinely invariant), Gamma itself, and the punctured a-axis.
  C5. The "1-parameter family of triple fibers" is the tau-orbit of
      (-1/4, 0, 0), i.e. the punctured a-axis; the special fiber has exactly
      3 distinct points.
  C6. Refutation probes of "image = C^3 \\ Gamma": exact Groebner (=Nullstellensatz)
      certificates that rational points ON Gamma have EMPTY fiber, and that
      sample points OFF Gamma (generic and on V(D)) have NONEMPTY fiber.
Additional checks:
  B1. The x-eliminant relation D(F) x^3 + (4 - 3 F2 F3) x - 2 F3 == 0 holds
      identically (depressed cubic: no x^2 term).
  B2. disc_x(D x^3 + (4-3bc) x - 2c) factors as -4 D S^2 -- so the cubic's
      discriminant class in C(a,b,c)^* / squares is exactly -4D.
  B3. D restricted to {bc = 4/3} is a perfect square vanishing exactly on
      Gamma:  V(D) intersect {bc=4/3} = Gamma.
"""
import itertools
import sys

import sympy as sp
from sympy import Rational as R

x, y, z, t, s = sp.symbols('x y z t s')
a, b, c = sp.symbols('a b c')

F1 = (1 + x*y)**3 * z + y**2 * (1 + x*y) * (4 + 3*x*y)
F2 = y + 3*x*(1 + x*y)**2 * z + 3*x*y**2 * (4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
FF = (F1, F2, F3)

D = 27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2
S = 27*a*c**2 - 9*b*c + 8

failures = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"{tag} {name}" + (f"  [{detail}]" if detail else ""))
    if not ok:
        failures.append(name)


# ---------------------------------------------------------------- sanity
check("S0 det J == -2 identically",
      sp.expand(sp.Matrix([F1, F2, F3]).jacobian([x, y, z]).det()) == -2)

# ---------------------------------------------------------------- C1
sub = {x: t*x, y: y/t, z: z/t**2}
check("C1 F1(sigma_t) == t^-2 F1", sp.expand(t**2 * F1.subs(sub) - F1) == 0)
check("C1 F2(sigma_t) == t^-1 F2", sp.expand(t * F2.subs(sub) - F2) == 0)
check("C1 F3(sigma_t) == t    F3", sp.expand(F3.subs(sub) - t * F3) == 0)

# ---------------------------------------------------------------- C2
# Quasi-homogeneity of each component w.r.t. unknown weights (p,q,r) forces a
# linear system: within each component all monomial exponent vectors must have
# equal weight. Nullspace must be 1-dimensional and spanned by (1,-1,-2).
rows = []
for Fm in FF:
    exps = [sp.Matrix(e) for e, _ in sp.Poly(Fm, x, y, z).terms()]
    e0 = exps[0]
    rows += [(e - e0).T for e in exps[1:]]
M = sp.Matrix.vstack(*rows)
ns = M.nullspace()
w = None
if len(ns) == 1:
    v = ns[0]
    # normalize so first entry is 1
    w = sp.simplify(v / v[0])
check("C2 weight vector unique up to scaling (nullspace dim == 1)",
      len(ns) == 1, f"nullspace dim = {len(ns)}")
check("C2 normalized weight vector == (1,-1,-2)",
      w is not None and list(w) == [1, -1, -2], f"w = {list(w) if w is not None else None}")
wt = (1, -1, -2)
dvec = []
for Fm in FF:
    e = sp.Poly(Fm, x, y, z).terms()[0][0]
    dvec.append(sum(ei * wi for ei, wi in zip(e, wt)))
check("C2 induced target weights == (-2,-1,1)", dvec == [-2, -1, 1], f"d = {dvec}")

# ---------------------------------------------------------------- C3
def Gam(u):
    return (R(4, 27) / u**2, R(4, 3) / u, u)

def tau(P, u):
    return (P[0] / u**2, P[1] / u, u * P[2])

ts = sp.Symbol('ts')
orb = tau(Gam(t), s)
tgt = Gam(s * t)
check("C3 tau_s(Gamma(t)) == Gamma(s t)  (Gamma is one orbit)",
      all(sp.cancel(orb[i] - tgt[i]) == 0 for i in range(3)))
base = tau(Gam(1), s)
check("C3 orbit of Gamma(1) sweeps Gamma exactly (c-coord = s => injective param)",
      all(sp.cancel(base[i] - Gam(s)[i]) == 0 for i in range(3)))
# stabilizer of Gamma(1) trivial: c-coordinate of tau_s(Gamma(1)) is s
check("C3 stabilizer trivial (tau_s fixes Gamma(1) only for s == 1)",
      sp.solve(sp.Eq(base[2], 1), s) == [1])

# Closedness: Gamma == V(bc - 4/3, 12a - b^2) as SETS.
g1 = b*c - R(4, 3)
g2 = 12*a - b**2
gsub = {a: Gam(t)[0], b: Gam(t)[1], c: Gam(t)[2]}
check("C3 Gamma subset V(bc-4/3, 12a-b^2)",
      sp.cancel(g1.subs(gsub)) == 0 and sp.cancel(g2.subs(gsub)) == 0)
sols = sp.solve([g1, g2], [a, b], dict=True)
ok_conv = (len(sols) == 1
           and sp.cancel(sols[0][b] - R(4, 3)/c) == 0
           and sp.cancel(sols[0][a] - R(4, 27)/c**2) == 0)
check("C3 V(bc-4/3, 12a-b^2) subset Gamma (unique solution == Gamma(c))",
      ok_conv, f"solutions: {sols}")
check("C3 no c=0 points in V (t->0 adds no closure point in C^3)",
      sp.solve([g1.subs(c, 0), g2.subs(c, 0)], [a, b]) == [])
# t -> oo: c-coordinate = t -> oo, leaves every bounded set; algebraically the
# variety V above IS the Zariski closure and equals Gamma, so Gamma is closed.

# ---------------------------------------------------------------- C4
tsub = {a: a/t**2, b: b/t, c: t*c}
check("C4 D quasi-invariant: D(tau_t) == t^-2 D  => V(D) invariant",
      sp.expand(t**2 * D.subs(tsub) - D) == 0)
check("C4 S genuinely invariant: S(tau_t) == S",
      sp.expand(S.subs(tsub) - S) == 0)
check("C4 Gamma subset V(D)", sp.cancel(D.subs(gsub)) == 0)
check("C4 Gamma subset V(S)", sp.cancel(S.subs(gsub)) == 0)
# punctured a-axis: tau_t(a,0,0) = (a/t^2, 0, 0) -- invariant by inspection.
check("C4 a-axis invariant", tau((a, 0, 0), t)[1] == 0 and tau((a, 0, 0), t)[2] == 0)

# ---------------------------------------------------------------- C5
axis = tau((R(-1, 4), 0, 0), t)
check("C5 tau-orbit of (-1/4,0,0) lies in punctured a-axis",
      axis[1] == 0 and axis[2] == 0)
check("C5 orbit sweeps ALL of a-axis minus origin (solve -1/(4t^2) = s)",
      len(sp.solve(sp.Eq(R(-1, 4)/t**2, s), t)) > 0)

# ---------------------------------------------------------------- fibers
def fiber_data(pt):
    """Exact fiber analysis over QQ. Returns (count_with_multiplicity, sols).
    Since det J == -2 != 0 everywhere, every fiber is reduced (etale), so the
    quotient dimension IS the number of distinct C-points."""
    polys = [sp.expand(F1 - pt[0]), sp.expand(F2 - pt[1]), sp.expand(F3 - pt[2])]
    G = sp.groebner(polys, x, y, z, order='grevlex')
    exprs = list(G.exprs)
    if exprs == [1]:
        return 0, []
    lms = [tuple(g.monoms(order='grevlex')[0]) for g in G.polys]
    # zero-dimensionality: a pure power of each variable must lead some element
    bounds = []
    for v in range(3):
        pures = [m[v] for m in lms if all(m[u] == 0 for u in range(3) if u != v)]
        if not pures:
            return None, None  # not zero-dimensional -- would itself be a finding
        bounds.append(min(pures))
    cnt = 0
    for e in itertools.product(*(range(bd) for bd in bounds)):
        if not any(all(e[i] >= m[i] for i in range(3)) for m in lms):
            cnt += 1
    try:
        sols = sp.solve(polys, [x, y, z], dict=True)
    except Exception:
        sols = None
    return cnt, sols

# special fiber
cnt, sols = fiber_data((R(-1, 4), 0, 0))
check("C5 fiber over (-1/4,0,0): quotient dim == 3 (reduced => 3 distinct pts)",
      cnt == 3, f"count = {cnt}")
ok_pts = (sols is not None and len(sols) == 3 and
          all(all(sp.simplify(Fi.subs(so) - v) == 0
                  for Fi, v in zip(FF, (R(-1, 4), 0, 0))) for so in sols))
check("C5 the 3 solutions verified to map to (-1/4,0,0)",
      ok_pts, f"sols = {sols}")

# another point of the SAME tau-orbit family (a=1 on the a-axis): still 3
cnt, _ = fiber_data((1, 0, 0))
check("C5 fiber over (1,0,0) (same family) has 3 points", cnt == 3, f"count = {cnt}")

# ---------------------------------------------------------------- C6
cnt, _ = fiber_data(Gam(1))
check("C6 fiber over Gamma(1)=(4/27,4/3,1) EMPTY (Groebner == [1], exact)",
      cnt == 0, f"count = {cnt}")
cnt, _ = fiber_data(Gam(2))
check("C6 fiber over Gamma(2)=(1/27,2/3,2) EMPTY (independent 2nd Gamma point)",
      cnt == 0, f"count = {cnt}")
cnt, _ = fiber_data(Gam(R(-1, 3)))
check("C6 fiber over Gamma(-1/3)=(4/3,-4,-1/3) EMPTY (3rd Gamma point)",
      cnt == 0, f"count = {cnt}")

# OFF Gamma probes: generic (D != 0) expect 3; on V(D)\Gamma expect >= 1.
for pt, expect, label in [
        ((2, 1, 1), 3, "generic (2,1,1), D=104"),
        ((1, 1, 1), 3, "generic (1,1,1), D=25"),
        ((R(-16, 27), 0, 1), 1, "V(D)\\Gamma (-16/27,0,1)"),
        ((0, 0, 1), 1, "V(D)\\Gamma (0,0,1)"),
        ((0, 0, 0), 1, "origin, D=0"),
]:
    Dval = D.subs({a: pt[0], b: pt[1], c: pt[2]})
    cnt, _ = fiber_data(pt)
    check(f"C6 fiber over {label} has {expect} point(s)",
          cnt == expect, f"count = {cnt}, D = {Dval}")

# ---------------------------------------------------------------- B1
E = D*x**3 + (4 - 3*b*c)*x - 2*c
Esub = E.subs({a: F1, b: F2, c: F3})
check("B1 eliminant relation E(F1,F2,F3; x) == 0 identically",
      sp.expand(Esub) == 0)
check("B1 E is depressed (coeff of x^2 is 0)",
      sp.Poly(E, x).all_coeffs()[1] == 0)

# ---------------------------------------------------------------- B2
disc = sp.discriminant(E, x)
q = sp.cancel(disc / (-4*D))
check("B2 disc_x(E) == -4 * D * S^2", sp.expand(q - S**2) == 0,
      f"disc/( -4D ) - S^2 simplifies to {sp.simplify(q - S**2)}")

# ---------------------------------------------------------------- B3
Drest = sp.cancel(D.subs(b, R(4, 3)/c) * 27 * c**2)
check("B3 27c^2 * D|_{bc=4/3} == (27 a c^2 - 4)^2  => V(D) cap {bc=4/3} == Gamma",
      sp.expand(Drest - (27*a*c**2 - 4)**2) == 0)

# ---------------------------------------------------------------- summary
print()
if failures:
    print(f"RESULT: {len(failures)} FAILURE(S): {failures}")
    sys.exit(1)
print("RESULT: ALL CHECKS PASS")
