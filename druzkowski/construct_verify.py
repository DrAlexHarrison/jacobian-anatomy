#!/usr/bin/env /usr/bin/python3
"""
Stage A of the classical reduction pipeline, made explicit on Alpoge's
Jacobian-Conjecture counterexample F: C^3 -> C^3 (degree 7, det JF = -2).

Following Bass-Connell-Wright, "The Jacobian conjecture: reduction of degree
and formal expansion of the inverse", Bull. AMS 7 (1982) 287-330:
  - Prop (3.1): degree reduction to <= 3 by adjoining variable pairs,
        F' = G o F[2] o H,  H = (X, X_{n+1}+P, X_{n+2}+Q),
        G  = (..., X_i - X_{n+1} X_{n+2}, ...),
    for a chosen top-degree term aM = P*Q of component i.
  - Section 4 (Step 2, "doubling"): F' = (X + F_(2)(X) + Y, Y - F_(3)(X))
    on C^{2n}; equals G(1) o F[n] o H(1) with G(1) = (X+Y, Y),
    H(1) = (X, Y - F_(3)(X)).  This makes J(N) nilpotent, N = F' - id.
  - Section 4 (Step 3, homogenization): adjoin T,
        L = (X + T^2 Y + T F_(2)(X), Y - F_(3)(X), T),
    identity + cubic homogeneous H with J(H) nilpotent; Keller exactly.

Every step is an explicit composition  A o (previous x id) o B  with A, B
elementary (or affine-triangular) automorphisms of determinant 1, so
det J is preserved EXACTLY and collisions transport EXPLICITLY:
    new colliding inputs  = B^{-1}(old point, 0)
    new common image      = A(old image, 0).
The script rebuilds everything from scratch and asserts:
  (0) ground truth: det JF = -2, three points collide;
  (1) efficient exact SYMBOLIC proof that det J = 1 for the degree<=3 map
      (column-elimination replay of the gadgets down to the 3x3 base case);
  (2) det J_L = 1 exactly-symbolically via Schur/block elimination + the
      same replay at scaled argument TX;
  (3) J(H) nilpotent for the final map, via det(I + s*JH) == 1 with s a
      fresh symbol (same elimination machinery at argument sTX), which
      forces char(JH) = lambda^N and hence JH^N = 0 (Cayley-Hamilton);
  (4) H homogeneous cubic, linear part = identity, F(0) = 0;
  (5) exact rational spot-checks of det J = 1 after EVERY gadget step;
  (6) the three collision points carried through EVERY step, with equal
      images asserted after EVERY step.

Run:  /usr/bin/python3 construct_verify.py
"""

import os
HERE = os.path.dirname(os.path.abspath(__file__))
import sys
import time
import random
from fractions import Fraction

import sympy as sp
from sympy import Rational

T0 = time.time()


def log(msg):
    print(f"[{time.time()-T0:7.1f}s] {msg}")
    sys.stdout.flush()


# ----------------------------------------------------------------------
# Stage 0: base map and ground truth
# ----------------------------------------------------------------------
x, y, z = sp.symbols('v1 v2 v3')

F1 = (1 + x*y)**3 * z + y**2 * (1 + x*y) * (4 + 3*x*y)
F2 = y + 3*x*(1 + x*y)**2 * z + 3*x*y**2*(4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z

F_base = [sp.expand(F1), sp.expand(F2), sp.expand(F3)]
base_vars = [x, y, z]

P0 = [Rational(0), Rational(0), Rational(-1, 4)]
Q0 = [Rational(1), Rational(-3, 2), Rational(13, 2)]
R0 = [Rational(-1), Rational(3, 2), Rational(13, 2)]

J_base = sp.Matrix(F_base).jacobian(base_vars)
detJ = sp.expand(J_base.det())
assert detJ == -2, f"ground truth det failed: {detJ}"

imgs = []
for pt in (P0, Q0, R0):
    sub = dict(zip(base_vars, pt))
    imgs.append(tuple(sp.nsimplify(f.subs(sub)) for f in F_base))
assert imgs[0] == imgs[1] == imgs[2] == (Rational(-1, 4), 0, 0), imgs
assert [f.subs(dict(zip(base_vars, [0, 0, 0]))) for f in F_base] == [0, 0, 0]
log("Stage 0 OK: det JF = -2 (symbolic); three points collide to (-1/4,0,0); F(0)=0")

# ----------------------------------------------------------------------
# Stage 1: linear normalization  Ftil = L0^{-1} o F
#   L0 = JF(0); then Ftil = X + H, H order >= 2, det J Ftil = 1.
# ----------------------------------------------------------------------
L0 = J_base.subs({x: 0, y: 0, z: 0})
assert L0 == sp.Matrix([[0, 0, 1], [0, 1, 0], [2, 0, 0]]), L0
L0inv = L0.inv()

Fv = sp.Matrix(F_base)
Ftil = list(sp.expand(L0inv * Fv))
# linear part must be the identity now
Jtil0 = sp.Matrix(Ftil).jacobian(base_vars).subs({x: 0, y: 0, z: 0})
assert Jtil0 == sp.eye(3)
detJtil = sp.expand(sp.Matrix(Ftil).jacobian(base_vars).det())
assert detJtil == 1, detJtil

common_image = list(L0inv * sp.Matrix(imgs[0]))          # = (0, 0, -1/4)
assert common_image == [0, 0, Rational(-1, 4)]
log(f"Stage 1 OK: Ftil = L0^-1 o F has identity linear part, det J = 1 (symbolic); "
    f"common image now {tuple(common_image)}")

# ----------------------------------------------------------------------
# Stage 2: BCW Prop (3.1) degree reduction to <= 3
# ----------------------------------------------------------------------
# state
VARS = list(base_vars)                       # v1..vn
FMAP = list(Ftil)                            # current map, expanded
POINTS = [list(P0), list(Q0), list(R0)]      # colliding inputs, transported
IMAGE = list(common_image)                   # common image, transported
STEPS = []                                   # (i, P, Q, ia, ib) for replay proof

RNG = random.Random(20260720)


def rand_point(n):
    return [Rational(RNG.randint(-9, 9), RNG.randint(1, 4)) for _ in range(n)]


def keller_spot_check(tag, npts=2):
    """Exact rational det J at random points; must equal 1."""
    Jm = sp.Matrix(FMAP).jacobian(VARS)
    for _ in range(npts):
        pt = rand_point(len(VARS))
        sub = dict(zip(VARS, pt))
        d = Jm.subs(sub).det()
        assert d == 1, f"{tag}: det J = {d} at {pt}"


def collision_check(tag):
    im = None
    for pt in POINTS:
        sub = dict(zip(VARS, pt))
        v = tuple(sp.nsimplify(f.subs(sub)) for f in FMAP)
        if im is None:
            im = v
        else:
            assert v == im, f"{tag}: images differ"
    assert tuple(im) == tuple(IMAGE), f"{tag}: image mismatch {im} vs {IMAGE}"
    # distinctness
    assert POINTS[0][:3] != POINTS[1][:3] != POINTS[2][:3]


def h_terms():
    """All (i, coeff, monom) of H = FMAP - id, as Poly monomials over VARS."""
    out = []
    for i, f in enumerate(FMAP):
        h = sp.expand(f - VARS[i])
        if h == 0:
            continue
        p = sp.Poly(h, *VARS)
        for c, m in zip(p.coeffs(), p.monoms()):
            out.append((i, c, m))
    return out


def monom_expr(m):
    e = sp.Integer(1)
    for v, k in zip(VARS, m):
        if k:
            e *= v**k
    return e


def split_monom(m, target):
    """Split exponent vector m into (mP, mQ), deg mP = target."""
    mP = [0]*len(m)
    left = target
    for j, e in enumerate(m):
        take = min(e, left)
        mP[j] = take
        left -= take
        if left == 0:
            break
    assert left == 0
    mQ = [e - p for e, p in zip(m, mP)]
    return tuple(mP), tuple(mQ)


def apply_gadget(i, coeff, m):
    """BCW (3.1) elementary reduction of term coeff*monom(m) in component i."""
    d = sum(m)
    # optimal split target for deg P (cost analysis: 4->2, 5->2, 6->3, 7->3)
    target = 2 if d in (4, 5) else 3
    mP, mQ = split_monom(m, target)
    P = sp.expand(coeff * monom_expr(mP))    # coefficient carried by P
    Q = monom_expr(mQ)

    n = len(VARS)
    va = sp.Symbol(f'v{n+1}')
    vb = sp.Symbol(f'v{n+2}')

    # transport collision inputs FIRST (P, Q involve only old vars)
    for pt in POINTS:
        sub = dict(zip(VARS, pt))
        pt.append(sp.nsimplify(-P.subs(sub)))
        pt.append(sp.nsimplify(-Q.subs(sub)))
    IMAGE.extend([Rational(0), Rational(0)])

    FMAP[i] = sp.expand(FMAP[i] - (va + P) * (vb + Q))
    FMAP.append(va + P)
    FMAP.append(vb + Q)
    VARS.extend([va, vb])
    STEPS.append((i, P, Q, n, n + 1))


log("Stage 2: degree reduction to <= 3 ...")
step = 0
while True:
    terms = [(i, c, m) for (i, c, m) in h_terms() if sum(m) >= 4]
    if not terms:
        break
    # deterministic: highest degree first, then component, then monomial
    terms.sort(key=lambda t: (-sum(t[2]), t[0], t[2]))
    i, c, m = terms[0]
    step += 1
    apply_gadget(i, c, m)
    collision_check(f"step {step}")
    keller_spot_check(f"step {step}", npts=1)
    log(f"  gadget {step:2d}: comp {i+1}, killed deg-{sum(m)} term, n = {len(VARS)}")

n1 = len(VARS)
# invariants of the degree <= 3 stage
degs = sorted({sum(m) for (_, _, m) in h_terms()})
assert set(degs) <= {2, 3}, degs
J0 = sp.Matrix(FMAP).jacobian(VARS).subs({v: 0 for v in VARS})
assert J0 == sp.eye(n1)
keller_spot_check("degree-3 stage", npts=3)
collision_check("degree-3 stage")
log(f"Stage 2 OK: {step} gadgets, n1 = {n1}; H-term degrees now {degs}; "
    f"identity linear part; Keller spot-checks pass")

F_deg3 = list(FMAP)          # frozen degree<=3 map on C^{n1}
VARS_deg3 = list(VARS)
PTS_deg3 = [list(p) for p in POINTS]
IMG_deg3 = list(IMAGE)


# ----------------------------------------------------------------------
# Symbolic det proof machinery: replay the gadgets as exact column
# eliminations.  For each gadget (i, P, Q, ia, ib), in the Jacobian of the
# CURRENT map (optionally at scaled argument X -> scale*X):
#     col_j -= dP/dv_j * col_ia + dQ/dv_j * col_ib     (all j < ia)
# turns rows ia, ib into unit rows (asserted), leaving det unchanged;
# Laplace along those unit rows deletes rows/cols ia, ib.  After replaying
# all gadgets backwards we reach J of Ftil (at scaled argument), whose
# 3x3 det sympy computes directly.  Column ops preserve det EXACTLY, so
# this is a complete symbolic proof that det J(F_deg3) == 1.
# ----------------------------------------------------------------------
def symbolic_det_via_replay(scale=None, tag=""):
    subsmap = {v: scale * v for v in VARS_deg3} if scale is not None else {}
    M = sp.Matrix(F_deg3).jacobian(VARS_deg3)
    if subsmap:
        M = M.subs(subsmap, simultaneous=True)
    M = sp.Matrix(M).applyfunc(sp.expand)
    for (i, P, Q, ia, ib) in reversed(STEPS):
        for j in range(ia):
            c1 = sp.diff(P, VARS_deg3[j])
            c2 = sp.diff(Q, VARS_deg3[j])
            if subsmap:
                c1 = c1.subs(subsmap, simultaneous=True)
                c2 = c2.subs(subsmap, simultaneous=True)
            if c1 == 0 and c2 == 0:
                continue
            M[:, j] = (M[:, j] - c1 * M[:, ia] - c2 * M[:, ib]).applyfunc(sp.expand)
        # rows ia, ib must now be unit rows e_ia, e_ib
        for k in range(M.shape[1]):
            assert sp.expand(M[ia, k] - (1 if k == ia else 0)) == 0, (tag, ia, k)
            assert sp.expand(M[ib, k] - (1 if k == ib else 0)) == 0, (tag, ib, k)
        keep = [r for r in range(M.shape[0]) if r not in (ia, ib)]
        M = sp.Matrix(M[keep, keep])
    assert M.shape == (3, 3)
    d = sp.expand(M.det())
    return d


log("Symbolic proof: det J(F_deg3) == 1 via elimination replay ...")
d_plain = symbolic_det_via_replay(scale=None, tag="plain")
assert d_plain == 1, d_plain
log("  det J(F_deg3) == 1  PROVED symbolically (exact, all-argument)")

# ----------------------------------------------------------------------
# Stage 3: BCW doubling  F' = (X + F2(X) + Y, Y - F3(X)) on C^{2 n1}
#   = G(1) o F_deg3[n1] o H(1),  G(1) = (X+Y, Y),  H(1) = (X, Y - F3(X)).
# ----------------------------------------------------------------------
H_of = [sp.expand(f - v) for f, v in zip(F_deg3, VARS_deg3)]
F2part, F3part = [], []
for h in H_of:
    p = sp.Poly(h, *VARS_deg3)
    f2 = sp.Integer(0)
    f3 = sp.Integer(0)
    for c, m in zip(p.coeffs(), p.monoms()):
        if sum(m) == 2:
            f2 += c * sp.prod([v**k for v, k in zip(VARS_deg3, m)])
        elif sum(m) == 3:
            f3 += c * sp.prod([v**k for v, k in zip(VARS_deg3, m)])
        else:
            raise AssertionError("non-quadratic/cubic H term")
    F2part.append(sp.expand(f2))
    F3part.append(sp.expand(f3))
for f, f2, f3, v in zip(F_deg3, F2part, F3part, VARS_deg3):
    assert sp.expand(f - v - f2 - f3) == 0

Wvars = [sp.Symbol(f'w{k+1}') for k in range(n1)]
F_dbl = ([sp.expand(v + f2 + w) for v, f2, w in zip(VARS_deg3, F2part, Wvars)]
         + [sp.expand(w - f3) for w, f3 in zip(Wvars, F3part)])
VARS_dbl = VARS_deg3 + Wvars

# collision transport: input (p, F3(p)), image (Fdeg3(p), 0)
PTS_dbl = []
for pt in PTS_deg3:
    sub = dict(zip(VARS_deg3, pt))
    PTS_dbl.append(list(pt) + [sp.nsimplify(f3.subs(sub)) for f3 in F3part])
IMG_dbl = list(IMG_deg3) + [Rational(0)] * n1

imset = set()
for pt in PTS_dbl:
    sub = dict(zip(VARS_dbl, pt))
    v = tuple(sp.nsimplify(f.subs(sub)) for f in F_dbl)
    imset.add(v)
assert len(imset) == 1 and next(iter(imset)) == tuple(IMG_dbl)
log(f"Stage 3 OK: doubled to C^{2*n1}; collisions transported, common image checked")

# ----------------------------------------------------------------------
# Stage 4: homogenization  L = (X + T^2 Y + T F2(X), Y - F3(X), T)
# ----------------------------------------------------------------------
t = sp.Symbol('t')
L_map = ([sp.expand(v + t**2 * w + t * f2)
          for v, w, f2 in zip(VARS_deg3, Wvars, F2part)]
         + [sp.expand(w - f3) for w, f3 in zip(Wvars, F3part)]
         + [t])
VARS_L = VARS_deg3 + Wvars + [t]
N_FINAL = len(VARS_L)

PTS_L = [list(pt) + [Rational(1)] for pt in PTS_dbl]
IMG_L = list(IMG_dbl) + [Rational(1)]

# (4a) collisions, distinctness
imset = set()
for pt in PTS_L:
    sub = dict(zip(VARS_L, pt))
    v = tuple(sp.nsimplify(f.subs(sub)) for f in L_map)
    imset.add(v)
assert len(imset) == 1 and next(iter(imset)) == tuple(IMG_L)
assert PTS_L[0] != PTS_L[1] and PTS_L[0] != PTS_L[2] and PTS_L[1] != PTS_L[2]

# (4b) identity-plus-cubic-homogeneous form
Hfin = [sp.expand(f - v) for f, v in zip(L_map, VARS_L)]
for h in Hfin:
    if h == 0:
        continue
    p = sp.Poly(h, *VARS_L)
    assert all(sum(m) == 3 for m in p.monoms()), "H not cubic homogeneous"
log(f"Stage 4 OK: L = id + H on C^{N_FINAL}, H homogeneous cubic; "
    f"collisions at t = 1 slice verified")

# ----------------------------------------------------------------------
# Stage 5: exact symbolic Keller proof for L, and nilpotency of JH
# ----------------------------------------------------------------------
# J_L block structure (X-rows | Y-rows | t-row):
#   [ I + t J2(X)   t^2 I    * ]
#   [   -J3(X)        I      0 ]
#   [     0           0      1 ]
# Row ops  X-row_i -= t^2 * Y-row_i  (det-preserving) yield
#   [ I + t J2 + t^2 J3      0      * ]
#   [   -J3                  I      0 ]
#   [    0                   0      1 ]
# Laplace along unit Y-columns and the unit t-row leaves
#   det J_L = det( I + t J2(X) + t^2 J3(X) ) = det J_{F_deg3}(t X),
# which the elimination replay evaluates symbolically to 1.
JL = sp.Matrix(L_map).jacobian(VARS_L)
M = sp.Matrix(JL).applyfunc(sp.expand)
for i in range(n1):
    M[i, :] = (M[i, :] - t**2 * M[n1 + i, :]).applyfunc(sp.expand)
# check Y-columns of X-rows vanished and Y-columns are unit columns
for i in range(n1):
    for j in range(n1, 2 * n1):
        assert M[i, j] == 0
    for r in list(range(n1)) + [2 * n1]:
        assert M[r, n1 + i] == 0 if r != n1 + i else True
# t-row is unit
assert all(M[2 * n1, k] == (1 if k == 2 * n1 else 0) for k in range(N_FINAL))
Mtop = sp.Matrix(M[:n1, :n1])
# entrywise identity: Mtop == J_{F_deg3}(t X)
Jd3 = sp.Matrix(F_deg3).jacobian(VARS_deg3)
Jd3_scaled = Jd3.subs({v: t * v for v in VARS_deg3}, simultaneous=True)
assert (Mtop - Jd3_scaled).applyfunc(sp.expand) == sp.zeros(n1, n1), "Schur mismatch"
log("Stage 5a: det J_L reduced symbolically to det J_{F_deg3}(tX) (exact block elimination)")

d_scaled = symbolic_det_via_replay(scale=t, tag="tX")
assert d_scaled == 1, d_scaled
log("  det J_{F_deg3}(tX) == 1 PROVED symbolically  =>  det J_L == 1 identically. KELLER.")

# Nilpotency of JH: det(I + s*JH) with s fresh reduces the same way to
# det J_{F_deg3}(s t X) == 1; so char(JH) = lambda^N and JH^N = 0.
s = sp.Symbol('s')
JH = sp.Matrix(JL) - sp.eye(N_FINAL)
M = (sp.eye(N_FINAL) + s * JH).applyfunc(sp.expand)
for i in range(n1):
    M[i, :] = (M[i, :] - s * t**2 * M[n1 + i, :]).applyfunc(sp.expand)
for i in range(n1):
    for j in range(n1, 2 * n1):
        assert M[i, j] == 0
Mtop = sp.Matrix(M[:n1, :n1])
Jd3_sscaled = Jd3.subs({v: s * t * v for v in VARS_deg3}, simultaneous=True)
assert (Mtop - Jd3_sscaled).applyfunc(sp.expand) == sp.zeros(n1, n1)
d_ss = symbolic_det_via_replay(scale=s * t, tag="stX")
assert d_ss == 1, d_ss
log("Stage 5b: det(I + s JH) == 1 identically  =>  JH nilpotent (Cayley-Hamilton).")

# belt-and-braces: numeric nilpotency index at a random integer point
sub = {v: Rational(RNG.randint(-5, 5)) for v in VARS_L}
Jn = JH.subs(sub)
Pw = sp.eye(N_FINAL)
idx = None
for k in range(1, N_FINAL + 1):
    Pw = Pw * Jn
    if Pw == sp.zeros(N_FINAL, N_FINAL):
        idx = k
        break
assert idx is not None, "JH not nilpotent at sample point?!"
log(f"  numeric check: JH^{idx} = 0 at a random integer point (index {idx})")

# belt-and-braces: exact det J_L at random integer points
for _ in range(5):
    sub = {v: Rational(RNG.randint(-7, 7)) for v in VARS_L}
    d = JL.subs(sub).det()
    assert d == 1, d
log("  det J_L == 1 at 5 random integer points (exact)")

# ----------------------------------------------------------------------
# Stage 6: outputs
# ----------------------------------------------------------------------
log("")
log("=" * 70)
log(f"RESULT: cubic-homogeneous Keller counterexample G = id + H on C^{N_FINAL}")
log(f"  dimension N = {N_FINAL}  (3 base + 2*{step} gadget + {n1} doubling + 1 homogenizing)")
log(f"  det J_G == 1 identically (proved symbolically);  J_H nilpotent;  H homogeneous cubic")
log(f"  gadget steps: {step}")
log("")

with open(os.path.join(HERE, 'G_map.txt'), 'w') as fh:
    fh.write(f"# Cubic-homogeneous Keller counterexample G = id + H on C^{N_FINAL}\n")
    fh.write(f"# Variables: v1..v{n1} (v1,v2,v3 = x,y,z of Alpoge's F), "
             f"w1..w{n1} (doubling partners), t (homogenizing)\n")
    fh.write(f"# G_i = var_i + H_i, H homogeneous cubic, det JG = 1, JH nilpotent.\n\n")
    for v, f in zip(VARS_L, L_map):
        fh.write(f"G[{v}] = {f}\n")
    fh.write("\n# Colliding points (all three share one image):\n")
    for name, pt in zip("PQR", PTS_L):
        fh.write(f"\n{name} = {[str(c) for c in pt]}\n")
    fh.write(f"\ncommon image = {[str(c) for c in IMG_L]}\n")

with open(os.path.join(HERE, 'collisions.txt'), 'w') as fh:
    for name, pt in zip("PQR", PTS_L):
        fh.write(f"{name}:\n")
        for v, c in zip(VARS_L, pt):
            fh.write(f"  {v} = {c}\n")
    fh.write("common image:\n")
    for v, c in zip(VARS_L, IMG_L):
        fh.write(f"  {v} = {c}\n")

print("\nColliding points (first 12 coords shown; full lists in G_map.txt):")
for name, pt in zip("PQR", PTS_L):
    print(f"  {name} = {[str(c) for c in pt[:12]]} ...")
print(f"  common image = {[str(c) for c in IMG_L[:12]]} ...")
log("STAGE A ALL CHECKS GREEN")


# ======================================================================
# STAGE B: Druzkowski cubic-linear form, via the Gorni-Zampieri pairing
# (Gorni & Zampieri, "On cubic-linear polynomial mappings", Prop 2.1 --
# specialized; our A = D0*B0 needs no matrix inversion at all).
#
# 1. Waring: write the Stage-A cubic homogeneous H as
#        H(x) = -B0 (D0 x)^{*3}
#    using the classical identities
#        a b^2 = ((a+b)^3 + (a-b)^3 - 2 a^3)/6
#        a b c = ((a+b+c)^3 + (a-b-c)^3 - (a+b-c)^3 - (a-b+c)^3)/24 ,
#    D0 (r0 x nA) rows = linear forms (deduped up to sign),
#    B0 (nA x r0) = coefficient matrix.  Verified SYMBOLICALLY.
# 2. GZ padding (B0 has zero rows -- e.g. the t-component and every
#    w-partner of a purely quadratic gadget component -- so B0 alone has
#    no right inverse):
#        B := [B0 | I_nA]   (nA x N),   D := [D0 ; 0]  (N x nA),
#        N := r0 + nA,      C := [0 ; I_nA]  (N x nA),  BC = I.
#    Then A := D*B = [[D0*B0, D0],[0, 0]]  (N x N)  -- no inversion --
#    satisfies A C = D, A M = 0 for M = [I ; -B0] spanning ker B, and
#        B F(C x) = x - B0 (D0 x)^{*3} = G(x)      (the GZ pairing),
#    where  F(X) = X - (A X)^{*3}  on C^N, i.e. writing X = (p, q):
#        F(p, q) = ( p - (D0 (B0 p + q))^{*3},  q ).
#    Each component is  X_k + (l_k(X))^3,  l_k = -(row k of [D0B0 | D0]),
#    with l_k = 0 on the padded block: literal Druzkowski form.
# 3. Keller: J_F = I - 3 diag((AX)^{*2}) D B; Sylvester's identity
#    det(I - UV) = det(I - VU) with U = 3 diag((D B X)^{*2}) D, V = B
#    gives det J_F(X) = det(I - 3 B0 diag((D0 y)^{*2}) D0) at y = B X,
#    i.e. det J_F(X) = det J_G(B X) == 1 by the Stage-A symbolic proof.
#    (rank D0 = nA is asserted => ker A = ker B as in GZ Prop 2.1.)
#    Same chain with a fresh s proves J_{H_F} nilpotent.
# 4. Collisions (GZ Prop 3.1 direction "G not injective => F not
#    injective", made explicit): for G(x1) = G(x2),
#        X1  = C x1 = (0, x1),
#        X2' = C x2 + F(C x1) - F(C x2) = ((D0 x2)^{*3} - (D0 x1)^{*3}, x1),
#    then F(X2') = F(X1) = (-(D0 x1)^{*3}, x1)  [uses B0 ((D0 x2)^{*3} -
#    (D0 x1)^{*3}) = x2 - x1, i.e. G(x1) = G(x2)], and
#    B X1 = x1 != x2 = B X2' certifies X1 != X2'.
# ======================================================================
log("")
log("STAGE B: Druzkowski cubic-linear form via Gorni-Zampieri pairing ...")

nA = N_FINAL
VARS_A = VARS_L

forms = {}          # normalized coeff-tuple -> index
B0 = {}             # (i, j) -> Fraction ;  H_i = - sum_j B0[i,j] * (l_j . x)^3


def form_index(vec):
    """vec: dict var-index -> int coeff.  Normalize sign; return (idx, sign)."""
    items = tuple(sorted((k, c) for k, c in vec.items() if c != 0))
    first = items[0][1]
    sg = 1
    if first < 0:
        sg = -1
        items = tuple((k, -c) for k, c in items)
    j = forms.setdefault(items, len(forms))
    return j, sg


def add_cube(i, coeff, vec):
    """Record H_i += coeff * (l.x)^3 ; since H_i = -sum B0[i,j] l_j^3,
    B0[i,j] -= coeff (sign-normalized: l = sg*l', l^3 = sg*l'^3)."""
    j, sg = form_index(vec)
    c = coeff * sg
    B0[(i, j)] = B0.get((i, j), Fraction(0)) - c


for i, h in enumerate(Hfin):
    if h == 0:
        continue
    p = sp.Poly(h, *VARS_A)
    for co, m in zip(p.coeffs(), p.monoms()):
        c = Fraction(str(co))
        sup = [(k, e) for k, e in enumerate(m) if e]
        if len(sup) == 1:
            k, e = sup[0]
            assert e == 3
            add_cube(i, c, {k: 1})
        elif len(sup) == 2:
            (k1, e1), (k2, e2) = sup
            u, v = (k1, k2) if e1 == 2 else (k2, k1)     # term = c * u^2 v
            add_cube(i, c / 6, {v: 1, u: 1})
            add_cube(i, c / 6, {v: 1, u: -1})
            add_cube(i, -c / 3, {v: 1})
        else:
            (k1, _), (k2, _), (k3, _) = sup              # term = c * u v w
            add_cube(i, c / 24, {k1: 1, k2: 1, k3: 1})
            add_cube(i, c / 24, {k1: 1, k2: -1, k3: -1})
            add_cube(i, -c / 24, {k1: 1, k2: 1, k3: -1})
            add_cube(i, -c / 24, {k1: 1, k2: -1, k3: 1})

# drop exact-zero B0 entries
B0 = {k: v for k, v in B0.items() if v != 0}
r0 = len(forms)
form_rows = [None] * r0
for items, j in forms.items():
    form_rows[j] = dict(items)
log(f"  Waring: r0 = {r0} distinct linear forms; {len(B0)} nonzero B0 entries")

# --- symbolic verification of H(x) = -B0 (D0 x)^{*3} --------------------
form_exprs = []
for j in range(r0):
    e = sp.Integer(0)
    for k, c in form_rows[j].items():
        e += c * VARS_A[k]
    form_exprs.append(e)

Brows = {}
for (i, j), b in B0.items():
    Brows.setdefault(i, []).append((j, b))
for i, h in enumerate(Hfin):
    expr = h
    for j, b in Brows.get(i, []):
        expr = expr + sp.Rational(b.numerator, b.denominator) * form_exprs[j]**3
    assert sp.expand(expr) == 0, f"Waring identity fails in component {i}"
log("  H(x) == -B0 (D0 x)^{*3} verified SYMBOLICALLY, all components")

# calculus closure for the Sylvester chain: the Jacobian of H in Waring form
JH_rhs = sp.zeros(nA, nA)
for (i, j), b in B0.items():
    bb = sp.Rational(b.numerator, b.denominator)
    lj2 = form_exprs[j]**2
    for k, c in form_rows[j].items():
        JH_rhs[i, k] += -3 * bb * c * lj2
assert (JH - JH_rhs).applyfunc(sp.expand) == sp.zeros(nA, nA)
log("  J_H == -3 B0 diag((D0 x)^{*2}) D0 verified SYMBOLICALLY")

# --- rank conditions ----------------------------------------------------
# rank D0 = nA  (=> ker A = ker B0, needed for the pairing / Prop 2.2)
Gram = [[Fraction(0)] * nA for _ in range(nA)]
for j in range(r0):
    row = form_rows[j]
    for k1, c1 in row.items():
        for k2, c2 in row.items():
            Gram[k1][k2] += c1 * c2
GramM = sp.Matrix(nA, nA, lambda a, b: sp.Rational(Gram[a][b].numerator,
                                                   Gram[a][b].denominator))
assert GramM.det() != 0, "rank D0 < nA; append missing coordinate forms"
log(f"  rank D0 = {nA} (Gram determinant nonzero): ker A = ker B0")

# --- the padded Druzkowski map -----------------------------------------
# N = r0 + nA;  X = (p, q);  F(p, q) = (p - (D0 (B0 p + q))^{*3}, q)
NB = r0 + nA
log(f"  GZ padding: B = [B0 | I], D = [D0 ; 0];  N = r0 + {nA} = {NB}")


def matvec_B(X):
    """B . X = B0 p + q  for X = (p, q), length N."""
    out = list(X[r0:])
    for (i, j), b in B0.items():
        if X[j]:
            out[i] += b * X[j]
    return out


def matvec_D0(ynA):
    """D0 . y for y a length-nA list."""
    out = []
    for j in range(r0):
        acc = Fraction(0)
        for k, c in form_rows[j].items():
            acc += c * ynA[k]
        out.append(acc)
    return out


def F_dru(X):
    AX = matvec_D0(matvec_B(X)) + [Fraction(0)] * nA     # A X = D (B X)
    return [xk - axk**3 for xk, axk in zip(X, AX)]


# --- collision transport ------------------------------------------------
PTS_A_frac = [[Fraction(str(c)) for c in pt] for pt in PTS_L]
D0cubes = [[a**3 for a in matvec_D0(pt)] for pt in PTS_A_frac]

X_P = [Fraction(0)] * r0 + PTS_A_frac[0]
X_Q = [d2 - d1 for d2, d1 in zip(D0cubes[1], D0cubes[0])] + PTS_A_frac[0]
X_R = [d3 - d1 for d3, d1 in zip(D0cubes[2], D0cubes[0])] + PTS_A_frac[0]

FP, FQ, FR = F_dru(X_P), F_dru(X_Q), F_dru(X_R)
assert FP == FQ == FR, "Druzkowski collision failed"
assert X_P != X_Q and X_P != X_R and X_Q != X_R
# distinctness certificate: B X_P = x1, B X_Q = x2, B X_R = x3
assert matvec_B(X_P) == PTS_A_frac[0]
assert matvec_B(X_Q) == PTS_A_frac[1]
assert matvec_B(X_R) == PTS_A_frac[2]
# common image sanity: F(X_P) = (-(D0 x1)^{*3}, x1)
assert FP == [-d for d in D0cubes[0]] + PTS_A_frac[0]
log(f"  COLLISION on C^{NB}: three distinct points, one common F-image (exact)")

# --- float spot-check of Keller (belt; the proof is the Sylvester chain)
try:
    import numpy as np
    rng = np.random.default_rng(20260720)
    Bnp = np.zeros((nA, NB))
    for (i, j), b in B0.items():
        Bnp[i, j] = float(b)
    Bnp[:, r0:] = np.eye(nA)
    Dnp = np.zeros((NB, nA))
    for j in range(r0):
        for k, c in form_rows[j].items():
            Dnp[j, k] = float(c)
    Anp = Dnp @ Bnp
    for _ in range(3):
        X = rng.uniform(-0.5, 0.5, NB)
        JF = np.eye(NB) - 3 * np.diag((Anp @ X) ** 2) @ Anp
        dv = np.linalg.det(JF)
        assert abs(dv - 1) < 1e-6, dv
    log("  float spot-check: det J_F = 1.000000 at 3 random points")
except ImportError:
    log("  (numpy unavailable; skipped float spot-check -- Sylvester chain stands)")

# --- outputs ------------------------------------------------------------
with open(os.path.join(HERE, 'F_druzkowski.txt'), 'w') as fh:
    fh.write(f"# Druzkowski cubic-linear counterexample F(X) = X - (A X)^{{*3}} "
             f"on C^{NB}\n")
    fh.write(f"# N = r0 + nA = {r0} + {nA} = {NB};  X = (p_1..p_{r0}, "
             f"q_1..q_{nA})\n")
    fh.write(f"# A = [[D0*B0, D0],[0,0]];  F(p,q) = (p - (D0(B0 p + q))^{{*3}}, q)\n")
    fh.write(f"# i.e. F_k = X_k + (l_k(X))^3, l_k = -(row k of [D0B0 | D0]), "
             f"l_k = 0 for k > r0.\n")
    fh.write(f"# Rows of D0 (the r0 = {r0} linear forms in the {nA} Stage-A "
             f"variables),\n# sparse format 'var_index:coeff':\n\n")
    for j in range(r0):
        row = " ".join(f"{k}:{form_rows[j][k]}" for k in sorted(form_rows[j]))
        fh.write(f"D0[{j}] = {row}\n")
    fh.write(f"\n# Nonzero entries of B0 ({nA} x {r0}), format (i, j): value\n\n")
    for (i, j), b in sorted(B0.items()):
        fh.write(f"B0[{i},{j}] = {b}\n")
    fh.write(f"\n# Colliding points on C^{NB} (X_P, X_Q, X_R; common image FP):\n\n")
    for name, vec in (("X_P", X_P), ("X_Q", X_Q), ("X_R", X_R), ("F_image", FP)):
        fh.write(f"{name} = [{', '.join(str(c) for c in vec)}]\n\n")

log("")
log("=" * 70)
log(f"STAGE B RESULT: Druzkowski map F(X) = X - (A X)^{{*3}} on C^{NB}")
log(f"  A = [[D0*B0, D0],[0,0]] explicit, no inversion;")
log(f"  Keller via Sylvester + Stage-A symbolic proof;")
log(f"  J_H nilpotent (same chain with fresh s);")
log(f"  three explicit distinct colliding points verified exactly.")
log("ALL STAGE A + STAGE B CHECKS GREEN")
