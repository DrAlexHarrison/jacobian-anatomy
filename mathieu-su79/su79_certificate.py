#!/usr/bin/env python3
"""Checks for the SU(79) Mathieu result (see README.md for the argument).

Verifies: H is cubic-homogeneous; det JL = 1 at sample rational points;
the three-point collision, in exact arithmetic; nonzero homogeneous
components of the formal inverse at odd degrees up to the bound (default
11, argv[1] to change); and a terminating control map. Needs sympy (for
parsing G_map.txt). Writes a resumable checkpoint json.
"""
import sys, os, re, json, time
from fractions import Fraction
from itertools import combinations_with_replacement

HERE = os.path.dirname(os.path.abspath(__file__))
# G_map.txt lives at <repo>/druzkowski/; this script ships both as
# mathieu-su79/su79_certificate.py (public tree) and
# programs/p3-mathieu/su79_certificate.py (working tree), so search upward.
_GMAP_CANDIDATES = [os.path.join(HERE, "..", "druzkowski", "G_map.txt"),
                    os.path.join(HERE, "..", "..", "druzkowski", "G_map.txt")]
GMAP = next((p for p in _GMAP_CANDIDATES if os.path.exists(p)),
            _GMAP_CANDIDATES[0])
CKPT = os.path.join(HERE, "su79_inverse_ckpt.json")
t0 = time.time()
def log(m): print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

# ----------------------------------------------------------------------
# Parse G_map.txt into: variable order, and H_i as list of (coeff, (idx,idx,idx))
# where the triple lists the three variable indices of the cubic monomial
# (with multiplicity). Zero H_i -> empty list.
# ----------------------------------------------------------------------
import sympy as sp

def parse_gmap(path):
    txt = open(path).read()
    # collect "G[var] = expr" up to the blank-line/comment block
    rows = {}
    for m in re.finditer(r'^G\[([a-zA-Z]\w*)\]\s*=\s*(.+)$', txt, re.M):
        rows[m.group(1)] = m.group(2).strip()
    # variable order: v1..v39, w1..w39, t  (as they appear)
    order = [v for v in rows.keys()]
    idx = {v: i for i, v in enumerate(order)}
    n = len(order)
    syms = {v: sp.Symbol(v) for v in order}
    H = []  # H[i] = list of (Fraction coeff, (i,j,k) sorted triple)
    for v in order:
        expr = sp.expand(sp.sympify(rows[v], locals=syms) - syms[v])  # strip +var
        terms = []
        if expr != 0:
            poly = sp.Poly(expr, *[syms[u] for u in order])
            for monom, coeff in poly.terms():
                deg = sum(monom)
                assert deg == 3, f"H[{v}] has a degree-{deg} term — not cubic-homog!"
                triple = []
                for j, e in enumerate(monom):
                    triple += [j] * e
                terms.append((Fraction(int(coeff.p), int(coeff.q))
                              if coeff.is_Rational else Fraction(str(coeff)),
                              tuple(triple)))
        H.append(terms)
    return order, idx, n, H

order, VIDX, N, H = parse_gmap(GMAP)
nz = sum(1 for h in H if h)
log(f"check 1: L = id + H on C^{N}; {nz}/{N} components with nonzero H; all H-terms degree 3")
assert N == 79, N

# ----------------------------------------------------------------------
# Exact evaluation L(X) = X + H(X) at a rational point (vector of Fraction)
# ----------------------------------------------------------------------
def evalL(X):
    Y = list(X)
    for i, terms in enumerate(H):
        s = Fraction(0)
        for c, (a, b, d) in terms:
            s += c * X[a] * X[b] * X[d]
        Y[i] += s
    return Y

# ----------------------------------------------------------------------
# CHECK 3: the three collision points (exact rationals) from G_map.txt header
# ----------------------------------------------------------------------
def parse_points(path):
    txt = open(path).read()
    pts = {}
    for name in ("P", "Q", "R"):
        m = re.search(rf'^{name}\s*=\s*\[(.*?)\]', txt, re.S | re.M)
        if not m:
            continue
        vals = [Fraction(s.strip().strip("'\"")) for s in m.group(1).split(",")]
        pts[name] = vals
    return pts

pts = parse_points(GMAP)
assert set(pts) >= {"P", "Q", "R"}, f"missing collision pts: {set(pts)}"
imgs = {k: evalL(v) for k, v in pts.items()}
same = imgs["P"] == imgs["Q"] == imgs["R"]
distinct = len({tuple(pts["P"]), tuple(pts["Q"]), tuple(pts["R"])}) == 3
img = imgs["P"]
# image should be the common target (0,0,-1/4,0,...,0, t=1)
log(f"check 3: L(P)=L(Q)=L(R): {same}; P,Q,R distinct: {distinct}")
log(f"         common image L(P) = (v1,v2,v3,t) = "
    f"({img[0]},{img[1]},{img[2]},...,{img[VIDX['t']]}); "
    f"rest zero: {all(x==0 for i,x in enumerate(img) if i not in (2, VIDX['t']))}")
assert same and distinct, "COLLISION CHECK FAILED"
log("check 3: collision verified in exact arithmetic")

# ----------------------------------------------------------------------
# CHECK 2: det J L = 1 — exact rational spot-checks (full symbolic proof is in
# druzkowski/construct_verify.py). J L = I + J H;  det at several points.
# ----------------------------------------------------------------------
import random
def jacobian_det_at(X):
    # J H[i][j] = d H_i / d x_j ; build dense 79x79 rational matrix, det via
    # fraction-free Bareiss.
    M = [[Fraction(0)]*N for _ in range(N)]
    for i in range(N):
        M[i][i] = Fraction(1)
    for i, terms in enumerate(H):
        for c, tr in terms:
            for pos in range(3):
                j = tr[pos]
                others = [tr[k] for k in range(3) if k != pos]
                M[i][j] += c * X[others[0]] * X[others[1]]
    # Bareiss determinant (exact)
    A = [row[:] for row in M]; sign = 1; prev = Fraction(1)
    for k in range(N):
        if A[k][k] == 0:
            sw = next((r for r in range(k+1, N) if A[r][k] != 0), None)
            if sw is None:
                return Fraction(0)
            A[k], A[sw] = A[sw], A[k]; sign = -sign
        for i in range(k+1, N):
            for j in range(k+1, N):
                A[i][j] = (A[i][j]*A[k][k] - A[i][k]*A[k][j]) / prev
        prev = A[k][k]
    return sign * A[N-1][N-1]

rng = random.Random(20260720)
ok = True
for _ in range(4):
    X = [Fraction(rng.randint(-4, 4), rng.randint(1, 4)) for _ in range(N)]
    d = jacobian_det_at(X)
    ok &= (d == 1)
d_at_P = jacobian_det_at(pts["P"])
log(f"check 2: det JL = 1 at 4 random points: {ok}; at P: {d_at_P == 1} (symbolic proof: druzkowski/)")
assert ok and d_at_P == 1

# ----------------------------------------------------------------------
# Sparse polynomial engine for the formal inverse.
# Poly = dict{ exponent-tuple(len N) -> Fraction }.  Truncate at total degree D.
# ----------------------------------------------------------------------
def padd(p, q):
    r = dict(p)
    for m, c in q.items():
        r[m] = r.get(m, Fraction(0)) + c
        if r[m] == 0: del r[m]
    return r

def pmul(p, q, D):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(a+b for a, b in zip(m1, m2))
            if sum(m) > D:
                continue
            v = r.get(m, Fraction(0)) + c1*c2
            if v == 0: r.pop(m, None)
            else: r[m] = v
    return r

def var_poly(i):
    e = [0]*N; e[i] = 1; return {tuple(e): Fraction(1)}

def deg_part(p, d):
    return {m: c for m, c in p.items() if sum(m) == d}

def compose_H_i(terms, F, D):
    """H_i(F) truncated to degree D, given F = list of polys."""
    acc = {}
    for c, (a, b, d) in terms:
        prod = pmul(pmul(F[a], F[b], D), F[d], D)
        for m, cc in prod.items():
            v = acc.get(m, Fraction(0)) + c*cc
            if v == 0: acc.pop(m, None)
            else: acc[m] = v
    return acc

def formal_inverse(H, N, D, label, ckpt=None):
    """Solve F = Y - H(F) degree by degree (odd degrees only survive)."""
    F = [var_poly(i) for i in range(N)]     # degree-1 part
    nonzero_degrees = {}
    start_deg = 3
    if ckpt and os.path.exists(ckpt):
        try:
            data = json.load(open(ckpt))
            if data.get("label") == label and data.get("N") == N:
                F = [{tuple(map(int, k.split(","))): Fraction(v)
                      for k, v in comp.items()} for comp in data["F"]]
                start_deg = data["next_deg"]
                nonzero_degrees = {int(k): v for k, v in data["nz"].items()}
                log(f"  [{label}] resumed from checkpoint at degree {start_deg}")
        except Exception as e:
            log(f"  [{label}] checkpoint ignored ({e})")
    for d in range(start_deg, D+1, 2):
        newF = [dict(F[i]) for i in range(N)]
        found = None
        for i in range(N):
            if not H[i]:
                continue
            hi = compose_H_i(H[i], F, d)         # H_i(F) up to degree d
            part = deg_part(hi, d)               # degree-d slice
            if part:
                for m, c in part.items():
                    newF[i][m] = newF[i].get(m, Fraction(0)) - c
                    if newF[i][m] == 0: del newF[i][m]
                if found is None:
                    # record an explicit witness monomial for this degree
                    m, c = next(iter(part.items()))
                    wit = "*".join(f"{order[j]}^{e}" for j, e in enumerate(m) if e)
                    found = (i, order[i], wit, str(-c))
        F = newF
        terms_total = sum(len(c) for c in F)
        if found:
            nonzero_degrees[d] = found[1:]
            log(f"  [{label}] F^({d}) nonzero: component {found[1]}, term ({found[3]})*{found[2]} ({terms_total} stored terms)")
        else:
            log(f"  [{label}] F^({d}) = 0 ({terms_total} stored terms)")
        if ckpt:
            json.dump({"label": label, "N": N, "next_deg": d+2,
                       "nz": {str(k): v for k, v in nonzero_degrees.items()},
                       "F": [{",".join(map(str, m)): str(c) for m, c in comp.items()}
                             for comp in F]}, open(ckpt, "w"))
    return F, nonzero_degrees

# ----------------------------------------------------------------------
# ILLUSTRATION: low-degree inverse of L (our counterexample)
# ----------------------------------------------------------------------
D = int(sys.argv[1]) if len(sys.argv) > 1 else 11  # NOTE section 7 cites degrees 3-11
log(f"Computing formal inverse of L to total degree {D} (sparse, checkpointed)...")
_, nzL = formal_inverse(H, N, D, "L", ckpt=CKPT)
log(f"inverse of L nonzero at degrees {sorted(nzL)} (non-termination follows from check 3, not this finite scan)")

# ----------------------------------------------------------------------
# NEGATIVE CONTROL: a cubic-homogeneous *triangular automorphism* on C^3,  [R6]
# H_i depends only on strictly-later variables => nilpotent => inverse
# TERMINATES.  Shows the machinery reports termination correctly.
# ----------------------------------------------------------------------
def neg_control():
    # Strictly-triangular cubic automorphism on C^3 whose inverse TERMINATES at
    # degree 9 (well below the degree-15 scan), so termination is visible:
    #   y0 = x0 + x1^3,  y1 = x1 + x2^3,  y2 = x2.
    # Inverse: x2=y2, x1=y1-y2^3, x0 = y0-(y1-y2^3)^3  -> highest term y2^9.
    n = 3
    Hn = [[] for _ in range(n)]
    Hn[0] = [(Fraction(1), (1, 1, 1))]      # H_0 = x1^3
    Hn[1] = [(Fraction(1), (2, 2, 2))]      # H_1 = x2^3
    # H_2 = 0  (last variable fixed)  -> whole map is a triangular automorphism
    global N, order
    savedN, savedorder = N, order
    N, order = n, [f"x{i}" for i in range(n)]
    _, nz = formal_inverse(Hn, n, 15, "AUTO", ckpt=None)
    N, order = savedN, savedorder
    return nz

log("control: triangular cubic automorphism on C^3")
nzA = neg_control()
top = max(nzA) if nzA else 0
log(f"control: inverse nonzero only at degrees {sorted(nzA)}; zero for odd d > {top} up to 15 (terminates)")

# ----------------------------------------------------------------------
print()
print(f"checks complete: H cubic-homogeneous; det JL = 1 at sample points;")
print(f"three distinct points collide; inverse nonzero at degrees {sorted(nzL)};")
print(f"control inverse terminates at degree {top}.")
print("See README.md for the statement this supports and audit/su79-review.md")
print("for the citation chain.")
