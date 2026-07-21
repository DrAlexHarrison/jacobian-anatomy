#!/usr/bin/env python3
"""
P3-Mathieu, Phase 2 — self-contained runnable certificate that the Mathieu
conjecture is FALSE on SU(79).

Logical chain (each link machine-checked below):

  (A) L = id + H on C^79 (druzkowski/G_map.txt) is cubic-homogeneous:
      every H_i is homogeneous of degree exactly 3, or zero.            [CHECK 1]
  (B) L is Keller: det J L = 1.  (Full symbolic proof lives in
      druzkowski/construct_verify.py; here: exact rational spot-checks
      at random points + at the collision points.)                     [CHECK 2]
  (C) L is NOT injective: three DISTINCT points P, Q, R in C^79 (exact
      rationals, druzkowski/G_map.txt) satisfy L(P) = L(Q) = L(R).      [CHECK 3]
      => the formal inverse F of L is NOT a polynomial: if it were, both
      formal identities L∘F = id and F∘L = id would be polynomial
      identities, and F∘L = id forces injectivity of L, contradicting
      CHECK 3. A power series with finitely many nonzero homogeneous
      parts is a polynomial, so F^{(2k+1)} != 0 for infinitely many k.  [R4]
  (D) Mathieu's reduction (PRIMARY: Mathieu 1997 — Prop 2.2(ii), Cor 1.3/1.7,
      Props 3.3/3.4, Formula 4.4. Zwart arXiv:2511.16561 is expository-
      secondary: its Thm 4.23 has a proof gap at the fixed-ξ step; the sound
      route is Mathieu's Baire-generic ξ. See audit/referee2-su79.md.)  [R1]:
      with Q_elt = sum_i H_i ⊗ ∂_i  in S^3 C^79 ⊗ (C^79)*,
        (Q_elt^k)_{2k·ω1}            = 0  for all k   <==  det J L = 1  (CHECK 2;
        only this direction is proven, and only it is needed)  [D3],
        (Q_elt^k)_{(2k+1)ω1 + ω_{78}} = k!·Σ_i F_i^{(2k+1)} ⊗ ∂_i =: Ψ(Q_elt^k).
      Non-termination of F (from C) => Ψ(Q_elt^k) != 0 for infinitely many k.
      This refutes Mathieu Prop 2.2(ii)'s conclusion for the graded algebra
      A' ∗ A(2ω₁), so by Cor 1.3/1.7 (per-group at proof level; the Peter-Weyl
      dictionary is Mathieu Lemma 1.2 = Zwart Lemmas 4.20-4.22 in v2 numbering)
      MC(SU(79)) is FALSE: a finite-type pair (f, g) on SU(79) with
      ∫ f^n = 0 ∀n but ∫ f^n g != 0 for infinitely many n EXISTS.
      (Existence, not an explicit pair: an explicit (f,g) additionally
      requires Mathieu's Cartan twist by a Baire-generic ξ — future work.) [R2]

  ILLUSTRATION: we compute the low-degree homogeneous parts F^{(3)}, F^{(5)},
  ... of the formal inverse in the SAME sparse engine and exhibit explicit
  nonzero components (concrete Ψ(Q_elt^k) != 0), plus a NEGATIVE CONTROL — a
  cubic-homogeneous *automorphism* whose inverse provably terminates — so the
  machinery demonstrably can report both outcomes.

Sign convention: druzkowski writes L = X + H; Mathieu/Zwart ingest X - H.
Only H' = -H changes, same graded class; nothing downstream depends on it.

Sparse, checkpointed, exact (fractions.Fraction). No CAS needed for the inverse;
sympy used only to parse G_map.txt robustly.
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
log(f"CHECK 1: parsed L = id + H on C^{N}; {nz}/{N} components have nonzero H; "
    f"every H-term verified homogeneous of degree 3.  PASS")
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
log(f"CHECK 3: L(P)=L(Q)=L(R): {same}; P,Q,R distinct: {distinct}")
log(f"         common image L(P) = (v1,v2,v3,t) = "
    f"({img[0]},{img[1]},{img[2]},...,{img[VIDX['t']]}); "
    f"rest zero: {all(x==0 for i,x in enumerate(img) if i not in (2, VIDX['t']))}")
assert same and distinct, "COLLISION CHECK FAILED"
log("CHECK 3: non-injectivity CERTIFIED (exact rational arithmetic).  PASS")

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
log(f"CHECK 2: det J L == 1 at 4 random rational points: {ok}; "
    f"at collision P: {d_at_P == 1}.  PASS (full symbolic proof: druzkowski/)")
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
            log(f"  [{label}] F^({d}) NONZERO: component {found[1]} has term "
                f"({found[3]})·{found[2]} ...  (total stored terms: {terms_total})")
        else:
            log(f"  [{label}] F^({d}) == 0  (total stored terms: {terms_total})")
        if ckpt:
            json.dump({"label": label, "N": N, "next_deg": d+2,
                       "nz": {str(k): v for k, v in nonzero_degrees.items()},
                       "F": [{",".join(map(str, m)): str(c) for m, c in comp.items()}
                             for comp in F]}, open(ckpt, "w"))
    return F, nonzero_degrees

# ----------------------------------------------------------------------
# ILLUSTRATION: low-degree inverse of L (our counterexample)
# ----------------------------------------------------------------------
D = int(sys.argv[1]) if len(sys.argv) > 1 else 11  # 11 = the refereed run; NOTE §7 cites degrees 3-11
log(f"Computing formal inverse of L to total degree {D} (sparse, checkpointed)...")
_, nzL = formal_inverse(H, N, D, "L", ckpt=CKPT)
log(f"ILLUSTRATION: nonzero inverse components of L at degrees "
    f"{sorted(nzL)} (all odd; each an explicit Ψ(Q_elt^k)!=0 witness). "
    f"Non-termination past any finite D is guaranteed by CHECK 3, not by this "
    f"finite scan.")

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

log("NEGATIVE CONTROL: cubic-homogeneous triangular automorphism on C^3 ...")
nzA = neg_control()
top = max(nzA) if nzA else 0
log(f"NEGATIVE CONTROL: inverse nonzero only at degrees {sorted(nzA)}; "
    f"F^(d)==0 for all odd d>{top} up to 15 — inverse TERMINATES => automorphism. "
    f"Machinery distinguishes automorphism (terminates) from L (does not).")

# ----------------------------------------------------------------------
print()
print("="*72)
print("CERTIFICATE: the Mathieu conjecture is FALSE on SU(79).")
print("  CHECK 1  L = id + H on C^79 is cubic-homogeneous.               PASS")
print("  CHECK 2  det J L = 1 (Keller) — exact spot-checks + druzkowski. PASS")
print("  CHECK 3  L not injective: 3 distinct exact preimages collide.   PASS")
print("  => formal inverse of L is non-terminating (no polynomial inverse)")
print("  => Ψ(Q_elt^k) != 0 for infinitely many k, while (Q_elt^k)_{2kω1}=0 ∀k")
print("  => MC(SU(79)) FALSE: a finite-type pair (f,g) with ∫f^n=0 ∀n, ∫f^n g≠0")
print("     ∞-often EXISTS (explicit pair = future work via the Cartan twist).")
print(f"  Illustration: inverse nonzero at degrees {sorted(nzL)} (Ψ witnesses).")
print(f"  Negative control terminates at degree {top}: machinery is falsifiable.")
print("="*72)
