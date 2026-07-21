#!/usr/bin/env python3
"""
REFEREE-2 independent verification of the SU(79) Mathieu-certificate compute
layer (public copy: mathieu-su79/su79_certificate.py). Written from scratch
per the referee independence mandate: own script, exact arithmetic, no code
shared with the certificate. Deliberate implementation differences:

  - parser: explicit expected variable order v1..v39,w1..w39,t, asserted
    against file order (the certificate trusts file order silently);
  - monomials: sorted index-tuples (certificate: 79-length exponent vectors);
  - determinant: plain Gaussian elimination over Fraction (certificate:
    fraction-free Bareiss);
  - inverse engine: degree-sliced recursion F^(d) = -[H(F)]^(d) (certificate:
    full-dict recompute per degree);
  - NEW CHECK A: full term-by-term comparison of my inverse against their
    checkpoint su79_inverse_ckpt.json (two independent engines, exact match);
  - NEW CHECK B: non-tautological two-sided composition test
    F(L(t*u)) == t*u mod t^13 along random rational directions u;
  - NEW CHECK C: exact-coefficient extraction of the five witness monomials
    from su79_run.log (not just "some nonzero term");
  - NEW CONTROL D: 79-variable cubic-homogeneous automorphism (v,w,t) ->
    (v + t^2 w, w, t): engine must terminate at degree 3 in the SAME ambient
    dimension as L (the certificate's control lives in dim 3);
  - NEW CONTROL E: closed-form check of the C^3 triangular control inverse;
  - NEW CONTROL F: broken-coefficient L~ must FAIL the det==1 spot-check
    (proves the Keller spot-check has teeth).
"""
import json, os, re, sys, time
from fractions import Fraction
from random import Random

HERE = os.path.dirname(os.path.abspath(__file__))
JC = os.path.dirname(HERE)  # repo root; this file lives in audit/
GMAP = f"{JC}/druzkowski/G_map.txt"
# The checkpoint is produced by running su79_certificate.py; in the public
# tree that is mathieu-su79/ (run the certificate first), in the working
# tree programs/p3-mathieu/.
_CKPTS = [f"{JC}/mathieu-su79/su79_inverse_ckpt.json",
          f"{JC}/programs/p3-mathieu/su79_inverse_ckpt.json"]
CKPT = next((p for p in _CKPTS if os.path.exists(p)), _CKPTS[0])
D = 11
t0 = time.time()
def log(m): print(f"[{time.time()-t0:8.1f}s] {m}", flush=True)
FAILURES = []
def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    log(f"{tag}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok: FAILURES.append(name)

# ---------------------------------------------------------------- parse
EXPECTED = [f"v{i}" for i in range(1, 40)] + [f"w{i}" for i in range(1, 40)] + ["t"]
VI = {v: i for i, v in enumerate(EXPECTED)}
N = 79

def parse_gmap():
    txt = open(GMAP).read()
    keys, H = [], [[] for _ in range(N)]
    for m in re.finditer(r'^G\[(\w+)\]\s*=\s*(.+)$', txt, re.M):
        var, expr = m.group(1), m.group(2)
        keys.append(var)
        i = VI[var]
        # tokenize into +/- separated terms; my own parser, no sympy
        expr = expr.replace("- ", "+ -").replace("+ ", "|")
        for term in [s.strip() for s in expr.split("|") if s.strip()]:
            neg = term.startswith("-")
            term = term.lstrip("-").strip().replace("**", "^")
            coeff = Fraction(1)
            idxs = []
            # grammar of this file: factor * factor * ... , factor is
            # int | name | name^int, any factor may carry a trailing /int
            for factor in term.split("*"):
                factor = factor.strip()
                if "/" in factor:
                    factor, dv = factor.split("/")
                    coeff /= int(dv)
                if re.fullmatch(r'\d+', factor):
                    coeff *= int(factor)
                elif "^" in factor:
                    name, e = factor.split("^")
                    idxs += [VI[name]] * int(e)
                else:
                    idxs.append(VI[factor])
            if neg: coeff = -coeff
            if idxs == [i] and coeff == 1:
                continue                      # the identity part var_i
            H[i].append((coeff, tuple(sorted(idxs))))
    return keys, H

keys, H = parse_gmap()
check("parse: variable order is v1..v39,w1..w39,t", keys == EXPECTED)
allmono = [(i, c, m) for i in range(N) for c, m in H[i]]
check("CHECK 1 (indep): every H-term is degree exactly 3",
      all(len(m) == 3 for _, _, m in allmono), f"{len(allmono)} cubic terms")
nz = sum(1 for h in H if h)
check("CHECK 1 (indep): 51/79 components nonzero", nz == 51, f"got {nz}")

# cross-validate my parser against sympy's on random points (catches
# parser bugs in EITHER direction; sympy is the certificate's parser)
import sympy as sp
syms = {v: sp.Symbol(v) for v in EXPECTED}
rows = {m.group(1): m.group(2) for m in
        re.finditer(r'^G\[(\w+)\]\s*=\s*(.+)$', open(GMAP).read(), re.M)}
rng = Random(987654321)
ok = True
for trial in range(3):
    X = [Fraction(rng.randint(-9, 9), rng.randint(1, 7)) for _ in range(N)]
    sub = {syms[v]: sp.Rational(X[VI[v]].numerator, X[VI[v]].denominator)
           for v in EXPECTED}
    for v in EXPECTED:
        mine = X[VI[v]] + sum(c * X[m[0]] * X[m[1]] * X[m[2]] for c, m in H[VI[v]])
        ref = sp.Rational(sp.sympify(rows[v], locals=syms).subs(sub))
        ok &= (Fraction(ref.p, ref.q) == mine)
check("parser cross-validation vs sympy at 3 random points", ok)

# ---------------------------------------------------------------- collisions
def parse_pts():
    txt = open(GMAP).read()
    out = {}
    for name in ("P", "Q", "R"):
        m = re.search(rf"^{name} = \[(.*?)\]", txt, re.M | re.S)
        out[name] = [Fraction(s.strip().strip("'")) for s in m.group(1).split(",")]
    return out

def evL(X):
    return [X[i] + sum(c * X[a] * X[b] * X[d] for c, (a, b, d) in H[i])
            for i in range(N)]

pts = parse_pts()
P, Q, R = pts["P"], pts["Q"], pts["R"]
check("CHECK 3 (indep): P,Q,R have 79 coords each",
      all(len(p) == N for p in (P, Q, R)))
check("CHECK 3 (indep): P,Q,R pairwise distinct",
      P != Q and P != R and Q != R)
iP, iQ, iR = evL(P), evL(Q), evL(R)
check("CHECK 3 (indep): L(P)=L(Q)=L(R) exactly", iP == iQ == iR)
check("CHECK 3 (indep): common image = (0,0,-1/4,0,...,0,t=1)",
      iP[2] == Fraction(-1, 4) and iP[N-1] == 1 and
      all(iP[j] == 0 for j in range(N) if j not in (2, N-1)))

# ---------------------------------------------------------------- det J L
def detJL_at(X, Hloc=None):
    Hloc = H if Hloc is None else Hloc
    M = [[Fraction(0)] * N for _ in range(N)]
    for i in range(N):
        M[i][i] = Fraction(1)
    for i in range(N):
        for c, tr in Hloc[i]:
            for pos in range(3):
                j = tr[pos]
                o = [tr[k] for k in range(3) if k != pos]
                M[i][j] += c * X[o[0]] * X[o[1]]
    # plain Gaussian elimination with pivoting, exact Fractions
    det = Fraction(1)
    for k in range(N):
        piv = next((r for r in range(k, N) if M[r][k] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != k:
            M[k], M[piv] = M[piv], M[k]
            det = -det
        det *= M[k][k]
        inv = 1 / M[k][k]
        for r in range(k + 1, N):
            if M[r][k] != 0:
                f = M[r][k] * inv
                M[r] = [M[r][j] - f * M[k][j] for j in range(N)]
    return det

rng = Random(31415926)
ok = True
for trial in range(10):
    X = [Fraction(rng.randint(-12, 12), rng.randint(1, 9)) for _ in range(N)]
    ok &= (detJL_at(X) == 1)
check("CHECK 2 (indep): det J L = 1 at 10 fresh random rational points", ok)
check("CHECK 2 (indep): det J L = 1 at P, Q, R",
      detJL_at(P) == 1 and detJL_at(Q) == 1 and detJL_at(R) == 1)
log("        (identical det==1 on all of C^79: druzkowski symbolic proof, "
    "VETTER post-hoc PASS; spot-checks here are the consumption audit)")

# CONTROL F: break one coefficient; the spot-check must catch it
Hbad = [list(h) for h in H]
c0, m0 = Hbad[0][0]
Hbad[0] = [(c0 * 2, m0)] + Hbad[0][1:]
X = [Fraction(rng.randint(-5, 5), rng.randint(1, 4)) for _ in range(N)]
check("CONTROL F: broken coefficient detected (det != 1)",
      detJL_at(X, Hbad) != 1)

# ---------------------------------------------------------------- inverse
# degree-sliced engine: F[i] = {deg: {sorted-index-tuple: Fraction}}
def mono_mul(m1, m2):
    return tuple(sorted(m1 + m2))

def slice_mul(s1, s2):
    out = {}
    for m1, c1 in s1.items():
        for m2, c2 in s2.items():
            m = mono_mul(m1, m2)
            v = out.get(m, Fraction(0)) + c1 * c2
            if v == 0: out.pop(m, None)
            else: out[m] = v
    return out

def odd_triples(d):
    for a in range(1, d + 1, 2):
        for b in range(1, d - a + 1, 2):
            c = d - a - b
            if c >= 1 and c % 2 == 1:
                yield a, b, c

def invert(Hloc, n, Dmax):
    F = [{1: {(i,): Fraction(1)}} for i in range(n)]
    for d in range(3, Dmax + 1, 2):
        for i in range(n):
            acc = {}
            for c, (a, b, e) in Hloc[i]:
                for da, db, dc in odd_triples(d):
                    sa, sbv, sc = (F[a].get(da), F[b].get(db), F[e].get(dc))
                    if not sa or not sbv or not sc:
                        continue
                    prod = slice_mul(slice_mul(sa, sbv), sc)
                    for m, cc in prod.items():
                        v = acc.get(m, Fraction(0)) - c * cc
                        if v == 0: acc.pop(m, None)
                        else: acc[m] = v
            if acc:
                F[i][d] = acc
        tot = sum(len(s) for f in F for s in f.values())
        nzd = sorted({dd for f in F for dd in f if dd > 1 and f[dd]})
        log(f"  inverse degree {d} done; nonzero degrees so far {nzd}; "
            f"total terms {tot}")
    return F

log("computing formal inverse of L to degree 11, independent engine ...")
F = invert(H, N, D)
nzdegs = sorted({d for f in F for d in f if d > 1 and f[d]})
check("inverse nonzero at every odd degree 3..11", nzdegs == [3, 5, 7, 9, 11])

# CHECK C: the five logged witness coefficients, exact
V = lambda name: VI[name]
wits = [
    (3,  (V("v16"), V("v17"), V("t")),                        Fraction(1)),
    (5,  (V("v16"), V("w17")) + (V("t"),) * 3,                Fraction(-1)),
    (7,  (V("v1"),) * 3 + (V("v3"),) + (V("t"),) * 3,         Fraction(-4)),
    (9,  (V("v1"),) * 3 + (V("v2"),) * 2 + (V("t"),) * 4,     Fraction(5, 2)),
    (11, (V("v1"),) * 4 + (V("v2"), V("v3")) + (V("t"),) * 5, Fraction(-15)),
]
ok = True
for d, mono, cexp in wits:
    got = F[0].get(d, {}).get(tuple(sorted(mono)), Fraction(0))
    if got != cexp:
        ok = False
        log(f"  MISMATCH deg {d}: expected {cexp}, got {got}")
check("CHECK C: all 5 logged witness coefficients match exactly (component v1)", ok)

# CHECK A: full comparison against their checkpoint
if os.path.exists(CKPT):
    data = json.load(open(CKPT))
    ok = data.get("label") == "L" and data.get("N") == N
    mismatches = 0
    if ok:
        theirs = data["F"]
        for i in range(N):
            mine_all = {}
            for d, s in F[i].items():
                for m, c in s.items():
                    ex = [0] * N
                    for j in m: ex[j] += 1
                    mine_all[",".join(map(str, ex))] = c
            their_i = {k: Fraction(v) for k, v in theirs[i].items()
                       if sum(map(int, k.split(","))) <= D}
            if mine_all != their_i:
                mismatches += 1
    check("CHECK A: term-by-term match with su79_inverse_ckpt.json (deg<=11)",
          ok and mismatches == 0, f"ckpt next_deg={data.get('next_deg')}")
else:
    check("CHECK A: checkpoint file present", False, "missing")

# CHECK B: two-sided composition F(L(t*u)) == t*u mod t^13, 3 random directions
T = 12
ok = True
for trial in range(3):
    u = [Fraction(rng.randint(-6, 6), rng.randint(1, 5)) for _ in range(N)]
    Lu = evL(u)                          # L(t*u)_i = t*u_i + t^3*H_i(u)
    coord = []                           # sparse univariate polys in t
    for i in range(N):
        pc = {1: u[i]}
        h = Lu[i] - u[i]
        if h != 0: pc[3] = h
        coord.append(pc)
    for i in range(N):
        acc = {}
        for d, s in F[i].items():
            for m, c in s.items():
                prod = {0: c}
                for j in m:
                    nxt = {}
                    for a, ca in prod.items():
                        for b, cb in coord[j].items():
                            if a + b > T: continue
                            v = nxt.get(a + b, Fraction(0)) + ca * cb
                            if v == 0: nxt.pop(a + b, None)
                            else: nxt[a + b] = v
                    prod = nxt
                    if not prod: break
                for a, ca in prod.items():
                    v = acc.get(a, Fraction(0)) + ca
                    if v == 0: acc.pop(a, None)
                    else: acc[a] = v
        good = (acc.get(1, Fraction(0)) == u[i]) and \
               all(a == 1 for a in acc)
        ok &= good
check("CHECK B: F(L(t*u)) == t*u mod t^13 along 3 random directions", ok)

# ---------------------------------------------------------------- controls
# CONTROL D: 79-var cubic-homogeneous automorphism (v,w,t)->(v+t^2 w,w,t)
Haut = [[] for _ in range(N)]
for i in range(39):
    Haut[i] = [(Fraction(1), tuple(sorted((39 + i, N - 1, N - 1))))]
Faut = invert(Haut, N, 9)
nzA = sorted({d for f in Faut for d in f if d > 1 and f[d]})
check("CONTROL D: 79-var automorphism inverse terminates at degree 3",
      nzA == [3], f"nonzero degrees {nzA}")

# CONTROL E: C^3 triangular x0+=x1^3, x1+=x2^3; closed form known
H3 = [[(Fraction(1), (1, 1, 1))], [(Fraction(1), (2, 2, 2))], []]
F3 = invert(H3, 3, 15)
nz3 = sorted({d for f in F3 for d in f if d > 1 and f[d]})
check("CONTROL E: C^3 triangular inverse terminates at degree 9, zeros 11-15",
      nz3 == [3, 5, 7, 9], f"nonzero degrees {nz3}")
# closed form: x1 = y1 - y2^3 ; x0 = y0 - (y1 - y2^3)^3
y0, y1, y2 = sp.symbols("y0 y1 y2")
closed = sp.expand(y0 - (y1 - y2**3)**3)
mine = sp.Integer(0)
for d, s in F3[0].items():
    for m, c in s.items():
        term = sp.Rational(c.numerator, c.denominator)
        for j in m:
            term *= (y0, y1, y2)[j]
        mine += term
check("CONTROL E: engine inverse == closed-form inverse (component x0)",
      sp.expand(mine - closed) == 0)

# ---------------------------------------------------------------- verdict
print()
if FAILURES:
    print("REFEREE-2 INDEPENDENT VERIFICATION: FAILURES:", FAILURES)
    sys.exit(1)
print("REFEREE-2 INDEPENDENT VERIFICATION: ALL CHECKS PASS")
print("  scope: compute layer only (parse, homogeneity, Keller spot-checks,")
print("  collision, inverse to deg 11, checkpoint agreement, composition,")
print("  controls). The representation-theoretic bridge is refereed separately")
print("  in audit/referee2-su79.md (report), not by this script.")
