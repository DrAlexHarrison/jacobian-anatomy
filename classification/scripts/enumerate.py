#!/usr/bin/env python3
import os
"""P2-EQUIVARIANT enumeration (T2 of PROGRAM.md).

Canonical form (T1): v = w; component i is supported on
S_i(w) = {alpha in Delta_D : <w, alpha - e_i> = 0}, Delta_D = {alpha >= 0, |alpha| <= D}.

Complete list of maximal families: U = {primitive(d x d') : d, d' in R independent},
R = {alpha - e_i} \\ {0}.  Every equivariant Keller map of degree <= D lives (up to
source/target coordinate permutation) inside family(u) for some u in U   [T2].

Outputs families.json: deduped families with sign-class verdict routing, plus
machine-checked completeness certificates:
  - cert_b: for every primitive direction d0 of R, an explicit u in U with
    chain-family(d0) componentwise contained in family(u);
  - cert_identity: e_i in S_i for every family;
  - cert_incomparable is a theorem (rank-2 relation sets force parallel), not checked here.

Usage: /usr/bin/python3 enumerate.py [D]   (default D = 6)
"""
import itertools, json, math, sys

D = int(sys.argv[1]) if len(sys.argv) > 1 else 6
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"families_deg{D}.json")

E = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def prim(d):
    g = math.gcd(math.gcd(abs(d[0]), abs(d[1])), abs(d[2]))
    if g == 0: return None
    d = (d[0]//g, d[1]//g, d[2]//g)
    for c in d:
        if c > 0: return d
        if c < 0: return (-d[0], -d[1], -d[2])
    return None

Delta = [(a, b, c) for a in range(D+1) for b in range(D+1) for c in range(D+1)
         if a+b+c <= D]
assert len(Delta) == (D+1)*(D+2)*(D+3)//6

R = set()
for al in Delta:
    for e in E:
        d = sub(al, e)
        if d != (0, 0, 0):
            R.add(d)
Rdirs = sorted({prim(d) for d in R})          # primitive up-to-sign directions of R

# ---- the complete direction list U (case (c) of T2) ----------------------
U = set()
for d, dp in itertools.combinations(Rdirs, 2):
    u = prim(cross(d, dp))
    if u is not None:
        U.add(u)
U = sorted(U)

def family(u):
    return [sorted(al for al in Delta if dot(u, sub(al, e)) == 0) for e in E]

def classify(u, S):
    maxdeg = max(sum(al) for Si in S for al in Si)
    if maxdeg <= 1: return "LINEAR"
    if maxdeg <= 2: return "WANG"
    pos = sum(1 for c in u if c > 0); neg = sum(1 for c in u if c < 0)
    zer = 3 - pos - neg
    if zer == 0 and (pos == 3 or neg == 3): return "PROPER"
    if zer == 1 and (pos == 2 or neg == 2): return "PLANAR-SPLIT"
    if zer == 2: return "AXIS"
    return "MIXED"

# ---- canonical dedupe under coordinate permutations ----------------------
PERMS = list(itertools.permutations(range(3)))
def act(p, vec): return tuple(vec[p.index(i)] for i in range(3))
# act(p, .): new coordinate i carries old coordinate p^{-1}(i); consistent for
# weights and exponent vectors alike.  Component relabeling follows the same p.
def canon(u, S):
    # new component j = old component p^{-1}(j), monomials transformed by p;
    # verified: transformed triple == family(act(p, u)).
    best = None
    for p in PERMS:
        pu = act(p, u)
        pinv = [p.index(j) for j in range(3)]
        pS = [tuple(sorted(act(p, al) for al in S[pinv[j]])) for j in range(3)]
        key = (pu, tuple(pS))
        if best is None or key < best:
            best = key
    return best

fams = {}
for u in U:
    S = family(u)
    key = canon(u, S)
    if key not in fams:
        fams[key] = {
            "w": list(u),
            "S": [[list(a) for a in Si] for Si in S],
            "sizes": [len(Si) for Si in S],
            "maxdeg": max(sum(al) for Si in S for al in Si),
            "class": classify(u, S),
        }

# identity certificate: e_i in S_i always
for f in fams.values():
    for i in range(3):
        assert list(E[i]) in f["S"][i], ("identity cert failed", f["w"])

# ---- completeness certificate for case (b): chain families embed ---------
cert_b_fail = []
Uset = set(U)
for d0 in Rdirs:
    chain = [sorted(al for al in Delta
                    if (lambda dd: dd == (0,0,0) or prim(dd) == d0)(sub(al, e)))
             for e in E]
    ok = False
    for dp in Rdirs:
        if dp == d0: continue
        u = prim(cross(d0, dp))
        if u is None: continue
        S = family(u)
        if all(set(chain[i]) <= set(S[i]) for i in range(3)):
            ok = True
            break
    if not ok:
        cert_b_fail.append(d0)
assert not cert_b_fail, f"case-(b) embedding failed for {cert_b_fail[:5]}"

from collections import Counter
cnt = Counter(f["class"] for f in fams.values())
report = {
    "degree_cap": D,
    "n_Delta": len(Delta),
    "n_R_dirs": len(Rdirs),
    "n_cross_directions": len(U),
    "n_families_after_dedupe": len(fams),
    "class_counts": dict(cnt),
    "cert_identity": "PASS",
    "cert_b_chain_embedding": f"PASS ({len(Rdirs)} directions)",
    "families": sorted(fams.values(), key=lambda f: (f["class"] != "MIXED",
                                                     sum(f["sizes"]))),
}
with open(OUT, "w") as fh:
    json.dump(report, fh, indent=1)

print(f"Delta_{D}: {len(Delta)} pts | R dirs: {len(Rdirs)} | cross dirs: {len(U)}"
      f" | families: {len(fams)}")
print("classes:", dict(cnt))
mixed = [f for f in fams.values() if f["class"] == "MIXED"]
mixed.sort(key=lambda f: sum(f["sizes"]))
print(f"MIXED systems: {len(mixed)}; coefficient counts "
      f"min={sum(mixed[0]['sizes']) if mixed else 0} "
      f"max={sum(mixed[-1]['sizes']) if mixed else 0}")
for f in mixed[-10:]:
    print("  big:", f["w"], f["sizes"], "maxdeg", f["maxdeg"])
print("wrote", OUT)
