#!/usr/bin/env python3
import os
# 2026-07-20: isomorphism audit of the 88-system hard
# family under the T1/T2 equivalence group (tied source+target coordinate
# permutation; sign flip acts trivially on S-triples). Prints every
# isomorphic group WITH the witnessing permutation, verified by exact
# S-triple equality. Run: python3 iso_audit.py
import json
from itertools import permutations

BASE = os.path.dirname(os.path.abspath(__file__))
rep = json.load(open(f'{BASE}/families_deg6.json'))
mixed = [f for f in rep['families'] if f['class'] == 'MIXED']
mixed.sort(key=lambda f: sum(f['sizes']))

hard = [630,631,632,635,638,639,641,645,646,648,649,650,651,655,657,658,659,
        664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,
        681,682,683,684,685,687,688,689,691,692,694,695,696,698,699,700,702,
        703,704] + list(range(706,741))

def transported(f, perm):
    # apply tied permutation: component i gets support of component perm[i],
    # with exponent vectors permuted by the same perm
    S = [set(map(tuple, s)) for s in f['S']]
    return tuple(frozenset(tuple(al[perm[j]] for j in range(3)) for al in S[perm[i]])
                 for i in range(3))

def canon(f):
    best, bestp = None, None
    for perm in permutations(range(3)):
        sig = tuple(tuple(sorted(s)) for s in transported(f, perm))
        if best is None or sig < best:
            best, bestp = sig, perm
    return best

sigs = {}
for fid in hard:
    sigs.setdefault(canon(mixed[fid]), []).append(fid)

groups = sorted(v for v in sigs.values() if len(v) > 1)
singles = sorted(v[0] for v in sigs.values() if len(v) == 1)
print(f'{len(sigs)} distinct classes among {len(hard)} systems; '
      f'{len(groups)} nontrivial groups; singletons: {singles}')
for g in groups:
    base = mixed[g[0]]
    for other in g[1:]:
        fo = mixed[other]
        wit = None
        target = tuple(frozenset(map(tuple, s)) for s in base['S'])
        for perm in permutations(range(3)):
            if transported(fo, perm) == target:
                wit = perm
                break
        assert wit is not None, (g[0], other)
        # verify weight transport: w'[perm[j]] pattern up to sign
        wo, wb = fo['w'], base['w']
        ok_w = all(wo[wit[j]] == wb[j] for j in range(3)) or \
               all(wo[wit[j]] == -wb[j] for j in range(3))
        print(f'  f{other:04d} {tuple(fo["w"])} == f{g[0]:04d} {tuple(wb)} '
              f'via perm {wit} (S-triples EXACTLY equal; w transport '
              f'{"OK" if ok_w else "sign-mixed (S equality is the certificate)"})')
