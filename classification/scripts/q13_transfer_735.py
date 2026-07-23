#!/usr/bin/env python3
import os
# 2026-07-22: EMPTY-BY-TRANSFER witness for f0735 from f0734,
# ground truth = the JOB FILES' actual monomial supports (certificate-consumed
# artifacts), NOT enumerate.py.
# w(734) = (1,-1,-1); w(735) = (1,-1,1).
# Step 1 (sign flip, no transfer needed): S(w) = S(-w) exactly, so family 735
#   equals the family of -w735 = (-1,1,-1) verbatim (same supports).
# Step 2 (permutation): sigma = transposition (1 2) maps -w735 = (-1,1,-1) to
#   (1,-1,-1) = w734. Verify the support bijection: transporting f0735's
#   supports by sigma (simultaneous source-var + target-component + weight
#   permutation) yields exactly f0734's supports.
import re

def parse_coll(path):
    """Extract per-component monomial supports (exponent triples in x,y,z)
    from the lowercase half of each COLL generator in a job file."""
    txt = open(path).read()
    m = re.search(r'ideal COLL = (.*?);', txt, re.S)
    gens = m.group(1).split(',')
    comps = []
    for g in gens:
        supp = set()
        for term in re.finditer(r'([+-]?)\s*(q\d+)((?:\*[xyz](?:\^\d+)?)*)', g):
            mono = term.group(3)
            if not mono:
                continue
            e = {'x': 0, 'y': 0, 'z': 0}
            for v in re.finditer(r'\*([xyz])(?:\^(\d+))?', mono):
                e[v.group(1)] += int(v.group(2) or 1)
            supp.add((e['x'], e['y'], e['z']))
        comps.append(supp)
    return comps

J = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'jobs')
S734 = parse_coll(f'{J}/f0734.sing')
S735 = parse_coll(f'{J}/f0735.sing')
print('sizes 734:', [len(s) for s in S734], ' 735:', [len(s) for s in S735])

# sigma = (1 2): source vars x<->y, target components 1<->2, weights permuted.
# Transport of a support: component i of the transported family :=
#   { alpha o sigma^{-1} : alpha in S_735[sigma(i)] }, sigma = (0 1) 0-indexed.
sig = [1, 0, 2]
T = [None]*3
for i in range(3):
    T[i] = {(a[sig[0]], a[sig[1]], a[sig[2]]) for a in S735[sig[i]]}
ok = all(T[i] == S734[i] for i in range(3))
print('SUPPORT BIJECTION sigma=(1 2) [S735 transported == S734]:',
      'PASS' if ok else 'FAIL')

# weight check: sigma applied to -w735 = (-1,1,-1) gives (1,-1,-1) = w734
w735 = (1, -1, 1)
neg = tuple(-w for w in w735)
perm = (neg[sig[0]], neg[sig[1]], neg[sig[2]])
print('WEIGHT CHECK sigma(-w735) == w734:', 'PASS' if perm == (1, -1, -1) else 'FAIL')

# negative control: the identity permutation must NOT match
Tid = S735
print('NEGATIVE CONTROL (identity perm rejected):',
      'PASS' if not all(Tid[i] == S734[i] for i in range(3)) else 'FAIL')
# negative control 2: transposition (2 3) must NOT match
sig2 = [0, 2, 1]
T2 = [{(a[sig2[0]], a[sig2[1]], a[sig2[2]]) for a in S735[sig2[i]]} for i in range(3)]
print('NEGATIVE CONTROL (perm (2 3) rejected):',
      'PASS' if not all(T2[i] == S734[i] for i in range(3)) else 'FAIL')
