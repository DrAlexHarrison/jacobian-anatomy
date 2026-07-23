#!/usr/bin/env python3
import os
# 2026-07-21: transfer-of-verdict machinery
# 00:10 ruling. Items delivered:
#  (1) transfer-lemma machine spot-check on load-bearing pairs (random
#      symbolic member of family B, transported by sigma: membership in A's
#      support + exact det J transport + det0 equality);
#  (2) per-pair witness checked against GROUND TRUTH job files (the .sing the
#      certificates consumed: scratch a/bNNNN.sing for the decomposed stream,
#      programs jobs/fNNNN.sing for the monolith stream) — NOT enumerate.py
#      output;
#  (3) negative controls: a wrong sigma on a true pair must FAIL; a non-pair
#      must admit NO sigma.
# Run: /usr/bin/python3 < iso_transfer.py
import re, os, random
from itertools import permutations
import sympy as sp

SCRATCH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'hardfamily')
JOBS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'jobs')

def source_file(fid):
    for p in (f'{SCRATCH}/b{fid:04d}.sing', f'{SCRATCH}/a{fid:04d}.sing',
              f'{JOBS}/f{fid:04d}.sing'):
        if os.path.exists(p):
            return p
    raise FileNotFoundError(fid)

def parse_supports(fid):
    """Extract the 3 components' lowercase monomial supports from the COLL
    ideal of the job file the certificates consumed."""
    txt = open(source_file(fid)).read()
    m = re.search(r'ideal COLL = (.*?);', txt, re.S)
    comps = []
    depth_split = m.group(1).split(',')
    assert len(depth_split) == 3, (fid, len(depth_split))
    for comp in depth_split:
        supp = set()
        # lowercase family terms: q<j>*<monomial in x,y,z> or <mono>*q? format is qJ last? observed: q0*z etc. and x-mono after q
        for tok in re.findall(r'q\d+(?:\*[a-z][a-z0-9^*]*)?', comp):
            expo = [0, 0, 0]
            for var, p in re.findall(r'([xyz])(?:\^(\d+))?', tok):
                expo['xyz'.index(var)] = int(p) if p else 1
            if tok.count('*'):  # has a monomial part
                supp.add(tuple(expo))
            else:
                supp.add((0, 0, 0))
        comps.append(frozenset(supp))
    return comps

def transported(S, perm):
    return tuple(frozenset(tuple(al[perm[j]] for j in range(3)) for al in S[perm[i]])
                 for i in range(3))

def find_sigma(fid_a, fid_b):
    """Return the sigma matching B's ground-truth supports onto A's, or None."""
    SA = tuple(parse_supports(fid_a))
    SB = parse_supports(fid_b)
    for perm in permutations(range(3)):
        if transported(SB, perm) == SA:
            return perm
    return None

def check_pair(fid_a, fid_b, expect=True):
    sigma = find_sigma(fid_a, fid_b)
    ok = (sigma is not None) == expect
    print(f'  pair (f{fid_a:04d}, f{fid_b:04d}): sigma={sigma} '
          f'[sources: {os.path.basename(source_file(fid_a))}, '
          f'{os.path.basename(source_file(fid_b))}] '
          f'{"OK" if ok else "*** UNEXPECTED ***"}')
    return sigma if ok else 'FAIL'

def lemma_spotcheck(fid_a, fid_b):
    """Random symbolic member of B, transported: membership in A + det J
    transport + det0 equality, exact."""
    sigma = find_sigma(fid_a, fid_b)
    assert sigma is not None
    SA = parse_supports(fid_a)
    SB = parse_supports(fid_b)
    xs = sp.symbols('x y z')
    rnd = random.Random(720721)
    FB = []
    for i in range(3):
        FB.append(sum(rnd.randint(1, 9) * xs[0]**a * xs[1]**b * xs[2]**c
                      for (a, b, c) in sorted(SB[i])))
    # G_i(x) = F_B[sigma[i]] with monomial x^alpha |-> prod_j x_j^{alpha_{sigma(j)}};
    # realized by substitution x_k -> x_{sigma^{-1}(k)} (inverse, NOT sigma:
    # a direction bug here passes silently on involutions and fails on 3-cycles).
    isig = [0, 0, 0]
    for j in range(3):
        isig[sigma[j]] = j
    sub = {xs[k]: xs[isig[k]] for k in range(3)}
    G = [sp.expand(FB[sigma[i]].subs(sub, simultaneous=True)) for i in range(3)]
    # membership: supp(G_i) subseteq SA[i]
    memb = True
    for i in range(3):
        for mono in sp.Poly(G[i], *xs).monoms():
            if tuple(mono) not in SA[i]:
                memb = False
    # det transport: det JG(x) == det JF_B evaluated at permuted variables
    dJB = sp.Matrix(FB).jacobian(xs).det()
    dJG = sp.Matrix(G).jacobian(xs).det()
    ok_det = sp.expand(dJG - dJB.subs(sub, simultaneous=True)) == 0
    # det0 equality
    ok_d0 = sp.simplify(dJG.subs({v: 0 for v in xs}) - dJB.subs({v: 0 for v in xs})) == 0
    print(f'  lemma spot-check (f{fid_a:04d} <- f{fid_b:04d}, sigma={sigma}): '
          f'membership {"PASS" if memb else "FAIL"}, detJ transport '
          f'{"PASS" if ok_det else "FAIL"}, det0 equal {"PASS" if ok_d0 else "FAIL"}')
    return memb and ok_det and ok_d0

print('(2) per-pair witnesses vs ground-truth job files, the load-bearing set:')
pairs = [(720, 723), (721, 722), (728, 729), (734, 735), (736, 737),
         (738, 739), (641, 633), (646, 644), (635, 654), (726, 727)]
sigmas = {}
allok = True
for a, b in pairs:
    s = check_pair(a, b)
    allok &= (s != 'FAIL')
    sigmas[(a, b)] = s

print('(1) transfer-lemma spot-checks (GL2 straggler / zero-weight / easy-hard crossing):')
for a, b in [(720, 723), (736, 737), (641, 633)]:
    allok &= lemma_spotcheck(a, b)

print('(3) negative controls:')
# The matching-sigma set must be exactly ONE coset of Aut(A): equivalently,
# (a) at least one sigma FAILS (checker can say no), and (b) #matches = #Aut(A).
# (Families with equal repeated weights carry genuine self-symmetries, e.g.
# Aut(720) = {id, x<->z}; demanding a UNIQUE sigma would be wrong.)
SA = parse_supports(720); SB = parse_supports(723)
match = [p for p in permutations(range(3)) if transported(SB, p) == tuple(SA)]
auts = [p for p in permutations(range(3)) if transported(SA, p) == tuple(SA)]
neg1 = (0 < len(match) < 6) and (len(match) == len(auts))
print(f'  (720,723): matches {match} = one Aut-coset (|Aut|={len(auts)}), '
      f'{6-len(match)}/6 sigmas rejected: {"PASS" if neg1 else "FAIL"}')
# non-isomorphic pair admits no sigma
neg2 = find_sigma(630, 631) is None
print(f'  non-pair (630,631) rejected: {"PASS" if neg2 else "FAIL"}')
allok &= neg1 and neg2
print('ALL TRANSFER CHECKS PASS' if allok else '*** FAILURES ***')
