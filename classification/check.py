#!/usr/bin/env python3
"""Quick independent check of the classification, needing only sympy.

Two checks, neither using Singular or any output stored in this
repository:
  1. For a sample of systems, rebuild the three collision ideals from
     jobs/fNNNN.sing and recompute their Groebner bases in sympy; the
     basis must be [1] (the system's family has no non-injective Keller
     member).
  2. Recompute the 741 -> 371 isomorphism-class fold from the job files'
     monomial supports.

Usage:
  python3 check.py            # default sample (12 systems) + the fold
  python3 check.py 26 210     # check specific systems instead
  python3 check.py --n 30     # larger sample

Systems marked Groebner-hard in LEDGER.md will not terminate here any
more than they did in Singular; the sample avoids them. Total runtime
for the default sample is under a minute on ordinary hardware.
"""
import sys, os, re, glob, time, itertools
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
JOBS = os.path.join(HERE, 'jobs')

HARD = set()
m = re.findall(r'f(\d{4})', open(os.path.join(HERE, 'LEDGER.md')).read())
# the ledger lists the hard systems by name; anything it names is excluded
HARD = {int(x) for x in m}

DEFAULT = [0, 17, 34, 51, 68, 85, 102, 119, 136, 153, 170, 187]

args = [a for a in sys.argv[1:]]
n = None
if '--n' in args:
    i = args.index('--n'); n = int(args[i + 1]); del args[i:i + 2]
fids = [int(a) for a in args] if args else DEFAULT
if n:
    pool = sorted(int(re.search(r'f(\d+)', p).group(1))
                  for p in glob.glob(os.path.join(JOBS, 'f*.sing')))
    fids = [f for f in pool if f not in HARD][::max(1, len(pool) // n)][:n]

def check_system(fid):
    job = open(os.path.join(JOBS, f'f{fid:04d}.sing')).read()
    vars_ = [v.strip() for v in
             re.search(r'ring r = 0, \(([^)]*)\), dp;', job).group(1).split(',')]
    syms = {v: sp.Symbol(v) for v in vars_}
    gens_syms = [syms[v] for v in vars_]
    def P(s): return sp.sympify(s.replace('^', '**'), locals=syms)
    ke = [g.strip() for g in re.search(r'ideal KE = (.*?);', job, re.S).group(1).split(',')]
    det0 = re.search(r'poly det0 = (.*?);', job, re.S).group(1).strip()
    coll = [g.strip() for g in re.search(r'ideal COLL = (.*?);', job, re.S).group(1).split(',')]
    for k, rab in enumerate(['t*(x-X)-1', 't*(y-Y)-1', 't*(z-Z)-1']):
        I = [P(g) for g in ke + coll] + [P(rab), P(f's*({det0})-1')]
        G = sp.groebner(I, *gens_syms, order='grevlex')
        if list(G.exprs) != [sp.Integer(1)]:
            return False
    return True

t0 = time.time()
bad = []
for fid in fids:
    ok = check_system(fid)
    print(f'f{fid:04d}: {"empty (Groebner basis = [1] for all three collision ideals)" if ok else "CHECK FAILED"}', flush=True)
    if not ok:
        bad.append(fid)
print(f'{len(fids) - len(bad)}/{len(fids)} sampled systems verified in {time.time()-t0:.0f}s')

# the class fold, from the job files themselves
def parse_supports(path):
    txt = open(path).read()
    comps = []
    for g in re.search(r'ideal COLL = (.*?);', txt, re.S).group(1).split(','):
        supp = set()
        for term in re.finditer(r'([+-]?)\s*(q\d+)((?:\*[xyz](?:\^\d+)?)*)', g):
            mono = term.group(3)
            if not mono:
                continue
            e = {'x': 0, 'y': 0, 'z': 0}
            for v in re.finditer(r'\*([xyz])(?:\^(\d+))?', mono):
                e[v.group(1)] += int(v.group(2) or 1)
            supp.add((e['x'], e['y'], e['z']))
        comps.append(frozenset(supp))
    return comps

def canon(sup):
    best = None
    for sig in itertools.permutations(range(3)):
        t = tuple(frozenset(tuple(a[sig[j]] for j in range(3)) for a in sup[sig[i]])
                  for i in range(3))
        key = tuple(tuple(sorted(s)) for s in t)
        if best is None or key < best:
            best = key
    return best

classes = {}
files = sorted(glob.glob(os.path.join(JOBS, 'f*.sing')))
for p in files:
    classes.setdefault(canon(parse_supports(p)), []).append(p)
print(f'class fold: {len(files)} systems -> {len(classes)} isomorphism classes')

if bad or len(classes) != 371 or len(files) != 741:
    print('RESULT: MISMATCH', bad)
    sys.exit(1)
print('RESULT: all checks pass')
