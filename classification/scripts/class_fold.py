#!/usr/bin/env python3
# Fold the 741 systems into isomorphism classes under simultaneous
# source/target coordinate permutation, from the job files' own supports.
import os, re, glob, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
JOBS = os.path.join(HERE, '..', 'jobs')

def parse_coll(path):
    txt = open(path).read()
    m = re.search(r'ideal COLL = (.*?);', txt, re.S)
    comps = []
    for g in m.group(1).split(','):
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
    fid = int(re.search(r'f(\d+)\.sing', p).group(1))
    classes.setdefault(canon(parse_coll(p)), []).append(fid)
sizes = {}
for members in classes.values():
    sizes[len(members)] = sizes.get(len(members), 0) + 1
print(f"systems: {len(files)}")
print(f"classes: {len(classes)}")
print(f"class sizes (size: count): {sizes}")
selfpaired = [m for m in classes.values() if len(m) == 1]
print(f"singletons: {sorted(x for m in selfpaired for x in m)[:5]}")
