#!/usr/bin/env python3
"""Collect sweep results into verdicts.json + a markdown verdict table.

Reads jobs*/manifest.json + results*/f*.out. Statuses:
  EMPTY     - all three unit certificates: no equivariant Keller counterexample
  NONEMPTY  - some run non-unit: WITNESS HUNT REQUIRED (shout!)
  TIMEOUT   - no VERDICT line: needs second pass
Usage: /usr/bin/python3 collect_results.py [jobsdir resultsdir outprefix]
"""
import glob, json, os, sys

jobsdir = sys.argv[1] if len(sys.argv) > 1 else 'jobs'
resdir = sys.argv[2] if len(sys.argv) > 2 else 'results'
outpfx = sys.argv[3] if len(sys.argv) > 3 else 'verdicts'
import os
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
man = {m['fid']: m for m in json.load(open(f'{base}/{jobsdir}/manifest.json'))}

verdicts = {}
for fid, m in man.items():
    out = f'{base}/{resdir}/f{fid:04d}.out'
    if not os.path.exists(out):
        verdicts[fid] = 'PENDING'
        continue
    txt = open(out).read()
    if 'EMPTY (no equivariant Keller counterexample' in txt:
        verdicts[fid] = 'EMPTY'
    elif 'NONEMPTY' in txt or 'NONUNIT' in txt:
        verdicts[fid] = 'NONEMPTY'
    else:
        verdicts[fid] = 'TIMEOUT'

from collections import Counter
cnt = Counter(verdicts.values())
rows = []
for fid in sorted(man):
    m = man[fid]
    rows.append({'fid': fid, 'w': m['w'], 'sizes': m['sizes'],
                 'maxdeg': m['maxdeg'], 'ncoef': m['ncoef'],
                 'verdict': verdicts[fid]})
json.dump({'counts': dict(cnt), 'rows': rows},
          open(f'{base}/{outpfx}.json', 'w'), indent=1)
print('counts:', dict(cnt))
bad = [r for r in rows if r['verdict'] == 'NONEMPTY']
for r in bad:
    print('!! NONEMPTY:', r)
to = [r for r in rows if r['verdict'] in ('TIMEOUT', 'PENDING')]
print(f'{len(to)} timeout/pending:', [(r['fid'], r['w'], r['ncoef']) for r in to[:20]])
