#!/usr/bin/env python3
"""Per-MIXED-system machinery: generic family -> Keller equations -> Singular jobs.

For family (w, S1,S2,S3):  F_i = sum_{alpha in S_i} q_j x^alpha  (q_j fresh).
KellerEqs = all (x,y,z)-coefficients of det J - det J|_{x=0}   (polys in q).
det0      = det J|_{x=0} = det(linear part)                     (poly in q).
Coll      = F(x,y,z) - F(X,Y,Z).
Three Singular runs (k = x,y,z):  I_k = KellerEqs + Coll + <t*(p_k - q_k) - 1,
s*det0 - 1>.  reduce(1, groebner(I_k)) == 0 for all k  <=>  no C*-equivariant
Keller counterexample with support in this system (Nullstellensatz over C).

Self-check per system: det J is quasi-homogeneous of tau-weight 0 (every
monomial x^beta in det J has <w,beta> = 0) -- validates the family build.

Usage: /usr/bin/python3 gen_singular.py families_deg6.json jobsdir
"""
import json, os, sys
import sympy as sp

fam_file = sys.argv[1]
jobsdir = sys.argv[2]
os.makedirs(jobsdir, exist_ok=True)

x, y, z, X, Y, Z, t, s = sp.symbols('x y z X Y Z t s')
P = (x, y, z); Q = (X, Y, Z)

rep = json.load(open(fam_file))
mixed = [f for f in rep['families'] if f['class'] == 'MIXED']
mixed.sort(key=lambda f: sum(f['sizes']))

manifest = []
for fid, f in enumerate(mixed):
    w = f['w']; S = f['S']
    qs = []
    F = []
    Fq = []
    k = 0
    for i in range(3):
        comp = 0; compq = 0
        for al in S[i]:
            q = sp.Symbol(f'q{k}'); k += 1; qs.append(q)
            mono = x**al[0] * y**al[1] * z**al[2]
            comp += q * mono
            compq += q * X**al[0] * Y**al[1] * Z**al[2]
        F.append(sp.expand(comp)); Fq.append(sp.expand(compq))
    detJ = sp.expand(sp.Matrix(F).jacobian(P).det())
    # self-check: quasi-homogeneous of weight 0
    pd = sp.Poly(detJ, x, y, z)
    for mono in pd.monoms():
        wt = w[0]*mono[0] + w[1]*mono[1] + w[2]*mono[2]
        assert wt == 0, (f['w'], mono, wt, "detJ not weight-0: family build broken")
    det0 = detJ.subs({x: 0, y: 0, z: 0})
    eqs = [c for m, c in zip(pd.monoms(), pd.coeffs()) if m != (0, 0, 0)]
    coll = [sp.expand(F[i] - Fq[i]) for i in range(3)]

    def sing(e):
        return str(e).replace('**', '^')

    lines = []
    nv = len(qs)
    lines.append(f'// system {fid}: w={w} sizes={f["sizes"]} maxdeg={f["maxdeg"]}')
    qnames = ', '.join(f'q{j}' for j in range(nv))
    lines.append(f'ring r = 0, ({qnames}, x, y, z, X, Y, Z, t, s), dp;')
    lines.append('ideal KE = ' + ','.join(sing(e) for e in eqs) + ';')
    lines.append('poly det0 = ' + sing(det0) + ';')
    lines.append('ideal COLL = ' + ','.join(sing(e) for e in coll) + ';')
    lines.append('int allunit = 1;')
    for kk, (pk, qk) in enumerate(zip(P, Q)):
        lines.append(f'ideal I{kk} = KE, COLL, t*({pk}-{qk})-1, s*det0-1;')
        lines.append(f'ideal G{kk} = groebner(I{kk});')
        lines.append(f'if (reduce(1, G{kk}) != 0) {{ allunit = 0; '
                     f'"NONUNIT k={kk} system {fid}"; }}')
    lines.append(f'if (allunit == 1) {{ "VERDICT system {fid}: EMPTY '
                 f'(no equivariant Keller counterexample in this family)"; }}')
    lines.append(f'else {{ "VERDICT system {fid}: NONEMPTY -- WITNESS HUNT REQUIRED"; }}')
    lines.append('exit;')
    with open(f'{jobsdir}/f{fid:04d}.sing', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    manifest.append({'fid': fid, 'w': w, 'sizes': f['sizes'],
                     'maxdeg': f['maxdeg'], 'ncoef': nv, 'neqs': len(eqs)})

with open(f'{jobsdir}/manifest.json', 'w') as fh:
    json.dump(manifest, fh, indent=1)
print(f'wrote {len(manifest)} jobs to {jobsdir}; '
      f'largest ncoef={max(m["ncoef"] for m in manifest)}, '
      f'neqs max={max(m["neqs"] for m in manifest)}')
