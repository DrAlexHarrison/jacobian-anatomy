#!/usr/bin/env python3
# 2026-07-22: branch beta (Bc=0, Ac !== 0) certificates + sub-cases,
# plus direct certs for gamma (Cc = rho*Bc) and a3 (Ac=0, Cc=rho*Bc).
# Usage: python3 < q11_beta.py [inline sub-case certs]
#        python3 < q11_beta.py FULL  [the three big certs; run detached]
import sympy as sp, subprocess, sys

r, u0 = sp.symbols('r u0')
al, rho = sp.Symbol('al'), sp.Symbol('rho')
N = sp.sympify(open('q8_layers.srepr').read())
S = sp.Symbol

def subs_zero(names):
    d = {}
    for nm, deg in names:
        for j in range(deg+1):
            d[S(f'{nm}{j}')] = 0
    return d

def cert(tag, eqs, vars_, hyps, targets, cap=600):
    coeffs = []
    for e in eqs:
        e = sp.expand(e)
        if e == 0:
            continue
        for x in sp.Poly(e, r).all_coeffs():
            if x == 0:
                continue
            num, _ = sp.fraction(sp.together(x))
            coeffs.append(sp.expand(num))
    svars = [f's{i}' for i in range(len(hyps))]
    L = [f'ring rr = 0, ({", ".join(vars_ + svars + ["tt"])}), dp;']
    gens = [str(x).replace('**', '^') for x in coeffs]
    gens += [f'{svars[i]}*({str(hp).replace("**", "^")})-1' for i, hp in enumerate(hyps)]
    L.append('ideal LD = ' + ','.join(gens) + ';')
    for name, expr in targets.items():
        es = str(sp.expand(expr)).replace('**', '^')
        L.append(f'ideal J{name} = LD, tt*({es})-1;')
        L.append(f'if (reduce(1, groebner(J{name})) == 0) {{ "UNIT {name}"; }} else {{ "NONUNIT {name}"; }}')
    L.append('exit;')
    try:
        p = subprocess.run(['/usr/bin/Singular', '-q'], input='\n'.join(L) + '\n',
                           capture_output=True, text=True, timeout=cap)
        out = ' | '.join(l for l in p.stdout.strip().splitlines() if l.strip())
    except subprocess.TimeoutExpired:
        out = f'TIMEOUT {cap}s (UNDECIDED)'
    print(f'CERT {tag}: {out}', flush=True)

zbc = subs_zero([('bc', 3)])
AC = {f'ac{j}': S(f'ac{j}') for j in range(3)}

if len(sys.argv) > 1 and sys.argv[1] == 'FULL':
    # beta monolithic: Bc=0, targets ac0..2 (14 vars)
    eqs = [N[m].subs(zbc) for m in range(2, 9)]
    cert('b-FULL', eqs,
         ['bh0', 'bh1', 'bh2', 'ch0', 'ch1', 'ch2', 'cc0', 'cc1', 'cc2', 'cc3',
          'ac0', 'ac1', 'ac2', 'al'], [al], AC, cap=10800)
    # gamma direct: Cc = rho*Bc, targets ac0..2 (15 vars)
    zg = {S(f'cc{j}'): rho*S(f'bc{j}') for j in range(4)}
    eqs = [N[m].subs(zg) for m in range(2, 9)]
    cert('g-RHO', eqs,
         ['bh0', 'bh1', 'bh2', 'ch0', 'ch1', 'ch2', 'bc0', 'bc1', 'bc2', 'bc3',
          'ac0', 'ac1', 'ac2', 'rho', 'al'], [al], AC, cap=10800)
    # a3 direct: Ac=0, Cc = rho*Bc, targets P0, P1 (12 vars)
    zac = subs_zero([('ac', 2)])
    eqs = [N[m].subs(zac).subs(zg) for m in range(2, 9)]
    cert('a3-RHO', eqs,
         ['bh0', 'bh1', 'bh2', 'ch0', 'ch1', 'ch2', 'bc0', 'bc1', 'bc2', 'bc3',
          'rho', 'al'], [al],
         {'P0': 2*S('bh0') + S('ch1'), 'P1': S('bh1') + 2*S('ch2')}, cap=10800)
    sys.exit(0)

# ---- inline sub-case certs for beta
zbh = subs_zero([('bh', 2)])
zcc = subs_zero([('cc', 3)])
# b-h0: Bh=0
eqs = [N[m].subs(zbc).subs(zbh) for m in range(2, 9)]
cert('b-h0', eqs, ['ch0', 'ch1', 'ch2', 'cc0', 'cc1', 'cc2', 'cc3',
                   'ac0', 'ac1', 'ac2', 'al'], [al], AC)
# b-K0: Cc=0
eqs = [N[m].subs(zbc).subs(zcc) for m in range(2, 9)]
cert('b-K0', eqs, ['bh0', 'bh1', 'bh2', 'ch0', 'ch1', 'ch2',
                   'ac0', 'ac1', 'ac2', 'al'], [al], AC)
# b-consts: h=bh0, K=cc0, Ac=ac0 (nonzero consts), c generic
zc = {S('bh1'): 0, S('bh2'): 0, S('cc1'): 0, S('cc2'): 0, S('cc3'): 0,
      S('ac1'): 0, S('ac2'): 0}
eqs = [N[m].subs(zbc).subs(zc) for m in range(2, 9)]
cert('b-consts', eqs, ['bh0', 'ch0', 'ch1', 'ch2', 'cc0', 'ac0', 'al'],
     [al, S('bh0'), S('cc0')], {'ac0': S('ac0')})
# b-coinc: h=h2 r^2, Ac=a2 r^2, K=k3 r^3 (h2,k3 != 0), c generic
zk = {S('bh0'): 0, S('bh1'): 0, S('cc0'): 0, S('cc1'): 0, S('cc2'): 0,
      S('ac0'): 0, S('ac1'): 0}
eqs = [N[m].subs(zbc).subs(zk) for m in range(2, 9)]
cert('b-coinc', eqs, ['bh2', 'ch0', 'ch1', 'ch2', 'cc3', 'ac2', 'al'],
     [al, S('bh2'), S('cc3')], {'ac2': S('ac2')})
