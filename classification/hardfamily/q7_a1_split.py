#!/usr/bin/env python3
# 2026-07-22: {734,735} normalized + N1-eliminated + N8-split certs.
# Ah is enslaved by N1: Ah = -al*(2*Bh + Ch' - r*Bh')/2. Substituted in.
# Q2-a1 target: P := 2*Bh + Ch' - r*Bh' == 0 (3 coeffs) and Ac == 0 (3 coeffs).
# N8 factors 5*Ac*Bc*sigma(Bc,Cc): branches T1 (Bc == 0), T2 (Ac == 0,
# targets P only), T3 (Cc = rho*Bc, sound lift).
import sympy as sp, subprocess

r, u0 = sp.symbols('r u0')
al, rho = sp.symbols('al rho')
def poly(name, d):
    return sum(sp.Symbol(f'{name}{j}')*r**j for j in range(d+1))

def build(branch):
    Bh, Ch = poly('bh', 2), poly('ch', 2)
    Ac = poly('ac', 2)
    if branch == 'T1':
        Bc = sp.S(0); Cc = poly('cc', 3)
    elif branch == 'T2':
        Bc, Cc, Ac = poly('bc', 3), poly('cc', 3), sp.S(0)
    else:  # T3
        Bc = poly('bc', 3); Cc = rho*Bc
    Ah = sp.expand(-al*(2*Bh + sp.diff(Ch, r) - r*sp.diff(Bh, r))/2)
    A = al + u0*Ah + u0**2*Ac
    B = 1 + u0*Bh + u0**2*Bc
    C = r + u0*Ch + u0**2*Cc
    I = sp.expand((B*sp.diff(C,r) - C*sp.diff(B,r))*(A*B + u0*sp.diff(A,u0)*B + u0*A*sp.diff(B,u0))
                  - (B*sp.diff(C,u0) - C*sp.diff(B,u0))*u0*(sp.diff(A,r)*B + A*sp.diff(B,r))
                  - al*B)
    P = sp.Poly(sp.expand(2*I), u0, r)   # clear the 1/2
    eqs = sorted({sp.expand(c) for c in P.coeffs()} - {sp.S(0)}, key=sp.default_sort_key)
    Pcomb = sp.expand(2*Bh + sp.diff(Ch, r) - r*sp.diff(Bh, r))
    pc = [sp.expand(c) for c in sp.Poly(Pcomb, r).all_coeffs()]
    return eqs, pc

def cert(branch, cap=1700):
    eqs, pc = build(branch)
    targets = {}
    if branch != 'T2':
        for j in range(3):
            targets[f'ac{j}'] = sp.Symbol(f'ac{j}')
    for i, c in enumerate(pc):
        targets[f'P{i}'] = c
    base_vars = ['bh0','bh1','bh2','ch0','ch1','ch2','al','t','t2']
    if branch == 'T1':
        vars_ = ['cc0','cc1','cc2','cc3','ac0','ac1','ac2'] + base_vars
    elif branch == 'T2':
        vars_ = ['bc0','bc1','bc2','bc3','cc0','cc1','cc2','cc3'] + base_vars
    else:
        vars_ = ['bc0','bc1','bc2','bc3','ac0','ac1','ac2','rho'] + base_vars
    L = [f'ring rr = 0, ({", ".join(vars_)}), dp;']
    L.append('ideal LD = ' + ','.join(str(e).replace('**','^') for e in eqs) + ', t*al-1;')
    outs_names = []
    for name, expr in targets.items():
        es = str(sp.expand(expr)).replace('**','^')
        L.append(f'ideal J{name} = LD, t2*({es})-1;')
        L.append(f'if (reduce(1, groebner(J{name})) == 0) {{ "UNIT {name}"; }} else {{ "NONUNIT {name}"; }}')
        outs_names.append(name)
    L.append('exit;')
    src = '\n'.join(L) + '\n'
    try:
        p = subprocess.run(['/usr/bin/Singular','-q'], input=src, capture_output=True,
                           text=True, timeout=cap)
        outs = [l for l in p.stdout.strip().splitlines() if l.strip()]
        nonunit = [l for l in outs if l.startswith('NONUNIT')]
        v = ('ALL UNIT' if (outs and not nonunit and len(outs) == len(outs_names))
             else ('NONUNIT: ' + ' '.join(nonunit) if nonunit else f'PARTIAL {len(outs)}/{len(outs_names)}: ' + ' '.join(outs)))
    except subprocess.TimeoutExpired:
        v = f'TIMEOUT {cap}s'
    print(f'branch {branch}: {v}', flush=True)
    return v

summary = {}
for br in ('T1', 'T2', 'T3'):
    summary[br] = cert(br)
print('A1-SPLIT SUMMARY:', summary, flush=True)
