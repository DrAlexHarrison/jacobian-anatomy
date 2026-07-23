#!/usr/bin/env python3
# 2026-07-20: the ladder FORCES P = S = 0
# (Q1-YES, d in {3,4}), via per-coefficient Rabinowitsch unit certificates in
# Singular. Normalized shape (justified by the one-line u^0 argument in
# notebook Entry 10 step (2)): A = a0(1+P u), B = b0(1+S u),
# f = c0 + phi x + f1(x) u (+ f2 u^2, d=3 only), phi*a0*b0 != 0.
# Run: /usr/bin/python3 < q1_forcing_check.py
# Expected: ALL UNIT x4, then "Q1 d in {3,4} MACHINE-FORCED".
# Verbatim the script run in-session; results logged in notebook Entry 10.
import sympy as sp, subprocess

x, u = sp.symbols('x u')

def run(cb):
    c_, b_ = cb
    d = c_ + b_
    degPS = 5 - d
    P = sum(sp.Symbol(f'p{i}')*x**i for i in range(degPS+1))
    S = sum(sp.Symbol(f's{i}')*x**i for i in range(degPS+1))
    f1 = sum(sp.Symbol(f'g{i}')*x**i for i in range(6-d+1))
    phi, f2c = sp.symbols('phi f2c')
    I1 = sp.expand(phi*((1+c_)*P + (1+b_)*S) + sp.diff(f1,x))
    I2 = sp.expand(phi*(1+b_+c_)*P*S + sp.diff(f1,x)*((1+c_)*P+(1+b_)*S)
                   - f1*(c_*sp.diff(P,x) + b_*sp.diff(S,x)))
    I3 = sp.expand(sp.diff(f1,x)*(1+b_+c_)*P*S
                   - f1*(c_*sp.diff(P,x)*S + b_*P*sp.diff(S,x))
                   - (2*f2c*(c_*sp.diff(P,x) + b_*sp.diff(S,x)) if d==3 else 0))
    eqs = []
    for E in ([I1, I2, I3] + ([sp.expand(2*f2c*(c_*sp.diff(P,x)*S + b_*P*sp.diff(S,x)))] if d==3 else [])):
        pe = sp.Poly(E, x)
        eqs += [sp.expand(co) for co in pe.all_coeffs() if co != 0]
    targets = [f'p{i}' for i in range(degPS+1)] + [f's{i}' for i in range(degPS+1)]
    gvars = [f'g{i}' for i in range(6-d+1)]
    vs = targets + gvars + ['phi'] + (['f2c'] if d==3 else []) + ['t', 't2']
    lines = [f'ring r = 0, ({", ".join(vs)}), dp;']
    lines.append('ideal L = ' + ','.join(str(e).replace('**','^') for e in eqs) + ', t*phi-1;')
    for tg in targets:
        lines.append(f'ideal J{tg} = L, t2*{tg}-1;')
        lines.append(f'if (reduce(1, groebner(J{tg})) == 0) {{ "UNIT {tg}"; }} else {{ "NONUNIT {tg}"; }}')
    lines.append('exit;')
    src = '\n'.join(lines) + '\n'
    p = subprocess.run(['/usr/bin/Singular','-q'], input=src, capture_output=True, text=True, timeout=280)
    outs = p.stdout.strip().splitlines()
    nonunit = [l for l in outs if l.startswith('NONUNIT')]
    print(f'(c,b)=({c_},{b_}) d={d}: {len(outs)} checks, '
          + ('ALL UNIT => P=S=0 forced (Q1-YES for this shape)' if not nonunit
         else 'NONUNIT!!! ' + ' '.join(nonunit)))
    return not nonunit

ok = all(run(cb) for cb in [(2,1),(1,2),(3,1),(1,3)])
print('Q1 d in {3,4} MACHINE-FORCED' if ok else 'JACKPOT-WATCH: forcing FAILED somewhere')
