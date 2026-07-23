#!/usr/bin/env python3
# 2026-07-21: S2-rho certificates. Branch S2 of (L4) has
# W(Bhat,Chat) = 0; with Bhat != 0 (Bhat == 0 is branch S1, already ALL
# UNIT at (3,1)) linear dependence gives Chat = rho*Bhat for a scalar rho.
# Soundness: any ladder point with lhat != 0, W2 = 0, Bhat != 0 lifts to a
# solution of the rho-substituted system; so UNIT certificates on the
# substituted system + t2*a_j - 1 empty the branch.
import sympy as sp, subprocess

r, u0 = sp.symbols('r u0')

def cert(a_, e_, cap=1700):
    al, rho = sp.symbols('al rho')
    b1, b2, g1, g2 = sp.symbols('b1 b2 g1 g2')
    ajs = sp.symbols(f'a0:{a_+1}')
    Bj = [sp.Symbol(f'Bb{j}') for j in range(a_+2)]
    lhat = sum(ajs[j]*r**j for j in range(a_+1))
    Bhat = sum(Bj[j]*r**j for j in range(a_+2))
    Chat = rho*Bhat
    A = al + u0*lhat
    B = b1 + b2*r + u0*Bhat
    C = g1 + g2*r + u0*Chat
    k = al*(b1*g2 - b2*g1)
    lhs = sp.expand((B*sp.diff(C, r) - C*sp.diff(B, r)) *
                    (A*B + u0*(e_*lhat*B + a_*A*Bhat))
                    - u0*(B*Chat - C*Bhat) *
                    (e_*u0*sp.diff(lhat, r)*B + a_*A*(b2 + u0*sp.diff(Bhat, r)))
                    - k*B)
    poly = sp.Poly(lhs, u0, r)
    eqs = sorted({sp.expand(c) for c in poly.coeffs()} - {sp.S(0)},
                 key=sp.default_sort_key)
    targets = [str(a) for a in ajs]
    vars_ = ([str(v) for v in ajs] + [str(v) for v in Bj] +
             ['rho', 'al', 'b1', 'b2', 'g1', 'g2', 't', 't2'])
    L = [f'ring rr = 0, ({", ".join(vars_)}), dp;']
    L.append('ideal LD = ' + ','.join(str(e).replace('**', '^') for e in eqs) +
             ', t*al*b1*g2-t*al*b2*g1-1;')
    for tg in targets:
        L.append(f'ideal J{tg} = LD, t2*{tg}-1;')
        L.append(f'if (reduce(1, groebner(J{tg})) == 0) {{ "UNIT {tg}"; }} '
                 f'else {{ "NONUNIT {tg}"; }}')
    L.append('exit;')
    src = '\n'.join(L) + '\n'
    try:
        p = subprocess.run(['/usr/bin/Singular', '-q'], input=src,
                           capture_output=True, text=True, timeout=cap)
        outs = [l for l in p.stdout.strip().splitlines() if l.strip()]
        nonunit = [l for l in outs if l.startswith('NONUNIT')]
        verdict = ('ALL UNIT' if (outs and not nonunit) else
                   ('NONUNIT: ' + ' '.join(nonunit) if nonunit else 'NO-OUTPUT'))
    except subprocess.TimeoutExpired:
        verdict = f'TIMEOUT {cap}s'
    print(f'(a,e)=({a_},{e_}) branch S2-rho: {verdict}', flush=True)
    return verdict

summary = {}
for pair in [(3, 1), (3, 2), (4, 1)]:
    summary[pair] = cert(*pair)
print('S2RHO SUMMARY:', summary, flush=True)
