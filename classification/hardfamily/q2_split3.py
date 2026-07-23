#!/usr/bin/env python3
# 2026-07-21: gauged S2 certificates. GAUGE LEMMA: post-composing
# (F2,F3) by the unimodular constant matrix [[1,0],[-rho,1]] stays in the
# family (weight -e components mix freely), preserves Keller, det0, and
# injectivity, and maps an S2 member (Chat = rho*Bhat) to one with
# Chat == 0. Hence S2 empty of nonconstant-A Keller members iff the
# Chat == 0 stratum is. Certificates: ladder with Chat = 0 substituted +
# t*k - 1 + t2*a_j - 1, per lhat coefficient.
import sympy as sp, subprocess

r, u0 = sp.symbols('r u0')

def cert(a_, e_, cap=1700):
    al = sp.Symbol('al')
    b1, b2, g1, g2 = sp.symbols('b1 b2 g1 g2')
    ajs = sp.symbols(f'a0:{a_+1}')
    Bj = [sp.Symbol(f'Bb{j}') for j in range(a_+2)]
    lhat = sum(ajs[j]*r**j for j in range(a_+1))
    Bhat = sum(Bj[j]*r**j for j in range(a_+2))
    A = al + u0*lhat
    B = b1 + b2*r + u0*Bhat
    C = g1 + g2*r
    k = al*(b1*g2 - b2*g1)
    lhs = sp.expand((B*sp.diff(C, r) - C*sp.diff(B, r)) *
                    (A*B + u0*(e_*lhat*B + a_*A*Bhat))
                    - u0*(-C*Bhat) *
                    (e_*u0*sp.diff(lhat, r)*B + a_*A*(b2 + u0*sp.diff(Bhat, r)))
                    - k*B)
    poly = sp.Poly(lhs, u0, r)
    eqs = sorted({sp.expand(c) for c in poly.coeffs()} - {sp.S(0)},
                 key=sp.default_sort_key)
    targets = [str(a) for a in ajs]
    vars_ = ([str(v) for v in ajs] + [str(v) for v in Bj] +
             ['al', 'b1', 'b2', 'g1', 'g2', 't', 't2'])
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
    print(f'(a,e)=({a_},{e_}) gauged-S2 (Chat=0): {verdict}', flush=True)
    return verdict

summary = {}
for pair in [(3, 1), (3, 2), (4, 1)]:
    summary[pair] = cert(*pair)
print('GAUGED-S2 SUMMARY:', summary, flush=True)
