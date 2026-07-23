#!/usr/bin/env python3
# 2026-07-21: branch-split Q2 certificates. (L4) reads
# lhat*W(Bhat,Chat)*Bhat = 0 in C[r] (domain), so with a_j != 0:
#   branch S1: Bhat == 0 (coefficients added), or
#   branch S2: W(Bhat,Chat) == 0 (its r-coefficients added).
# Both branch certs UNIT for every lhat coefficient a_j  =>  A const forced
# (Q2-YES) for the class  =>  EMPTY via slice-Moh (Entry 24 Branch I).
import sympy as sp, subprocess

r, u0 = sp.symbols('r u0')

def build(a_, e_):
    al = sp.Symbol('al')
    b1, b2, g1, g2 = sp.symbols('b1 b2 g1 g2')
    ajs = sp.symbols(f'a0:{a_+1}')
    Bj = [sp.Symbol(f'Bb{j}') for j in range(a_+2)]
    Cj = [sp.Symbol(f'Cc{j}') for j in range(a_+2)]
    lhat = sum(ajs[j]*r**j for j in range(a_+1))
    Bhat = sum(Bj[j]*r**j for j in range(a_+2))
    Chat = sum(Cj[j]*r**j for j in range(a_+2))
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
    W2 = sp.expand(Bhat*sp.diff(Chat, r) - Chat*sp.diff(Bhat, r))
    w2c = [sp.expand(c) for c in sp.Poly(W2, r).all_coeffs() if c != 0]
    return eqs, w2c, ajs, Bj, Cj

def cert(a_, e_, branch, cap=1800):
    eqs, w2c, ajs, Bj, Cj = build(a_, e_)
    targets = [str(a) for a in ajs]
    vars_ = ([str(v) for v in ajs] + [str(v) for v in Bj] + [str(v) for v in Cj] +
             ['al', 'b1', 'b2', 'g1', 'g2', 't', 't2'])
    L = [f'ring rr = 0, ({", ".join(vars_)}), dp;']
    extra = [str(v) for v in Bj] if branch == 'S1' else \
            [str(e).replace('**', '^') for e in w2c]
    L.append('ideal LD = ' +
             ','.join(str(e).replace('**', '^') for e in eqs) + ',' +
             ','.join(extra) + ', t*al*b1*g2-t*al*b2*g1-1;')
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
    print(f'(a,e)=({a_},{e_}) branch {branch}: {verdict}', flush=True)
    return verdict

summary = {}
for pair in [(3, 1), (3, 2), (4, 1)]:
    for br in ('S1', 'S2'):
        summary[(pair, br)] = cert(*pair, br)
print('SPLIT SUMMARY:', summary, flush=True)
