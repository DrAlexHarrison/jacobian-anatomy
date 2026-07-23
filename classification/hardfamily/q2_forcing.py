#!/usr/bin/env python3
# 2026-07-21: Q2 forcing certificates, affine GL(2) classes.
# For (a,e) in {(3,1),(3,2),(4,1)}: ladder = u0-coefficients of (GL2-ZW)
# with fully generic coefficients; per-coefficient Rabinowitsch unit
# certificates on every lhat coefficient. ALL UNIT => A const forced =>
# class EMPTY via the slice-Moh endgame (notebook Entry 24).
import sympy as sp, subprocess

r, u0 = sp.symbols('r u0')

def run(a_, e_):
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
    targets = [str(a) for a in ajs]
    allvars = ([str(v) for v in ajs] + [f'Bb{j}' for j in range(a_+2)] +
               [f'Cc{j}' for j in range(a_+2)] +
               ['al', 'b1', 'b2', 'g1', 'g2', 't', 't2'])
    L = [f'ring rr = 0, ({", ".join(allvars)}), dp;']
    L.append('ideal LD = ' + ','.join(str(e).replace('**','^') for e in eqs) +
             ', t*al*(b1*g2-b2*g1)-1;')
    for tg in targets:
        L.append(f'ideal J{tg} = LD, t2*{tg}-1;')
        L.append(f'if (reduce(1, groebner(J{tg})) == 0) {{ "UNIT {tg}"; }} '
                 f'else {{ "NONUNIT {tg}"; }}')
    L.append('exit;')
    src = '\n'.join(L) + '\n'
    p = subprocess.run(['/usr/bin/Singular','-q'], input=src, capture_output=True,
                       text=True, timeout=3500)
    outs = [l for l in p.stdout.strip().splitlines() if l.strip()]
    nonunit = [l for l in outs if l.startswith('NONUNIT')]
    print(f'(a,e)=({a_},{e_}): {len(eqs)} ladder eqs, {len(targets)} lhat coeffs: '
          + ('ALL UNIT => A const forced => class EMPTY via slice-Moh'
             if not nonunit else 'NOT FORCED: ' + ' '.join(nonunit)), flush=True)
    return not nonunit

res = {}
for pair in [(3,1), (3,2), (4,1)]:
    try:
        res[pair] = run(*pair)
    except subprocess.TimeoutExpired:
        print(f'(a,e)={pair}: TIMEOUT at 3500s (UNDECIDED)', flush=True)
        res[pair] = None
print('SUMMARY:', res, flush=True)
