#!/usr/bin/env python3
# 2026-07-22: {734,735} normalized-ladder forcing certificates.
# NORMALIZATION LEMMA (hand, Entry 30): source-GL(2) on (y,z) (Mobius on r,
# u0-twist) + target-GL(2) on (F2,F3) + scalings preserve the family, Keller,
# det0 != 0, and injectivity, and bring (B0, C0) -> (1, r), detM0 -> 1
# (possible since det0 = alpha*detM0 != 0 makes B0, C0 independent).
# Normalized ladder: A = al + u0*Ah + u0^2*Ac, B = 1 + u0*Bh + u0^2*Bc,
# C = r + u0*Ch + u0^2*Cc; k = al. Certs: Rabinowitsch unit per Ah/Ac coeff.
import sympy as sp, subprocess

r, u0 = sp.symbols('r u0')
al = sp.Symbol('al')
def poly(name, d):
    return sum(sp.Symbol(f'{name}{j}')*r**j for j in range(d+1))
Ah, Ac = poly('ah', 1), poly('ac', 2)
Bh, Bc = poly('bh', 2), poly('bc', 3)
Ch, Cc = poly('ch', 2), poly('cc', 3)
A = al + u0*Ah + u0**2*Ac
B = 1 + u0*Bh + u0**2*Bc
C = r + u0*Ch + u0**2*Cc
k = al
I = sp.expand((B*sp.diff(C,r) - C*sp.diff(B,r))*(A*B + u0*sp.diff(A,u0)*B + u0*A*sp.diff(B,u0))
              - (B*sp.diff(C,u0) - C*sp.diff(B,u0))*u0*(sp.diff(A,r)*B + A*sp.diff(B,r))
              - k*B)
P = sp.Poly(I, u0, r)
eqs = sorted({sp.expand(c) for c in P.coeffs()} - {sp.S(0)}, key=sp.default_sort_key)
targets = ['ah0','ah1','ac0','ac1','ac2']
vars_ = (['ah0','ah1','ac0','ac1','ac2'] +
         [f'bh{j}' for j in range(3)] + [f'bc{j}' for j in range(4)] +
         [f'ch{j}' for j in range(3)] + [f'cc{j}' for j in range(4)] +
         ['al','t','t2'])
L = [f'ring rr = 0, ({", ".join(vars_)}), dp;']
L.append('ideal LD = ' + ','.join(str(e).replace('**','^') for e in eqs) + ', t*al-1;')
res = []
for tg in targets:
    L.append(f'ideal J{tg} = LD, t2*{tg}-1;')
    L.append(f'if (reduce(1, groebner(J{tg})) == 0) {{ "UNIT {tg}"; }} else {{ "NONUNIT {tg}"; }}')
L.append('exit;')
src = '\n'.join(L) + '\n'
print(f'{len(eqs)} normalized ladder eqs, {len(vars_)} vars', flush=True)
try:
    p = subprocess.run(['/usr/bin/Singular','-q'], input=src, capture_output=True,
                       text=True, timeout=3500)
    outs = [l for l in p.stdout.strip().splitlines() if l.strip()]
    print('\n'.join(outs), flush=True)
    nonunit = [l for l in outs if l.startswith('NONUNIT')]
    print('Q2-a1: ' + ('ALL UNIT => A const forced => {734,735} EMPTY via Reduction 12.1'
          if (outs and not nonunit) else 'NOT ALL UNIT / INCOMPLETE'), flush=True)
except subprocess.TimeoutExpired:
    print('TIMEOUT 3500s (UNDECIDED)', flush=True)
