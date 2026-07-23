#!/usr/bin/env python3
# Probe: (3,1) ladder + k!=0 + a3!=0 unit check via msolve (F4, char 0).
# GOTCHA: msolve silently misparses parens; emit fully expanded.
import sympy as sp, subprocess

r, u0 = sp.symbols('r u0')
a_, e_ = 3, 1
al = sp.Symbol('al')
b1, b2, g1, g2 = sp.symbols('b1 b2 g1 g2')
ajs = sp.symbols('a0:4')
Bj = [sp.Symbol(f'Bb{j}') for j in range(5)]
Cj = [sp.Symbol(f'Cc{j}') for j in range(5)]
lhat = sum(ajs[j]*r**j for j in range(4))
Bhat = sum(Bj[j]*r**j for j in range(5))
Chat = sum(Cj[j]*r**j for j in range(5))
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
eqs = sorted({sp.expand(c) for c in poly.coeffs()} - {sp.S(0)}, key=sp.default_sort_key)
t, t2 = sp.symbols('t t2')
eqs.append(sp.expand(t*al*(b1*g2 - b2*g1) - 1))
eqs.append(sp.expand(t2*sp.Symbol('a3') - 1))
vars_ = [str(v) for v in ajs] + [f'Bb{j}' for j in range(5)] + \
        [f'Cc{j}' for j in range(5)] + ['al','b1','b2','g1','g2','t','t2']
with open('q2_ms_31_a3.ms', 'w') as f:
    f.write(', '.join(vars_) + '\n0\n')
    f.write(',\n'.join(str(e).replace('**','^') for e in eqs) + '\n')
print(f'wrote q2_ms_31_a3.ms ({len(eqs)} polys, {len(vars_)} vars)')
