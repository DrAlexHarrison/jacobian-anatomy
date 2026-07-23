#!/usr/bin/env python3
# 2026-07-21: a = e = 1 class ({734,735}) ladder extraction.
# Chart (r, u0) = (z/y, xy); A = al + u0*Ah(r) + u0^2*Ac(r) (deg 1, 2),
# B = B0 + u0*Bh + u0^2*Bc (deg 1, 2, 3), C likewise. Identity:
# (B*C_r - C*B_r)*(AB + u0*A_u0*B + u0*A*B_u0)
#   - (B*C_u0 - C*B_u0)*u0*(A_r*B + A*B_r) = k*B,  k = al*detM0.
# Emit each u0-layer; FACTOR the top layers to reveal branch structure.
import sympy as sp
r, u0 = sp.symbols('r u0')
al = sp.Symbol('al')
b1, b2, g1, g2 = sp.symbols('b1 b2 g1 g2')
def poly(name, d):
    return sum(sp.Symbol(f'{name}{j}')*r**j for j in range(d+1))
Ah, Ac = poly('ah', 1), poly('ac', 2)
Bh, Bc = poly('bh', 2), poly('bc', 3)
Ch, Cc = poly('ch', 2), poly('cc', 3)
A = al + u0*Ah + u0**2*Ac
B = (b1 + b2*r) + u0*Bh + u0**2*Bc
C = (g1 + g2*r) + u0*Ch + u0**2*Cc
k = al*(b1*g2 - b2*g1)
I = sp.expand((B*sp.diff(C,r) - C*sp.diff(B,r))*(A*B + u0*sp.diff(A,u0)*B + u0*A*sp.diff(B,u0))
              - (B*sp.diff(C,u0) - C*sp.diff(B,u0))*u0*(sp.diff(A,r)*B + A*sp.diff(B,r))
              - k*B)
P = sp.Poly(I, u0)
layers = {int(m[0]): sp.expand(c) for m, c in zip(P.monoms(), P.coeffs())}
print('u0-layers present:', sorted(layers))
for m in sorted(layers, reverse=True)[:3]:
    f = sp.factor(layers[m])
    print(f'--- u0^{m} FACTORED:')
    print(sp.sstr(f)[:600])
