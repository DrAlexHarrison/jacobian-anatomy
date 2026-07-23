#!/usr/bin/env python3
# 2026-07-21: GL(2)-stratum master identity check,
# class rep (a,e) = (3,1) (f0723 shape), fully generic coefficients.
# Chart (r, u0) = (z/y, x^e y^a); Phi = (C/B, u0 A^e B^a);
# Jac_{(r,u0)}(Phi) * B^2 = detJF * A^{e-1} * B^a;  det0 = alpha*detM0.
# Run: /usr/bin/python3 < q3_gl2_master_check.py   (expect three PASS lines)
import sympy as sp
a_, e_ = 3, 1
x, y, z, r, u0 = sp.symbols('x y z r u0')
al = sp.Symbol('alpha')
aj = sp.symbols('a0:4')
b1, b2 = sp.symbols('beta1 beta2')
g1, g2 = sp.symbols('gamma1 gamma2')
bj = sp.symbols('b0:5')
cj = sp.symbols('c0:5')
lhat = sum(aj[j]*r**j for j in range(4))
Bhat = sum(bj[j]*r**j for j in range(5))
Chat = sum(cj[j]*r**j for j in range(5))
uband = [x*y**(3-j)*z**j for j in range(4)]
F1 = sp.expand(x*(al + sum(aj[j]*uband[j] for j in range(4))))
F2 = sp.expand(b1*y + b2*z + sum(bj[j]*x*y**(4-j)*z**j for j in range(5)))
F3 = sp.expand(g1*y + g2*z + sum(cj[j]*x*y**(4-j)*z**j for j in range(5)))
detJ = sp.Matrix([F1, F2, F3]).jacobian([x, y, z]).det()
chart = {x: u0/y**3, z: r*y}
D = sp.simplify(sp.expand(detJ.subs(chart)))
assert not D.free_symbols & {y}
A = al + u0*lhat
B = b1 + b2*r + u0*Bhat
C = g1 + g2*r + u0*Chat
rp = sp.simplify((F3/F2).subs(chart))
assert sp.simplify(rp - C/B) == 0
u0p = sp.simplify((F1*F2**3).subs(chart))
assert sp.simplify(u0p - u0*A*B**3) == 0
print('Phi = (C/B, u0*A*B^3) verified in-chart')
Phi1, Phi2 = C/B, u0*A**e_*B**a_
Jac = sp.simplify(sp.diff(Phi1, r)*sp.diff(Phi2, u0) - sp.diff(Phi1, u0)*sp.diff(Phi2, r))
lhs = sp.simplify(Jac*B**2)
rhs = sp.simplify(D*A**(e_-1)*B**a_)
print('Identity Jac(Phi)*B^2 = detJF*A^{e-1}*B^a:',
      'PASS' if sp.simplify(sp.expand(lhs - rhs)) == 0 else 'FAIL')
det0 = detJ.subs({x: 0, y: 0, z: 0})
print('det0 = alpha*detM0:',
      'PASS' if sp.simplify(det0 - al*(b1*g2 - b2*g1)) == 0 else 'FAIL')
