#!/usr/bin/env python3
# 2026-07-20: GL(2)-block stratum quotient
# identity, w = (1,-1,-1) frame (system 734). All four checks must PASS.
# Run: /usr/bin/python3 < q2_gl2_identity_check.py
import sympy as sp
x, y, z, u, v = sp.symbols('x y z u v')

def gen(name, d):
    cs = []
    e = 0
    for m in range(d+1):
        for n in range(d+1-m):
            c = sp.Symbol(f'{name}{m}{n}')
            cs.append(c); e += c*u**m*v**n
    return e

A  = gen('a', 2)
B1 = gen('b', 2); B2 = gen('c', 2)
C1 = gen('d', 2); C2 = gen('e', 2)
U, V = x*y, x*z
F1 = sp.expand(x*A.subs({u:U, v:V}))
F2 = sp.expand(y*B1.subs({u:U, v:V}) + z*B2.subs({u:U, v:V}))
F3 = sp.expand(y*C1.subs({u:U, v:V}) + z*C2.subs({u:U, v:V}))
detJ = sp.expand(sp.Matrix([F1,F2,F3]).jacobian([x,y,z]).det())
L1 = u*B1 + v*B2
L2 = u*C1 + v*C2
Phi1 = A*L1; Phi2 = A*L2
JacPhi = sp.expand(sp.diff(Phi1,u)*sp.diff(Phi2,v) - sp.diff(Phi1,v)*sp.diff(Phi2,u))
lhs = sp.expand(JacPhi.subs({u:U, v:V}))
rhs = sp.expand(detJ * A.subs({u:U, v:V}))
print('Identity Jac(Phi) = detJF * A :', 'PASS' if sp.expand(lhs-rhs)==0 else 'FAIL')
def jac2(f,g): return sp.diff(f,u)*sp.diff(g,v) - sp.diff(f,v)*sp.diff(g,u)
div = sp.expand(A*jac2(L1,L2) + L2*jac2(L1,A) + L1*jac2(A,L2))
print('Identity divided-bracket = detJF:', 'PASS' if sp.expand(div.subs({u:U,v:V}) - detJ)==0 else 'FAIL')
det0 = detJ.subs({x:0,y:0,z:0})
pred = (A*(B1*C2-B2*C1)).subs({u:0,v:0})
print('det0 structure:', 'PASS' if sp.expand(det0-pred)==0 else 'FAIL')
Fs = [x, sp.expand(y + z*sp.Symbol('q')*(x*z)), z]
dJs = sp.Matrix(Fs).jacobian([x,y,z]).det()
print('shear automorphism det == 1:', 'PASS' if sp.simplify(dJs-1)==0 else 'FAIL')
