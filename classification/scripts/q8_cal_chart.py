#!/usr/bin/env python3
# 2026-07-22: CAL-A for the {734,735} hand kill.
# Verify, on the RAW f0734 family (exact monomial supports from jobs/f0734.sing,
# generic q0..q23, no normalization):
#   (1) det JF is weight-0: det JF = Delta(u,v) with u=xy, v=xz, exactly.
#   (2) chart bijection: A = q0 + u0*(q2+q1*r) + u0^2*(q5+q4*r+q3*r^2),
#       B = (q7+q6*r) + u0*(q10+q9*r+q8*r^2) + u0^2*(q14+q13*r+q12*r^2+q11*r^3),
#       C = (q16+q15*r) + u0*(q19+q18*r+q17*r^2) + u0^2*(q23+q22*r+q21*r^2+q20*r^3),
#       reproduces F exactly: F1 = x*A, F2 = y*B, F3 = y*C at r=z/y, u0=xy.
#   (3) LADDER IDENTITY: Delta(u0, u0*r)*B = D*G - u0*E*(AB)_r as polynomials
#       in (r, u0), where D = B*C_r - C*B_r, E = B*C_u0 - C*B_u0,
#       G = (u0*A*B)_u0.  [This is Jac(Phi)*B^2 = detJF*B for Phi=(C/B, u0*A*B).]
#   (4) det0 = q0*(q7*q15 - q6*q16) = A(0)*detM0.
import sympy as sp

x, y, z, r, u0 = sp.symbols('x y z r u0')
q = sp.symbols('q0:24')

F1 = q[0]*x + q[1]*x**2*z + q[2]*x**2*y + q[3]*x**3*z**2 + q[4]*x**3*y*z + q[5]*x**3*y**2
F2 = (q[6]*z + q[7]*y + q[8]*x*z**2 + q[9]*x*y*z + q[10]*x*y**2
      + q[11]*x**2*z**3 + q[12]*x**2*y*z**2 + q[13]*x**2*y**2*z + q[14]*x**2*y**3)
F3 = (q[15]*z + q[16]*y + q[17]*x*z**2 + q[18]*x*y*z + q[19]*x*y**2
      + q[20]*x**2*z**3 + q[21]*x**2*y*z**2 + q[22]*x**2*y**2*z + q[23]*x**2*y**3)

J = sp.Matrix(3, 3, lambda i, j: sp.diff([F1, F2, F3][i], [x, y, z][j]))
Delta = sp.expand(J.det())

# (1) weight-0 check + extraction of Delta as poly in (u,v)
u, v = sp.symbols('u v')
duv = {}
ok_w = True
for mono, c in sp.Poly(Delta, x, y, z).terms():
    a, b, cc = mono
    if a != b + cc:
        ok_w = False
        break
    duv[(b, cc)] = duv.get((b, cc), 0) + c
print('CHECK-1 weight-0:', 'PASS' if ok_w else 'FAIL')
Duv = sum(c*u**i*v**j for (i, j), c in duv.items())
print('CHECK-1b reconstruct:', 'PASS' if sp.expand(
    Duv.subs({u: x*y, v: x*z}) - Delta) == 0 else 'FAIL')

# (2) chart blocks
A = q[0] + u0*(q[2] + q[1]*r) + u0**2*(q[5] + q[4]*r + q[3]*r**2)
B = (q[7] + q[6]*r) + u0*(q[10] + q[9]*r + q[8]*r**2) + u0**2*(q[14] + q[13]*r + q[12]*r**2 + q[11]*r**3)
C = (q[16] + q[15]*r) + u0*(q[19] + q[18]*r + q[17]*r**2) + u0**2*(q[23] + q[22]*r + q[21]*r**2 + q[20]*r**3)
sub = {r: z/y, u0: x*y}
ok2 = (sp.simplify(x*A.subs(sub) - F1) == 0 and
       sp.simplify(y*B.subs(sub) - F2) == 0 and
       sp.simplify(y*C.subs(sub) - F3) == 0)
print('CHECK-2 chart bijection:', 'PASS' if ok2 else 'FAIL')

# (3) ladder identity
D = sp.expand(B*sp.diff(C, r) - C*sp.diff(B, r))
E = sp.expand(B*sp.diff(C, u0) - C*sp.diff(B, u0))
G = sp.expand(sp.diff(u0*A*B, u0))
LHS = sp.expand(Duv.subs({u: u0, v: u0*r})*B)
RHS = sp.expand(D*G - u0*E*sp.diff(A*B, r))
print('CHECK-3 ladder identity det JF * B = D*G - u0*E*(AB)_r:',
      'PASS' if sp.expand(LHS - RHS) == 0 else 'FAIL')

# (4) det0
det0 = q[0]*(q[7]*q[15] - q[6]*q[16])
print('CHECK-4 det0 = A(0)*detM0:',
      'PASS' if sp.expand(Delta.subs({x: 0, y: 0, z: 0}) - det0) == 0 else 'FAIL')
