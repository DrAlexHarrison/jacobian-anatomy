#!/usr/bin/env python3
# 2026-07-20: (ZW) determinant identity +
# u-ladder extraction, fully generic capped coefficients, (c,b) in
# {(2,1),(1,2),(3,1),(1,3)}. Run: /usr/bin/python3 < q1_identity_check.py
# Expected: identity PASS + ladder PASS for all four, then ALL PASS.
# Verbatim the script run in-session; results logged in notebook Entry 10.
import sympy as sp

x, y, z, u = sp.symbols('x y z u')
k = sp.Symbol('k')

def caps_poly(name, degx, degu_terms):
    coeffs = []
    expr = 0
    for m, dm in degu_terms:
        for i in range(dm+1):
            c = sp.Symbol(f'{name}{m}_{i}')
            coeffs.append(c)
            expr += c * x**i * u**m
    return expr, coeffs

ok_all = True
for (c_, b_) in [(2,1),(1,2),(3,1),(1,3)]:
    d = b_ + c_
    fterms = [(m, 6 - d*m) for m in range(0, 6//d + 1)]
    aterms = [(m, 5 - d*m) for m in range(0, 5//d + 1)]
    f, fc = caps_poly('f', 6, fterms)
    A, ac = caps_poly('a', 5, aterms)
    B, bc = caps_poly('b', 5, aterms)
    U = y**c_ * z**b_
    F1 = f.subs(u, U); F2 = sp.expand(y*A.subs(u, U)); F3 = sp.expand(z*B.subs(u, U))
    det = sp.expand(sp.Matrix([F1,F2,F3]).jacobian([x,y,z]).det())
    fx, fu = sp.diff(f,x), sp.diff(f,u)
    Ax, Au = sp.diff(A,x), sp.diff(A,u)
    Bx, Bu = sp.diff(B,x), sp.diff(B,u)
    rhs = fx*(A*B + u*(c_*B*Au + b_*A*Bu)) - u*fu*(c_*Ax*B + b_*A*Bx)
    diff = sp.expand(det - rhs.subs(u, U))
    ident_ok = (diff == 0)
    ok_all &= ident_ok
    print(f'(c,b)=({c_},{b_}) d={d}: identity {"PASS" if ident_ok else "FAIL"} '
          f'(#coeffs f/A/B = {len(fc)}/{len(ac)}/{len(bc)})')

    if d in (3,4):
        a0, b0, phi, f2c = sp.symbols('a0 b0 phi f2c')
        degPS = 5 - d
        P = sum(sp.Symbol(f'p{i}')*x**i for i in range(degPS+1))
        S = sum(sp.Symbol(f's{i}')*x**i for i in range(degPS+1))
        f1 = sum(sp.Symbol(f'g{i}')*x**i for i in range(6-d+1))
        An = a0*(1 + P*u)
        Bn = b0*(1 + S*u)
        fn = sp.Symbol('c0') + phi*x + f1*u + (f2c*u**2 if d == 3 else 0)
        fxn, fun = sp.diff(fn,x), sp.diff(fn,u)
        dJ = sp.expand(fxn*(An*Bn + u*(c_*Bn*sp.diff(An,u) + b_*An*sp.diff(Bn,u)))
                       - u*fun*(c_*sp.diff(An,x)*Bn + b_*An*sp.diff(Bn,x)))
        poly_u = sp.Poly(dJ - k, u)
        eqs = {int(m[0]): sp.expand(co) for m, co in zip(poly_u.monoms(), poly_u.coeffs())}
        I0 = phi*a0*b0 - k
        I1 = sp.expand(phi*a0*b0*((1+c_)*P + (1+b_)*S) + sp.diff(f1,x)*a0*b0)
        I2 = sp.expand(phi*a0*b0*(1+b_+c_)*P*S + sp.diff(f1,x)*a0*b0*((1+c_)*P+(1+b_)*S)
                       - f1*a0*b0*(c_*sp.diff(P,x) + b_*sp.diff(S,x)))
        I3 = sp.expand(sp.diff(f1,x)*a0*b0*(1+b_+c_)*P*S
                       - f1*a0*b0*(c_*sp.diff(P,x)*S + b_*P*sp.diff(S,x))
                       - (2*f2c*a0*b0*(c_*sp.diff(P,x) + b_*sp.diff(S,x)) if d==3 else 0))
        ladder = {0: I0, 1: I1, 2: I2, 3: I3}
        if d == 3:
            ladder[4] = sp.expand(2*f2c*a0*b0*(c_*sp.diff(P,x)*S + b_*P*sp.diff(S,x)))
        lad_ok = True
        for m, e in eqs.items():
            hand = ladder.get(m, sp.Integer(0))
            if sp.expand(e - hand) != 0 and sp.expand(e + hand) != 0:
                lad_ok = False
                print(f'   u^{m} MISMATCH: extracted {e} vs hand {hand}')
        ok_all &= lad_ok
        print(f'   ladder (u^0..u^{max(eqs)}) vs hand (I)-(IV): {"PASS" if lad_ok else "FAIL"}')
print('ALL PASS' if ok_all else 'FAILURES PRESENT')
