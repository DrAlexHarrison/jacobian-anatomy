#!/usr/bin/env python3
# 2026-07-21: ladder-extraction checks backing Thm 26.1
# (M1g/M3g on the Chat==0 stratum; S1-1/S1-2 on the Bhat==0 stratum), fully
# generic, (a,e) in {(3,1),(3,2),(4,1)}. Run: /usr/bin/python3 < this. Expect
# 12x PASS + ALL PASS. Verbatim the script run in-session (Entry 26).
import sympy as sp
r, u0 = sp.symbols('r u0')
ok_all = True
for (a_, e_) in [(3,1),(3,2),(4,1)]:
    al = sp.Symbol('al'); b1,b2,g1,g2 = sp.symbols('b1 b2 g1 g2')
    ajs = sp.symbols(f'a0:{a_+1}')
    Bj = [sp.Symbol(f'Bb{j}') for j in range(a_+2)]
    Cj = [sp.Symbol(f'Cc{j}') for j in range(a_+2)]
    lhat = sum(ajs[j]*r**j for j in range(a_+1))
    Bhat = sum(Bj[j]*r**j for j in range(a_+2))
    Chat = sum(Cj[j]*r**j for j in range(a_+2))
    detM0 = b1*g2 - b2*g1
    B0 = b1 + b2*r; C0 = g1 + g2*r
    k = al*detM0
    def zw(Bh, Ch):
        A = al + u0*lhat
        B = B0 + u0*Bh
        C = C0 + u0*Ch
        return sp.expand((B*sp.diff(C,r) - C*sp.diff(B,r)) *
                         (A*B + u0*(e_*lhat*B + a_*A*Bh))
                         - u0*(B*Ch - C*Bh) *
                         (e_*u0*sp.diff(lhat,r)*B + a_*A*(b2 + u0*sp.diff(Bh,r)))
                         - k*B)
    L = sp.Poly(zw(Bhat, sp.S(0)), u0)
    cofs = {int(m[0]): sp.expand(c) for m,c in zip(L.monoms(), L.coeffs())}
    M1g = sp.expand(detM0*((1+e_)*lhat*B0 + a_*al*Bhat) + al*B0*(g2*Bhat - C0*sp.diff(Bhat,r)) + a_*al*b2*C0*Bhat)
    M3g = sp.expand((g2*Bhat - C0*sp.diff(Bhat,r))*(1+e_+a_)*lhat*Bhat + C0*Bhat*(e_*sp.diff(lhat,r)*Bhat + a_*lhat*sp.diff(Bhat,r)))
    ok1 = sp.expand(cofs.get(1, sp.S(0)) - M1g) == 0
    ok3 = sp.expand(cofs.get(3, sp.S(0)) - M3g) == 0
    Ls = sp.Poly(sp.expand(zw(sp.S(0), Chat)/(B0)), u0)
    cs = {int(m[0]): sp.expand(c) for m,c in zip(Ls.monoms(), Ls.coeffs())}
    S11 = sp.expand((1+e_)*detM0*lhat + al*B0*sp.diff(Chat,r) - (1+a_)*al*b2*Chat)
    S12 = sp.expand(B0*((1+e_)*lhat*sp.diff(Chat,r) - e_*sp.diff(lhat,r)*Chat) - (1+e_+a_)*b2*lhat*Chat)
    ok4 = sp.expand(cs.get(1, sp.S(0)) - S11) == 0
    ok5 = sp.expand(cs.get(2, sp.S(0)) - S12) == 0
    ok_all &= ok1 and ok3 and ok4 and ok5
    print(f'(a,e)=({a_},{e_}): M1g {"PASS" if ok1 else "FAIL"}, M3g {"PASS" if ok3 else "FAIL"}, '
          f'S1-1 {"PASS" if ok4 else "FAIL"}, S1-2 {"PASS" if ok5 else "FAIL"}')
print('ALL PASS' if ok_all else 'FAILURES')
