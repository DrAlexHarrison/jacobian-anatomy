#!/usr/bin/env python3
# 2026-07-22: GAUGE LEMMA machine checks for the {734,735} hand kill.
# The family (raw f0734 supports, 24 coefficients q0..q23):
#   F1 = x*A(u,v), F2 = y*B1 + z*B2, F3 = y*C1 + z*C2  (blocks deg <= 2 in u,v)
# Actions:
#  (I) General source x target: G = h_t o F o g_s with
#      g_s(x,y,z) = (x, p*y + q*z, s*y + t*z),  h_t(X,Y,Z) = (X, a*Y+b*Z, c*Y+d*Z)
#      (x-scalings omitted: not needed for normalization of (M0, detM0) up to
#       the constants we track; A(0) = q0 is untouched by these).
#      CHECK: G stays in the family (support containment), det JG = det JF o g_s
#      times det(h23)*det(g23), M0_G = h23*M0*g23^T-form, injectivity trivial.
#  (II) Conjugation instance (the residual gauge on normalized members):
#      g generic, h = g^{-1}: verify chart-block transformation formulas
#        Ah_G = (p+q*r)*Ah(mob), Ac_G = (p+q*r)^2*Ac(mob),
#        B_G,m = (p+q*r)^(m+1)/Del * (t*B_m - q*C_m)(mob),
#        C_G,m = (p+q*r)^(m+1)/Del * (-s*B_m + p*C_m)(mob),  mob = (s+t*r)/(p+q*r)
#      and that (B0,C0) = (1, r) is preserved.
#  (III) The two instances used in the tree: shear (p=t=1, q=0) and
#      swap (g = [[0,1],[1,0]]), exact block maps.
import sympy as sp

x, y, z, r, u0 = sp.symbols('x y z r u0')
q = sp.symbols('q0:24')
p_, q_, s_, t_ = sp.symbols('p q s t')
a_, b_, c_, d_ = sp.symbols('a b c d')

F1 = q[0]*x + q[1]*x**2*z + q[2]*x**2*y + q[3]*x**3*z**2 + q[4]*x**3*y*z + q[5]*x**3*y**2
F2 = (q[6]*z + q[7]*y + q[8]*x*z**2 + q[9]*x*y*z + q[10]*x*y**2
      + q[11]*x**2*z**3 + q[12]*x**2*y*z**2 + q[13]*x**2*y**2*z + q[14]*x**2*y**3)
F3 = (q[15]*z + q[16]*y + q[17]*x*z**2 + q[18]*x*y*z + q[19]*x*y**2
      + q[20]*x**2*z**3 + q[21]*x**2*y*z**2 + q[22]*x**2*y**2*z + q[23]*x**2*y**3)
F = [F1, F2, F3]
SUPP = [set(sp.Poly(Fi, x, y, z).monoms()) for Fi in F]

def act(F, gs, ht):
    Fg = [sp.expand(Fi.subs({y: gs[0], z: gs[1]}, simultaneous=True)) for Fi in F]
    return [Fg[0], sp.expand(ht[0](Fg[1], Fg[2])), sp.expand(ht[1](Fg[1], Fg[2]))]

# ---- (I) general action
gs = (p_*y + q_*z, s_*y + t_*z)
ht = (lambda Y, Z: a_*Y + b_*Z, lambda Y, Z: c_*Y + d_*Z)
G = act(F, gs, ht)
ok = all(set(sp.Poly(G[i], x, y, z).monoms()) <= SUPP[i] for i in range(3))
print('CHECK I-support (G in family):', 'PASS' if ok else 'FAIL')
J = lambda H: sp.Matrix(3, 3, lambda i, j: sp.diff(H[i], [x, y, z][j]))
dG = sp.expand(J(G).det())
dF = sp.expand(J(F).det())
dFg = sp.expand(dF.subs({y: gs[0], z: gs[1]}, simultaneous=True))
lhs = sp.expand(dG - (a_*d_ - b_*c_)*(p_*t_ - q_*s_)*dFg)
print('CHECK I-detJ transport:', 'PASS' if lhs == 0 else 'FAIL')
M0G = sp.Matrix([[sp.diff(G[1], y), sp.diff(G[1], z)],
                 [sp.diff(G[2], y), sp.diff(G[2], z)]]).subs({x: 0, y: 0, z: 0})
M0 = sp.Matrix([[q[7], q[6]], [q[16], q[15]]])
h23 = sp.Matrix([[a_, b_], [c_, d_]])
g23 = sp.Matrix([[p_, q_], [s_, t_]])
print('CHECK I-M0 transport (M0_G = h23*M0*g23):',
      'PASS' if sp.simplify(M0G - h23*M0*g23) == sp.zeros(2, 2) else 'FAIL')
print('CHECK I-A(0) untouched (al_G = al):',
      'PASS' if sp.expand(G[0].coeff(x, 1).subs({y: 0, z: 0}) - q[0]) == 0 else 'FAIL')

# ---- (II) conjugation on chart blocks
Del = p_*t_ - q_*s_
hinv = (lambda Y, Z: (t_*Y - q_*Z)/Del, lambda Y, Z: (-s_*Y + p_*Z)/Del)
G = act(F, gs, hinv)
G = [sp.expand(sp.cancel(Gi)) for Gi in G]
ok = all(set(sp.Poly(Gi, x, y, z).monoms()) <= SUPP[i] for i, Gi in enumerate(G))
print('CHECK II-support:', 'PASS' if ok else 'FAIL')

def chart_blocks(H):
    """extract (A_m), (B_m), (C_m) chart blocks of a family member"""
    Aq = sp.Poly(H[0], x, y, z)
    out = {}
    A = sp.expand(H[0]/x)
    # A in u,v: substitute u->u0, v->u0*r via x*A form: A(x,y,z) with u=xy,v=xz
    Bl = sp.expand(H[1])
    Cl = sp.expand(H[2])
    # chart: B = F2/y at z = r*y, u0 = x*y
    sub = {z: r*y}
    Bch = sp.expand(sp.cancel(Bl.subs(sub)/y))
    Cch = sp.expand(sp.cancel(Cl.subs(sub)/y))
    Ach = sp.expand(sp.cancel(H[0].subs(sub)/x))
    # now in terms of x*y = u0: replace x*y -> u0 by poly in (x*y, r)
    res = []
    for expr in (Ach, Bch, Cch):
        pe = sp.Poly(expr, x, y)
        acc = 0
        for (i, j), cc in zip(pe.monoms(), pe.coeffs()):
            assert i == j, (i, j, expr)
            acc += cc*u0**i
        res.append(sp.expand(acc))
    return res  # each a poly in (u0, r)

AF, BF, CF = chart_blocks(F)
AG, BG, CG = chart_blocks(G)
mob = (s_ + t_*r)/(p_ + q_*r)
w = p_ + q_*r
okII = True
for m in range(0, 3):
    Am = AF.coeff(u0, m)
    Bm, Cm = BF.coeff(u0, m), CF.coeff(u0, m)
    predA = sp.expand(sp.cancel(w**m * Am.subs(r, mob)))
    predB = sp.expand(sp.cancel(w**(m+1)/Del * (t_*Bm - q_*Cm).subs(r, mob)))
    predC = sp.expand(sp.cancel(w**(m+1)/Del * (-s_*Bm + p_*Cm).subs(r, mob)))
    for got, pred, nm in ((AG.coeff(u0, m), predA, f'A{m}'),
                          (BG.coeff(u0, m), predB, f'B{m}'),
                          (CG.coeff(u0, m), predC, f'C{m}')):
        if sp.simplify(got - pred) != 0:
            okII = False
            print(f'  MISMATCH {nm}')
print('CHECK II-block formulas (all levels A,B,C):', 'PASS' if okII else 'FAIL')

# normalized preservation: set B0 = 1, C0 = r i.e. q7=1,q6=0,q16=0,q15=1
norm = {q[7]: 1, q[6]: 0, q[16]: 0, q[15]: 1}
B0G = BG.coeff(u0, 0).subs(norm)
C0G = CG.coeff(u0, 0).subs(norm)
print('CHECK II-normalization preserved (B0_G,C0_G)=(1,r):',
      'PASS' if sp.simplify(B0G - 1) == 0 and sp.simplify(C0G - r) == 0 else 'FAIL')

# ---- (III) exact instances
inst = {'shear': {p_: 1, q_: 0, t_: 1}, 'swap': {p_: 0, q_: 1, s_: 1, t_: 0}}
sh = inst['shear']
BGs = sp.expand(BG.subs(sh))
CGs = sp.expand(CG.subs(sh))
okS = (sp.expand(BGs - BF.subs(r, r + s_)) == 0 and
       sp.expand(CGs - (CF.subs(r, r + s_) - s_*BF.subs(r, r + s_))) == 0 and
       sp.expand(AG.subs(sh) - AF.subs(r, r + s_)) == 0)
print('CHECK III-shear (B->B(r+s), C->(C-sB)(r+s), A->A(r+s)):',
      'PASS' if okS else 'FAIL')
sw = inst['swap']
BGw = sp.expand(sp.cancel(BG.subs(sw)))
CGw = sp.expand(sp.cancel(CG.subs(sw)))
AGw = sp.expand(sp.cancel(AG.subs(sw)))
# swap: Del=-1, w=r, mob=1/r: B_G,m = r^(m+1)/(-1)*(0*B - 1*C)(1/r) = r^(m+1)*C_m(1/r)
okW = True
for m in range(0, 3):
    pB = sp.expand(sp.cancel(r**(m+1)*CF.coeff(u0, m).subs(r, 1/r)))
    pC = sp.expand(sp.cancel(r**(m+1)*BF.coeff(u0, m).subs(r, 1/r)))
    pA = sp.expand(sp.cancel(r**m*AF.coeff(u0, m).subs(r, 1/r)))
    if (sp.simplify(BGw.coeff(u0, m) - pB) != 0 or
            sp.simplify(CGw.coeff(u0, m) - pC) != 0 or
            sp.simplify(AGw.coeff(u0, m) - pA) != 0):
        okW = False
print('CHECK III-swap (B_m <-> reversed C_m, A_m reversed):',
      'PASS' if okW else 'FAIL')
