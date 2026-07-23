#!/usr/bin/env python3
# 2026-07-22: machine verification of the intermediate HAND steps
# quoted in the notebook entries.
import sympy as sp

r, u0 = sp.symbols('r u0')
al = sp.Symbol('al')
N = sp.sympify(open('q8_layers.srepr').read())
S = sp.Symbol

def dr(f): return sp.diff(f, r)
W = lambda f, g: sp.expand(f*dr(g) - dr(f)*g)

def specialize(expr, hh, cc_, KK, QQ, AA):
    sub = {}
    for j in range(3):
        sub[S(f'bh{j}')] = sp.expand(hh).coeff(r, j)
        sub[S(f'ch{j}')] = sp.expand(cc_).coeff(r, j)
        sub[S(f'ac{j}')] = sp.expand(AA).coeff(r, j)
    for j in range(4):
        sub[S(f'cc{j}')] = sp.expand(KK).coeff(r, j)
        sub[S(f'bc{j}')] = sp.expand(QQ).coeff(r, j)
    return sp.expand(expr.subs(sub))

# ---- A. a2b-i: Ac=Bc=K=0, c = tau*h; N2* forces h0 = h1 = 0 (shifted), then P = 0
tau, h0, h1, h2 = sp.symbols('tau h0 h1 h2')
u = r - tau
h = h0 + h1*u + h2*u**2
c = tau*h
N2s = specialize(N[2], h, c, 0, 0, 0)*2/al
sols = sp.solve([sp.expand(N2s).coeff(r, j) for j in range(4)], [h0, h1], dict=True)
print('A. a2b-i N2* forces h0=h1=0:',
      'PASS' if all(s.get(h0, 0) == 0 and s.get(h1, 0) == 0 for s in sols) and sols else sols)
P = sp.expand(2*h + dr(c) - r*dr(h))
Pforced = sp.expand(P.subs({h0: 0, h1: 0}))
print('A2. then P == 0 identically:', 'PASS' if Pforced == 0 else f'FAIL: {Pforced}')

# ---- B. a2b-ii coincidence at xi=0: N4 == (al/2)*(6*p1*h2^2*c0*r^4 - 4*p1*h2*k3*r^5)
p1, k3, c0 = sp.symbols('p1 k3 c0')
h = h2*r**2
K = k3*r**3
c = c0 + (p1/2)*r**2
N4s = specialize(N[4], h, c, K, 0, 0)
pred = sp.expand(al/2*(6*p1*h2**2*c0*r**4 - 4*p1*h2*k3*r**5))
print('B. a2b-ii N4 form:', 'PASS' if sp.expand(N4s - pred) == 0 else 'FAIL')
# sanity: N5 vanishes identically under the shapes (root-matched)
N5s = specialize(N[5], h, c, K, 0, 0)
print('B2. a2b-ii N5 == 0 under shapes:', 'PASS' if N5s == 0 else f'note: N5 = {sp.factor(N5s)}')

# ---- C. beta-K0: Bc=Cc=0, c=tau*h; N4 (Eq A) forces h = h2*u^2, Ac = a2*u^2;
#         then N2 (Eq C) forces a2 = 0
a0, a1, a2s = sp.symbols('a0 a1 a2s')
h = h0 + h1*u + h2*u**2
Ac = a0 + a1*u + a2s*u**2
c = tau*h
N4s = specialize(N[4], h, c, 0, 0, Ac)
eqs = [sp.expand(N4s).coeff(r, j) for j in range(6)]
sols = sp.solve(eqs, [h0, h1, a0, a1], dict=True)
ok = sols and all(all(s.get(v, 0) == 0 for v in (h0, h1, a0, a1)) or
                  (s.get(h2) == 0 if h2 in s else False) for s in sols)
print('C. beta-K0 Eq A solutions:', sols if not ok else 'PASS (h0=h1=a0=a1=0 or degenerate)')
N2s = specialize(N[2], h2*u**2, tau*h2*u**2, 0, 0, a2s*u**2)*2
N2red = sp.expand(N2s)
print('C2. beta-K0 Eq C reduces to 6*a2*u^2 == 0:',
      'PASS' if sp.expand(N2red - 6*a2s*u**2) == 0 else f'FAIL: {sp.factor(N2red)}')

# ---- D. beta-coinc at xi=0: the N5/N4 coefficient relations (*), (**)
c1, c2 = sp.symbols('c1 c2')
h = h2*r**2
Ac = a2s*r**2
K = k3*r**3
c = c0 + c1*r + c2*r**2
N5s = specialize(N[5], h, c, K, 0, Ac)
f5 = sp.factor(N5s)
print('D. beta-coinc N5 =', str(f5)[:300])
# claimed: N5 = h2*r^5*(-8*a2*h2*c0 + (5*a2*k3 - 4*a2*h2*c1 - (al/2)*h2*k3*c1)*r) (up to overall const)
pred5 = sp.expand(h2*r**5*(-8*a2s*h2*c0 + (5*a2s*k3 - 4*a2s*h2*c1 - sp.Rational(1, 2)*al*h2*k3*c1)*r))
print('D1. N5 matches claimed form:', 'PASS' if sp.expand(N5s - pred5) == 0 else 'FAIL')
N4s = specialize(N[4], h, sp.expand(c.subs(c0, 0)), K, 0, Ac)
r4 = sp.expand(N4s).coeff(r, 4)
r5 = sp.expand(N4s).coeff(r, 5)
r6 = sp.expand(N4s).coeff(r, 6)
pred_r5 = sp.expand(h2*c2*(4*a2s + 3*al*h2*c1 - 4*al*k3))
pred_r4 = sp.expand(-3*a2s*h2*c1 + 5*a2s*k3 + sp.Rational(3, 2)*al*h2**2*c1**2 - sp.Rational(7, 2)*al*h2*k3*c1)
print('D2. N4 r^6 == 0:', 'PASS' if r6 == 0 else f'FAIL: {r6}')
print('D3. N4 r^5 matches (**r5):', 'PASS' if sp.expand(r5 - pred_r5) == 0 else f'FAIL: {sp.factor(r5)}')
print('D4. N4 r^4 matches (**):', 'PASS' if sp.expand(r4 - pred_r4) == 0 else f'FAIL: {sp.factor(r4)}')

# ---- E. beta-consts: h=b, K=k, Ac=a consts; N5 => c2=c1=0; then N4 = 4ab^2
b_, k_, a_ = sp.symbols('b_ k_ a_')
c = c0 + c1*r + c2*r**2
N5s = specialize(N[5], b_, c, k_, 0, a_)
pred5 = sp.expand(b_**2*(4*a_*(c1 + 2*c2*r) + sp.Rational(1, 2)*al*k_*2*c2))
print('E. beta-consts N5 =', sp.factor(N5s))
N4s = specialize(N[4], b_, c0, k_, 0, a_)
print('E2. beta-consts N4 (c=c0) == 4*a*b^2:',
      'PASS' if sp.expand(N4s - 4*a_*b_**2) == 0 else f'FAIL: {sp.factor(N4s)}')

# ---- F. beta-h0 (K const nonzero): N3 => 3*a*c' + al*k*c'' = 0 form; N2 => 6a = 0
ch0v, ch1v, ch2v = sp.symbols('ch0v ch1v ch2v')
c = ch0v + ch1v*r + ch2v*r**2
N3s = specialize(N[3], 0, c, k_, 0, a_)*2
print('F. beta-h0 N3x2 =', sp.factor(N3s))
N2s = specialize(N[2], 0, ch0v, k_, 0, a_)*2
print('F2. beta-h0 N2x2 (c const) == 6a:',
      'PASS' if sp.expand(N2s - 6*a_) == 0 else f'FAIL: {sp.factor(N2s)}')

# ---- G. a2a: Ac=Bc=Bh=0: N3* == 2*(c''K - c'K') and N2* == 2K' - 2c'^2 + c*c''
K = sum(S(f'k{j}')*r**j for j in range(4))
N3s = specialize(N[3], 0, c, K, 0, 0)*2/al
pred = sp.expand(2*(dr(dr(c))*K - dr(c)*dr(K)))
print('G. a2a N3* == 2(c\'\'K - c\'K\'):', 'PASS' if sp.expand(N3s - pred) == 0 else 'FAIL')
N2s = specialize(N[2], 0, c, K, 0, 0)*2/al
pred = sp.expand(2*dr(K) - 2*dr(c)**2 + c*dr(dr(c)))
print('G2. a2a N2* == 2K\' - 2c\'^2 + c*c\'\':', 'PASS' if sp.expand(N2s - pred) == 0 else 'FAIL')
