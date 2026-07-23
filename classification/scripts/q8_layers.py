#!/usr/bin/env python3
# 2026-07-22: normalized {734,735} ladder — layer extraction + the
# structural convolution table, machine-verified.
# Normalized family (Entry 30 lemma): A = al + u0*Ah + u0^2*Ac,
# B = 1 + u0*Bh + u0^2*Bc, C = r + u0*Ch + u0^2*Cc; k = al.
# Keller <=> N_m = 0 for m = 1..8 where N_m = u0^m-layer of
#   D*G - u0*E*(AB)_r - al*B,  D = B*C_r - C*B_r, E = B*C_u0 - C*B_u0,
#   G = (u0*A*B)_u0.
# Convolution table (hand-derived, verified below):
#   D0 = 1, D1 = Bh + Ch' - r*Bh', D2 = Bc - r*Bc' + Cc' + W(Bh,Ch),
#   D3 = Bh*Cc' + Bc*Ch' - Ch*Bc' - Cc*Bh', D4 = W(Bc,Cc)
#   E0 = Ch - r*Bh, E1 = 2*(Cc - r*Bc), E2 = Bh*Cc - Bc*Ch, E3 = 0
#   (AB)_0 = al, (AB)_1 = al*Bh + Ah, (AB)_2 = al*Bc + Ah*Bh + Ac,
#   (AB)_3 = Ah*Bc + Ac*Bh, (AB)_4 = Ac*Bc;  G_m = (m+1)*(AB)_m
#   N_m = sum_{i+j=m} D_i*G_j - sum_{i+j=m-1} E_i*(AB)_j' - al*B_m
# N1 = 2*Ah + al*P, P := 2*Bh + Ch' - r*Bh'  (deg <= 1: r^2 cancels).
# Then Ah := -al*P/2 substituted; N2..N8 saved to q8_layers.srepr.
import sympy as sp

r, u0 = sp.symbols('r u0')
al = sp.Symbol('al')
def poly(name, d):
    return sum(sp.Symbol(f'{name}{j}')*r**j for j in range(d+1))
Ah_free = poly('ah', 1)
Ac = poly('ac', 2)
Bh, Bc = poly('bh', 2), poly('bc', 3)
Ch, Cc = poly('ch', 2), poly('cc', 3)

def dr(f): return sp.diff(f, r)
W = lambda f, g: sp.expand(f*dr(g) - dr(f)*g)

def brute(Ah):
    A = al + u0*Ah + u0**2*Ac
    B = 1 + u0*Bh + u0**2*Bc
    C = r + u0*Ch + u0**2*Cc
    I = sp.expand((B*dr(C) - C*dr(B))*sp.diff(u0*A*B, u0)
                  - u0*(B*sp.diff(C, u0) - C*sp.diff(B, u0))*dr(A*B) - al*B)
    Po = sp.Poly(I, u0)
    return {int(m[0]): sp.expand(c) for m, c in zip(Po.monoms(), Po.coeffs())}

def conv_layers(Ah):
    D = [sp.S(1), sp.expand(Bh + dr(Ch) - r*dr(Bh)),
         sp.expand(Bc - r*dr(Bc) + dr(Cc) + W(Bh, Ch)),
         sp.expand(Bh*dr(Cc) + Bc*dr(Ch) - Ch*dr(Bc) - Cc*dr(Bh)),
         W(Bc, Cc)]
    E = [sp.expand(Ch - r*Bh), sp.expand(2*(Cc - r*Bc)),
         sp.expand(Bh*Cc - Bc*Ch), sp.S(0)]
    AB = [al, sp.expand(al*Bh + Ah), sp.expand(al*Bc + Ah*Bh + Ac),
          sp.expand(Ah*Bc + Ac*Bh), sp.expand(Ac*Bc)]
    Blev = {0: sp.S(1), 1: Bh, 2: Bc}
    N = {}
    for m in range(0, 9):
        t = sum(D[i]*(j+1)*AB[j] for i in range(5) for j in range(5) if i + j == m)
        t -= sum(E[i]*dr(AB[j]) for i in range(4) for j in range(5) if i + j == m - 1)
        t -= al*Blev.get(m, 0)
        N[m] = sp.expand(t)
    return N

# --- 1. convolution table == brute expansion (Ah free)
b = brute(Ah_free)
c = conv_layers(Ah_free)
ok = all(sp.expand(b.get(m, 0) - c[m]) == 0 for m in range(0, 9))
print('CHECK conv-table == brute (all layers, Ah free):', 'PASS' if ok else 'FAIL')

# --- 2. N0, N1
P = sp.expand(2*Bh + dr(Ch) - r*dr(Bh))
print('CHECK N0 == 0:', 'PASS' if c[0] == 0 else 'FAIL')
print('CHECK deg P <= 1:', 'PASS' if sp.degree(P, r) <= 1 else 'FAIL')
print('CHECK N1 == 2*Ah + al*P:',
      'PASS' if sp.expand(c[1] - (2*Ah_free + al*P)) == 0 else 'FAIL')

# --- 3. substitute Ah = -al*P/2, save N2..N8
Ah_sub = sp.expand(-al*P/2)
N = conv_layers(Ah_sub)
print('CHECK N1 == 0 after substitution:', 'PASS' if sp.expand(N[1]) == 0 else 'FAIL')
for m in range(2, 9):
    dm = sp.degree(N[m], r) if N[m] != 0 else -1
    print(f'N{m}: r-deg {dm}, {len(sp.Poly(N[m], r).all_coeffs()) if N[m] != 0 else 0} coeff eqs')

with open('q8_layers.srepr', 'w') as f:
    f.write(sp.srepr({m: N[m] for m in range(2, 9)}))
print('saved q8_layers.srepr')

# --- 4. structural spot checks used by the hand tree
Wbc = W(Bc, Cc)
print('CHECK N8 == 5*Ac*Bc*W(Bc,Cc):',
      'PASS' if sp.expand(N[8] - 5*Ac*Bc*Wbc) == 0 else 'FAIL')
z_ac = {sp.Symbol(f'ac{j}'): 0 for j in range(3)}
N7a = sp.expand(N[7].subs(z_ac))
print('CHECK N7|_{Ac=0} == 4*Ah*Bc*W(Bc,Cc):',
      'PASS' if sp.expand(N7a - 4*Ah_sub*Bc*Wbc) == 0 else 'FAIL')
z_bc = {sp.Symbol(f'bc{j}'): 0 for j in range(4)}
print('CHECK N7|_{Bc=0} == 0:',
      'PASS' if sp.expand(N[7].subs(z_bc)) == 0 else 'FAIL')
N6b = sp.expand(N[6].subs(z_bc))
tgt = sp.expand(4*Ac*Bh*W(Bh, Cc) - Bh*Cc*dr(Ac*Bh))
print('CHECK N6|_{Bc=0} == 4*Ac*Bh*W(Bh,Cc) - Bh*Cc*(Ac*Bh)\':',
      'PASS' if sp.expand(N6b - tgt) == 0 else 'FAIL')
z_cc = {sp.Symbol(f'cc{j}'): 0 for j in range(4)}
N7g = sp.expand(N[7].subs(z_cc))
tgt7g = sp.expand(Bc*(5*Ac*Bc*dr(Ch) - 4*Ac*dr(Bc)*Ch + dr(Ac)*Bc*Ch))
print('CHECK N7|_{Cc=0} == Bc*(5*Ac*Bc*Ch\' - 4*Ac*Bc\'*Ch + Ac\'*Bc*Ch):',
      'PASS' if sp.expand(N7g - tgt7g) == 0 else 'FAIL')
N5ab = sp.expand(N[5].subs(z_ac).subs(z_bc))
tgt5 = sp.expand(Bh*(3*Ah_sub*W(Bh, Cc) - Cc*dr(Ah_sub*Bh)))
print('CHECK N5|_{Ac=Bc=0} == Bh*(3*Ah*W(Bh,Cc) - Cc*(Ah*Bh)\'):',
      'PASS' if sp.expand(N5ab - tgt5) == 0 else 'FAIL')
print('CHECK N6|_{Ac=Bc=0} == 0:',
      'PASS' if sp.expand(N[6].subs(z_ac).subs(z_bc)) == 0 else 'FAIL')
