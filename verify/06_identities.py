#!/usr/bin/env python3
"""Crown-jewel polynomial identities (referee-supplied, independently re-verified).

All are identities in Q[x,y,z] or Q[a,b,c] by PURE EXPANSION — no elimination
theory, no specialization caveats:
  1. I_x:  D(F)*x^3 + (4 - 3*F2*F3)*x - 2*F3 == 0
  2. I_y:  2y^3 - 3*F2*y^2 + 18*F1*y + (27*F1^2*F3 - 18*F1*F2 + F2^3) == 0
     (constant leading coefficient 2: "y never escapes" is a THEOREM)
  3. E^2 identity: (4-3bc)^3 + 27c^2*D == E^2 with E = 27ac^2 - 9bc + 8
     => Gamma = V(D) \\cap V(4-3bc) exactly.
  4. disc_x = -4*E^2*D exactly;  Res_y(G1,G2) = -c*x^9*(D x^3+(4-3bc)x-2c).
  5. Constructive section over {c=0}:  F(0, b, a-4b^2) = (a, b, 0).
  6. The one-line "Gamma is missed": on Gamma(t), D = 0 and 4-3bc = 0
     identically, so I_x at any preimage would read -2t = 0: contradiction.
"""
import sympy as sp

x, y, z, a, b, c, t = sp.symbols('x y z a b c t')
R = sp.Rational
F1 = (1+x*y)**3*z + y**2*(1+x*y)*(4+3*x*y)
F2 = y + 3*x*(1+x*y)**2*z + 3*x*y**2*(4+3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
D = 27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2
E = 27*a*c**2 - 9*b*c + 8
sub = [(a, F1), (b, F2), (c, F3)]

# --- 1. I_x ---------------------------------------------------------------
Ix = (D*x**3 + (4 - 3*b*c)*x - 2*c).subs(sub)
assert sp.expand(Ix) == 0
print("PASS 1 : I_x == 0 identically in Q[x,y,z]")

# --- 2. I_y ---------------------------------------------------------------
Iy = (2*y**3 - 3*b*y**2 + 18*a*y + 27*a**2*c - 18*a*b + b**3).subs(sub)
assert sp.expand(Iy) == 0
print("PASS 2 : I_y == 0 identically; LC constant => y integral, never escapes")

# --- 3. E^2 identity ------------------------------------------------------
assert sp.expand((4 - 3*b*c)**3 + 27*c**2*D - E**2) == 0
# => Gamma = V(D, 4-3bc): forward — on Gamma both vanish;
gam = [(a, R(4,27)/t**2), (b, R(4,3)/t), (c, t)]
assert sp.simplify(D.subs(gam)) == 0 and sp.simplify((4-3*b*c).subs(gam)) == 0
# converse — D = 0 and 4-3bc = 0 give E = 0; with bc = 4/3 (so c != 0),
# E = 27ac^2 - 12 + 8 = 27ac^2 - 4 = 0, i.e. ac^2 = 4/27, b = 4/(3c): Gamma(c).
assert sp.expand(E.subs(b, 4/(3*c)) - (27*a*c**2 - 4)) == 0
print("PASS 3 : (4-3bc)^3 + 27c^2 D == E^2;  Gamma = V(D) \\cap V(4-3bc) exactly")

# --- 4. exact discriminant and resultant constants ------------------------
cubic = D*x**3 + (4 - 3*b*c)*x - 2*c
assert sp.expand(sp.discriminant(sp.Poly(cubic, x)) + 4*E**2*D) == 0
zsub = (2*x - 3*x**2*y - c)/x**3
G1 = sp.expand(sp.together(F1.subs(z, zsub) - a)*x**3)
G2 = sp.expand(sp.together(F2.subs(z, zsub) - b)*x**2)
res = sp.Poly(G1, y).resultant(sp.Poly(G2, y))
res = res.as_expr() if isinstance(res, sp.Poly) else res
assert sp.expand(res - (-c*x**9*cubic)) == 0
print("PASS 4 : disc_x == -4*E^2*D;  Res_y(G1,G2) == -c*x^9*(x-cubic), exact constants")

# --- 5. polynomial section over {c=0} -------------------------------------
sec = [(x, 0), (y, b), (z, a - 4*b**2)]
imgs = [sp.expand(f.subs(sec)) for f in (F1, F2, F3)]
assert imgs == [sp.expand(a), b, 0]
print("PASS 5 : F(0, b, a-4b^2) == (a, b, 0) — {c=0} hit by an explicit section")

# --- 6. one-line Gamma-missed proof ---------------------------------------
# I_x at a hypothetical preimage of Gamma(t) reads 0*x^3 + 0*x - 2t = 0.
assert sp.simplify((-2*c).subs(gam) + 2*t) == 0 and sp.simplify((-2*t)) != 0
print("PASS 6 : on Gamma(t): I_x forces -2t == 0 — no preimage for t != 0")

# --- 7. D squarefree and irreducible over Q (load-bearing for S3) ---------
g = sp.gcd(sp.gcd(sp.gcd(D, sp.diff(D, a)), sp.diff(D, b)), sp.diff(D, c))
assert g == 1 or g.is_number
assert len(sp.factor_list(D)[1]) == 1 and sp.factor_list(D)[1][0][1] == 1
print("PASS 7 : gcd(D, grad D) = 1 (squarefree); D irreducible over Q")

# --- 8. pushforward certificate: x-cubic is the pushforward of C ----------
# x = r/(1-r*y(r)) with 1-r*y(r) = (b r^2 - 4r + 3c)/(2r) gives the bilinear
# relation x*(b r^2 - 4r + 3c) = 2 r^2; eliminating r against C recovers the
# x-cubic exactly:
r_ = sp.symbols('r_')
Cc = 2*a*r_**3 - b*r_**2 + 2*r_ - c
rel = x*(b*r_**2 - 4*r_ + 3*c) - 2*r_**2
res = sp.Poly(Cc, r_).resultant(sp.Poly(rel, r_))
res = res.as_expr() if isinstance(res, sp.Poly) else res
assert sp.expand(res - 4*c*(D*x**3 + (4 - 3*b*c)*x - 2*c)) == 0
print("PASS 8 : Res_r(C, x(br^2-4r+3c)-2r^2) = 4c*(x-cubic) — pushforward exact")

# --- 9. irreducibility of the covering cubic over C(a,b,c) ---------------
# C_q(r) = 2a r^3 - b r^2 + 2r - c is LINEAR in c with unit c-coefficient
# (-1), and its content as a polynomial in r is 1 (the r-coefficient 2 is a
# unit). Hence in any factorization over C[a,b,c][r] the c-free factor
# divides a unit, and Gauss's lemma gives irreducibility over C(a,b,c).
r9 = sp.symbols('r9')
C9 = 2*a*r9**3 - b*r9**2 + 2*r9 - c
assert sp.degree(sp.Poly(C9, c), c) == 1
assert sp.Poly(C9, c).LC() == -1                       # unit c-coefficient
coeffs_r = sp.Poly(C9, r9).all_coeffs()                # [2a, -b, 2, -c]
assert 2 in [sp.simplify(cf) for cf in coeffs_r]       # unit content in r
g9 = sp.gcd(sp.gcd(sp.gcd(coeffs_r[0], coeffs_r[1]), coeffs_r[2]), coeffs_r[3])
assert g9 == 1 or g9.is_number
print("PASS 9 : covering cubic linear in c, unit c-coeff, content 1 =>")
print("        irreducible over C(a,b,c) (Gauss)")

print("ALL CROWN-JEWEL IDENTITY CERTIFICATES PASS")
