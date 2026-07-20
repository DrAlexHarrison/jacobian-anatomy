#!/usr/bin/env python3
"""The second counterexample (GPT-5.6 Sol via @hari65535): G : C^4 -> C^4.

Verifies:
  1. det JG == 4 identically.
  2. The four distinct points (-2,0,1,0), (-1,0,1,-2), (1,-2,-7,14),
     (2,-1,0,3) all map to (-1,-4,8,-4).  Since G is etale (Keller), a
     four-point fiber forces generic degree >= 4.
  3. C*-equivariance with source weights (1,-1,-3,-2) -> target
     (-3,-2,-1,1): the same reversed-weight pattern as the 3D map.
     The equivariant ansatz is 2-for-2 across counterexamples.
"""
import sympy as sp

x, y, z, w, t = sp.symbols('x y z w t')
u = 1 + x*y
G = [-(u)**2*z - y**3*u,
     2*x*u*z + u**3*w + y**2*(7 + 12*x*y + 4*x**2*y**2),
     2*x**2*z + 3*x*u**2*w + 2*y*(1 + 10*x*y + 6*x**2*y**2),
     2*x - 4*x**2*y - x**3*w]
V = (x, y, z, w)

# --- 1. determinant -------------------------------------------------------
assert sp.expand(sp.Matrix(G).jacobian(V).det()) == 4
print("PASS 1 : det JG == 4 identically")

# --- 2. the four-point collision -----------------------------------------
pts = [(-2, 0, 1, 0), (-1, 0, 1, -2), (1, -2, -7, 14), (2, -1, 0, 3)]
assert len(set(pts)) == 4
for p in pts:
    img = tuple(sp.expand(g.subs(dict(zip(V, p)))) for g in G)
    assert img == (-1, -4, 8, -4), (p, img)
print("PASS 2 : four distinct points -> (-1,-4,8,-4); etale => generic degree >= 4")

# --- 3. equivariance ------------------------------------------------------
sub = {x: t*x, y: y/t, z: z/t**3, w: w/t**2}
sc = [sp.simplify(g.subs(sub)) for g in G]
assert sp.simplify(sc[0] - G[0]/t**3) == 0
assert sp.simplify(sc[1] - G[1]/t**2) == 0
assert sp.simplify(sc[2] - G[2]/t) == 0
assert sp.simplify(sc[3] - t*G[3]) == 0
print("PASS 3 : G(tx, y/t, z/t^3, w/t^2) = (G1/t^3, G2/t^2, G3/t, t*G4)")
print("        reversed-weight C* symmetry: the equivariant pattern is 2-for-2")

print("ALL SECOND-COUNTEREXAMPLE CERTIFICATES PASS")
