# The BCW / Drużkowski reduction, made explicit on Alpöge's counterexample

**Status: Stage A (cubic homogeneous, BCW) complete and machine-verified.
Stage B (cubic linear, Drużkowski) complete and machine-verified.**
Everything below is rebuilt from scratch and asserted by
`construct_verify.py` (runs green under `/usr/bin/python3`, sympy 1.12).

## 0. Ground truth (re-verified, not assumed)

The base counterexample F: C³ → C³,

    F₁ = (1+xy)³z + y²(1+xy)(4+3xy)
    F₂ = y + 3x(1+xy)²z + 3xy²(4+3xy)
    F₃ = 2x − 3x²y − x³z

satisfies (symbolically, exact): det JF ≡ −2, F(0,0,0) = 0, and the three
points

    p = (0, 0, −1/4),   q = (1, −3/2, 13/2),   r = (−1, 3/2, 13/2)

all map to (−1/4, 0, 0). So F is a Keller map that is not injective: the
Jacobian Conjecture is false, and by Bass–Connell–Wright (Bull. AMS 7,
1982; "BCW") and Drużkowski (Math. Ann. 264, 1983) the failure must
propagate to a cubic-homogeneous counterexample and to a cubic-linear one
in higher dimension. This directory carries those existence statements
out **explicitly**, with every intermediate map verified Keller and the
collision (non-injectivity witness) transported through **every**
elementary step.

## 1. The pipeline

Every step below is of the form  new = A ∘ (old × id_k) ∘ B  with A, B
explicit (affine-)elementary automorphisms of Jacobian determinant 1.
Hence, structurally:

* det J(new) = det J(old)  (as polynomials — chain rule);
* old(p) = old(q)  ⟹  new(B⁻¹(p,w)) = A(old(p),w) = new(B⁻¹(q,w)),
  so collisions transport explicitly (we always take w = 0).

### Stage 1 — linear normalization (dimension 3)

L₀ := JF(0) = [[0,0,1],[0,1,0],[2,0,0]], det L₀ = −2.
F̃ := L₀⁻¹ ∘ F = X + H,  H of order ≥ 2, degree 7, **det J F̃ ≡ 1**
(3×3, verified symbolically). Collisions: same three points, common image
now (0, 0, −1/4).

### Stage 2 — BCW Prop. (3.1): degree reduction to ≤ 3 (dim 3 → 39)

While some component F_i contains a monomial aM of degree d ≥ 4: write
aM = P·Q with 2 ≤ deg P, deg Q ≤ d−2 (we use the cost-optimal split
deg P = 2 for d ∈ {4,5}, 3 for d ∈ {6,7}), adjoin two variables and set

    B (BCW's "H") = (X, X_{n+1} + P, X_{n+2} + Q)          (elementary)
    A (BCW's "G") = (…, X_i − X_{n+1}·X_{n+2}, …)          (elementary)

so that the new map is

    F′ = ( …, (F_i − aM) − X_{n+1}Q − P X_{n+2} − X_{n+1}X_{n+2}, …,
           X_{n+1} + P,  X_{n+2} + Q ).

The pair (d, #top-degree monomials) drops lexicographically ⟹
termination. All created terms have degree in [2, d−1], so the map stays
X + (order ≥ 2) with identity linear part.

Collision transport per step: point p ↦ (p, −P(p), −Q(p)); common image
↦ (image, 0, 0).

**Outcome: 18 gadget steps, n₁ = 3 + 2·18 = 39.** The degree-≤3 map
F₁ = X + H₂ + H₃ on C³⁹ has identity linear part and H-terms only of
degrees 2 and 3.

**Symbolic Keller proof** (not just spot-checks): det J F₁ ≡ 1 is proved
mechanically by replaying the 18 gadgets backwards as exact column
eliminations: for each step, col_j −= (∂P/∂X_j)·col_a + (∂Q/∂X_j)·col_b
turns rows a, b into unit rows (asserted symbolically; the product rule
cancels exactly), Laplace deletes them, and after 18 replays the 3×3
det J F̃ = 1 finishes the proof. Column operations preserve determinants
exactly, so this is a complete symbolic proof, executed in sympy.

### Stage 3 — BCW §4 "doubling" (dim 39 → 78)

Split H = H₂ + H₃ (quadratic/cubic parts) and set, on C^{2·39},

    F′ = ( X + H₂(X) + Y,  Y − H₃(X) )
       = G(1) ∘ F₁[39] ∘ H(1),   G(1) = (X+Y, Y),  H(1) = (X, Y − H₃(X)).

This is the step that makes J(N) **nilpotent** (N = F′ − id): BCW run it
with a parameter T and conclude nilpotency from invertibility of
I + T·J(N) via their graded-ring Lemma (4.1). We do not assume this — we
prove the consequences symbolically in Stage 5.

Collision transport: p ↦ (p, H₃(p)); common image ↦ (image, 0).

### Stage 4 — BCW §4 homogenization (dim 78 → 79)

Adjoin T as an honest coordinate:

    G := L = ( X + T²Y + T·H₂(X),  Y − H₃(X),  T )   on C⁷⁹.

Every H-term is degree-3 homogeneous in (X, Y, T): T²Y_i and T·H₂ᵢ and
H₃ᵢ are all cubic; the linear part is the identity (T²Y_i is cubic, not
linear). Collisions live on the T = 1 slice: p ↦ (p, H₃(p), 1);
common image ↦ (image, 0, 1).

## 2. Stage A result

**G = id + H on C⁷⁹, H homogeneous of degree 3, det JG ≡ 1 identically,
JH nilpotent, with three explicit distinct points sharing one image.**

Dimension bookkeeping: 79 = 3 (base) + 2·18 (gadget pairs) + 39
(doubling partners) + 1 (homogenizing T).

Verification ledger (everything asserted by `construct_verify.py`):

| Claim | How verified |
|---|---|
| det JG ≡ 1 (all of C⁷⁹) | **Symbolic proof**: block elimination X-rows −= T²·(Y-rows) reduces JG to J F₁(TX) (Schur, D-block = I); entrywise identity J F₁(TX) = I + T·JH₂ + T²·JH₃ asserted; then the 18-step elimination replay at scaled argument TX ends at det J F̃(T·) = 1. Belt: exact det = 1 at 5 random integer points; exact det = 1 spot-checks after **every** gadget step. |
| JH nilpotent | **Symbolic proof**: same elimination chain applied to I + s·JH (s fresh) gives det(I + s·JH) ≡ 1 identically ⟹ char(JH) = λ⁷⁹ ⟹ JH⁷⁹ = 0 (Cayley–Hamilton over ℚ[X,Y,T]). Belt: JH¹⁹ = 0 exactly at a random integer point (index 19). |
| H homogeneous cubic | every monomial of every H_i asserted to have total degree exactly 3. |
| Non-injectivity | the three lifted points evaluated through G exactly (rational arithmetic); images equal, points pairwise distinct. Checked after **every** intermediate step as well. |

Colliding points, complete and exact (coordinates ordered
v1..v39, w1..w39, t; also in `G_map.txt` / `collisions.txt`):

    P = (0, 0, −1/4, 0×36 | 0×39 | 1)      (only v3 = −1/4 and t = 1 nonzero)

    Q: v = (1, −3/2, 13/2, −1, 351/16, −3, −117/8, 9/2, 27/8, −9, 27/8,
            9/4, −39/4, −3, −117/8, 1/2, −13/2, 9/4, 39/2, 9/4, −27/2,
            −6, 39/4, 3, 117/8, −39/4, 9/4, 9/4, 39/2, 9/4, 27/4, 21/2,
            −9/4, 3, 81/16, 1, −351/16, −9/4, 39/4)
       w = (−17/4, −45/8, 99/16, 1, −351/8, 3, 117/8, −9/2, −27/8, 0,
            −27/8, 0, 39/4, 0, 117/8, 0×24)
       t = 1

    R: v = (−1, 3/2, 13/2, 1, −351/16, 3, −117/8, −9/2, −27/8, −9,
            −27/8, 9/4, −39/4, −3, −117/8, 1/2, 13/2, 9/4, −39/2, 9/4,
            27/2, −6, −39/4, 3, −117/8, 39/4, −9/4, 9/4, 39/2, 9/4,
            27/4, 21/2, −9/4, 3, 81/16, 1, −351/16, −9/4, −39/4)
       w = (17/4, 45/8, 99/16, −1, 351/8, −3, 117/8, 9/2, 27/8, 0,
            27/8, 0, 39/4, 0, 117/8, 0×24)
       t = 1

Common image: G(P) = G(Q) = G(R) = (0, 0, −1/4, 0×75, 1).

## 3. Stage B — Drużkowski cubic-linear form (via Gorni–Zampieri pairing)

Construction (Gorni & Zampieri, *On cubic-linear polynomial mappings*,
Prop. 2.1, specialized — our A needs **no matrix inversion**):

1. **Waring step.** Using a·b² = ((a+b)³ + (a−b)³ − 2a³)/6 and
   a·b·c = ((a+b+c)³ + (a−b−c)³ − (a+b−c)³ − (a−b+c)³)/24, write the
   Stage-A cubic H as H(x) = −B₀ (D₀x)^{*3} (componentwise cube), where
   the rows of D₀ are the **r₀ = 347** distinct linear forms (deduped up
   to sign) in the 79 variables and B₀ is the 79 × 347 coefficient
   matrix (411 nonzero entries). **Verified symbolically, component by
   component.** rank D₀ = 79 (asserted via Gram determinant).
2. **Padding** (GZ Prop 2.1; needed because B₀ has zero rows — the
   t-component and every w-partner of a purely quadratic gadget
   component have H₃ = 0):

       B := [B₀ | I₇₉],   D := [D₀ ; 0],   C := [0 ; I₇₉],
       N := r₀ + 79 = 426,   BC = I₇₉.

3. **The map.** A := D·B = [[D₀B₀, D₀],[0, 0]] — explicit, **no matrix
   inversion anywhere** — and, writing X = (p, q) ∈ C³⁴⁷ × C⁷⁹,

       F(X) = X − (A X)^{*3}   on C⁴²⁶,
       F(p, q) = ( p − (D₀(B₀p + q))^{*3},  q ).

   Each component is literally X_k + (ℓ_k(X))³ with
   ℓ_k = −(row k of [D₀B₀ | D₀]), and ℓ_k = 0 on the padded block:
   the Drużkowski cubic-linear form.
   Pairing identities: AC = D, AM = 0 for M = [I; −B₀] spanning ker B,
   and B·F(Cx) = x − B₀(D₀x)^{*3} = G(x). rank D₀ = 79 gives
   ker A = ker B (GZ's pairing condition).
4. **Keller.** J_F = I − 3·diag((AX)^{*2})·D·B. Sylvester's identity
   det(I − UV) = det(I − VU) with U = 3·diag((DBX)^{*2})·D, V = B
   gives det J_F(X) = det(I₇₉ − 3·B₀·diag((D₀y)^{*2})·D₀) at y = BX,
   which is det JG(BX) ≡ 1 — riding on the Stage-A symbolic proof.
   The same chain with a fresh s gives det(I + s·J_{H_F}) ≡ 1, so
   J_{H_F} is nilpotent. (Sylvester's identity is the one cited, not
   machine-derived, ingredient; float spot-checks of det J_F at random
   426-dim points are run as a belt.)
5. **Collisions** (GZ Prop 3.1 direction "G not injective ⟹ F not
   injective", made explicit). For G(x₁) = G(x₂):

       X₁  = C·x₁ = (0, x₁),
       X₂′ = C·x₂ + F(C·x₁) − F(C·x₂) = ( (D₀x₂)^{*3} − (D₀x₁)^{*3},  x₁ ),

   then F(X₂′) = F(X₁) = (−(D₀x₁)^{*3}, x₁) — the key cancellation is
   B₀[(D₀x₂)^{*3} − (D₀x₁)^{*3}] = x₂ − x₁, which is exactly
   G(x₁) = G(x₂) — and B·X₁ = x₁ ≠ x₂ = B·X₂′ certifies X₁ ≠ X₂′.
   All three Stage-A points lift; the three lifted points and their
   common image are **verified exactly** by direct rational evaluation
   of F on C⁴²⁶.

**Outcome: a Drużkowski map F(X) = X − (AX)^{*3} on C⁴²⁶** (full data:
D₀, B₀ and the three colliding points in `F_druzkowski.txt`), Keller,
J_{H_F} nilpotent, with three explicit distinct colliding points.

## 4. What is proof and what is evidence

* Symbolic, machine-executed, exact: base det = −2; det J F̃ = 1; the
  18 unit-row assertions per elimination replay (×3 replays: plain, TX,
  sTX); the Schur block reductions; H = −B₀(D₀x)^{*3} **and its
  Jacobian** J_H = −3·B₀·diag((D₀x)^{*2})·D₀ (so the Sylvester chain's
  calculus step is machine-checked too); homogeneity; all collision
  evaluations; rank certificates.
* Mathematical theorems cited (not re-derived in code): chain rule /
  multiplicativity of det under the explicit compositions;
  Cayley–Hamilton; Sylvester's determinant identity; BCW Lemma (4.1)
  is *not* needed by our verification route (we prove its consequence
  det(I + s·JH) ≡ 1 directly).
* Spot-checks (falsification belts, not load-bearing): exact rational
  det = 1 after every gadget; 5 exact integer-point dets of JG;
  JH¹⁹ = 0 at a sample point; float dets for the Drużkowski map.

## 5. Dimension notes (not minimal, deliberately)

79 (Stage A) and 426 (Stage B) are what the *plain* BCW/GZ pipeline
yields on this F; they are upper bounds, not optima. Known slack, left
on the table for correctness-first reasons: (a) gadget factor sharing
(two gadgets with a common factor P can share one adjoined variable)
would cut several of the 36 gadget coordinates; (b) BCW's "linear in
each variable" refinement is skipped (not needed for cubic-homogeneous
form); (c) partial doubling (only components with quadratic part
strictly need partners) could shrink the 39 doubling coordinates;
(d) the Waring step uses per-monomial decompositions — a global change
of coordinates would reduce r₀ = 347; (e) the 79 padding coordinates
exist only to give B a right inverse — padding only the zero rows of B₀
would shrink that block. Anyone racing for a *smaller* N can start from
`construct_verify.py` and tighten these; the maps here are the first
fully verified ones, not the smallest.

## Files

* `construct_verify.py` — rebuilds everything from scratch; asserts all
  of the above; green under `/usr/bin/python3`.
* `G_map.txt` — the full 79 components of G and the three colliding
  points with common image.
* `collisions.txt` — the Stage-A collision data, one coordinate per line.
* `F_druzkowski.txt` — D₀, B₀ (sparse), and the three colliding points
  of the cubic-linear F with common image.
