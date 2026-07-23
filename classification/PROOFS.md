# The hand theorems

The 88 Gröbner-hard systems of the sweep fall into strata by weight
structure. Three theorems and a transfer lemma cover the six classes that
no machine route decided. Notation: for a weight vector w the family
consists of all maps F = (F₁,F₂,F₃) with F ∘ λ_t = λ_t ∘ F for
λ_t = (t^{w₁}x, t^{w₂}y, t^{w₃}z), components supported on the enumerated
degree ≤ 6 monomials. A Keller member has det JF ≡ k ≠ 0 constant (then
k = det JF(0) =: det0). A family is EMPTY if it has no non-injective
Keller member. W(f,g) := fg' − f'g; primes are d/dr unless stated.
Computations cited as (scripts/…) are exact symbolic checks or Gröbner
unit certificates; the index is at the end.

**Quotient principle.** Differentiating equivariance at t = 1 gives
JF·E = E∘F for the Euler field E = (w₁x, w₂y, w₃z), hence
F*(ι_E Ω) = det JF · ι_E Ω with Ω = dx∧dy∧dz. The form ι_E Ω descends to
each stratum's invariant quotient, so det JF ≡ k says the induced map on
the quotient plane scales a fixed 2-form by k. Each theorem below is this
principle in a chart adapted to its stratum.

---

## 1. The zero-weight stratum, w = (0, b, −c)

Here gcd(b,c) = 1, b,c ≥ 1, d := b + c. Invariants: x and u = y^c z^b.
Equivariance forces F = (f(x,u), y·A(x,u), z·B(x,u)). Define
Φ = (f, u·A^c·B^b) on the (x,u)-plane.

**Lemma 1.1.** Jac_{(x,u)}(Φ) = A^{c−1} B^{b−1} · det JF.
(From du' = c(yA)^{c−1}(zB)^b d(yA) + b(yA)^c(zB)^{b−1} d(zB), pulling
back dx∧du; checked symbolically for (c,b) ∈ {(1,1),(2,1),(1,2),(3,2)}.)

**Lemma 1.2.** At a point with A = B = 0, every row of JF lies in the
span of (1,0,0) and (0,u_y,u_z), so det JF = 0. A Keller member therefore
has no common zero of A and B.

**Theorem 1.3 (d = 2, w = (0,1,−1)).** The family is EMPTY.

*Proof.* A collision F(P) = F(Q) descends to Φ(πP) = Φ(πQ), π = (x, yz).
By Lemma 1.1 at (c,b) = (1,1), Jac(Φ) = det JF = k, so Φ is a planar
Keller map of degree ≤ 6, injective by Moh; hence πP = πQ = (x₀,u₀). If
u₀ ≠ 0, P ≠ Q on one fiber forces A(x₀,u₀) = B(x₀,u₀) = 0, contradicting
Lemma 1.2. If u₀ = 0: restricting Jac(Φ) ≡ k to u = 0 gives
f_x(x,0)A(x,0)B(x,0) = k, a product of univariate polynomials equal to a
nonzero constant, so each factor is constant; F restricted to {y = 0},
{z = 0}, and the axis is then affine and injective, and the
cross-configurations are impossible. ∎

**Theorem 1.4 (d ≥ 3).** Every degree ≤ 6 Keller member with
w = (0,b,−c), d ≥ 3, has A, B constant and f = c₀ + φx + h(u), and is
injective. All stratum families are EMPTY.

*Proof.* The caps force f = f₀(x) + f₁(x)u + f₂u² (f₂ only at d = 3),
A = A₀(x) + A₁(x)u, B = B₀(x) + B₁(x)u. The determinant identity
    det JF = f_x·(AB + u(cBA_u + bAB_u)) − u·f_u·(cA_xB + bAB_x)
has u⁰-layer f₀'A₀B₀ = k, so f₀ = c₀ + φx, A₀ = a₀, B₀ = b₀,
φa₀b₀ = k. Normalize P := A₁/a₀, S := B₁/b₀, α := φ. The u-layers give,
with W := (1+c)P + (1+b)S, V := PS, T' := cP' + bS', L := cP'S + bPS':
    (I) f₁' = −αW   (II) α(1+d)V − αW² − f₁T' = 0
    (III) f₁'(1+d)V − f₁L − 2f₂T' = 0   (IV) f₂L = 0 (d = 3 only).
Forcing P = S = 0, by cases on d:

d = 6: A₁ = B₁ = 0 by the caps.

d = 5: P, S constant. If f₁' ≠ 0, (III) gives V = 0, then (II) gives
W = 0; W = V = 0 with positive weights kills P, S. If f₁' = 0, (I) and
(II) do the same.

d = 4 (f₂ absent): if f₁PS ≢ 0, dividing (III) by f₁V gives
(1+d)f₁'/f₁ = cP'/P + bS'/S, so f₁⁵ = κP^cS^b; the degree count forces
f₁ constant; then f₁' = 0, (I) gives W = 0, (III) gives L = 0;
substituting S = −(1+c)P/(1+b) into L yields P·P' = 0, so P' = 0, then
T' = 0 and (II) gives V = 0, contradicting PS ≢ 0. If f₁ ≡ 0, (I) and
(II) kill W and V. If S ≡ 0, P ≢ 0 (or symmetrically): (I) and (II) give
f₁ = −α(1+c)²P²/(cP'), and differentiating forces 2(1+c) = c,
impossible.

d = 3: if f₂ ≠ 0, (IV) gives (P^cS^b)' = 0; if PS ≢ 0 then P, S are
constant and (III), (I), (II) force V = 0 then P = S = 0; if S ≡ 0,
(III) gives P' = 0 and (II) kills P. If f₂ = 0 and f₁PS ≢ 0:
f₁⁴ = κP^cS^b with 4·deg f₁ ≤ 6 forces deg f₁ ≤ 1; constant f₁ dies by
(I), (II); f₁ = γ(x−ρ) forces P = p(x−ρ)^i, S = s(x−ρ)^j with
ci + bj = 4, and (I) (W constant) kills the top coefficient in each
admissible (i,j). If f₂ = 0, f₁ ≡ 0, PS ≢ 0: (I) gives W = 0 and (II)
gives V = PS = 0. If f₂ = 0, S ≡ 0, P ≢ 0: polynomiality of f₁ forces
P = p(x−ξ)², and (I) demands 3(1+c)/(2c) = 1, impossible.

So A = a₀, B = b₀, f_x ≡ φ: F = (c₀ + φx + h(u), a₀y, b₀z) with
φa₀b₀ ≠ 0, visibly injective. ∎

Each branch of the case analysis also carries a Gröbner unit certificate
on the ladder ideal with the relevant coefficient inverted
(scripts/q1_forcing_check.py, four (c,b) shapes).

---

## 2. The affine GL(2) stratum, w = (a, −e, −e), a ≥ 2

Here gcd(a,e) = 1. Both y and z have weight −e, so the family contains
genuine shear automorphisms; rigidity cannot mean "all coefficients
zero". The correct statement:

**Theorem 2.1.** For a ≥ 2 under the degree ≤ 6 caps (classes {720,723},
{721,722}, {728,729}): every Keller member has F₁ = αx with α constant,
and is injective. The three classes are EMPTY.

*Setup.* The caps make every multiplier affine in the band:
F₁ = x(α + Σaⱼuⱼ) with uⱼ = x^e y^{a−j}z^j, and (F₂,F₃) = linear part
plus band terms. In the chart (r, u₀) = (z/y, x^e y^a) write
A = α + u₀ℓ(r), B = β₁ + β₂r + u₀B̂(r), C = γ₁ + γ₂r + u₀Ĉ(r), with
deg ℓ ≤ a, deg B̂, Ĉ ≤ a+1 and detM₀ = β₁γ₂ − β₂γ₁ ≠ 0. The quotient
principle gives
    Jac_{(r,u₀)}(Φ)·B² = det JF·A^{e−1}·B^a,  Φ = (C/B, u₀A^eB^a)
(scripts/q4_affine_ladder_check.py), whose u₀-layers form a polynomial
ladder with top layer ℓ·W(B̂,Ĉ)·B̂ = 0.

*Gauge.* Post-composing (F₂,F₃) with a constant unimodular matrix stays
in the family and preserves Keller, det0 and injectivity; it reduces the
branch Ĉ = ρB̂ to Ĉ ≡ 0.

*Proof.* Three branches from the top layer.

Branch B̂ ≡ 0: the reduced ladder gives
(S1-1) (1+e)detM₀ℓ + αB₀Ĉ' − (1+a)αβ₂Ĉ = 0 and
(S1-2) B₀[(1+e)ℓĈ' − eℓ'Ĉ] − (1+e+a)β₂ℓĈ = 0.
If β₂ = 0, B₀ is constant and (S1-2) forces ℓ^e = κĈ^{1+e}; the degree
count with (S1-1) kills ℓ. If β₂ ≠ 0: at the root ξ of B₀, (S1-2) gives
ℓ(ξ)Ĉ(ξ) = 0 while (S1-1) makes ℓ(ξ) and Ĉ(ξ) proportional, so both
vanish; dividing by (r−ξ) decrements the ladder constants and repeats;
deg ℓ ≤ a exhausts: ℓ ≡ 0.

Branch W(B̂,Ĉ) = 0, B̂ ≢ 0: gauge to Ĉ ≡ 0. The ladder gives
(M1) detM₀[(1+e)ℓB₀ + aαB̂] + αB₀(γ₂B̂ − C₀B̂') + aαβ₂C₀B̂ = 0 and,
dividing the next layer by B̂,
(M3) e·ℓ'/ℓ − (1+e)·B̂'/B̂ = −(1+e+a)·γ₂/C₀ as rational functions.
If γ₂ ≠ 0: matching residues in (M3) at the root η of C₀ gives
mult_ℓ(η) < mult_B̂(η), using deg B̂ ≤ a+1 < 1+e+a; reading (M1) at order
mult_ℓ(η) in (r−η) leaves (1+e)detM₀·lc·B₀(η) = 0, forcing B₀(η) = 0,
impossible since Res(B₀,C₀) = detM₀ ≠ 0. If γ₂ = 0 (so
β₂γ₁ = −detM₀ ≠ 0): (M3) gives ℓ^e = κB̂^{1+e}; in (M1) the two B̂-blocks
cancel exactly, leaving ℓ ∝ B̂', and (B̂')^e ∝ B̂^{1+e} forces
deg B̂ = −e < 0: B̂ constant, ℓ ≡ 0.

With ℓ ≡ 0 in every branch: A ≡ α, F₁ = αx pins x, and
det JF = α·det ∂(F₂,F₃)/∂(y,z) ≡ k makes every x-slice a planar Keller
map of degree ≤ 6, injective by Moh; collisions across slices die on the
first coordinate. ∎

Ladder extractions are checked generically at (a,e) = (3,1), (3,2),
(4,1); the B̂ ≡ 0 branch carries unit certificates at (3,1) and (3,2),
and at (4,1) rests on the uniform descent above.

---

## 3. The class {734, 735}, w = (1, −1, −1)

**Theorem 3.1.** The w = (1,−1,−1) family contains no non-injective
Keller member. With Lemma 4.1, class {734,735} is EMPTY and the
degree ≤ 6 equivariant classification is complete.

*Chart.* With u = xy, v = xz the family is exactly F₁ = x·A, F₂ = y·B,
F₃ = y·C in the chart (r,u₀) = (z/y, xy):
A = α + u₀Ah(r) + u₀²Ac(r), B = B₀ + u₀Bh + u₀²Bc,
C = C₀ + u₀Ch + u₀²Cc, with degree profiles (1|2), (linear|2|3),
(linear|2|3). This is a coefficient bijection with the enumerated
24-coefficient family, and
    det JF · B = D·G − u₀·E·(AB)_r,
D = BC_r − CB_r, E = BC_{u₀} − CB_{u₀}, G = (u₀AB)_{u₀}, holds on the
raw family (scripts/q8_cal_chart.py); det0 = A(0)·detM₀.

*Reduction to constancy of A.* An injective Keller member is an
automorphism (Ax–Grothendieck; polynomial inverse by
Białynicki-Birula–Rosenlicht) with equivariant inverse, so the induced
quotient map Φ = (A·L₁, A·L₂), where L₁ = uB₁ + vB₂ and L₂ = uC₁ + vC₂,
is an automorphism of C² whose Jacobian det JF·A is then constant: A is
constant. Conversely if A ≡ α, then Jac(L₁,L₂) = k/α = detM₀ ≠ 0, so
(L₁,L₂) is a planar Keller pair of degree ≤ 3, injective by Moh; a
collision F(p) = F(p') pushes down to Φ, so π(p) = π(p'), F₁ = αx gives
x = x', and the fiber and boundary cases are immediate. **The theorem is
therefore equivalent to: the Keller equations force A constant.** A leaf
of the analysis below admitting a nonconstant-A solution would be an
explicit degree ≤ 6 counterexample; none does.

*Gauge.* Source GL(2) on (y,z) and target GL(2) on (F₂,F₃) preserve the
family, Keller, det0 ≠ 0, and injectivity, and act transitively on
invertible linear parts: normalize (B₀, C₀) = (1, r), detM₀ = 1, k = α.
The residual conjugations act on chart blocks by A_m ↦ w^m A_m(mob) and
(B,C)_m ↦ w^{m+1}/Δ·(mix)(B,C)_m(mob), mob = (s+tr)/(p+qr), w = p+qr.
Two instances are used below: the shear (r ↦ r+s; C-blocks pick up
−s·B-blocks) and the swap (B_m and C_m exchange with reversal
f(r) ↦ r^{deg}f(1/r)). Formulas: scripts/q12_gauge.py.

*Ladder.* Keller is equivalent to the vanishing of layers N₁..N₈ of
D·G − u₀E(AB)_r − αB (convolution table in scripts/q8_layers.py). Two
enslavements: N₁ = 2Ah + αP with P := 2Bh + Ch' − rBh' (the r² term
cancels, deg P ≤ 1), so Ah = −αP/2; and N₂ is linear in Ac with unit
coefficient. The top layers factor as
    N₈ = 5·Ac·Bc·W(Bc,Cc),   N₇|_{Ac=0} = −2αP·Bc·W(Bc,Cc),
and "A constant" is equivalent to P ≡ 0 and Ac ≡ 0.

*The tree.* N₈ = 0 and, on Ac = 0, N₇ = 0 partition all Keller
solutions:

(α₁) Ac ≡ 0, P ≡ 0: A ≡ α, injective by the reduction. (These are the
genuine automorphism strata, e.g. the shears (x, y + czv, z).)

(a2) Ac ≡ 0, Bc ≡ 0, P ≢ 0: contradiction. Write h := Bh, c := Ch,
K := Cc, T := rh' − c'. Layer N₅* = h(K(Ph)' − 3PW(h,K)). If h ≡ 0,
layers N₃* and N₂* give W(c',K) = 0 and then a coefficient
contradiction. If K ≡ 0, N₄* forces c = τh, and N₂* in the shifted
variable forces h = h₂(r−τ)², making P ≡ 0. Otherwise
K(Ph)' = 3PW(h,K) integrates to K³ = κPh⁴, whose degree and root
multiplicity matching under the caps admits only constant blocks or the
coincidence P = p₁(r−ξ), h = h₂(r−ξ)², K = k₃(r−ξ)³; constants die in
N₄* and N₂*; the coincidence dies at the r⁵ coefficient of N₄, which
reduces to p₁h₂k₃ = 0. The whole branch also carries a single Gröbner
certificate: the branch ideal plus the inverse of either coefficient of
P is the unit ideal (scripts/q10_leafchecks.py), and the coincidence
leaf a certificate with the root position ξ generic
(scripts/q11b_xi_certs.py).

(a3) Ac ≡ 0, Bc ≢ 0, W(Bc,Cc) = 0, P ≢ 0: then Cc = ρBc; the shear
with s = ρ makes Cc ≡ 0 and the swap produces an (a2) configuration,
closed by (a2). A direct certificate with ρ generic also closes the
branch (scripts/q11_beta.py).

(β) Ac ≢ 0, Bc ≡ 0: contradiction. Layer N₆ = h(4AcW(h,K) − K(Ach)').
If h ≡ 0: N₄ gives 3AcK' = 2KAc', so K³ = κAc² forces K, Ac constant
(or K ≡ 0), and N₃, N₂ then force Ac = 0 either way. If K ≡ 0: N₅
forces c = τh, N₄ forces h and Ac proportional to (r−τ)², and N₂
reduces to 6·Ac = 0. Otherwise K⁴ = κAc·h⁵ admits only constants
(killed by N₅ and then N₄ = 4Ac·h² ≠ 0) or the coincidence
Ac = a₂(r−ξ)², h = h₂(r−ξ)², K = k₃(r−ξ)³, which is closed by a
Gröbner certificate with ξ generic (scripts/q11_beta.py,
q11b_xi_certs.py); for this last sub-leaf the written chain stops at a
relation system and the certificate decides. Every other leaf has both
the written chain and a certificate.

(γ) Ac ≢ 0, Bc ≢ 0, W(Bc,Cc) = 0: the shear and swap produce a (β)
configuration, closed by (β).

Every leaf except (α₁) is contradictory, so Keller forces A constant,
and Theorem 3.1 follows. ∎

The displayed intermediate identities are checked in
scripts/q14_handsteps.py; branch-specialized layers in
scripts/q9_branches.py. No collision certificate in the family returned
NONUNIT.

---

## 4. Transfer between mirror systems

**Lemma 4.1.** Let σ ∈ S₃ act simultaneously on source coordinates,
target components, and the weight vector. The induced map on coefficient
families is a bijection sending Keller members to Keller members (det J
transported by conjugation, det0 preserved) and collisions to
collisions; EMPTY transfers along σ. A global sign flip of the weights
fixes each family, since the supports satisfy S(w) = S(−w) exactly.

System f0735 (w = (1,−1,1)): the sign flip identifies its family with
that of (−1,1,−1), and σ = (12) transports the latter onto f0734's
family; the support bijection is checked against the two job files with
negative controls (scripts/q13_transfer_735.py). Hence f0735 is EMPTY by
Theorem 3.1.

---

## Computation index

| Script | Checks |
|---|---|
| q1_identity_check.py, q1_forcing_check.py | §1 identity, ladder, forcing certificates |
| q4_affine_ladder_check.py | §2 ladder extractions |
| q5_a1_ladder.py | §3 layer structure before normalization |
| q8_cal_chart.py | §3 chart bijection, ladder identity, det0 |
| q8_layers.py, q8_layers.srepr | §3 layer table; N₁, N₂ enslavements; N₈, N₇ factorizations |
| q9_branches.py | §3 branch-specialized layers |
| q10_leafchecks.py | §3 branch master forms; (a2) certificates |
| q11_beta.py, q11b_xi_certs.py | §3 (β), (a3) certificates; generic-ξ coincidence certificates |
| q12_gauge.py | §3 gauge action, shear and swap instances |
| q13_transfer_735.py | §4 transfer witness, negative controls |
| q14_handsteps.py | §3 intermediate identities |
| class_fold.py | the 741 → 371 class fold from the job-file supports |
| iso_audit.py, iso_transfer.py | permutation witnesses and transfer checks for the hard systems |
