# Vetter report: C*-equivariance claim-set + sheet parametrization

**Verdict: CONFIRMED** (both claim-sets; one scoping caveat on the image
corollary, one wording requirement on the lifting law — see the addendum).

Scope note: the NOTE.md referee pass was reassigned to an external session by
the orchestrator (2026-07-20); this report covers job 1 only — the
equivariance claim-set and the later-added sheet-parametrization claims.

- Auditor: vetter (pre-publication referee agent), 2026-07-20.
- Script: `/home/alex/code/jc/audit/vetter_equivariance.py` — written from scratch,
  exact arithmetic only, zero reuse of `verify/03_new_facts.py` (the artifact under
  audit). 35/35 checks PASS on first run (sympy 1.12, /usr/bin/python3).
- Map formulas cross-checked against `verify/01_basic.py` to rule out relay typos.

## What was verified (all symbolic, all exact)

1. **Equivariance identity** `F(t·x, t⁻¹y, t⁻²z) = (t⁻²F₁, t⁻¹F₂, tF₃)` holds
   identically as a Laurent-polynomial identity in t — valid for every t ∈ C*,
   not just sampled values.

2. **The weights are forced.** Requiring each component of F to be
   quasi-homogeneous for unknown weights (p,q,r) yields a linear system on the
   monomial exponent vectors whose nullspace is exactly 1-dimensional, spanned by
   (1,−1,−2). Concretely, within F₁ the monomials z and y² force r = 2q and z,
   xyz force p+q = 0; everything else is consistent. Induced target weights are
   (−2,−1,1). So the C*-structure is intrinsic to F, not a choice.

3. **Γ is a single closed orbit.** τ_s(Γ(t)) = Γ(st) identically, the orbit of
   Γ(1) = (4/27, 4/3, 1) sweeps Γ bijectively (the c-coordinate is the
   parameter, so the parametrization is injective and the stabilizer is
   trivial: orbit ≅ C*). Degenerate-limit probe: Γ equals the affine variety
   V(bc − 4/3, 12a − b²) **as a set, both inclusions verified** (the system has
   the unique solution family (4/27c², 4/3c, c), and has no solution at c = 0).
   Hence Γ is Zariski-closed in C³ already — the t→0 and t→∞ limits leave every
   affine chart and add no closure points. "C³ minus one orbit" is a correct
   restatement of "C³ ∖ Γ".

4. **All named strata are τ-invariant.**
   - D = 27a²c² − 18abc + 16a + b³c − b² is quasi-homogeneous of weight −2
     (D∘τ_t = t⁻²D), so V(D) is invariant.
   - S = 27ac² − 9bc + 8 is genuinely invariant (S∘τ_t = S).
   - Γ is an orbit; the punctured a-axis is invariant by inspection.
   - Consistency: Γ ⊂ V(D) and Γ ⊂ V(S), both verified identically in t.

5. **The "1-parameter family of triple fibers" = τ-orbit of the special
   fiber.** τ_t(−1/4, 0, 0) = (−1/(4t²), 0, 0) sweeps exactly the punctured
   a-axis (t² ↦ −1/(4s) solvable for every s ≠ 0). The special fiber over
   (−1/4, 0, 0) is exactly {(−1, 3/2, 13/2), (0, 0, −1/4), (1, −3/2, 13/2)}:
   quotient dimension 3 by Gröbner, and since det J ≡ −2 the fiber scheme is
   reduced, so 3 is the exact point count, complete, no missed solutions.
   By equivariance F⁻¹(τ_t p) = σ_t(F⁻¹(p)) with σ_t an automorphism of C³, so
   every fiber over the punctured a-axis has exactly 3 points. ⚠️ One wording
   caution for the note: **triple fibers are not special** — generic fibers are
   also 3 points (verified at (2,1,1) and (1,1,1); consistent with 01_basic.py's
   generically-3:1 finding). What is special about the a-axis family is only
   that it is the announced fiber's orbit; the note should not imply 3-point
   fibers are rare.

## The image corollary — status upgraded, one boundary drawn

- **"Γ misses the image" is now PROVEN, not sampled.** Gröbner basis of
  ⟨F₁−4/27, F₂−4/3, F₃−1⟩ over Q is [1] — a Nullstellensatz unit certificate,
  exact over Q, hence no solutions over C. Equivariance transports this single
  certificate to every point of Γ: if F(v) = Γ(t) then F(σ_{1/t}(v)) = Γ(1),
  contradiction. Independent re-certification at Γ(2) and Γ(−1/3) also returned
  [1]. So Γ ∩ Im(F) = ∅ rigorously — the "finite sample points" worry is
  retired for this half of the image claim.
- **The other half ("everything off Γ is hit") is NOT established by
  equivariance** and must carry its own proof in the note (the elimination /
  degeneration argument). Equivariance only reduces it to one point per orbit —
  still a 2-parameter family of orbits. My probes are all consistent with it:
  fibers over generic points have 3 points, over V(D)∖Γ sample points
  ((−16/27,0,1), (0,0,1)) and the origin exactly 1 point, over Γ zero. If the
  note words the image characterization as a consequence of equivariance alone,
  that is an overclaim; as a claim proven elsewhere and *organized* by
  equivariance, it is fine.

## Bonus facts verified while arming the note referee pass

- **Eliminant relation:** D(F)·x³ + (4 − 3F₂F₃)·x − 2F₃ ≡ 0 identically in
  Q[x,y,z] — the load-bearing depressed cubic (no x² term) is a true identity;
  every fiber point's x-coordinate satisfies E(a,b,c;x) = D·x³ + (4−3bc)x − 2c
  at its image values.
- **disc_x(E) = −4·D·S².** So the discriminant class of the fiber cubic in
  C(a,b,c)*/(squares) is exactly −4D, and the surface S = 27ac²−9bc+8 enters
  precisely as the square factor — this is presumably why S is a stratum, and
  it is the exact fact the S₃-monodromy argument needs.
- **V(D) ∩ {bc = 4/3} = Γ exactly:** 27c²·D|_{b=4/(3c)} = (27ac² − 4)², a
  perfect square vanishing only at a = 4/(27c²). Together with E: on
  {bc = 4/3, c ≠ 0, D = 0} the cubic E degenerates to the nonzero constant
  −2c, which is the clean structural reason the fiber over Γ is empty.

## Refutation attempts that failed (i.e., the claim survived)

1. Sought a second nullspace direction in the weight system — none (rank is 2).
2. Sought extra components or c=0 points in V(bc−4/3, 12a−b²) that would make
   "Γ closed / one orbit" false — none.
3. Sought a τ_t direction convention mismatch (source vs target weights) — the
   identity as relayed is the one that holds.
4. Sought a nonzero D∘τ_t − t⁻²D or S∘τ_t − S residue — both vanish identically.
5. Sought a 4th point or non-reducedness in the special fiber — quotient
   dimension is exactly 3 and étaleness forces reducedness.
6. Sought a Γ-point with a nonempty fiber (would refute the announcement) at
   three independent rational parameters including a negative one — all empty.

---

# Addendum: sheet parametrization / wall identity (added on orchestrator request)

**Verdict: CONFIRMED — and upgraded to a complete proof of the lifting law.**
Script: `/home/alex/code/jc/audit/vetter_wall_identity.py` (23/23 PASS, exact
symbolic, written from scratch; no reuse of verify/03 or verify/04).

## The claims, as verified

1. **Covering identity** (P0): with r = x/(1+xy),
   2F₁r³ − F₂r² + 2r − F₃ ≡ 0 identically in Q(x,y,z). So the r-coordinate of
   every fiber point over q = (a,b,c) with 1+xy ≠ 0 is a root of the covering
   cubic C_q(r) = 2ar³ − br² + 2r − c.

2. **Wall identity** (P1): CONFIRMED in the stronger exact-cofactor form
   (1 − r·y(r)) − C′(r)/2 = −3·C(r)/(2r) identically, with
   y(r) = −(br² + 3c − 6r)/(2r²). On C(r) = 0, r ≠ 0 this collapses to the
   claimed 1 − r·y(r) = C′(r)/2. Consequence, proved not sampled: the
   section's x = r/(1 − r·y(r)) has a pole exactly at multiple roots — the
   3-lines-of-algebra derivation checks out.

3. **Section** (P2, P3): x(r) simplifies to 2r²/(br² − 4r + 3c); z(r) solved
   from F₃ = c. Stronger than relayed: **F₂ = b and F₃ = c hold EXACTLY on
   the section** (identically in Q(a,b,c)(r), no reduction needed); C_q(r) = 0
   is precisely the remaining condition F₁ = a (F₁(sec) − a is a unit multiple
   of C_q). r = x/(1+xy) and x = r/(1−ry) verified mutually inverse.

4. **disc_r(C_q) = −4D exactly** (P4). This is the clean home of the "−4D"
   in the monodromy argument; the −4D·S² I found earlier for the x-eliminant E
   is the same class mod squares (S² is the projection artifact).

5. **Γ = triple-root locus, both directions** (P5): over Γ(t) the cubic is
   the perfect cube (8/27t²)(r − 3t/2)³, and coefficient-matching
   C = 2a(r − ρ)³ solves uniquely to Γ(2ρ/3). So "fiber empty ⟺ all sheets
   walled simultaneously ⟺ q ∈ Γ" is now structural, not just sampled.

## The lifting law: exact statement, and a complete proof from the verified identities

Relayed form: "#fiber = #simple roots of C". As an **affine** statement this
is FALSE precisely on {a = 0, b ≠ 0} (degree drop): fiber over (0,1,0) has 3
points but the degree-2 cubic has only 2 simple affine roots. It is TRUE in
the **projective** reading — simple roots of the homogenized cubic in P¹,
where r = ∞ has multiplicity 3 − deg(C_q) — verified at a 14-point exact
battery spanning every stratum type (generic; V(D)∖Γ; Γ ×2; a-axis; a=0 with
b≠0, with bc=1, with b=0; c=0 generic and on V(D); origin). Any published
statement must either homogenize or carve out a = 0 with an explicit +[b≠0]
term.

The verified identities assemble into a full proof, valid at every q:

- Every fiber point has an r ∈ P¹: r = x/(1+xy) if 1+xy ≠ 0 (a root of C_q by
  P0); r = ∞ if 1+xy = 0, forcing a = 0, b = 2/x ≠ 0 (P8).
- **r = 0 sheet** ⟺ x = 0 ⟺ c = 0: F(0,y,z) = (z+4y², y, 0) is an
  isomorphism onto {c = 0} (P7) — exactly one such point, and r = 0 is always
  a simple root when it is a root (C′(0) = 2 ≠ 0 unconditionally).
- **r = ∞ sheet** ⟺ 1+xy = 0: exists iff a = 0, b ≠ 0 (unique point:
  x = 2/b, z solved from c; P8) iff ∞ is a simple root of the homogenization
  (multiplicity 3 − deg = 1 iff a = 0, b ≠ 0).
- **Affine r₀ ≠ 0, simple**: the section is defined (wall = C′(r₀)/2 ≠ 0,
  denominators nonzero) and lands in the fiber (P2), giving existence.
- **Uniqueness** (P9, the closing lemma): y ≡ y(r)|_{b=F₂, c=F₃, r=x/(1+xy)}
  holds identically on the source — so any fiber point with r ∉ {0,∞} IS
  section(r): at most one point per root.
- **Affine r₀ ≠ 0, multiple**: a fiber point there would satisfy
  1 − r₀y(r₀) = 0 by P9 + the wall identity, while x(1 − r₀y) = r₀ forces
  0 = r₀ ≠ 0 — contradiction. No point.

Hence #F⁻¹(q) = #{simple roots of C_q in P¹} for ALL q, proved from machine-
verified identities plus the four displayed one-line arguments. Combined with
P5 (triple root ⟺ Γ) and disc = −4D, this yields the full stratification
3/1/0 and the exact image C³ ∖ Γ as theorems — closing the "onto-half"
gap I flagged in the main report, provided the note carries this argument (or
an equivalent) rather than sample points alone.

## Refutation attempts that failed (addendum)

1. Sought a nonzero remainder in F₁∘section mod C_q — zero; and F₂, F₃ exact.
2. Sought a case where a simple root fails to lift (denominator vanishing at a
   simple root) — impossible: the denominator IS 2r·(wall) = r·C′(r).
3. Sought fiber points not covered by the three charts (section / x=0 /
   1+xy=0) — none exist: the charts partition by r ∈ C*∪{0}∪{∞}.
4. Sought a battery point where the projective count mismatches the exact
   Groebner fiber count — none in 14 stratum-spanning probes; the affine
   reading duly failed where predicted ((0,1,0), (0,2,3), (0,1,1)), confirming
   the machinery can detect miscounts.
