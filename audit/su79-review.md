# Review notes: the SU(79) Mathieu chain

**Reviewer:** a separate session that produced none of the work below. **Date:** 2026-07-20 ~18:30 MDT (stamped from `date`).
**Objects reviewed:**
`programs/p3-mathieu/{PROGRAM.md, su79_certificate.py, su79_run.log,
su79_inverse_ckpt.json}`; consumed chain `druzkowski/{REDUCTION.md, G_map.txt}`
(prior internal check relied on, consumption spot-checked); primary sources
`refs/mathieu1997.pdf` (all 17 pages, visual read) and `refs/zwart.txt` (all
load-bearing proofs read in full).

## Summary

**The mathematical result stands: the Mathieu conjecture
is FALSE for SU(79), and for SU(N) for every N ≥ 79.**

**Publishability grade: PUBLISHABLE AFTER REVISIONS (major revision of the
bridge section and one claim de-escalation; no change to the result).**
The compute layer is flawless and now double-verified by an independent
engine. The representation-theoretic bridge as WRITTEN cites a proof with a
genuine gap (Zwart Thm 4.23) and overclaims explicitness of the witness pair;
both are repairable by citing Mathieu's original argument, which I have
walked step-by-step and which is sound. Required revisions R1-R8 below.

---

## 1. What is PROVED (my independent reconstruction of the correct chain)

Let L = id + H on C^79 (druzkowski/G_map.txt), H homogeneous cubic.
Machine-verified facts, each now confirmed by TWO independent implementations
(their `su79_certificate.py` + my `audit/su79_crosscheck.py`, disjoint
code, exact rational arithmetic):

- (F1) Every H-term is homogeneous of degree exactly 3; 51/79 components
  nonzero; 124 cubic terms. [IN-SESSION VERIFIED]
- (F2) det J L = 1: identically, by druzkowski's symbolic elimination-replay
  proof (checked separately); spot-checked here at 10 fresh random rational
  points + P, Q, R, and a deliberately broken coefficient is caught by the
  same spot-check (control has teeth). [IN-SESSION VERIFIED consumption]
- (F3) L(P) = L(Q) = L(R) = (0,0,-1/4,0,...,0,1) with P, Q, R pairwise
  distinct exact rationals: L is NOT injective. [IN-SESSION VERIFIED]
- (F4) The formal inverse F of L has nonzero homogeneous components at every
  odd degree 3..11; my engine's F agrees with their checkpoint TERM-BY-TERM
  (77,934 stored terms), all five logged witness coefficients match exactly,
  and the non-tautological two-sided composition F(L(t·u)) = t·u mod t^13
  holds along 3 random directions. [IN-SESSION VERIFIED]

Mathematical chain from (F1)-(F4) to the headline, with the citations I
verified from the PRIMARY sources:

1. **Non-injective ⇒ non-polynomial formal inverse.** If the formal inverse F
   were polynomial, both formal identities L∘F = id and F∘L = id would be
   polynomial identities, and F∘L = id forces injectivity of L,
   contradicting (F3). Hence F^{(2k+1)} ≠ 0 for infinitely many k (a power
   series with finitely many nonzero homogeneous components is a polynomial).
2. **Keller ⇒ hypothesis side.** Abhyankar/BCW inversion formula
   [Mathieu 1997, Formula 4.4(i), p. 273; requires j = 1, i.e. (F2)]:
   1 = Σ_α ∂^(α)h^α with h = −H. Extracting the degree-2k homogeneous part:
   div(Q_elt^k) = 0 for all k ≥ 1, where Q_elt = Σ H_i ⊗ ∂_i ∈ S³C^79 ⊗
   (C^79)*, degree-one element of the graded SL(79,C)-algebra
   A' = ⊕_k S^{3k}C^79 ⊗ S^k(C^79)*. By multiplicity-one [Mathieu Prop 3.3,
   proof read: S^{3k}V ⊗ S^kV* is multiplicity free since S^kV* is weight
   multiplicity free, Lemma 3.2(i)]: **(Q_elt^k)_{2kω₁} = 0 for ALL k.**
3. **Witness side.** Formula 4.4(ii) + Prop 3.4(ii)-(iii) [Mathieu p. 272-274]:
   given div(Q_elt^k) = 0, the projection Φ satisfies Φ(Q_elt^k/k!) =
   Σ_i F_i^{(2k+1)} ⊗ ∂_i, and Φ vanishes on every isotypic component except
   L(2kω₁) and L(2kω₁ + θ), θ = ω₁ + ω₇₈. Hence F^{(2k+1)} ≠ 0 ⇒
   **(Q_elt^k)_{(2k+1)ω₁+ω₇₈} ≠ 0**, for infinitely many k by step 1.
   (Weight bookkeeping at N = 79 checked: τ = (d−1)ω₁ = 2ω₁, kτ + µ with
   µ = θ = ω₁ + ω₇₈ gives (2k+1)ω₁ + ω₇₈; θ = highest root, valid for
   N ≥ 3. Focus point (c): CORRECT.)
4. **The bridge (REPAIRED, see Finding B).** Steps 2-3 refute exactly the
   conclusion of Mathieu's Proposition 2.2(ii) [p. 268-269] for the graded
   G-algebra A' with f = Q_elt ∈ A'_1, τ = 2ω₁, µ = θ. Hence the conjecture
   C(A' ∗ A(2ω₁)) is FALSE for G = SL(79,C). By Corollaries 1.3 + 1.7
   [p. 266-267; per-group at proof level: for a FIXED G, both proofs consume
   the Main Conjecture only for K = a maximal compact subgroup of that same
   G], MC(SU(79)) implies C(B) for every SL(79,C)-algebra B. Contrapositive:
   **the Mathieu conjecture is false for K = SU(79).** By definition of MC
   being false, there EXIST K-finite f, g on SU(79) with ∫ f^n dk = 0 for
   all n ≥ 1 and ∫ f^n g dk ≠ 0 for infinitely many n.
5. **N ≥ 79.** Pad L to C^N, N > 79, by identity coordinates (H_i = 0 for
   i > 79): still cubic-homogeneous, Keller, non-injective (same P, Q, R
   padded by zeros). The whole chain reruns verbatim at N. (Currently
   recorded NOWHERE in the artifacts, revision R5.)

I verified Mathieu's Prop 2.2(ii) proof line by line (p. 269): the
Baire-category choice of ξ ∈ ∩_n U_n over the uncountable field C (each U_n
dense open because L(nτ*) is spanned by n-th powers and isotypic components
of submodules inject into ambient isotypic components), the Lemma 2.1
surjection giving a FIXED detectable type ν for all n ≥ N, and the trivial
component computation (L^n)_triv = [(f^n)_{nτ} ⊗ ξ^n]_triv. Sound. Cor 1.3's
finitely-generated-subfield embedding into C and Lemma 1.2's ∫ = triv
projection: sound. Focus point (b), "non-polynomial inverse ⇒ infinitely
many nonzero degrees", is trivially correct once STATED as in step 1;
revision R4 fixes the garbled in-file justification.

## 2. FINDINGS (defects; none overturn the result)

**FINDING A (process, serious): su79_run.log was produced by a DIFFERENT
script than the one on disk, and its negative control is vacuous.**
The logged control (log lines 14-21) is a C^5 triangular automorphism whose
inverse is nonzero at EVERY odd degree 3..15 through the scan ceiling, then
declares "F^(d)==0 for all odd d>15 up to 15: TERMINATES": an empty claim.
As logged, the control is indistinguishable from L. (A C^5 chain terminates
at degree 3^4 = 81, far past the 15-ceiling, the control was miscalibrated.)
The script on disk (mtime 03:37, newer than the log's 03:35) already contains
a rewritten C^3 control that genuinely terminates at degree 9 inside the
15-window, but no log of it existed. **Remediated during this review:**
fresh clean-room rerun of the on-disk script (no stale checkpoint), log
committed (now `../mathieu-su79/su79_rerun.log`); an independent
CONTROL D adds a 79-variable in-dimension automorphism (terminates at
degree 3) and CONTROL E checks the C^3 control against its closed-form
inverse. Focus point (d): the ORIGINAL logged control did NOT exercise the
failure mode; the remediated controls now do. Residual: script banner and
docstring still say "C^5" for the C^3 control, R6.

**FINDING B (math, load-bearing citation): the certificate's cited bridge,
Zwart Thm 4.16/4.23 [arXiv:2511.16561], has a genuine proof gap at the final
step of Thm 4.23, so the artifacts must repoint to Mathieu's original
Prop 2.2(ii) + Cor 1.3/1.7 (which are sound; §1 step 4).**
Zwart fixes ξ = a highest weight vector of V(τ*) once and for all and infers,
from vanishing of the µ-isotypic component of f^k ⊗ ξ^k, that
(f^k)_{µ+kτ} = 0 ("this must mean (f^k)_λ ⊗ (ξ^k)_{kτ*} = 0, since the sum
is direct", refs/zwart.txt ~line 3310). That inference needs injectivity of
w ↦ proj_µ(w ⊗ ξ^k) on the (µ+kτ)-isotypic component, which is FALSE in
general: for SL(2), λ = ω, σ = ω, µ = 0, proj_triv(e₁ ⊗ e₁) = 0 with
e₁ ≠ 0; the kernel is exactly the highest-weight line. Mathieu's original
proof avoids this precisely by choosing ξ generically (Baire over the
uncountable field, his Prop 2.2 proof, p. 269); Zwart's fixed-ξ
simplification is an unflagged deviation from Mathieu and does not close.
This does NOT affect the truth of "MC(SU(N)) ⇒ homogeneous-form JC_N"
(Mathieu's own proof stands) and hence does not affect our headline; it DOES
mean (i) the certificate docstring's citation "(Mathieu 1997; Zwart
arXiv:2511.16561, Thm 4.16)" must become Mathieu-primary, and (ii) this is
a second erratum-class item in that preprint, independent of the Thm 2.2
statement-level scope issue discussed in NOTE.md §7.

**FINDING C (overclaim): "explicit finite-type (f,g) on SU(79)" is not
delivered; existence is.** No pair (f,g) is written down anywhere in the
artifacts. Worse, PROGRAM.md §4 Step 2's construction sketch (f = the image
of Q_elt under an embedding into C[SU(79)], g = a fixed matrix coefficient of
type µ) is WRONG as written: the isotypic types occurring in Q_elt^k are
(2k+i)ω₁ + iω₇₈, all k-dependent, so for ANY fixed g the pairing
∫ f^k g dk vanishes for all large k automatically; the sketched untwisted
pair satisfies the MC conclusion vacuously and can never witness failure.
The k-independent detectable type exists only after the Cartan twist
L = f ⊗ ξ ∈ A' ∗ A(τ) with a Baire-generic ξ (Mathieu Prop 2.2), which
PROGRAM.md §4 omits entirely and §5's obstruction table sidesteps ("the
counterexample runs the contrapositive and bypasses it"); you cannot both
bypass the machinery and keep the word "explicit". Correct claim: MC(SU(79))
is false, hence a violating K-finite pair EXISTS; an explicit pair would
additionally require carrying the twist ξ (generic, any ξ outside countably
many proper closed sets works, so "explicit modulo one generic parameter" is
achievable but is NEW WORK, not done). R2 de-escalates the claim; the
"endgame deliverable" in PROGRAM.md §7 already points the right way.

**FINDING D (minor, accumulated):**
- D1. PROGRAM.md/certificate cite "Zwart Lemmas 4.18-4.20" for the
  Peter-Weyl dictionary; in the actual v2 text these are Lemmas 4.20, 4.21,
  4.22 + Thm 4.23. Citation drift throughout.
- D2. Certificate docstring justifies "no polynomial inverse" by "(a
  bijective Keller map is an automorphism)", inverted logic; the correct
  one-liner is §1 step 1. (R4)
- D3. Certificate line "(Q_elt^k)_{2k·ω1} = 0 for all k <=> det J L = 1":
  only ⇐ is proven and only ⇐ is needed; drop "<=>".
- D4. PROGRAM.md §0 point 1 still quotes Zwart Thm 2.2 AS STATED ("full JC on
  C^N, the same N") as "Mathieu's theorem", contradicting its own §0.1
  drift-confirmation eight lines later. Align §0 with §0.1 (per-N content =
  homogeneous-form only). Same overshoot in §0 point 2 / §2 table: "SU(2)
  governs planar JC₂"; via Mathieu's actual chain MC(SU(2)) delivers only
  homogeneous-form JC₂, which is classically true, so SU(2) is not wired to
  open JC₂ content by this route at all.
- D5. sweep/STATE references to the certificate as "explicit (f,g)" inherit
  Finding C's language; v1.1 text must not.

## 3. Revisions required for PUBLISHABLE

- **R1.** Rewrite the bridge as §1 steps 1-4: cite Mathieu 1997 Prop 2.2(ii),
  Cor 1.3, Cor 1.7, Lemma 1.2, Props 3.3/3.4, Formula 4.4, with the
  PER-GROUP refinement stated explicitly ("the proofs of Cor 1.3/1.7 consume
  the Main Conjecture only for the maximal compact of the fixed G"); Zwart
  becomes expository-secondary with the Thm 4.23 caveat.
- **R2.** Replace every "explicit finite-type (f,g)" by existence phrasing;
  optionally add the constructive path (Cartan twist, generic ξ) as
  future work.
- **R3.** Commit the fresh rerun log (done: ../mathieu-su79/su79_rerun.log); retire
  su79_run.log or annotate it as superseded (Finding A).
- **R4.** Fix the non-injectivity ⇒ non-polynomial-inverse justification
  (two-sided formal identity argument).
- **R5.** Add the one-line padding argument for N > 79 (§1 step 5) to
  PROGRAM.md §2 and to any v1.1 text claiming "SU(N), N ≥ 79".
- **R6.** Fix the control's "C^5" label to C^3 in banner + docstring;
  keep the C^3 control (it terminates inside the window).
- **R7.** Fix D1 lemma numbers, D3 "<=>", D4 internal inconsistencies.
- **R8.** For a standalone note: state that S^kV* weight-multiplicity-free ⇒
  S^{3k}V ⊗ S^kV* multiplicity free (one line, makes Prop 3.3's use
  self-contained), and remark reducedness/primality of A' per Mathieu
  Lemma 5.4's spectrum analysis (or route through Lemma 1.6 for arbitrary
  G-algebras, which Cor 1.7 already does).

## 4. Focus-point scorecard (ORDERS §B B1)

- (a) integral bridge proved-not-assumed: **DEFECT FOUND** (Findings B, C);
  repaired chain in §1; result stands.
- (b) non-poly inverse ⇒ ∞ many nonzero degrees: **CORRECT** once stated
  (R4 wording fix).
- (c) τ/τ+µ bookkeeping at N = 79: **CORRECT** (§1 step 3).
- (d) negative control exercises true failure mode: **FAILED in the logged
  run** (Finding A); remediated (fresh rerun + controls D/E/F).
- (e) independence: **DONE**: audit/su79_crosscheck.py, from scratch,
  exact arithmetic, disjoint algorithms (plain-Gaussian det vs Bareiss,
  degree-sliced index-tuple engine vs exponent-vector engine, own parser
  cross-validated against sympy), checkpoint agreement term-by-term,
  composition check, three controls. ALL PASS in 41s.

## 5. What I did NOT verify

- druzkowski's symbolic det/nilpotency proof internals (relied on the prior
  post-hoc PASS; consumption spot-checked per doctrine).
- The genuine Haar-integral route (PROGRAM.md §4 route (a), arXiv 2304.02648
  machinery): untested prose, not load-bearing.
- Any statement about N < 79 (N_min bracket [4,79] is [CITED]-grade:
  de Bondt-van den Essen 2005 and Hubbers 1994 were not re-read this pass).
- Zhang/GMC/jacobianfun items: separate jobs.
