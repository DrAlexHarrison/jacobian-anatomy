# Referee report: NOTE.md v0.9 (internal pre-referee pass)

**Verdict: PASS-WITH-EDITS.** No claim in the note is false. Every mathematical
statement I checked is true, every suite reruns green, and every citation I
could verify says what the note says it says. The edits below are: three
places where a [PROVED]/[COMPUTED] label points at a certificate chain that
does not fully cover the stated claim (the claims are true — I supply or
locate the missing certificates), one proof-wording gap, and a handful of
precision edits. Nothing requires structural revision.

- Referee: vetter (internal pre-referee; the external release gate is Alex's
  separate session — this report informs revisions, it does not open the gate).
- Date: 2026-07-20. Note version refereed: v0.9 banner, mtime 02:36.
- My scripts: `audit/vetter_equivariance.py` (35/35), `audit/vetter_wall_identity.py`
  (23/23), `audit/vetter_note_checks.py` (15/15, new this pass) — all exact,
  written independently of `verify/`.

## What I reran and re-verified myself

1. `./verify.sh` — suites 01–05: all PASS (rerun in this session). Suite 06
   (added mid-referee): all PASS.
2. `verify/cas/run-all.sh` — **36 PASS / 0 FAIL tokens** reproduced, including
   all 15 negative controls detecting their breaks. The note's "36/36" and
   negative-control meta-claims are accurate as of my rerun.
3. **`JC.lean` compiles and its axiom claim is TRUE** — I ran
   `lake env lean` on a copy with `#print axioms` appended (the repo had no
   artifact for this; the existing `audit/build.log` covers only Dean's repo).
   Result: exit 0;
   `jacobian_det`, `not_injective`, `jacobian_conjecture_false` each depend on
   exactly `[propext, Classical.choice, Quot.sound]`. Recommend committing
   this output as an artifact.
4. Independent note-checks (`vetter_note_checks.py`): D squarefree over Q
   (gcd with gradient = 1) and irreducible over Q; covering cubic absolutely
   irreducible (linear in c with unit coefficient + primitive in r — valid
   over C, not just Q); y-eliminant as a pure identity in Q[x,y,z] (x = 0
   included); rigidity re-derived by Groebner basis, not solve: the ideal of
   the (exactly ten) nonconstant det-coefficients has GB {t−4, s−3}, so
   (4,3) is the unique solution over C with multiplicity 1; the covering cubic
   itself (not just the x-cubic) at (3/7,−2/5,1/3) is irreducible/Q with
   non-square disc = −606772/18375 = −4D there; the char-2 rescaling G has
   integer coefficients, det ≡ 1, and G(0,1,0) ≡ G(1,1,0) mod 2 (diffs
   (−52,−60,−2)); escape arc, Siphon family, slowform333 identity, deg F = 7,
   origin fiber = {(0,0,0)}.
5. **Pushforward probe** (the one [PROVED] sentence with no certificate):
   Res_r(C_q, x·(br²−4r+3c) − 2r²) = **4c · (D·x³ + (4−3bc)x − 2c)** exactly.
   The §3 "x-eliminant is the pushforward of C_q under the section transform"
   sentence is true and now has a one-line certificate available.
6. Citation sweep (details in the table below).

## (a) The image characterization — SOUND, with one certificate gap

- **Γ-missed half: fully proved, twice over.** (i) Groebner unit certificate
  at Γ(1) + equivariance transport (+ re-certs at Γ(2), Γ(−1/3));
  (ii) the identity route now in `verify/06`: I_x ≡ 0 plus D|_Γ = 0 and
  (4−3bc)|_Γ = 0 forces −2t = 0 at any hypothetical preimage. Route (ii) is
  pure expansion — no Groebner trust needed. Airtight.
- **Everything-else-hit half: sound.** For a ≠ 0, non-perfect-cube ⇒ simple
  root ⇒ the section lifts it (denominators are r·C′(r) ≠ 0 at a simple
  root r ≠ 0; r = 0 is always simple and handled by the x = 0 isomorphism
  onto {c = 0}); perfect-cube locus = Γ exactly (coefficient matching, both
  directions verified). For a = 0: b ≠ 0 gives the explicit {1+xy=0}-sheet
  preimage; b = 0 gives the linear C with a simple root. Complete case split;
  the depressed x-eliminant and its degeneration behave exactly as claimed.
- **R1 (required).** The lifting law's **equality** (and hence the exact
  3/1/0 stratum counts) needs two directions the cited certificates
  (`verify/04`, certs 1–7) do not cover: (i) *at most one* fiber point per
  root — the forcing identity y ≡ y(r)|_{b=F₂,c=F₃,r=x/(1+xy)}, which holds
  identically (my P9, `audit/vetter_wall_identity.py`); (ii) *no* fiber point
  at a multiple affine root (from (i) + the wall identity: such a point would
  force x·0 = r ≠ 0). Note: for generic q the equality can instead be closed
  by "distinct simple roots lift to distinct points (r separates them) plus
  #fiber ≤ generic degree 3", but the exact count 1 on V(D)∖Γ genuinely needs
  (ii). Fix: add the P9 identity + the one-line exclusion to `verify/04`
  (or cite the audit file next to "certificates 1–7").
- **R4 (required, small).** "V(D) is the Jelonek non-properness set" asserts
  both inclusions but only ⊇ is evidenced (fiber-count drop + the escape
  arc). The ⊆ direction (F is proper over C³∖V(D)) is true and needs three
  lines: over U = C³∖V(D) the fiber count is constantly 3 = generic degree;
  take the normalization V̄ of U in C(x,y,z) — finite of degree 3 over the
  smooth U, containing F⁻¹(U) as an open subset; a point of V̄∖F⁻¹(U) over U
  would make some fiber exceed 3; so F⁻¹(U) = V̄ is finite over U, hence
  proper. Add this (it also uses R1's exact counts) or soften to "contains
  and, by the sheet count, equals".

## (b) The primality chain — each arrow checked

Chain rule ⇒ (det JG∘H)·det JH = −2 ⇒ both factors units in C[x,y,z] ⇒
constants ✓. Degrees multiply for dominant quasi-finite compositions ✓.
3 prime ⇒ one factor birational ✓. Birational + étale + normal target ⇒ open
immersion is the correct form of ZMT ✓. Open immersion ⇒ injective ✓.
Injective polynomial self-map ⇒ surjective is Ax–Grothendieck ✓. Bijective
open immersion onto C³ ⇒ automorphism ✓.

- **R3 (required, one clause).** "both determinants are units" elides a step:
  what is constant is det JH and (det JG)∘H. To conclude det JG itself is
  constant one needs H dominant — which holds because H is étale (det JH a
  nonzero constant) ⇒ open image ⇒ dense. Insert "(H is étale, hence
  dominant, so det JG is itself constant)".

## (c) The S₃ monodromy — sound; its input was uncertified

- The argument as written is correct and does NOT need D irreducible:
  degree of C(x,y,z)/C(a,b,c) is 3 (three distinct points in one fiber give
  ≥ 3 via the finiteness bound; the covering cubic gives ≤ 3, and the section
  formulas give C(a,b,c)(r) = C(x,y,z) ✓ — y(r), x(r), z(r) are rational in r
  over the base, and r = x/(1+xy) is rational in the fiber variables). The
  monic normalization of C_q is then the minimal polynomial, hence C_q is
  irreducible over C(a,b,c) (independently: C is linear in c with unit
  c-coefficient and primitive in r, so absolutely irreducible — my N3).
  disc = −4D, and −4D is a nonsquare in C(a,b,c) **iff D is squarefree** (a
  nonconstant squarefree element of the UFD C[a,b,c] has a prime factor to
  odd multiplicity; units of C are squares). Transitive + not in A₃ ⇒ S₃ ✓.
- **R2 (required).** "**[COMPUTED]** D squarefree (gcd with its gradient
  is 1) and irreducible over ℚ" (§4) has **no certificate anywhere in
  `verify/` or `verify/cas/`** — and squarefreeness is the load-bearing input
  to §5. It is true: certified now in `audit/vetter_note_checks.py` (N1, N2).
  Add the two one-liners to a verify script (06 is the natural home) or point
  the label at the audit file.
- **E2 (recommended).** The specialization sentence "a specialized Galois
  group embeds into the generic one": the embedding lands in the *arithmetic*
  group Gal over Q(a,b,c), not directly in the geometric group over C(a,b,c)
  that the main argument computes (the geometric group is the smaller one a
  priori). The main proof stands on its own; reword the confirmation to
  "confirms the arithmetic group is S₃" or drop the embedding clause. Also
  `verify/03` specializes the x-cubic while §5's subject is the covering
  cubic — same field extension, disc differs by the square E², so both are
  valid; my N6 now certifies the covering cubic itself.

## (d) Citations — every checkable one checked

| Citation | Verdict |
|---|---|
| arXiv 2512.23614 for "statement and history" | **Verified.** "On the origin of the Jacobian conjecture" (Rodríguez Díaz, C.R. Math. 364 (2026)) — a statement/history reference (Kraus 1884 vs Keller 1939). Appropriate. |
| Guccione–Guccione–Horruitiner–Valqui, arXiv 2204.14178 | **Verified exactly.** Abstract: all pairs with max(deg P, deg Q) < 125 discarded "except the pair (72,108) (and the symmetric pair (108,72))". The note's frontier sentence is precise. |
| van den Essen, *Polynomial Automorphisms*, **Thm 10.4.2** for DC_n ⇒ JC_n | **Verified** including the numbering (corroborated by the 2021 "New Equivalences" follow-up volume citing Thm 10.4.2 for exactly this implication). |
| Zheglov, arXiv 2410.06959 **v5 (Jan 2026)** | **Verified.** "The Conjecture of Dixmier for the first Weyl algebra is true", v5 submitted 19 Jan 2026, stated unconditionally; "claimed proven, under community review" is the right hedge. |
| Tsuchimoto / Belov-Kanel–Kontsevich JC_{2n} ⇒ DC_n; ladder logic | **Verified** (standard; the note's "JC₄ ⇒ DC₂ leg now vacuous" and "DC₂ ⇒ JC₂, JC₂ ⇒ DC₁ each live" are correct logic given ¬JC₃). DC₁-true and DC_{n≥3}-false are consistent; no tension. |
| Adjamagbo–van den Essen, Acta Math. Vietnam. 32 (2007) | **Consistent with training knowledge** (equivalence of Dixmier/Jacobian/Poisson at the family level); the note's "at the family level" hedge is exactly right. |
| Mathieu 1997; Duistermaat–van der Kallen abelian case | **Consistent** (Mathieu conjecture ⇒ JC; DvdK proved the abelian/Laurent case). |
| Zhao, TAMS 359 (2007) vanishing-conjecture equivalence | **Consistent** (Hessian-nilpotent formulation, JC ⟺ vanishing conjecture). |
| BCW Bull. AMS 7 (1982); Drużkowski Math. Ann. 264 (1983) | **Consistent**; "reductions run backwards" is the correct use (cubic-homogeneous / cubic-linear counterexamples now exist in higher dimension). |
| Jelonek, Ann. Polon. Math. 58 (1993) | **Consistent** (the non-properness set S_F). Cited for the concept; see R4 for the equality's proof obligation. |
| Pinchuk, Math. Z. 217 (1994); "Pinchuk-type maps generate degree-6 field extensions" | **Verified** — Campbell (arXiv 1208.5177, 1202.2949): Pinchuk maps are generically 2-to-1 as real maps yet induce **degree-6** function-field extensions, trivial automorphism group. The note's degree-6 parenthetical is accurate. |
| Shestakov–Umirbaev 2004 preemption | **Correct and valuable** — they proved wildness of the Nagata automorphism; no "2004 counterexample to JC" exists. |
| deancureton/jacobian | **Verified live**: repo exists, Lean 4, theorem names exactly `not_jacobianConjecture{,_all_char,_complex}`, README credits "Levent Alpöge/Fable 5's counterexample". The note's audit summary matches `audit/deancureton-audit.md`, whose build.log I inspected (8661 jobs, exit 0, three standard axioms for every theorem). |
| Paul-Lez, formal-conjectures PR #4474 | **Verified live**: "feat: add Jacobian disproof", authored by Paul-Lez, open. |
| Wang (degree 2 impossible), §8 | Correct (Wang 1980); currently uncited — optional cite since §8 is a questions section. |

- **E3 (recommended).** §7 real form: "**strictly stronger** than Pinchuk's
  counterexample" compares across dimensions — Pinchuk's lives in the harder
  n = 2 and is not dominated by an n = 3 example. The parenthetical already
  states the differences honestly; change "strictly stronger than" to e.g.
  "with a strictly stronger Jacobian condition than" or "a constant-Jacobian
  counterpart, in n = 3, to". The degree-6 and constant-vs-nonconstant facts
  themselves are verified.

## (e) The étale-properness remark

**PRESENT AND CORRECT in both occurrences** — §4 ("a **connected** finite
étale cover of the simply connected C³ would be trivial") and §7 ("its source
C^n is connected and its target C^n is simply connected, so the covering has
degree 1"). The four-line repaired-conjecture theorem is sound as written.

## New-content checks (since the equivariance pass)

- **Three presentations (§3):** v = y + 1/x = (1+xy)/x = 1/r verified; the
  reversed-cubic identity verified; the pushforward sentence now certified by
  Res_r(C, x(br²−4r+3c) − 2r²) = 4c·(x-cubic) — **E5 (recommended):** add
  that one-liner to verify/ so the sentence's [PROVED] label is self-contained.
- **Rigidity (§3):** verified independently at Groebner level (stronger than
  the suite's solve-based check): exactly ten nonconstant coefficients, ideal
  = ⟨t−4, s−3⟩, so (4,3) is unique over C with multiplicity 1, det = −2 there.
  The note's "ten nonconstant coefficients" is exactly right. Wording "isolated
  point of its own family" is scoped to the ansatz — fine as written.
- **Zheglov hedge (§7):** verified and properly hedged (see table).
- **Dual-CAS (provenance section):** reran; 36/36 with all negative controls
  detected; "different algorithms" claim is substantiated (Groebner vs
  iterated resultants; poldisc vs resultant identity; vdim vs numeric roots).
- **Lean paragraph (§1):** Dean-repo claims match the audit file and build
  log; JC.lean's own claims now verified by my compile+axioms run (item 3
  above). "First public ... as far as we can tell" is properly hedged. The
  char-2 trick arithmetic verified (N7). The honest caveat (char p > 0
  classically false) is present and correct.

## Required edits (before external referee)

1. **R1** — lifting-law uniqueness: add the y ≡ y(r) forcing identity and the
   multiple-root exclusion to `verify/04` (both available in
   `audit/vetter_wall_identity.py`, check P9), or cite that file where the
   note says "certificates 1–7".
2. **R2** — certify "D squarefree + irreducible over Q" inside `verify/`
   (one-liners; currently only in `audit/vetter_note_checks.py` N1/N2).
3. **R3** — §6: insert the H-dominance clause so det JG (not just det JG∘H)
   is constant.
4. **R4** — §4 non-properness: add the three-line properness-off-V(D)
   argument (or soften "is" to an evidenced form).
5. **Stale line** — "Reproducing everything" says "all sympy suites (01-05)";
   verify.sh now runs 01-06.

## Recommended edits

- **E1** — §4 stratification bullets: the parentheticals ("3 simple roots",
  "double + simple") silently assume a ≠ 0; the counts are correct on a = 0
  too, but add a pointer to the a = 0 clause of the lifting law.
- **E2** — §5 specialization wording (arithmetic vs geometric group), and
  optionally cite N6 (covering cubic specialization) beside verify/03's
  x-cubic one.
- **E3** — §7 "strictly stronger than Pinchuk" rewording.
- **E4** — §4 "as an ideal, Γ = V(bc−4/3, 12a−b²)": mixed ideal/variety
  language. The intended stronger claim is true — the ideal is prime (the
  quotient is ≅ C[c, c⁻¹]) hence radical, so I(Γ) = ⟨bc−4/3, 12a−b²⟩ — say
  "with prime defining ideal ⟨bc−4/3, 12a−b²⟩".
- **E5** — add the pushforward resultant certificate (4c constant) to verify/.
- **E6** — commit the JC.lean axiom-audit output as an artifact (currently
  only Dean's repo has a committed build log).

## What I tried to break and could not

Weight forcing; Γ closedness and orbit structure; both halves of the image
claim (including three Γ-points and 14 stratum-spanning fiber counts);
equality of disc classes between the r- and x-presentations; the rigidity
uniqueness (via Groebner, checking sympy's solver didn't silently drop
solutions); the char-2 collision arithmetic; the negative controls' ability
to actually detect breaks (reran, all detected); every arXiv/book citation I
could reach, including the two easiest to have gotten wrong — the (72,108)
frontier and the Pinchuk degree-6 parenthetical — both of which the note
states exactly correctly.
