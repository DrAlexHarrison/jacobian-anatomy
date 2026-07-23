# The equivariant classification at degree ≤ 6

**Theorem.** Every polynomial map F: C³ → C³ that (i) is equivariant for
a nontrivial C*-action with matched source and target weights (the
target weight vector a permutation of the source weights, the canonical
form established by the enumeration), (ii) has components of total
degree ≤ 6, and (iii) has constant nonzero Jacobian determinant, is
injective. Equivalently, no C*-equivariant Keller counterexample of
total degree ≤ 6 exists in C³.

**Consequence.** The Alpöge counterexample (degree 7, weights (1,−1,−2))
attains the minimal possible degree within its symmetry class.

**Scope.**
- Equivariant maps only. Whether a non-equivariant Keller counterexample
  of degree ≤ 6 exists in C³ is open and not addressed here.
- Over C.
- "Degree ≤ 6" means total degree of each component within the enumerated
  supports; the enumeration of weight systems and supports is part of the
  artifact (scripts/enumerate.py, scripts/families_deg6.json).

## Evidence

The enumeration yields 864 families. 123 fall to structural lemmas
(PROPER 113, PLANAR-SPLIT 6, AXIS 1, WANG 2, LINEAR 1), each class
injective. The remaining 741 systems form 371 isomorphism classes under
simultaneous source/target coordinate permutation (370 twin pairs and
one self-paired system; the fold, from the job files' own supports:
scripts/class_fold.py; permutation witnesses and transfer checks for the
Gröbner-hard systems: scripts/iso_audit.py, scripts/iso_transfer.py).

Each of the 371 classes contains no non-injective Keller member:

1. 365 classes by Gröbner unit certificates: triple-Rabinowitsch
   collision systems (jobs/fNNNN.sing → results/fNNNN.out), or
   decomposed per-component certificates for the Gröbner-hard systems
   (hardfamily/aNNNN.out, bNNNN.out). No collision certificate in any
   system returned NONUNIT.
2. 5 classes by the hand theorems of PROOFS.md §1–§2 (the zero-weight
   and affine GL(2) strata), with forcing certificates.
3. The class {734,735} by PROOFS.md §3–§4. One sub-case there rests on
   Gröbner certificates rather than a written chain (see Limitations);
   the second member of the class follows by the transfer lemma
   (scripts/q13_transfer_735.py).

## Negative control

The pipeline was tested against the known degree-7 counterexample
family, which it must never certify empty. With coefficients frozen at
the counterexample's values and k of the 26 coefficients freed:

- k = 7, 10, 16: both decomposition routes locate the counterexample in
  a live component and its three collision systems are non-unit
  (at k = 16 the minAssGTZ route on the point component).
- k = 19, 22: the decomposition completes and locates the counterexample
  in a live component; the Gröbner half of the collision check exceeded
  its time budget. At k = 22 the collision systems were instead shown
  non-unit directly: the known collision fiber is an explicit common
  zero of every generator, evaluated independently by two systems.
- k = 24 and k = 26 (the full family): the decomposition itself did not
  terminate within budget (24 h at k = 26) and these rungs are
  undecided.

Two a-priori facts supplement the ladder: the component-discard rule
cannot drop the witness component (det0 at the counterexample's
coefficients is nonzero), and the production ideal has an explicit zero,
so a unit certificate there is impossible. No full-generality control is
claimed.

## Dependencies

- Moh (1983): planar Keller maps of degree ≤ 100 are injective. Endgame
  of every hand theorem.
- Ax–Grothendieck and Białynicki-Birula–Rosenlicht: an injective
  polynomial self-map of Cⁿ is an automorphism with polynomial inverse.
  Used in PROOFS.md §3.
- A Gröbner unit certificate (1 ∈ ideal) is a Nullstellensatz proof of
  emptiness. The computations were carried out in Singular 4.x;
  cross-checks on sampled systems in msolve; reruns under reversed
  variable orders; the two routes (monolithic and decomposed) agree on
  every system both decided, and the 30+ twin pairs decided
  independently on both sides agree in every case.

## Limitations

- The control ladder is as stated above: collision non-unitness is
  Gröbner-certified through 16 free coefficients, shown by explicit zero
  at 22, and undecided at 24 and 26.
- 84 systems never received a monolithic sweep verdict at any budget;
  their classes are covered by the decomposed certificates, the twin
  certificates with the transfer lemma, or the hand theorems. Timeouts
  were never counted as verdicts in either direction
  (scripts/verdicts.json records the sweep; LEDGER.md the per-class
  totals).
- The enumeration's permutation dedupe retained twin pairs it was meant
  to merge. This produces redundant work, not gaps; counts here are
  stated at class level for this reason.
- In PROOFS.md §3, the sub-case with the three quadratic-level blocks
  nonzero and coincident roots is decided by Gröbner certificates; the
  written chain stops at a relation system.
- Two early runs failed from machine memory pressure and were re-run
  cleanly; the affected log lines in hardfamily/ are annotated
  RESOURCE-UNDECIDED.
- Three long-budget certificates for PROOFS.md §3 branches already
  closed by other routes were left running: one finished (unit), two
  exceeded a 10800 s budget and are undecided. Nothing rests on them.

## Contents

- CLASSIFICATION.md: this file.
- PROOFS.md: the hand theorems.
- LEDGER.md: per-class totals, recomputed from the raw outputs.
- REPRODUCE.md: how to re-run any certificate or the whole sweep.
- jobs/ (741 Singular jobs), results/ (741 outputs), hardfamily/
  (decomposed-route certificates, control-ladder outputs, and logs for
  the 88 Gröbner-hard systems, 4 of which later completed monolithically), scripts/ (enumeration, generation,
  isomorphism audit, and the computation scripts cited in PROOFS.md).
