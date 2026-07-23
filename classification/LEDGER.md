# Per-class ledger

Totals recomputed from the raw outputs (results/fNNNN.out verdict lines,
the hardfamily/ logs, and the decomposed-route outputs), not from any
summary. The recomputation output is hardfamily/final_ledger_data.txt.

## Headline count

**371 isomorphism classes of MIXED equivariant families (741 enumerated
systems = 370 tied-permutation twin pairs + the self-paired f0740).
371/371 covered: 365 classes by DIRECT machine certificate, 6 classes by
hand theorems, 0 uncovered. Zero NONUNIT appeared in any
collision certificate of any system across the entire operation.**

Structural classes (T3 lemmas, unchanged): PROPER 113, PLANAR-SPLIT 6,
AXIS 1, WANG 2, LINEAR 1 — INJECTIVE by the structural lemmas; no duplicates.

## Per-class grades

1. **DIRECT (365 classes; 730 members hold certificates).** At least one
   member (in all but a handful, both members) carries a machine EMPTY
   certificate: monolith triple-Rabinowitsch Gröbner unit certificates
   (657 members incl. the 7200s retry cracks f0632/f0641/f0648/f0671) or
   decomposed-route per-component certificates (77 members: minAssGTZ
   a-route and facstd b-route, with per-component DEAD/LIVE evidence and
   3x UNIT lines; aNNNN.out / bNNNN.out on disk). Twin pairs where both
   members finished agree EMPTY in every single case (30+ silent dual-run
   confirmations, supporting evidence only).

2. **HAND THEOREMS (6 classes):**
   - {720,723}, {721,722}, {728,729} — Theorem 26.1 (affine GL(2)
     stratum: gauge lemma, residue-matching kill, slice-Moh endgame;
     ladder extractions machine-verified 12/12; S1-branch certificates
     ALL UNIT at (3,1) and (3,2); (4,1) S1 rests on the uniform hand
     descent; see PROOFS.md §2).
   - {738,739} — Theorem 10.1 (zero-weight stratum; d = 3 case; per-
     coefficient forcing certificates 6/6 + 6/6 UNIT at the generic
     shape level; per-system z-certificates timed out, not load-bearing).
   - {740} — Theorem 4.2 = fable (7.5): two independent blind
     derivations; Moh-cited planar reduction.
   - {734,735} — Theorem 35.1 (the {734,735} proof, notebook Entries 32-35, commit
     2883a57): chart bijection CAL-A, Normalization Lemma machine-
     grounded 9/9, normalized-ladder hand tree with N1/N2 enslavements
     and N8/N7 dichotomies, all leaf certificates UNIT including a
     monolithic gauge-free confirmation. **f0735 member grade:
     empty by transfer** (sigma = (12) witness, ground-truth job-file
     supports, negative controls, per the blessed 00:10 format) — the
     only transfer-graded member in the ledger.

3. **UNCOVERED: none.**

## Honesty items (enumerated per the 04:38 order; every one carries into any write-up)

H1. **Deg-7 control, zero-witness route:** the 22-free gate rung passed by
    the EXPLICIT-ZERO route (Z5E: facstd decomposition + point-in-LIVE-
    component + all generators of the point component's three collision
    systems vanishing at the announced fiber pairs, dual-evaluated sympy +
    Singular). The Gröbner-unit route passed at ≤ 16-free only (both
    minAssGTZ and facstd code paths). State both routes distinctly.
    Gate SATISFIED per coordination 20:36 Jul 21, subject to review's
    re-run of the rung (condition 6).
H2. **Full-26 deg-7 control: UNDECIDED, permanently logged as such.**
    1800s, 7200s (both routes), and 86400s runs all died inside
    decomposition (ctrl3.log FULL26-24H rc=124 dt=86401s). The a-priori
    dead-discard note (det0 at the Alpöge vector = +2 ≠ 0) plus CAL-4
    close the false-negative channels except termination itself.
H3. **RESOURCE-UNDECIDED lines:** stragglers round-1 720/721 entries were
    a memory famine (~23:33 Jul 20), annotated in-log, not genuine
    attempts; the affected lines are annotated RESOURCE-UNDECIDED.
H4. **enumerate.py dedupe weaker than documented:** counts 864/741/88 are
    duplicate-inclusive; the fold is 371 classes (iso_audit.py, exact
    S-triple witnesses). SAFE direction (under-merging); review
    verifies the merged-side sample per the 00:10 ruling item 5b. The
    1750s easy/hard split is partly coordinate-artifact (three twin pairs
    straddle it).
H5. **Citation dependencies:** Moh 1983 (planar JC through degree 100)
    carries the slice/planar endgames of T3-PLANAR/AXIS, Thm 4.2,
    Thm 10.1 (d=2), Thm 26.1, Thm 35.1; Białynicki-Birula–Rosenlicht 1962
    carries T3-PROPER; Wang 1980 carries WANG. All [CITED], classical,
    not re-proved here.
H6. **Independence framing:** Entries 3-4 (Claim 3.1, Thm 4.2) are blind
    two-minds vs fable-hardfamily. Entry 10 (Thm 10.1's d ∈ {3,4}) is an
    INFORMED COMPLETION (author had read fable §7), backed by machine
    certificates — not a second blind mind. Thm 26.1 and Thm 35.1 are
    single-author-plus-machine.
H7. **the beta-coincidence sub-leaf:** one sub-leaf of the Entry-34
    tree is closed by CERTIFICATE ONLY (no hand argument); flagged by
    the {734,735} proof, must be stated in the write-up and walked by review.
H8. **Timeouts are UNDECIDED** everywhere and always were; no timeout
    contributed to any EMPTY grade. facstd branch counts are not
    isomorphism invariants (order-sensitivity; verdicts, not
    decompositions, are the comparable objects).

## What this ledger does NOT claim

No "complete classification" statement is made here. That claim remains
FORBIDDEN until coordination's declaration on this filing and review's
route walk (Thm 26.1 + Thm 35.1 + transfer lemma + the H-items) complete.
This document is the drain accounting, nothing more.

Filed 2026-07-22 04:38 MDT (stamped from date). — the ledger pass
