# Reproducing the classification

Requirements: Singular 4.x, Python 3 with sympy. Optional: msolve.

No Singular installed? `python3 check.py` re-proves a sample of the
certificates in sympy's own Groebner engine and recomputes the class
fold, with no dependency beyond sympy.

## Any single certificate
    Singular -q < jobs/f0001.sing
Compare with results/f0001.out; a job prints its verdict line when all
three collision certificates are Gröbner units. The 88 Gröbner-hard
systems (LEDGER.md; 4 of them later completed monolithically at a 7200 s
budget) otherwise time out monolithically; their decomposed
certificates are in hardfamily/ (aNNNN: minAssGTZ route, bNNNN:
factorized-Buchberger route), runnable the same way.
scripts/verdicts.json records the sweep verdicts; LEDGER.md the
per-class totals.

## The enumeration
    python3 scripts/enumerate.py
    python3 scripts/class_fold.py
    python3 scripts/iso_audit.py

## The computations cited in PROOFS.md
    cd scripts
    python3 q8_cal_chart.py
    python3 q8_layers.py
    python3 q12_gauge.py
    python3 q14_handsteps.py
    python3 q13_transfer_735.py
    python3 q10_leafchecks.py
    python3 q11_beta.py
    python3 q1_identity_check.py
    python3 q1_forcing_check.py
    python3 q4_affine_ladder_check.py

## Negative control
Outputs in hardfamily/: adeg7z*.out and ctrl*.log; the explicit-zero
run at 22 free coefficients is adeg7z5e.out. CLASSIFICATION.md states
what the control does and does not establish.
