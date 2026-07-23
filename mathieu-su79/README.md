# The Mathieu conjecture fails for SU(N), N ≥ 79

**Theorem.** For K = SU(79) there is a finite-type pair (f, g) with
∫ f^n dk = 0 for all n ≥ 1 and ∫ f^n g dk ≠ 0 for infinitely many n.
The same holds for every N ≥ 79 by padding with identity coordinates.
Existence is what is proved; writing an explicit pair requires Mathieu's
Cartan twist by a generic parameter, which is not done here.

**Proof.** `../druzkowski/G_map.txt` gives a cubic-homogeneous Keller map
L = id + H on C^79, lifted from the announced counterexample with its
collision transported (see `../druzkowski/`). Since L is not injective,
its formal inverse is not a polynomial, so the inverse has nonzero
homogeneous components in infinitely many degrees. By Mathieu (1997),
Prop 2.2(ii) with Cor 1.3/1.7, this refutes the conjecture for SU(79).
Notes on the citation chain, including a proof gap in the fixed-ξ step of
Zwart arXiv:2511.16561 Thm 4.23 and its repair via Mathieu's
Baire-generic argument, are in `../audit/su79-review.md`.

**Checks.** `python3 su79_certificate.py` (needs sympy) verifies:
homogeneity of H; det JL = 1 at sample points (symbolic proof in
`../druzkowski/construct_verify.py`); the three-point collision in exact
arithmetic; and nonzero homogeneous components of the formal inverse at
degrees 3–11, with a terminating control map. A second implementation
with different algorithms is `../audit/su79_crosscheck.py`. Reference
log: `su79_rerun.log`.
