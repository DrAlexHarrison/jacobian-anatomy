# The Mathieu conjecture is false for SU(79) — runnable certificate

**Theorem (refereed; report: `../audit/referee2-su79.md`).** The Mathieu
conjecture fails for K = SU(79), and for SU(N) for every N >= 79: a
finite-type pair (f, g) on SU(79) exists with ∫ f^n dk = 0 for all n >= 1
yet ∫ f^n g dk != 0 for infinitely many n. (Existence is what is certified;
an explicit pair additionally requires Mathieu's Cartan twist by a generic
parameter — open work.) For N > 79, pad L by identity coordinates: the
padded map is still cubic-homogeneous, Keller, and non-injective, and the
whole chain reruns verbatim (report §1, step 5).

**How it is proved.** `../druzkowski/G_map.txt` is the explicit
cubic-homogeneous Keller map L = id + H on C^79 lifted from the announced
counterexample with its collision transported (see `../druzkowski/`).
`su79_certificate.py` machine-checks, in exact rational arithmetic:
homogeneity (CHECK 1), det J L = 1 spot-checks (CHECK 2; full symbolic proof
in `../druzkowski/construct_verify.py`), the explicit 3-point collision
(CHECK 3), and nonzero homogeneous components of the formal inverse at
degrees 3-11 with a genuinely terminating negative control. Non-injectivity
forces the formal inverse to be non-polynomial, which through Mathieu's
original 1997 argument (Prop 2.2(ii) + Cor 1.3/1.7; citation chain audited
in the referee report, including a documented proof gap in the fixed-xi step
of Zwart arXiv:2511.16561 Thm 4.23 and its sound repair) refutes the Mathieu
conjecture for SU(79).

**Run it:** `python3 su79_certificate.py` (pure Python + sympy for parsing;
exact fractions throughout; ~minutes; writes a checkpoint json). The
authoritative reference log is `su79_rerun_referee3.log`. An INDEPENDENT
from-scratch verifier with disjoint algorithms is
`../audit/referee2_su79_verify.py` (term-by-term agreement with the
certificate's formal inverse; all checks pass in ~41s).
