# jacobian-anatomy

Machine-verified anatomy of the Alpöge–Mathew–Fable counterexample to the
Jacobian Conjecture (announced 2026-07-19), and the classification proving
its degree is minimal within its symmetry class.

- `NOTE.md`: the anatomy. C\*-equivariance, the covering cubic and lifting
  law, S₃ monodromy, non-properness quartic, exact image C³ ∖ Γ, rigidity,
  consequences (Dixmier, Mathieu, Zhao, real forms).
- `classification/` (v1.2): no C\*-equivariant Keller counterexample of
  total degree ≤ 6 exists in C³, so the announced degree-7 map is minimal
  in its symmetry class. 371 isomorphism classes; 741 rerunnable
  certificate jobs; four hand theorems with full proofs; negative
  controls. Scope and qualifiers: `classification/CLASSIFICATION.md`.
  Quick independent check with sympy alone: `python3 classification/check.py`.
- `verify.sh`: reproduces every computational claim in NOTE.md
  (needs python3 + sympy).
- `JC.lean`: Lean 4 / mathlib certificate (det ≡ −2 and non-injectivity
  over ℚ; axioms: `propext`, `Classical.choice`, `Quot.sound`). Second
  public formalization; the first is
  [deancureton/jacobian](https://github.com/deancureton/jacobian),
  audited clean in `audit/`.

**Authorship.** Researched, computed, and written by Claude (Fable 5),
directed by Alex Harrison. All computations are reproducible from this
repository (`./verify.sh`, `classification/REPRODUCE.md`).

**Errata.** Report errors by GitHub issue or
harrison.alexander.p@gmail.com. Confirmed errors are acknowledged at the
top of `NOTE.md`.

**Citing.** Cite the tagged release. DOI (all versions):
[10.5281/zenodo.21465361](https://doi.org/10.5281/zenodo.21465361).

**License.** Code (verify/, JC.lean, classification/ scripts and jobs):
MIT. Prose (NOTE.md, README, audits, classification/ documents): CC BY 4.0.
See `LICENSE`.
