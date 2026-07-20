# jacobian-anatomy

Independent verification and complete geometric anatomy of the counterexample
to the Jacobian Conjecture (the map of Alpöge–Mathew–Fable, announced
2026-07-19).

- **`NOTE.md`** — the anatomy note: C\*-equivariance, the covering cubic and
  lifting law, monodromy S₃, non-properness quartic, exact image C³ ∖ Γ,
  rigidity, and consequences (Dixmier, Mathieu, Zhao, real forms).
- **`verify.sh`** — one command reproduces every computational claim
  (`./verify.sh`; needs python3 + sympy).
- **`JC.lean`** — a Lean 4 / mathlib certificate (det ≡ −2 and non-injectivity
  over ℚ; axioms: `propext`, `Classical.choice`, `Quot.sound` only). This is an
  independent *second* formalization; the first public one, as far as we can
  tell, is [deancureton/jacobian](https://github.com/deancureton/jacobian),
  whose build-and-axiom audit by us is in **`audit/`** (verdict: clean).

**Who made this.** Researched, computed, and written by Claude (Fable 5),
directed by Alex Harrison (research physiologist and software developer). The
mathematics is machine-generated and machine-verified; the standards —
everything reproducible, everything adversarially checked, errors fixed
prominently — are the human contribution. Run `./verify.sh`; trust the
certificates, not us.

**Errata policy.** Report any error against a labeled claim (GitHub issue,
or email harrison.alexander.p@gmail.com); confirmed errors are acknowledged
prominently at the top of `NOTE.md`.

**Citing.** Please cite the tagged release (a DOI via Zenodo is planned at
v1.0).

**License.** Code (verify/, JC.lean): MIT. Prose (NOTE.md, README, audits):
CC BY 4.0. See `LICENSE`.
