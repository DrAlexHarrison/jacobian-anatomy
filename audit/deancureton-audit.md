# Audit: deancureton/jacobian (Lean 4 formalization of Alpöge's Jacobian Conjecture counterexample)

Auditor: Claude (Fable 5), Alex Harrison's workstation. Audit date: 2026-07-20.
Audited commit: `0d4a9212d874226ad81ce5a926becddfa94e6a88` (HEAD of `main`).
Local clone: `/home/alex/code/jc/audit/jacobian`

## 1. Repo metadata

- URL: https://github.com/deancureton/jacobian (clone succeeded on first try; no 404).
- GitHub `created_at`: 2026-07-20T04:30:27Z. `pushed_at`: 2026-07-20T05:06:46Z. 7 stars, 1 fork at audit time.
- Author: Dean Cureton (GitHub uid 89793485), noreply email; no GPG data checked.
- Commit history (`git log --format=fuller`, all times -0400 EDT):

```
0d4a9212 readme                AuthorDate/CommitDate: Mon Jul 20 01:06:19 2026 -0400
7197fc2d all characteristics   AuthorDate/CommitDate: Mon Jul 20 01:00:17 2026 -0400
66b4ec86 init                  AuthorDate/CommitDate: Mon Jul 20 00:28:38 2026 -0400
```

  In UTC: init 04:28:38Z, all-characteristics 05:00:17Z, readme 05:06:19Z, all in a 38-minute window early 2026-07-20, i.e., hours after Alpöge's 2026-07-19 announcement.
- File inventory (complete): `Jacobian/Counterexample.lean` (162 lines, the entire mathematical content), `Jacobian.lean` (one-line root import), `lakefile.toml`, `lake-manifest.json`, `lean-toolchain`, `README.md`, `.gitignore` (`/.lake`). No CI workflow, no other Lean files.
- Toolchain: `leanprover/lean4:v4.33.0-rc1`. Mathlib pin: tag `v4.33.0-rc1`, rev `79d0395a1825a6264ad5d269e35e60537518955e` (differs from our local mathlib pin `c732b96d05ef`, so their cache was fetched fresh via `lake exe cache get`).
- README title credits "Levent Alpöge/Fable 5's counterexample": the author publicly attributes the formalization assistance to a Claude Fable 5 model.

## 2. Build result: SUCCESS

Executed on this workstation 2026-07-20 (elan 4.2.3, toolchain auto-installed from their `lean-toolchain`):

```
cd /home/alex/code/jc/audit/jacobian
lake exe cache get     # 08:14:14Z -> 08:18:02Z, 8643/8643 cache files downloaded
lake build             # 08:18:02Z -> 08:18:27Z
```

Tail of build log (`/home/alex/code/jc/audit/build.log`), verbatim:

```
⚠ [8659/8661] Built Jacobian.Counterexample (13s)
warning: Jacobian/Counterexample.lean:1:1: * '-/':
Copyright too short!
Note: This linter can be disabled with `set_option linter.style.header false`
info: Jacobian/Counterexample.lean:159:0: 'not_jacobianConjecture' depends on axioms: [propext, Classical.choice, Quot.sound]
info: Jacobian/Counterexample.lean:160:0: 'not_jacobianConjecture_all_char' depends on axioms: [propext, Classical.choice, Quot.sound]
info: Jacobian/Counterexample.lean:161:0: 'not_jacobianConjecture_complex' depends on axioms: [propext, Classical.choice, Quot.sound]
✔ [8660/8661] Built Jacobian (5.5s)
Build completed successfully (8661 jobs).
=== build exit: 0 at 2026-07-20T08:18:27Z ===
```

**`lake build` exits 0 (8661 jobs).** The single warning is a mathlib style linter complaining the file lacks a copyright header: purely cosmetic, zero mathematical content. `Jacobian.Counterexample` itself elaborated in 13s.

## 3. Sorry/axiom sweep

`grep -rn -E 'sorry|admit|native_decide|axiom|unsafe|implemented_by|partial def|extern|opaque'` over all `.lean` sources:

- Only hits: lines 159-161 of `Counterexample.lean`, which are the author's own `#print axioms` commands for the three headline theorems. Zero occurrences of `sorry`, `admit`, `native_decide`, `axiom` declarations, `unsafe`, `@[implemented_by]`, `partial def`, `extern`, or `opaque`.
- No custom macros, no custom tactics, no `set_option` in the source file. `lakefile.toml` leanOptions are benign (`pp.unicode.fun`, `relaxedAutoImplicit false`, mathlib linter set, `maxSynthPendingDepth 3`).

## 4. #print axioms (gold standard): ALL CLEAN

Independent check: auditor-written scratch file `AxiomsAudit.lean` (local only, in the clone but never committed) importing `Jacobian.Counterexample`, run via `lake env lean AxiomsAudit.lean` against their completed build. Output verbatim:

```
@not_jacobianConjecture : ∀ {K : Type u_1} [inst : Field K],
  2 ≠ 0 →
    ¬∀ (F : Fin 3 → MvPolynomial (Fin 3) K),
        IsUnit (MvPolynomial.jacobianDet F) → Function.Injective (MvPolynomial.evalMap F)
not_jacobianConjecture_all_char : ∀ (K : Type u_1) [inst : Field K],
  ¬∀ (F : Fin 3 → MvPolynomial (Fin 3) K),
      IsUnit (MvPolynomial.jacobianDet F) → Function.Injective (MvPolynomial.evalMap F)
not_jacobianConjecture_complex : ¬∀ (F : Fin 3 → MvPolynomial (Fin 3) ℂ),
    IsUnit (MvPolynomial.jacobianDet F) → Function.Injective (MvPolynomial.evalMap F)
@MvPolynomial.jacobianDet : {R : Type u_1} →
  {σ : Type u_2} → [inst : CommRing R] → [Fintype σ] → [DecidableEq σ] → (σ → MvPolynomial σ R) → MvPolynomial σ R
@MvPolynomial.evalMap : {R : Type u_1} →
  {σ : Type u_2} → [inst : CommSemiring R] → (σ → MvPolynomial σ R) → (σ → R) → σ → R
'MvPolynomial.jacobianDet' depends on axioms: [propext, Classical.choice, Quot.sound]
'MvPolynomial.evalMap' depends on axioms: [propext, Classical.choice, Quot.sound]
'Jacobian.jacobianDet_F' depends on axioms: [propext, Classical.choice, Quot.sound]
'Jacobian.jacobianDet_G' depends on axioms: [propext, Classical.choice, Quot.sound]
'Jacobian.evalMap_F_p0' depends on axioms: [propext, Classical.choice, Quot.sound]
'Jacobian.evalMap_F_p1' depends on axioms: [propext, Classical.choice, Quot.sound]
'Jacobian.evalMap_F_p2' depends on axioms: [propext, Classical.choice, Quot.sound]
'Jacobian.evalMap_G_p0' depends on axioms: [propext, Classical.choice, Quot.sound]
'Jacobian.evalMap_G_p1' depends on axioms: [propext, Classical.choice, Quot.sound]
'Jacobian.evalMap_G_p2' depends on axioms: [propext, Classical.choice, Quot.sound]
'Jacobian.evalMap_G_char_two' depends on axioms: [propext, Classical.choice, Quot.sound]
'not_jacobianConjecture' depends on axioms: [propext, Classical.choice, Quot.sound]
'not_jacobianConjecture_all_char' depends on axioms: [propext, Classical.choice, Quot.sound]
'not_jacobianConjecture_complex' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Every definition, lemma, and headline theorem reports exactly `[propext, Classical.choice, Quot.sound]`, the three standard axioms of classical mathlib. No `sorryAx`, no `Lean.ofReduceBool`/`ofReduceNat` (i.e., no `native_decide`), no custom axioms anywhere in the dependency cone. The elaborated `#check` statements also confirm the theorems say exactly what the source reads: no hidden hypotheses beyond `[Field K]` (plus `2 ≠ 0` only on the det=−2 variant, as intended).

## 5. Faithfulness audit

### 5.1 Definitions, quoted verbatim from Counterexample.lean

```lean
def jacobianMatrix [CommSemiring R] [DecidableEq σ] (F : σ → MvPolynomial σ R) :
    Matrix σ σ (MvPolynomial σ R) :=
  Matrix.of fun i j => pderiv j (F i)

def jacobianDet [CommRing R] [Fintype σ] [DecidableEq σ] (F : σ → MvPolynomial σ R) :
    MvPolynomial σ R :=
  (jacobianMatrix F).det

def evalMap [CommSemiring R] (F : σ → MvPolynomial σ R) (p : σ → R) : σ → R :=
  fun i => eval p (F i)
```

- `pderiv` is mathlib's formal partial derivative on `MvPolynomial` (a `Derivation`); `Matrix.det` is mathlib's determinant. The Jacobian determinant is computed: `jacobianDet_F : jacobianDet (F K) = C (-2)` is *proved* by `simp` (expanding `det_fin_three` and the `pderiv` rules) followed by `ring`. Same for `jacobianDet_G : jacobianDet (G K) = 1`.
- `evalMap` is the genuine evaluation self-map of K³. The instance arguments (`DecidableEq (Fin 3)`, `Fintype (Fin 3)`) are canonical and carry no hidden content.

### 5.2 The map F is exactly Alpöge's map

Their definition (verbatim, X 0 = x, X 1 = y, X 2 = z):

```lean
def F : Fin 3 → MvPolynomial (Fin 3) K :=
  ![(1 + X 0 * X 1) ^ 3 * X 2 + X 1 ^ 2 * (1 + X 0 * X 1) * (C 4 + C 3 * (X 0 * X 1)),
    X 1 + C 3 * X 0 * (1 + X 0 * X 1) ^ 2 * X 2 + C 3 * X 0 * X 1 ^ 2 * (C 4 + C 3 * (X 0 * X 1)),
    C 2 * X 0 - C 3 * X 0 ^ 2 * X 1 - X 0 ^ 3 * X 2]
```

Character-by-character this is Alpöge's announced map:
F₁ = (1+xy)³z + y²(1+xy)(4+3xy); F₂ = y + 3x(1+xy)²z + 3xy²(4+3xy); F₃ = 2x − 3x²y − x³z.
Independently re-verified with sympy on this machine (script: scratchpad `faithfulness_check.py`):

```
F identical to Alpöge: True True True
det J(F) = -2
F (0, 0, -1/4)  = [-1/4, 0, 0]
F (1, -3/2, 13/2)  = [-1/4, 0, 0]
F (-1, 3/2, 13/2)  = [-1/4, 0, 0]
```

Witness points in Lean match Alpöge's exactly: `![0, 0, -(1/4)]`, `![1, -(3/2), 13/2]`, `![-1, 3/2, 13/2]`, all proven to map to `![-(1/4), 0, 0]` (the char≠2 theorems carry the necessary `(2 : K) ≠ 0` hypothesis so the 1/2- and 1/4-denominators make sense; distinctness of the two collision witnesses is extracted from first coordinates 1 ≠ −1, which is exactly char ≠ 2).

### 5.3 The det-1 form G (their original contribution beyond Alpöge's post)

```lean
def G : Fin 3 → MvPolynomial (Fin 3) K :=
  ![(1 + C 2 * X 0 * X 1) ^ 3 * X 2
      + C 4 * X 1 ^ 2 * (1 + C 2 * X 0 * X 1) * (C 2 + C 3 * (X 0 * X 1)),
    X 1 + C 3 * X 0 * (1 + C 2 * X 0 * X 1) ^ 2 * X 2
      + C 12 * X 0 * X 1 ^ 2 * (C 2 + C 3 * (X 0 * X 1)),
    -X 0 + C 3 * X 0 ^ 2 * X 1 + X 0 ^ 3 * X 2]
```

Claimed provenance: G = diag(1/2, 1/2, −1/2) ∘ F ∘ diag(1, 2, 2). Sympy-verified on this machine:

```
G matches diag(1/2,1/2,-1/2) o F o diag(1,2,2): True True True
det J(G) = 1
G (0, 0, -1/8)   = [-1/8, 0, 0]
G (1, -3/4, 13/4)  = [-1/8, 0, 0]
G (-1, 3/4, 13/4)  = [-1/8, 0, 0]
G(0,1,0) - G(1,1,0) = [-52, -60, -2] | all even: True   (char-2 collision)
G mod 2 = [z, x*z + y, x**3*z + x**2*y + x]              (matches their comment)
```

Although the conjugation uses 1/2, G itself has integer coefficients (2, 4, 3, 12), so it is well-defined over every field, including char 2; `jacobianDet_G = 1` is a polynomial identity over ℤ, hence holds in all characteristics (their `ring` closes it uniformly). In char 2 the collision is (0,1,0) vs (1,1,0), genuinely distinct in every field (0 ≠ 1), proved via `linear_combination` with coefficients (−26, −30, −1) on `(2 : K) = 0`, which matches the componentwise differences (−52, −60, −2) exactly.

### 5.4 Statement strength

```lean
theorem not_jacobianConjecture {K : Type*} [Field K] (h2 : (2 : K) ≠ 0) :
    ¬ ∀ F : Fin 3 → MvPolynomial (Fin 3) K, IsUnit (jacobianDet F) → Injective (evalMap F)

theorem not_jacobianConjecture_all_char (K : Type*) [Field K] :
    ¬ ∀ F : Fin 3 → MvPolynomial (Fin 3) K, IsUnit (jacobianDet F) → Injective (evalMap F)

theorem not_jacobianConjecture_complex :
    ¬ ∀ F : Fin 3 → MvPolynomial (Fin 3) ℂ, IsUnit (jacobianDet F) → Injective (evalMap F)
```

- The only hypothesis on the headline theorem is `[Field K]`: no characteristic assumption, no decidability instances, no smuggled side conditions. The statement genuinely says: for EVERY field K there exists a polynomial self-map of K³ with unit Jacobian determinant that is not injective. That is the correct formal negation of the (dimension-3 instance of the) Jacobian Conjecture property over K; a dim-3 counterexample refutes the general conjecture.
- `IsUnit (jacobianDet F)` lives in `MvPolynomial (Fin 3) K`, where units are the nonzero constants, the correct notion. For the det=−2 form the unit proof correctly requires 2 ≠ 0; the det=1 form uses `isUnit_one`.
- Caveat on mathematical (not formal) novelty: for char p > 0 the Jacobian Conjecture was already classically false (e.g., x − xᵖ in dimension 1), so the "all characteristics" breadth is not itself news; the mathematically new content is char 0 (and specifically ℂ), which the repo covers via `not_jacobianConjecture` and `not_jacobianConjecture_complex`. As a *formal statement*, though, all-char is proven and strictly broader.
- Comparison with ours (`/home/alex/code/jc/JC.lean`): ours proves the det=−2 map non-injective over ℚ with the same genuine pderiv-matrix determinant. Theirs subsumes ours in stated generality (arbitrary field, char-2 handled via G, plus the ℂ specialization).

## 6. Conclusion: clean

The repository at commit `0d4a921` is a sound, faithful, and complete Lean 4 formalization of Alpöge's counterexample to the Jacobian Conjecture, and it is stronger than a ℂ-only formalization. `lake build` on the pinned toolchain (Lean `v4.33.0-rc1`, mathlib `79d0395a`) completes successfully with zero errors (one cosmetic missing-copyright-header lint warning); the source contains no `sorry`, `admit`, `native_decide`, custom `axiom`, `unsafe`, `@[implemented_by]`, or `partial def`; and an independent `#print axioms` run over every definition, lemma, and headline theorem reports exactly `[propext, Classical.choice, Quot.sound]`. The map `F` is character-for-character Alpöge's announced map (sympy re-verified on this machine, including det J(F) = −2 from a `pderiv`/`Matrix.det` encoding), the witnesses are Alpöge's exact points, distinct, and the det-1 form `G` is exactly diag(1/2,1/2,−1/2)∘F∘diag(1,2,2) with integer coefficients, det ≡ 1 in all characteristics, and a valid char-2 collision (0,1,0) ≠ (1,1,0), making `not_jacobianConjecture_all_char (K) [Field K]` proven with no additional hypotheses. Our appropriate public posture: their repo is the first formalization and it is clean; ours (`/home/alex/code/jc/JC.lean`, over ℚ) is an independent second formalization, and this document is the build audit of theirs.
