/-
A machine-checked certificate that the Jacobian Conjecture is false.

The map (Alpöge, announced 2026-07-19, found with Claude Fable):
  F₁ = (1+xy)³z + y²(1+xy)(4+3xy)
  F₂ = y + 3x(1+xy)²z + 3xy²(4+3xy)
  F₃ = 2x - 3x²y - x³z

Verified here in Lean 4 / mathlib:
  1. `jacobian_det` : the formal Jacobian determinant of F is the constant -2,
     so F is a Keller map.
  2. `evalF_p₁ / p₂ / p₃` : the three distinct points (0,0,-1/4), (1,-3/2,13/2),
     (-1,3/2,13/2) all map to (-1/4,0,0).
  3. `not_injective` : the induced polynomial map ℚ³ → ℚ³ is not injective.
  4. `jacobian_conjecture_false` : the combined existential. Since ℚ has
     characteristic zero and the same identities hold verbatim over ℂ (the
     polynomials have rational coefficients and the witnesses are rational
     points), this falsifies the Jacobian Conjecture as classically stated.
-/
import Mathlib

open MvPolynomial

noncomputable section
namespace JCCounterexample

private abbrev x : MvPolynomial (Fin 3) ℚ := X 0
private abbrev y : MvPolynomial (Fin 3) ℚ := X 1
private abbrev z : MvPolynomial (Fin 3) ℚ := X 2

/-- Alpöge's map. -/
def F : Fin 3 → MvPolynomial (Fin 3) ℚ
  | 0 => (1 + x * y) ^ 3 * z + y ^ 2 * (1 + x * y) * (4 + 3 * x * y)
  | 1 => y + 3 * x * (1 + x * y) ^ 2 * z + 3 * x * y ^ 2 * (4 + 3 * x * y)
  | 2 => 2 * x - 3 * x ^ 2 * y - x ^ 3 * z

/-- The Jacobian matrix of formal partial derivatives. -/
def J : Matrix (Fin 3) (Fin 3) (MvPolynomial (Fin 3) ℚ) :=
  Matrix.of fun i j => pderiv j (F i)

private lemma pd_one (i : Fin 3) : pderiv i (1 : MvPolynomial (Fin 3) ℚ) = 0 := by
  rw [← map_one (C : ℚ →+* MvPolynomial (Fin 3) ℚ), pderiv_C]

private lemma pd00 : pderiv 0 (x) = 1 := by simp
private lemma pd01 : pderiv 1 (x) = 0 := by simp
private lemma pd02 : pderiv 2 (x) = 0 := by simp
private lemma pd10 : pderiv 0 (y) = 0 := by simp
private lemma pd11 : pderiv 1 (y) = 1 := by simp
private lemma pd12 : pderiv 2 (y) = 0 := by simp
private lemma pd20 : pderiv 0 (z) = 0 := by simp
private lemma pd21 : pderiv 1 (z) = 0 := by simp
private lemma pd22 : pderiv 2 (z) = 1 := by simp

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- The Jacobian determinant of F is the nonzero constant polynomial -2:
F is a Keller map. -/
theorem jacobian_det : J.det = -2 := by
  have h3 : ∀ p : MvPolynomial (Fin 3) ℚ, p ^ 3 = p * p * p := fun p => by ring
  have h2 : ∀ p : MvPolynomial (Fin 3) ℚ, p ^ 2 = p * p := fun p => by ring
  have e4 : (4 : MvPolynomial (Fin 3) ℚ) = 1 + 1 + 1 + 1 := by norm_num
  have e3 : (3 : MvPolynomial (Fin 3) ℚ) = 1 + 1 + 1 := by norm_num
  have e2 : (2 : MvPolynomial (Fin 3) ℚ) = 1 + 1 := by norm_num
  simp only [J, Matrix.det_fin_three, Matrix.of_apply, F, h3, h2, e4, e3, e2]
  simp only [pderiv_mul, map_add, map_sub, pd_one,
    pd00, pd01, pd02, pd10, pd11, pd12, pd20, pd21, pd22]
  simp only [mul_zero, zero_mul, mul_one, one_mul, add_zero, zero_add, sub_zero,
    zero_sub]
  ring

/-- Evaluation of F at a point of ℚ³. -/
def evalF (v : Fin 3 → ℚ) : Fin 3 → ℚ := fun i => eval v (F i)

def p₁ : Fin 3 → ℚ := ![0, 0, -1/4]
def p₂ : Fin 3 → ℚ := ![1, -3/2, 13/2]
def p₃ : Fin 3 → ℚ := ![-1, 3/2, 13/2]
def q  : Fin 3 → ℚ := ![-1/4, 0, 0]

theorem evalF_p₁ : evalF p₁ = q := by
  funext i; fin_cases i <;> (simp [evalF, F, p₁, q]; try norm_num)

theorem evalF_p₂ : evalF p₂ = q := by
  funext i; fin_cases i <;> (simp [evalF, F, p₂, q]; try norm_num)

theorem evalF_p₃ : evalF p₃ = q := by
  funext i; fin_cases i <;> (simp [evalF, F, p₃, q]; try norm_num)

theorem p₁_ne_p₂ : p₁ ≠ p₂ := by
  intro h
  have h0 := congrFun h 0
  simp [p₁, p₂] at h0

/-- The polynomial map F : ℚ³ → ℚ³ is not injective. -/
theorem not_injective : ¬ Function.Injective evalF := by
  intro h
  exact p₁_ne_p₂ (h (evalF_p₁.trans evalF_p₂.symm))

/-- **The Jacobian Conjecture is false**: an explicit polynomial map in three
variables over a field of characteristic zero whose Jacobian determinant is the
nonzero constant -2, yet whose induced map is not injective (hence has no
polynomial inverse). -/
theorem jacobian_conjecture_false :
    ∃ G : Fin 3 → MvPolynomial (Fin 3) ℚ,
      (Matrix.of fun i j => pderiv j (G i)).det = -2 ∧
      ¬ Function.Injective fun v : Fin 3 → ℚ => fun i => eval v (G i) :=
  ⟨F, jacobian_det, not_injective⟩

end JCCounterexample
