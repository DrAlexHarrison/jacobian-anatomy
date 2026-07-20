# CAS cross-verification of Alpöge's Jacobian-Conjecture counterexample

**Goal.** Re-derive every claim about the counterexample map in **two computer
algebra systems that are independent of each other and of the original sympy
verification** — [Singular](https://www.singular.uni-kl.de/) 4.3.2 and
[PARI/GP](https://pari.math.u-bordeaux.fr/) 2.15.4 — and guard every claim with a
**negative control** (a deliberately broken input that the same machinery must
reject). Sympy appears in **no** script here; independence is the whole point.

The map `F : C^3 -> C^3` (Alpöge, announced 2026-07-19, found with Claude Fable):

```
F1 = (1+xy)^3 z + y^2 (1+xy)(4+3xy)
F2 = y + 3x(1+xy)^2 z + 3 x y^2 (4+3xy)
F3 = 2x - 3 x^2 y - x^3 z
```

**Bottom line: 36/36 checks PASS, 0 FAIL.** Every claim reproduced in Singular
and (where feasible) PARI; every negative control detected its break.
Run everything with:

```
bash /home/alex/code/jc/verify/cas/run-all.sh
```

Independence is genuine at the algorithm level, not just the tool level:
determinants are taken by cofactor expansion; elimination is done **two
different ways** — Gröbner-basis elimination (Singular) *and* iterated
resultants (PARI); the discriminant is taken by `poldisc` (PARI) and
cross-checked by the resultant identity `Res(C,C') = -a_n·disc` (Singular);
fibers are counted by Gröbner `vdim` (Singular) *and* by explicit numeric
root-solving (PARI).

---

## Results table

| # | Claim | Singular | PARI/GP | Negative control |
|---|-------|:--------:|:-------:|------------------|
| 1 | `det Jacobian(F) ≡ −2` (Keller map) | **PASS** | **PASS** | `(4+3xy)→(4+2xy)` ⇒ det is a degree-13 polynomial, **not** constant — **detected** |
| 2 | Fiber over `(−1/4,0,0)`: 3 witnesses map there; complete fiber ideal is 0-dim with `vdim = 3` | **PASS** | **PASS** | witness `z=13/2→13/3` ⇒ `F(1,−3/2,13/3)=(1/48,−13/8,13/6)≠target` — **detected** |
| 3 | Generic fiber count `= 3` at 3 random targets (`vdim = 3`) | **PASS** | **PASS** | dropping `F1` ⇒ `⟨F2,F3⟩` has `dim = 1` (not 0-dim) — **detected** |
| 4 | Eliminants: `elim(y,z)=D·x³+(4−3bc)x−2c`, `elim(x,z)=2y³−3by²+18ay+(27a²c−18ab+b³)` | **PASS** | **PASS** | perturbed map ⇒ `elim(y,z)` is a **degree-7** poly, not the cubic; resultant not divisible by claimed cubic — **detected** |
| 5 | `disc_x(cubic) = −4·(27ac²−9bc+8)²·D` (square factor) | **PASS** | **PASS** | `D→D+1` ⇒ discriminant loses the square-factor identity — **detected** |
| 6 | Empty fibers on Γ = unit ideal; single point has `vdim = 1` | **PASS** | **PASS** | generic target ⇒ **not** unit ideal (`vdim = 3`); generic eliminant is degree 3, not a constant — **detected** |
| 7 | `G_m`-equivariance `F(tx,y/t,z/t²)=(F1/t²,F2/t,t·F3)` | **PASS** | **PASS** | wrong weight `z→z·t` ⇒ identity residual `=(t³−1)·z(1+xy)³ ≠ 0` — **detected** |
| 8 | Punt covering identity `2F1·r³−F2·r²+2r−F3=0`, `r=x/(1+xy)` | **PASS** | **PASS** | coefficient `2F1→3F1` ⇒ numerator `≠ 0` — **detected** |

**PASS tokens: 36 · FAIL tokens: 0.** (Singular positives 11, Singular negatives
8, PARI positives 10, PARI negatives 7.)

---

## Per-script detail, tool, and wall time

All wall times are from a single `run-all.sh` pass on `alex-desktop`
(warm binaries). Singular runs one script per claim; the PARI positives and
negatives are each one bundled script.

| Script | Tool | Wall | Outcome |
|--------|------|-----:|---------|
| `claim1_det.sing` | Singular | 0.014 s | `det = −2` exactly; `d == -2` ⇒ PASS |
| `claim2_fiber.sing` | Singular | 0.010 s | 3 witnesses → `(−1/4,0,0)`; fiber ideal `dim 0, vdim 3` ⇒ PASS |
| `claim3_generic_fiber.sing` | Singular | 0.011 s | `dim 0, vdim 3` at all 3 generic targets ⇒ PASS |
| `claim4_eliminants.sing` | Singular | 0.025 s | eliminants returned **exactly** as claimed cubics (see below) ⇒ PASS |
| `claim5_discriminant.sing` | Singular | 0.009 s | `Res(C,C') − 4·D²·(27ac²−9bc+8)² = 0` ⇒ PASS |
| `claim6_special_fibers.sing` | Singular | 0.010 s | 3 empty fibers = unit ideal; single point `vdim 1` ⇒ PASS |
| `claim7_equivariance.sing` | Singular | 0.008 s | all three weighted residuals `= 0` ⇒ PASS |
| `claim8_punt.sing` | Singular | 0.011 s | cleared-denominator numerator `= 0` ⇒ PASS |
| `neg_controls.sing` | Singular | 0.707 s | 8/8 breaks detected (heaviest: degree-7 broken eliminant) |
| `pari_claims.gp` | PARI/GP | 0.011 s | claims 1–8, 10 PASS tokens |
| `pari_negatives.gp` | PARI/GP | 0.010 s | 7/7 breaks detected |

---

## What each tool actually computed (key artifacts)

**Claim 1 — Jacobian determinant.** Both CASs return the constant `−2` for
`det[∂F_i/∂x_j]` over `Q[x,y,z]`. Constant nonzero determinant ⇒ Keller map.

**Claim 4 — eliminants.** The mission allowed a match "up to factors `x^k`, `c`".
Singular's Gröbner elimination in fact returns the generators **exactly**, with
no extraneous factors, so we test *exact equality* (strictly stronger):

```
eliminate(y,z):  (27a²c²−18abc+16a+b³c−b²)·x³ + (4−3bc)·x − 2c
eliminate(x,z):  2y³ − 3b·y² + 18a·y + (27a²c − 18ab + b³)
```

The `x`-eliminant's leading coefficient is exactly `D = 27a²c²−18abc+16a+b³c−b²`;
the `y`-eliminant's leading coefficient is exactly `2`. PARI independently
derives the same `x`-eliminant by **iterated resultants** (eliminate `z` from the
pairs, then `y`) and confirms the claimed cubic divides the resultant — a
completely different elimination algorithm reaching the same locus.

**Claim 5 — discriminant.** PARI's `poldisc` factors the `x`-cubic's discriminant
as

```
factor(disc) = [ 27ac² − 9bc + 8 , 2 ;      <- squared factor
                 27a²c² − 18abc + 16a + b³c − b² , 1 ]   <- = D
```

and the exact rational constant is `−4`: verified as the polynomial identity
`disc = −4·(27ac²−9bc+8)²·D` (`disc − target == 0`). Singular corroborates via
`Res(C,C') = 4·D²·(27ac²−9bc+8)²`, i.e. `disc = −Res/D = −4·(27ac²−9bc+8)²·D`.

**Claim 6 — special fibers.** Singular is authoritative: at `(4/27,4/3,1)`,
`(1/27,2/3,2)`, `(4/243,−4/9,−3)` the ideal `⟨F−target⟩` is the **unit ideal**
(`reduce(1,std I)=0`, so `1 ∈ I`, so the fiber is empty); at `(−16/27,0,1)` it is
0-dimensional with `vdim = 1` (one point). PARI gives the geometric cross-view:
at the three empty targets the `x`-eliminant degenerates to a **nonzero constant**
(`−2`, `−4`, `6`) — no `x`-solution, hence no fiber point — while at the
single-point target it degenerates to the **linear** `4x − 2` (unique root
`x = 1/2`).

**Claim 3 — fiber counting, independent method.** Singular counts via Gröbner
`vdim = 3`. PARI counts by *solving*: it takes the three `x`-roots of the cubic,
solves `F3=c` for `z` linearly, substitutes into `F2=b` to get the `y`-roots,
and keeps only triples that also satisfy `F1=a` to machine precision — yielding
exactly 3 distinct points at every generic target.

---

## Why the negative controls matter

A verification suite that cannot fail is worthless. Each control perturbs one
ingredient and confirms the **same code path** rejects it:

* **det** — breaking the `(4+3xy)` factor makes the Jacobian non-constant
  (a degree-13 polynomial), so the Keller property genuinely depends on the
  exact coefficients.
* **fiber evaluation** — nudging a witness coordinate makes it miss the target,
  so the witness check is not vacuously true.
* **zero-dimensionality / vdim** — an underdetermined system is caught as
  positive-dimensional (`dim = 1`), so `vdim = 3` is a real finiteness statement.
* **eliminant** — the perturbed map yields a *degree-7* eliminant, so "the
  eliminant is this cubic" is a real constraint.
* **discriminant** — shifting `D` by 1 destroys the square-factor pattern, so the
  `(27ac²−9bc+8)²` structure is special to the true `D`.
* **empty vs. generic** — a generic target is *not* the unit ideal, so the
  emptiness test discriminates.
* **equivariance** — the wrong `G_m` weight leaves a nonzero residual
  `(t³−1)·z(1+xy)³`, so the weights `(−2,−1,+1)` are exactly right.
* **punt identity** — changing the leading coefficient leaves a nonzero
  numerator, so the covering relation is exact.

---

## Files

```
/home/alex/code/jc/verify/cas/
├── run-all.sh                 # orchestrator — absolute binary paths only
├── claim1_det.sing            # Jacobian det = -2
├── claim2_fiber.sing          # fiber over (-1/4,0,0): witnesses + vdim=3
├── claim3_generic_fiber.sing  # generic fiber count = 3
├── claim4_eliminants.sing     # x- and y-eliminants (Groebner)
├── claim5_discriminant.sing   # discriminant square factor (resultant)
├── claim6_special_fibers.sing # empty fibers (unit ideal) + single point
├── claim7_equivariance.sing   # G_m equivariance, weights (-2,-1,+1)
├── claim8_punt.sing           # Punt covering identity
├── neg_controls.sing          # 8 Singular negative controls
├── pari_claims.gp             # claims 1-8 in PARI/GP (independent CAS)
├── pari_negatives.gp          # 7 PARI negative controls
└── RESULTS.md                 # this file
```

Environment: Singular 4.3.2 (GMP 6.3.0 / NTL 11.5.1 / FLINT 3.0.1),
PARI/GP 2.15.4, on `alex-desktop` (Linux). House rules honored: absolute binary
paths (`/usr/bin/Singular`, `/usr/bin/gp`, `/bin/*`) throughout; no git
operations; nothing outside `~/code/jc/verify/cas/` touched.
