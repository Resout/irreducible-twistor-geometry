# Formalization Roadmap
## Irreducible Twistor Geometry — J. R. Manuel (2026)

Every theorem from the paper, sorted by what it needs.
Check each box when the Lean proof compiles or the Python script runs.

**Status: 157/~200 theorems verified. 7 Lean files, 1227 lines, zero sorry.**

---

## LAYER A: VERIFIED (138 theorems, done)

These compile in `lean_verification/TwistorVerified/`. Zero sorry. Zero warnings.

### SelfConsistency.lean (68 theorems)
- [x] Self-consistency: (H−1)² = H+1 at H=3
- [x] Self-consistency unique: H=3 is the only integer ≥ 2 satisfying (H−1)² = H+1
- [x] Quartic: n(n−1)(n+2)(n+3) = 0 has unique positive root n=1
- [x] H from n=1: (n+1)²−1 = 3
- [x] K* = 7/30 as rational (norm_num over ℚ)
- [x] η = 4/27 as rational
- [x] Conservation law at H=3: K*(H²+1) − η·H² = 1 (over ℚ)
- [x] Conservation law universal: holds for all H ≠ 0 over ℚ (field_simp; ring)
- [x] Conservation numerator: (H²−H+1) − (H−1)² = H for all H over ℤ (ring)
- [x] Born floor: 1/H³ = η/(H+1) at H=3
- [x] Born floor value: 1/27
- [x] Born denominator: H³−1 = 26
- [x] Sym² dimension: (H+1)(H+2)/2 = 10 at H=3
- [x] Self-consistency fails at H=2,4,5,...,10 (8 theorems)
- [x] Quartic nonzero at n=2,3,...,10 (9 theorems)
- [x] Quadratic factor: u²−5u+4=0 ↔ u=1 or u=4
- [x] Monotonicity: (H−1)² > H+1 for H ≥ 4 (over ℤ, nlinarith)
- [x] All dimensional identities, instanton 810/7, 5!=120, Koide 2/3, S₃, angle 1/30
- [x] Gaussian gap 30/23, Koide angle, mass ratio exponents

### DSCombination.lean (15 theorems)
- [x] MassFn structure defined over ℚ
- [x] DS combination rule defined (ds_pre, ds_combine)
- [x] Structural filter: pre-norm total = 1−K (algebraic identity)
- [x] Conflict symmetric: K(m,e) = K(e,m)
- [x] Agreement symmetric
- [x] DS pre-normalisation commutative
- [x] Conflict decomposition: K = (1−θ)(1−φ) − agreement
- [x] K* = 7/30 bounds, filter rate 23/30, contraction
- [x] (23/30)^n < 1 for n ≥ 1 (pow_le_pow_of_le_one)
- [x] Structural decay: pre-norm total < 1 when K > 0
- [x] Pre-norm positive when K < 1

### BornFloor.lean (10 theorems)
- [x] Born probability numerator/denominator defined
- [x] Floor equality condition: s₁²+s₂²+s₃² = 26θ²
- [x] det(M) = −25θ²/2 at floor equality
- [x] det(M) ≠ 0 at floor equality when θ ≠ 0
- [x] det(M) < 0 at floor (timelike)
- [x] Floor preserves L₁
- [x] K < 1 for non-negative masses with θ,φ > 0 (conflict_lt_one)
- [x] DS preserves L₁ (pre-norm / (1-K) = 1)
- [x] θ contracts when φ+K < 1

### Contraction.lean (19 theorems)
- [x] Singleton ratio multiplicativity (no-ignorance case)
- [x] (23/30)^10 < 1/10, ^20 < 1/100, ^50 < 1/10000
- [x] Contraction rate κ < 1 at equilibrium (explicit bound)
- [x] κ > 0
- [x] κ^100 < 0.012, κ^200 < 0.0002
- [x] Spectral radius bound: ρ ≤ 23/30 < 1
- [x] Mass gap exists: 0 < 1 − 23/30
- [x] Gaussian gap approximation: 30/23 within 3.2% of 1.263
- [x] Approximate equilibrium L₁ checks (m* and e*)
- [x] Approximate K within 0.001 of 7/30
- [x] Complete mass gap chain theorem (composite)

### PauliEmbedding.lean (20 theorems)
- [x] Minkowski quadratic Q = θ² − Σsᵢ² defined
- [x] Q ≤ 0 when θ=0
- [x] Q = −25θ² at Born floor
- [x] Q ≠ 0 at floor when θ ≠ 0
- [x] Q < 0 at floor (timelike)
- [x] DS component formulas: θ'' = θφ, s''ᵢ = sᵢeᵢ + sᵢφ + θeᵢ
- [x] DS-Matrix identity: s''ᵢ = θeᵢ + φsᵢ + sᵢeᵢ (ring)
- [x] Commutator [M,E] absent from DS output (proven component-wise)
- [x] K = symmetric off-diagonal sum (K_is_symmetric_offdiag)
- [x] Cross product defined: (s×e)ₖ = sᵢeⱼ − sⱼeᵢ
- [x] Cross product antisymmetric: s×e = −(e×s)
- [x] Self-commutator vanishes: s×s = 0
- [x] Off-diagonal split: sᵢeⱼ = (symmetric + antisymmetric)/2
- [x] Q self-combination formula

### NonHolomorphic.lean (6 theorems)
- [x] DS combination is polynomial (holomorphic)
- [x] Born condition cleared: 26θ² ≥ Σsᵢ²
- [x] Born quadratic: t²(26S²−Sq) + 2tSq − Sq = 0
- [x] Leading coefficient 26S²−Sq > 0 for positive masses
- [x] θ drops after DS step (verified at equilibrium values)
- [x] Floor fires at equilibrium: 26θ² < (1−θ)²

### BudgetMatching.lean (19 theorems)
- [x] Budget cubic cross-multiplied form
- [x] Cubic form: 7H³−30H²+37H−30 = 0
- [x] Cubic factors: (H−3)(7H²−9H+10)
- [x] Discriminant: 81−280 = −199 < 0
- [x] Quadratic 7H²−9H+10 > 0 for all H (nlinarith)
- [x] Budget unique: H=3 is the only rational root
- [x] K_cons(3) = 7/30 (verification)
- [x] Sym² equivalence: (H+1)(H+2)/2 = H²+1 iff H²−3H = 0
- [x] H²−3H = 0 iff H=0 or H=3
- [x] Sym² at H=3: 4·5/2 = 10 = 9+1
- [x] Partial fraction: K* = 1/H − 1/(H²+1) universal over ℚ
- [x] K partial at 3: 1/3 − 1/10 = 7/30
- [x] η(2) = 1/8, η(3) = 4/27, η(4) = 9/64
- [x] η(3) > η(2): 4/27 > 1/8
- [x] η(3) > η(4): 4/27 > 9/64
- [x] Integer peak: η(3) > η(2) ∧ η(3) > η(4)
- [x] All four routes to H=3 machine-verified

---

## LAYER B: ALGEBRA (doable now with Mathlib)

Concrete algebra and computation on ℚ⁴. No differential geometry.

### Part I — Uniqueness (partially done)

- [ ] **Complex encoding uniqueness** (thm:complex_unique)
  - ℂ is unique normed division algebra with ordering-sensitive encoding.
  - Needs: Hurwitz classification (axiomatize), Skolem-Noether for ℍ (in Mathlib).
  - Paper: §3.3, lines 334–348.

- [ ] **Unique product theorem** (thm:unique_product)
  - DS rule is unique bilinear product under 5 axioms (0 free params).
  - Proof: parameter elimination from 10 → 3 → 2 → 0.
  - Paper: §3.4, lines 360–384.

- [x] **Budget matching cubic** (thm:budget) ✓ BudgetMatching.lean
  - K_cons(H) = 7/30 iff H=3. Cubic factors, disc = −199 < 0.

- [x] **Sym² characterization** (thm:sym2) ✓ BudgetMatching.lean
  - dim(Sym²(ℂ^(H+1))) = H²+1 iff H²−3H=0 iff H=3.

- [x] **Efficiency peak** (thm:optimal) ✓ BudgetMatching.lean
  - η(3) > η(2) and η(3) > η(4). Integer optimum verified.

- [ ] **Efficiency robustness** (thm:robust)
  - Integer optimum H=3 for β ∈ (0.82, 1.42).
  - Paper: §3.2, lines 244–250.

### Part A — DS Dynamics

- [ ] **Phase sensitivity** (thm:phase)
  - Explicit example: two mass vectors with same Born, different DS outputs.
  - Easy: `norm_num` on concrete ℚ values.
  - Paper: §3.3, lines 316–322.

- [ ] **Self-sorting** (thm:selfsort)
  - Ratios m_i/m_j evolve as (e_i/e_j)^n. Partially done (ratio_multiplicative).
  - Paper: §4.2, lines 450–456.

- [x] **Partial fraction of K*** ✓ BudgetMatching.lean
  - K* = 1/H − 1/(H²+1). Universal over ℚ.

### Part C — det(M) Protection

- [ ] **Light cone repulsion** (thm:lightcone)
  - det(M)=0 implies det(M'') < 0 after DS step. R < 0 for positive masses.
  - Algebraic: expand Q_pre, verify all terms in R have negative coefficient.
  - Paper: §4.4, lines 958–976.

- [ ] **Self-entanglement excludes det=0** (thm:selfentangle)
  - Fixed-point equations at det=0 force contradiction.
  - Case analysis: θ_A ≠ 0 and θ_A = 0 both lead to det ≠ 0.
  - Paper: §4.4, lines 1033–1045.

### Part D — Conformal Breaking

- [ ] **Conformal symmetry breaking** (thm:conformal_break)
  - Explicit counterexample: m=(0.8,0.05,0.05,0.1), swap s₁↔θ.
  - Floor activates for m but not for swapped m.
  - Easy: `norm_num` on concrete values.
  - Paper: §8.1, lines 1429–1443.

- [ ] **Condensate coefficient** (thm:condensate)
  - C·det(M*)² = 8332/625 where 8332 = 4(3·26²+2·26+3).
  - Paper: §5.3, lines 1371–1387.

### Part D — Spin Decomposition

- [ ] **SO(4) decomposition of End(ℝ⁴)**
  - (2,2)⊗(2,2) = (1,1)⊕(3,1)⊕(1,3)⊕(3,3). Dimensions 1+3+3+9=16.
  - Trace + antisymmetric + symmetric traceless = complete.
  - Could state as dimensional identity: 1+3+3+9=16.
  - Paper: §8.3, lines 1454–1474.

---

## LAYER C: GEOMETRY (axioms + consequences)

### External Theorems — Declare as `axiom` in Lean

- [ ] Hurwitz's theorem (1898)
- [ ] Gleason's theorem (1957)
- [ ] Petz classification (1996)
- [ ] BMU reconstruction (2014)
- [ ] Ward correspondence (1977)
- [ ] Birkhoff-Grothendieck (1957)
- [ ] Newlander-Nirenberg
- [ ] Mason's framework (2005)
- [ ] Popov confirmation (2021)
- [ ] Hitchin minitwistor (1982)
- [ ] Beale-Kato-Majda criterion (1984)
- [ ] Osterwalder-Schrader reconstruction (1973/75)
- [ ] Maldacena reduction (2011)
- [ ] Adamo-Mason (2014)
- [ ] Perron-Frobenius for positive operators

### Consequences from Axioms (once stubs exist)

- [ ] From Mason axiom + verified ∂̄Φ≠0: F⁺ ≠ 0
- [ ] From Perron-Frobenius axiom + DS positivity: OS2
- [ ] From OS axiom stubs + verified Δ>0: QFT with mass gap
- [ ] From Maldacena axiom + OS2: Einstein gravity
- [ ] From BKM axiom + Born floor bound: conditional NS regularity

### Geometry to Build (long-term, collaborative)

- [ ] Phase 1: ℂPⁿ, line bundles, section dimensions
- [ ] Phase 2: Holomorphic vector bundles, ∂̄-operator
- [ ] Phase 3: Twistor space, double fibration
- [ ] Phase 4: Ward correspondence (full proof)
- [ ] Phase 6: Wirtinger calculus on ℂ⁴
- [ ] Phase 7: Almost complex structures, Nijenhuis, Chern-Simons
- [ ] Phase 8: Penrose transform, sheaf cohomology
- [ ] Phase 9: Koopman operator, spectral theory on L²(B)
- [ ] Phase 10: OS axioms as Lean definitions
- [ ] Phase 11: OS reconstruction theorem

---

## LAYER D: COMPUTATIONAL (Python, 287 scripts)

Verified numerically. Scripts at `computations/` on GitHub.

- [ ] λ₀ = 0.2829 to 500 digits (exact analytical Jacobian)
- [ ] λ₁ = 0.2813 to 500 digits
- [ ] det(I−J) = 0.5154 to 50 digits
- [ ] Λ ≈ 7.76 × 10⁻⁵¹
- [ ] PSLQ: no minimal polynomial for λ₀ up to degree 30
- [ ] PSLQ: degree-19 candidate refuted at 500 digits
- [ ] Equilibrium degree-24 polynomial (Gröbner basis)
- [ ] K* = 0.2320 ± 0.0033 across 105 lattice analyses
- [ ] SU(2) and SU(3) lattice verification
- [ ] U(1) control (same K*, different mechanism)
- [ ] det(M) protection: min |det| = 0.00084 across 10⁶ steps
- [ ] Rank-2 vs rank-1 coupling comparison
- [ ] Penrose residue ‖ρ₋₁‖ = 0.638 at 120 digits
- [ ] |F⁺|² = 0.407
- [ ] Non-abelian response ‖[R₀,R₁]‖ = 0.386
- [ ] Correlator decay Δ_eff = 1.2625
- [ ] SO(4) fractions: 23%/26%/51%
- [ ] DS-Matrix identity to 10⁻¹⁵ across 50,000 trials
- [ ] Non-triviality: kurtosis 9σ, G₄ 14.5σ
- [ ] ADE classification: 36 groups, all Δ > 0
- [ ] Folding invariance: 10 pairs to 10⁻⁹
- [ ] Descended equation: 0% Levi-Civita content
- [ ] Penrose gauntlet: 10/14 pass
- [ ] All 287 scripts

---

## LAYER E: PHYSICAL IDENTIFICATIONS (not formalizable)

Modeling choices. Outside Lean permanently.

- [ ] One DS step = one lattice spacing
- [ ] Born floor = viscosity
- [ ] (3,1)⊕(1,3) sector = vorticity
- [ ] DS conservation = incompressibility
- [ ] Spectral gap Δ = Yang-Mills mass gap
- [ ] Λ = physical cosmological constant
- [ ] Equilibrium = YM vacuum

---

## PRIORITY ORDER

### Do first (biggest impact, least effort)
1. Layer B: Budget matching cubic, Sym² characterization, partial fraction
2. Layer B: Conformal symmetry breaking (concrete counterexample)
3. Layer B: Light cone repulsion, self-entanglement exclusion
4. Layer B: Unique product theorem (parameter elimination)
5. Layer B: Phase sensitivity (concrete computation)

### Do second (moderate effort)
6. Layer B: Efficiency peak and robustness
7. Layer B: Complex encoding uniqueness (needs Hurwitz axiom)
8. Layer B: Condensate coefficient, Frobenius identity
9. Layer B: Spin decomposition (1+3+3+9=16)

### Do third (axiom file)
10. Write GeometricAxioms.lean with all 15 external declarations
11. Write MasonConsequences.lean: axioms + algebra → F⁺ ≠ 0
12. Write SpectralConsequences.lean: axioms + contraction → mass gap

### Do fourth (Python)
13. Verify all 287 scripts run, add README
14. Ensure reproducibility (versions, seeds, precision)

### Long game (years)
15. Phases 1-3, 6, 9 (projective geometry, Wirtinger, spectral)
16. Everything else

---

## FILE STRUCTURE (current → target)

```
TwistorVerified/
├── SelfConsistency.lean          ✓ done (68)
├── DSCombination.lean            ✓ done (15)
├── BornFloor.lean                ✓ done (10)
├── Contraction.lean              ✓ done (19)
├── PauliEmbedding.lean           ✓ done (20)
├── NonHolomorphic.lean           ✓ done (6)
├── BudgetMatching.lean           ✓ done (19)
├── UniqueProduct.lean            → 5-axiom uniqueness, parameter elimination
├── DetProtection.lean            → light cone repulsion, self-entanglement
├── ConformalBreaking.lean        → explicit counterexample, DS holomorphic
├── GeometricAxioms.lean          → Ward, Mason, Gleason, etc. (15 axioms)
├── MasonConsequences.lean        → from axioms: F⁺ ≠ 0, YM
├── SpectralConsequences.lean     → from axioms: OS, mass gap
└── COMPLETE_VERIFICATION.lean    ✓ master file (1072 lines)
```

---

*Last updated: 2026-04-07*
*Lean: 4.30.0-rc1 + Mathlib*
*Verified: 157 theorems, zero sorry, 7 files*
*Paper: 80 proved statements, 88 pages*
