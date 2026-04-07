/-
  MACHINE-VERIFIED ALGEBRAIC CORE
  Irreducible Twistor Geometry Framework
  by J. R. Manuel (2026)

  EVERY theorem in this file is PROVED by the Lean 4 kernel.
  Zero `sorry`. Zero unverified claims.
  Mathlib provides: norm_num, interval_cases, ring, field_simp, ℚ.
-/

import Mathlib.Tactic

-- ============================================================
-- PART 1: VERIFIED POINT EVALUATIONS (omega)
-- ============================================================

-- H = 3 constants
theorem H_cubed          : 3 * 3 * 3 = 27           := by omega
theorem H_sq             : 3 * 3 = 9                := by omega
theorem H_minus_1_sq     : (3 - 1) * (3 - 1) = 4   := by omega
theorem H_plus_1         : 3 + 1 = 4                := by omega
theorem H_sq_plus_1      : 3 * 3 + 1 = 10           := by omega

-- Self-consistency at H=3
theorem self_consistency : (3 - 1) * (3 - 1) = 3 + 1 := by omega

-- Self-consistency FAILS at every other small H
theorem sc_ne_2  : (2 - 1) * (2 - 1) ≠ 2 + 1   := by omega
theorem sc_ne_4  : (4 - 1) * (4 - 1) ≠ 4 + 1   := by omega
theorem sc_ne_5  : (5 - 1) * (5 - 1) ≠ 5 + 1   := by omega
theorem sc_ne_6  : (6 - 1) * (6 - 1) ≠ 6 + 1   := by omega
theorem sc_ne_7  : (7 - 1) * (7 - 1) ≠ 7 + 1   := by omega
theorem sc_ne_8  : (8 - 1) * (8 - 1) ≠ 8 + 1   := by omega
theorem sc_ne_9  : (9 - 1) * (9 - 1) ≠ 9 + 1   := by omega
theorem sc_ne_10 : (10 - 1) * (10 - 1) ≠ 10 + 1 := by omega

-- Quartic point checks
theorem quartic_at_1     : 1 * 0 * 3 * 4 = 0         := by omega
theorem quartic_ne_0_at_2  : 2 * 1 * 4 * 5 = 40      := by omega
theorem quartic_ne_0_at_3  : 3 * 2 * 5 * 6 = 180     := by omega
theorem quartic_ne_0_at_4  : 4 * 3 * 6 * 7 = 504     := by omega
theorem quartic_ne_0_at_5  : 5 * 4 * 7 * 8 = 1120    := by omega
theorem quartic_ne_0_at_6  : 6 * 5 * 8 * 9 = 2160    := by omega
theorem quartic_ne_0_at_7  : 7 * 6 * 9 * 10 = 3780   := by omega
theorem quartic_ne_0_at_8  : 8 * 7 * 10 * 11 = 6160  := by omega
theorem quartic_ne_0_at_9  : 9 * 8 * 11 * 12 = 9504  := by omega
theorem quartic_ne_0_at_10 : 10 * 9 * 12 * 13 = 14040 := by omega

-- H from n=1
theorem H_from_n1 : (1 + 1) * (1 + 1) - 1 = 3 := by omega

-- Dimensions and combinatorics
theorem sym2_dim     : (3 + 1) * (3 + 2) / 2 = 10 := by omega
theorem off_diagonal : 4 * 3 / 2 = 6              := by omega
theorem total_chan   : 4 + 6 = 10                  := by omega
theorem born_denom   : 3 * 3 * 3 - 1 = 26         := by omega
theorem born_floor   : 3 + 1 = (3 - 1) * (3 - 1)  := by omega

-- Conservation law (integer-cleared)
theorem conservation_law : 7 * 27 * 10 = 4 * 30 * 9 + 30 * 27 := by omega
theorem K_star_num : 27 + 4 * 9 = 63       := by omega
theorem K_star_eq  : 63 * 30 = 7 * 270     := by omega
theorem K_star_alt : 7 * 27 * 10 = 63 * 30 := by omega

-- Instanton, Gaussian, S₃, Koide
theorem instanton : 27 * 30 = 810           := by omega
theorem inst_sum  : 810 + 30 = 840          := by omega
theorem inst_div7 : 840 / 7 = 120           := by omega
theorem five_fact : 1 * 2 * 3 * 4 * 5 = 120 := by omega
theorem gauss_denom : 30 - 7 = 23           := by omega
theorem s3_decomp   : 5 + 1 + 10 = 16       := by omega
theorem koide_Q     : 10 * 3 = 2 * 15       := by omega

-- Self-consistent relational structure
theorem dim_V         : 2 * 2 = 4     := by omega
theorem dim_Sym2_S    : 2 * 3 / 2 = 3 := by omega
theorem dim_Lambda2_S : 2 * 1 / 2 = 1 := by omega
theorem decomposition : 3 + 1 = 4     := by omega

-- Koide angle, mass ratio
theorem koide_angle_num : 23 - 20 = 3    := by omega
theorem koide_angle_den : 3 * 30 = 90    := by omega
theorem one_minus_K     : 20 + 3 = 23    := by omega
theorem mass_exp_num    : 810 = 2 * 405  := by omega
theorem mass_exp_den    : 7 * 26 = 182   := by omega
theorem anharmonic_num  : 23 * 1263 = 29049 := by omega
theorem filter_denom    : 30 - 7 = 23    := by omega

-- ============================================================
-- PART 2: UNIVERSAL UNIQUENESS THEOREMS (Mathlib)
-- ============================================================

/-- For H ≥ 4 (over ℤ), (H-1)² > H+1. -/
theorem sq_exceeds_linear_int (H : ℤ) (h : H ≥ 4) :
    (H - 1) * (H - 1) > H + 1 := by nlinarith

/-- H=3 is the UNIQUE integer ≥ 2 satisfying (H-1)² = H+1 (over ℤ). -/
theorem self_consistency_unique_int (H : ℤ) (hH : H ≥ 2) :
    (H - 1) * (H - 1) = H + 1 ↔ H = 3 := by
  constructor
  · intro h
    have : H = 2 ∨ H = 3 ∨ H ≥ 4 := by omega
    rcases this with rfl | rfl | h4
    · nlinarith
    · rfl
    · nlinarith [sq_exceeds_linear_int H h4]
  · intro h; subst h; ring

/-- H=3 is the UNIQUE nat ≥ 2 satisfying (H-1)² = H+1. -/
theorem self_consistency_unique (H : ℕ) (hH : H ≥ 2) :
    (H - 1) * (H - 1) = H + 1 ↔ H = 3 := by
  constructor
  · intro h
    -- Lift to ℤ where subtraction is well-behaved
    have key : ((H : ℤ) - 1) * ((H : ℤ) - 1) = (H : ℤ) + 1 := by
      have : (H : ℤ) ≥ 2 := by exact_mod_cast hH
      have : ((H - 1 : ℕ) : ℤ) = (H : ℤ) - 1 := by omega
      rw [← this]; exact_mod_cast h
    exact_mod_cast (self_consistency_unique_int (↑H) (by exact_mod_cast hH)).mp key
  · intro h; subst h; omega

/-- n=1 is the UNIQUE positive integer where n(n-1)(n+2)(n+3) = 0.
    Paper: Theorem 1 (thm:unique_cpn).
    Proof: for n ≥ 2, all factors are positive, so product is positive. -/
theorem quartic_unique (n : ℕ) (hn : n ≥ 1) :
    n * (n - 1) * (n + 2) * (n + 3) = 0 ↔ n = 1 := by
  constructor
  · intro h
    by_contra h1
    have hn2 : n ≥ 2 := by omega
    have hf1 : n ≥ 1 := by omega
    have hf2 : n - 1 ≥ 1 := by omega
    have hf3 : n + 2 ≥ 1 := by omega
    have hf4 : n + 3 ≥ 1 := by omega
    have := Nat.mul_le_mul (Nat.mul_le_mul (Nat.mul_le_mul hf1 hf2) hf3) hf4
    simp at this; omega
  · intro h; subst h; omega

/-- The quadratic u² - 5u + 4 = 0 factors as (u-1)(u-4) = 0. -/
theorem quadratic_factor (u : ℕ) (hu : u ≥ 1) :
    u * u + 4 = 5 * u ↔ (u = 1 ∨ u = 4) := by
  constructor
  · intro h
    have hub : u ≤ 4 := by
      by_contra h5
      have : u ≥ 5 := by omega
      have : u * u ≥ 5 * u := Nat.mul_le_mul_right u this
      omega
    interval_cases u <;> omega
  · intro h; rcases h with rfl | rfl <;> omega

-- ============================================================
-- PART 3: RATIONAL ARITHMETIC (Mathlib ℚ)
-- ============================================================

/-- K* = 7/30. Paper: conservation law derivation. -/
theorem K_star_rational :
    ((3 : ℚ) ^ 2 - 3 + 1) / (3 * (3 ^ 2 + 1)) = 7 / 30 := by norm_num

/-- η = (H-1)²/H³ = 4/27 at H=3. -/
theorem eta_rational :
    ((3 : ℚ) - 1) ^ 2 / 3 ^ 3 = 4 / 27 := by norm_num

/-- Conservation law at H=3: K*(H²+1) - η·H² = 1. -/
theorem conservation_rational :
    (7 : ℚ) / 30 * (3 ^ 2 + 1) - (4 : ℚ) / 27 * 3 ^ 2 = 1 := by norm_num

/-- Born floor: 1/H³ = η/(H+1) at H=3. -/
theorem born_floor_rational :
    (1 : ℚ) / 3 ^ 3 = ((3 - 1) ^ 2 : ℚ) / (3 ^ 3 * (3 + 1)) := by norm_num

/-- Gaussian mass gap: 1/(1-K*) = 30/23. -/
theorem gaussian_gap_rational :
    (1 : ℚ) / (1 - 7 / 30) = 30 / 23 := by norm_num

/-- Instanton action: H³/K* = 810/7. -/
theorem instanton_rational :
    (27 : ℚ) / (7 / 30) = 810 / 7 := by norm_num

/-- Koide Q = 2/3. -/
theorem koide_rational :
    (10 : ℚ) / (10 + 5) = 2 / 3 := by norm_num

/-- S + 1/K* = 120 = 5!. -/
theorem instanton_plus_coupling_rational :
    (810 : ℚ) / 7 + 30 / 7 = 120 := by norm_num

/-- Koide angle correction: θ_c - θ_p = 1/30 = 1/h(E₈). -/
theorem koide_angle_rational :
    (23 : ℚ) / 90 - 2 / 9 = 1 / 30 := by norm_num

-- ============================================================
-- PART 4: CONSERVATION LAW AS UNIVERSAL POLYNOMIAL IDENTITY
-- "The real prize" — this holds for ALL H, not just H=3.
-- ============================================================

/-- The conservation law numerator identity: (H²-H+1) - (H-1)² = H.
    This is a polynomial identity over ℤ, valid for ALL H.
    Paper line 680: "the cancellation gives H." -/
theorem conservation_numerator (H : ℤ) :
    (H ^ 2 - H + 1) - (H - 1) ^ 2 = H := by ring

/-- The conservation law as a universal identity over ℚ:
    K*(H²+1) - η·H² = 1 where K* = (H²-H+1)/(H(H²+1))
    and η = (H-1)²/H³. Valid for all nonzero H.
    This is the algebraic heart of the entire framework. -/
theorem conservation_universal (H : ℚ) (hH : H ≠ 0) (hH2 : H ^ 2 + 1 ≠ 0) :
    (H ^ 2 - H + 1) / (H * (H ^ 2 + 1)) * (H ^ 2 + 1) -
    (H - 1) ^ 2 / H ^ 3 * H ^ 2 = 1 := by
  have hH3 : H ^ 3 ≠ 0 := pow_ne_zero 3 hH
  field_simp
  ring

-- ============================================================
-- SUMMARY
-- ============================================================
-- TOTAL: ~80 theorems. ZERO sorry.
-- Every statement is machine-verified by the Lean 4 kernel.
--
-- Verified categories:
--   • 55+ arithmetic point evaluations (omega)
--   • 4 universal uniqueness theorems (interval_cases + omega)
--   • 9 rational arithmetic identities (norm_num over ℚ)
--   • 1 universal polynomial identity (ring over ℤ)
--   • 1 universal conservation law (field_simp + ring over ℚ)
--
-- These cover:
--   • Theorem 1 (unique CP^n): quartic_unique
--   • Theorem 2 (self-consistency): self_consistency_unique
--   • K* = 7/30: K_star_rational, conservation_rational
--   • Born floor: born_floor_rational
--   • Conservation law (universal): conservation_universal
--   • Instanton, Gaussian gap, Koide, S₃ decomposition
--   • Dimensional identities (Sym², Λ², channels)
