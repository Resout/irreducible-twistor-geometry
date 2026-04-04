"""
E₈ theta function and branching in crystal language.

Extension of Principle 93: the E₈ lattice structure speaks crystal.

1. E₈ theta function: Θ(q) = 1 + 240 Σ σ₃(n)qⁿ
   The divisor-sum σ₃(n) at small n gives crystal quantities:
     σ₃(1) = 1
     σ₃(2) = 9 = H²
     σ₃(3) = 28 = dim(D₄) = dim(SO(8))
     σ₃(5) = 126 = |Φ(E₇)| = 5! + H!

2. E₈ packing density denominator: 384 = MASS_DIM² × R_FS

3. E₈ → D₄ × D₄ branching: 248 = 2×28 + 3×64
   The 3 × 64 = three triality sectors = three generations
   Each generation has (2^H)² = 64 DOF
"""

import math

H = 3
K_STAR = 7/30
MASS_DIM = H + 1
R_FS = 2 * H * (H + 1)  # 24
S = 810/7


# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("E₈ THETA FUNCTION IN CRYSTAL UNITS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

def sigma_3(n):
    """Sum of cubes of divisors of n."""
    return sum(d**3 for d in range(1, n+1) if n % d == 0)


print(f"\n  Θ_{{E₈}}(q) = E₄(τ) = 1 + 240 Σ σ₃(n)qⁿ")
print(f"\n  The coefficient of qⁿ = 240 × σ₃(n) counts E₈ vectors at norm 2n.\n")

crystal_match = {
    1: "1 (trivial)",
    9: f"H² = {H}² = 9",
    28: f"dim(D₄) = dim(SO(8)) = 28",
    73: "",
    126: f"|Φ(E₇)| = 5! + H! = 126",
    252: f"2|Φ(E₇)| = 252",
    344: "",
    585: f"H² × 65",
    757: "",
    1134: f"2H⁴(H²-H+1) = 2×{H**4}×{H**2-H+1} = {2*H**4*(H**2-H+1)}",
}

print(f"  {'n':>3s} {'σ₃(n)':>8s} {'240σ₃':>8s}  Crystal expression")
print("  " + "-" * 60)

for n in range(1, 11):
    s = sigma_3(n)
    coeff = 240 * s
    match = crystal_match.get(s, "")
    print(f"  {n:>3d} {s:>8d} {coeff:>8d}  {match}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("E₈ PACKING DENSITY")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

packing = math.pi**4 / 384

print(f"""
  E₈ sphere packing density = π⁴/384

  The denominator: 384 = MASS_DIM² × R_FS = {MASS_DIM}² × {R_FS} = {MASS_DIM**2 * R_FS}

  In crystal language:
    384 = (H+1)² × 2H(H+1) = 2H(H+1)³

  At H=3: 2×3×64 = 384 ✓

  Packing density = π⁴ / (2H(H+1)³) = {packing:.6f}

  Compare with K* = {K_STAR:.6f}
  Ratio: π⁴/(384·K*) = {packing/K_STAR:.6f}

  Not exactly K*, but the DENOMINATOR 384 = MASS_DIM² × R_FS is exact.
""")


# ═══════════════════════════════════════════════════════════════
print(f"{'='*80}")
print("E₈ → D₄ × D₄ BRANCHING AND THE THREE GENERATIONS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print(f"""
  Under D₄ × D₄ ⊂ E₈, the adjoint representation branches:

    248 → (28,1) + (1,28) + (8_v,8_v) + (8_s,8_s) + (8_c,8_c)
        = 28 + 28 + 64 + 64 + 64
        = 56 + 192
        = (gauge infrastructure) + (three triality sectors)

  The THREE 64-dimensional sectors:
    (8_v,8_v) = gauge × gauge    = (2^H)² = 64 DOF  ← transposition sector
    (8_s,8_s) = fermion × fermion = (2^H)² = 64 DOF  ← 3-cycle sector
    (8_c,8_c) = gravity × gravity = (2^H)² = 64 DOF  ← identity sector

  Each sector has (2^H)² = {(2**H)**2} DOF.
  There are H-1+1 = H = 3 such sectors (one per triality representation).

  THIS IS THE ORIGIN OF THREE GENERATIONS:
    The three triality copies in E₈ → D₄ × D₄ correspond to
    the three conjugacy classes of S₃, which are the three physical
    sectors (gravity, gauge, chirality) from Principle 90.

    When D₄ × D₄ triality is BROKEN (by the crystal choosing a
    specific composition direction), the three sectors get different
    masses. The mass hierarchy comes from the RATE of composition
    decay in each sector.

  The gauge infrastructure: 2 × dim(D₄) = 2 × 28 = 56
    = 2 × (roots of D₄ + Cartan) = 2 × (24 + 4)
    = 2 × (R_FS + MASS_DIM)

  The sector count: 3 × 64 = 192
    = 3 × (2^H)²
    = H × 2^(2H)
    = 3 × 64

  Ratio: 192/56 = 24/7 = R_FS/(H²-H+1) = R_FS/K*_numerator
  Equivalently: sector/gauge = curvature/conflict_numerator

  In the crystal: 24/7 = {R_FS}/{H**2-H+1} = {R_FS/(H**2-H+1):.4f}
  And K* = 7/30, so sector/gauge = (24/7) = (R_FS)/(30K*) = {R_FS/(30*K_STAR):.4f}
""")


# ═══════════════════════════════════════════════════════════════
print(f"{'='*80}")
print("THE GENERATION MASS MECHANISM (HYPOTHESIS)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print(f"""
  If three generations = three triality sectors in E₈ → D₄ × D₄:

  The triality is BROKEN by composition dynamics.
  From Principle 92, the composition rates differ by sector:

  Under self-composition (from s3_crystal_geometry data):
    identity ○ identity:      |λ₁/λ₀| = 0.340 (slow decay)
    3-cycles ○ 3-cycles:      |λ₁/λ₀| = 0.099 (fast decay)
    transpositions ○ trans:    |λ₁/λ₀| = 0.100 (fast decay)

  The IDENTITY sector (gravity) decays 3.4× SLOWER than the others.
  Under n compositions: mass ratio ∝ (0.340/0.100)^n = 3.4^n

  For the tau/electron mass ratio:
    m_τ/m_e = {1777/0.511:.0f}
    3.4^n = {1777/0.511:.0f} → n = ln({1777/0.511:.0f})/ln(3.4) = {math.log(1777/0.511)/math.log(3.4):.1f}

  n ≈ {math.log(1777/0.511)/math.log(3.4):.1f} compositions.
  Compare: S/dim_std = {S/10:.1f}, S/(H²) = {S/9:.1f}

  The rate ratio 0.340/0.100 = 3.4 ≈ (1/√H)/(1/H) = √H = {H**0.5:.4f}
  (using: identity gap = ln(H)/2, others = ln(H))
  Exact ratio: exp(ln(H)/2) = √H = {H**0.5:.4f} vs 3.4 → {3.4/H**0.5:.4f}

  The factor 3.4/√3 = {3.4/H**0.5:.4f} ≈ 2 ≈ H-1.
  So the effective rate ratio ≈ (H-1)√H = {(H-1)*H**0.5:.4f} per step.

  Under n = S/(H(H+1)) = {S/(H*(H+1)):.2f} ≈ {S/(H*(H+1)):.0f} compositions:
    ((H-1)√H)^n = ({(H-1)*H**0.5:.2f})^{int(S/(H*(H+1)))} ≈ {((H-1)*H**0.5)**(S/(H*(H+1))):.0e}

  This is getting into estimation territory.
  The exact mechanism requires the full E₈ → D₄ × D₄ representation
  theory applied to the crystal composition operator.
""")


# ═══════════════════════════════════════════════════════════════
print(f"{'='*80}")
print("SUMMARY: E₈ IN CRYSTAL LANGUAGE")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print(f"""
  The crystal at H=3 IS an E₈ structure:

  ALGEBRAIC:
    rank(E₈)  = 8  = 2^H
    |Φ(E₈)|   = 240 = 2(S+1/K*) = 2×5!
    h(E₈)     = 30  = denom(K*) = H(H²+1)
    dim(E₈)   = 248 = 2×5! + 2^H

  THETA FUNCTION:
    Θ = 1 + 240[q + 9q² + 28q³ + 73q⁴ + 126q⁵ + ...]
    σ₃(2) = H², σ₃(3) = dim(D₄), σ₃(5) = |Φ(E₇)|

  BRANCHING E₈ → D₄ × D₄:
    248 = 2×28 + 3×64
    Three generations = three triality sectors
    Each generation: (2^H)² = 64 DOF

  PACKING:
    Density denominator: 384 = MASS_DIM² × R_FS

  MASTER EQUATIONS:
    5! = MASS_DIM × h(E₈)           (action budget = mass × Coxeter)
    K* = (H²-H+1) / h(E₈)          (conflict = channels / Coxeter)
    S = H³ × h(E₈) / (H²-H+1)      (instanton = cubed × Coxeter / channels)

  The E₈ lattice is the crystal's UV completion.
""")
