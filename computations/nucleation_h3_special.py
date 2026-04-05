"""
At what H does crystallization happen BEFORE equilibrium?

K_peak = K*/2 + 1/(H²+1)
K* = (H²-H+1)/(H(H²+1))

K_peak < K* iff K*/2 + 1/(H²+1) < K*
             iff 1/(H²+1) < K*/2
             iff 2/(H²+1) < K*
             iff 2/(H²+1) < (H²-H+1)/(H(H²+1))
             iff 2H < H²-H+1
             iff 0 < H²-3H+1
             iff H > (3+√5)/2 ≈ 2.618

So K_peak < K* iff H ≥ 3. The smallest integer is H=3.

This means: at H=2, the crystal doesn't crystallize until AFTER
reaching equilibrium conflict. The entropy peak is a POST-equilibrium
phenomenon. At H≥3, crystallization starts BEFORE equilibrium — the
entropy peak is a genuine nucleation barrier that precedes the steady state.

At H=3 (the self-consistent dimension), the crystal begins organizing
before the conflict has fully established. Structure precipitates
from sub-equilibrium conflict. The crystal doesn't wait — it nucleates.

The critical value (3+√5)/2 is the golden ratio + 1 = φ + 1.
"""

import math

phi = (1 + math.sqrt(5)) / 2

print("=" * 80)
print("CRYSTALLIZATION BEFORE EQUILIBRIUM: H ≥ 3")
print("=" * 80)

print(f"\n  K_peak < K* iff H > (3+√5)/2 = {(3+math.sqrt(5))/2:.6f} = φ+1 = {phi+1:.6f}")
print(f"  Golden ratio φ = {phi:.6f}")
print(f"  Critical H = φ + 1 = {phi+1:.6f}")
print(f"  Smallest integer H: 3")

print(f"\n  {'H':>3s}  {'K_peak':>10s}  {'K*':>10s}  {'K_peak < K*?':>15s}  {'K_peak/K*':>10s}")
for H in range(2, 10):
    K_star = (H**2 - H + 1) / (H * (H**2 + 1))
    K_peak = (H**2 + H + 1) / (2 * H * (H**2 + 1))
    before = K_peak < K_star
    ratio = K_peak / K_star
    print(f"  {H:3d}  {K_peak:>10.6f}  {K_star:>10.6f}  {'YES (nucleates)' if before else 'NO (post-equil)':>15s}  {ratio:>10.4f}")

# The critical condition H²-3H+1 = 0 at H = (3+√5)/2
print(f"\n\n{'='*80}")
print("THE GOLDEN RATIO IN THE NUCLEATION CONDITION")
print("="*80)

print(f"""
  The condition K_peak = K* (nucleation at equilibrium) requires:
    H² - 3H + 1 = 0
    H = (3 ± √5) / 2

  The positive root is:
    H_crit = (3 + √5)/2 = φ + 1 = 1/φ + 2

  where φ = (1+√5)/2 is the golden ratio.

  Below H_crit: entropy peaks AFTER K reaches K* (post-equilibrium)
  Above H_crit: entropy peaks BEFORE K reaches K* (nucleation)

  H=3 is the smallest integer where the crystal nucleates before
  reaching equilibrium conflict. This is ANOTHER special property
  of H=3, alongside (H-1)²=H+1.

  Are these related? At H=3:
    (H-1)² = H+1       → 4 = 4     (self-consistency)
    H² - 3H + 1 = -1   → just barely past the critical point

  The self-consistency equation (H-1)²=H+1 gives H=3 or H=(3+√5)/2.
  Wait — let me check:
""")

# (H-1)² = H+1 → H²-2H+1 = H+1 → H²-3H = 0 → H(H-3) = 0
# Solutions: H=0 or H=3. NOT related to golden ratio.

# BUT the nucleation condition H²-3H+1 = 0 has roots at golden ratio.
# And the self-consistency H²-3H = 0 has root at H=3.

# The DIFFERENCE between the two equations:
# Self-consistency: H²-3H = 0 → H=3
# Nucleation critical: H²-3H+1 = 0 → H = (3+√5)/2 ≈ 2.618

# At H=3: H²-3H+1 = 9-9+1 = 1. Just above zero.
# The crystal at H=3 is ONE STEP past the critical nucleation point.

print(f"  Self-consistency: H²-3H = 0  → H = 3")
print(f"  Nucleation:       H²-3H+1 = 0 → H = (3+√5)/2 ≈ {(3+math.sqrt(5))/2:.4f}")
print(f"")
print(f"  At H=3: H²-3H+1 = {3**2 - 3*3 + 1}")
print(f"  At H=φ+1: H²-3H+1 = {(phi+1)**2 - 3*(phi+1) + 1:.10f}")
print(f"")
print(f"  The self-consistency condition (H²-3H=0) and the nucleation")
print(f"  condition (H²-3H+1=0) differ by exactly 1.")
print(f"")
print(f"  H=3 satisfies H²-3H = 0 (self-consistency) and gives")
print(f"  H²-3H+1 = 1 > 0 (nucleation before equilibrium).")
print(f"")
print(f"  The crystal at H=3 nucleates ONE UNIT past the critical point.")
print(f"  The golden ratio lives in the boundary between pre- and post-")
print(f"  equilibrium crystallization. H=3 is the first integer past it.")


# Connection: K_peak/K* at H=3
print(f"\n\n{'='*80}")
print("THE NUCLEATION RATIO AT H=3")
print("="*80)

H = 3
K_star = (H**2 - H + 1) / (H * (H**2 + 1))
K_peak = (H**2 + H + 1) / (2 * H * (H**2 + 1))
ratio = K_peak / K_star
print(f"\n  K_peak/K* = {ratio:.6f} = 13/14")
print(f"  1 - K_peak/K* = 1/14 = 1/(2(H²-H+1))")
print(f"  The crystal nucleates when it has accumulated 13/14 = 92.9% of K*.")
print(f"  The last 1/14 of the conflict journey (K_peak → K*) is the")
print(f"  crystallization regime where entropy FALLS.")
print(f"")
print(f"  The nucleation gap: K* - K_peak = 1/60 = 1/(2H(H²+1))")
print(f"  This is HALF the Born floor: 1/(2H³) = 1/54 ≠ 1/60")
print(f"  Actually: 1/60 = K*/(H²-H+1) = K*/7 × ... no:")
print(f"  K* - K_peak = K*/2 - 1/(H²+1) = ... let me compute:")

gap = K_star - K_peak
print(f"  K* - K_peak = {gap:.10f}")
print(f"  1/(2H(H²+1)) = {1/(2*H*(H**2+1)):.10f}")
print(f"  Match: {abs(gap - 1/(2*H*(H**2+1))) < 1e-12}")
print(f"")
print(f"  So K* - K_peak = 1/(2H(H²+1)) = K*/(2(H²-H+1))")
print(f"  At H=3: 1/60 = 7/(30×14) = K*/14")
