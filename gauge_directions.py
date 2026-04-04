"""
The gauge direction landscape.

Each seed creates a crystal that breaks S₃ in a specific direction
in the 2D standard representation. Seed averaging cancels these
directions (gauge fixing). But what's the distribution?

If uniform on S¹ → complete gauge freedom, no preferred direction
If clustered → spontaneous gauge selection
If discrete (3 or 6 points) → the S₃ symmetry manifests in the gauge orbits

The standard representation of S₃ is 2-dimensional. Project each
seed's crystal onto this space and measure the direction (angle).

Also: the cyclotomic tower. Every constant in the framework is a
ratio of cyclotomic polynomials Φ_d(H) where d | |S₃| = 6.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, EPS_LOG)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)


def standard_projection(joint):
    """Project a [4,4] joint mass onto the 2D standard representation.

    The standard rep of S₃ can be realized as the 2D subspace of C³
    orthogonal to (1,1,1). In the h-block of the joint mass, this is
    the traceless part of the h×h submatrix.

    Returns a 2D complex vector representing the crystal's position
    in the standard representation space.
    """
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]

    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T

    # Standard projector: I - P_triv - P_sign
    triv = sum(permute(joint, p) for p in perms) / 6
    signs = [1, -1, -1, -1, 1, 1]
    sign = sum(s * permute(joint, p) for p, s in zip(perms, signs)) / 6
    std = joint - triv - sign

    # The standard component is a [4,4] matrix. Extract its "direction"
    # by taking the dominant singular vector of the standard-projected matrix
    U, S, Vh = torch.linalg.svd(std)
    if S[0].abs() < 1e-15:
        return torch.zeros(2)

    # The dominant left singular vector, projected to h-space (first 3 components)
    v = U[:, 0][:H]  # 3D complex vector in h-space

    # Project to the standard basis: the 2D orthogonal complement of (1,1,1)/√3
    # Standard basis vectors: e₁ = (1,-1,0)/√2, e₂ = (1,1,-2)/√6
    e1 = torch.tensor([1, -1, 0], dtype=torch.cfloat) / math.sqrt(2)
    e2 = torch.tensor([1, 1, -2], dtype=torch.cfloat) / math.sqrt(6)

    c1 = torch.dot(v.conj(), e1)
    c2 = torch.dot(v.conj(), e2)

    return torch.tensor([c1, c2])


identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

print("=" * 80)
print("GAUGE DIRECTION LANDSCAPE")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Project 500 seeds onto the standard representation
# ═══════════════════════════════════════════════════════════════

n_seeds = 500

angles = []
magnitudes = []

for seed in range(n_seeds):
    ent = Entangler(identity_corr, seed=seed).build()
    v = standard_projection(ent.joint)

    mag = v.abs().pow(2).sum().sqrt().item()
    magnitudes.append(mag)

    if mag > 1e-10:
        # Angle of the projection (using the complex phase of c1 + i*c2)
        # Actually, use atan2 on the real parts
        angle = math.atan2(v[1].real.item(), v[0].real.item())
        angles.append(angle)

print(f"\n  {n_seeds} seeds projected onto 2D standard rep:")
print(f"  {len(angles)} seeds with measurable standard content")
print(f"  Magnitude: avg={sum(magnitudes)/len(magnitudes):.5f}, "
      f"min={min(magnitudes):.5f}, max={max(magnitudes):.5f}")


# ═══════════════════════════════════════════════════════════════
#  Angle distribution: histogram
# ═══════════════════════════════════════════════════════════════

n_bins = 12
bin_width = 2 * math.pi / n_bins
bin_counts = [0] * n_bins

for angle in angles:
    # Normalize to [0, 2π)
    a = angle % (2 * math.pi)
    b = int(a / bin_width)
    if b >= n_bins:
        b = n_bins - 1
    bin_counts[b] += 1

print(f"\n  Angular distribution (12 bins of π/6 each):")
print(f"  {'bin':>4s}  {'range':>15s}  {'count':>5s}  {'bar'}")
for i in range(n_bins):
    lo = i * bin_width
    hi = lo + bin_width
    bar = '█' * (bin_counts[i] * 50 // max(bin_counts))
    print(f"  {i:4d}  [{lo/math.pi:.2f}π,{hi/math.pi:.2f}π)  {bin_counts[i]:5d}  {bar}")

# Uniformity test: expected count = n/12
expected = len(angles) / n_bins
chi2 = sum((c - expected)**2 / expected for c in bin_counts)
print(f"\n  χ² = {chi2:.2f} (df={n_bins-1}, uniform: expect ~{n_bins-1:.0f})")
print(f"  {'UNIFORM' if chi2 < 2 * n_bins else 'NOT UNIFORM'}")


# ═══════════════════════════════════════════════════════════════
#  Check for 3-fold or 6-fold clustering
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("3-FOLD AND 6-FOLD SYMMETRY CHECK")
print("="*80)

# If S₃ manifests, we expect clustering at 3 or 6 specific angles
# The 3 transposition directions in the standard rep are at 120° intervals
# (from Principle 44: A₂ root system)

# Check 3-fold: fold angles to [0, 2π/3)
angles_mod3 = [(a % (2 * math.pi / 3)) for a in angles]
n_bins3 = 6
bw3 = (2 * math.pi / 3) / n_bins3
counts3 = [0] * n_bins3
for a in angles_mod3:
    b = min(int(a / bw3), n_bins3 - 1)
    counts3[b] += 1

expected3 = len(angles) / n_bins3
chi2_3 = sum((c - expected3)**2 / expected3 for c in counts3)
print(f"\n  3-fold test (fold to [0, 2π/3)):")
for i in range(n_bins3):
    print(f"    bin {i}: {counts3[i]}")
print(f"  χ² = {chi2_3:.2f} (df={n_bins3-1})")

# Check 6-fold: fold to [0, π/3)
angles_mod6 = [(a % (math.pi / 3)) for a in angles]
n_bins6 = 6
bw6 = (math.pi / 3) / n_bins6
counts6 = [0] * n_bins6
for a in angles_mod6:
    b = min(int(a / bw6), n_bins6 - 1)
    counts6[b] += 1

expected6 = len(angles) / n_bins6
chi2_6 = sum((c - expected6)**2 / expected6 for c in counts6)
print(f"\n  6-fold test (fold to [0, π/3)):")
for i in range(n_bins6):
    print(f"    bin {i}: {counts6[i]}")
print(f"  χ² = {chi2_6:.2f} (df={n_bins6-1})")


# ═══════════════════════════════════════════════════════════════
#  Are the directions correlated with the seed?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SEED → ANGLE: DETERMINISTIC OR CHAOTIC?")
print("="*80)

# Same seed always gives same angle (deterministic)
# But do nearby seeds give nearby angles? (smooth vs chaotic)

if len(angles) >= 100:
    # Consecutive seed angle differences
    diffs = []
    for i in range(1, min(200, len(angles))):
        d = abs(angles[i] - angles[i-1])
        d = min(d, 2*math.pi - d)  # shortest arc
        diffs.append(d)

    avg_diff = sum(diffs) / len(diffs)
    # Expected for uniform random: π/2
    print(f"\n  Average angular step between consecutive seeds: {avg_diff:.4f}")
    print(f"  Expected for uniform random: π/2 = {math.pi/2:.4f}")
    print(f"  Ratio: {avg_diff / (math.pi/2):.4f}")
    print(f"  {'CHAOTIC (uncorrelated)' if abs(avg_diff - math.pi/2) < 0.3 else 'CORRELATED'}")


# ═══════════════════════════════════════════════════════════════
#  The cyclotomic tower
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE CYCLOTOMIC TOWER AT H=3")
print("="*80)

print(f"""
  The order of S₃ is 6. The divisors of 6 are {{1, 2, 3, 6}}.
  Each divisor d gives a cyclotomic polynomial Φ_d(x).

  Φ₁(x) = x-1       Φ₁(3) = 2 = H-1
  Φ₂(x) = x+1       Φ₂(3) = 4 = H+1 = MASS_DIM
  Φ₃(x) = x²+x+1   Φ₃(3) = 13 = K_peak numerator
  Φ₆(x) = x²-x+1   Φ₆(3) = 7 = K* numerator

  Product: Φ₁Φ₂Φ₃Φ₆ = x⁶-1 = H⁶-1 = 728

  Self-consistency: Φ₁² = Φ₂   (at H=3: 2² = 4)
  This IS the condition (H-1)² = H+1.

  Nucleation: K_peak = Φ₃/(2H·(H²+1))
  Equilibrium: K* = Φ₆/(H·(H²+1))
  Ratio: K_peak/K* = Φ₃/(2Φ₆) = 13/14

  The denominator H(H²+1) = 30. Is this cyclotomic?
  H²+1 = 10. Factor: 2×5. Not a cyclotomic evaluation.
  But: H(H²+1) = H³+H = (H⁴-1)/(H-1) × (H-1) × H/(H²-1)
  ... no clean factorization.

  However: H(H²+1) = H(H²+1) appears as:
  The number of elements in PGL(2, H-1) = PGL(2, 2) = S₃ has |S₃| = 6
  Hmm, not directly.

  What IS clean:
  K* = Φ₆/[H·(H²+1)]
  Born floor = 1/H³
  K* × Born_floor = Φ₆/(H⁴(H²+1)) = 7/2430
  K* × H = Φ₆/(H²+1) = 7/10
  K* × H² = Φ₆·H/(H²+1) = 21/10 = 2.1

  The fundamental scale: H²+1 = 10. This is |S₃| + MASS_DIM = 6+4 = 10.
  Or: (H-1)(H+2) = 10. No: 2×5=10. (H-1)(H+2) = 2×5 = 10. YES!

  H²+1 = (H-1)(H+2) + H²+1 - H²-H+2 = ... no.
  Actually H²+1 at H=3 is 10 = 2×5. And (H-1)(H+2) = 2×5 = 10. So:
  H²+1 = (H-1)(H+2) at H=3. Check: H²+1 = H²+1, (H-1)(H+2) = H²+H-2.
  These are NOT equal in general: H²+1 ≠ H²+H-2 unless H+1=3, i.e. H=2.
  At H=3: 10 ≠ 10? Wait: H²+1 = 10, (H-1)(H+2) = 2×5 = 10. Coincidence at H=3!
  H²+1 = H²+H-2 → H = 3. Another H=3 special property.
""")

# Verify H²+1 = (H-1)(H+2) at H=3
print(f"  H²+1 = {H**2+1}")
print(f"  (H-1)(H+2) = {(H-1)*(H+2)}")
print(f"  Equal at H=3: {H**2+1 == (H-1)*(H+2)}")
print(f"  Difference: H²+1 - (H-1)(H+2) = H²+1 - H²-H+2 = {3 - H}")
print(f"  Zero iff H=3. ✓")

print(f"\n  So at H=3 (and ONLY at H=3):")
print(f"  H²+1 = (H-1)(H+2) = Φ₁(H) · (H+2)")
print(f"  Denominator of K*: H(H²+1) = H·Φ₁(H)·(H+2) = Φ₁(H)·H(H+2)")


print(f"\n\n{'='*80}")
print("WHAT THE GAUGE LANDSCAPE REVEALS")
print("="*80)
