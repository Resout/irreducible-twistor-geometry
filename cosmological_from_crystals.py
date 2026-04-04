"""
Cosmological constant from crystal constants.

The cosmological constant Λ in the paper arises from the conformal-to-Einstein
reduction: Bach-flat conformal gravity → Einstein gravity + Λ.

Maldacena showed: Λ = (scalar sector fraction) × (curvature of the moduli space).

In the crystal framework:
  Scalar fraction = 1 - K* = 23/30 (Principle 43)
  Curvature = κ ≈ 1.3 (Principle 22, positive Gaussian curvature)
  Born floor = 1/27 (compactness of the mass manifold)
  Fixed-point ignorance = 10/33 (Principle 48)
  Fisher-Rao = 2 × Fubini-Study (Principle 39)

The question: can we compute a dimensionless ratio that equals Λ/m²
where m is the mass gap?

Λ/m² should be a pure number built from H=3 constants.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities, born_fidelity,
                            sym2_fingerprint, sym2_distance)
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)


print("=" * 80)
print("COSMOLOGICAL CONSTANT FROM CRYSTAL CONSTANTS")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  All the constants
# ═══════════════════════════════════════════════════════════════

print(f"\n--- The crystal constants (all from H=3) ---\n")

constants = {
    'H': H,
    'K*': K_STAR,
    'Δ': DELTA,
    '1-K*': 1-K_STAR,
    'Born_floor': BORN_FLOOR,
    'p*': math.exp(DELTA/2),
    'Born(θ)_fp': 10/33,
    'Schmidt_peak': 3.272,  # measured
    'composition_rate': 299/405,
    '2Δ': 2*DELTA,
}

for name, val in constants.items():
    print(f"  {name:>20s} = {val:.6f}")


# ═══════════════════════════════════════════════════════════════
#  Candidate Λ/m² ratios
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Candidate dimensionless ratios for Λ/m² ---\n")

# The mass gap m² = Δ² (in natural units where the lattice spacing = 1)
# The cosmological constant Λ has dimensions of 1/length²

# Candidates:
candidates = {
    # From the scalar/graviton fractions
    '(1-K*) × K*': (1-K_STAR) * K_STAR,
    '(1-K*) × Born_floor': (1-K_STAR) * BORN_FLOOR,
    '(1-K*) × Δ': (1-K_STAR) * DELTA,
    '(1-K*) / H²': (1-K_STAR) / H**2,
    'K* × Born_floor': K_STAR * BORN_FLOOR,
    'K* / (H(H²+1))': K_STAR / (H*(H**2+1)),

    # From the fixed-point chain
    'Born(θ)_fp × K*': (10/33) * K_STAR,
    'Born(θ)_fp - Born_floor': 10/33 - BORN_FLOOR,
    'Born(θ)_fp × (1-K*)': (10/33) * (1-K_STAR),
    '(Born(θ)_fp)²': (10/33)**2,

    # From composition
    '1 - 299/405': 1 - 299/405,
    'Δ²': DELTA**2,
    'Δ × Born_floor': DELTA * BORN_FLOOR,
    'Δ / (4π)': DELTA / (4*math.pi),
    'Δ² / (4π)': DELTA**2 / (4*math.pi),

    # Curvature-based
    '4Δ² / H²': 4 * DELTA**2 / H**2,
    '(1-K*)²/H³': (1-K_STAR)**2 / H**3,
    'K*/H³': K_STAR / H**3,

    # The Penrose integral normalization
    # In twistor theory, Λ = 3/(twistor radius)²
    '3/Schmidt_peak²': 3 / 3.272**2,
    '3 × Born_floor': 3 * BORN_FLOOR,
    'H × Born_floor': H * BORN_FLOOR,

    # Simple fractions
    '7/(30×27)': 7/(30*27),
    '1/(H²(H²+1))': 1/(H**2 * (H**2+1)),
    '1/(H(H²+1))': 1/(H * (H**2+1)),
    '(H²-H+1)/(H³(H²+1))': (H**2-H+1)/(H**3*(H**2+1)),
}

# Sort by value
for name, val in sorted(candidates.items(), key=lambda x: x[1]):
    # Check if it's a nice fraction
    frac = ""
    for num in range(1, 100):
        for den in range(1, 1000):
            if abs(num/den - val) < 1e-6:
                frac = f"= {num}/{den}"
                break
        if frac:
            break
    print(f"  {name:>30s} = {val:.8f}  {frac}")


# ═══════════════════════════════════════════════════════════════
#  The actual curvature of crystal space
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CURVATURE OF CRYSTAL SPACE (SECTIONAL CURVATURE)")
print("="*80)

# Build a grid of crystals and measure curvature
# Using Sym² fingerprints in 9D space

n_seeds = 30

# Reference crystals
S3_CORRS = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

# Build all crystals with seed averaging
crystals = {}
for name, corr in S3_CORRS.items():
    s = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        s += Entangler(corr, seed=seed).build().joint
    crystals[name] = s / n_seeds

# Sym² fingerprints
fps = {n: sym2_fingerprint(c) for n, c in crystals.items()}

# Fubini-Study distances
print(f"\n  Fubini-Study distances between S₃ crystals:")
names = list(S3_CORRS.keys())
for i in range(len(names)):
    for j in range(i+1, len(names)):
        a, b = names[i], names[j]
        # FS distance
        z1, z2 = crystals[a].reshape(-1), crystals[b].reshape(-1)
        inner = (z1.conj() @ z2).abs().item()
        n1 = z1.abs().pow(2).sum().sqrt().item()
        n2 = z2.abs().pow(2).sum().sqrt().item()
        cos_d = min(inner / (n1 * n2), 1.0)
        fs = math.acos(cos_d)

        # Sym² distance
        sd = sym2_distance(fps[a], fps[b])

        print(f"    {a:>6s} ↔ {b:>6s}: FS={fs:.4f}, Sym²={sd:.6f}")

# Mean FS distance
fs_vals = []
for i in range(len(names)):
    for j in range(i+1, len(names)):
        z1, z2 = crystals[names[i]].reshape(-1), crystals[names[j]].reshape(-1)
        inner = (z1.conj() @ z2).abs().item()
        n1 = z1.abs().pow(2).sum().sqrt().item()
        n2 = z2.abs().pow(2).sum().sqrt().item()
        cos_d = min(inner / (n1 * n2), 1.0)
        fs_vals.append(math.acos(cos_d))

mean_fs = sum(fs_vals) / len(fs_vals)
print(f"\n  Mean FS distance: {mean_fs:.4f}")
print(f"  FS diameter (max): {max(fs_vals):.4f}")
print(f"  FS diameter / π: {max(fs_vals)/math.pi:.4f}")


# ═══════════════════════════════════════════════════════════════
#  The scalar/graviton connection to Λ
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SCALAR/GRAVITON RATIO AND Λ")
print("="*80)

scalar_frac = 1 - K_STAR  # 23/30
graviton_frac = K_STAR     # 7/30 (conflicting fraction)

print(f"\n  Scalar fraction (1-K*) = {scalar_frac:.6f} = 23/30")
print(f"  Graviton fraction K* = {graviton_frac:.6f} = 7/30")
print(f"  Ratio scalar/graviton = {scalar_frac/graviton_frac:.6f} = 23/7 = {23/7:.6f}")
print(f"  23/7 = {23/7:.10f}")
print(f"  This is NOT π (3.14159...) but close: {23/7 - math.pi:.6f} difference")

# In conformal gravity: Λ = (Weyl²)_scalar / (Weyl²)_total × R²
# where R is the curvature radius
# The scalar fraction = surviving fraction after DS equilibrium
# The graviton fraction = conflicting fraction

# What if Λ/m² = K* × (1-K*)?
# = 7/30 × 23/30 = 161/900
lambda_candidate = K_STAR * (1-K_STAR)
print(f"\n  K* × (1-K*) = {lambda_candidate:.6f} = 161/900")
print(f"  = {161/900:.6f}")

# What about Λ/m² = K*² / (1-K*)?
# = 49/900 / (23/30) = 49/900 × 30/23 = 49/690
lambda2 = K_STAR**2 / (1-K_STAR)
print(f"  K*² / (1-K*) = {lambda2:.6f} = 49/690 = {49/690:.6f}")

# What about using the composition rate?
# Λ/m² = 1 - (composition rate) = 1 - 299/405 = 106/405
lambda3 = 1 - 299/405
print(f"  1 - 299/405 = {lambda3:.6f} = 106/405 = {106/405:.6f}")

# Check: is 106/405 = something nice?
# 106 = 2 × 53
# 405 = 5 × 81 = 5 × 3⁴
print(f"  106/405 = 2×53 / (5×3⁴) — not obviously simple")

# The difference between the entangling rate (Δ) and the composition rate (2Δ)
# gives the NET information gain per cycle:
# In one cycle: gain Δ/4 from entangling, lose 2Δ from composition
# Net: -7Δ/4 per cycle
# Λ might be related to this net loss
net_loss = 7 * DELTA / 4
print(f"\n  Net loss per cycle = 7Δ/4 = {net_loss:.6f}")
print(f"  7Δ/4 / (4π) = {net_loss/(4*math.pi):.6f}")

# Actually, the most natural candidate:
# Λ = Δ² because Δ sets the scale and Λ has dimensions of 1/length²
print(f"\n  Δ² = {DELTA**2:.6f}")
print(f"  Δ² × (1-K*) = {DELTA**2 * (1-K_STAR):.6f}")
print(f"  Δ² × K* = {DELTA**2 * K_STAR:.6f}")
print(f"  Δ² / (4π) = {DELTA**2 / (4*math.pi):.6f}")
print(f"  Δ² / 4 = {DELTA**2 / 4:.6f}")
print(f"  Δ² × Born_floor = {DELTA**2 * BORN_FLOOR:.6f}")


# ═══════════════════════════════════════════════════════════════
#  Volume of crystal space
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("VOLUME OF CRYSTAL SPACE")
print("="*80)

# The crystal space is CP³ with Born floor S⁵.
# Volume of CP³ = π³/6
# But the Born floor restricts to a compact subset.
# The effective volume might be the relevant normalization.

vol_CP3 = math.pi**3 / 6
print(f"\n  Vol(CP³) = π³/6 = {vol_CP3:.6f}")
print(f"  Vol(S⁵) = π³ = {math.pi**3:.6f}")
print(f"  Born(θ) = 1/27 restricts to a sphere of radius...")

# The Born floor defines |m_h|² / (|m_h|² + |m_θ|²) ≥ 1 - 1/27 = 26/27
# So |m_θ|² / |m_h|² ≤ 1/26
# This is a constraint on the ratio, defining a cap on S⁵

# Effective volume fraction ≈ Born_floor (fraction of probability space
# that's "locked" into ignorance)
print(f"  Effective volume fraction ≈ Born_floor = {BORN_FLOOR:.6f}")
print(f"  Λ ∝ 1/V_eff ∝ 1/Born_floor = {1/BORN_FLOOR:.1f} = H³ = {H**3}")


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
print(f"""
  The cosmological constant computation requires identifying which
  dimensionless ratio equals Λ/m². The candidates are all built from
  K* = 7/30 and H = 3 (zero free parameters).

  The most natural candidates:
    K*(1-K*) = 161/900 ≈ 0.179  (scalar × graviton fraction)
    Δ²       = 0.0706           (mass gap squared)
    1-299/405 = 106/405 ≈ 0.262 (composition deficit)

  What's missing: the geometric normalization from the Penrose integral.
  In twistor theory, Λ = 3/(twistor radius)². The twistor radius in
  the crystal is the FS diameter of the S₃ orbit, which needs to be
  related to the paper's continuous twistor space.
""")
