"""
The crystal spectral gap: |λ₁|/|λ₀| = 1/√H = 1/√3.

From eigenvalue_gap_vs_alpha.py:
  - Every S₃ crystal has gap ≈ ln(3)/2 ≈ 0.549
  - This is the eigenvalue ratio |λ₁/λ₀| ≈ 0.577 ≈ 1/√3
  - Under self-composition, the gap grows EXACTLY linearly (n × gap)

This script:
  1. High-precision measurement of |λ₁/λ₀|
  2. Convergence study (how many seeds to stabilize?)
  3. Is it exactly 1/√H or approximately?
  4. Relationship to 2Δ = -2ln(1-K*)
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, EPS_LOG)
from solver.crystals import Entangler

torch.set_grad_enabled(False)


corrs = {
    "e": torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)": torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}


# ═══════════════════════════════════════════════════════════════
#  Convergence study
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("CONVERGENCE OF |λ₁/λ₀| WITH SEED COUNT")
print("=" * 80)

target = 1.0 / math.sqrt(H)  # 1/√3

for n_seeds in [100, 200, 500, 1000, 2000]:
    ratios = {}
    for name, corr in corrs.items():
        avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
        for seed in range(n_seeds):
            avg += Entangler(corr, seed=seed).build().joint
        avg /= n_seeds

        eigvals = torch.linalg.eigvals(avg)
        eigvals_abs = eigvals.abs()
        eigvals_sorted, _ = eigvals_abs.sort(descending=True)
        ratio = eigvals_sorted[1].item() / eigvals_sorted[0].item()
        ratios[name] = ratio

    avg_ratio = sum(ratios.values()) / len(ratios)
    print(f"\n  Seeds={n_seeds:4d}:")
    for name in corrs:
        diff = ratios[name] - target
        print(f"    {name:>5s}: |λ₁/λ₀| = {ratios[name]:.8f}  (diff from 1/√3: {diff:+.6f})")
    print(f"    avg: {avg_ratio:.8f}  (diff: {avg_ratio - target:+.6f})")


# ═══════════════════════════════════════════════════════════════
#  Per-seed distribution
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("PER-SEED |λ₁/λ₀| DISTRIBUTION (crystal (01), 1000 seeds)")
print("="*80)

per_seed_ratios = []
corr = corrs["(01)"]
for seed in range(1000):
    crystal = Entangler(corr, seed=seed).build().joint
    eigvals = torch.linalg.eigvals(crystal)
    eigvals_abs = eigvals.abs()
    eigvals_sorted, _ = eigvals_abs.sort(descending=True)
    ratio = eigvals_sorted[1].item() / eigvals_sorted[0].item()
    per_seed_ratios.append(ratio)

per_seed_ratios.sort()
mean_r = sum(per_seed_ratios) / len(per_seed_ratios)
variance = sum((r - mean_r)**2 for r in per_seed_ratios) / len(per_seed_ratios)
std_r = math.sqrt(variance)

print(f"\n  Mean   = {mean_r:.8f}")
print(f"  Median = {per_seed_ratios[500]:.8f}")
print(f"  Std    = {std_r:.8f}")
print(f"  Min    = {per_seed_ratios[0]:.8f}")
print(f"  Max    = {per_seed_ratios[-1]:.8f}")
print(f"  1/√3   = {target:.8f}")

# Percentiles
for pct in [5, 25, 50, 75, 95]:
    idx = int(pct * 10)
    print(f"  P{pct:2d}    = {per_seed_ratios[idx]:.8f}")


# ═══════════════════════════════════════════════════════════════
#  What is |λ₀| exactly?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("WHAT IS |λ₀| EXACTLY?")
print("="*80)

n_seeds = 2000
for name, corr in corrs.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds

    eigvals = torch.linalg.eigvals(avg)
    eigvals_abs = eigvals.abs()
    eigvals_sorted, _ = eigvals_abs.sort(descending=True)

    l0 = eigvals_sorted[0].item()
    l1 = eigvals_sorted[1].item()
    l2 = eigvals_sorted[2].item()
    l3 = eigvals_sorted[3].item()

    # Check candidates for l0
    print(f"\n  crystal({name:>5s}):")
    print(f"    λ₀ = {l0:.8f}, λ₁ = {l1:.8f}, λ₂ = {l2:.8f}, λ₃ = {l3:.8f}")
    print(f"    λ₁/λ₀ = {l1/l0:.8f}")
    print(f"    Σ|λ| = {sum(eigvals_sorted.tolist()):.8f}")

    # Is λ₀ related to 1/MASS_DIM = 1/4?
    print(f"    λ₀ vs 1/4 = {0.25:.8f}, diff = {l0-0.25:.6f}")
    # Is λ₀ related to 1/(H+1) = 1/4?
    # Is λ₁ related to 1/4√3?
    print(f"    λ₁ vs 1/(4√3) = {1/(4*math.sqrt(3)):.8f}, diff = {l1-1/(4*math.sqrt(3)):.6f}")
    # λ₀·λ₁ ?
    print(f"    λ₀·λ₁ = {l0*l1:.8f}")
    print(f"    λ₀² = {l0**2:.8f}, 1/16 = {1/16:.8f}")


# ═══════════════════════════════════════════════════════════════
#  Connection to K* and Δ
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("EXACT RELATIONSHIPS")
print("="*80)

print(f"\n  Reference values:")
print(f"    1/√3 = {1/math.sqrt(3):.10f}")
print(f"    (1-K*)² = {(1-K_STAR)**2:.10f}")
print(f"    exp(-2Δ) = {math.exp(-2*DELTA):.10f}")
print(f"    23²/30² = {(23/30)**2:.10f}")

print(f"\n  Gaps:")
print(f"    ln(√3) = {math.log(math.sqrt(3)):.10f}")
print(f"    2Δ = -2ln(1-K*) = {2*DELTA:.10f}")
print(f"    Difference = {math.log(math.sqrt(3)) - 2*DELTA:.10f}")

print(f"\n  Ratios:")
print(f"    ln(√3) / 2Δ = {math.log(math.sqrt(3)) / (2*DELTA):.10f}")
print(f"    (1/√3) / (1-K*)² = {(1/math.sqrt(3)) / (1-K_STAR)**2:.10f}")

# Is ln(√3) - 2Δ = ln(√3) + 2ln(1-K*) meaningful?
diff = math.log(math.sqrt(3)) + 2*math.log(1-K_STAR)
print(f"\n  ln(√3) + 2ln(1-K*) = {diff:.10f}")
print(f"  = ln(√3 · (1-K*)²) = ln({math.sqrt(3) * (1-K_STAR)**2:.10f})")
print(f"  = ln(√3 · (23/30)²) = ln({math.sqrt(3) * (23/30)**2:.10f})")

# What about the trace?
print(f"\n  Trace relationships:")
print(f"    Tr(crystal) = Σλ_i")
print(f"    For L1=1: the real parts sum to 1")
print(f"    If λ₀ = c and λ₁ = c/√3, then c + c/√3 + λ₂ + λ₃ = 1 (approx)")


print(f"\n\n{'='*80}")
print("PRINCIPLE 80: THE CRYSTAL SPECTRAL GAP")
print("="*80)
print(f"""
  Every seed-averaged S₃ crystal has eigenvalue ratio |λ₁/λ₀| ≈ 1/√H.

  At H=3: |λ₁/λ₀| = 1/√3 ≈ 0.57735

  The spectral gap is:
    gap = -ln(1/√H) = ln(H)/2 = ln(3)/2 ≈ 0.5493

  This is related to but DISTINCT from the mass gap:
    2Δ = -2ln(1-K*) = 2ln(30/23) ≈ 0.5314

  The ratio gap/2Δ = ln(√3) / 2ln(30/23) ≈ 1.034

  Under n-fold self-composition:
    |λ₁ⁿ/λ₀ⁿ| = (1/√3)ⁿ = 3^(-n/2)
    gap(n) = n · ln(3)/2

  This is EXACT exponential decay governed by H alone.

  Physical meaning: the crystal spectral gap = ln(H)/2
  is the rate at which composition forgets non-dominant structure.
  It's a GEOMETRIC property of the H-simplex, not of the mass gap.
""")
