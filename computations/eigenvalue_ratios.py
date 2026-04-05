"""
Eigenvalue ratios of the crystal matrix C.

The dominant eigenvalue of C determines the composition fixed point.
The RATIO λ₁/λ₀ determines the convergence rate to the fixed point.

Preliminary observation: λ₁/λ₀ ≈ 0.618 ≈ 1/φ (golden ratio).
Is this real or seed-averaging artifact?

Check:
1. Per-seed eigenvalue ratios (not averaged)
2. Across all S₃ elements
3. For non-bijection crystals
4. At different numbers of entangler steps
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities)
from solver.crystals import Entangler

torch.set_grad_enabled(False)

PHI = (1 + math.sqrt(5)) / 2  # golden ratio ≈ 1.618
INV_PHI = 1 / PHI             # ≈ 0.618


S3_CORRS = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}


print("=" * 80)
print("EIGENVALUE RATIOS OF THE CRYSTAL MATRIX")
print(f"1/φ = {INV_PHI:.10f}")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Per-seed eigenvalue ratios for identity crystal
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Per-seed eigenvalue ratios (identity crystal) ---\n")

ratios_12 = []  # λ₁/λ₀
ratios_23 = []  # λ₂/λ₀
ratios_34 = []  # λ₃/λ₀

for seed in range(50):
    c = Entangler(S3_CORRS["e"], seed=seed).build().joint
    eigvals = torch.linalg.eigvals(c)
    mags = eigvals.abs()
    idx = mags.argsort(descending=True)
    sorted_mags = mags[idx]

    r12 = sorted_mags[1].item() / sorted_mags[0].item()
    r23 = sorted_mags[2].item() / sorted_mags[0].item()
    r34 = sorted_mags[3].item() / sorted_mags[0].item()

    ratios_12.append(r12)
    ratios_23.append(r23)
    ratios_34.append(r34)

    if seed < 5:
        print(f"  seed {seed}: λ = [{', '.join(f'{sorted_mags[i].item():.4f}' for i in range(4))}]"
              f"  λ₁/λ₀ = {r12:.4f}")

mean_12 = sum(ratios_12) / len(ratios_12)
std_12 = (sum((r - mean_12)**2 for r in ratios_12) / len(ratios_12))**0.5
mean_23 = sum(ratios_23) / len(ratios_23)
mean_34 = sum(ratios_34) / len(ratios_34)

print(f"\n  Over 50 seeds:")
print(f"    λ₁/λ₀ = {mean_12:.6f} ± {std_12:.6f}")
print(f"    λ₂/λ₀ = {mean_23:.6f}")
print(f"    λ₃/λ₀ = {mean_34:.6f}")
print(f"    1/φ    = {INV_PHI:.6f}")
print(f"    λ₁/λ₀ vs 1/φ: {abs(mean_12 - INV_PHI)/INV_PHI*100:.2f}% off")


# ═══════════════════════════════════════════════════════════════
#  Averaged crystal eigenvalue ratios for all S₃ elements
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- 50-seed averaged crystal eigenvalue ratios ---\n")

for name, corr in S3_CORRS.items():
    crystal = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(50):
        crystal += Entangler(corr, seed=seed).build().joint
    crystal /= 50

    eigvals = torch.linalg.eigvals(crystal)
    mags = eigvals.abs()
    idx = mags.argsort(descending=True)
    sorted_mags = mags[idx]

    r12 = sorted_mags[1].item() / sorted_mags[0].item()

    print(f"  {name:>6s}: λ = [{', '.join(f'{sorted_mags[i].item():.5f}' for i in range(4))}]"
          f"  λ₁/λ₀ = {r12:.6f}")

print(f"\n  1/φ = {INV_PHI:.6f}")


# ═══════════════════════════════════════════════════════════════
#  Non-bijection crystal eigenvalue ratios
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Non-bijection crystal eigenvalue ratios ---\n")

OTHER_CORRS = {
    "modular":      torch.tensor([[.4,.3,.3],[.3,.4,.3],[.3,.3,.4]], dtype=torch.float32),
    "exponential":  torch.tensor([[1,0,0],[.6,.4,0],[0,0,1]], dtype=torch.float32),
    "quadratic":    torch.tensor([[1,0,0],[0,.3,.7],[0,0,1]], dtype=torch.float32),
    "logarithmic":  torch.tensor([[.8,.2,0],[0,.6,.4],[0,.1,.9]], dtype=torch.float32),
    "constant_0":   torch.tensor([[1,0,0],[1,0,0],[1,0,0]], dtype=torch.float32),
    "two_to_one":   torch.tensor([[1,0,0],[1,0,0],[0,1,0]], dtype=torch.float32),
    "uniform":      torch.tensor([[1,1,1],[1,1,1],[1,1,1]], dtype=torch.float32) / 3,
}

for name, corr in OTHER_CORRS.items():
    crystal = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(50):
        crystal += Entangler(corr, seed=seed).build().joint
    crystal /= 50

    eigvals = torch.linalg.eigvals(crystal)
    mags = eigvals.abs()
    idx = mags.argsort(descending=True)
    sorted_mags = mags[idx]

    r12 = sorted_mags[1].item() / sorted_mags[0].item()

    print(f"  {name:>14s}: λ = [{', '.join(f'{sorted_mags[i].item():.5f}' for i in range(4))}]"
          f"  λ₁/λ₀ = {r12:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Eigenvalue ratios at different entangler steps
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Eigenvalue ratio vs entangler steps ---\n")

for n_steps in [10, 20, 30, 40, 50, 80, 100]:
    crystal = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(50):
        crystal += Entangler(S3_CORRS["e"], seed=seed).build(n_steps=n_steps).joint
    crystal /= 50

    eigvals = torch.linalg.eigvals(crystal)
    mags = eigvals.abs()
    idx = mags.argsort(descending=True)
    sorted_mags = mags[idx]

    r12 = sorted_mags[1].item() / sorted_mags[0].item()
    print(f"  {n_steps:3d} steps: λ₁/λ₀ = {r12:.6f}")

print(f"\n  1/φ = {INV_PHI:.6f}")


# ═══════════════════════════════════════════════════════════════
#  Eigenvalue ratios at different discount factors
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Eigenvalue ratio vs discount factor ---\n")

for discount in [0.1, 0.2, 0.3, 0.4, 0.5]:
    crystal = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(50):
        crystal += Entangler(S3_CORRS["e"], seed=seed).build(discount=discount).joint
    crystal /= 50

    eigvals = torch.linalg.eigvals(crystal)
    mags = eigvals.abs()
    idx = mags.argsort(descending=True)
    sorted_mags = mags[idx]

    r12 = sorted_mags[1].item() / sorted_mags[0].item()
    print(f"  d={discount:.1f}: λ₁/λ₀ = {r12:.6f}, λ₀ = {sorted_mags[0].item():.6f}")


# ═══════════════════════════════════════════════════════════════
#  Is the ratio related to the eigenvalues of the CORRELATION?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Eigenvalues of the correlation matrix ---\n")

for name in ["e", "modular", "exponential"]:
    if name in S3_CORRS:
        corr = S3_CORRS[name]
    else:
        corr = OTHER_CORRS[name]
    ev = torch.linalg.eigvals(corr.to(torch.cfloat))
    mags = ev.abs()
    idx = mags.argsort(descending=True)
    print(f"  {name:>14s}: corr eigenvalues = [{', '.join(f'{mags[idx[i]].item():.4f}' for i in range(H))}]")


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
