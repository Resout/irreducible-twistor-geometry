"""
Does the composition spectral gap scale with n_steps?

If the gap comes from the bilinear structure of the entangler
(each step contributes Δ of contraction, 50 steps → 50Δ?), then
the gap should scale linearly with n_steps.

But the gap is ≈ 2Δ regardless of n_steps (at optimal discount).
This suggests the gap is a PROPERTY of the crystal at equilibrium,
not proportional to the number of construction steps.

Test: measure the gap at different n_steps (adjusting discount to
maintain the optimal total evidence budget n_steps × discount ≈ 20).
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, schmidt_number
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)

corr = torch.eye(H, dtype=torch.float32)


def eigendata(joint):
    ev, _ = torch.linalg.eig(joint)
    idx = ev.abs().argsort(descending=True)
    return ev[idx]


def composition_gap(joint):
    """Spectral gap of joint-as-matrix = -ln(|λ₁/λ₀|)."""
    ev = eigendata(joint)
    if ev[0].abs().item() > 1e-12 and ev[1].abs().item() > 1e-12:
        return math.log(ev[0].abs().item() / ev[1].abs().item())
    return 0


# ═══════════════════════════════════════════════════════════════
#  Gap vs n_steps at fixed discount (d=0.3)
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("SPECTRAL GAP vs n_steps (discount=0.3)")
print(f"Δ = {DELTA:.4f}, 2Δ = {2*DELTA:.4f}")
print("=" * 80)

print(f"\n  {'n_steps':>8s} {'gap':>8s} {'gap/Δ':>8s} {'Schmidt':>8s} {'|λ₀|':>8s} {'|λ₁|':>8s}")
print("  " + "-" * 55)

for n_steps in [5, 10, 15, 20, 30, 40, 50, 75, 100, 150, 200]:
    gaps = []
    schmidts = []
    for seed in range(20):
        ent = Entangler(corr, seed=seed)
        ent.build(n_steps=n_steps, discount=0.3)
        g = composition_gap(ent.joint)
        s = schmidt_number(ent.joint)
        gaps.append(g)
        schmidts.append(s)

    avg_gap = sum(gaps) / len(gaps)
    avg_sc = sum(schmidts) / len(schmidts)

    # Get eigenvalues for one seed to show structure
    ent = Entangler(corr, seed=42)
    ent.build(n_steps=n_steps, discount=0.3)
    ev = eigendata(ent.joint)

    print(f"  {n_steps:>8d} {avg_gap:>8.4f} {avg_gap/DELTA:>8.3f} {avg_sc:>8.3f} "
          f"{ev[0].abs().item():>8.4f} {ev[1].abs().item():>8.4f}")


# ═══════════════════════════════════════════════════════════════
#  Gap vs n_steps at constant total evidence (n*d ≈ 15)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SPECTRAL GAP vs n_steps at constant total evidence (n*d ≈ 15)")
print("="*80)

total_evidence = 15.0

print(f"\n  {'n_steps':>8s} {'discount':>10s} {'gap':>8s} {'gap/Δ':>8s} {'Schmidt':>8s}")
print("  " + "-" * 50)

for n_steps in [10, 15, 20, 30, 50, 75, 100, 150]:
    discount = total_evidence / n_steps
    if discount > 0.8:
        continue  # too unstable

    gaps = []
    for seed in range(20):
        try:
            ent = Entangler(corr, seed=seed)
            ent.build(n_steps=n_steps, discount=discount)
            g = composition_gap(ent.joint)
            s = schmidt_number(ent.joint)
            if not math.isnan(g) and abs(g) < 100:
                gaps.append(g)
        except:
            pass

    if gaps:
        avg_gap = sum(gaps) / len(gaps)
        print(f"  {n_steps:>8d} {discount:>10.4f} {avg_gap:>8.4f} {avg_gap/DELTA:>8.3f} "
              f"{sum(schmidts)/len(schmidts) if schmidts else 0:>8.3f}")


# ═══════════════════════════════════════════════════════════════
#  The eigenvalue structure vs n_steps
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("EIGENVALUE RATIOS vs n_steps (seed 42, d=0.3)")
print("="*80)

print(f"\n  {'n':>4s} {'|λ₀|':>8s} {'|λ₁/λ₀|':>10s} {'|λ₂/λ₀|':>10s} {'|λ₃/λ₀|':>10s} {'gap/Δ':>8s}")
print("  " + "-" * 55)

for n_steps in [5, 10, 20, 30, 50, 75, 100, 200, 500]:
    ent = Entangler(corr, seed=42)
    ent.build(n_steps=n_steps, discount=0.3)
    ev = eigendata(ent.joint)

    l0 = ev[0].abs().item()
    r1 = ev[1].abs().item() / l0 if l0 > 1e-12 else 0
    r2 = ev[2].abs().item() / l0 if l0 > 1e-12 else 0
    r3 = ev[3].abs().item() / l0 if l0 > 1e-12 else 0
    gap = -math.log(r1) if r1 > 1e-12 else 0

    print(f"  {n_steps:>4d} {l0:>8.4f} {r1:>10.4f} {r2:>10.4f} {r3:>10.4f} {gap/DELTA:>8.3f}")


# ═══════════════════════════════════════════════════════════════
#  Critical question: does the gap CONVERGE to a limit?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("DOES THE GAP CONVERGE? (many seeds, d=0.3)")
print("="*80)

for n_steps in [10, 20, 50, 100, 200]:
    gaps = []
    for seed in range(50):
        ent = Entangler(corr, seed=seed)
        ent.build(n_steps=n_steps, discount=0.3)
        g = composition_gap(ent.joint)
        gaps.append(g)

    avg = sum(gaps) / len(gaps)
    std = (sum((g - avg)**2 for g in gaps) / len(gaps)) ** 0.5
    print(f"  n={n_steps:>3d}: gap/Δ = {avg/DELTA:.4f} ± {std/DELTA:.4f}")


# ═══════════════════════════════════════════════════════════════
#  What if we build a crystal from EXACTLY the DS fixed point?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CRYSTAL FROM THE DS FIXED POINT")
print("="*80)

# At the DS fixed point: Born(θ) = 1/27, dominant singleton ≈ 0.963
# Build a correlation matrix that IS the fixed point
# Identity crystal at fixed point: correlation = identity, Born(θ) = 1/27

# What eigenvalue ratio does this crystal have?
# The crystal's eigenvalue structure is set by the entangler dynamics.
# At the DS fixed point (n_steps → ∞, d → 0), the crystal should
# have a specific eigenvalue structure.

# Build with many steps, small discount (approaching fixed point)
print(f"\n  Approaching DS fixed point (increasing n, decreasing d):")
for n, d in [(50, 0.3), (100, 0.15), (200, 0.075), (500, 0.03), (1000, 0.015)]:
    gaps = []
    for seed in range(10):
        try:
            ent = Entangler(corr, seed=seed)
            ent.build(n_steps=n, discount=d)
            g = composition_gap(ent.joint)
            if not math.isnan(g) and abs(g) < 100:
                gaps.append(g)
        except:
            pass
    if gaps:
        avg = sum(gaps) / len(gaps)
        print(f"    n={n:>4d}, d={d:.4f}, n×d={n*d:>5.1f}: gap/Δ = {avg/DELTA:.4f}")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
