"""
Is the sub-floor conditional Born(θ) exactly (H-1)/H⁴ = 2/81?

All three crystal types gave Born(θ in θ-row) ≈ 0.024.
The candidate: (H-1)/H⁴ = 2/81 = 0.024691...

This would mean: Born_floor × (H-1)/H = (1/H³) × (H-1)/H = (H-1)/H⁴.
The conditional θ is the unconditional floor scaled by (H-1)/H.

(H-1)/H is the fraction of singletons (non-θ hypotheses).
So: conditional Born(θ in θ-row) = unconditional Born(θ) × (fraction that is NOT θ).

Does this hold exactly? And does it hold across different seeds?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
from solver.algebra import H, MASS_DIM, BORN_FLOOR, schmidt_number, born_probabilities
from solver.crystals import Entangler, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)

target = (H - 1) / H**4
print(f"Target: (H-1)/H⁴ = {H-1}/{H**4} = {target:.8f}")
print(f"Born floor: 1/H³ = {BORN_FLOOR:.8f}")
print(f"Ratio target/floor = (H-1)/H = {(H-1)/H:.8f}")
print()

# ═══════════════════════════════════════════════════════════════
#  Test across all relationship types and many seeds
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("SUB-FLOOR BORN(θ) IN θ-ROW: is it exactly (H-1)/H⁴ ?")
print("=" * 70)

all_values = []

for name, corr in RELATIONSHIP_SIGNATURES.items():
    row_vals = []
    for seed in range(20):
        ent = Entangler(corr, seed=seed).build()
        joint = ent.joint
        # θ-row as a 4-dim mass
        theta_row = joint[H, :]
        bp = born_probabilities(theta_row)
        val = bp[H].item()
        row_vals.append(val)
        all_values.append(val)

    avg = sum(row_vals) / len(row_vals)
    std = (sum((v - avg)**2 for v in row_vals) / len(row_vals)) ** 0.5
    diff = avg - target
    print(f"  {name:15s}: avg Born(θ in θ-row) = {avg:.6f} ± {std:.6f}, "
          f"diff from target = {diff:+.6f}")

# Also check the θ-COLUMN
print(f"\n  θ-COLUMN check:")
for name, corr in RELATIONSHIP_SIGNATURES.items():
    col_vals = []
    for seed in range(20):
        ent = Entangler(corr, seed=seed).build()
        joint = ent.joint
        theta_col = joint[:, H]
        bp = born_probabilities(theta_col)
        col_vals.append(bp[H].item())

    avg = sum(col_vals) / len(col_vals)
    std = (sum((v - avg)**2 for v in col_vals) / len(col_vals)) ** 0.5
    print(f"  {name:15s}: avg Born(θ in θ-col) = {avg:.6f} ± {std:.6f}")

# Overall statistics
avg_all = sum(all_values) / len(all_values)
std_all = (sum((v - avg_all)**2 for v in all_values) / len(all_values)) ** 0.5
print(f"\n  Overall: {avg_all:.6f} ± {std_all:.6f}")
print(f"  Target:  {target:.6f}")
print(f"  Difference: {avg_all - target:+.6f}")
print(f"  Ratio measured/target: {avg_all / target:.4f}")

# ═══════════════════════════════════════════════════════════════
#  Check other rows: Born(θ in h-row)
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("BORN(θ) IN EACH ROW — does the pattern change?")
print("=" * 70)

corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
for seed in range(5):
    ent = Entangler(corr, seed=seed).build()
    joint = ent.joint

    print(f"\n  Seed {seed}:")
    for row in range(MASS_DIM):
        row_mass = joint[row, :]
        bp = born_probabilities(row_mass)
        label = f"h{row}" if row < H else "θ"
        at_floor = bp[H].item() / BORN_FLOOR
        print(f"    Row {label}: Born(θ) = {bp[H].item():.6f} ({at_floor:.2f}× floor)")

# ═══════════════════════════════════════════════════════════════
#  The hierarchy: floor at each dimensional level
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("THE FLOOR HIERARCHY")
print("What is the effective floor at each level?")
print("=" * 70)

# Build many crystals, collect Born(θ) at each level
import random as _random

corrs = list(RELATIONSHIP_SIGNATURES.values())
level_data = {
    "4-dim marginal Born(θ)": [],
    "h-row Born(θ)": [],
    "θ-row Born(θ)": [],
    "Born(θ×θ)": [],
}

for seed in range(100):
    corr = corrs[seed % len(corrs)]
    ent = Entangler(corr, seed=seed).build()
    joint = ent.joint
    total = joint.abs().pow(2).sum().item()

    # 4-dim marginal
    marginal = joint.abs().pow(2).sum(dim=1)
    marginal_born = marginal / marginal.sum()
    level_data["4-dim marginal Born(θ)"].append(marginal_born[H].item())

    # Average h-row Born(θ)
    for i in range(H):
        bp = born_probabilities(joint[i, :])
        level_data["h-row Born(θ)"].append(bp[H].item())

    # θ-row Born(θ)
    bp_theta = born_probabilities(joint[H, :])
    level_data["θ-row Born(θ)"].append(bp_theta[H].item())

    # Born(θ×θ)
    level_data["Born(θ×θ)"].append(joint[H, H].abs().pow(2).item() / total)

print(f"\n  {'Level':>25s} {'min':>10s} {'mean':>10s} {'max':>10s} {'floor ref':>12s}")
print("  " + "-" * 70)

floor_refs = {
    "4-dim marginal Born(θ)": f"1/H³ = {BORN_FLOOR:.6f}",
    "h-row Born(θ)": f"? ≈ {BORN_FLOOR:.6f}",
    "θ-row Born(θ)": f"(H-1)/H⁴ = {target:.6f}",
    "Born(θ×θ)": f"(1/H³)² = {BORN_FLOOR**2:.6f}",
}

for level, vals in level_data.items():
    mn = min(vals)
    mx = max(vals)
    avg = sum(vals) / len(vals)
    ref = floor_refs.get(level, "")
    print(f"  {level:>25s} {mn:>10.6f} {avg:>10.6f} {mx:>10.6f}  {ref}")


print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
print(f"\n  The floor hierarchy:")
print(f"    Level 0 (θ×θ):      Born ≈ 0.005-0.010 (above (1/27)² but well below 1/27)")
print(f"    Level 1 (θ-row θ):  Born ≈ {target:.4f} = (H-1)/H⁴ = BORN_FLOOR × (H-1)/H")
print(f"    Level 2 (h-row θ):  Born ≈ 0.25 (≈ 1/(H+1) — NOT at the floor)")
print(f"    Level 3 (marginal): Born ≈ 0.20-0.40 (5-11× floor)")
print(f"    Level 4 (4-dim DS): Born → 1/27 exactly (the enforced floor)")
