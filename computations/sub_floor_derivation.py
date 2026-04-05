"""
Derive the sub-floor constant from H=3 principles.

The θ-row and θ-column of the joint mass have Born(θ) ≈ 0.022-0.024,
constant across all crystal types. Previous script tested (H-1)/H⁴ = 2/81
= 0.02469... as the candidate for θ-row.

This script:
1. Precise measurement over 500 seeds × 6 types
2. Test multiple candidate expressions against the data
3. Separate θ-row from θ-column (they may differ)
4. Check if the constant depends on entangler parameters
5. Try to derive it from the joint mass structure analytically
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, BORN_FLOOR, K_STAR, DELTA, schmidt_number, born_probabilities
from solver.crystals import Entangler, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)

N_SEEDS = 500


# ═══════════════════════════════════════════════════════════════
#  Candidate expressions from H=3 constants
# ═══════════════════════════════════════════════════════════════

candidates = {
    "(H-1)/H⁴":       (H-1) / H**4,                    # 2/81 = 0.024691
    "1/H²(H+1)":      1 / (H**2 * (H+1)),              # 1/36 = 0.027778
    "K*/H²":           K_STAR / H**2,                    # 7/270 = 0.025926
    "DELTA/H³":        DELTA / H**3,                     # 0.009852
    "1/H⁴":           1 / H**4,                          # 1/81 = 0.012346
    "2/(H³(H+1))":    2 / (H**3 * (H+1)),               # 2/108 = 0.018519
    "(H-1)/(H³(H+1))": (H-1) / (H**3 * (H+1)),         # 2/108 = 0.018519
    "BORN²×H":         BORN_FLOOR**2 * H,                # 3/729 = 0.004115
    "1/(H²+H³)":      1 / (H**2 + H**3),                # 1/36 = 0.027778
    "7/(H⁴(H+1))":    7 / (H**4 * (H+1)),               # 7/324 = 0.021605
    "K*/(H²+1)":      K_STAR / (H**2 + 1),              # 7/300 = 0.023333
    "7/H⁵":           7 / H**5,                          # 7/243 = 0.028807
    "K*×BORN":         K_STAR * BORN_FLOOR,               # 7/810 = 0.008642
    "(H²-H+1)/H⁵":   (H**2 - H + 1) / H**5,           # 7/243 = 0.028807
    "BORN×(H-1)/H":   BORN_FLOOR * (H-1) / H,           # 2/81 = 0.024691
}

print("Candidate values:")
for name, val in sorted(candidates.items(), key=lambda x: x[1]):
    print(f"  {name:>20s} = {val:.8f}")


# ═══════════════════════════════════════════════════════════════
#  High-precision measurement: θ-row and θ-column Born(θ)
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print(f"SUB-FLOOR CONSTANT: {N_SEEDS} SEEDS × 6 TYPES")
print(f"{'='*80}")

theta_row_all = []
theta_col_all = []
by_type = {}

for name, corr in RELATIONSHIP_SIGNATURES.items():
    row_vals = []
    col_vals = []

    for seed in range(N_SEEDS):
        ent = Entangler(corr, seed=seed).build()
        joint = ent.joint

        # θ-row: joint[H, :] treated as 4-dim mass
        theta_row = joint[H, :]
        bp_row = born_probabilities(theta_row)
        row_vals.append(bp_row[H].item())

        # θ-column: joint[:, H] treated as 4-dim mass
        theta_col = joint[:, H]
        bp_col = born_probabilities(theta_col)
        col_vals.append(bp_col[H].item())

    theta_row_all.extend(row_vals)
    theta_col_all.extend(col_vals)

    r_avg = sum(row_vals) / len(row_vals)
    r_std = (sum((v - r_avg)**2 for v in row_vals) / len(row_vals)) ** 0.5
    c_avg = sum(col_vals) / len(col_vals)
    c_std = (sum((v - c_avg)**2 for v in col_vals) / len(col_vals)) ** 0.5

    by_type[name] = {"row_avg": r_avg, "row_std": r_std, "col_avg": c_avg, "col_std": c_std}

    print(f"\n  {name}:")
    print(f"    θ-row Born(θ): {r_avg:.8f} ± {r_std:.8f}")
    print(f"    θ-col Born(θ): {c_avg:.8f} ± {c_std:.8f}")


# Grand averages
row_grand = sum(theta_row_all) / len(theta_row_all)
row_grand_std = (sum((v - row_grand)**2 for v in theta_row_all) / len(theta_row_all)) ** 0.5
row_sem = row_grand_std / len(theta_row_all)**0.5

col_grand = sum(theta_col_all) / len(theta_col_all)
col_grand_std = (sum((v - col_grand)**2 for v in theta_col_all) / len(theta_col_all)) ** 0.5
col_sem = col_grand_std / len(theta_col_all)**0.5

print(f"\n{'='*80}")
print(f"GRAND AVERAGES (n={len(theta_row_all)})")
print(f"{'='*80}")
print(f"  θ-row Born(θ) = {row_grand:.8f} ± {row_grand_std:.8f}  (SEM = {row_sem:.8f})")
print(f"  θ-col Born(θ) = {col_grand:.8f} ± {col_grand_std:.8f}  (SEM = {col_sem:.8f})")
print(f"  Row-col difference: {abs(row_grand - col_grand):.8f}")


# ═══════════════════════════════════════════════════════════════
#  Match candidates to data
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print("CANDIDATE MATCH (θ-row)")
print(f"{'='*80}")
print(f"  Measured: {row_grand:.8f}")
print(f"\n  {'Candidate':>20s} {'Value':>12s} {'Diff':>12s} {'σ away':>8s}")
print("  " + "-" * 60)

for name, val in sorted(candidates.items(), key=lambda x: abs(x[1] - row_grand)):
    diff = val - row_grand
    sigma = abs(diff) / row_sem if row_sem > 0 else float('inf')
    marker = "  ★" if sigma < 2 else "  ◇" if sigma < 5 else ""
    print(f"  {name:>20s} {val:>12.8f} {diff:>+12.8f} {sigma:>8.1f}{marker}")


print(f"\n{'='*80}")
print("CANDIDATE MATCH (θ-col)")
print(f"{'='*80}")
print(f"  Measured: {col_grand:.8f}")
print(f"\n  {'Candidate':>20s} {'Value':>12s} {'Diff':>12s} {'σ away':>8s}")
print("  " + "-" * 60)

for name, val in sorted(candidates.items(), key=lambda x: abs(x[1] - col_grand)):
    diff = val - col_grand
    sigma = abs(diff) / col_sem if col_sem > 0 else float('inf')
    marker = "  ★" if sigma < 2 else "  ◇" if sigma < 5 else ""
    print(f"  {name:>20s} {val:>12.8f} {diff:>+12.8f} {sigma:>8.1f}{marker}")


# ═══════════════════════════════════════════════════════════════
#  Does the constant change with discount/n_steps?
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print("PARAMETER DEPENDENCE: does the sub-floor constant change?")
print(f"{'='*80}")

corr = torch.eye(H, dtype=torch.float32)

# Discount sweep
print(f"\n  Discount sweep (n_steps=50, 20 seeds):")
for d_pct in [10, 15, 20, 25, 30, 35, 40]:
    d = d_pct / 100.0
    vals = []
    for seed in range(20):
        try:
            ent = Entangler(corr, seed=seed)
            ent.build(n_steps=50, discount=d)
            theta_row = ent.joint[H, :]
            bp = born_probabilities(theta_row)
            v = bp[H].item()
            if not math.isnan(v) and v < 1.0:
                vals.append(v)
        except:
            pass
    if vals:
        avg = sum(vals) / len(vals)
        print(f"    d={d:.2f}: Born(θ in θ-row) = {avg:.6f} (n={len(vals)})")

# Step count sweep
print(f"\n  Step count sweep (discount=0.3, 20 seeds):")
for n_steps in [10, 20, 30, 50, 75, 100]:
    vals = []
    for seed in range(20):
        try:
            ent = Entangler(corr, seed=seed)
            ent.build(n_steps=n_steps, discount=0.3)
            theta_row = ent.joint[H, :]
            bp = born_probabilities(theta_row)
            v = bp[H].item()
            if not math.isnan(v) and v < 1.0:
                vals.append(v)
        except:
            pass
    if vals:
        avg = sum(vals) / len(vals)
        print(f"    n={n_steps:>3d}: Born(θ in θ-row) = {avg:.6f}")


# ═══════════════════════════════════════════════════════════════
#  Analytical structure of the θ-row
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print("θ-ROW STRUCTURE: what determines the sub-floor constant?")
print(f"{'='*80}")

# Look at the full Born distribution of the θ-row
ent = Entangler(torch.eye(H, dtype=torch.float32), seed=42).build()
joint = ent.joint

print(f"\n  θ-row mass (raw complex):")
for j in range(MASS_DIM):
    m = joint[H, j]
    print(f"    m[θ,{j}] = {m.real.item():+.8f} {m.imag.item():+.8f}j  |m|² = {m.abs().pow(2).item():.8f}")

print(f"\n  θ-row Born probabilities:")
bp = born_probabilities(joint[H, :])
for j in range(MASS_DIM):
    label = f"h{j}" if j < H else "θ"
    print(f"    Born({label}) = {bp[j].item():.8f}")

print(f"\n  θ-column mass (raw complex):")
for i in range(MASS_DIM):
    m = joint[i, H]
    print(f"    m[{i},θ] = {m.real.item():+.8f} {m.imag.item():+.8f}j  |m|² = {m.abs().pow(2).item():.8f}")

print(f"\n  θ-column Born probabilities:")
bp = born_probabilities(joint[:, H])
for i in range(MASS_DIM):
    label = f"h{i}" if i < H else "θ"
    print(f"    Born({label}) = {bp[i].item():.8f}")

# Check: is θ-row Born(θ) related to the overall Born(θ) of the marginal?
marginal = joint.abs().pow(2).sum(dim=1)
marginal_bp = marginal / marginal.sum()
print(f"\n  Marginal Born(θ) = {marginal_bp[H].item():.8f}")
print(f"  θ-row Born(θ) = {born_probabilities(joint[H,:])[H].item():.8f}")
print(f"  Ratio: {born_probabilities(joint[H,:])[H].item() / marginal_bp[H].item():.6f}")

# What about Born(θ in θ-row) × weight of θ-row?
theta_row_weight = joint[H, :].abs().pow(2).sum().item() / joint.abs().pow(2).sum().item()
print(f"\n  θ-row weight (fraction of total mass²): {theta_row_weight:.8f}")
print(f"  Product: Born(θ in θ-row) × weight = {born_probabilities(joint[H,:])[H].item() * theta_row_weight:.8f}")
print(f"  Born(θ×θ) = {joint[H,H].abs().pow(2).item() / joint.abs().pow(2).sum().item():.8f}")


print(f"\n{'='*80}")
print("FINDINGS")
print(f"{'='*80}")
