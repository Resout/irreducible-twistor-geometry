"""
The product state has ~8% standard sector under composition.
Is this real or an artifact?

A product state m(i,j) = m(i)·m(j) has EXACT analytical form.
If m = [α, α, α, β] (S₃-symmetric in h), the product state
is already in the trivial sector → f_std = 0 exactly.

If m = [α₁, α₂, α₃, β] (NOT S₃-symmetric), the product state
has nonzero standard content proportional to the asymmetry.

Question: does the composition fixed point have equal h-marginals?
If not, the standard sector residual measures the h-asymmetry
of the product state — which reflects the specific starting crystal.

Also: the p* = √(30/23) formula gives Born(θ)_fp = 10/33 for
bijections. But the h-level distribution at the fixed point
might NOT be uniform. If it reflects the starting permutation,
the product state is a TYPE-SPECIFIC projector.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities, EPS_LOG)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)


def s3_fractions(joint):
    """S₃ irrep fractions."""
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]
    signs = [1, -1, -1, -1, 1, 1]

    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M.real.float() @ P.T

    triv = sum(permute(joint, p) for p in perms) / 6
    sign = sum(s * permute(joint, p) for p, s in zip(perms, signs)) / 6

    bp = joint.abs().pow(2)
    total = bp.sum().item()
    if total < 1e-15:
        return 1.0, 0.0, 0.0

    f_triv = triv.to(joint.dtype).abs().pow(2).sum().item() / total
    f_sign = sign.to(joint.dtype).abs().pow(2).sum().item() / total
    f_std = max(0, 1 - f_triv - f_sign)
    return f_triv, f_sign, f_std


# ═══════════════════════════════════════════════════════════════
#  Analytical product states
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("PRODUCT STATE IRREP DECOMPOSITION")
print("=" * 80)


# 1. Perfect S₃-symmetric product state
print(f"\n--- 1. S₃-symmetric product state: m = [α, α, α, β] ---\n")

for p_star in [1.0, 1.1421, 1.333]:
    # Born(h) = 1/(H+p²), Born(θ) = p²/(H+p²)
    born_h = 1 / (H + p_star**2)
    born_theta = p_star**2 / (H + p_star**2)

    # Mass amplitudes (Born = |m|², so m = √Born)
    alpha = math.sqrt(born_h)
    beta = math.sqrt(born_theta)

    # Product state: m(i,j) = m(i) * m(j)
    m = torch.tensor([alpha, alpha, alpha, beta], dtype=torch.cfloat)
    joint = torch.einsum("i,j->ij", m, m)

    f_t, f_s, f_st = s3_fractions(joint)
    sn = schmidt_number(joint)

    print(f"  p* = {p_star:.4f}: Born(θ) = {born_theta:.4f}, "
          f"trivial = {f_t:.6f}, sign = {f_s:.6f}, standard = {f_st:.6f}, Schmidt = {sn:.4f}")


# 2. Asymmetric product state (h-values differ)
print(f"\n--- 2. Asymmetric product states ---\n")

test_marginals = [
    ("uniform-h", [0.30, 0.30, 0.30, 0.10]),
    ("slight-asym", [0.32, 0.30, 0.28, 0.10]),
    ("moderate-asym", [0.35, 0.30, 0.25, 0.10]),
    ("strong-asym", [0.50, 0.30, 0.10, 0.10]),
    ("extreme-asym", [0.80, 0.10, 0.05, 0.05]),
]

print(f"  {'name':>16s}  {'trivial':>8s}  {'sign':>8s}  {'standard':>8s}")
for name, born in test_marginals:
    m = torch.tensor([math.sqrt(b) for b in born], dtype=torch.cfloat)
    joint = torch.einsum("i,j->ij", m, m)
    f_t, f_s, f_st = s3_fractions(joint)
    print(f"  {name:>16s}  {f_t:>8.5f}  {f_s:>8.5f}  {f_st:>8.5f}")


# ═══════════════════════════════════════════════════════════════
#  Actual composition fixed points: h-marginal distribution
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMPOSITION FIXED POINT: h-MARGINALS")
print("="*80)

identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
inverse_corr = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)
cycle_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)

corrs = {
    "identity": identity_corr,
    "inverse": inverse_corr,
    "cycle": cycle_corr,
}

n_seeds = 50
n_compositions = 30

for name, corr in corrs.items():
    # Build crystal, then compose many times
    avg_joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        ent = Entangler(corr, seed=seed).build(n_steps=50, discount=0.3)
        current = ent.joint.clone()
        original = ent.joint.clone()
        for _ in range(n_compositions):
            current = compose(current, original)
        avg_joint += current
    avg_joint /= n_seeds

    # Marginals
    bp = born_probabilities(avg_joint)
    bp2d = bp.reshape(MASS_DIM, MASS_DIM)

    # Row marginal (sum over columns)
    row_marg = bp2d.sum(dim=1)
    # Column marginal (sum over rows)
    col_marg = bp2d.sum(dim=0)

    f_t, f_s, f_st = s3_fractions(avg_joint)
    sn = schmidt_number(avg_joint)

    # h-asymmetry: max - min of h-marginals
    h_vals = row_marg[:H].tolist()
    h_asym = max(h_vals) - min(h_vals)

    print(f"\n  {name}:")
    print(f"    Row marginal Born: [{', '.join(f'{v:.5f}' for v in row_marg.tolist())}]")
    print(f"    Col marginal Born: [{', '.join(f'{v:.5f}' for v in col_marg.tolist())}]")
    print(f"    h-asymmetry: {h_asym:.6f}")
    print(f"    Irreps: trivial={f_t:.5f}, sign={f_s:.5f}, standard={f_st:.5f}")
    print(f"    Schmidt: {sn:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Does seed-averaging eliminate the standard sector?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SEED-AVERAGE THEN COMPOSE vs COMPOSE THEN SEED-AVERAGE")
print("="*80)

print(f"\n  Strategy A: Average 50 seeds FIRST, then compose 30 times")
print(f"  Strategy B: Compose each seed 30 times, then average")

for name, corr in corrs.items():
    # Strategy A: average first
    avg_crystal = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        ent = Entangler(corr, seed=seed).build()
        avg_crystal += ent.joint
    avg_crystal /= n_seeds

    current_A = avg_crystal.clone()
    for _ in range(n_compositions):
        current_A = compose(current_A, avg_crystal)

    f_t_A, f_s_A, f_st_A = s3_fractions(current_A)

    # Strategy B: compose first (already done above, recompute)
    avg_composed = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        ent = Entangler(corr, seed=seed).build()
        current = ent.joint.clone()
        original = ent.joint.clone()
        for _ in range(n_compositions):
            current = compose(current, original)
        avg_composed += current
    avg_composed /= n_seeds

    f_t_B, f_s_B, f_st_B = s3_fractions(avg_composed)

    print(f"\n  {name}:")
    print(f"    Strategy A (avg→compose): triv={f_t_A:.5f}, sign={f_s_A:.5f}, std={f_st_A:.5f}")
    print(f"    Strategy B (compose→avg): triv={f_t_B:.5f}, sign={f_s_B:.5f}, std={f_st_B:.5f}")


# ═══════════════════════════════════════════════════════════════
#  Single seed: track h-marginals through composition
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SINGLE SEED: h-MARGINALS THROUGH COMPOSITION")
print("="*80)

ent = Entangler(identity_corr, seed=42).build()
current = ent.joint.clone()
original = ent.joint.clone()

print(f"\n  Identity crystal, seed 42:")
print(f"  {'comp':>4s}  {'h0':>8s}  {'h1':>8s}  {'h2':>8s}  {'θ':>8s}  {'h_asym':>8s}  {'f_std':>8s}  {'Schmidt':>8s}")

for step in range(20):
    bp = born_probabilities(current).reshape(MASS_DIM, MASS_DIM)
    row_marg = bp.sum(dim=1)
    h_vals = row_marg[:H].tolist()
    h_asym = max(h_vals) - min(h_vals)
    f_t, f_s, f_st = s3_fractions(current)
    sn = schmidt_number(current)

    print(f"  {step:4d}  {h_vals[0]:>8.5f}  {h_vals[1]:>8.5f}  {h_vals[2]:>8.5f}  "
          f"{row_marg[H].item():>8.5f}  {h_asym:>8.5f}  {f_st:>8.5f}  {sn:>8.4f}")

    current = compose(current, original)


# ═══════════════════════════════════════════════════════════════
#  Correlation between h-asymmetry and standard fraction
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("h-ASYMMETRY vs STANDARD FRACTION")
print("="*80)

asym_vals = []
std_vals = []

for seed in range(100):
    ent = Entangler(identity_corr, seed=seed).build()
    current = ent.joint.clone()
    original = ent.joint.clone()

    # Compose 15 times (near the fixed point)
    for _ in range(15):
        current = compose(current, original)

    bp = born_probabilities(current).reshape(MASS_DIM, MASS_DIM)
    row_marg = bp.sum(dim=1)
    h_vals = row_marg[:H].tolist()
    h_asym = max(h_vals) - min(h_vals)

    f_t, f_s, f_st = s3_fractions(current)

    asym_vals.append(h_asym)
    std_vals.append(f_st)

# Correlation
mean_a = sum(asym_vals) / len(asym_vals)
mean_s = sum(std_vals) / len(std_vals)
cov = sum((a-mean_a)*(s-mean_s) for a, s in zip(asym_vals, std_vals)) / len(asym_vals)
var_a = sum((a-mean_a)**2 for a in asym_vals) / len(asym_vals)
var_s = sum((s-mean_s)**2 for s in std_vals) / len(std_vals)
corr = cov / (var_a * var_s) ** 0.5 if var_a > 0 and var_s > 0 else 0

print(f"\n  Over 100 seeds (identity crystal, 15 compositions):")
print(f"  avg h-asymmetry = {mean_a:.5f}")
print(f"  avg f_std = {mean_s:.5f}")
print(f"  correlation(h_asym, f_std) = {corr:.4f}")

# Distribution of f_std
sorted_std = sorted(std_vals)
print(f"\n  f_std distribution:")
print(f"    min = {sorted_std[0]:.5f}")
print(f"    25% = {sorted_std[25]:.5f}")
print(f"    50% = {sorted_std[50]:.5f}")
print(f"    75% = {sorted_std[75]:.5f}")
print(f"    max = {sorted_std[99]:.5f}")


print(f"\n\n{'='*80}")
print("WHAT THE PRODUCT STATE REVEALS")
print("="*80)
