"""
The standard sector residual at the composition fixed point comes from
COLUMN asymmetry, not row asymmetry. Rows are always uniform (composition
is row-wise). Columns carry the crystal's specific projection bias.

For a single seed: f_std → ~10%. The column marginal Born distribution
at the fixed point encodes WHICH crystal performed the composition.

Questions:
1. Is 10% a universal value or seed/type-dependent?
2. Does the column asymmetry equal some function of H?
3. What sets the specific column marginals?

The column marginal at the fixed point should be the DOMINANT EIGENVECTOR
of the crystal's column operator. Composition A@B: columns of the result
are linear combinations of columns of B, weighted by rows of A. At the
fixed point, the column structure IS the right eigenvector of the crystal.
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
    return f_triv, f_sign, max(0, 1 - f_triv - f_sign)


identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

print("=" * 80)
print("COLUMN ASYMMETRY AT THE COMPOSITION FIXED POINT")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  f_std at fixed point across many seeds
# ═══════════════════════════════════════════════════════════════

n_compositions = 25
n_seeds = 200

std_vals = []
col_asym_vals = []
col_marginals = []

print(f"\n  Computing {n_seeds} seeds × {n_compositions} compositions...")

for seed in range(n_seeds):
    ent = Entangler(identity_corr, seed=seed).build()
    current = ent.joint.clone()
    original = ent.joint.clone()
    for _ in range(n_compositions):
        current = compose(current, original)

    f_t, f_s, f_st = s3_fractions(current)
    std_vals.append(f_st)

    # Column marginal
    bp = current.abs().pow(2)
    bp_norm = bp / bp.sum().clamp(min=1e-15)
    col_marg = bp_norm.sum(dim=0)
    col_marginals.append(col_marg[:H].tolist())

    h_col = col_marg[:H].tolist()
    col_asym_vals.append(max(h_col) - min(h_col))


avg_std = sum(std_vals) / len(std_vals)
std_std = (sum((v - avg_std)**2 for v in std_vals) / len(std_vals)) ** 0.5

avg_col_asym = sum(col_asym_vals) / len(col_asym_vals)
std_col_asym = (sum((v - avg_col_asym)**2 for v in col_asym_vals) / len(col_asym_vals)) ** 0.5

print(f"\n  f_std at fixed point: {avg_std:.5f} ± {std_std:.5f}")
print(f"  Column h-asymmetry:  {avg_col_asym:.5f} ± {std_col_asym:.5f}")

# Distribution
sorted_std = sorted(std_vals)
print(f"\n  f_std distribution (200 seeds):")
print(f"    min = {sorted_std[0]:.5f}")
print(f"    10% = {sorted_std[20]:.5f}")
print(f"    25% = {sorted_std[50]:.5f}")
print(f"    50% = {sorted_std[100]:.5f}")
print(f"    75% = {sorted_std[150]:.5f}")
print(f"    90% = {sorted_std[180]:.5f}")
print(f"    max = {sorted_std[199]:.5f}")


# ═══════════════════════════════════════════════════════════════
#  The column structure IS the crystal's right eigenvector
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COLUMN MARGINAL = RIGHT EIGENVECTOR")
print("="*80)

# For a single seed, compare the composition fixed point's column
# marginal with the crystal's column SVD

for seed in [0, 1, 42, 99]:
    ent = Entangler(identity_corr, seed=seed).build()
    crystal = ent.joint.clone()

    # Right singular vector of the crystal (as a matrix)
    U, S, Vh = torch.linalg.svd(crystal)
    right_sv = Vh[0]  # dominant right singular vector

    # Column marginal of crystal
    bp = crystal.abs().pow(2)
    col_marg_crystal = bp.sum(dim=0) / bp.sum()

    # Composition fixed point column marginal
    current = crystal.clone()
    for _ in range(25):
        current = compose(current, crystal)

    bp_fp = current.abs().pow(2)
    col_marg_fp = bp_fp.sum(dim=0) / bp_fp.sum()

    # Born probabilities of right singular vector
    rv_born = right_sv.abs().pow(2)
    rv_born = rv_born / rv_born.sum()

    print(f"\n  seed {seed}:")
    print(f"    Crystal col marginal:    [{', '.join(f'{v:.4f}' for v in col_marg_crystal.tolist())}]")
    print(f"    Right SV Born:           [{', '.join(f'{v:.4f}' for v in rv_born.tolist())}]")
    print(f"    Fixed point col marginal:[{', '.join(f'{v:.4f}' for v in col_marg_fp.tolist())}]")


# ═══════════════════════════════════════════════════════════════
#  Does S₃-symmetrizing the crystal eliminate f_std at fixed point?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("S₃-SYMMETRIZED CRYSTAL: FIXED POINT f_std")
print("="*80)

perms_list = [
    [0,1,2,3], [1,0,2,3], [2,1,0,3],
    [0,2,1,3], [1,2,0,3], [2,0,1,3],
]

for seed in [0, 42, 99]:
    ent = Entangler(identity_corr, seed=seed).build()
    crystal = ent.joint.clone()

    # S₃-symmetrize
    sym_crystal = torch.zeros_like(crystal)
    for p in perms_list:
        P = torch.zeros(MASS_DIM, MASS_DIM)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        sym_crystal += (P.to(crystal.dtype) @ crystal @ P.T.to(crystal.dtype))
    sym_crystal /= 6

    # Compose
    current = sym_crystal.clone()
    for _ in range(25):
        current = compose(current, sym_crystal)

    f_t, f_s, f_st = s3_fractions(current)

    bp_fp = current.abs().pow(2)
    col_marg = bp_fp.sum(dim=0) / bp_fp.sum()

    print(f"\n  seed {seed} (symmetrized):")
    print(f"    f_std at fixed point: {f_st:.6f}")
    print(f"    Col marginal: [{', '.join(f'{v:.5f}' for v in col_marg.tolist())}]")


# ═══════════════════════════════════════════════════════════════
#  Correlation between column asymmetry and standard fraction
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COLUMN ASYMMETRY DETERMINES STANDARD FRACTION")
print("="*80)

# Compute correlation
mean_a = sum(col_asym_vals) / len(col_asym_vals)
mean_s = sum(std_vals) / len(std_vals)
cov = sum((a-mean_a)*(s-mean_s) for a, s in zip(col_asym_vals, std_vals)) / len(std_vals)
var_a = sum((a-mean_a)**2 for a in col_asym_vals) / len(col_asym_vals)
var_s = sum((s-mean_s)**2 for s in std_vals) / len(std_vals)
corr = cov / (var_a * var_s) ** 0.5 if var_a > 0 and var_s > 0 else 0

print(f"\n  correlation(col_asymmetry, f_std) = {corr:.4f}")
print(f"  (1.0 = perfectly determined, 0.0 = independent)")


print(f"\n\n{'='*80}")
print("WHAT THE COLUMN ASYMMETRY REVEALS")
print("="*80)
