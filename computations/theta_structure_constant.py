"""
The θ-channel contribution to the structure constant.

N = √(MASS_DIM/H) = √(4/3) = 2/√3 ≈ 1.155

Hypothesis: the θ dimension inflates the structure constant.
If we restrict to the h×h block (3×3 submatrix), the structure
constant should be N_hh = 1 (Chevalley normalization).

The full 4×4 structure constant N = N_hh × √(MASS_DIM/H).

This would mean: θ participates in the Lie bracket with weight
(H+1)/H relative to the h-sector. The ignorance dimension
contributes to the gauge field strength.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, EPS_LOG)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)


def compose_block(A, B, block_size):
    """Compose restricted to a subblock."""
    A_sub = A[:block_size, :block_size]
    B_sub = B[:block_size, :block_size]
    result = torch.einsum("ab,bc->ac", A_sub, B_sub)
    re_sum = result.reshape(-1).real.sum()
    if abs(re_sum) > 1e-8:
        result = result / re_sum
    return result


def s3_sign_project_block(joint, block_size):
    """Sign projection restricted to a subblock."""
    perms_h = [
        [0,1,2], [1,0,2], [2,1,0],
        [0,2,1], [1,2,0], [2,0,1],
    ]
    signs = [1, -1, -1, -1, 1, 1]

    sub = joint[:block_size, :block_size]

    def permute(M, p):
        P = torch.zeros(block_size, block_size, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T

    sign_proj = sum(s * permute(sub, p) for p, s in zip(perms_h, signs)) / 6
    return sign_proj


n_seeds = 200

# Build seed-averaged transposition and 3-cycle crystals
corrs = {
    "(01)": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)": torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

crystals = {}
for name, corr in corrs.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    crystals[name] = avg


print("=" * 80)
print("THE θ CONTRIBUTION TO THE STRUCTURE CONSTANT")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Full 4×4 structure constant (recap)
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Full [4,4] commutator ---\n")

trans_pairs = [("(01)", "(02)"), ("(01)", "(12)"), ("(02)", "(12)")]
cycle_sign_norms_full = []

for name in ["(012)", "(021)"]:
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]
    signs = [1, -1, -1, -1, 1, 1]
    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T
    sign_proj = sum(s * permute(crystals[name], p) for p, s in zip(perms, signs)) / 6
    cycle_sign_norms_full.append(sign_proj.abs().pow(2).sum().sqrt().item())

avg_cycle_sign_full = sum(cycle_sign_norms_full) / 2

N_full_list = []
for name_a, name_b in trans_pairs:
    comm = compose(crystals[name_a], crystals[name_b]) - compose(crystals[name_b], crystals[name_a])
    # Sign project the full commutator
    perms = [[0,1,2,3], [1,0,2,3], [2,1,0,3], [0,2,1,3], [1,2,0,3], [2,0,1,3]]
    signs = [1, -1, -1, -1, 1, 1]
    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T
    sign_comm = sum(s * permute(comm, p) for p, s in zip(perms, signs)) / 6
    comm_sign_norm = sign_comm.abs().pow(2).sum().sqrt().item()
    N = comm_sign_norm / avg_cycle_sign_full
    N_full_list.append(N)
    print(f"  [{name_a}, {name_b}]: N_full = {N:.6f}")

N_full_avg = sum(N_full_list) / len(N_full_list)
print(f"  Average N_full = {N_full_avg:.6f}")
print(f"  √(4/3) = {math.sqrt(4/3):.6f}")
print(f"  Ratio N/√(4/3) = {N_full_avg/math.sqrt(4/3):.6f}")


# ═══════════════════════════════════════════════════════════════
#  h×h block (3×3) structure constant
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("h×h BLOCK [3,3] COMMUTATOR (no θ)")
print("="*80)

# Sign norms of 3-cycle crystals in h×h block
cycle_sign_norms_hh = []
for name in ["(012)", "(021)"]:
    sign_proj = s3_sign_project_block(crystals[name], H)
    norm = sign_proj.abs().pow(2).sum().sqrt().item()
    cycle_sign_norms_hh.append(norm)
    print(f"\n  {name} h×h sign norm: {norm:.6f}")

avg_cycle_sign_hh = sum(cycle_sign_norms_hh) / 2

N_hh_list = []
for name_a, name_b in trans_pairs:
    comm_hh = compose_block(crystals[name_a], crystals[name_b], H) - \
              compose_block(crystals[name_b], crystals[name_a], H)

    sign_comm_hh = s3_sign_project_block(
        torch.nn.functional.pad(comm_hh, (0, 1, 0, 1)),  # pad to 4×4 for projection
        H
    )
    # Actually, just project the 3×3 directly
    perms_h = [[0,1,2], [1,0,2], [2,1,0], [0,2,1], [1,2,0], [2,0,1]]
    signs = [1, -1, -1, -1, 1, 1]
    def permute3(M, p):
        P = torch.zeros(H, H, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T

    sign_hh = sum(s * permute3(comm_hh, p) for p, s in zip(perms_h, signs)) / 6
    comm_sign_norm_hh = sign_hh.abs().pow(2).sum().sqrt().item()

    N_hh = comm_sign_norm_hh / avg_cycle_sign_hh if avg_cycle_sign_hh > 1e-10 else 0
    N_hh_list.append(N_hh)
    print(f"\n  [{name_a}, {name_b}] h×h:")
    print(f"    ||sign(comm_hh)|| = {comm_sign_norm_hh:.6f}")
    print(f"    N_hh = {N_hh:.6f}")

N_hh_avg = sum(N_hh_list) / len(N_hh_list)


# ═══════════════════════════════════════════════════════════════
#  Compare: does the ratio = √(MASS_DIM/H)?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE θ FACTOR")
print("="*80)

print(f"\n  N_full = {N_full_avg:.6f}")
print(f"  N_hh   = {N_hh_avg:.6f}")
if N_hh_avg > 1e-10:
    ratio = N_full_avg / N_hh_avg
    print(f"  Ratio N_full/N_hh = {ratio:.6f}")
    print(f"  √(MASS_DIM/H) = √(4/3) = {math.sqrt(MASS_DIM/H):.6f}")
    print(f"  √((H+1)/H) = {math.sqrt((H+1)/H):.6f}")
    print(f"  Match: {abs(ratio - math.sqrt(MASS_DIM/H)) < 0.05}")

# Also check: block decomposition of the commutator norm
print(f"\n  Block decomposition of commutator:")
for name_a, name_b in trans_pairs[:1]:  # just first pair
    comm = compose(crystals[name_a], crystals[name_b]) - compose(crystals[name_b], crystals[name_a])
    bp = comm.abs().pow(2)
    hh_mass = bp[:H, :H].sum().item()
    htheta_mass = bp[:H, H].sum().item() + bp[H, :H].sum().item()
    theta_theta = bp[H, H].item()
    total = bp.sum().item()

    print(f"\n  [{name_a}, {name_b}] Born mass distribution:")
    print(f"    h×h block:     {hh_mass/total:.5f} ({hh_mass/total*100:.1f}%)")
    print(f"    h×θ + θ×h:     {htheta_mass/total:.5f} ({htheta_mass/total*100:.1f}%)")
    print(f"    θ×θ:           {theta_theta/total:.5f} ({theta_theta/total*100:.1f}%)")
    print(f"    h-total:       {hh_mass/total:.5f}")
    print(f"    θ-total:       {(htheta_mass+theta_theta)/total:.5f}")

    # If θ contributes (H+1)/H times as much per dimension:
    # h×h has H²=9 entries, θ-related has 2H+1=7 entries
    # Ratio (H²+2H+1)/(H²) = (H+1)²/H² = 16/9

    print(f"\n    H² = {H**2} entries in h×h, {2*H+1} entries with θ")
    print(f"    Per-entry h×h: {hh_mass/(H**2):.6f}")
    per_theta = (htheta_mass + theta_theta) / (2*H + 1) if (2*H+1) > 0 else 0
    print(f"    Per-entry θ:   {per_theta:.6f}")
    if hh_mass > 0:
        per_entry_ratio = per_theta / (hh_mass / H**2)
        print(f"    Per-entry ratio θ/h = {per_entry_ratio:.4f}")


print(f"\n\n{'='*80}")
print("WHAT THE θ STRUCTURE CONSTANT REVEALS")
print("="*80)
