"""
Structure constants of A₂ from crystal composition.

In the A₂ Lie algebra:
  [e_α₁, e_α₂] = N_{12} · e_{α₁+α₂}

where α₁, α₂ are simple roots (transpositions) and α₁+α₂ is the
highest root (3-cycle in S₃).

The crystal commutator is 99.8% sign. The 3-cycle crystal is ~26% sign.
The commutator should be proportional to the sign-projected 3-cycle.

Extract:
1. Born fidelity of commutator with 3-cycle crystals
2. The structure constant N = ||commutator|| / ||sign(3-cycle)||
3. Whether N is universal across all root pairs (A₂: should be ±1)
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_fidelity, EPS_LOG)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)


def s3_project(joint):
    """Return (trivial, sign, standard) projected matrices."""
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

    triv = sum(permute(joint, p) for p in perms) / 6
    sign = sum(s * permute(joint, p) for p, s in zip(perms, signs)) / 6
    std = joint - triv - sign
    return triv, sign, std


n_seeds = 200

# Build all S₃ element crystals
corrs = {
    "e": torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
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
print("STRUCTURE CONSTANTS OF A₂ FROM CRYSTAL COMPOSITION")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Sign projections of all elements
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Sign sector content and norm ---\n")
sign_norms = {}
for name, crystal in crystals.items():
    triv, sign, std = s3_project(crystal)
    sign_norm = sign.abs().pow(2).sum().sqrt().item()
    total_norm = crystal.abs().pow(2).sum().sqrt().item()
    frac = sign.abs().pow(2).sum().item() / crystal.abs().pow(2).sum().item()
    sign_norms[name] = sign_norm
    print(f"  {name:>5s}: sign_norm = {sign_norm:.6f}, total_norm = {total_norm:.6f}, "
          f"sign_frac = {frac:.5f}")


# ═══════════════════════════════════════════════════════════════
#  Commutator fidelity with 3-cycle crystals
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMMUTATOR FIDELITY WITH 3-CYCLE CRYSTALS")
print("="*80)

trans_pairs = [("(01)", "(02)"), ("(01)", "(12)"), ("(02)", "(12)")]

for name_a, name_b in trans_pairs:
    A = crystals[name_a]
    B = crystals[name_b]

    comm = compose(A, B) - compose(B, A)
    comm_norm = comm.abs().pow(2).sum().sqrt().item()

    # Fidelity with each S₃ element
    print(f"\n  [{name_a}, {name_b}] (norm = {comm_norm:.6f}):")
    for ref_name, ref_crystal in crystals.items():
        fid = born_fidelity(comm, ref_crystal)
        print(f"    fid({ref_name:>5s}) = {fid:.5f}")

    # Fidelity with SIGN-PROJECTED 3-cycles
    _, sign_012, _ = s3_project(crystals["(012)"])
    _, sign_021, _ = s3_project(crystals["(021)"])
    _, sign_comm, _ = s3_project(comm)

    # Inner product of sign projections
    inner_012 = torch.sum(sign_comm.conj() * sign_012).item()
    inner_021 = torch.sum(sign_comm.conj() * sign_021).item()

    print(f"    ⟨sign(comm), sign((012))⟩ = {inner_012:.6f}")
    print(f"    ⟨sign(comm), sign((021))⟩ = {inner_021:.6f}")


# ═══════════════════════════════════════════════════════════════
#  The group algebra product: (01)∘(02) = (012), (02)∘(01) = (021)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMPOSITION = GROUP MULTIPLICATION")
print("="*80)

for name_a, name_b in trans_pairs:
    A = crystals[name_a]
    B = crystals[name_b]

    AB = compose(A, B)
    BA = compose(B, A)

    # Which S₃ element does AB most resemble?
    print(f"\n  compose({name_a}, {name_b}):")
    best_fid = 0
    best_name = ""
    for ref_name, ref_crystal in crystals.items():
        fid = born_fidelity(AB, ref_crystal)
        if fid > best_fid:
            best_fid = fid
            best_name = ref_name
        print(f"    fid({ref_name:>5s}) = {fid:.5f}")
    print(f"    → best match: {best_name} (fid = {best_fid:.5f})")

    print(f"\n  compose({name_b}, {name_a}):")
    best_fid = 0
    best_name = ""
    for ref_name, ref_crystal in crystals.items():
        fid = born_fidelity(BA, ref_crystal)
        if fid > best_fid:
            best_fid = fid
            best_name = ref_name
        print(f"    fid({ref_name:>5s}) = {fid:.5f}")
    print(f"    → best match: {best_name} (fid = {best_fid:.5f})")


# ═══════════════════════════════════════════════════════════════
#  Structure constant extraction
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("STRUCTURE CONSTANT N")
print("="*80)

# N = ||sign(commutator)|| / ||sign(3-cycle)||
# Should be universal across all transposition pairs

for name_a, name_b in trans_pairs:
    A = crystals[name_a]
    B = crystals[name_b]
    comm = compose(A, B) - compose(B, A)

    _, sign_comm, _ = s3_project(comm)
    comm_sign_norm = sign_comm.abs().pow(2).sum().sqrt().item()

    # Use average of (012) and (021) sign norms
    avg_cycle_sign = (sign_norms["(012)"] + sign_norms["(021)"]) / 2

    N = comm_sign_norm / avg_cycle_sign if avg_cycle_sign > 1e-10 else 0

    print(f"\n  [{name_a}, {name_b}]:")
    print(f"    ||sign(comm)|| = {comm_sign_norm:.6f}")
    print(f"    ||sign(cycle)|| = {avg_cycle_sign:.6f}")
    print(f"    N = {N:.6f}")

# What should N be?
print(f"\n  For A₂ with standard normalization: |N_{{α₁,α₂}}| = 1")
print(f"  The crystal's N should be related to the composition")
print(f"  contraction rate 299/405 = {299/405:.6f}")


print(f"\n\n{'='*80}")
print("WHAT THE STRUCTURE CONSTANTS REVEAL")
print("="*80)
