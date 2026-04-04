"""
The A₂ Lie bracket from crystal composition.

The Lie bracket [A,B] = AB - BA is the fundamental operation of a Lie algebra.
For A₂ = su(3), the bracket of two simple roots gives the third root:
  [α₁, α₂] ∝ α₁ + α₂ (the third positive root)

In the crystal framework:
  compose(A,B) - compose(B,A) = the crystal commutator

If this commutator projects onto the third transposition's direction
in the standard representation, the A₂ Lie algebra IS the crystal's
composition algebra.

Use the corrected complex S₃ projector throughout.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_fidelity, EPS_LOG)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)


def s3_fractions(joint):
    """Correct complex S₃ projector."""
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
    bp = joint.abs().pow(2)
    total = bp.sum().item()
    if total < 1e-15:
        return 1.0, 0.0, 0.0, triv, sign, std
    f_triv = triv.abs().pow(2).sum().item() / total
    f_sign = sign.abs().pow(2).sum().item() / total
    f_std = std.abs().pow(2).sum().item() / total
    return f_triv, f_sign, f_std, triv, sign, std


def standard_angle(joint):
    """Get angle in standard rep."""
    _, _, _, _, _, std = s3_fractions(joint)
    U, S, Vh = torch.linalg.svd(std)
    if S[0].abs() < 1e-15:
        return None
    v = U[:, 0][:H]
    e1 = torch.tensor([1, -1, 0], dtype=torch.cfloat) / math.sqrt(2)
    e2 = torch.tensor([1, 1, -2], dtype=torch.cfloat) / math.sqrt(6)
    c1 = torch.dot(v.conj(), e1).real.item()
    c2 = torch.dot(v.conj(), e2).real.item()
    return math.atan2(c2, c1)


# Build seed-averaged crystals for all 3 transpositions
swap01_corr = torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32)
swap02_corr = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)
swap12_corr = torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32)
identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

n_seeds = 200

transpositions = {}
for name, corr in [("(01)", swap01_corr), ("(02)", swap02_corr), ("(12)", swap12_corr)]:
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    transpositions[name] = avg

# Also build identity crystal
avg_id = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for seed in range(n_seeds):
    avg_id += Entangler(identity_corr, seed=seed).build().joint
avg_id /= n_seeds


print("=" * 80)
print("THE A₂ LIE BRACKET FROM CRYSTAL COMPOSITION")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Standard sector angles of each transposition
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Transposition standard sector angles ---\n")
for name, crystal in transpositions.items():
    angle = standard_angle(crystal)
    f_t, f_s, f_st, _, _, _ = s3_fractions(crystal)
    print(f"  {name}: angle = {angle/math.pi:.4f}π, f_std = {f_st:.5f}")


# ═══════════════════════════════════════════════════════════════
#  The commutator: compose(A,B) - compose(B,A)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CRYSTAL COMMUTATORS [A, B] = compose(A,B) - compose(B,A)")
print("="*80)

pairs = [
    ("(01)", "(02)"),
    ("(01)", "(12)"),
    ("(02)", "(12)"),
]

for name_a, name_b in pairs:
    A = transpositions[name_a]
    B = transpositions[name_b]

    AB = compose(A, B)
    BA = compose(B, A)
    commutator = AB - BA

    # Normalize commutator
    comm_norm = commutator.abs().pow(2).sum().sqrt().item()

    # Irrep decomposition of commutator
    f_t, f_s, f_st, triv_c, sign_c, std_c = s3_fractions(commutator)

    # Angle in standard rep
    angle_comm = standard_angle(commutator)

    # Also measure AB and BA separately
    angle_AB = standard_angle(AB)
    angle_BA = standard_angle(BA)

    # The THIRD transposition (the one not in the pair)
    remaining = [n for n in ["(01)", "(02)", "(12)"] if n != name_a and n != name_b][0]
    angle_remaining = standard_angle(transpositions[remaining])

    print(f"\n  [{name_a}, {name_b}]:")
    print(f"    ||[A,B]|| = {comm_norm:.6f}")
    print(f"    Irreps: triv={f_t:.5f}, sign={f_s:.5f}, std={f_st:.5f}")
    print(f"    Commutator angle: {angle_comm/math.pi:.4f}π" if angle_comm else "    Commutator angle: None")
    print(f"    compose(A,B) angle: {angle_AB/math.pi:.4f}π" if angle_AB else "    compose(A,B) angle: None")
    print(f"    compose(B,A) angle: {angle_BA/math.pi:.4f}π" if angle_BA else "    compose(B,A) angle: None")
    print(f"    Third transposition ({remaining}) angle: {angle_remaining/math.pi:.4f}π")

    if angle_comm is not None and angle_remaining is not None:
        da = abs(angle_comm - angle_remaining)
        da = min(da, 2*math.pi - da)
        print(f"    Angular distance to third: {da/math.pi:.4f}π")


# ═══════════════════════════════════════════════════════════════
#  Commutator with identity (should be zero)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMMUTATOR WITH IDENTITY (should be near-zero)")
print("="*80)

for name in transpositions:
    comm = compose(avg_id, transpositions[name]) - compose(transpositions[name], avg_id)
    comm_norm = comm.abs().pow(2).sum().sqrt().item()
    print(f"  [id, {name}]: ||comm|| = {comm_norm:.6f}")


# ═══════════════════════════════════════════════════════════════
#  The anticommutator: compose(A,B) + compose(B,A)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CRYSTAL ANTICOMMUTATORS {{A, B}} = compose(A,B) + compose(B,A)")
print("="*80)

for name_a, name_b in pairs:
    A = transpositions[name_a]
    B = transpositions[name_b]

    AB = compose(A, B)
    BA = compose(B, A)
    anticomm = AB + BA

    # Normalize
    anticomm_norm = anticomm.abs().pow(2).sum().sqrt().item()

    f_t, f_s, f_st, _, _, _ = s3_fractions(anticomm)

    # Fidelity to identity crystal
    fid_id = born_fidelity(anticomm, avg_id)

    # Fidelity to each transposition
    fids = {}
    for name, crystal in transpositions.items():
        fids[name] = born_fidelity(anticomm, crystal)

    print(f"\n  {{{name_a}, {name_b}}}:")
    print(f"    Irreps: triv={f_t:.5f}, sign={f_s:.5f}, std={f_st:.5f}")
    print(f"    Fidelity to identity: {fid_id:.5f}")
    for name, fid in fids.items():
        print(f"    Fidelity to {name}: {fid:.5f}")


# ═══════════════════════════════════════════════════════════════
#  Sign content of commutator vs anticommutator
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMMUTATOR → SIGN, ANTICOMMUTATOR → TRIVIAL?")
print("="*80)

print(f"\n  If the Lie bracket structure is clean:")
print(f"    Commutator [A,B] should be antisymmetric → sign sector")
print(f"    Anticommutator {{A,B}} should be symmetric → trivial sector")

for name_a, name_b in pairs:
    A = transpositions[name_a]
    B = transpositions[name_b]
    AB = compose(A, B)
    BA = compose(B, A)

    comm = AB - BA
    anti = AB + BA

    _, f_s_c, f_st_c, _, _, _ = s3_fractions(comm)
    f_t_a, f_s_a, f_st_a, _, _, _ = s3_fractions(anti)

    print(f"\n  [{name_a},{name_b}]: sign={f_s_c:.5f}, std={f_st_c:.5f}")
    print(f"  {{{name_a},{name_b}}}: triv={f_t_a:.5f}, sign={f_s_a:.5f}, std={f_st_a:.5f}")


print(f"\n\n{'='*80}")
print("WHAT THE LIE BRACKET REVEALS")
print("="*80)
