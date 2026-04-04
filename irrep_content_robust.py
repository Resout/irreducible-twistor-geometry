"""
S₃ irrep content by conjugacy class: robustness check and precision.

Principle: The three conjugacy classes of S₃ map to the three sectors
of physics through their irrep decomposition:
  Identity:      pure trivial (gravity)
  Transpositions: trivial + standard (gravity + gauge)
  3-cycles:      trivial + sign (gravity + chirality)

This script: verify this is seed-independent, find exact fractions,
and trace the origin of each component.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import numpy as np
from solver.algebra import H, MASS_DIM, K_STAR, BORN_FLOOR
from solver.crystals import Entangler

torch.set_grad_enabled(False)

# S₃ permutations
S3_PERMS = {
    "e":     [0, 1, 2],
    "(01)":  [1, 0, 2],
    "(02)":  [2, 1, 0],
    "(12)":  [0, 2, 1],
    "(012)": [1, 2, 0],
    "(021)": [2, 0, 1],
}

S3_CORR = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

CONJUGACY_CLASSES = {
    "identity":      ["e"],
    "transpositions": ["(01)", "(02)", "(12)"],
    "3-cycles":      ["(012)", "(021)"],
}


# ═══════════════════════════════════════════════════════════════
#  Build S₃ representation on 16D joint mass space
# ═══════════════════════════════════════════════════════════════

dim_sq = MASS_DIM ** 2  # 16

def permutation_matrix(perm):
    """4×4 permutation matrix from permutation of {0,1,2}. θ index fixed."""
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return P

def group_rep_16(perm):
    """16×16 representation: σ acts on vec(M) as (P_σ ⊗ P_σ)vec(M)."""
    P = permutation_matrix(perm)
    return torch.kron(P, P)

# Build projectors
reps = {name: group_rep_16(perm) for name, perm in S3_PERMS.items()}

chi_trivial  = {"e": 1, "(01)": 1, "(02)": 1, "(12)": 1, "(012)": 1, "(021)": 1}
chi_sign     = {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1}
chi_standard = {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1}

def build_projector(chi, dim_chi):
    P = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for name in S3_PERMS:
        P += chi[name] * reps[name]
    P *= dim_chi / 6.0
    return P

P_triv = build_projector(chi_trivial, 1)
P_sign = build_projector(chi_sign, 1)
P_std  = build_projector(chi_standard, 2)


def irrep_content(joint):
    """Returns (frac_trivial, frac_sign, frac_standard) of a [4,4] crystal."""
    v = joint.reshape(dim_sq)
    norm = v.abs().pow(2).sum().item()
    if norm < 1e-20:
        return (1/3, 1/3, 1/3)
    ft = (P_triv @ v).abs().pow(2).sum().item() / norm
    fs = (P_sign @ v).abs().pow(2).sum().item() / norm
    fd = (P_std @ v).abs().pow(2).sum().item() / norm
    return (ft, fs, fd)


# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("PART 1: IRREP CONTENT — SINGLE SEED (seed=42) vs SEED-AVERAGED")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

N_SEEDS = 200

for class_name, members in CONJUGACY_CLASSES.items():
    # Single seed
    single_fracs = []
    for name in members:
        j = Entangler(S3_CORR[name], seed=42).build().joint
        single_fracs.append(irrep_content(j))

    avg_single = tuple(np.mean([f[i] for f in single_fracs]) for i in range(3))

    # Seed-averaged crystal
    avg_joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    count = 0
    for name in members:
        for seed in range(N_SEEDS):
            avg_joint += Entangler(S3_CORR[name], seed=seed).build().joint
            count += 1
    avg_joint /= count
    avg_joint /= avg_joint.reshape(-1).real.sum()
    avg_fracs = irrep_content(avg_joint)

    # Per-seed statistics
    all_fracs = []
    for name in members:
        for seed in range(N_SEEDS):
            j = Entangler(S3_CORR[name], seed=seed).build().joint
            all_fracs.append(irrep_content(j))

    per_seed_avg = tuple(np.mean([f[i] for f in all_fracs]) for i in range(3))
    per_seed_std = tuple(np.std([f[i] for f in all_fracs]) for i in range(3))

    print(f"\n  {class_name}:")
    print(f"    Single (seed=42): triv={avg_single[0]:.4f}, sign={avg_single[1]:.4f}, std={avg_single[2]:.4f}")
    print(f"    Seed-averaged:    triv={avg_fracs[0]:.4f}, sign={avg_fracs[1]:.4f}, std={avg_fracs[2]:.4f}")
    print(f"    Per-seed mean:    triv={per_seed_avg[0]:.4f}±{per_seed_std[0]:.4f}, "
          f"sign={per_seed_avg[1]:.4f}±{per_seed_std[1]:.4f}, "
          f"std={per_seed_avg[2]:.4f}±{per_seed_std[2]:.4f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 2: DEPENDENCE ON ENTANGLER PARAMETERS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Test how n_steps and discount affect the irrep content
corr_012 = S3_CORR["(012)"]
corr_01  = S3_CORR["(01)"]
corr_e   = S3_CORR["e"]

print("\n  3-cycle (012) — varying n_steps (discount=0.3, seed=42):")
for n_steps in [10, 20, 50, 100, 200, 500]:
    j = Entangler(corr_012, seed=42).build(n_steps=n_steps)
    ft, fs, fd = irrep_content(j.joint)
    print(f"    n={n_steps:>4d}: triv={ft:.4f}, sign={fs:.4f}, std={fd:.4f}")

print("\n  3-cycle (012) — varying discount (n_steps=50, seed=42):")
for discount in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8]:
    j = Entangler(corr_012, seed=42).build(discount=discount)
    ft, fs, fd = irrep_content(j.joint)
    print(f"    d={discount:.1f}: triv={ft:.4f}, sign={fs:.4f}, std={fd:.4f}")

print("\n  Transposition (01) — varying n_steps (discount=0.3, seed=42):")
for n_steps in [10, 20, 50, 100, 200, 500]:
    j = Entangler(corr_01, seed=42).build(n_steps=n_steps)
    ft, fs, fd = irrep_content(j.joint)
    print(f"    n={n_steps:>4d}: triv={ft:.4f}, sign={fs:.4f}, std={fd:.4f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 3: THE EXACT IRREP CONTENT OF PURE PERMUTATION MATRICES")
print("(without entangler — the algebraic prediction)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# The PURE permutation matrix P_σ (4×4) has exact irrep content
# that we can compute algebraically

for name, perm in S3_PERMS.items():
    P = permutation_matrix(perm).real  # The 4×4 permutation matrix
    v = P.reshape(dim_sq).to(torch.cfloat)
    ft, fs, fd = irrep_content(v.reshape(MASS_DIM, MASS_DIM))
    print(f"  {name:>6s}: triv={ft:.6f}, sign={fs:.6f}, std={fd:.6f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 4: ALGEBRAIC ANALYSIS — WHY THESE FRACTIONS?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# For the identity permutation P_e = I₄:
# v(e) = vec(I₄) = sum over i of e_i ⊗ e_i
# This is the "identity vector" in the 16D space.
# Under S₃: P_triv(v(e)) projects onto the S₃-invariant part of I₄.
# Since I₄ is S₃-invariant (permutations fix I₄ because they only permute
# among h0,h1,h2 and these appear equally on the diagonal), we expect
# v(e) to be 100% trivial.

# For a transposition, say (01): P_(01) swaps rows and columns 0,1.
# This is NOT S₃-invariant: it treats pair (0,1) differently from (0,2), (1,2).
# The standard representation captures this directionality.

# For a 3-cycle (012): det(P) = +1 (even permutation).
# The sign representation detects parity. But the sign rep on C⁴
# acts as: (P_σ ⊗ P_σ) on vec(M), and the sign character is sgn(σ).
# So 3-cycles have sgn = +1 and project differently than transpositions.

# Let's verify: decompose each permutation matrix into irreps ALGEBRAICALLY
print("  Algebraic decomposition of P_σ ⊗ P_σ action on vec(P_σ):")
print()

for name, perm in S3_PERMS.items():
    P = permutation_matrix(perm)  # 4×4
    v = P.reshape(dim_sq)  # 16D vector

    # Apply each group element and decompose
    orbit = torch.zeros(6, dim_sq, dtype=torch.cfloat)
    for i, (gname, gperm) in enumerate(S3_PERMS.items()):
        G = group_rep_16(gperm)
        orbit[i] = G @ v

    # The orbit of v under S₃ tells us the representation content
    # Project onto characters:
    triv_proj = orbit.mean(0)  # average = trivial projection
    sign_proj = sum((-1 if S3_PERMS[gname] in [[1,0,2],[2,1,0],[0,2,1]] else 1) * orbit[i]
                    for i, gname in enumerate(S3_PERMS)) / 6.0

    # Check: is orbit[0] = orbit[1] = ... ? (i.e., is v S₃-invariant?)
    spread = torch.stack([orbit[i] - orbit[0] for i in range(1, 6)]).abs().max().item()

    print(f"  {name:>6s}: orbit spread = {spread:.6f} "
          f"(0 = fully invariant)")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 5: THE TRIVIAL-SIGN SPLIT FOR 3-CYCLES")
print("Does sign fraction → 1/4 = 1/MASS_DIM as n_steps → ∞?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Track sign fraction for (012) as entangler converges
corr = S3_CORR["(012)"]
print("\n  (012) sign fraction vs n_steps (200 seeds):")
for n_steps in [10, 20, 50, 100, 200, 500, 1000]:
    sign_fracs = []
    for seed in range(200):
        j = Entangler(corr, seed=seed).build(n_steps=n_steps)
        _, fs, _ = irrep_content(j.joint)
        sign_fracs.append(fs)
    mean_s = np.mean(sign_fracs)
    std_s = np.std(sign_fracs)
    print(f"    n={n_steps:>5d}: sign = {mean_s:.6f} ± {std_s:.6f}  "
          f"(1/4 = {0.25:.6f}, 1/3 = {0.3333:.4f})")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 6: THE STANDARD FRACTION FOR TRANSPOSITIONS")
print("Does standard fraction → 1/3 as n_steps → ∞?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

corr = S3_CORR["(01)"]
print("\n  (01) standard fraction vs n_steps (200 seeds):")
for n_steps in [10, 20, 50, 100, 200, 500, 1000]:
    std_fracs = []
    for seed in range(200):
        j = Entangler(corr, seed=seed).build(n_steps=n_steps)
        _, _, fd = irrep_content(j.joint)
        std_fracs.append(fd)
    mean_d = np.mean(std_fracs)
    std_d = np.std(std_fracs)
    print(f"    n={n_steps:>5d}: standard = {mean_d:.6f} ± {std_d:.6f}  "
          f"(1/3 = {1/3:.6f}, K* = {K_STAR:.6f})")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("PART 7: COLLINEARITY VERIFICATION")
print("d(e,t) + d(t,3) = d(e,3)?")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

from solver.algebra import born_fidelity

def fs_dist(j1, j2):
    fid = born_fidelity(j1, j2)
    return math.acos(min(max(fid, -1), 1))

# Build seed-averaged class crystals
class_crystals = {}
for class_name, members in CONJUGACY_CLASSES.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    count = 0
    for name in members:
        for seed in range(N_SEEDS):
            avg += Entangler(S3_CORR[name], seed=seed).build().joint
            count += 1
    avg /= count
    avg /= avg.reshape(-1).real.sum()
    class_crystals[class_name] = avg

d_et = fs_dist(class_crystals["identity"], class_crystals["transpositions"])
d_e3 = fs_dist(class_crystals["identity"], class_crystals["3-cycles"])
d_t3 = fs_dist(class_crystals["transpositions"], class_crystals["3-cycles"])

print(f"\n  d(identity, transpositions) = {d_et:.6f}")
print(f"  d(identity, 3-cycles)       = {d_e3:.6f}")
print(f"  d(transpositions, 3-cycles) = {d_t3:.6f}")
print(f"\n  d(e,t) + d(t,3) = {d_et + d_t3:.6f}")
print(f"  d(e,3)           = {d_e3:.6f}")
print(f"  Difference:        {abs(d_e3 - d_et - d_t3):.6f} ({abs(d_e3 - d_et - d_t3)/d_e3*100:.3f}%)")

# Ratio test: is transpositions at (H-1)/H of the way?
ratio_et = d_et / d_e3
print(f"\n  d(e,t)/d(e,3) = {ratio_et:.6f}  ((H-1)/H = {(H-1)/H:.6f})")
print(f"  d(t,3)/d(e,3) = {d_t3/d_e3:.6f}  (1/H = {1/H:.6f})")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("FINDINGS")
print("=" * 80)
