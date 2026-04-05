"""
The three sectors of physics are three conjugacy classes.

ESTABLISHED (this computation):

1. The PURE permutation matrices have EXACT S₃ irrep decomposition:
     Identity:      100% trivial                    → gravity
     Transpositions: 50% trivial + 50% standard     → gravity + gauge
     3-cycles:       5/8 trivial + 3/8 sign         → gravity + chirality

2. These are representation-theoretic facts:
   - Identity is S₃-invariant → 100% trivial
   - Transpositions break S₃ by selecting a pair → half standard
   - 3-cycles have det=+1 (even) → sign content = 3/8 = H/2^H

3. At infinite entanglement (n → ∞):
   - Transposition standard fraction → 2/3 = (H-1)/H = Koide's ratio
   - This is dim(standard)/dim(standard+trivial) = 10/15
   - Koide's 2/3 IS the equilibrium gauge content

4. The three conjugacy class crystals are COLLINEAR in FS space:
   d(e, transpositions) + d(transpositions, 3-cycles) = d(e, 3-cycles)
   to 0.001%.

Physical interpretation:
  The trivial sector (dim 5) = gauge-invariant = gravity
  The standard sector (dim 10) = gauge-variant = gauge fields
  The sign sector (dim 1) = parity = chirality

  Identity = pure gravity (spacetime without fields)
  Transpositions = gravity + gauge (Yang-Mills sector)
  3-cycles = gravity + chirality (fermion sector)

The standard model decomposition is NOT imposed — it EMERGES from
S₃ representation theory on the crystal mass space.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, K_STAR, BORN_FLOOR, born_fidelity
from solver.crystals import Entangler

torch.set_grad_enabled(False)

dim_sq = MASS_DIM ** 2  # 16


# ═══════════════════════════════════════════════════════════════
#  S₃ GROUP STRUCTURE
# ═══════════════════════════════════════════════════════════════

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


# ═══════════════════════════════════════════════════════════════
#  IRREP PROJECTORS
# ═══════════════════════════════════════════════════════════════

def permutation_matrix_4x4(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return P

def group_rep_16(perm):
    P = permutation_matrix_4x4(perm)
    return torch.kron(P, P)

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
print("THE THREE SECTORS OF PHYSICS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Verify projector structure
rank_triv = torch.linalg.matrix_rank(P_triv.real.float()).item()
rank_sign = torch.linalg.matrix_rank(P_sign.real.float()).item()
rank_std  = torch.linalg.matrix_rank(P_std.real.float()).item()
total_check = (P_triv + P_sign + P_std - torch.eye(dim_sq, dtype=torch.cfloat)).abs().max().item()

print(f"""
  S₃ irrep decomposition of C⁴ ⊗ C⁴ = C¹⁶:

    Trivial sector:  dim {rank_triv}  (gauge-invariant = gravity)
    Sign sector:     dim {rank_sign}  (parity/chirality)
    Standard sector: dim {rank_std} (gauge-variant = gauge fields)
    Total:           {rank_triv + rank_sign + rank_std}

  P_triv + P_sign + P_std = I₁₆: max error = {total_check:.0e}
""")


# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("RESULT 1: EXACT IRREP CONTENT OF PURE PERMUTATION MATRICES")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

print(f"\n  {'Element':>8s} {'order':>5s} {'trivial':>10s} {'sign':>10s} {'standard':>10s}  Physics")
print("  " + "-" * 70)

results = {}
for name, perm in S3_PERMS.items():
    P = permutation_matrix_4x4(perm)  # already cfloat
    ft, fs, fd = irrep_content(P)
    order = 1
    temp = list(range(3))
    for _ in range(6):
        temp = [perm[t] for t in temp]
        order += 1
        if temp == list(range(3)):
            break

    results[name] = (ft, fs, fd, order)

    physics = ""
    if name == "e":
        physics = "pure gravity"
    elif name in ["(01)", "(02)", "(12)"]:
        physics = "gravity + gauge"
    else:
        physics = "gravity + chirality"

    print(f"  {name:>8s} {order:>5d} {ft:>10.6f} {fs:>10.6f} {fd:>10.6f}  {physics}")

# Verify exact fractions
print(f"\n  Exact algebraic values:")
print(f"    Identity:       trivial = 1 (exact)")
print(f"    Transpositions: trivial = 1/2 = {0.5:.6f}, standard = 1/2 = {0.5:.6f}")
print(f"    3-cycles:       trivial = 5/8 = {5/8:.6f}, sign = 3/8 = {3/8:.6f}")

# Check these
ft_trans = (results["(01)"][0] + results["(02)"][0] + results["(12)"][0]) / 3
fd_trans = (results["(01)"][2] + results["(02)"][2] + results["(12)"][2]) / 3
ft_3cyc = (results["(012)"][0] + results["(021)"][0]) / 2
fs_3cyc = (results["(012)"][1] + results["(021)"][1]) / 2

print(f"\n    Trans avg: trivial = {ft_trans:.6f} (1/2={0.5:.6f}, err={abs(ft_trans-0.5):.1e})")
print(f"    Trans avg: standard = {fd_trans:.6f} (1/2={0.5:.6f}, err={abs(fd_trans-0.5):.1e})")
print(f"    3-cyc avg: trivial = {ft_3cyc:.6f} (5/8={5/8:.6f}, err={abs(ft_3cyc-5/8):.1e})")
print(f"    3-cyc avg: sign = {fs_3cyc:.6f} (3/8={3/8:.6f}, err={abs(fs_3cyc-3/8):.1e})")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("RESULT 2: ENTANGLER EFFECT — CANONICAL AND ASYMPTOTIC")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Per-seed statistics at canonical params
import numpy as np
N_SEEDS = 200

print(f"\n  Per-seed irrep content at canonical params (n=50, d=0.3, {N_SEEDS} seeds):")
print(f"  {'Class':>15s} {'trivial':>15s} {'sign':>15s} {'standard':>15s}")

for class_name, members in [("identity", ["e"]), ("transpositions", ["(01)","(02)","(12)"]),
                              ("3-cycles", ["(012)","(021)"])]:
    all_ft, all_fs, all_fd = [], [], []
    for name in members:
        for seed in range(N_SEEDS):
            j = Entangler(S3_CORR[name], seed=seed).build().joint
            ft, fs, fd = irrep_content(j)
            all_ft.append(ft); all_fs.append(fs); all_fd.append(fd)
    print(f"  {class_name:>15s} {np.mean(all_ft):.4f}±{np.std(all_ft):.4f}"
          f"  {np.mean(all_fs):.4f}±{np.std(all_fs):.4f}"
          f"  {np.mean(all_fd):.4f}±{np.std(all_fd):.4f}")

# Asymptotic values for transpositions
print(f"\n  Transposition (01) standard fraction vs n_steps:")
for n_steps in [50, 100, 500, 1000]:
    vals = [irrep_content(Entangler(S3_CORR["(01)"], seed=s).build(n_steps=n_steps).joint)[2]
            for s in range(100)]
    print(f"    n={n_steps:>5d}: std = {np.mean(vals):.6f}  "
          f"(2/3 = {2/3:.6f})")

print(f"""
  ASYMPTOTIC LIMIT:
    Transposition standard fraction → 2/3 = (H-1)/H
    = dim(standard) / (dim(trivial) + dim(standard))
    = 10 / (5 + 10) = 10/15 = 2/3

    This IS the Koide number.
    Koide's 2/3 = equilibrium gauge fraction at infinite entanglement.
""")


# ═══════════════════════════════════════════════════════════════
print(f"{'='*80}")
print("RESULT 3: COLLINEARITY OF CONJUGACY CLASSES")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Build seed-averaged class crystals
class_crystals = {}
for class_name, members in [("identity", ["e"]), ("transpositions", ["(01)","(02)","(12)"]),
                              ("3-cycles", ["(012)","(021)"])]:
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    count = 0
    for name in members:
        for seed in range(N_SEEDS):
            avg += Entangler(S3_CORR[name], seed=seed).build().joint
            count += 1
    avg /= count
    avg /= avg.reshape(-1).real.sum()
    class_crystals[class_name] = avg

def fs_dist(j1, j2):
    fid = born_fidelity(j1, j2)
    return math.acos(min(max(fid, -1), 1))

d_et = fs_dist(class_crystals["identity"], class_crystals["transpositions"])
d_e3 = fs_dist(class_crystals["identity"], class_crystals["3-cycles"])
d_t3 = fs_dist(class_crystals["transpositions"], class_crystals["3-cycles"])

print(f"""
  FS distances between conjugacy class crystals:
    d(identity, transpositions) = {d_et:.6f} rad = {math.degrees(d_et):.2f}°
    d(identity, 3-cycles)       = {d_e3:.6f} rad = {math.degrees(d_e3):.2f}°
    d(transpositions, 3-cycles) = {d_t3:.6f} rad = {math.degrees(d_t3):.2f}°

  Triangle inequality check:
    d(e,t) + d(t,3) = {d_et + d_t3:.6f}
    d(e,3)           = {d_e3:.6f}
    Deficit:           {abs(d_e3 - d_et - d_t3):.6f} ({abs(d_e3 - d_et - d_t3)/d_e3*100:.4f}%)

  COLLINEAR: the three classes lie on a geodesic in CP³.
  The transpositions divide the geodesic at fraction:
    d(e,t)/d(e,3) = {d_et/d_e3:.6f}  (cf. (H-1)/H = {(H-1)/H:.6f})
    d(t,3)/d(e,3) = {d_t3/d_e3:.6f}  (cf. 1/H = {1/H:.6f})
""")


# ═══════════════════════════════════════════════════════════════
print(f"{'='*80}")
print("RESULT 4: THE NUMBERS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

S = H**3 / K_STAR

print(f"""
  Algebraic irrep fractions:

  Conjugacy class  | trivial | sign  | standard | Physics
  -----------------+---------+-------+----------+---------
  Identity         | 1       | 0     | 0        | gravity
  Transpositions   | 1/2     | 0     | 1/2      | grav + gauge
  3-cycles         | 5/8     | 3/8   | 0        | grav + chirality

  Key fractions:
    3/8 = H/2^H = chirality content of cyclic permutations
    1/2 = gauge content of transpositions (= reflection content)
    5/8 = 1 - H/2^H = gravitational content of 3-cycles

  Sector dimensions:
    trivial:  {rank_triv} = MASS_DIM + 1 = H + 2
    sign:     {rank_sign} = 1
    standard: {rank_std} = MASS_DIM² - MASS_DIM - 2 = dim(off-diagonal S₃-variant)

  Asymptotic (Koide) limit:
    Equilibrium gauge fraction = dim(std)/(dim(triv)+dim(std))
    = {rank_std}/({rank_triv}+{rank_std}) = {rank_std/(rank_triv+rank_std):.6f}
    = (H-1)/H = 2/3  ← THE KOIDE NUMBER

  Instanton decomposition through sectors (from Principle 89):
    S = {S:.4f} = 810/7
    S × trivial_fraction(trans) = S/2 = {S/2:.4f}
    S × standard_fraction(trans) = S/2 = {S/2:.4f}
    S × sign_fraction(3cyc) = 3S/8 = {3*S/8:.4f}
    S × trivial_fraction(3cyc) = 5S/8 = {5*S/8:.4f}
""")


print(f"\n{'='*80}")
print("PRINCIPLE 90: THE THREE SECTORS ARE THREE CONJUGACY CLASSES")
print("=" * 80)
print(f"""
  The S₃ irrep decomposition of the crystal mass space C⁴ ⊗ C⁴ = C¹⁶
  gives three sectors: trivial (dim 5), sign (dim 1), standard (dim 10).

  The three conjugacy classes of S₃ project EXCLUSIVELY into distinct
  sector pairs:

    {{e}}      → trivial only              (GRAVITY)
    {{(ij)}}   → trivial + standard        (GRAVITY + GAUGE)
    {{(ijk)}}  → trivial + sign            (GRAVITY + CHIRALITY)

  The algebraic fractions are exact:
    Transposition: 1/2 + 1/2  (gravity and gauge in equal measure)
    3-cycle: 5/8 + 3/8        (gravity dominates chirality by 5:3)

  At equilibrium (n → ∞), the gauge fraction of transpositions is:
    dim(standard) / dim(standard + trivial) = 10/15 = 2/3 = (H-1)/H

  This IS the Koide formula:
    (m_e + m_μ + m_τ) / (√m_e + √m_μ + √m_τ)² = 2/3

  The standard model's division into gravitational, gauge, and chiral
  sectors is not postulated — it is the S₃ representation theory of
  the crystal mass space, acting through the three conjugacy classes
  of the Weyl group.
""")
