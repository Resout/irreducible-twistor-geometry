"""
The 5 standard doublets and quark masses.

The standard sector has 10 dimensions = 5 copies of V_std (the 2D standard
representation of S₃). Under composition, these 5 doublets have nearly
degenerate eigenvalues that split at ~3% per step (Principle 92).

Questions:
1. How do the 5 doublet eigenvalues split under iterated composition?
2. Do the split eigenvalues, raised to appropriate powers, give quark mass ratios?
3. Does the doublet structure map to up-type/down-type quarks (3 + 2 or 2 + 3)?
4. At n = 35/9 (the Koide exponent), do the 5 eigenvalue pairs give Koide-like
   relations for quark triplets?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import numpy as np
from solver.algebra import H, MASS_DIM, K_STAR, BORN_FLOOR, DELTA
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)

dim_sq = MASS_DIM ** 2

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

def perm_mat_4x4(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return P

reps = {n: torch.kron(perm_mat_4x4(p), perm_mat_4x4(p)) for n, p in S3_PERMS.items()}
chi_d = {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1}

P_std = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
for n in S3_PERMS:
    P_std += chi_d[n] * reps[n]
P_std *= 2 / 6.0


def build_T_B(B_joint):
    T = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for col in range(dim_sq):
        A = torch.zeros(dim_sq, dtype=torch.cfloat)
        A[col] = 1.0
        result = compose(A.reshape(MASS_DIM, MASS_DIM), B_joint)
        T[:, col] = result.reshape(dim_sq)
    return T


# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("STANDARD SECTOR EIGENVALUES: SEED-AVERAGED")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

N_SEEDS = 100

# For each conjugacy class, build seed-averaged T_B and extract standard spectrum
for class_name, members in [("identity", ["e"]),
                             ("transpositions", ["(01)", "(02)", "(12)"]),
                             ("3-cycles", ["(012)", "(021)"])]:
    all_evals = []

    for name in members:
        T_avg = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
        for seed in range(N_SEEDS):
            ent = Entangler(S3_CORR[name], seed=seed).build()
            T_avg += build_T_B(ent.joint)
        T_avg /= N_SEEDS

        restricted = P_std @ T_avg @ P_std
        evals = torch.linalg.eigvals(restricted)
        nonzero = evals[evals.abs() > 1e-8]
        mags = nonzero.abs().sort(descending=True).values
        all_evals.append(mags)

    # Average across elements in class
    min_len = min(len(e) for e in all_evals)
    avg_evals = torch.stack([e[:min_len] for e in all_evals]).mean(dim=0)

    print(f"\n  {class_name} ({len(members)} elements):")
    print(f"  {'idx':>4s} {'|λ|':>10s} {'|λ|/|λ₀|':>10s} {'gap':>10s}")

    prev = None
    for i in range(min(10, len(avg_evals))):
        mag = avg_evals[i].item()
        ratio = mag / avg_evals[0].item()
        gap = f"{-math.log(mag/prev):.4f}" if prev and mag > 0 else ""
        prev = mag
        print(f"  {i:>4d} {mag:>10.6f} {ratio:>10.6f} {gap:>10s}")

    # Identify doublets (pairs with < 1% splitting)
    print(f"\n  Doublet structure:")
    i = 0
    doublet_idx = 0
    while i < min(10, len(avg_evals)):
        mag_i = avg_evals[i].item()
        if i + 1 < len(avg_evals):
            mag_next = avg_evals[i+1].item()
            split = abs(mag_i - mag_next) / mag_i * 100
            if split < 5:
                avg_doublet = (mag_i + mag_next) / 2
                print(f"    Doublet {doublet_idx}: |λ| = {avg_doublet:.6f} (split = {split:.2f}%)")
                i += 2
                doublet_idx += 1
                continue
        print(f"    Singlet: |λ| = {mag_i:.6f}")
        i += 1
        doublet_idx += 1


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("DOUBLET SPLITTING UNDER ITERATED COMPOSITION")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Build seed-averaged identity crystal
id_joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for seed in range(N_SEEDS):
    id_joint += Entangler(S3_CORR["e"], seed=seed).build().joint
id_joint /= N_SEEDS

# Iterate self-composition and track standard sector spectrum
current = id_joint.clone()
print(f"\n  {'step':>5s}  {'λ₀':>8s} {'λ₁':>8s} {'λ₂':>8s} {'λ₃':>8s} {'λ₄':>8s}"
      f"  {'λ₅':>8s} {'λ₆':>8s} {'λ₇':>8s} {'λ₈':>8s} {'λ₉':>8s}")

for step in range(8):
    T = build_T_B(current)
    restricted = P_std @ T @ P_std
    evals = torch.linalg.eigvals(restricted)
    nonzero = evals[evals.abs() > 1e-8]
    mags = nonzero.abs().sort(descending=True).values

    vals = " ".join(f"{mags[i].item():>8.5f}" if i < len(mags) else f"{'':>8s}"
                    for i in range(10))
    print(f"  {step:>5d}  {vals}")

    if step < 7:
        current = compose(current, id_joint)


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("DOUBLET RATIOS AT n = 35/9 (KOIDE EXPONENT)")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Get the standard sector eigenvalues at step 0
T = build_T_B(id_joint)
restricted = P_std @ T @ P_std
evals = torch.linalg.eigvals(restricted)
nonzero = evals[evals.abs() > 1e-8]
mags = nonzero.abs().sort(descending=True).values

n_koide = 35/9

# Group into doublets
doublets = []
i = 0
while i < len(mags):
    if i + 1 < len(mags) and abs(mags[i].item() - mags[i+1].item()) / mags[i].item() < 0.05:
        doublets.append((mags[i].item() + mags[i+1].item()) / 2)
        i += 2
    else:
        doublets.append(mags[i].item())
        i += 1

print(f"\n  {len(doublets)} doublet averages:")
for d_idx, d_val in enumerate(doublets):
    d_mass = d_val ** n_koide
    print(f"    D{d_idx}: |λ| = {d_val:.6f}, |λ|^(35/9) = {d_mass:.6e}")

# Mass ratios between doublets
if len(doublets) >= 3:
    print(f"\n  Doublet mass ratios at n = 35/9:")
    for i in range(len(doublets)):
        for j in range(i+1, len(doublets)):
            ratio = doublets[i]**n_koide / doublets[j]**n_koide
            print(f"    D{i}/D{j} = {ratio:.2f}")

# Compare with quark mass ratios
print(f"\n  Physical quark mass ratios:")
quarks = {"u": 2.16, "d": 4.67, "s": 93.4, "c": 1270, "b": 4180, "t": 172500}
for pair in [("t", "b"), ("t", "c"), ("b", "c"), ("b", "s"), ("c", "s"), ("s", "d"), ("d", "u")]:
    q1, q2 = pair
    print(f"    m_{q1}/m_{q2} = {quarks[q1]/quarks[q2]:.2f}")


# ═══════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("KOIDE FOR DOUBLET TRIPLETS")
print("=" * 80)
# ═══════════════════════════════════════════════════════════════

# Take triplets from the 5 doublets and check Koide
if len(doublets) >= 3:
    from itertools import combinations

    print(f"\n  Koide Q for all triplets of doublet masses (at n = 35/9):")
    for combo in combinations(range(len(doublets)), 3):
        m = [doublets[i]**n_koide for i in combo]
        Q = sum(m) / (sum(math.sqrt(mi) for mi in m))**2
        label = f"D{''.join(str(c) for c in combo)}"
        marker = " ←" if abs(Q - 2/3) < 0.05 else ""
        print(f"    {label}: Q = {Q:.6f}{marker}")


print(f"\n\n{'='*80}")
print("SUMMARY")
print("=" * 80)
