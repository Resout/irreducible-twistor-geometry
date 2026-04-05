"""
Cartan matrix from DS conflict between crystals.

The Cartan matrix A_{ij} = 2⟨α_i,α_j⟩/⟨α_j,α_j⟩ encodes root angles.
For A₂: A = [[2,-1],[-1,2]], meaning cos(angle) = -1/2 (120°).

Principle 44 derived A₁₂ = -0.97 from standard-sector inner products.
Can we get the same from the CONFLICT K between crystals?

When two crystals are DS-combined, the conflict K measures their
disagreement. For two identical crystals: K = 0 (perfect agreement).
For orthogonal crystals: K = K_max. For anti-correlated: K = ???

If K between transposition crystals encodes their root angle,
then: K_{(01),(12)} should give cos(120°) = -1/2.

The conflict K lives in [0, 1). Mapping K → cos(angle):
  cos(angle) = 1 - 2K/K_max  (maps [0, K_max] → [1, -1])
or some other mapping.

Questions:
1. What is K between each pair of S₃ crystals?
2. Is there a mapping K → cos(angle) that gives the Cartan matrix?
3. Does K respect conjugacy classes?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            born_probabilities, born_fidelity, ds_combine,
                            schmidt_number)
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)

S3_CORRS = {
    "e":     torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

n_seeds = 50

# Build seed-averaged crystals
refs = {}
for name, corr in S3_CORRS.items():
    s = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        s += Entangler(corr, seed=seed).build().joint
    refs[name] = s / n_seeds


print("=" * 80)
print("CARTAN MATRIX FROM DS CONFLICT")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Conflict K between all pairs of S₃ crystals
# ═══════════════════════════════════════════════════════════════

print(f"\n--- DS conflict K between crystal marginals ---\n")

# Use the row marginal (4D mass) for DS combination
marginals = {}
for name, joint in refs.items():
    # Row marginal: sum |m_{ij}|² over j, then take square root for mass
    bp = born_probabilities(joint).reshape(MASS_DIM, MASS_DIM)
    marg = bp.sum(dim=1)
    # Convert to complex mass with L1=1
    m = marg.sqrt().to(torch.cfloat)
    m = m / m.real.sum()
    marginals[name] = m

# Compute K for all pairs
names = list(S3_CORRS.keys())
print(f"  {'':>8s}", end="")
for j in names:
    print(f" {j:>7s}", end="")
print()

K_matrix = {}
for i in names:
    print(f"  {i:>8s}", end="")
    for j in names:
        _, K = ds_combine(marginals[i].unsqueeze(0), marginals[j].unsqueeze(0))
        k_val = K.abs().item() if K.is_complex() else K.item()
        K_matrix[(i,j)] = k_val
        print(f" {k_val:>7.4f}", end="")
    print()


# ═══════════════════════════════════════════════════════════════
#  Conflict K from JOINT mass (16D DS combination)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- DS conflict K between JOINT masses (16D) ---\n")

print(f"  {'':>8s}", end="")
for j in names:
    print(f" {j:>7s}", end="")
print()

K_joint = {}
for i in names:
    print(f"  {i:>8s}", end="")
    for j in names:
        flat_i = refs[i].reshape(-1).to(torch.cfloat)
        flat_j = refs[j].reshape(-1).to(torch.cfloat)
        # Normalize to proper mass functions
        flat_i = flat_i / flat_i.real.sum()
        flat_j = flat_j / flat_j.real.sum()
        _, K = ds_combine(flat_i.unsqueeze(0), flat_j.unsqueeze(0))
        k_val = K.abs().item() if K.is_complex() else K.item()
        K_joint[(i,j)] = k_val
        print(f" {k_val:>7.4f}", end="")
    print()


# ═══════════════════════════════════════════════════════════════
#  Born fidelity between all pairs
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Born fidelity between joint masses ---\n")

print(f"  {'':>8s}", end="")
for j in names:
    print(f" {j:>7s}", end="")
print()

BF_matrix = {}
for i in names:
    print(f"  {i:>8s}", end="")
    for j in names:
        bf = born_fidelity(refs[i], refs[j])
        BF_matrix[(i,j)] = bf
        print(f" {bf:>7.4f}", end="")
    print()


# ═══════════════════════════════════════════════════════════════
#  Extract Cartan-like matrix from K or BF
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CARTAN-LIKE MATRICES FROM DIFFERENT INNER PRODUCTS")
print("="*80)

# For the Cartan matrix, we need:
# A_{ij} = 2⟨α_i, α_j⟩/⟨α_j, α_j⟩ = -1 for i≠j (A₂)
# So ⟨α_i, α_j⟩/⟨α_j, α_j⟩ = -1/2 for i≠j (cos 120°)

trans = ["(01)", "(02)", "(12)"]

# Method 1: K as "distance" → convert to inner product
# cos(angle) = 1 - 2K/K_self or cos(angle) = 1 - K/K_max
print(f"\n  From DS conflict K (4D marginals):")
for i in trans:
    for j in trans:
        if i >= j:
            continue
        K_ij = K_matrix[(i,j)]
        K_ii = K_matrix[(i,i)]
        K_jj = K_matrix[(j,j)]
        print(f"    K({i},{j}) = {K_ij:.6f}, K({i},{i}) = {K_ii:.6f}, K({j},{j}) = {K_jj:.6f}")

print(f"\n  From DS conflict K (16D joints):")
for i in trans:
    for j in trans:
        if i >= j:
            continue
        K_ij = K_joint[(i,j)]
        K_ii = K_joint[(i,i)]
        K_jj = K_joint[(j,j)]
        print(f"    K({i},{j}) = {K_ij:.6f}, K({i},{i}) = {K_ii:.6f}, K({j},{j}) = {K_jj:.6f}")

# Method 2: Born fidelity as inner product
# BF ∈ [0,1], with BF=1 for identical crystals
# cos(angle) = (2*BF - 1) maps BF to [-1, 1]
print(f"\n  Born fidelity → Cartan:")
print(f"    Transposition pairs:")
for i in range(len(trans)):
    for j in range(i+1, len(trans)):
        bf = BF_matrix[(trans[i], trans[j])]
        bf_self_i = BF_matrix[(trans[i], trans[i])]
        bf_self_j = BF_matrix[(trans[j], trans[j])]
        # Normalized inner product
        ip = bf / (bf_self_i * bf_self_j)**0.5
        cartan = 2 * ip
        print(f"    BF({trans[i]},{trans[j]}) = {bf:.4f}, normalized = {ip:.4f}, "
              f"A_ij = {cartan:.4f} (target: -1 impossible from BF > 0)")

# Method 3: composition Schmidt ratio as inner product
print(f"\n  Composition Schmidt ratio:")
for i in range(len(trans)):
    for j in range(i+1, len(trans)):
        comp = compose(refs[trans[i]], refs[trans[j]])
        sn = schmidt_number(comp)
        sn_i = schmidt_number(refs[trans[i]])
        sn_j = schmidt_number(refs[trans[j]])
        ratio = sn / (sn_i * sn_j)**0.5
        print(f"    Schmidt({trans[i]}∘{trans[j]}) / √(S_i S_j) = {ratio:.4f}")

# Method 4: 1 - 2*BF as inner product (maps to [-1, 1])
print(f"\n  1 - 2*BF (maps BF to [-1, 1]):")
for i in range(len(trans)):
    for j in range(i+1, len(trans)):
        bf = BF_matrix[(trans[i], trans[j])]
        val = 1 - 2*bf
        print(f"    1 - 2*BF({trans[i]},{trans[j]}) = {val:.4f} (target cos(120°) = -0.500)")


# ═══════════════════════════════════════════════════════════════
#  K between identity and each S₃ element
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- K(e, g) for each g ∈ S₃ (distance from identity) ---\n")
for name in names:
    k4 = K_matrix[("e", name)]
    k16 = K_joint[("e", name)]
    bf = BF_matrix[("e", name)]
    print(f"  {name:>6s}: K_4D = {k4:.4f}, K_16D = {k16:.4f}, BF = {bf:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Does K respect conjugacy classes?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Does K respect conjugacy classes? ---\n")

# Within-class K:
print(f"  Within transpositions:")
for i in range(len(trans)):
    for j in range(i+1, len(trans)):
        print(f"    K({trans[i]}, {trans[j]}) = {K_joint[(trans[i], trans[j])]:.4f}")

cycles = ["(012)", "(021)"]
print(f"\n  Within 3-cycles:")
print(f"    K((012), (021)) = {K_joint[('(012)', '(021)')]:.4f}")

print(f"\n  Cross-class (trans → 3-cycle):")
for t in trans:
    for c in cycles:
        print(f"    K({t}, {c}) = {K_joint[(t, c)]:.4f}")


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
