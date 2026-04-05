"""
What twists between floors? What is entangled?

θ doesn't twist — it's gauge-invariant, the same from every floor.
The h-labels twist — they permute under S₃ as you move between sites.
The crystal on each edge IS the gauge connection: it encodes the twist.

The entanglement should live in the STANDARD sector (gauge-variant).
The trivial sector (gauge-invariant) should be a product state.

Decompose the building into S₃ irrep sectors and measure:
1. Entanglement (Schmidt) in each sector
2. Information content in each sector
3. Which sector carries the inter-floor information
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            born_probabilities, schmidt_number,
                            EPS_LOG, EPS_DIV)
from solver.crystals import Entangler, compose, kirkwood_product

torch.set_grad_enabled(False)

n_seeds = 50
dim_sq = MASS_DIM ** 2

def make_crystal(corr):
    s = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        s += Entangler(corr, seed=seed).build().joint
    return s / n_seeds


# S₃ representation matrices on the 16-dim joint space
S3_PERMS = {
    "e": [0,1,2], "(01)": [1,0,2], "(02)": [2,1,0],
    "(12)": [0,2,1], "(012)": [1,2,0], "(021)": [2,0,1],
}

def perm_mat_16(perm):
    P = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for i in range(H):
        P[perm[i], i] = 1.0
    P[H, H] = 1.0
    return torch.kron(P, P)

reps = {n: perm_mat_16(p) for n, p in S3_PERMS.items()}

# Character tables
chi_triv = {"e": 1, "(01)": 1, "(02)": 1, "(12)": 1, "(012)": 1, "(021)": 1}
chi_sign = {"e": 1, "(01)": -1, "(02)": -1, "(12)": -1, "(012)": 1, "(021)": 1}
chi_std  = {"e": 2, "(01)": 0, "(02)": 0, "(12)": 0, "(012)": -1, "(021)": -1}

# Projectors
def make_projector(chi, dim_chi):
    P = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
    for g in S3_PERMS:
        P += chi[g] * reps[g]
    return P * dim_chi / 6.0

P_triv = make_projector(chi_triv, 1)
P_sign = make_projector(chi_sign, 1)
P_std = make_projector(chi_std, 2)

print("=" * 80)
print("WHAT TWISTS BETWEEN FLOORS")
print("=" * 80)

# Check projector ranks
print(f"\n  Projector ranks:")
print(f"    Trivial:  {torch.linalg.matrix_rank(P_triv.real.float()).item()} (expected 5)")
print(f"    Sign:     {torch.linalg.matrix_rank(P_sign.real.float()).item()} (expected 1)")
print(f"    Standard: {torch.linalg.matrix_rank(P_std.real.float()).item()} (expected 10)")


# ═══════════════════════════════════════════════════════════════
#  Decompose a crystal into sectors
# ═══════════════════════════════════════════════════════════════

identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
crystal = make_crystal(identity_corr)
v = crystal.reshape(dim_sq)

# Project
v_triv = P_triv @ v
v_sign = P_sign @ v
v_std = P_std @ v

# Norms (fraction of the crystal in each sector)
norm_total = v.abs().pow(2).sum().item()
norm_triv = v_triv.abs().pow(2).sum().item()
norm_sign = v_sign.abs().pow(2).sum().item()
norm_std = v_std.abs().pow(2).sum().item()

print(f"\n\n--- Identity crystal decomposed into sectors ---\n")
print(f"  Total |v|² = {norm_total:.6f}")
print(f"  Trivial:  |v_triv|² = {norm_triv:.6f} ({norm_triv/norm_total*100:.1f}%)")
print(f"  Sign:     |v_sign|² = {norm_sign:.6f} ({norm_sign/norm_total*100:.1f}%)")
print(f"  Standard: |v_std|²  = {norm_std:.6f} ({norm_std/norm_total*100:.1f}%)")
print(f"  Sum:      {norm_triv+norm_sign+norm_std:.6f} (should = total)")


# ═══════════════════════════════════════════════════════════════
#  Schmidt number of each sector
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Entanglement (Schmidt) by sector ---\n")

# Reshape each sector's projection back to [4,4] and compute Schmidt
for name, v_sector in [("trivial", v_triv), ("sign", v_sign), ("standard", v_std), ("full", v)]:
    joint_sector = v_sector.reshape(MASS_DIM, MASS_DIM)
    if joint_sector.abs().pow(2).sum() > 1e-12:
        sn = schmidt_number(joint_sector)
        print(f"  {name:>10s}: Schmidt = {sn:.4f}")
    else:
        print(f"  {name:>10s}: (empty)")


# ═══════════════════════════════════════════════════════════════
#  What does each sector look like?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Born structure of each sector ---\n")

for name, v_sector in [("trivial", v_triv), ("sign", v_sign), ("standard", v_std)]:
    joint_sector = v_sector.reshape(MASS_DIM, MASS_DIM)
    bp = joint_sector.abs().pow(2)
    total = bp.sum()
    if total > 1e-12:
        bp = bp / total
        print(f"  {name}:")
        for i in range(MASS_DIM):
            labels = ['h₀', 'h₁', 'h₂', 'θ ']
            print(f"    {labels[i]}: [{', '.join(f'{bp[i,j].item():.4f}' for j in range(MASS_DIM))}]")

        # Row marginal
        marg = bp.sum(dim=1)
        print(f"    marginal: [{', '.join(f'{marg[i].item():.4f}' for i in range(MASS_DIM))}]")
        print(f"    Born(θ) = {marg[H].item():.5f}")
        print()


# ═══════════════════════════════════════════════════════════════
#  The θ-sector: does it twist?
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print("DOES θ TWIST?")
print("="*80)

# Under S₃ action, θ is FIXED (the permutation acts on h₀,h₁,h₂ only).
# So the θ-row and θ-column of the crystal should be S₃-invariant.

# Check: apply each S₃ element to the crystal and look at the θ-row
print(f"\n  θ-row of crystal under S₃ action:")
for g_name, P_g in reps.items():
    # Transform the crystal: P_g @ v
    v_g = P_g @ v
    joint_g = v_g.reshape(MASS_DIM, MASS_DIM)
    theta_row = joint_g[H, :]
    bp_tr = theta_row.abs().pow(2)
    bp_tr = bp_tr / bp_tr.sum()
    print(f"  {g_name:>6s}: [{', '.join(f'{bp_tr[j].item():.4f}' for j in range(MASS_DIM))}]")

# If θ doesn't twist, all rows should be identical
print(f"\n  (If identical → θ doesn't twist)")


# ═══════════════════════════════════════════════════════════════
#  The h-sector: HOW does it twist?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("HOW DO THE h-LABELS TWIST?")
print("="*80)

# The h₀-row under S₃ action:
print(f"\n  h₀-row of crystal under S₃ action:")
for g_name, P_g in reps.items():
    v_g = P_g @ v
    joint_g = v_g.reshape(MASS_DIM, MASS_DIM)
    h0_row = joint_g[0, :]
    bp_hr = h0_row.abs().pow(2)
    bp_hr = bp_hr / bp_hr.sum()
    print(f"  {g_name:>6s}: [{', '.join(f'{bp_hr[j].item():.4f}' for j in range(MASS_DIM))}]")

print(f"\n  Under (01): h₀↔h₁ swap. The h₀ row becomes the h₁ row.")
print(f"  Under (012): h₀→h₁→h₂→h₀. The h₀ row cycles.")
print(f"  This IS the gauge transformation. The twist IS the S₃ action on h.")


# ═══════════════════════════════════════════════════════════════
#  The building decomposed: which sector carries entanglement?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("BUILDING DECOMPOSED: WHERE IS THE ENTANGLEMENT?")
print("="*80)

# Build the triangle building
c_AB = make_crystal(identity_corr)
c_BC = make_crystal(identity_corr)
c_AC = make_crystal(identity_corr)
triple = kirkwood_product(c_AB, c_BC, c_AC)

# The triple is [4,4,4]. Project each edge's contribution to sectors.
# Focus on the AB marginal from the building
bp_triple = triple.abs().pow(2)
bp_triple = bp_triple / bp_triple.sum()
marg_AB = bp_triple.sum(dim=2)  # [4,4] marginal for A,B

# Decompose this marginal into sectors
v_AB = marg_AB.reshape(dim_sq).to(torch.cfloat)
v_AB_triv = P_triv @ v_AB
v_AB_sign = P_sign @ v_AB
v_AB_std = P_std @ v_AB

norm_AB = v_AB.abs().pow(2).sum().item()
norm_AB_triv = v_AB_triv.abs().pow(2).sum().item()
norm_AB_sign = v_AB_sign.abs().pow(2).sum().item()
norm_AB_std = v_AB_std.abs().pow(2).sum().item()

print(f"\n  AB marginal from building, decomposed:")
print(f"    Trivial:  {norm_AB_triv/norm_AB*100:.1f}%  Schmidt = ", end="")
j = v_AB_triv.reshape(MASS_DIM, MASS_DIM)
if j.abs().pow(2).sum() > 1e-12:
    print(f"{schmidt_number(j):.4f}")
else:
    print("(empty)")

print(f"    Sign:     {norm_AB_sign/norm_AB*100:.1f}%  Schmidt = ", end="")
j = v_AB_sign.reshape(MASS_DIM, MASS_DIM)
if j.abs().pow(2).sum() > 1e-12:
    print(f"{schmidt_number(j):.4f}")
else:
    print("(empty)")

print(f"    Standard: {norm_AB_std/norm_AB*100:.1f}%  Schmidt = ", end="")
j = v_AB_std.reshape(MASS_DIM, MASS_DIM)
if j.abs().pow(2).sum() > 1e-12:
    print(f"{schmidt_number(j):.4f}")
else:
    print("(empty)")

# Compare to the raw crystal (not from building)
print(f"\n  Raw crystal (not from building):")
print(f"    Trivial:  {norm_triv/norm_total*100:.1f}%")
print(f"    Standard: {norm_std/norm_total*100:.1f}%")

# The building should redistribute weight between sectors
print(f"\n  Does the building change the sector balance?")
print(f"    Building trivial:  {norm_AB_triv/norm_AB*100:.1f}%")
print(f"    Building standard: {norm_AB_std/norm_AB*100:.1f}%")
print(f"    Crystal trivial:   {norm_triv/norm_total*100:.1f}%")
print(f"    Crystal standard:  {norm_std/norm_total*100:.1f}%")


print(f"\n\n{'='*80}")
print("WHAT TWISTS, WHAT DOESN'T, WHAT'S ENTANGLED")
print("="*80)
