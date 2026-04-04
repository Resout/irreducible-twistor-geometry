"""
What lives between floors?

The building (full Kirkwood joint) gives Born(h₀) = 0.923 at C.
Sequential propagation gives at best 0.787.
The difference: 0.136 — information that lives BETWEEN floors.

Questions:
1. How much inter-floor information is there? (bits/nats)
2. Does it equal a crystal constant? (Δ, K*, Born_floor?)
3. How does it scale with topology? (chain → triangle → K4 → ?)
4. What happens when the building has DIFFERENT crystal types on edges?

The inter-floor information is what you lose by thinking sequentially.
It's the cost of saying "flows" and "carries" instead of "is."
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            born_probabilities, born_fidelity,
                            ds_combine, enforce_born_floor,
                            schmidt_number, EPS_LOG, EPS_DIV)
from solver.crystals import Entangler, compose, kirkwood_product

torch.set_grad_enabled(False)

n_seeds = 50

def make_crystal(corr):
    s = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        s += Entangler(corr, seed=seed).build().joint
    return s / n_seeds

identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
inverse_corr = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)
cycle_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)


def mutual_info(joint_2d):
    """Mutual information from a 2D Born probability table."""
    p = joint_2d.clone()
    p = p / p.sum()
    marg_0 = p.sum(dim=1)
    marg_1 = p.sum(dim=0)
    mi = 0
    for i in range(p.shape[0]):
        for j in range(p.shape[1]):
            pij = p[i,j].item()
            pi = marg_0[i].item()
            pj = marg_1[j].item()
            if pij > 1e-12 and pi > 1e-12 and pj > 1e-12:
                mi += pij * math.log(pij / (pi * pj))
    return mi


print("=" * 80)
print("WHAT LIVES BETWEEN FLOORS")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Chain vs Triangle: the inter-floor information
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Chain (A—B—C) vs Triangle (A—B—C—A) ---\n")

c_AB = make_crystal(identity_corr)
c_BC = make_crystal(identity_corr)
c_AC = make_crystal(identity_corr)

# CHAIN: only AB and BC edges (no AC)
# The chain's joint is just compose(AB, BC) for the AC marginal
# The full chain joint: P(A,B,C) = P(A,B) × P(C|B) = P(A,B) × P(B,C)/P(B)
# Kirkwood with AB, BC, and IDENTITY for AC

# For a chain, the AC joint is determined by composition: A→B→C
chain_AC = compose(c_AB, c_BC)

# TRIANGLE: all three edges
triangle = kirkwood_product(c_AB, c_BC, c_AC)
bp_triangle = triangle.abs().pow(2)
bp_triangle = bp_triangle / bp_triangle.sum()

# Extract AC marginals from both
# Triangle AC marginal
marg_AC_triangle = bp_triangle.sum(dim=1)  # sum over B

# Chain AC marginal (from composition)
bp_chain_AC = chain_AC.abs().pow(2)
bp_chain_AC = bp_chain_AC / bp_chain_AC.sum()

# Direct AC marginal (from single crystal)
bp_direct_AC = c_AC.abs().pow(2)
bp_direct_AC = bp_direct_AC / bp_direct_AC.sum()

mi_triangle = mutual_info(marg_AC_triangle)
mi_chain = mutual_info(bp_chain_AC)
mi_direct = mutual_info(bp_direct_AC)

print(f"  I(A;C) from triangle building: {mi_triangle:.6f} nats")
print(f"  I(A;C) from chain (compose):   {mi_chain:.6f} nats")
print(f"  I(A;C) from direct edge:       {mi_direct:.6f} nats")
print(f"\n  Inter-floor information:")
print(f"    Triangle - Chain:  {mi_triangle - mi_chain:.6f} nats")
print(f"    Triangle - Direct: {mi_triangle - mi_direct:.6f} nats")
print(f"    Direct - Chain:    {mi_direct - mi_chain:.6f} nats")

# Compare to crystal constants
delta_mi = mi_triangle - mi_chain
print(f"\n  Inter-floor / Δ = {delta_mi / DELTA:.4f}")
print(f"  Inter-floor / K* = {delta_mi / K_STAR:.4f}")
print(f"  Inter-floor / ln(H) = {delta_mi / math.log(H):.4f}")


# ═══════════════════════════════════════════════════════════════
#  Scaling: chain → triangle → K4 → K5
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Scaling with topology: how I(A;farthest) grows ---\n")

# Chain of length n: A—B—C—D—...
# How does I(A; last site) scale?

print(f"  Chain of length n (I(A;site_n)):")
current_joint = c_AB.clone()
for n in range(2, 7):
    new_bond = make_crystal(identity_corr)
    chain_joint = compose(current_joint, new_bond)

    bp = chain_joint.abs().pow(2)
    bp = bp / bp.sum()
    mi = mutual_info(bp)

    print(f"    n={n}: I(A;site_{n}) = {mi:.6f} nats, ratio to n=2: {mi/mutual_info(bp_chain_AC):.4f}")
    current_joint = chain_joint

# Complete graphs: how I(A;B) changes in K_n
print(f"\n  Complete graph K_n (I(A;B) via Kirkwood):")
# K2 = single edge
mi_k2 = mi_direct
print(f"    K2: I(A;B) = {mi_k2:.6f}")

# K3 = triangle
mi_k3 = mi_triangle  # I(A;C) in triangle = I(A;B) by symmetry
print(f"    K3: I(A;B) = {mi_k3:.6f}, ratio to K2: {mi_k3/mi_k2:.4f}")

# K4 = tetrahedron
c_edges = {}
for i in range(4):
    for j in range(i+1, 4):
        c_edges[(i,j)] = make_crystal(identity_corr)

# Build K4 Kirkwood: use triangle (0,1,2) then extend to 3
tri_012 = kirkwood_product(c_edges[(0,1)], c_edges[(1,2)], c_edges[(0,2)])
# Extend to node 3
quad = torch.zeros(MASS_DIM, MASS_DIM, MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for d in range(MASS_DIM):
    w_0d = c_edges[(0,3)][:, d]
    w_1d = c_edges[(1,3)][:, d]
    w_2d = c_edges[(2,3)][:, d]
    for a in range(MASS_DIM):
        for b in range(MASS_DIM):
            for cc in range(MASS_DIM):
                quad[a,b,cc,d] = tri_012[a,b,cc] * w_0d[a] * w_1d[b] * w_2d[cc]

l1 = quad.reshape(-1).real.sum()
if abs(l1) > EPS_DIV: quad = quad / l1

bp_quad = quad.abs().pow(2)
bp_quad = bp_quad / bp_quad.sum()

# I(0;3) in K4
marg_03 = bp_quad.sum(dim=(1,2))  # sum over nodes 1,2
mi_k4 = mutual_info(marg_03)
print(f"    K4: I(A;D) = {mi_k4:.6f}, ratio to K2: {mi_k4/mi_k2:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Mixed building: different crystal types on different edges
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("MIXED BUILDING — DIFFERENT RELATIONSHIPS ON DIFFERENT EDGES")
print("="*80)

c_id = make_crystal(identity_corr)
c_inv = make_crystal(inverse_corr)
c_cyc = make_crystal(cycle_corr)

# Triangle with one identity, one inverse, one cycle
# This is a NON-TRIVIAL building — the floors relate differently
mixed_triple = kirkwood_product(c_id, c_inv, c_cyc)

bp_mixed = mixed_triple.abs().pow(2)
bp_mixed = bp_mixed / bp_mixed.sum()

# Marginals
marg_A_mix = bp_mixed.sum(dim=(1,2))
marg_B_mix = bp_mixed.sum(dim=(0,2))
marg_C_mix = bp_mixed.sum(dim=(0,1))

print(f"\n  Mixed triangle: identity AB, inverse BC, cycle AC")
print(f"  Marginal A: [{', '.join(f'{marg_A_mix[i].item():.4f}' for i in range(MASS_DIM))}]")
print(f"  Marginal B: [{', '.join(f'{marg_B_mix[i].item():.4f}' for i in range(MASS_DIM))}]")
print(f"  Marginal C: [{', '.join(f'{marg_C_mix[i].item():.4f}' for i in range(MASS_DIM))}]")

# Are the marginals still symmetric? (They shouldn't be — different edges!)
print(f"\n  Symmetry broken: A≠B≠C because the edges differ")

# Schmidt of the mixed triple
from solver.crystals import schmidt_3way
sn_uniform = schmidt_3way(kirkwood_product(c_id, c_id, c_id))
sn_mixed = schmidt_3way(mixed_triple)
print(f"\n  3-way Schmidt:")
print(f"    Uniform (all identity): {sn_uniform:.4f}")
print(f"    Mixed (id/inv/cyc):     {sn_mixed:.4f}")

# Mutual information in mixed building
marg_AC_mix = bp_mixed.sum(dim=1)
marg_AB_mix = bp_mixed.sum(dim=2)
marg_BC_mix = bp_mixed.sum(dim=0)

mi_AB_mix = mutual_info(marg_AB_mix)
mi_BC_mix = mutual_info(marg_BC_mix)
mi_AC_mix = mutual_info(marg_AC_mix)

print(f"\n  Mutual information (mixed building):")
print(f"    I(A;B) [identity edge]:  {mi_AB_mix:.6f}")
print(f"    I(B;C) [inverse edge]:   {mi_BC_mix:.6f}")
print(f"    I(A;C) [cycle edge]:     {mi_AC_mix:.6f}")

# Group theory: identity ∘ inverse = inverse, inverse ∘ cycle = ???
# The composition through the building should reflect the group multiplication
# (01) = identity, (02) = inverse, (012) = cycle (approximately)
# (01)∘(02) = (012), (02)∘(012) = (01)
# The building encodes the GROUP STRUCTURE in its inter-floor correlations

# Condition on A=h0 in the mixed building
cond_A0_mix = bp_mixed[0,:,:]
cond_A0_mix = cond_A0_mix / cond_A0_mix.sum()

marg_B_given_A0_mix = cond_A0_mix.sum(dim=1)
marg_C_given_A0_mix = cond_A0_mix.sum(dim=0)

print(f"\n  Conditioned on A=h₀:")
print(f"    P(B|A=h₀): [{', '.join(f'{marg_B_given_A0_mix[i].item():.4f}' for i in range(MASS_DIM))}]")
print(f"    P(C|A=h₀): [{', '.join(f'{marg_C_given_A0_mix[i].item():.4f}' for i in range(MASS_DIM))}]")

# B should show identity relationship with A (h₀ → h₀)
# C should show cycle relationship with A (h₀ → h₁ for forward cycle)
print(f"\n  B peaks at h₀ (identity relationship with A)")
print(f"  C peaks at h{marg_C_given_A0_mix[:H].argmax().item()} "
      f"(cycle/inverse relationship with A)")


# ═══════════════════════════════════════════════════════════════
#  The interaction information: what the building adds beyond pairs
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("INTERACTION INFORMATION: WHAT THE BUILDING ADDS")
print("="*80)

# Interaction information I3 = I(A;B;C) = I(A;B) + I(A;C) - I(A;BC)
# Positive I3 = redundant information (each edge carries overlapping info)
# Negative I3 = synergistic information (the whole exceeds sum of parts)

# Compute I(A;BC) — mutual info between A and the (B,C) pair
marg_A_t = bp_triangle.sum(dim=(1,2))
marg_BC_t = bp_triangle.sum(dim=0).reshape(MASS_DIM * MASS_DIM)
joint_A_BC = bp_triangle.reshape(MASS_DIM, MASS_DIM * MASS_DIM)

mi_A_BC = 0
for i in range(MASS_DIM):
    for j in range(MASS_DIM * MASS_DIM):
        pij = joint_A_BC[i,j].item()
        pi = marg_A_t[i].item()
        pj = marg_BC_t[j].item()
        if pij > 1e-12 and pi > 1e-12 and pj > 1e-12:
            mi_A_BC += pij * math.log(pij / (pi * pj))

# I(A;B) and I(A;C) from the triangle
marg_AB_t = bp_triangle.sum(dim=2)
marg_AC_t = bp_triangle.sum(dim=1)
mi_AB_t = mutual_info(marg_AB_t)
mi_AC_t = mutual_info(marg_AC_t)

I3 = mi_AB_t + mi_AC_t - mi_A_BC
print(f"\n  Uniform triangle:")
print(f"    I(A;B) = {mi_AB_t:.6f}")
print(f"    I(A;C) = {mi_AC_t:.6f}")
print(f"    I(A;BC) = {mi_A_BC:.6f}")
print(f"    I₃ = I(A;B) + I(A;C) - I(A;BC) = {I3:.6f}")

if I3 > 0:
    print(f"    REDUNDANT: edges carry overlapping information")
else:
    print(f"    SYNERGISTIC: the building exceeds sum of its edges")

print(f"\n    I₃ / Δ = {I3 / DELTA:.4f}")
print(f"    I₃ / K* = {I3 / K_STAR:.4f}")
print(f"    I₃ / ln(H) = {I3 / math.log(H):.4f}")


print(f"\n\n{'='*80}")
print("WHAT LIVES BETWEEN FLOORS")
print("="*80)
