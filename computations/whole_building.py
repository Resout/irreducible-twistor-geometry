"""
The whole building exists at once.

Don't propagate. Build the full joint state. Then look at it
from different floors.

A triangle of 3 sites with crystals on each edge.
The Kirkwood product gives the full [4,4,4] joint state.
Marginalizing to any single site should automatically contain
the multi-path information — without sequential propagation.

The sequential propagation (Principle 50) was looking at the map
one floor at a time. This computation looks at the whole building.
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


print("=" * 80)
print("THE WHOLE BUILDING — FULL JOINT STATE")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Triangle: 3 sites, 3 edges
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Triangle: sites A, B, C ---")

c_AB = make_crystal(identity_corr)
c_BC = make_crystal(identity_corr)
c_AC = make_crystal(identity_corr)

# The Kirkwood product: full [4,4,4] joint from 3 pairwise [4,4] joints
triple = kirkwood_product(c_AB, c_BC, c_AC)

print(f"\n  Triple joint shape: {triple.shape}")
print(f"  Triple L1: {triple.reshape(-1).real.sum().item():.6f}")

# Born probabilities of the full triple
bp_triple = triple.abs().pow(2)
bp_triple = bp_triple / bp_triple.sum()

# Marginals: what each site looks like when the building is complete
marg_A = bp_triple.sum(dim=(1,2))  # sum over B and C
marg_B = bp_triple.sum(dim=(0,2))  # sum over A and C
marg_C = bp_triple.sum(dim=(0,1))  # sum over A and B

print(f"\n  Marginal at A: [{', '.join(f'{marg_A[i].item():.4f}' for i in range(MASS_DIM))}]")
print(f"  Marginal at B: [{', '.join(f'{marg_B[i].item():.4f}' for i in range(MASS_DIM))}]")
print(f"  Marginal at C: [{', '.join(f'{marg_C[i].item():.4f}' for i in range(MASS_DIM))}]")

# These should be identical by symmetry (all crystals are the same)
print(f"\n  A-B difference: {(marg_A - marg_B).abs().max().item():.6f}")
print(f"  A-C difference: {(marg_A - marg_C).abs().max().item():.6f}")


# ═══════════════════════════════════════════════════════════════
#  Condition on A = h0: what does B look like?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Conditioning: if A = h₀, what is B? ---")

# Condition on A = h0: slice the triple at A=0, then marginal over C
conditioned_A0 = bp_triple[0, :, :]  # [4,4]: B × C given A=h0
conditioned_A0 = conditioned_A0 / conditioned_A0.sum()
marg_B_given_A0 = conditioned_A0.sum(dim=1)  # sum over C

print(f"  P(B | A=h₀): [{', '.join(f'{marg_B_given_A0[i].item():.4f}' for i in range(MASS_DIM))}]")

# Compare to what sequential propagation gives
signal_A0 = torch.zeros(MASS_DIM, dtype=torch.cfloat)
signal_A0[0] = 0.85 + 0j; signal_A0[1] = 0.05 + 0j
signal_A0[2] = 0.05 + 0j; signal_A0[H] = 0.05 + 0j

weighted = c_AB * signal_A0.unsqueeze(1)
marginal_seq = weighted.sum(dim=0)
re_sum = marginal_seq.real.sum()
if abs(re_sum) > EPS_DIV:
    marginal_seq = marginal_seq / re_sum
bp_seq = born_probabilities(marginal_seq)

print(f"  Sequential A→B:   [{', '.join(f'{bp_seq[i].item():.4f}' for i in range(MASS_DIM))}]")

# The Kirkwood conditioning includes info from ALL paths (A→B and A→C→B)
# The sequential only uses A→B
print(f"\n  The building sees more: conditioning includes all paths automatically")


# ═══════════════════════════════════════════════════════════════
#  Condition on A = h0: what is C?
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Conditioning: if A = h₀, what is C? ---")

# Two perspectives on C given A:
# 1. Direct (A→C): one floor
# 2. From the building: all floors at once

# Building:
conditioned_A0_C = bp_triple[0, :, :]  # [4,4]: B × C given A=h0
conditioned_A0_C = conditioned_A0_C / conditioned_A0_C.sum()
marg_C_given_A0 = conditioned_A0_C.sum(dim=0)  # sum over B

print(f"  Building P(C|A=h₀): [{', '.join(f'{marg_C_given_A0[i].item():.4f}' for i in range(MASS_DIM))}]")

# Sequential via direct A→C
weighted = c_AC * signal_A0.unsqueeze(1)
marginal_direct = weighted.sum(dim=0)
re_sum = marginal_direct.real.sum()
if abs(re_sum) > EPS_DIV:
    marginal_direct = marginal_direct / re_sum
bp_direct = born_probabilities(marginal_direct)
print(f"  Direct A→C:         [{', '.join(f'{bp_direct[i].item():.4f}' for i in range(MASS_DIM))}]")

# Sequential via A→B→C (2 hops)
weighted1 = c_AB * signal_A0.unsqueeze(1)
m1 = weighted1.sum(dim=0)
re1 = m1.real.sum()
if abs(re1) > EPS_DIV: m1 = m1 / re1
m1 = enforce_born_floor(m1.unsqueeze(0)).squeeze(0)

weighted2 = c_BC * m1.unsqueeze(1)
m2 = weighted2.sum(dim=0)
re2 = m2.real.sum()
if abs(re2) > EPS_DIV: m2 = m2 / re2
bp_2hop = born_probabilities(m2)
print(f"  Sequential A→B→C:   [{', '.join(f'{bp_2hop[i].item():.4f}' for i in range(MASS_DIM))}]")

# Multi-path DS combine
combined_seq, K = ds_combine(marginal_direct.unsqueeze(0), m2.unsqueeze(0))
bp_multi_seq = born_probabilities(combined_seq.squeeze(0))
print(f"  Multi-path (DS):    [{', '.join(f'{bp_multi_seq[i].item():.4f}' for i in range(MASS_DIM))}]")

# The building should see more than any sequential approach
print(f"\n  Born(h₀) comparison:")
print(f"    Building:        {marg_C_given_A0[0].item():.5f}")
print(f"    Direct:          {bp_direct[0].item():.5f}")
print(f"    2-hop:           {bp_2hop[0].item():.5f}")
print(f"    Multi-path DS:   {bp_multi_seq[0].item():.5f}")


# ═══════════════════════════════════════════════════════════════
#  The Schmidt number of the triple
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Entanglement structure of the building ---")

from solver.crystals import schmidt_3way
sn3 = schmidt_3way(triple)
print(f"  3-way Schmidt (A vs BC): {sn3:.4f}")

# 2-way Schmidt of each edge
print(f"  2-way Schmidt AB: {schmidt_number(c_AB):.4f}")
print(f"  2-way Schmidt BC: {schmidt_number(c_BC):.4f}")
print(f"  2-way Schmidt AC: {schmidt_number(c_AC):.4f}")


# ═══════════════════════════════════════════════════════════════
#  Tetrahedron: 4 sites, 6 edges — the complete building
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("TETRAHEDRON: 4 SITES, 6 EDGES")
print("="*80)

# Build all 6 crystals
c = {}
pairs = [('A','B'), ('A','C'), ('A','D'), ('B','C'), ('B','D'), ('C','D')]
for p in pairs:
    c[p] = make_crystal(identity_corr)

# Build the full [4,4,4,4] joint via iterated Kirkwood
# triple_ABC = kirkwood(AB, BC, AC)
# Then extend to D via the three edges BD, CD, AD

triple_ABC = kirkwood_product(c[('A','B')], c[('B','C')], c[('A','C')])

# For each d-index, weight the triple by the three D-edges
quad = torch.zeros(MASS_DIM, MASS_DIM, MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for d in range(MASS_DIM):
    # Weight from A→D, B→D, C→D
    w_AD = c[('A','D')][:, d]  # [4]: A-component given D=d
    w_BD = c[('B','D')][:, d]  # [4]: B-component given D=d
    w_CD = c[('C','D')][:, d]  # [4]: C-component given D=d

    for a in range(MASS_DIM):
        for b in range(MASS_DIM):
            for cc_idx in range(MASS_DIM):
                quad[a, b, cc_idx, d] = triple_ABC[a, b, cc_idx] * w_AD[a] * w_BD[b] * w_CD[cc_idx]

# Normalize
l1 = quad.reshape(-1).real.sum()
if abs(l1) > EPS_DIV:
    quad = quad / l1

bp_quad = quad.abs().pow(2)
bp_quad = bp_quad / bp_quad.sum()

# Marginals
marg_A4 = bp_quad.sum(dim=(1,2,3))
marg_D4 = bp_quad.sum(dim=(0,1,2))

print(f"\n  Marginal at A: [{', '.join(f'{marg_A4[i].item():.4f}' for i in range(MASS_DIM))}]")
print(f"  Marginal at D: [{', '.join(f'{marg_D4[i].item():.4f}' for i in range(MASS_DIM))}]")

# Condition on A = h0
cond_A0_quad = bp_quad[0, :, :, :]
cond_A0_quad = cond_A0_quad / cond_A0_quad.sum()
marg_D_given_A0 = cond_A0_quad.sum(dim=(0,1))  # sum over B, C

print(f"\n  Building P(D|A=h₀): [{', '.join(f'{marg_D_given_A0[i].item():.4f}' for i in range(MASS_DIM))}]")

# Compare to sequential
# Direct A→D
weighted = c[('A','D')] * signal_A0.unsqueeze(1)
m_direct = weighted.sum(dim=0)
re = m_direct.real.sum()
if abs(re) > EPS_DIV: m_direct = m_direct / re
bp_d_direct = born_probabilities(m_direct)
print(f"  Direct A→D:         [{', '.join(f'{bp_d_direct[i].item():.4f}' for i in range(MASS_DIM))}]")

print(f"\n  Born(h₀) at D:")
print(f"    Building:  {marg_D_given_A0[0].item():.5f}")
print(f"    Direct:    {bp_d_direct[0].item():.5f}")
print(f"    Ratio:     {marg_D_given_A0[0].item() / bp_d_direct[0].item():.4f}")


# ═══════════════════════════════════════════════════════════════
#  The key question: does the building contain more than the sum?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("DOES THE BUILDING CONTAIN MORE THAN THE SUM OF ITS FLOORS?")
print("="*80)

# If the user's insight is right, the full joint state (the building)
# should contain information that NO sequential propagation can recover.
# The building IS the state. Sequential propagation only APPROXIMATES it.

# Test: mutual information between non-adjacent sites in the building
# vs what propagation gives

# In the triangle: conditional mutual information I(A;C|B)
# Building: compute from the triple directly
# Sequential: zero (Markov chain A→B→C has I(A;C|B) = 0)

# P(A,C) from building
marg_AC = bp_triple.sum(dim=1)  # sum over B
marg_A_only = bp_triple.sum(dim=(1,2))
marg_C_only = bp_triple.sum(dim=(0,1))

# Mutual information I(A;C)
mi_AC = 0
for a in range(MASS_DIM):
    for cc in range(MASS_DIM):
        p_ac = marg_AC[a, cc].item()
        p_a = marg_A_only[a].item()
        p_c = marg_C_only[cc].item()
        if p_ac > 1e-10 and p_a > 1e-10 and p_c > 1e-10:
            mi_AC += p_ac * math.log(p_ac / (p_a * p_c))

# For comparison: what is I(A;B)?
marg_AB = bp_triple.sum(dim=2)
marg_B_only = bp_triple.sum(dim=(0,2))
mi_AB = 0
for a in range(MASS_DIM):
    for b in range(MASS_DIM):
        p_ab = marg_AB[a, b].item()
        p_a = marg_A_only[a].item()
        p_b = marg_B_only[b].item()
        if p_ab > 1e-10 and p_a > 1e-10 and p_b > 1e-10:
            mi_AB += p_ab * math.log(p_ab / (p_a * p_b))

print(f"\n  Mutual information in the triangle building:")
print(f"    I(A;B) = {mi_AB:.6f} nats (direct neighbors)")
print(f"    I(A;C) = {mi_AC:.6f} nats (connected via B)")
print(f"    I(A;C)/I(A;B) = {mi_AC/mi_AB:.4f}")
print(f"\n  In a Markov chain A→B→C, I(A;C|B) = 0 by data processing.")
print(f"  In the building, A and C are connected BOTH directly AND via B.")
print(f"  The building has I(A;C) ≈ I(A;B): the unrepresented floors persist.")


print(f"\n\n{'='*80}")
print("WHAT THE BUILDING SHOWS")
print("="*80)
