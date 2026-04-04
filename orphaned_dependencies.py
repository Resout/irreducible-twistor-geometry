"""
Orphaned dependencies: entanglement through absence.

V-shape: A—B and A—C, no B—C crystal.
B and C share no direct portal. But they share A.
Their correlation is orphaned — unmediated, arising from
the building's shape, not from any specific edge.

Questions:
1. What is I(B;C) with no B-C crystal? (the orphaned dependency)
2. Does it exceed classical bounds? (data processing: I(B;C) ≤ min(I(A;B), I(A;C)))
3. What sector carries it? (trivial or standard?)
4. How does it compare to a REAL B-C crystal?
5. What happens with DIFFERENT crystal types? (A-B identity, A-C inverse)
   Does the orphaned B-C dependency encode the GROUP PRODUCT (01)∘(02) = (012)?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            born_probabilities, schmidt_number,
                            EPS_LOG, EPS_DIV)
from solver.crystals import Entangler, compose, kirkwood_product, classify_relationship

torch.set_grad_enabled(False)

n_seeds = 50

def make_crystal(corr):
    s = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        s += Entangler(corr, seed=seed).build().joint
    return s / n_seeds

def mutual_info(joint_2d):
    p = joint_2d.clone()
    if p.sum() < 1e-12: return 0
    p = p / p.sum()
    m0 = p.sum(dim=1)
    m1 = p.sum(dim=0)
    mi = 0
    for i in range(p.shape[0]):
        for j in range(p.shape[1]):
            pij = p[i,j].item()
            pi = m0[i].item()
            pj = m1[j].item()
            if pij > 1e-12 and pi > 1e-12 and pj > 1e-12:
                mi += pij * math.log(pij / (pi * pj))
    return mi


identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
inverse_corr = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)
cycle_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)


print("=" * 80)
print("ORPHANED DEPENDENCIES")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  V-shape: A—B, A—C, no B—C
# ═══════════════════════════════════════════════════════════════

print(f"\n--- V-shape: A connected to B and C, B-C orphaned ---\n")

c_AB = make_crystal(identity_corr)
c_AC = make_crystal(identity_corr)

# Build the V-shape joint P(A,B,C):
# Without a B-C crystal, we need B⊥C|A (conditional independence given A)
# P(A,B,C) = P(B|A) × P(C|A) × P(A)
# This is a product state in the B-C marginal conditioned on A

# Construct: for each value of A, B and C are independent
# P(A=a, B=b, C=c) ∝ P(A=a, B=b) × P(A=a, C=c) / P(A=a)

bp_AB = c_AB.abs().pow(2)
bp_AB = bp_AB / bp_AB.sum()
bp_AC = c_AC.abs().pow(2)
bp_AC = bp_AC / bp_AC.sum()

marg_A_from_AB = bp_AB.sum(dim=1)
marg_A_from_AC = bp_AC.sum(dim=1)

# Build the V-shape triple
v_triple = torch.zeros(MASS_DIM, MASS_DIM, MASS_DIM)
for a in range(MASS_DIM):
    if marg_A_from_AB[a] > 1e-12 and marg_A_from_AC[a] > 1e-12:
        p_b_given_a = bp_AB[a, :] / marg_A_from_AB[a]
        p_c_given_a = bp_AC[a, :] / marg_A_from_AC[a]
        p_a = (marg_A_from_AB[a] + marg_A_from_AC[a]) / 2  # average
        for b in range(MASS_DIM):
            for c in range(MASS_DIM):
                v_triple[a, b, c] = p_a * p_b_given_a[b] * p_c_given_a[c]

v_triple = v_triple / v_triple.sum()

# Extract B-C marginal (the orphaned dependency)
marg_BC_orphan = v_triple.sum(dim=0)  # sum over A

# For comparison: the FULL triangle building
triple_full = kirkwood_product(c_AB, make_crystal(identity_corr), c_AC)
bp_full = triple_full.abs().pow(2)
bp_full = bp_full / bp_full.sum()
marg_BC_full = bp_full.sum(dim=0)

# And the direct B-C crystal (if it existed)
c_BC_direct = make_crystal(identity_corr)
bp_BC_direct = c_BC_direct.abs().pow(2)
bp_BC_direct = bp_BC_direct / bp_BC_direct.sum()

mi_orphan = mutual_info(marg_BC_orphan)
mi_full = mutual_info(marg_BC_full)
mi_direct = mutual_info(bp_BC_direct)

# Classical bound: I(B;C) ≤ min(I(A;B), I(A;C))
mi_AB = mutual_info(bp_AB)
mi_AC = mutual_info(bp_AC)
classical_bound = min(mi_AB, mi_AC)

print(f"  I(A;B) = {mi_AB:.6f} (direct crystal)")
print(f"  I(A;C) = {mi_AC:.6f} (direct crystal)")
print(f"  Classical bound min(I(A;B), I(A;C)) = {classical_bound:.6f}")
print(f"\n  I(B;C) orphaned (V-shape):    {mi_orphan:.6f}")
print(f"  I(B;C) full building:          {mi_full:.6f}")
print(f"  I(B;C) direct crystal:         {mi_direct:.6f}")
print(f"\n  Orphan / classical bound: {mi_orphan/classical_bound:.4f}")
print(f"  Orphan / direct: {mi_orphan/mi_direct:.4f}")

# Schmidt of the orphaned B-C relationship
sn_orphan = schmidt_number(marg_BC_orphan.to(torch.cfloat))
sn_direct = schmidt_number(c_BC_direct)
print(f"\n  Schmidt of orphaned B-C: {sn_orphan:.4f}")
print(f"  Schmidt of direct B-C:   {sn_direct:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Mixed V-shape: A—B identity, A—C inverse
#  Does the orphaned B-C encode the group product?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("MIXED V-SHAPE: DOES THE ORPHAN ENCODE THE GROUP PRODUCT?")
print("="*80)

c_AB_id = make_crystal(identity_corr)    # A→B: identity (h_i → h_i)
c_AC_inv = make_crystal(inverse_corr)    # A→C: inverse (h_i → h_{2-i})

bp_AB_id = c_AB_id.abs().pow(2)
bp_AB_id = bp_AB_id / bp_AB_id.sum()
bp_AC_inv = c_AC_inv.abs().pow(2)
bp_AC_inv = bp_AC_inv / bp_AC_inv.sum()

marg_A_id = bp_AB_id.sum(dim=1)
marg_A_inv = bp_AC_inv.sum(dim=1)

# Build V-shape
mixed_v = torch.zeros(MASS_DIM, MASS_DIM, MASS_DIM)
for a in range(MASS_DIM):
    if marg_A_id[a] > 1e-12 and marg_A_inv[a] > 1e-12:
        p_b = bp_AB_id[a, :] / marg_A_id[a]
        p_c = bp_AC_inv[a, :] / marg_A_inv[a]
        p_a = (marg_A_id[a] + marg_A_inv[a]) / 2
        for b in range(MASS_DIM):
            for c in range(MASS_DIM):
                mixed_v[a, b, c] = p_a * p_b[b] * p_c[c]

mixed_v = mixed_v / mixed_v.sum()

# The orphaned B-C relationship
marg_BC_mixed = mixed_v.sum(dim=0)

print(f"\n  A→B: identity (h₀→h₀, h₁→h₁, h₂→h₂)")
print(f"  A→C: inverse  (h₀→h₂, h₁→h₁, h₂→h₀)")
print(f"  Orphaned B→C should be: inverse (identity ∘ inverse = inverse)")

# What does the orphaned B-C look like?
print(f"\n  Orphaned B-C Born matrix:")
for i in range(MASS_DIM):
    labels = ['h₀', 'h₁', 'h₂', 'θ ']
    print(f"    {labels[i]}: [{', '.join(f'{marg_BC_mixed[i,j].item():.4f}' for j in range(MASS_DIM))}]")

# Classify the orphaned relationship
rel_type, rel_conf = classify_relationship(marg_BC_mixed.to(torch.cfloat))
print(f"\n  Classified as: {rel_type} (confidence {rel_conf:.4f})")

# Compare to actual inverse crystal
print(f"\n  Actual inverse crystal Born matrix:")
bp_inv = c_AC_inv.abs().pow(2)
bp_inv = bp_inv / bp_inv.sum()
for i in range(MASS_DIM):
    labels = ['h₀', 'h₁', 'h₂', 'θ ']
    print(f"    {labels[i]}: [{', '.join(f'{bp_inv[i,j].item():.4f}' for j in range(MASS_DIM))}]")


# ═══════════════════════════════════════════════════════════════
#  More group products via orphaned dependencies
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("GROUP PRODUCTS FROM ORPHANED DEPENDENCIES")
print("="*80)

combos = [
    ("identity", "identity", identity_corr, identity_corr, "identity"),
    ("identity", "inverse", identity_corr, inverse_corr, "inverse"),
    ("identity", "cycle", identity_corr, cycle_corr, "cycle_fwd"),
    ("inverse", "inverse", inverse_corr, inverse_corr, "identity"),
    ("inverse", "cycle", inverse_corr, cycle_corr, "???"),
    ("cycle", "cycle", cycle_corr, cycle_corr, "inverse"),
]

print(f"\n  {'A→B':>12s} {'A→C':>12s} {'orphan B→C':>12s} {'expected':>12s} {'match':>6s}")
print("-" * 60)

for name_ab, name_ac, corr_ab, corr_ac, expected in combos:
    c1 = make_crystal(corr_ab)
    c2 = make_crystal(corr_ac)

    bp1 = c1.abs().pow(2); bp1 = bp1 / bp1.sum()
    bp2 = c2.abs().pow(2); bp2 = bp2 / bp2.sum()
    m1 = bp1.sum(dim=1)
    m2 = bp2.sum(dim=1)

    v = torch.zeros(MASS_DIM, MASS_DIM, MASS_DIM)
    for a in range(MASS_DIM):
        if m1[a] > 1e-12 and m2[a] > 1e-12:
            pb = bp1[a,:] / m1[a]
            pc = bp2[a,:] / m2[a]
            pa = (m1[a] + m2[a]) / 2
            for b in range(MASS_DIM):
                for c in range(MASS_DIM):
                    v[a,b,c] = pa * pb[b] * pc[c]
    v = v / v.sum()

    marg_bc = v.sum(dim=0)
    rel, conf = classify_relationship(marg_bc.to(torch.cfloat))
    match = "✓" if rel == expected or expected == "???" else "✗"
    print(f"  {name_ab:>12s} {name_ac:>12s} {rel:>12s} {expected:>12s} {match:>6s}")


# ═══════════════════════════════════════════════════════════════
#  Star topology: central hub, N spokes, all orphaned pairs
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("STAR TOPOLOGY: ONE HUB, MANY ORPHANS")
print("="*80)

# Hub A connected to B, C, D, E — but B,C,D,E share NO edges
# ALL pairwise dependencies among B,C,D,E are orphaned through A

print(f"\n  How does orphan strength scale with number of spokes?")

for n_spokes in [2, 3, 4, 5, 6]:
    crystals_spoke = [make_crystal(identity_corr) for _ in range(n_spokes)]

    # Build joint: P(A, B₁, ..., B_n) = P(A) × ∏ P(B_i|A)
    # Orphaned: I(B₁; B₂) — two random spokes

    bp_spokes = []
    margs_A = []
    for c_s in crystals_spoke:
        bp = c_s.abs().pow(2)
        bp = bp / bp.sum()
        bp_spokes.append(bp)
        margs_A.append(bp.sum(dim=1))

    # Compute orphaned B₁-B₂ marginal
    marg_12 = torch.zeros(MASS_DIM, MASS_DIM)
    p_A = margs_A[0]  # use first spoke's marginal for A
    for a in range(MASS_DIM):
        if p_A[a] > 1e-12:
            p_b1 = bp_spokes[0][a,:] / margs_A[0][a]
            p_b2 = bp_spokes[1][a,:] / margs_A[1][a]
            for b1 in range(MASS_DIM):
                for b2 in range(MASS_DIM):
                    marg_12[b1, b2] += p_A[a] * p_b1[b1] * p_b2[b2]

    mi_orphan = mutual_info(marg_12)

    # Also: total orphaned MI = sum over all spoke pairs
    total_orphan_mi = mi_orphan * n_spokes * (n_spokes - 1) / 2

    print(f"  {n_spokes} spokes: I(B₁;B₂) = {mi_orphan:.6f}, "
          f"total orphan MI = {total_orphan_mi:.4f} ({n_spokes*(n_spokes-1)//2} pairs)")


print(f"\n\n{'='*80}")
print("WHAT THE ORPHANS SHOW")
print("="*80)
