"""
Wilson loops in the crystal lattice.

W(C) = tr(compose(c₁, c₂, ..., cₙ)) around a closed path C.

The trace measures the return amplitude: how much of h₀ comes
back as h₀ after traversing the loop. It's gauge-invariant.

Confinement: |W| decays with AREA of the loop (area law)
Freedom: |W| decays with PERIMETER (perimeter law)

The structural filter (299/405 per step) gives perimeter decay.
But the BUILDING might show something different — if the full
joint state preserves what sequential composition loses, the
Wilson loop computed from the building might not decay at all.

Two Wilson loops:
1. Sequential: compose crystals around the path, take trace
2. Building: construct full joint, trace the diagonal
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities,
                            EPS_LOG, EPS_DIV)
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


def wilson_loop_sequential(crystals):
    """Compose crystals around a closed path, return tr(result)."""
    result = crystals[0]
    for c in crystals[1:]:
        result = compose(result, c)
    # Trace: sum of diagonal
    tr = sum(result[i, i] for i in range(MASS_DIM))
    return tr, result


def wilson_loop_value(crystals):
    """Just the magnitude of the Wilson loop."""
    tr, _ = wilson_loop_sequential(crystals)
    return tr.abs().item()


print("=" * 80)
print("WILSON LOOPS IN THE CRYSTAL LATTICE")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Simple loops: triangle, square, pentagon
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Wilson loops of different sizes (all identity bonds) ---\n")

bond = make_crystal(identity_corr)

print(f"  {'perimeter':>9s} {'|W|':>10s} {'ln|W|':>10s} {'-ln|W|/P':>10s}")
print("-" * 45)

for n_sides in range(3, 10):
    crystals = [bond] * n_sides
    w = wilson_loop_value(crystals)
    ln_w = math.log(max(w, 1e-15))
    rate = -ln_w / n_sides

    print(f"  {n_sides:9d} {w:>10.6f} {ln_w:>10.4f} {rate:>10.5f}")

print(f"\n  If perimeter law: -ln|W|/P = constant")
print(f"  Δ = {DELTA:.5f}")
print(f"  2Δ = {2*DELTA:.5f}")


# ═══════════════════════════════════════════════════════════════
#  Rectangle loops: area vs perimeter
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Rectangular loops: R×T (area vs perimeter) ---\n")

print(f"  {'R×T':>5s} {'area':>5s} {'perim':>5s} {'|W|':>10s} {'ln|W|':>10s} "
      f"{'-ln/area':>10s} {'-ln/perim':>10s}")
print("-" * 65)

for R in range(1, 5):
    for T in range(R, 5):
        # Rectangle: R bonds right, T bonds up, R bonds left, T bonds down
        perimeter = 2 * (R + T)
        area = R * T
        crystals = [bond] * perimeter
        w = wilson_loop_value(crystals)
        ln_w = math.log(max(w, 1e-15))

        print(f"  {R}×{T:>2d} {area:>5d} {perimeter:>5d} {w:>10.6f} {ln_w:>10.4f} "
              f"{-ln_w/area:>10.5f} {-ln_w/perimeter:>10.5f}")


# ═══════════════════════════════════════════════════════════════
#  Wilson loop from the BUILDING (Kirkwood) vs sequential
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("WILSON LOOP: BUILDING vs SEQUENTIAL")
print("="*80)

# Triangle: sequential W = tr(compose(AB, BC, CA))
c_AB = make_crystal(identity_corr)
c_BC = make_crystal(identity_corr)
c_CA = make_crystal(identity_corr)

w_seq, comp_seq = wilson_loop_sequential([c_AB, c_BC, c_CA])
print(f"\n  Triangle (sequential):")
print(f"    W = tr(C_AB ∘ C_BC ∘ C_CA) = {w_seq.abs().item():.6f}")

# Triangle: building W = sum_a P(A=a, B=a, C=a) ???
# The Wilson loop in the building is: the probability that all
# sites agree on the same hypothesis
triple = kirkwood_product(c_AB, c_BC, c_CA)
bp = triple.abs().pow(2)
bp = bp / bp.sum()

# Diagonal: P(A=a, B=a, C=a) for each a
w_building = 0
for a in range(MASS_DIM):
    w_building += bp[a, a, a].item()

print(f"\n  Triangle (building):")
print(f"    W = Σ_a P(A=a, B=a, C=a) = {w_building:.6f}")

# Also: the "return amplitude" — condition on A=a, what's P(it returns to a)?
print(f"\n  Return amplitude P(C=a|A=a):")
for a in range(H):
    # Condition on A=a
    slice_a = bp[a, :, :]
    if slice_a.sum() > 1e-12:
        slice_a = slice_a / slice_a.sum()
        # Trace over B: P(C=a|A=a) = Σ_b P(B=b, C=a | A=a)
        p_return = sum(slice_a[b, a].item() for b in range(MASS_DIM))
        print(f"    P(C=h{a}|A=h{a}) = {p_return:.5f}")

# Compare: sequential P(C=a|A=a)
bp_seq = born_probabilities(comp_seq)
for a in range(H):
    print(f"    Sequential: (C_AB∘C_BC∘C_CA)[h{a},h{a}] Born = {bp_seq.reshape(MASS_DIM,MASS_DIM)[a,a].item():.5f}")


# ═══════════════════════════════════════════════════════════════
#  Mixed Wilson loops: non-trivial flux
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("MIXED WILSON LOOPS — NON-TRIVIAL FLUX")
print("="*80)

# Triangle with identity, inverse, cycle bonds
# Group product: cycle ∘ inverse ∘ identity = ???
# (012) ∘ (02) ∘ e = (012) ∘ (02) = (01)
# So the flux is the transposition (01)

c_id = make_crystal(identity_corr)
c_inv = make_crystal(inverse_corr)
c_cyc = make_crystal(cycle_corr)

print(f"\n  Identity → Inverse → Cycle triangle:")
w_mixed, comp_mixed = wilson_loop_sequential([c_id, c_inv, c_cyc])
print(f"    |W| = {w_mixed.abs().item():.6f}")
print(f"    Flux = group product = cycle ∘ inverse ∘ identity")

# The diagonal of the composition tells the flux type
bp_mixed = born_probabilities(comp_mixed).reshape(MASS_DIM, MASS_DIM)
print(f"    Born diagonal: [{', '.join(f'{bp_mixed[i,i].item():.4f}' for i in range(MASS_DIM))}]")

# Pure identity triangle for comparison
w_pure, _ = wilson_loop_sequential([c_id, c_id, c_id])
print(f"\n  Identity × 3:")
print(f"    |W| = {w_pure.abs().item():.6f}")

# Pure inverse triangle: inverse³ = inverse (since inverse² = identity)
w_inv3, _ = wilson_loop_sequential([c_inv, c_inv, c_inv])
print(f"\n  Inverse × 3:")
print(f"    |W| = {w_inv3.abs().item():.6f}")

# Pure cycle triangle: cycle³ = identity
w_cyc3, _ = wilson_loop_sequential([c_cyc, c_cyc, c_cyc])
print(f"\n  Cycle × 3:")
print(f"    |W| = {w_cyc3.abs().item():.6f}")


# ═══════════════════════════════════════════════════════════════
#  The orphan IS the virtual Wilson loop
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE ORPHAN IS THE VIRTUAL WILSON LOOP")
print("="*80)

# V-shape A→B, A→C (identity bonds)
# The orphaned B-C dependency: what happens if we "close" the loop
# through the orphan?

# The Wilson loop A→B→(orphan)→C→A is:
# compose(c_AB, orphan_BC, c_CA)
# where orphan_BC is the "virtual crystal" from the V-shape

# Build the orphan crystal: it's the B-C marginal from the V-shape
c_AB2 = make_crystal(identity_corr)
c_AC2 = make_crystal(identity_corr)

bp_AB2 = c_AB2.abs().pow(2); bp_AB2 = bp_AB2 / bp_AB2.sum()
bp_AC2 = c_AC2.abs().pow(2); bp_AC2 = bp_AC2 / bp_AC2.sum()
mA_ab = bp_AB2.sum(dim=1)
mA_ac = bp_AC2.sum(dim=1)

orphan_BC = torch.zeros(MASS_DIM, MASS_DIM)
for a in range(MASS_DIM):
    if mA_ab[a] > 1e-12 and mA_ac[a] > 1e-12:
        pb = bp_AB2[a,:] / mA_ab[a]
        pc = bp_AC2[a,:] / mA_ac[a]
        for b in range(MASS_DIM):
            for c in range(MASS_DIM):
                orphan_BC[b,c] += (mA_ab[a] + mA_ac[a])/2 * pb[b] * pc[c]
orphan_BC = orphan_BC / orphan_BC.sum()

# The "virtual Wilson loop" through the orphan
# But first: what is the trace of the orphan itself?
tr_orphan = sum(orphan_BC[i,i].item() for i in range(MASS_DIM))
print(f"\n  Trace of orphan B-C: {tr_orphan:.6f}")
print(f"  Trace of direct B-C: {sum(bp_BC_direct[i,i].item() for i in range(MASS_DIM)):.6f}"
      if 'bp_BC_direct' in dir() else "")

# The Wilson loop through real edges
w_real = wilson_loop_value([c_AB2, make_crystal(identity_corr), c_AC2.T])
# Note: c_AC goes from A to C, but we need C to A, so transpose
print(f"\n  W(real triangle): |W| = {w_real:.6f}")

# The orphan trace tells us: how much of the loop's flux is
# captured by the orphan alone?
print(f"\n  The orphan's trace ({tr_orphan:.4f}) captures the loop's")
print(f"  diagonal structure without any explicit B-C crystal.")
print(f"  The virtual Wilson loop IS the orphaned dependency.")


print(f"\n\n{'='*80}")
print("WHAT THE LOOPS SHOW")
print("="*80)
