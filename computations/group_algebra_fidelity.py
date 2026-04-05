"""
The 8% fidelity deficit in the S₃ multiplication table.

Crystal composition gives (01)∘(01) ≈ e with fidelity 0.92, not 1.0.
Why? Three candidates:

1. Seed noise: the entangler is stochastic. Seed-averaging might fix it.
2. Structural filter: composition decays Schmidt by 2Δ per step.
   (01)∘(01) has Schmidt < identity, losing correlation.
3. The Born floor: it enforces ignorance (θ ≥ 1/27), which mixes
   every crystal toward the center. Two compositions compound this.

If the deficit is structural (not seed noise), it tells us something
about how the continuous group algebra sits inside the discrete crystal.
The structural filter means the group algebra is CONTRACTIVE — repeated
multiplication converges to the identity, not away from it.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities, born_fidelity)
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

# Group multiplication table for S₃
S3_MULT = {
    ("e","e"): "e", ("e","(01)"): "(01)", ("e","(02)"): "(02)",
    ("e","(12)"): "(12)", ("e","(012)"): "(012)", ("e","(021)"): "(021)",
    ("(01)","e"): "(01)", ("(01)","(01)"): "e", ("(01)","(02)"): "(012)",
    ("(01)","(12)"): "(021)", ("(01)","(012)"): "(02)", ("(01)","(021)"): "(12)",
    ("(02)","e"): "(02)", ("(02)","(01)"): "(021)", ("(02)","(02)"): "e",
    ("(02)","(12)"): "(012)", ("(02)","(012)"): "(12)", ("(02)","(021)"): "(01)",
    ("(12)","e"): "(12)", ("(12)","(01)"): "(012)", ("(12)","(02)"): "(021)",
    ("(12)","(12)"): "e", ("(12)","(012)"): "(01)", ("(12)","(021)"): "(02)",
    ("(012)","e"): "(012)", ("(012)","(01)"): "(12)", ("(012)","(02)"): "(01)",
    ("(012)","(12)"): "(02)", ("(012)","(012)"): "(021)", ("(012)","(021)"): "e",
    ("(021)","e"): "(021)", ("(021)","(01)"): "(02)", ("(021)","(02)"): "(12)",
    ("(021)","(12)"): "(01)", ("(021)","(012)"): "e", ("(021)","(021)"): "(012)",
}


def make_crystal(corr, seed=42):
    return Entangler(corr, seed=seed).build().joint


print("=" * 80)
print("THE FIDELITY DEFICIT IN THE S₃ GROUP ALGEBRA")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Test 1: Fidelity vs number of seeds
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Test 1: Does seed-averaging improve fidelity? ---\n")

for n_seeds in [1, 5, 10, 20, 50, 100]:
    # Build seed-averaged crystals
    refs = {}
    for name, corr in S3_CORRS.items():
        s = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
        for seed in range(n_seeds):
            s += make_crystal(corr, seed=seed)
        refs[name] = s / n_seeds

    # Test: (01)∘(01) should = e
    comp = compose(refs["(01)"], refs["(01)"])
    fid = born_fidelity(comp, refs["e"])
    sn_comp = schmidt_number(comp)
    sn_e = schmidt_number(refs["e"])

    # Also test a cross product: (01)∘(02) should = (012)
    comp2 = compose(refs["(01)"], refs["(02)"])
    fid2 = born_fidelity(comp2, refs["(012)"])

    print(f"  {n_seeds:3d} seeds: (01)²→e fid={fid:.4f}  (01)∘(02)→(012) fid={fid2:.4f}"
          f"  Schmidt(comp)={sn_comp:.3f} vs Schmidt(e)={sn_e:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 2: Schmidt decay through composition
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Test 2: Schmidt decay from composition (structural filter) ---\n")

n_seeds = 50
refs = {}
for name, corr in S3_CORRS.items():
    s = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        s += make_crystal(corr, seed=seed)
    refs[name] = s / n_seeds

for g1 in ["(01)", "(02)", "(12)", "(012)", "(021)"]:
    for g2 in ["(01)", "(02)", "(12)", "(012)", "(021)"]:
        expected = S3_MULT[(g1, g2)]
        comp = compose(refs[g1], refs[g2])
        sn_comp = schmidt_number(comp)
        sn_expected = schmidt_number(refs[expected])
        ratio = sn_comp / sn_expected
        fid = born_fidelity(comp, refs[expected])
        if g2 == "(01)":  # print one column
            print(f"  {g1:>6s}∘{g2:>6s} = {expected:>6s}: "
                  f"Schmidt {sn_comp:.3f}/{sn_expected:.3f} = {ratio:.3f}, fid={fid:.4f}")

# Average ratio
all_ratios = []
all_fids = []
for g1 in S3_CORRS:
    for g2 in S3_CORRS:
        expected = S3_MULT[(g1, g2)]
        comp = compose(refs[g1], refs[g2])
        all_ratios.append(schmidt_number(comp) / schmidt_number(refs[expected]))
        all_fids.append(born_fidelity(comp, refs[expected]))

print(f"\n  Across all 36 products:")
print(f"    Mean Schmidt ratio: {sum(all_ratios)/len(all_ratios):.4f}")
print(f"    Mean fidelity:      {sum(all_fids)/len(all_fids):.4f}")
print(f"    Min fidelity:       {min(all_fids):.4f}")
print(f"    Max fidelity:       {max(all_fids):.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 3: Iterated composition — does it converge?
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Test 3: Iterated composition — does the group algebra contract? ---\n")

# Start with (01), compose with itself repeatedly
# (01)^1 = (01), (01)^2 = e, (01)^3 = (01), (01)^4 = e, ...
crystal = refs["(01)"]
for step in range(1, 13):
    crystal = compose(crystal, refs["(01)"])
    expected = "e" if step % 2 == 1 else "(01)"
    fid = born_fidelity(crystal, refs[expected])
    sn = schmidt_number(crystal)
    print(f"  (01)^{step+1:2d}: Schmidt={sn:.4f}, fid({expected})={fid:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 4: What does the composition converge TO?
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Test 4: What is the fixed point of iterated composition? ---\n")

# After many compositions, what does the crystal look like?
crystal = refs["(01)"]
for _ in range(20):
    crystal = compose(crystal, refs["(01)"])

bp = born_probabilities(crystal).reshape(MASS_DIM, MASS_DIM)
print(f"  Born matrix of (01)^21:")
for i in range(MASS_DIM):
    print(f"    [{', '.join(f'{bp[i,j].item():.4f}' for j in range(MASS_DIM))}]")

print(f"\n  Schmidt: {schmidt_number(crystal):.4f}")
print(f"  Is this the product state? (q = p²)")

# Check: is each row proportional to the same vector?
rows = [bp[i] / bp[i].sum() for i in range(MASS_DIM)]
for i in range(MASS_DIM):
    print(f"    Row {i} (normalized): [{', '.join(f'{rows[i][j].item():.4f}' for j in range(MASS_DIM))}]")


# ═══════════════════════════════════════════════════════════════
#  Test 5: The deficit is 2Δ per composition — verify
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Test 5: Is the fidelity deficit exactly exp(-2Δ)? ---\n")

# If composition applies the structural filter at rate 2Δ,
# then one composition should reduce fidelity by exp(-2Δ) ≈ 0.588
# But fidelity is 0.92, not 0.59. So the deficit is NOT 2Δ.
# The 2Δ applies to Schmidt ratio, not to Born fidelity.

exp_2delta = math.exp(-2 * DELTA)
print(f"  exp(-2Δ) = {exp_2delta:.4f}")
print(f"  Mean fidelity = {sum(all_fids)/len(all_fids):.4f}")
print(f"  Mean Schmidt ratio = {sum(all_ratios)/len(all_ratios):.4f}")
print(f"  exp(-2Δ) * initial_Schmidt/composed_Schmidt?")

# Check if 1 - fidelity scales with something
mean_fid = sum(all_fids)/len(all_fids)
deficit = 1 - mean_fid
print(f"\n  Deficit = 1 - fidelity = {deficit:.4f}")
print(f"  Deficit / Δ = {deficit / DELTA:.4f}")
print(f"  Deficit / K* = {deficit / K_STAR:.4f}")
print(f"  Deficit / Born(θ) = {deficit / (1/27):.4f}")

# Is the deficit the Born floor contribution?
# Every crystal has Born(θ) ≥ 1/27. The "wrong" part of the
# composed crystal could be the θ×θ cross-term.
print(f"\n  Born floor = {1/27:.6f}")
print(f"  2 × Born floor = {2/27:.6f}")
print(f"  Deficit = {deficit:.6f}")


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
