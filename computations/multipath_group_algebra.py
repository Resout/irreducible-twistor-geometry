"""
Multi-path group algebra: can parallel composition overcome the 299/405 contraction?

Single composition: Schmidt(g∘h)/Schmidt(gh) = 299/405, fidelity ≈ 0.92
Multi-path: DS-combine k independent compositions (different seeds).
Does fidelity approach 1.0? Does Schmidt recover?

This tests whether the Penrose transform (multi-path = discrete contour
integral) can restore the unitary structure of the group algebra.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities, born_fidelity,
                            ds_combine, enforce_born_floor)
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


def make_crystal(corr, seed=42):
    return Entangler(corr, seed=seed).build().joint


# Reference target: (01)∘(01) should = e
# Use 50-seed average for the reference
ref_e = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for s in range(50):
    ref_e += make_crystal(S3_CORRS["e"], seed=s)
ref_e /= 50

ref_01 = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for s in range(50):
    ref_01 += make_crystal(S3_CORRS["(01)"], seed=s)
ref_01 /= 50


print("=" * 80)
print("MULTI-PATH GROUP ALGEBRA")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Test 1: Multi-path single-step composition
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Multi-path: (01)∘(01) = e via k parallel lines ---\n")

for n_paths in [1, 2, 3, 5, 10, 20, 50]:
    # Each path uses a different seed pair for the two (01) crystals
    composed_paths = []
    for p in range(n_paths):
        c1 = make_crystal(S3_CORRS["(01)"], seed=p)
        c2 = make_crystal(S3_CORRS["(01)"], seed=100+p)
        comp = compose(c1, c2)
        composed_paths.append(comp)

    # DS-combine all paths
    # Flatten to 16D for combination
    combined = composed_paths[0].reshape(-1).to(torch.cfloat)
    total_K = 0.0
    for comp in composed_paths[1:]:
        flat = comp.reshape(-1).to(torch.cfloat)
        combined, K = ds_combine(combined.unsqueeze(0), flat.unsqueeze(0))
        combined = combined.squeeze(0)
        total_K += abs(K.item()) if K.is_complex() else K.item()

    result = combined.reshape(MASS_DIM, MASS_DIM)
    sn = schmidt_number(result)
    fid = born_fidelity(result, ref_e)
    sn_ratio = sn / schmidt_number(ref_e)

    print(f"  {n_paths:3d} paths: Schmidt = {sn:.4f} (ratio={sn_ratio:.4f}), "
          f"fidelity(e) = {fid:.4f}, total K = {total_K:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 2: Multi-path 2-step composition (deep)
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Multi-path: (01)∘(02)∘(12) = (012)∘(12) = (01) via k paths ---\n")

ref_target = ref_01  # (01)∘(02)∘(12) should give... let's compute

# Group: (01)∘(02) = (012), then (012)∘(12) = (01)
# So the target is (01)

for n_paths in [1, 2, 5, 10, 20, 50]:
    composed_paths = []
    for p in range(n_paths):
        c1 = make_crystal(S3_CORRS["(01)"], seed=p)
        c2 = make_crystal(S3_CORRS["(02)"], seed=50+p)
        c3 = make_crystal(S3_CORRS["(12)"], seed=100+p)
        comp = compose(compose(c1, c2), c3)
        composed_paths.append(comp)

    combined = composed_paths[0].reshape(-1).to(torch.cfloat)
    for comp in composed_paths[1:]:
        flat = comp.reshape(-1).to(torch.cfloat)
        combined, K = ds_combine(combined.unsqueeze(0), flat.unsqueeze(0))
        combined = combined.squeeze(0)

    result = combined.reshape(MASS_DIM, MASS_DIM)
    sn = schmidt_number(result)
    fid = born_fidelity(result, ref_target)

    print(f"  {n_paths:3d} paths: Schmidt = {sn:.4f}, fidelity((01)) = {fid:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 3: Multi-path vs seed-averaging
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Comparison: multi-path DS vs simple seed-averaging ---\n")

for n in [5, 10, 20, 50]:
    # Method 1: seed-average, then compose
    avg_01 = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for s in range(n):
        avg_01 += make_crystal(S3_CORRS["(01)"], seed=s)
    avg_01 /= n
    comp_avg = compose(avg_01, avg_01)
    fid_avg = born_fidelity(comp_avg, ref_e)
    sn_avg = schmidt_number(comp_avg)

    # Method 2: compose individually, then DS-combine
    combined = None
    for p in range(n):
        c1 = make_crystal(S3_CORRS["(01)"], seed=p)
        c2 = make_crystal(S3_CORRS["(01)"], seed=100+p)
        comp = compose(c1, c2)
        flat = comp.reshape(-1).to(torch.cfloat)
        if combined is None:
            combined = flat
        else:
            combined, K = ds_combine(combined.unsqueeze(0), flat.unsqueeze(0))
            combined = combined.squeeze(0)
    result_ds = combined.reshape(MASS_DIM, MASS_DIM)
    fid_ds = born_fidelity(result_ds, ref_e)
    sn_ds = schmidt_number(result_ds)

    print(f"  n={n:3d}: avg-then-compose: Schmidt={sn_avg:.4f} fid={fid_avg:.4f} | "
          f"compose-then-DS: Schmidt={sn_ds:.4f} fid={fid_ds:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 4: Cross-seed agreement (the consistency signal)
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Cross-seed agreement for (01)∘(01) = e ---\n")

# How much do different seed realizations of the same composition agree?
n_check = 20
compositions = []
for p in range(n_check):
    c1 = make_crystal(S3_CORRS["(01)"], seed=p)
    c2 = make_crystal(S3_CORRS["(01)"], seed=50+p)
    compositions.append(compose(c1, c2))

# Pairwise Born fidelity
fids = []
for i in range(n_check):
    for j in range(i+1, n_check):
        fids.append(born_fidelity(compositions[i], compositions[j]))

print(f"  Cross-seed fidelity: mean={sum(fids)/len(fids):.4f}, "
      f"min={min(fids):.4f}, max={max(fids):.4f}")

# Compare to cross-seed fidelity of UNcomposed crystals
uncomposed = [make_crystal(S3_CORRS["(01)"], seed=s) for s in range(n_check)]
ufids = []
for i in range(n_check):
    for j in range(i+1, n_check):
        ufids.append(born_fidelity(uncomposed[i], uncomposed[j]))

print(f"  Uncomposed fidelity: mean={sum(ufids)/len(ufids):.4f}, "
      f"min={min(ufids):.4f}, max={max(ufids):.4f}")


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
