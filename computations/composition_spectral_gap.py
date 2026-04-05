"""
The composition operator's spectral gap.

The 16×16 composition operator T_B has:
  - 4-fold degenerate leading eigenvalue ≈ 0.975
  - Second cluster ≈ 0.57
  - Gap: 0.975/0.57 ≈ 1.71

Questions:
1. Why is the leading eigenvalue 4-fold degenerate? (= MASS_DIM)
2. Is the gap related to Δ? Gap ≈ -ln(0.57/0.975) ≈ 0.536 ≈ 2Δ
3. Does T_B have a rank-1 fixed point? What is it?
4. How does the spectrum change across S₃ elements?
5. Connection to the 4-dim eigenvalue spectrum (r₁ = 2.03-2.21 Δ)
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import numpy as np
from solver.algebra import H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, schmidt_number
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)


def build_composition_operator(B_joint):
    """16×16 linear map: A → compose(A, B)."""
    dim = MASS_DIM * MASS_DIM
    T = torch.zeros(dim, dim, dtype=torch.cfloat)
    for col in range(dim):
        A = torch.zeros(dim, dtype=torch.cfloat)
        A[col] = 1.0 + 0j
        A = A.reshape(MASS_DIM, MASS_DIM)
        result = compose(A, B_joint)
        T[:, col] = result.reshape(dim)
    return T


# ═══════════════════════════════════════════════════════════════
#  The spectral gap of T_B for each crystal type
# ═══════════════════════════════════════════════════════════════

from solver.crystals import RELATIONSHIP_SIGNATURES

print("=" * 80)
print("COMPOSITION OPERATOR SPECTRAL GAP")
print("=" * 80)

corr_types = dict(RELATIONSHIP_SIGNATURES)
# Add S₃ elements
corr_types["identity_perm"] = torch.eye(H, dtype=torch.float32)
corr_types["swap(0,2)"] = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)
corr_types["cycle(012)"] = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)

for name, corr in corr_types.items():
    ent = Entangler(corr, seed=42).build()
    T = build_composition_operator(ent.joint)

    eigenvalues = torch.linalg.eigvals(T)
    mags = eigenvalues.abs()
    idx = mags.argsort(descending=True)
    eigenvalues = eigenvalues[idx]
    mags_sorted = mags[idx]

    # Find the spectral gap: where does the first significant drop occur?
    # Group eigenvalues by magnitude (within 1% tolerance)
    clusters = []
    current_cluster = [mags_sorted[0].item()]
    for i in range(1, len(mags_sorted)):
        m = mags_sorted[i].item()
        if m > 0.99 * current_cluster[-1] and m > 1e-8:
            current_cluster.append(m)
        elif m > 1e-8:
            clusters.append(current_cluster)
            current_cluster = [m]
    if current_cluster and current_cluster[0] > 1e-8:
        clusters.append(current_cluster)

    leading_mag = clusters[0][0] if clusters else 0
    leading_mult = len(clusters[0]) if clusters else 0
    second_mag = clusters[1][0] if len(clusters) > 1 else 0
    second_mult = len(clusters[1]) if len(clusters) > 1 else 0
    gap_ratio = leading_mag / second_mag if second_mag > 0 else float('inf')
    gap_log = math.log(gap_ratio) if gap_ratio > 0 and gap_ratio < float('inf') else 0

    print(f"\n  {name}:")
    print(f"    λ₀ = {leading_mag:.6f} (mult ≈ {leading_mult})")
    print(f"    λ₁ = {second_mag:.6f} (mult ≈ {second_mult})")
    print(f"    gap ratio = {gap_ratio:.4f}")
    print(f"    gap (log) = {gap_log:.4f} (Δ = {DELTA:.4f}, 2Δ = {2*DELTA:.4f})")
    print(f"    gap/Δ = {gap_log/DELTA:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Why 4-fold degeneracy?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("WHY 4-FOLD DEGENERACY?")
print("="*80)

# Hypothesis: the leading eigenvalue corresponds to the 4 rows
# of the joint mass being independently preserved.
# If T_B acts row-by-row: T_B(A)_{i,:} = f_i(A_{i,:}) for some f_i,
# then each row gives one leading eigenvalue.

# Test: build T_B for identity and check if it's block-diagonal in rows
ent_id = Entangler(torch.eye(H, dtype=torch.float32), seed=42).build()
T_id = build_composition_operator(ent_id.joint)

# The flattening is row-major: index = i * MASS_DIM + j
# Row i corresponds to indices [i*4, i*4+1, i*4+2, i*4+3]
print("\n  Is T_B block-diagonal in rows?")
for i in range(MASS_DIM):
    for k in range(MASS_DIM):
        if k == i:
            continue
        # Check if row-i block of T has any connection to row-k block
        row_i = slice(i * MASS_DIM, (i+1) * MASS_DIM)
        row_k = slice(k * MASS_DIM, (k+1) * MASS_DIM)
        cross = T_id[row_i, row_k].abs().max().item()
        if cross > 0.01:
            print(f"    Block ({i},{k}): max = {cross:.4f} (NOT zero)")
            break
    else:
        continue
    break

# Actually compute: for each row block, what's the leading eigenvalue?
print("\n  Row-block analysis of T_B (identity crystal):")
for i in range(MASS_DIM):
    row_block = slice(i * MASS_DIM, (i+1) * MASS_DIM)
    # T_B restricted to row i → row i
    T_block = T_id[row_block, row_block]
    ev = torch.linalg.eigvals(T_block)
    mags = ev.abs()
    idx = mags.argsort(descending=True)
    ev = ev[idx]

    label = f"h{i}" if i < H else "θ"
    print(f"    Row {label}: eigenvalues = ", end="")
    for j in range(min(4, len(ev))):
        print(f"|λ|={ev[j].abs().item():.4f} ", end="")
    print()


# ═══════════════════════════════════════════════════════════════
#  The fixed point of T_B
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("FIXED POINT OF COMPOSITION OPERATOR")
print("="*80)

# T_B^n(A) → fixed point as n → ∞
# This should be the rank-1 product state from the compositional fixed point discovery

A = ent_id.joint.clone()
B = ent_id.joint.clone()

for step in range(30):
    A = compose(A, B)

print(f"  After 30 self-compositions:")
print(f"    Schmidt = {schmidt_number(A):.6f}")

# Check if A is an eigenvector of T_B with eigenvalue λ₀
T_A = compose(A, B)
A_flat = A.reshape(-1)
T_A_flat = T_A.reshape(-1)
valid = A_flat.abs() > 1e-8
if valid.any():
    ratios = T_A_flat[valid] / A_flat[valid]
    print(f"    T_B(A_∞)/A_∞: mean |ratio| = {ratios.abs().mean().item():.6f}, std = {ratios.abs().std().item():.6f}")
    print(f"    This should be ≈ λ₀ = 0.975")


# ═══════════════════════════════════════════════════════════════
#  Connection: 16-dim gap → 4-dim gap
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("16-DIM vs 4-DIM SPECTRAL GAP")
print("="*80)

# The 4×4 joint mass eigenvalues under self-composition:
# We measured r₁/Δ ≈ 2.0-2.2 for bijections
# The 16×16 operator gap is -ln(λ₁/λ₀) / Δ

# Are these the same thing?

ent = Entangler(torch.eye(H, dtype=torch.float32), seed=42).build()
joint = ent.joint

# 4-dim eigenvalues
ev4 = torch.linalg.eigvals(joint)
idx4 = ev4.abs().argsort(descending=True)
ev4 = ev4[idx4]

# After one composition
composed = compose(joint, joint)
ev4c = torch.linalg.eigvals(composed)
idx4c = ev4c.abs().argsort(descending=True)
ev4c = ev4c[idx4c]

# Decay rate of |λ₁/λ₀| per composition
for i in range(1, MASS_DIM):
    r_orig = ev4[i].abs().item() / ev4[0].abs().item()
    r_comp = ev4c[i].abs().item() / ev4c[0].abs().item()
    if r_orig > 1e-12 and r_comp > 1e-12:
        rate = math.log(r_orig / r_comp)
        print(f"  4-dim: decay of |λ_{i}/λ₀| per step = {rate:.4f} = {rate/DELTA:.2f}Δ")

# 16-dim spectral gap
T = build_composition_operator(joint)
ev16 = torch.linalg.eigvals(T)
mags16 = ev16.abs()
idx16 = mags16.argsort(descending=True)
mags16_sorted = mags16[idx16]

# Leading and second cluster
leading = mags16_sorted[0].item()
# Find first eigenvalue significantly below leading
for i in range(1, len(mags16_sorted)):
    if mags16_sorted[i].item() < 0.95 * leading:
        second = mags16_sorted[i].item()
        print(f"\n  16-dim: λ₀ = {leading:.6f}, λ_gap = {second:.6f}")
        gap16 = math.log(leading / second)
        print(f"  16-dim: spectral gap = {gap16:.4f} = {gap16/DELTA:.2f}Δ")
        break

print(f"\n  Comparison:")
print(f"    4-dim r₁ decay rate ≈ 2Δ (composition of 4×4 eigenvalues)")
print(f"    16-dim spectral gap ≈ ?Δ (eigenvalues of 16×16 operator)")
print(f"    Δ = {DELTA:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Seed-averaged 16-dim spectral gap
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SEED-AVERAGED 16-DIM SPECTRAL GAP (identity crystal)")
print("="*80)

gaps = []
for seed in range(50):
    ent = Entangler(torch.eye(H, dtype=torch.float32), seed=seed).build()
    T = build_composition_operator(ent.joint)
    ev = torch.linalg.eigvals(T)
    mags = ev.abs()
    idx = mags.argsort(descending=True)
    mags = mags[idx]

    # Leading eigenvalue
    l0 = mags[0].item()
    # Find spectral gap
    for i in range(1, len(mags)):
        if mags[i].item() < 0.95 * l0:
            gap = math.log(l0 / mags[i].item())
            gaps.append(gap)
            break

avg_gap = sum(gaps) / len(gaps)
std_gap = (sum((g - avg_gap)**2 for g in gaps) / len(gaps)) ** 0.5
print(f"  Gap = {avg_gap:.4f} ± {std_gap:.4f}")
print(f"  Gap/Δ = {avg_gap/DELTA:.3f} ± {std_gap/DELTA:.3f}")
print(f"  2Δ = {2*DELTA:.4f}")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
