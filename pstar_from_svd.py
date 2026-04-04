"""
Why p* = exp(Δ/2)? — Deriving the fixed point from the SVD.

The composition operator T_C(A) = A @ C maps [4,4] → [4,4].
For a product state u⊗v: T_C(u⊗v) = u ⊗ (C^T v).
So the fixed point is u⊗v where v is the dominant eigenvector of C^T.

The marginal Born(θ) of the product state u⊗v is just Born_u(θ).
And u is determined by the row-marginal of the crystal.

Questions:
1. What is Born(θ) of the dominant right singular vector of C?
2. What is Born(θ) of the row marginal of C?
3. Do these give 10/33?
4. Can we derive 10/33 from K* = 7/30?
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

n_seeds = 50


print("=" * 80)
print("WHY p* = exp(Δ/2)? — THE SVD STRUCTURE")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  SVD of each S₃ crystal
# ═══════════════════════════════════════════════════════════════

print(f"\n--- SVD of S₃ crystals (50-seed averaged) ---\n")

for name, corr in S3_CORRS.items():
    # Build seed-averaged crystal
    crystal = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        crystal += Entangler(corr, seed=seed).build().joint
    crystal /= n_seeds

    # SVD
    U, S, Vh = torch.linalg.svd(crystal)

    # Singular values
    S_norm = S / S.sum()
    schmidt = 1.0 / (S_norm**2).sum().item()

    print(f"  {name:>6s}: σ = [{', '.join(f'{s:.4f}' for s in S.tolist())}]")
    print(f"         σ_norm = [{', '.join(f'{s:.4f}' for s in S_norm.tolist())}]")

    # Dominant left singular vector (column of U)
    u0 = U[:, 0]
    bp_u0 = u0.abs().pow(2) / u0.abs().pow(2).sum()
    print(f"         u₀ Born = [{', '.join(f'{b:.4f}' for b in bp_u0.tolist())}], Born(θ)={bp_u0[H].item():.5f}")

    # Dominant right singular vector (row of Vh)
    v0 = Vh[0, :]
    bp_v0 = v0.abs().pow(2) / v0.abs().pow(2).sum()
    print(f"         v₀ Born = [{', '.join(f'{b:.4f}' for b in bp_v0.tolist())}], Born(θ)={bp_v0[H].item():.5f}")

    # Row marginal of crystal: sum over columns
    row_marg = crystal.abs().pow(2).sum(dim=1)
    bp_row = row_marg / row_marg.sum()
    print(f"         row_marg Born = [{', '.join(f'{b:.4f}' for b in bp_row.tolist())}], Born(θ)={bp_row[H].item():.5f}")

    # Column marginal: sum over rows
    col_marg = crystal.abs().pow(2).sum(dim=0)
    bp_col = col_marg / col_marg.sum()
    print(f"         col_marg Born = [{', '.join(f'{b:.4f}' for b in bp_col.tolist())}], Born(θ)={bp_col[H].item():.5f}")

    # The actual fixed point from iterated composition
    fp = crystal.clone()
    for _ in range(40):
        fp = compose(fp, crystal)
    bp_fp = born_probabilities(fp).reshape(MASS_DIM, MASS_DIM)
    print(f"         fixed_pt Born(θ) = {bp_fp[0, H].item():.5f}")
    print()


# ═══════════════════════════════════════════════════════════════
#  Eigenvalues of the composition operator
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print("EIGENVALUES OF THE COMPOSITION OPERATOR T_C")
print("="*80)

# Build T_C for the identity crystal
crystal_e = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for seed in range(n_seeds):
    crystal_e += Entangler(S3_CORRS["e"], seed=seed).build().joint
crystal_e /= n_seeds

dim_sq = MASS_DIM * MASS_DIM
T = torch.zeros(dim_sq, dim_sq, dtype=torch.cfloat)
for col in range(dim_sq):
    A = torch.zeros(dim_sq, dtype=torch.cfloat)
    A[col] = 1.0 + 0j
    result = compose(A.reshape(MASS_DIM, MASS_DIM), crystal_e)
    T[:, col] = result.reshape(dim_sq)

eigvals = torch.linalg.eigvals(T)
idx = eigvals.abs().argsort(descending=True)
eigvals_sorted = eigvals[idx]

print(f"\n  Top 8 eigenvalues of T_e (by magnitude):")
for i in range(8):
    ev = eigvals_sorted[i]
    print(f"    λ_{i} = {ev.abs().item():.6f} (phase {torch.angle(ev).item()/math.pi:.4f}π)")

# The ratio λ₁/λ₀ gives the contraction rate per step
if eigvals_sorted[0].abs().item() > 1e-10:
    ratio_10 = eigvals_sorted[1].abs().item() / eigvals_sorted[0].abs().item()
    ratio_40 = eigvals_sorted[4].abs().item() / eigvals_sorted[0].abs().item()
    print(f"\n  λ₁/λ₀ = {ratio_10:.6f}")
    print(f"  λ₄/λ₀ = {ratio_40:.6f}")
    print(f"  exp(-Δ) = {math.exp(-DELTA):.6f}")
    print(f"  exp(-2Δ) = {math.exp(-2*DELTA):.6f}")
    print(f"  (1-K*) = {1-K_STAR:.6f}")
    print(f"  299/405 = {299/405:.6f}")


# ═══════════════════════════════════════════════════════════════
#  The dominant eigenvector of T_C
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("DOMINANT EIGENVECTOR OF T_C")
print("="*80)

eigvals_full, eigvecs = torch.linalg.eig(T)
idx = eigvals_full.abs().argsort(descending=True)

# Dominant eigenvector (the fixed point of T_C)
dom_vec = eigvecs[:, idx[0]]
dom_joint = dom_vec.reshape(MASS_DIM, MASS_DIM)

# Normalize to real L1 = 1
re_sum = dom_joint.real.sum()
if abs(re_sum) > 1e-10:
    dom_joint = dom_joint / re_sum

bp_dom = dom_joint.abs().pow(2)
bp_dom = bp_dom / bp_dom.sum()

print(f"\n  Born matrix of dominant eigenvector:")
for i in range(MASS_DIM):
    print(f"    [{', '.join(f'{bp_dom[i,j].item():.5f}' for j in range(MASS_DIM))}]")

sn_dom = schmidt_number(dom_joint)
print(f"\n  Schmidt = {sn_dom:.4f} (1.0 = product state)")

# Marginal Born(θ)
marg = bp_dom.sum(dim=1)
print(f"  Marginal = [{', '.join(f'{marg[i].item():.5f}' for i in range(MASS_DIM))}]")
print(f"  Born(θ) = {marg[H].item():.5f}")
print(f"  10/33   = {10/33:.5f}")


# ═══════════════════════════════════════════════════════════════
#  The second eigenvector (contraction direction)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SECOND EIGENVECTOR (CONTRACTION DIRECTION)")
print("="*80)

sec_vec = eigvecs[:, idx[1]]
sec_joint = sec_vec.reshape(MASS_DIM, MASS_DIM)

bp_sec = sec_joint.abs().pow(2)
bp_sec = bp_sec / bp_sec.sum()

print(f"\n  Born matrix of second eigenvector:")
for i in range(MASS_DIM):
    print(f"    [{', '.join(f'{bp_sec[i,j].item():.5f}' for j in range(MASS_DIM))}]")

sn_sec = schmidt_number(sec_joint)
print(f"\n  Schmidt = {sn_sec:.4f}")
print(f"  Eigenvalue ratio λ₁/λ₀ = {eigvals_full[idx[1]].abs().item() / eigvals_full[idx[0]].abs().item():.6f}")


# ═══════════════════════════════════════════════════════════════
#  Key test: is the dominant eigenvector a product state?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("IS THE DOMINANT EIGENVECTOR A PRODUCT STATE?")
print("="*80)

# If it's a product state, all rows should be proportional
rows = [bp_dom[i] / bp_dom[i].sum() for i in range(MASS_DIM)]
print(f"\n  Row-normalized Born matrix:")
for i in range(MASS_DIM):
    print(f"    Row {i}: [{', '.join(f'{rows[i][j].item():.5f}' for j in range(MASS_DIM))}]")

# Maximum pairwise difference between rows
max_diff = 0
for i in range(MASS_DIM):
    for j in range(i+1, MASS_DIM):
        diff = (rows[i] - rows[j]).abs().max().item()
        max_diff = max(max_diff, diff)
print(f"\n  Max row difference: {max_diff:.6f} (0 = perfect product state)")


# ═══════════════════════════════════════════════════════════════
#  Connecting the eigenvalue structure to p*
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE EIGENVALUE-p* CONNECTION")
print("="*80)

lambda_0 = eigvals_full[idx[0]].abs().item()
lambda_1 = eigvals_full[idx[1]].abs().item()

print(f"\n  λ₀ = {lambda_0:.6f}")
print(f"  λ₁ = {lambda_1:.6f}")
print(f"  λ₁/λ₀ = {lambda_1/lambda_0:.6f}")
print(f"  -ln(λ₁/λ₀) = {-math.log(lambda_1/lambda_0):.6f}")
print(f"  Δ = {DELTA:.6f}")
print(f"  2Δ = {2*DELTA:.6f}")

# The contraction rate of T_C in the standard sector should be
# related to 2Δ (composition gap from Principle 30).
# And p* = exp(Δ/2) connects the gap to the fixed point.

# Check: does λ₀ = 1/(H+p*²)?
# The dominant eigenvalue should be the "weight" of the product state
p_star = math.exp(DELTA/2)
predicted_lambda0 = 1 / (H + p_star**2)
print(f"\n  Predicted λ₀ = 1/(H+p*²) = 1/(3+30/23) = 23/99 = {predicted_lambda0:.6f}")
print(f"  Measured λ₀ = {lambda_0:.6f}")
print(f"  Ratio: {lambda_0 / predicted_lambda0:.4f}")

# Check if λ₀ relates to the composition rate
print(f"\n  λ₀ × (H+p*²) = {lambda_0 * (H + p_star**2):.6f}")
print(f"  λ₀ × 99/23 = {lambda_0 * 99/23:.6f}")
print(f"  λ₀ × H = {lambda_0 * H:.6f}")
print(f"  λ₀ × H² = {lambda_0 * H**2:.6f}")
print(f"  λ₀ × H³ = {lambda_0 * H**3:.6f}")

# The singular values of the crystal
U, S, Vh = torch.linalg.svd(crystal_e)
print(f"\n  Crystal singular values: [{', '.join(f'{s:.6f}' for s in S.tolist())}]")
print(f"  σ₀² = {S[0].item()**2:.6f}")
print(f"  σ₀ = {S[0].item():.6f}")
print(f"  σ₀ / (1-K*) = {S[0].item() / (1-K_STAR):.6f}")
print(f"  σ₀² / (1-K*) = {S[0].item()**2 / (1-K_STAR):.6f}")


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
