"""
The eigenvalue gap of the composition operator as a function of α.

The composition operator L_α(X) = compose_α(A, X) acts on the space
of 4×4 matrices. Its eigenvalue gap determines the rate of decorrelation.

At α=0: L₀(X) = A·X (raw multiplication), eigenvalues are products of A's eigenvalues
At α=1: L₁(X) = normalize(A·X), the structural filter

The mass gap 2Δ should appear as the eigenvalue gap of L₁.
Question: how does this gap depend on α?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, EPS_LOG)
from solver.crystals import Entangler

torch.set_grad_enabled(False)

n_seeds = 500

corrs = {
    "A": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "B": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
}

crystals = {}
for name, corr in corrs.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    crystals[name] = avg

A = crystals["A"]


def compose_alpha(X, Y, alpha):
    raw = torch.einsum("ab,bc->ac", X, Y)
    L1 = raw.reshape(-1).real.sum()
    if abs(L1) < 1e-15:
        return raw
    return raw / (abs(L1) ** alpha)


def composition_operator_matrix(crystal, alpha, dim=MASS_DIM):
    """Build the (dim²×dim²) matrix of the linear map X ↦ compose_α(crystal, X).

    Since compose_α involves L1(crystal·X) which depends on X,
    the map is NOT linear for α>0. But we can linearize around
    a reference point (like the identity matrix).

    For α=0 it IS linear: L₀(X) = A·X.
    """
    N = dim * dim
    L = torch.zeros(N, N, dtype=torch.cfloat)

    for j in range(N):
        # Basis matrix e_j
        ej = torch.zeros(dim, dim, dtype=torch.cfloat)
        ej[j // dim, j % dim] = 1.0

        # Apply composition
        result = compose_alpha(crystal, ej, alpha)

        # Flatten result
        L[:, j] = result.reshape(-1)

    return L


# ═══════════════════════════════════════════════════════════════
#  Eigenvalue spectrum at α=0 (raw multiplication)
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("COMPOSITION OPERATOR EIGENVALUES")
print("=" * 80)

L0 = composition_operator_matrix(A, 0.0)
eigvals_0 = torch.linalg.eigvals(L0)
eigvals_0_abs = eigvals_0.abs()
eigvals_0_sorted, _ = eigvals_0_abs.sort(descending=True)

print(f"\n  α=0 (raw multiplication, X ↦ A·X):")
print(f"    Top eigenvalues (absolute):")
for i in range(min(8, len(eigvals_0_sorted))):
    print(f"      λ{i} = {eigvals_0_sorted[i].item():.8f}")

gap_0 = -math.log(eigvals_0_sorted[1].item() / eigvals_0_sorted[0].item()) \
    if eigvals_0_sorted[0].item() > 1e-10 and eigvals_0_sorted[1].item() > 1e-10 else 0
print(f"    -ln(λ₁/λ₀) = {gap_0:.6f}")


# ═══════════════════════════════════════════════════════════════
#  Eigenvalue gap vs α
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("EIGENVALUE GAP vs α")
print("="*80)

# Note: for α>0, the composition X ↦ A·X/||A·X||₁^α is NOT linear in X
# because the normalization depends on X. So the "composition operator"
# is only meaningful at α=0.
#
# However, we can look at the ITERATED composition spectrum:
# How fast does compose_α(A, compose_α(A, ... X ...)) converge?
# This is measured by the contraction rate of the iterated map.

print(f"\n  The composition operator is linear ONLY at α=0.")
print(f"  For α>0, we measure the contraction rate empirically.")
print(f"\n  At α=0: the eigenvalues of A (the crystal matrix) determine the map.")

# Eigenvalues of A itself
eigvals_A = torch.linalg.eigvals(A)
eigvals_A_abs = eigvals_A.abs()
eigvals_A_sorted, _ = eigvals_A_abs.sort(descending=True)

print(f"\n  Eigenvalues of crystal A:")
for i in range(MASS_DIM):
    ev = eigvals_A[i]
    print(f"    λ{i} = {ev.real.item():.6f} + {ev.imag.item():.6f}i  (|λ| = {abs(ev.item()):.6f})")

gap_A = -math.log(eigvals_A_sorted[1].item() / eigvals_A_sorted[0].item()) \
    if eigvals_A_sorted[0].item() > 1e-10 and eigvals_A_sorted[1].item() > 1e-10 else 0
print(f"\n  Crystal eigenvalue gap: -ln(|λ₁|/|λ₀|) = {gap_A:.6f}")
print(f"  Δ = {DELTA:.6f}")
print(f"  2Δ = {2*DELTA:.6f}")
print(f"  Gap/Δ = {gap_A/DELTA:.6f}")
print(f"  Gap/2Δ = {gap_A/(2*DELTA):.6f}")


# ═══════════════════════════════════════════════════════════════
#  Contraction rate of iterated composition
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CONTRACTION RATE OF ITERATED COMPOSITION")
print("="*80)

# Start with B, repeatedly compose with A
B = crystals["B"]

from solver.crystals import compose

def norm(X):
    return X.abs().pow(2).sum().sqrt().item()


# At α=1 (full normalization, the actual compose())
print(f"\n  α=1.0 (full normalization):")
X = B.clone()
prev_diff = None
for step in range(10):
    X_new = compose(A, X)
    diff = norm(X_new - compose(A, compose(A, X)))  # second-order convergence check

    # Track how X evolves
    eigvals_X = torch.linalg.eigvals(X_new)
    eigvals_X_abs = eigvals_X.abs()
    eigvals_X_sorted, _ = eigvals_X_abs.sort(descending=True)

    schmidt = 1.0 / (eigvals_X_sorted / eigvals_X_sorted.sum()).pow(2).sum().item() \
        if eigvals_X_sorted.sum().item() > 1e-10 else 1.0

    if prev_diff is not None and prev_diff > 1e-10:
        contraction = diff / prev_diff
        ln_contraction = -math.log(contraction) if contraction > 1e-10 else float('inf')
    else:
        contraction = 0
        ln_contraction = 0

    print(f"    step {step}: schmidt={schmidt:.4f}, "
          f"|λ₁/λ₀|={eigvals_X_sorted[1].item()/eigvals_X_sorted[0].item():.6f}")

    prev_diff = diff
    X = X_new


# ═══════════════════════════════════════════════════════════════
#  Direct: eigenvalue gap of self-composition
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("EIGENVALUE GAP UNDER SELF-COMPOSITION")
print("="*80)

# A, A², A³, ...  (using normalized composition)
print(f"\n  Self-composition of crystal A:")
X = A.clone()
for step in range(8):
    eigvals_X = torch.linalg.eigvals(X)
    eigvals_X_abs = eigvals_X.abs()
    eigvals_X_sorted, _ = eigvals_X_abs.sort(descending=True)
    gap = eigvals_X_sorted[1].item() / eigvals_X_sorted[0].item() \
        if eigvals_X_sorted[0].item() > 1e-10 else 0
    ln_gap = -math.log(gap) if gap > 1e-10 else float('inf')

    print(f"    A^{step+1}: |λ₁/λ₀| = {gap:.8f}, -ln(|λ₁/λ₀|) = {ln_gap:.6f}, "
          f"ratio to 2Δ = {ln_gap/(2*DELTA):.4f}")

    X_new = compose(X, A)
    X = X_new


# ═══════════════════════════════════════════════════════════════
#  All S₃ element eigenvalue gaps
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("EIGENVALUE GAPS OF ALL S₃ ELEMENT CRYSTALS")
print("="*80)

all_corrs = {
    "e": torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)": torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

all_crystals = {}
for name, corr in all_corrs.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    all_crystals[name] = avg

print(f"\n  {'Name':>6s}  {'|λ₀|':>10s}  {'|λ₁|':>10s}  {'|λ₁/λ₀|':>10s}  "
      f"{'gap':>8s}  {'gap/Δ':>8s}  {'gap/2Δ':>8s}")
print("-" * 80)

for name, crystal in all_crystals.items():
    eigvals_c = torch.linalg.eigvals(crystal)
    eigvals_c_abs = eigvals_c.abs()
    eigvals_c_sorted, _ = eigvals_c_abs.sort(descending=True)

    ratio = eigvals_c_sorted[1].item() / eigvals_c_sorted[0].item() \
        if eigvals_c_sorted[0].item() > 1e-10 else 0
    gap = -math.log(ratio) if ratio > 1e-10 else float('inf')

    print(f"  {name:>6s}  {eigvals_c_sorted[0].item():>10.6f}  "
          f"{eigvals_c_sorted[1].item():>10.6f}  {ratio:>10.6f}  "
          f"{gap:>8.4f}  {gap/DELTA:>8.4f}  {gap/(2*DELTA):>8.4f}")


# ═══════════════════════════════════════════════════════════════
#  The composition operator at α=0 eigenvalue structure
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMPOSITION OPERATOR L₀(X) = A·X: FULL 16×16 SPECTRUM")
print("="*80)

L0 = composition_operator_matrix(A, 0.0)
eigvals_full = torch.linalg.eigvals(L0)
eigvals_full_abs = eigvals_full.abs()
eigvals_full_sorted, idx = eigvals_full_abs.sort(descending=True)

print(f"\n  16 eigenvalues of L₀ (grouped by magnitude):")
for i in range(16):
    ev = eigvals_full[idx[i]]
    print(f"    λ{i:2d} = {ev.real.item():>10.6f} + {ev.imag.item():>10.6f}i  "
          f"(|λ| = {eigvals_full_sorted[i].item():.8f})")

# The eigenvalues of L₀(X) = A·X are just the eigenvalues of A,
# each with multiplicity 4 (one for each column of X)
print(f"\n  Expected: eigenvalues of A (each ×4 multiplicity)")
print(f"  Eigenvalues of A: {[f'{x:.6f}' for x in eigvals_A_sorted.tolist()]}")


print(f"\n\n{'='*80}")
print("WHAT THE EIGENVALUE GAP REVEALS")
print("="*80)
