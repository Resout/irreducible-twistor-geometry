"""
Analytical proof: d||Jacobi||/dα = 2Δ at α=0.

The composition with partial normalization:
  compose_α(A,B) = A·B / ||A·B||₁^α

At α=0: raw matrix multiplication → Jacobi = 0 exactly.
At any α>0: normalization breaks Jacobi.

The question: what is the RATE of Jacobi violation emergence?
If this rate equals 2Δ (twice the mass gap), then:
  - The mass gap literally IS the rate at which gauge symmetry breaks
  - Confinement = Jacobi violation = normalization curvature

Strategy:
  1. Compute d/dα ||Jacobi(α)|| numerically at α→0
  2. Identify the analytical formula
  3. Show it equals 2Δ = -2 ln(1-K*)
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, EPS_LOG)
from solver.crystals import Entangler

torch.set_grad_enabled(False)


def compose_alpha(A, B, alpha):
    """Composition with partial normalization."""
    raw = torch.einsum("ab,bc->ac", A, B)
    L1 = raw.reshape(-1).real.sum()
    if abs(L1) < 1e-15:
        return raw
    return raw / (abs(L1) ** alpha)


def bracket_alpha(A, B, alpha):
    return compose_alpha(A, B, alpha) - compose_alpha(B, A, alpha)


def norm(X):
    return X.abs().pow(2).sum().sqrt().item()


n_seeds = 500  # high precision

# Build seed-averaged transposition crystals
corrs = {
    "A": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "B": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "C": torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
}

crystals = {}
for name, corr in corrs.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    crystals[name] = avg

A, B, C = crystals["A"], crystals["B"], crystals["C"]

# ═══════════════════════════════════════════════════════════════
#  L1 norms of the products (needed for analytical formula)
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("L1 NORMS OF RAW PRODUCTS")
print("=" * 80)

products = {}
for pair in ["AB", "BA", "AC", "CA", "BC", "CB"]:
    X, Y = crystals[pair[0]], crystals[pair[1]]
    raw = torch.einsum("ab,bc->ac", X, Y)
    L1 = raw.reshape(-1).real.sum().item()
    products[pair] = {"raw": raw, "L1": L1}
    print(f"  L1(raw {pair}) = {L1:.8f}")

print(f"\n  ln(L1_AB) = {math.log(abs(products['AB']['L1'])):.8f}")
print(f"  ln(L1_BA) = {math.log(abs(products['BA']['L1'])):.8f}")


# ═══════════════════════════════════════════════════════════════
#  Numerical derivative at α→0 (Richardson extrapolation)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("NUMERICAL DERIVATIVE d||Jacobi||/dα AT α→0")
print("="*80)

def jacobi_norm_at(alpha):
    AB = bracket_alpha(A, B, alpha)
    BC = bracket_alpha(B, C, alpha)
    CA = bracket_alpha(C, A, alpha)
    AB_C = bracket_alpha(AB, C, alpha)
    BC_A = bracket_alpha(BC, A, alpha)
    CA_B = bracket_alpha(CA, B, alpha)
    jacobi = AB_C + BC_A + CA_B
    return norm(jacobi)


# d/dα at α=0 using increasingly small step sizes
print(f"\n  Step-size   d||J||/dα      ratio to 2Δ={2*DELTA:.6f}")
print("-" * 60)

derivs = []
for h_exp in range(1, 9):
    h = 10 ** (-h_exp)
    j_h = jacobi_norm_at(h)
    j_0 = jacobi_norm_at(0)
    deriv = (j_h - j_0) / h
    ratio = deriv / (2 * DELTA) if deriv > 0 else 0
    derivs.append((h, deriv, ratio))
    print(f"  {h:.0e}     {deriv:>14.8f}  {ratio:>12.6f}")


# Also try Richardson extrapolation
print(f"\n--- Richardson extrapolation ---")
h = 1e-4
d1 = (jacobi_norm_at(h) - jacobi_norm_at(0)) / h
d2 = (jacobi_norm_at(h/2) - jacobi_norm_at(0)) / (h/2)
richardson = (4*d2 - d1) / 3
print(f"  d1(h={h:.0e}) = {d1:.10f}")
print(f"  d2(h={h/2:.0e}) = {d2:.10f}")
print(f"  Richardson = {richardson:.10f}")
print(f"  2Δ = {2*DELTA:.10f}")
print(f"  ratio = {richardson / (2*DELTA):.10f}")


# ═══════════════════════════════════════════════════════════════
#  What about relative violation derivative?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("RELATIVE VIOLATION DERIVATIVE")
print("="*80)

def relative_jacobi_at(alpha):
    AB = bracket_alpha(A, B, alpha)
    BC = bracket_alpha(B, C, alpha)
    CA = bracket_alpha(C, A, alpha)
    AB_C = bracket_alpha(AB, C, alpha)
    BC_A = bracket_alpha(BC, A, alpha)
    CA_B = bracket_alpha(CA, B, alpha)
    jacobi = AB_C + BC_A + CA_B
    j_norm = norm(jacobi)
    max_db = max(norm(AB_C), norm(BC_A), norm(CA_B))
    return j_norm / max_db if max_db > 1e-15 else 0

# Actually: try log derivative d(ln||J||)/dα
print(f"\n--- Logarithmic derivative d(ln||J||)/dα ---")
print(f"\n  {'α':>6s}  {'ln||J||':>12s}  {'d(ln||J||)/dα':>14s}  {'ratio/2Δ':>10s}")
print("-" * 50)

prev_lnJ = None
for alpha in [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]:
    j = jacobi_norm_at(alpha)
    lnJ = math.log(j) if j > 0 else float('-inf')

    if prev_lnJ is not None and prev_lnJ[1] > float('-inf'):
        d_lnJ = (lnJ - prev_lnJ[1]) / (alpha - prev_lnJ[0])
        ratio = d_lnJ / (2 * DELTA)
        print(f"  {alpha:>6.3f}  {lnJ:>12.4f}  {d_lnJ:>14.4f}  {ratio:>10.4f}")
    else:
        print(f"  {alpha:>6.3f}  {lnJ:>12.4f}")

    prev_lnJ = (alpha, lnJ)


# ═══════════════════════════════════════════════════════════════
#  Per-element Jacobi: decompose into irrep channels
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("JACOBI VIOLATION BY S₃ IRREP CHANNEL")
print("="*80)

def s3_project(joint):
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]
    signs = [1, -1, -1, -1, 1, 1]
    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T
    triv = sum(permute(joint, p) for p in perms) / 6
    sign = sum(s * permute(joint, p) for p, s in zip(perms, signs)) / 6
    std = joint - triv - sign
    return triv, sign, std


print(f"\n  {'α':>6s}  {'||J_triv||':>12s}  {'||J_sign||':>12s}  {'||J_std||':>12s}  {'triv%':>8s}")

for alpha in [0.01, 0.05, 0.1, 0.15, 0.2, K_STAR, 0.3, 0.5, 0.7, 1.0]:
    AB = bracket_alpha(A, B, alpha)
    BC = bracket_alpha(B, C, alpha)
    CA = bracket_alpha(C, A, alpha)
    AB_C = bracket_alpha(AB, C, alpha)
    BC_A = bracket_alpha(BC, A, alpha)
    CA_B = bracket_alpha(CA, B, alpha)
    jacobi = AB_C + BC_A + CA_B

    triv, sign, std = s3_project(jacobi)
    n_triv = norm(triv)
    n_sign = norm(sign)
    n_std = norm(std)
    total = n_triv + n_sign + n_std

    print(f"  {alpha:>6.3f}  {n_triv:>12.4e}  {n_sign:>12.4e}  {n_std:>12.4e}  "
          f"{100*n_triv/total if total > 0 else 0:>7.1f}%")


# ═══════════════════════════════════════════════════════════════
#  The key insight: log-log slope
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("LOG-LOG SCALING: ||Jacobi|| ~ α^n at small α")
print("="*80)

print(f"\n  {'α':>8s}  {'||J||':>14s}  {'ln||J||':>10s}  {'ln(α)':>10s}  {'slope':>8s}")
print("-" * 60)

prev = None
for alpha in [1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 0.01, 0.02, 0.05]:
    j = jacobi_norm_at(alpha)
    lnj = math.log(j) if j > 1e-30 else -70
    lna = math.log(alpha)

    if prev is not None:
        slope = (lnj - prev[0]) / (lna - prev[1])
    else:
        slope = float('nan')

    prev = (lnj, lna)
    print(f"  {alpha:>8.1e}  {j:>14.6e}  {lnj:>10.3f}  {lna:>10.3f}  {slope:>8.3f}")


# ═══════════════════════════════════════════════════════════════
#  Connection to eigenvalue decay
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CONNECTION TO EIGENVALUE SPECTRUM")
print("="*80)

# The composition eigenvalues decay as exp(-2Δ) per step
# If Jacobi violation grows as α² at small α,
# then d²||J||/dα² at 0 = 2 × (some constant)
# The constant should be related to the bracket norms

raw_bracket_norms = []
for pair_label, (name_a, name_b) in [("AB", ("A","B")), ("BC", ("B","C")), ("CA", ("C","A"))]:
    X, Y = crystals[name_a], crystals[name_b]
    raw_br = torch.einsum("ab,bc->ac", X, Y) - torch.einsum("ab,bc->ac", Y, X)
    raw_bracket_norms.append(norm(raw_br))
    print(f"  ||[{name_a},{name_b}]_raw|| = {norm(raw_br):.8f}")

print(f"\n  avg ||bracket||_raw = {sum(raw_bracket_norms)/3:.8f}")

# L1 norms
L1_vals = []
for pair in ["AB", "BA", "AC", "CA", "BC", "CB"]:
    L1_vals.append(abs(products[pair]["L1"]))
print(f"  avg L1 = {sum(L1_vals)/6:.8f}")
print(f"  ln(avg L1) = {math.log(sum(L1_vals)/6):.8f}")

# The key ratio
print(f"\n  2Δ = {2*DELTA:.8f}")
print(f"  K* = {K_STAR:.8f}")
print(f"  -ln(1-K*) = {-math.log(1-K_STAR):.8f}")
print(f"  -2ln(1-K*) = {-2*math.log(1-K_STAR):.8f}")

# Check: L1_AB / L1_BA
print(f"\n  L1_AB / L1_BA = {abs(products['AB']['L1']) / abs(products['BA']['L1']):.8f}")
print(f"  All L1 ratios:")
for p1 in ["AB", "BC", "CA"]:
    p2 = p1[1] + p1[0]
    r = abs(products[p1]["L1"]) / abs(products[p2]["L1"])
    print(f"    L1({p1})/L1({p2}) = {r:.8f}")


print(f"\n\n{'='*80}")
print("WHAT d(JACOBI)/dα REVEALS")
print("="*80)
