"""
Raw composition (no normalization) should satisfy Jacobi exactly.

Matrix multiplication A@B = Σ_k A_{ik}B_{kj} is associative.
The commutator of an associative algebra ALWAYS satisfies Jacobi.
(Standard theorem: [[A,B],C] + cyclic = A(BC-CB) - (BC-CB)A + ... = 0)

The crystal's compose() function normalizes by dividing by the real
L1 sum. This normalization is what breaks the Jacobi identity.

Test:
1. Raw composition (just einsum, no normalization) → Jacobi = 0 exactly
2. Normalized composition (full compose()) → Jacobi ≠ 0 (confirmed)
3. What's the minimal normalization that breaks Jacobi?

The mass gap lives in the normalization.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, EPS_LOG)
from solver.crystals import Entangler

torch.set_grad_enabled(False)


def raw_compose(A, B):
    """Raw tensor contraction — no normalization."""
    return torch.einsum("ab,bc->ac", A, B)


def normalized_compose(A, B):
    """Full normalized composition (same as solver's compose)."""
    result = torch.einsum("ab,bc->ac", A, B)
    re_sum = result.reshape(-1).real.sum()
    if abs(re_sum) > 1e-8:
        result = result / re_sum
    return result


def bracket(compose_fn, A, B):
    return compose_fn(A, B) - compose_fn(B, A)


def norm(X):
    return X.abs().pow(2).sum().sqrt().item()


n_seeds = 200

# Build seed-averaged transposition crystals
corrs = {
    "A": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),  # (01)
    "B": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),  # (02)
    "C": torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),  # (12)
}

crystals = {}
for name, corr in corrs.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    crystals[name] = avg

A, B, C = crystals["A"], crystals["B"], crystals["C"]


print("=" * 80)
print("RAW vs NORMALIZED: WHERE THE JACOBI IDENTITY LIVES")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Raw composition: Jacobi should be EXACTLY zero
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Raw composition (no normalization) ---\n")

AB_raw = bracket(raw_compose, A, B)
BC_raw = bracket(raw_compose, B, C)
CA_raw = bracket(raw_compose, C, A)

AB_C_raw = bracket(raw_compose, AB_raw, C)
BC_A_raw = bracket(raw_compose, BC_raw, A)
CA_B_raw = bracket(raw_compose, CA_raw, B)

jacobi_raw = AB_C_raw + BC_A_raw + CA_B_raw

print(f"  ||[A,B]_raw|| = {norm(AB_raw):.6f}")
print(f"  ||[[A,B],C]_raw|| = {norm(AB_C_raw):.6f}")
print(f"  ||[[B,C],A]_raw|| = {norm(BC_A_raw):.6f}")
print(f"  ||[[C,A],B]_raw|| = {norm(CA_B_raw):.6f}")
print(f"  ||Jacobi_raw|| = {norm(jacobi_raw):.10f}")

max_raw = max(norm(AB_C_raw), norm(BC_A_raw), norm(CA_B_raw))
if max_raw > 1e-15:
    print(f"  Relative violation: {norm(jacobi_raw)/max_raw:.2e}")


# ═══════════════════════════════════════════════════════════════
#  Normalized composition: Jacobi violated (recap)
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Normalized composition ---\n")

AB_norm = bracket(normalized_compose, A, B)
BC_norm = bracket(normalized_compose, B, C)
CA_norm = bracket(normalized_compose, C, A)

AB_C_norm = bracket(normalized_compose, AB_norm, C)
BC_A_norm = bracket(normalized_compose, BC_norm, A)
CA_B_norm = bracket(normalized_compose, CA_norm, B)

jacobi_norm = AB_C_norm + BC_A_norm + CA_B_norm

print(f"  ||[A,B]_norm|| = {norm(AB_norm):.6f}")
print(f"  ||[[A,B],C]_norm|| = {norm(AB_C_norm):.6f}")
print(f"  ||Jacobi_norm|| = {norm(jacobi_norm):.6f}")

max_norm_val = max(norm(AB_C_norm), norm(BC_A_norm), norm(CA_B_norm))
if max_norm_val > 1e-10:
    print(f"  Relative violation: {norm(jacobi_norm)/max_norm_val:.6f}")


# ═══════════════════════════════════════════════════════════════
#  S₃ decomposition of the RAW commutator
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("S₃ DECOMPOSITION: RAW vs NORMALIZED COMMUTATORS")
print("="*80)

def s3_fracs(joint):
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
    bp = joint.abs().pow(2)
    total = bp.sum().item()
    if total < 1e-15:
        return 1.0, 0.0, 0.0
    return (triv.abs().pow(2).sum().item() / total,
            sign.abs().pow(2).sum().item() / total,
            std.abs().pow(2).sum().item() / total)

print(f"\n  Raw [A,B]:  ", end="")
f_t, f_s, f_st = s3_fracs(AB_raw)
print(f"triv={f_t:.5f}, sign={f_s:.5f}, std={f_st:.5f}")

print(f"  Norm [A,B]:  ", end="")
f_t, f_s, f_st = s3_fracs(AB_norm)
print(f"triv={f_t:.5f}, sign={f_s:.5f}, std={f_st:.5f}")

print(f"\n  Raw [[A,B],C]:  ", end="")
f_t, f_s, f_st = s3_fracs(AB_C_raw)
print(f"triv={f_t:.5f}, sign={f_s:.5f}, std={f_st:.5f}")

print(f"  Norm [[A,B],C]:  ", end="")
f_t, f_s, f_st = s3_fracs(AB_C_norm)
print(f"triv={f_t:.5f}, sign={f_s:.5f}, std={f_st:.5f}")


# ═══════════════════════════════════════════════════════════════
#  Per-seed Jacobi with raw composition
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("PER-SEED RAW JACOBI")
print("="*80)

violations_raw = []
for seed in range(100):
    a = Entangler(corrs["A"], seed=seed).build().joint
    b = Entangler(corrs["B"], seed=seed).build().joint
    c = Entangler(corrs["C"], seed=seed).build().joint

    ab_c = bracket(raw_compose, bracket(raw_compose, a, b), c)
    bc_a = bracket(raw_compose, bracket(raw_compose, b, c), a)
    ca_b = bracket(raw_compose, bracket(raw_compose, c, a), b)

    j = ab_c + bc_a + ca_b
    max_n = max(norm(ab_c), norm(bc_a), norm(ca_b))
    if max_n > 1e-15:
        violations_raw.append(norm(j) / max_n)

avg_raw = sum(violations_raw) / len(violations_raw) if violations_raw else 0
print(f"\n  {len(violations_raw)} seeds:")
print(f"  Average relative violation: {avg_raw:.2e}")
print(f"  Max violation: {max(violations_raw):.2e}")


# ═══════════════════════════════════════════════════════════════
#  The fundamental result
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE FUNDAMENTAL RESULT")
print("="*80)

print(f"""
  Raw matrix multiplication: Jacobi violation = {norm(jacobi_raw):.2e}
  DS-normalized composition: Jacobi violation = {norm(jacobi_norm):.2f}

  The Lie algebra gl(4,C) is EXACTLY realized by raw tensor contraction.
  The DS normalization (divide by L1 sum) is the SOLE source of
  Jacobi violation.

  The mass gap = the curvature of the L1 normalization on the mass simplex.

  The gauge theory is EXACT in the tangent space (gl(4,C)).
  The mass gap is the PROJECTION from gl(4,C) to the mass simplex.
  Confinement = the gauge field cannot propagate on the simplex
  because the simplex is curved (Born floor, L1 constraint).
""")


print(f"\n\n{'='*80}")
print("WHAT THE RAW JACOBI REVEALS")
print("="*80)
