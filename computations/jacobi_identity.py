"""
The Jacobi identity: the final test of a Lie algebra.

[[A,B],C] + [[B,C],A] + [[C,A],B] = 0

For the three transposition crystals A=(01), B=(02), C=(12):
  [A,B] = compose(A,B) - compose(B,A)
  [[A,B],C] = compose([A,B],C) - compose(C,[A,B])

If the Jacobi identity holds, the crystal composition algebra IS
a Lie algebra — not approximately, but structurally.

Also test: does [A,[B,C]] = [[A,B],C] + [B,[A,C]]? (Derivation property)
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_fidelity, EPS_LOG)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)


n_seeds = 200

# Build seed-averaged transposition crystals
corrs = {
    "A=(01)": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "B=(02)": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "C=(12)": torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
}

crystals = {}
for name, corr in corrs.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    crystals[name] = avg

A = crystals["A=(01)"]
B = crystals["B=(02)"]
C = crystals["C=(12)"]


def bracket(X, Y):
    """Crystal Lie bracket [X, Y] = compose(X,Y) - compose(Y,X)."""
    return compose(X, Y) - compose(Y, X)


def norm(X):
    """Frobenius norm."""
    return X.abs().pow(2).sum().sqrt().item()


print("=" * 80)
print("THE JACOBI IDENTITY")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Compute the three double brackets
# ═══════════════════════════════════════════════════════════════

AB = bracket(A, B)
BC = bracket(B, C)
CA = bracket(C, A)

AB_C = bracket(AB, C)
BC_A = bracket(BC, A)
CA_B = bracket(CA, B)

jacobi = AB_C + BC_A + CA_B

print(f"\n--- Individual terms ---\n")
print(f"  ||[A,B]|| = {norm(AB):.6f}")
print(f"  ||[B,C]|| = {norm(BC):.6f}")
print(f"  ||[C,A]|| = {norm(CA):.6f}")
print()
print(f"  ||[[A,B],C]|| = {norm(AB_C):.6f}")
print(f"  ||[[B,C],A]|| = {norm(BC_A):.6f}")
print(f"  ||[[C,A],B]|| = {norm(CA_B):.6f}")
print()
print(f"  ||[[A,B],C] + [[B,C],A] + [[C,A],B]|| = {norm(jacobi):.6f}")

# Relative Jacobi violation
max_double = max(norm(AB_C), norm(BC_A), norm(CA_B))
if max_double > 1e-10:
    relative = norm(jacobi) / max_double
    print(f"\n  Relative violation: {relative:.6f} ({relative*100:.2f}%)")
    print(f"  {'JACOBI HOLDS' if relative < 0.05 else 'JACOBI VIOLATED'} (threshold 5%)")
else:
    print(f"\n  All double brackets near zero.")


# ═══════════════════════════════════════════════════════════════
#  S₃ decomposition of the Jacobi sum
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("S₃ DECOMPOSITION OF JACOBI SUM")
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

for label, mat in [("[[A,B],C]", AB_C), ("[[B,C],A]", BC_A),
                    ("[[C,A],B]", CA_B), ("Jacobi sum", jacobi)]:
    f_t, f_s, f_st = s3_fracs(mat)
    print(f"  {label:>15s}: triv={f_t:.5f}, sign={f_s:.5f}, std={f_st:.5f}, norm={norm(mat):.6f}")


# ═══════════════════════════════════════════════════════════════
#  Does the Jacobi sum vanish in each sector separately?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("PER-SECTOR JACOBI VIOLATION")
print("="*80)

perms = [
    [0,1,2,3], [1,0,2,3], [2,1,0,3],
    [0,2,1,3], [1,2,0,3], [2,0,1,3],
]
signs_list = [1, -1, -1, -1, 1, 1]

def project_sign(M):
    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T
    return sum(s * permute(M, p) for p, s in zip(perms, signs_list)) / 6

def project_triv(M):
    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T
    return sum(permute(M, p) for p in perms) / 6

# Sign sector of each double bracket and Jacobi sum
sign_AB_C = project_sign(AB_C)
sign_BC_A = project_sign(BC_A)
sign_CA_B = project_sign(CA_B)
sign_jacobi = sign_AB_C + sign_BC_A + sign_CA_B

print(f"\n  Sign sector:")
print(f"    ||sign([[A,B],C])|| = {norm(sign_AB_C):.6f}")
print(f"    ||sign([[B,C],A])|| = {norm(sign_BC_A):.6f}")
print(f"    ||sign([[C,A],B])|| = {norm(sign_CA_B):.6f}")
print(f"    ||sign(Jacobi)|| = {norm(sign_jacobi):.6f}")
if norm(sign_AB_C) > 1e-10:
    print(f"    Relative sign violation: {norm(sign_jacobi)/max(norm(sign_AB_C), norm(sign_BC_A), norm(sign_CA_B)):.6f}")

# Trivial sector
triv_AB_C = project_triv(AB_C)
triv_BC_A = project_triv(BC_A)
triv_CA_B = project_triv(CA_B)
triv_jacobi = triv_AB_C + triv_BC_A + triv_CA_B

print(f"\n  Trivial sector:")
print(f"    ||triv([[A,B],C])|| = {norm(triv_AB_C):.6f}")
print(f"    ||triv(Jacobi)|| = {norm(triv_jacobi):.6f}")


# ═══════════════════════════════════════════════════════════════
#  Per-seed Jacobi violation
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("PER-SEED JACOBI VIOLATION")
print("="*80)

violations = []
for seed in range(min(100, n_seeds)):
    a = Entangler(corrs["A=(01)"], seed=seed).build().joint
    b = Entangler(corrs["B=(02)"], seed=seed).build().joint
    c = Entangler(corrs["C=(12)"], seed=seed).build().joint

    ab_c = bracket(bracket(a, b), c)
    bc_a = bracket(bracket(b, c), a)
    ca_b = bracket(bracket(c, a), b)

    j = ab_c + bc_a + ca_b
    j_norm = norm(j)
    max_norm = max(norm(ab_c), norm(bc_a), norm(ca_b))

    if max_norm > 1e-10:
        violations.append(j_norm / max_norm)

if violations:
    avg_v = sum(violations) / len(violations)
    std_v = (sum((v - avg_v)**2 for v in violations) / len(violations)) ** 0.5
    sorted_v = sorted(violations)

    print(f"\n  {len(violations)} seeds with measurable double brackets:")
    print(f"  Relative violation: {avg_v:.6f} ± {std_v:.6f}")
    print(f"  Median: {sorted_v[len(sorted_v)//2]:.6f}")
    print(f"  Min: {sorted_v[0]:.6f}, Max: {sorted_v[-1]:.6f}")
    print(f"  10%: {sorted_v[len(sorted_v)//10]:.6f}, 90%: {sorted_v[9*len(sorted_v)//10]:.6f}")


print(f"\n\n{'='*80}")
print("WHAT THE JACOBI IDENTITY REVEALS")
print("="*80)
