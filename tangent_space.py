"""
The linearized Lie algebra: tangent space at the identity crystal.

The identity crystal e has joint mass concentrated on the diagonal.
The tangent space T_e(Crystal) should be the space of infinitesimal
deformations e + εX that remain valid crystals.

Key questions:
1. What is the dimension of the tangent space?
2. Does it close under the bracket [X,Y] = XY - YX?
3. What is its structure (gl(4,C)? sl(4,C)? su(3)? a₂?)
4. How do the S₃ irrep sectors decompose the tangent space?

Strategy:
  - Compute (crystal - identity) / ε for each S₃ element
  - These are the "directions" in crystal space
  - Check if their brackets form a closed algebra
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, EPS_LOG,
                            schmidt_number, born_fidelity)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)


n_seeds = 500

# Build all S₃ element crystals
corrs = {
    "e": torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "(01)": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "(02)": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),
    "(12)": torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

crystals = {}
for name, corr in corrs.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    crystals[name] = avg

E = crystals["e"]


# ═══════════════════════════════════════════════════════════════
#  Tangent vectors: X_g = crystal_g - identity
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("TANGENT VECTORS AT IDENTITY")
print("=" * 80)

tangents = {}
for name in ["(01)", "(02)", "(12)", "(012)", "(021)"]:
    X = crystals[name] - E
    tangents[name] = X
    print(f"\n  X_{name} = crystal({name}) - e")
    print(f"    ||X|| = {X.abs().pow(2).sum().sqrt().item():.8f}")
    print(f"    real L1 = {X.reshape(-1).real.sum().item():.8f}")
    print(f"    imag L1 = {X.reshape(-1).imag.sum().item():.8f}")

    # SVD of tangent vector
    _, S, _ = torch.linalg.svd(X)
    S_abs = S.abs()
    S_norm = S_abs / S_abs.sum().clamp(min=1e-10)
    eff_rank = 1 / (S_norm**2).sum().item() if (S_norm**2).sum().item() > 1e-10 else 0
    print(f"    effective rank: {eff_rank:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Gram matrix: inner products of tangent vectors
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("GRAM MATRIX ⟨X_g, X_h⟩ (Frobenius inner product)")
print("="*80)

names = ["(01)", "(02)", "(12)", "(012)", "(021)"]
gram = torch.zeros(5, 5)

print(f"\n  {'':>6s}", end="")
for n in names:
    print(f"  {n:>8s}", end="")
print()

for i, ni in enumerate(names):
    print(f"  {ni:>6s}", end="")
    for j, nj in enumerate(names):
        ip = torch.sum(tangents[ni].conj() * tangents[nj]).real.item()
        gram[i, j] = ip
        print(f"  {ip:>8.5f}", end="")
    print()

# Eigenvalues of Gram matrix
eigvals = torch.linalg.eigvalsh(gram)
print(f"\n  Gram eigenvalues: {eigvals.tolist()}")
print(f"  Effective dimension: {(eigvals > 1e-6).sum().item()}")


# ═══════════════════════════════════════════════════════════════
#  Do tangent vectors close under raw bracket?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CLOSURE UNDER RAW BRACKET [X,Y] = XY - YX")
print("="*80)

def raw_bracket(X, Y):
    return torch.einsum("ab,bc->ac", X, Y) - torch.einsum("ab,bc->ac", Y, X)


# Bracket of tangent vectors
for ni, nj in [("(01)","(02)"), ("(01)","(12)"), ("(02)","(12)"),
                ("(01)","(012)"), ("(012)","(021)")]:
    Xi, Xj = tangents[ni], tangents[nj]
    bracket = raw_bracket(Xi, Xj)
    br_norm = bracket.abs().pow(2).sum().sqrt().item()

    # Project bracket onto tangent space
    # Coefficients c_k = ⟨bracket, X_k⟩ / ⟨X_k, X_k⟩
    # Better: use Gram matrix
    projections = []
    for k, nk in enumerate(names):
        coeff = torch.sum(bracket.conj() * tangents[nk]).real.item() / \
                torch.sum(tangents[nk].conj() * tangents[nk]).real.item()
        projections.append((nk, coeff))

    # Residual
    reconstructed = sum(c * tangents[n] for n, c in projections)
    residual = bracket - reconstructed
    res_norm = residual.abs().pow(2).sum().sqrt().item()

    print(f"\n  [{ni}, {nj}]: ||bracket|| = {br_norm:.6f}, ||residual|| = {res_norm:.6f}")
    print(f"    residual/bracket = {res_norm/br_norm if br_norm > 1e-10 else 0:.6f}")
    for nk, c in projections:
        if abs(c) > 0.001:
            print(f"    coefficient of X_{nk} = {c:.6f}")


# ═══════════════════════════════════════════════════════════════
#  S₃ decomposition of tangent vectors
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("S₃ DECOMPOSITION OF TANGENT SPACE")
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

for name in names:
    X = tangents[name]
    triv, sign, std = s3_project(X)
    nt = triv.abs().pow(2).sum().item()
    ns = sign.abs().pow(2).sum().item()
    nst = std.abs().pow(2).sum().item()
    total = nt + ns + nst
    print(f"\n  X_{name}: triv={100*nt/total:.1f}%, sign={100*ns/total:.1f}%, std={100*nst/total:.1f}%")


# ═══════════════════════════════════════════════════════════════
#  The identity crystal itself
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE IDENTITY CRYSTAL")
print("="*80)

print(f"\n  e (identity crystal, {MASS_DIM}×{MASS_DIM}):")
for i in range(MASS_DIM):
    row = ""
    for j in range(MASS_DIM):
        val = E[i, j]
        row += f"  {val.real.item():>8.5f}"
        if abs(val.imag.item()) > 1e-5:
            row += f"+{val.imag.item():.3f}i"
    print(f"    [{row} ]")

_, S_e, _ = torch.linalg.svd(E)
print(f"\n  SVD of e: {S_e.abs().tolist()}")
print(f"  Schmidt(e) = {schmidt_number(E):.4f}")

triv_e, sign_e, std_e = s3_project(E)
nt = triv_e.abs().pow(2).sum().item()
ns = sign_e.abs().pow(2).sum().item()
nst = std_e.abs().pow(2).sum().item()
total = nt + ns + nst
print(f"  S₃ decomposition: triv={100*nt/total:.1f}%, sign={100*ns/total:.1f}%, std={100*nst/total:.1f}%")


# ═══════════════════════════════════════════════════════════════
#  Composition acts like matrix product near identity
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("LINEARIZATION: compose(e+εX, e+εY) ≈ e + ε(X+Y) + ε²XY + O(ε³)")
print("="*80)

# Check: compose(e, crystal_g) ≈ crystal_g ?
for name in ["(01)", "(02)", "(012)"]:
    comp = compose(E, crystals[name])
    diff = (comp - crystals[name]).abs().pow(2).sum().sqrt().item()
    print(f"\n  compose(e, crystal({name})) - crystal({name}):")
    print(f"    ||diff|| = {diff:.8f}")
    fid = born_fidelity(comp, crystals[name])
    print(f"    fidelity = {fid:.8f}")

# Check: compose(crystal_g, e) ≈ crystal_g ?
for name in ["(01)", "(02)", "(012)"]:
    comp = compose(crystals[name], E)
    diff = (comp - crystals[name]).abs().pow(2).sum().sqrt().item()
    print(f"\n  compose(crystal({name}), e) - crystal({name}):")
    print(f"    ||diff|| = {diff:.8f}")
    fid = born_fidelity(comp, crystals[name])
    print(f"    fidelity = {fid:.8f}")


print(f"\n\n{'='*80}")
print("WHAT THE TANGENT SPACE REVEALS")
print("="*80)
