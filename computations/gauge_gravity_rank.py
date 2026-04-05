"""
Rank of gauge (commutator) vs gravity (anticommutator).

Gauge is rank 2 = dim(root space). What is gravity?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import torch
from solver.algebra import H, MASS_DIM
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)

n_seeds = 300
corrs = [
    ("(01)", torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32)),
    ("(02)", torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)),
]
crystals = {}
for name, corr in corrs:
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    crystals[name] = avg

A, B = crystals["(01)"], crystals["(02)"]

AB = compose(A, B)
BA = compose(B, A)
comm = AB - BA
anti = AB + BA

def raw(X, Y):
    return torch.einsum("ab,bc->ac", X, Y)

print("=" * 60)
print("RANK: GAUGE vs GRAVITY")
print("=" * 60)

for label, M in [("Crystal A", A), ("Product AB", AB),
                  ("Commutator [A,B]", comm), ("Anticommutator {A,B}", anti)]:
    _, S, _ = torch.linalg.svd(M)
    S_abs = S.abs()
    S_norm = S_abs / S_abs.sum()
    eff_rank = 1 / (S_norm**2).sum().item()
    print(f"\n  {label}:")
    for i in range(MASS_DIM):
        bar = "#" * int(S_norm[i].item() * 50)
        print(f"    σ{i} = {S_abs[i].item():.6f} ({S_norm[i].item():.4f}) {bar}")
    print(f"    Effective rank: {eff_rank:.3f}")

print("\n\n--- Raw (un-normalized) ---")
for label, M in [("Raw AB", raw(A,B)), ("Raw [A,B]", raw(A,B)-raw(B,A)),
                  ("Raw {A,B}", raw(A,B)+raw(B,A))]:
    _, S, _ = torch.linalg.svd(M)
    S_abs = S.abs()
    S_norm = S_abs / S_abs.sum()
    eff_rank = 1 / (S_norm**2).sum().item()
    print(f"\n  {label}:")
    for i in range(MASS_DIM):
        print(f"    σ{i} = {S_abs[i].item():.6f} ({S_norm[i].item():.4f})")
    print(f"    Effective rank: {eff_rank:.3f}")

print("\n\nSUMMARY:")
for label, M in [("Gauge [A,B]", comm), ("Gravity {A,B}", anti)]:
    _, S, _ = torch.linalg.svd(M)
    S_norm = S.abs() / S.abs().sum()
    eff_rank = 1 / (S_norm**2).sum().item()
    print(f"  {label}: rank {eff_rank:.1f}")

print(f"\n  Root space dim = 2 (A₂)")
print(f"  Mass dim = {MASS_DIM}")
print(f"  Joint dim = {MASS_DIM**2}")
