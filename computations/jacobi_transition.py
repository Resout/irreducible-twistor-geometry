"""
The Jacobi phase transition: from Lie algebra to mass gap.

Interpolate the normalization strength α from 0 to 1:
  compose_α(A,B) = raw(A,B) × (1-α) + normalized(A,B) × α

Or more physically: partial L1 normalization
  compose_α(A,B) = raw(A,B) / L1^α   where L1 = ||raw(A,B)||_L1

At α=0: raw (Jacobi holds, free gauge field)
At α=1: full normalization (Jacobi fails, confined gauge field)

The transition point α_crit should be where the Jacobi violation
first becomes significant. If α_crit = K* = 7/30, the mass gap
IS the critical normalization strength.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, EPS_LOG)
from solver.crystals import Entangler

torch.set_grad_enabled(False)


def compose_alpha(A, B, alpha):
    """Composition with partial normalization.

    alpha=0: raw tensor contraction
    alpha=1: full L1 normalization
    """
    raw = torch.einsum("ab,bc->ac", A, B)
    L1 = raw.reshape(-1).real.sum()
    if abs(L1) < 1e-15:
        return raw
    # Partial normalization: divide by L1^alpha
    return raw / (abs(L1) ** alpha)


def bracket_alpha(A, B, alpha):
    return compose_alpha(A, B, alpha) - compose_alpha(B, A, alpha)


def norm(X):
    return X.abs().pow(2).sum().sqrt().item()


n_seeds = 200

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


print("=" * 80)
print("THE JACOBI PHASE TRANSITION")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Jacobi violation as a function of α
# ═══════════════════════════════════════════════════════════════

alphas = [0, 0.05, 0.1, 0.15, 0.2, 0.233, 0.25, 0.3, 0.35, 0.4,
          0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

print(f"\n  {'α':>6s}  {'||Jacobi||':>12s}  {'max||[[,],]||':>14s}  {'relative':>10s}  {'note':>10s}")
print("-" * 60)

for alpha in alphas:
    AB = bracket_alpha(A, B, alpha)
    BC = bracket_alpha(B, C, alpha)
    CA = bracket_alpha(C, A, alpha)

    AB_C = bracket_alpha(AB, C, alpha)
    BC_A = bracket_alpha(BC, A, alpha)
    CA_B = bracket_alpha(CA, B, alpha)

    jacobi = AB_C + BC_A + CA_B
    j_norm = norm(jacobi)
    max_db = max(norm(AB_C), norm(BC_A), norm(CA_B))
    relative = j_norm / max_db if max_db > 1e-15 else 0

    note = ""
    if abs(alpha - K_STAR) < 0.01:
        note = "≈ K*"
    elif abs(alpha - DELTA) < 0.01:
        note = "≈ Δ"

    print(f"  {alpha:>6.3f}  {j_norm:>12.4e}  {max_db:>14.4e}  {relative:>10.6f}  {note:>10s}")


# ═══════════════════════════════════════════════════════════════
#  Fine-grained scan around the transition
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("FINE-GRAINED SCAN: WHERE DOES RELATIVE VIOLATION REACH 1%?")
print("="*80)

fine_alphas = [i * 0.01 for i in range(0, 101)]
transitions = []

for alpha in fine_alphas:
    AB = bracket_alpha(A, B, alpha)
    BC = bracket_alpha(B, C, alpha)
    CA = bracket_alpha(C, A, alpha)

    AB_C = bracket_alpha(AB, C, alpha)
    BC_A = bracket_alpha(BC, A, alpha)
    CA_B = bracket_alpha(CA, B, alpha)

    jacobi = AB_C + BC_A + CA_B
    j_norm = norm(jacobi)
    max_db = max(norm(AB_C), norm(BC_A), norm(CA_B))
    relative = j_norm / max_db if max_db > 1e-15 else 0

    transitions.append((alpha, relative, j_norm))

# Find where relative first exceeds various thresholds
for threshold in [0.001, 0.01, 0.05, 0.1, 0.5]:
    for alpha, rel, jn in transitions:
        if rel > threshold:
            print(f"  Relative > {threshold:.3f} at α = {alpha:.2f}")
            break


# ═══════════════════════════════════════════════════════════════
#  The transition shape
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("TRANSITION SHAPE")
print("="*80)

print(f"\n  Reference values: K* = {K_STAR:.4f}, Δ = {DELTA:.4f}, 1/(H²+1) = {1/(H**2+1):.4f}")
print()
print(f"  {'α':>6s}  {'relative':>10s}")
for alpha, rel, jn in transitions:
    if alpha <= 0.50 and int(alpha * 100) % 2 == 0:
        bar = '█' * min(int(rel * 50), 50)
        marker = " ← K*" if abs(alpha - K_STAR) < 0.005 else ""
        marker = " ← Δ" if abs(alpha - DELTA) < 0.005 else marker
        marker = " ← 1/(H²+1)" if abs(alpha - 0.10) < 0.005 else marker
        print(f"  {alpha:>6.2f}  {rel:>10.6f}  {bar}{marker}")


# ═══════════════════════════════════════════════════════════════
#  What's at the critical point?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE CRITICAL NORMALIZATION STRENGTH")
print("="*80)

# Find the α where relative violation = 50% (midpoint of transition)
for i in range(len(transitions) - 1):
    a1, r1, _ = transitions[i]
    a2, r2, _ = transitions[i+1]
    if r1 < 0.5 and r2 >= 0.5:
        # Linear interpolation
        alpha_crit = a1 + (0.5 - r1) * (a2 - a1) / (r2 - r1)
        print(f"\n  α_crit (50% relative violation) ≈ {alpha_crit:.4f}")
        print(f"\n  Candidates:")
        print(f"    K* = {K_STAR:.4f}  diff = {alpha_crit - K_STAR:+.4f}")
        print(f"    Δ = {DELTA:.4f}  diff = {alpha_crit - DELTA:+.4f}")
        print(f"    1/MASS_DIM = {1/MASS_DIM:.4f}  diff = {alpha_crit - 1/MASS_DIM:+.4f}")
        print(f"    1/H = {1/H:.4f}  diff = {alpha_crit - 1/H:+.4f}")
        print(f"    BORN_FLOOR = {BORN_FLOOR:.4f}  diff = {alpha_crit - BORN_FLOOR:+.4f}")
        break


print(f"\n\n{'='*80}")
print("WHAT THE JACOBI TRANSITION REVEALS")
print("="*80)
