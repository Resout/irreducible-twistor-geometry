"""
Born(θ) = p*²/(H + p*²) at the composition fixed point.

At the product state with a = b, c = p*a, d = p*²a:
  Born(h_i) = (3a² + c²) / (9a² + 6c² + d²)  [each of 3 h's]
  Born(θ) = (3c² + d²) / (9a² + 6c² + d²)

Substituting c = pa, d = p²a:
  Born(θ) = a²(3p² + p⁴) / a²(9 + 6p² + p⁴)
           = p²(3 + p²) / (3 + p²)²
           = p² / (3 + p²)
           = p²/(H + p²)

This is EXACT (no approximation). Let me verify numerically and
explore its implications.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR, schmidt_number, born_probabilities
from solver.crystals import Entangler, compose, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)


S3_PERMS = [
    [0, 1, 2], [1, 0, 2], [2, 1, 0],
    [0, 2, 1], [1, 2, 0], [2, 0, 1],
]

def apply_perm(joint, perm):
    inv_perm = [0] * H
    for i in range(H):
        inv_perm[perm[i]] = i
    result = torch.zeros_like(joint)
    for i in range(MASS_DIM):
        for j in range(MASS_DIM):
            ii = inv_perm[i] if i < H else H
            jj = inv_perm[j] if j < H else H
            result[i, j] = joint[ii, jj]
    return result

def symmetrize(joint):
    result = torch.zeros_like(joint)
    for perm in S3_PERMS:
        result += apply_perm(joint, perm)
    return result / 6.0


# ═══════════════════════════════════════════════════════════════
#  Verify Born(θ) = p²/(H + p²) at fixed points
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("Born(θ) = p*²/(H + p*²) AT COMPOSITION FIXED POINT")
print("=" * 80)

print(f"\n  {'Type':>15s} {'p*':>8s} {'Born(θ) pred':>14s} {'Born(θ) meas':>14s} {'diff':>10s}")
print("  " + "-" * 65)

for name, corr in RELATIONSHIP_SIGNATURES.items():
    p_vals = []
    born_vals = []

    for seed in range(100):
        ent = Entangler(corr, seed=seed).build()
        M = symmetrize(ent.joint)
        current = M.clone()
        for _ in range(40):
            current = compose(current, M)

        a = current[0, 0]
        c = current[0, H]
        d = current[H, H]

        if a.abs().item() > 1e-10:
            p = (c / a).real.item()
            p_vals.append(p)

            # Measured Born(θ): use the actual Born probabilities
            marginal = current.abs().pow(2).sum(dim=1)
            bp = marginal / marginal.sum()
            born_vals.append(bp[H].item())

    p_avg = sum(p_vals) / len(p_vals)
    born_avg = sum(born_vals) / len(born_vals)

    # Predicted Born(θ) from p*
    born_pred = p_avg**2 / (H + p_avg**2)

    diff = abs(born_avg - born_pred)
    rel_diff = diff / born_avg * 100

    print(f"  {name:>15s} {p_avg:>8.4f} {born_pred:>14.6f} {born_avg:>14.6f} {rel_diff:>9.2f}%")


# ═══════════════════════════════════════════════════════════════
#  The inverse: p* = √(H·Born(θ)/(1 - Born(θ)))
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("INVERSE: p* = √(H·Born(θ)/(1 - Born(θ)))")
print("="*80)

print(f"\n  Special values:")
print(f"    Born(θ) = 1/H³ = {BORN_FLOOR:.6f} (floor):")
p_floor = math.sqrt(H * BORN_FLOOR / (1 - BORN_FLOOR))
print(f"      p* = {p_floor:.6f}")

print(f"    Born(θ) = 1/(H+1) = {1/(H+1):.6f} (uniform):")
p_unif = math.sqrt(H * (1/(H+1)) / (1 - 1/(H+1)))
print(f"      p* = {p_unif:.6f} (should be 1.0)")

print(f"    Born(θ) = 1/H = {1/H:.6f}:")
p_third = math.sqrt(H * (1/H) / (1 - 1/H))
print(f"      p* = {p_third:.6f} = √(H/(H-1)) = √(3/2)")

print(f"    Born(θ) = 1/2:")
p_half = math.sqrt(H * 0.5 / 0.5)
print(f"      p* = {p_half:.6f} = √H = √3")


# ═══════════════════════════════════════════════════════════════
#  Key identity: Born(h_i) = (3 + p²)/(3(3 + p²)) × 3 = 1/(3+p²) × 3?
#  Let me compute Born(h_i) carefully.
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("BORN DISTRIBUTION AT FIXED POINT")
print("="*80)

# At the product state with a=b, c=pa, d=p²a:
# M = [[a, a, a, pa],
#      [a, a, a, pa],
#      [a, a, a, pa],
#      [pa, pa, pa, p²a]]
#
# Born probabilities from |m|²:
# Total |m|² = 9|a|² + 6|pa|² + |p²a|² = |a|²(9 + 6p² + p⁴) = |a|²(3+p²)²
#
# Born(h_i) for i<H: sum of row i = |a|²(3 + p²) → Born = 3(3+p²)/(3+p²)² = 3/(3+p²)
# Wait... Born is from marginal, not row.
#
# Marginal: m_marginal[i] = Σ_j |m[i,j]|²
# For i < H: = 3|a|² + |pa|² = |a|²(3 + p²)
# For i = θ: = 3|pa|² + |p²a|² = |a|²p²(3 + p²)
#
# Total marginal = 3|a|²(3+p²) + |a|²p²(3+p²) = |a|²(3+p²)(3+p²) = |a|²(3+p²)²
#
# Born(h_i) = |a|²(3+p²) / [|a|²(3+p²)²] = 1/(3+p²) = 1/(H+p²)
# Born(θ) = |a|²p²(3+p²) / [|a|²(3+p²)²] = p²/(3+p²) = p²/(H+p²)
#
# Check: 3·Born(h) + Born(θ) = 3/(H+p²) + p²/(H+p²) = (3+p²)/(H+p²) = 1 ✓

print(f"\n  At product state (a=b, c=pa, d=p²a):")
print(f"    Born(h_i) = 1/(H + p²)   for each i = 0,1,2")
print(f"    Born(θ)   = p²/(H + p²)")
print(f"    Sum = H/(H+p²) + p²/(H+p²) = 1  ✓")

print(f"\n  Table of values:")
print(f"    {'p':>8s} {'Born(h_i)':>10s} {'Born(θ)':>10s} {'Born floor?':>12s}")
print("    " + "-" * 45)

for p in [0.0, 0.34, 0.5, 1.0, 1.10, 1.33, 1.73, 3.0]:
    bh = 1 / (H + p**2)
    bt = p**2 / (H + p**2)
    floor = "AT FLOOR" if abs(bt - BORN_FLOOR) < 0.001 else "above" if bt > BORN_FLOOR else "BELOW"
    print(f"    {p:>8.3f} {bh:>10.6f} {bt:>10.6f} {floor:>12s}")


# ═══════════════════════════════════════════════════════════════
#  The Born(θ) of the ORIGINAL (non-symmetrized) crystal
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("BORN(θ) TRAJECTORY: original → composition → fixed point")
print("="*80)

for name in ["proportional", "modular"]:
    corr = RELATIONSHIP_SIGNATURES[name]
    ent = Entangler(corr, seed=42).build()
    M = symmetrize(ent.joint)

    print(f"\n  {name}:")
    current = M.clone()
    for step in range(15):
        marginal = current.abs().pow(2).sum(dim=1)
        bp = marginal / marginal.sum()
        born_theta = bp[H].item()

        a = current[0, 0]
        c = current[0, H]
        p = (c/a).real.item() if a.abs().item() > 1e-10 else 0
        pred = p**2 / (H + p**2) if (H + p**2) > 0 else 0

        print(f"    Step {step:>2d}: Born(θ) = {born_theta:.6f}, p²/(H+p²) = {pred:.6f}, diff = {abs(born_theta-pred):.6f}")
        current = compose(current, M)


# ═══════════════════════════════════════════════════════════════
#  SCHMIDT at the product state fixed point
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SCHMIDT AT THE PRODUCT STATE")
print("="*80)

# Schmidt = 1/Σ Born(i)² = 1/(H·Born(h)² + Born(θ)²)
# = 1/(H/(H+p²)² + p⁴/(H+p²)²)
# = (H+p²)²/(H + p⁴)

print(f"  Schmidt = (H + p²)² / (H + p⁴)")
print(f"\n  {'p':>8s} {'Schmidt':>10s}")
print("  " + "-" * 22)

for p in [0.0, 0.5, 1.0, 1.10, 1.33, 1.73, 3.0, 10.0]:
    sc = (H + p**2)**2 / (H + p**4)
    print(f"  {p:>8.3f} {sc:>10.4f}")

print(f"\n  At p=0: Schmidt = H²/H = H = {H}")
print(f"  At p=1: Schmidt = (H+1)²/(H+1) = H+1 = {H+1}")
print(f"  At p→∞: Schmidt → p⁴/p⁴ = 1")
print(f"  Maximum at p = ? (differentiate and set to 0)")

# d/dp[(H+p²)²/(H+p⁴)] = 0
# Let u = p². d/du[(H+u)²/(H+u²)] = 0
# [2(H+u)(H+u²) - (H+u)²·2u] / (H+u²)² = 0
# 2(H+u)[(H+u²) - u(H+u)] = 0
# (H+u²) - u(H+u) = 0
# H + u² - uH - u² = 0
# H(1 - u) = 0
# u = 1 → p = 1!

print(f"  Maximum at p = 1: Schmidt = {(H+1)**2/(H+1)} = H+1 = {H+1}")
print(f"  The product state with p=1 (uniform Born) has maximum Schmidt!")

# Fixed point Schmidts
print(f"\n  Fixed point Schmidts:")
for name_p, p_val in [("proportional", 1.10), ("modular", 1.33)]:
    sc = (H + p_val**2)**2 / (H + p_val**4)
    print(f"    {name_p}: p*={p_val:.2f}, Schmidt = {sc:.4f}")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
print(f"""
  EXACT RESULTS at the composition fixed point (product state):

  Born(h_i) = 1/(H + p*²)
  Born(θ)   = p*²/(H + p*²)
  Schmidt   = (H + p*²)²/(H + p*⁴)

  These are EXACT (analytical derivation from the product state form).

  The parameter p* = c*/a* encodes:
    p* = 0:  pure knowledge (Born(θ) = 0, impossible by floor)
    p* = 1:  uniform distribution (Born = 1/(H+1) for all)
    p* > 1:  ignorance-dominated (all actual fixed points)
    p* → ∞: pure ignorance (Born(θ) → 1, Schmidt → 1)

  All actual composition fixed points have p* > 1:
    proportional: p* ≈ 1.10, Born(θ) ≈ 0.29
    modular:      p* ≈ 1.33, Born(θ) ≈ 0.37

  The Born floor constrains p* ≥ √(H·BORN/(1-BORN)) = √(3/26) ≈ 0.34.
  All actual p* are well above this minimum.

  Schmidt is maximized at p* = 1 (uniform) where Schmidt = H+1 = 4.
  All product states have Schmidt < H+1. This is the information-
  theoretic maximum for a rank-1 state.
""")
