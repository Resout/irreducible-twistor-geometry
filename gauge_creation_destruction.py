"""
Gauge field lifecycle: creation in the entangler, destruction in composition.

The entangler creates standard sector content (gauge-variant structure).
Composition destroys it at rate 2Δ. At the K* equilibrium crystal,
f_std ≈ 0.20. Under composition, the fixed point is f_std ≈ 0.02.

The physical gauge field content should be the balance between these
two processes. If we compose the equilibrium crystal once and measure
the standard sector change, we get the destruction rate PER COMPOSITION
STEP. If we entangle one more step and measure the creation rate,
we get the creation rate PER ENTANGLING STEP.

The ratio creation/destruction should equal the entangler/composition
step ratio (50 build : 1 compose? Or something else?).

Also: the sign sector differentiates identity from cycle. Does it
survive composition?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import random as _random
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, EPS_LOG, EPS_DIV, EPS_NORM)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)


def s3_fractions(joint):
    """S₃ irrep fractions."""
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]
    signs = [1, -1, -1, -1, 1, 1]

    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M.real.float() @ P.T

    triv = sum(permute(joint, p) for p in perms) / 6
    sign = sum(s * permute(joint, p) for p, s in zip(perms, signs)) / 6

    bp = joint.abs().pow(2)
    total = bp.sum().item()
    if total < 1e-15:
        return 1.0, 0.0, 0.0

    f_triv = triv.to(joint.dtype).abs().pow(2).sum().item() / total
    f_sign = sign.to(joint.dtype).abs().pow(2).sum().item() / total
    f_std = max(0, 1 - f_triv - f_sign)
    return f_triv, f_sign, f_std


# Build crystals at the K* equilibrium
identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
cycle_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)
inverse_corr = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)

corrs = {
    "identity": identity_corr,
    "inverse": inverse_corr,
    "cycle (012)": cycle_corr,
}

print("=" * 80)
print("GAUGE CREATION AND DESTRUCTION")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Starting crystals and their irrep content
# ═══════════════════════════════════════════════════════════════

n_seeds = 50
crystals = {}

for name, corr in corrs.items():
    avg_joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        ent = Entangler(corr, seed=seed).build(n_steps=50, discount=0.3)
        avg_joint += ent.joint
    avg_joint /= n_seeds
    crystals[name] = avg_joint

print(f"\n--- Starting crystals (50-seed average) ---\n")
print(f"  {'name':>15s}  {'trivial':>8s}  {'sign':>8s}  {'standard':>8s}  {'Schmidt':>8s}")
for name, joint in crystals.items():
    f_t, f_s, f_st = s3_fractions(joint)
    sn = schmidt_number(joint)
    print(f"  {name:>15s}  {f_t:>8.4f}  {f_s:>8.4f}  {f_st:>8.4f}  {sn:>8.3f}")


# ═══════════════════════════════════════════════════════════════
#  Self-composition cascade: how does the standard sector decay?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("STANDARD SECTOR UNDER SELF-COMPOSITION")
print("="*80)

n_comp_steps = 15

for name in corrs:
    print(f"\n  {name}:")
    print(f"  {'comp':>4s}  {'trivial':>8s}  {'sign':>8s}  {'standard':>8s}  {'Schmidt':>8s}  {'ln(f_std)':>10s}")

    current = crystals[name].clone()
    original = crystals[name].clone()

    for comp_step in range(n_comp_steps):
        f_t, f_s, f_st = s3_fractions(current)
        sn = schmidt_number(current)
        log_std = math.log(f_st) if f_st > 1e-10 else float('-inf')

        print(f"  {comp_step:4d}  {f_t:>8.4f}  {f_s:>8.4f}  {f_st:>8.4f}  {sn:>8.3f}  {log_std:>10.4f}")

        current = compose(current, original)


# ═══════════════════════════════════════════════════════════════
#  Standard sector decay rate under composition
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("STANDARD SECTOR DECAY RATE UNDER COMPOSITION")
print("="*80)

for name in corrs:
    current = crystals[name].clone()
    original = crystals[name].clone()
    prev_std = s3_fractions(current)[2]

    rates = []
    for comp_step in range(1, 8):
        current = compose(current, original)
        f_t, f_s, f_st = s3_fractions(current)
        if prev_std > 1e-10 and f_st > 1e-10:
            rate = -math.log(f_st / prev_std)
            rates.append(rate)
        prev_std = f_st

    if rates:
        avg_rate = sum(rates) / len(rates)
        print(f"\n  {name}: avg decay rate = {avg_rate:.4f} per composition")
        print(f"    = {avg_rate/DELTA:.3f}Δ  (compare: overall gap = 2Δ)")
        print(f"    Rates per step: {', '.join(f'{r:.4f}' for r in rates)}")


# ═══════════════════════════════════════════════════════════════
#  Cross-composition: does sign survive?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SIGN SECTOR UNDER COMPOSITION")
print("="*80)

# Compose identity × cycle, cycle × identity, cycle × cycle
combos = [
    ("id ∘ id", crystals["identity"], crystals["identity"]),
    ("id ∘ cyc", crystals["identity"], crystals["cycle (012)"]),
    ("cyc ∘ id", crystals["cycle (012)"], crystals["identity"]),
    ("cyc ∘ cyc", crystals["cycle (012)"], crystals["cycle (012)"]),
    ("inv ∘ cyc", crystals["inverse"], crystals["cycle (012)"]),
]

print(f"\n  {'product':>10s}  {'trivial':>8s}  {'sign':>8s}  {'standard':>8s}  {'Schmidt':>8s}")
for label, a, b in combos:
    result = compose(a, b)
    f_t, f_s, f_st = s3_fractions(result)
    sn = schmidt_number(result)
    print(f"  {label:>10s}  {f_t:>8.4f}  {f_s:>8.4f}  {f_st:>8.4f}  {sn:>8.3f}")


# ═══════════════════════════════════════════════════════════════
#  The sign/standard ratio
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SIGN-TO-STANDARD RATIO")
print("="*80)

print(f"\n  For the cycle crystal:")
for comp_step in range(8):
    if comp_step == 0:
        current = crystals["cycle (012)"].clone()
    else:
        current = compose(current, crystals["cycle (012)"])
    f_t, f_s, f_st = s3_fractions(current)
    ratio = f_s / f_st if f_st > 1e-10 else float('inf')
    print(f"    comp {comp_step}: sign={f_s:.5f}, std={f_st:.5f}, sign/std={ratio:.4f}")


print(f"\n\n{'='*80}")
print("WHAT GAUGE CREATION AND DESTRUCTION SHOW")
print("="*80)
