"""
The sign sector oscillation period = group element order.

For cycles (order 3): sign vanishes at steps 2, 5, 8, ... (every 3rd)
For transpositions (order 2): sign should vanish at steps 1, 3, 5, ... (every 2nd)

Because σⁿ = e at multiples of ord(σ), and the identity has zero sign.
The envelope decays from the structural filter (rate ~3Δ).

The GROUP ALGEBRA is alive inside the composition cascade.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, EPS_LOG)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)


def s3_fractions(joint):
    """Correct complex S₃ projector."""
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


identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
inverse_corr = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)  # transposition (02)
cycle_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)    # 3-cycle (012)
swap01_corr = torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32)   # transposition (01)

corrs = {
    "identity (e)": identity_corr,
    "trans (02)": inverse_corr,
    "trans (01)": swap01_corr,
    "cycle (012)": cycle_corr,
}

print("=" * 80)
print("SIGN OSCILLATION PERIOD = GROUP ELEMENT ORDER")
print("=" * 80)

n_seeds = 200
n_comp = 12

for name, corr in corrs.items():
    # Build seed-averaged crystal
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds

    print(f"\n  {name} ({n_seeds}-seed avg):")
    print(f"  {'comp':>4s}  {'f_sign':>12s}  {'f_std':>12s}  {'f_triv':>12s}  {'note':>15s}")

    current = avg.clone()
    for comp in range(n_comp):
        f_t, f_s, f_st = s3_fractions(current)

        # Determine algebraic product
        # For cycle (order 3): comp%3==2 → e
        # For transposition (order 2): comp%2==1 → e
        note = ""
        if "cycle" in name:
            if comp % 3 == 2:
                note = "≡ e (mod 3)"
            elif comp % 3 == 0:
                note = f"≡ σ"
            else:
                note = f"≡ σ²"
        elif "trans" in name:
            if comp % 2 == 1:
                note = "≡ e (mod 2)"
            elif comp % 2 == 0:
                note = f"≡ σ"
        elif "identity" in name:
            note = "always e"

        print(f"  {comp:4d}  {f_s:>12.8f}  {f_st:>12.8f}  {f_t:>12.8f}  {note:>15s}")

        current = compose(current, avg)


# ═══════════════════════════════════════════════════════════════
#  Envelope decay rate (on the non-zero steps)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("ENVELOPE DECAY RATE OF SIGN SECTOR")
print("="*80)

# For cycle: non-zero at steps 0, 1, 3, 4, 6, 7, ...
# Compare step 0 vs 3, 1 vs 4, etc. (same position in period)

print(f"\n  Cycle (012):")
avg_cyc = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for seed in range(n_seeds):
    avg_cyc += Entangler(cycle_corr, seed=seed).build().joint
avg_cyc /= n_seeds

sign_vals = []
current = avg_cyc.clone()
for comp in range(n_comp):
    f_t, f_s, f_st = s3_fractions(current)
    sign_vals.append(f_s)
    current = compose(current, avg_cyc)

# Period-3: compare step n with step n+3
print(f"  {'step':>4s}  {'f_sign':>12s}  {'step+3':>6s}  {'f_sign(+3)':>12s}  {'ratio':>10s}  {'rate_per_3':>10s}")
for n in range(min(7, n_comp - 3)):
    if sign_vals[n] > 1e-8 and sign_vals[n+3] > 1e-8:
        ratio = sign_vals[n+3] / sign_vals[n]
        rate = -math.log(ratio) / 3
        print(f"  {n:4d}  {sign_vals[n]:>12.8f}  {n+3:6d}  {sign_vals[n+3]:>12.8f}  {ratio:>10.6f}  {rate:>10.5f}")


# For transposition: non-zero at steps 0, 2, 4, ...
print(f"\n  Transposition (02):")
avg_trans = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for seed in range(n_seeds):
    avg_trans += Entangler(inverse_corr, seed=seed).build().joint
avg_trans /= n_seeds

sign_trans = []
current = avg_trans.clone()
for comp in range(n_comp):
    f_t, f_s, f_st = s3_fractions(current)
    sign_trans.append(f_s)
    current = compose(current, avg_trans)

# Period-2: compare step n with step n+2
print(f"  {'step':>4s}  {'f_sign':>12s}  {'step+2':>6s}  {'f_sign(+2)':>12s}  {'ratio':>10s}  {'rate_per_2':>10s}")
for n in range(min(8, n_comp - 2)):
    if sign_trans[n] > 1e-8 and sign_trans[n+2] > 1e-8:
        ratio = sign_trans[n+2] / sign_trans[n]
        rate = -math.log(ratio) / 2
        print(f"  {n:4d}  {sign_trans[n]:>12.8f}  {n+2:6d}  {sign_trans[n+2]:>12.8f}  {ratio:>10.6f}  {rate:>10.5f}")


# ═══════════════════════════════════════════════════════════════
#  The standard sector: does it also oscillate?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("STANDARD SECTOR OSCILLATION?")
print("="*80)

print(f"\n  Does f_std oscillate with the group order too?")
print(f"\n  Cycle (012), f_std:")
current = avg_cyc.clone()
for comp in range(n_comp):
    f_t, f_s, f_st = s3_fractions(current)
    note = "≡ e" if comp % 3 == 2 else f"≡ σ{'²' if comp%3==1 else ''}"
    print(f"    comp {comp}: f_std = {f_st:.8f}  ({note})")
    current = compose(current, avg_cyc)


# ═══════════════════════════════════════════════════════════════
#  What about the FIDELITY to the expected group product?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("FIDELITY TO EXPECTED GROUP PRODUCT")
print("="*80)

from solver.algebra import born_fidelity

# Build reference crystals
refs = {}
for name, corr in corrs.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    refs[name] = avg

# Cycle composition should cycle through: (012), (021), e, (012), ...
# (012)^1 = (012), (012)^2 = (021) ≈ inverse of cycle
# Actually (012)^2 = (021). At H=3, (021) maps 0→2→1→0,
# which is [[0,0,1],[1,0,0],[0,1,0]] — let me build this

cycle021_corr = torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32)
avg_021 = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for seed in range(n_seeds):
    avg_021 += Entangler(cycle021_corr, seed=seed).build().joint
avg_021 /= n_seeds

expected_products = [
    ("(012)", refs["cycle (012)"]),
    ("(021)", avg_021),
    ("e", refs["identity (e)"]),
]

print(f"\n  Composition of cycle (012) with itself:")
print(f"  {'comp':>4s}  {'expected':>8s}  {'fid→(012)':>10s}  {'fid→(021)':>10s}  {'fid→e':>10s}  {'best':>8s}")

current = avg_cyc.clone()
for comp in range(9):
    expected_idx = comp % 3
    fids = []
    for label, ref in expected_products:
        fid = born_fidelity(current, ref)
        fids.append((label, fid))

    best = max(fids, key=lambda x: x[1])
    exp_label = expected_products[expected_idx][0]

    print(f"  {comp:4d}  {exp_label:>8s}  {fids[0][1]:>10.5f}  {fids[1][1]:>10.5f}  {fids[2][1]:>10.5f}  {best[0]:>8s}")

    current = compose(current, avg_cyc)


print(f"\n\n{'='*80}")
print("WHAT THE GROUP ORDER OSCILLATION REVEALS")
print("="*80)
