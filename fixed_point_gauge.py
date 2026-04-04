"""
Is f_std at the composition fixed point exactly 1/(H³+1) = 1/28?

The corrected projector gives f_std ≈ 3.57% for the seed-averaged
identity crystal after 15 compositions. Need higher precision:
- More seeds (200)
- More compositions (50)
- Check convergence to see if it's still moving

Also: the sign sector decay rate in the cycle crystal.
0.263 → 0.127 → 0.000. Ratio ≈ 0.48 ≈ 1/2.
Is the sign decay rate exactly ln(2) per composition step?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities, EPS_LOG)
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

    f_triv = triv.abs().pow(2).sum().item() / total
    f_sign = sign.abs().pow(2).sum().item() / total
    f_std = std.abs().pow(2).sum().item() / total
    return f_triv, f_sign, f_std


identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
cycle_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)
inverse_corr = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)

print("=" * 80)
print("FIXED POINT GAUGE CONTENT: HIGH PRECISION")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  High-precision f_std at composition fixed point
# ═══════════════════════════════════════════════════════════════

n_seeds = 200
n_compositions = 40

corrs = {
    "identity": identity_corr,
    "inverse": inverse_corr,
    "cycle": cycle_corr,
}

for name, corr in corrs.items():
    # Build seed-averaged crystal
    avg_crystal = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        ent = Entangler(corr, seed=seed).build()
        avg_crystal += ent.joint
    avg_crystal /= n_seeds

    # Track convergence through composition
    current = avg_crystal.clone()
    print(f"\n  {name} ({n_seeds}-seed average):")
    print(f"  {'comp':>4s}  {'f_std':>10s}  {'f_sign':>10s}  {'f_triv':>10s}  {'Schmidt':>8s}")

    for comp in range(n_compositions):
        f_t, f_s, f_st = s3_fractions(current)
        sn = schmidt_number(current)
        if comp < 10 or comp % 5 == 0:
            print(f"  {comp:4d}  {f_st:>10.7f}  {f_s:>10.7f}  {f_t:>10.7f}  {sn:>8.4f}")
        current = compose(current, avg_crystal)


# ═══════════════════════════════════════════════════════════════
#  Compare to candidate expressions
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CANDIDATE EXPRESSIONS FOR f_std AT FIXED POINT")
print("="*80)

# Get the converged value for identity
avg_id = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for seed in range(n_seeds):
    ent = Entangler(identity_corr, seed=seed).build()
    avg_id += ent.joint
avg_id /= n_seeds

current = avg_id.clone()
for _ in range(n_compositions):
    current = compose(current, avg_id)
f_t_final, f_s_final, f_std_final = s3_fractions(current)

print(f"\n  Converged f_std (identity, {n_seeds} seeds, {n_compositions} compositions):")
print(f"  f_std = {f_std_final:.8f}")

candidates = {
    "1/(H³+1) = 1/28": 1 / (H**3 + 1),
    "1/H³ = BORN_FLOOR": BORN_FLOOR,
    "BORN_FLOOR²": BORN_FLOOR**2,
    "K*/(H²+1)": K_STAR / (H**2 + 1),
    "K*/H(H+1)": K_STAR / (H * (H + 1)),
    "Δ/(H²+1)": DELTA / (H**2 + 1),
    "1/(H²(H+1))": 1 / (H**2 * (H + 1)),
    "1/(2H(H²+1))": 1 / (2 * H * (H**2 + 1)),
    "(H-1)/(H(H²+1))": (H-1) / (H * (H**2 + 1)),
    "1/H⁴": 1 / H**4,
    "2/(H²(H+1))": 2 / (H**2 * (H + 1)),
    "K*²": K_STAR**2,
    "Δ²": DELTA**2,
    "1/(H(H²+H+1))": 1 / (H * (H**2 + H + 1)),
    "BORN_FLOOR × K*": BORN_FLOOR * K_STAR,
    "1/30": 1/30,
    "BORN_FLOOR/H": BORN_FLOOR / H,
}

print(f"\n  {'expression':>25s}  {'value':>10s}  {'diff':>12s}")
for expr, val in sorted(candidates.items(), key=lambda x: abs(x[1] - f_std_final)):
    diff = val - f_std_final
    print(f"  {expr:>25s}  {val:>10.7f}  {diff:>+12.7f}")


# ═══════════════════════════════════════════════════════════════
#  Sign sector decay rate in the cycle crystal
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SIGN SECTOR DECAY RATE")
print("="*80)

avg_cyc = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
for seed in range(n_seeds):
    ent = Entangler(cycle_corr, seed=seed).build()
    avg_cyc += ent.joint
avg_cyc /= n_seeds

current = avg_cyc.clone()
prev_sign = None

print(f"\n  Cycle crystal ({n_seeds}-seed average):")
print(f"  {'comp':>4s}  {'f_sign':>10s}  {'ratio':>10s}  {'rate':>10s}  {'rate/Δ':>8s}")

for comp in range(10):
    f_t, f_s, f_st = s3_fractions(current)

    if prev_sign is not None and prev_sign > 1e-6 and f_s > 1e-6:
        ratio = f_s / prev_sign
        rate = -math.log(ratio)
        rate_delta = rate / DELTA
    elif prev_sign is not None and f_s < 1e-6:
        ratio = 0
        rate = float('inf')
        rate_delta = float('inf')
    else:
        ratio = None
        rate = None
        rate_delta = None

    ratio_str = f"{ratio:.6f}" if ratio is not None else "---"
    rate_str = f"{rate:.5f}" if rate is not None and rate != float('inf') else "---"
    delta_str = f"{rate_delta:.3f}" if rate_delta is not None and rate_delta != float('inf') else "---"

    print(f"  {comp:4d}  {f_s:>10.7f}  {ratio_str:>10s}  {rate_str:>10s}  {delta_str:>8s}")
    prev_sign = f_s
    current = compose(current, avg_cyc)


# ═══════════════════════════════════════════════════════════════
#  Per-seed sign decay to get statistics
# ═══════════════════════════════════════════════════════════════

print(f"\n\n  Per-seed sign decay (first composition step):")

first_step_ratios = []
for seed in range(n_seeds):
    ent = Entangler(cycle_corr, seed=seed).build()
    crystal = ent.joint.clone()

    f_t0, f_s0, f_st0 = s3_fractions(crystal)
    composed = compose(crystal, crystal)
    f_t1, f_s1, f_st1 = s3_fractions(composed)

    if f_s0 > 1e-4 and f_s1 > 1e-6:
        first_step_ratios.append(f_s1 / f_s0)

if first_step_ratios:
    avg_ratio = sum(first_step_ratios) / len(first_step_ratios)
    std_ratio = (sum((r - avg_ratio)**2 for r in first_step_ratios) / len(first_step_ratios)) ** 0.5
    avg_rate = -math.log(avg_ratio) if avg_ratio > 0 else float('inf')

    print(f"    n = {len(first_step_ratios)} seeds with measurable sign")
    print(f"    avg ratio = {avg_ratio:.6f} ± {std_ratio:.6f}")
    print(f"    avg rate = {avg_rate:.5f} = {avg_rate/DELTA:.3f}Δ")
    print(f"    1/2 = {0.5:.6f}  (if rate = ln(2))")
    print(f"    ln(2) = {math.log(2):.5f} = {math.log(2)/DELTA:.3f}Δ")


# ═══════════════════════════════════════════════════════════════
#  Single-seed convergence to check f_std scatter
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("PER-SEED f_std AT FIXED POINT")
print("="*80)

seed_fstd = []
for seed in range(n_seeds):
    ent = Entangler(identity_corr, seed=seed).build()
    current = ent.joint.clone()
    original = ent.joint.clone()
    for _ in range(30):
        current = compose(current, original)
    f_t, f_s, f_st = s3_fractions(current)
    seed_fstd.append(f_st)

sorted_fstd = sorted(seed_fstd)
avg_fstd = sum(seed_fstd) / len(seed_fstd)
std_fstd = (sum((v - avg_fstd)**2 for v in seed_fstd) / len(seed_fstd)) ** 0.5

print(f"\n  {n_seeds} seeds, 30 compositions each:")
print(f"  avg f_std = {avg_fstd:.6f} ± {std_fstd:.6f}")
print(f"  median = {sorted_fstd[n_seeds//2]:.6f}")
print(f"  min = {sorted_fstd[0]:.6f}, max = {sorted_fstd[-1]:.6f}")
print(f"  10% = {sorted_fstd[n_seeds//10]:.6f}, 90% = {sorted_fstd[9*n_seeds//10]:.6f}")

print(f"\n  Reference: 1/(H³+1) = {1/(H**3+1):.6f}")
print(f"  Reference: BORN_FLOOR = {BORN_FLOOR:.6f}")


print(f"\n\n{'='*80}")
print("WHAT THE FIXED POINT GAUGE CONTENT REVEALS")
print("="*80)
