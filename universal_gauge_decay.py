"""
Is the gauge-variant decay rate universal per ORDER-CYCLE?

Sign of cycle: decays at ~1.0 per 3 steps (order 3)
Standard of transposition: decays at ~1.0 per 2 steps (order 2)

Hypothesis: the gauge-variant content loses exactly e⁻¹ per
trip around the group element's cyclic subgroup.

If true: Δ_gauge = 1/ord(σ) per composition step.
Or equivalently: Δ_gauge × ord(σ) = 1 for all σ.

This would be a UNIVERSAL decay rate per order-cycle, connecting
the mass gap to the group order through the gauge structure.

Also check: do the K* numerator cyclotomic polynomials appear?
K* = Φ₆(H)/H(H²+1), K_peak = Φ₃(H)/2H(H²+1)
where Φ₃(x) = x²+x+1, Φ₆(x) = x²-x+1.
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


# Build 200-seed averaged crystals for all S₃ elements
identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
swap01_corr = torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32)   # (01)
swap02_corr = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)   # (02)
swap12_corr = torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32)   # (12)
cycle012_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32) # (012)
cycle021_corr = torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32) # (021)

elements = {
    "e": (identity_corr, 1),
    "(01)": (swap01_corr, 2),
    "(02)": (swap02_corr, 2),
    "(12)": (swap12_corr, 2),
    "(012)": (cycle012_corr, 3),
    "(021)": (cycle021_corr, 3),
}

n_seeds = 200
n_comp = 15

print("=" * 80)
print("UNIVERSAL GAUGE DECAY: ONE e-FOLD PER ORDER-CYCLE?")
print("=" * 80)

# Build seed-averaged crystals
crystals = {}
for name, (corr, order) in elements.items():
    avg = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(n_seeds):
        avg += Entangler(corr, seed=seed).build().joint
    avg /= n_seeds
    crystals[name] = avg


# ═══════════════════════════════════════════════════════════════
#  Track the DOMINANT gauge-variant sector for each element
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Gauge-variant content through composition ---\n")

decay_data = {}  # name → list of (step, gauge_content) at same-phase steps

for name, (corr, order) in elements.items():
    crystal = crystals[name]
    current = crystal.clone()

    gauge_at_phase = []  # gauge content at steps ≡ 0 mod order

    all_steps = []
    for comp in range(n_comp):
        f_t, f_s, f_st = s3_fractions(current)

        # The "gauge content" is whichever of sign/standard is dominant for this element
        if order == 3:
            gauge = f_s  # 3-cycles live in sign
        elif order == 2:
            gauge = f_st  # transpositions live in standard
        else:
            gauge = f_st  # identity: standard (tiny)

        all_steps.append((comp, f_t, f_s, f_st, gauge))

        # Record at same-phase steps (comp ≡ 0 mod order)
        if comp % order == 0:
            gauge_at_phase.append((comp, gauge))

        current = compose(current, crystal)

    decay_data[name] = (order, gauge_at_phase, all_steps)


# ═══════════════════════════════════════════════════════════════
#  Compute decay rates per ORDER-CYCLE
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print("DECAY RATES PER ORDER-CYCLE")
print("="*80)

print(f"\n  For each element, compare gauge content at steps 0, ord, 2·ord, ...")
print(f"  Rate per order-cycle = -ln(gauge(n·ord) / gauge((n-1)·ord))")
print(f"  If universal: rate per cycle = 1 (= e-fold per trip)")
print(f"  Δ = {DELTA:.5f}")

for name, (order, gauge_at_phase, all_steps) in decay_data.items():
    if order == 1:
        continue  # identity has no gauge content

    print(f"\n  {name} (order {order}):")
    print(f"  {'cycle':>6s}  {'step':>4s}  {'gauge':>12s}  {'ratio':>10s}  {'rate/cycle':>12s}  {'rate/step':>12s}  {'rate/Δ':>8s}")

    for i in range(len(gauge_at_phase)):
        step, gauge = gauge_at_phase[i]
        if i > 0:
            prev_gauge = gauge_at_phase[i-1][1]
            if prev_gauge > 1e-10 and gauge > 1e-10:
                ratio = gauge / prev_gauge
                rate_cycle = -math.log(ratio)
                rate_step = rate_cycle / order
                rate_delta = rate_step / DELTA
            else:
                ratio = 0
                rate_cycle = float('inf')
                rate_step = float('inf')
                rate_delta = float('inf')
            r_str = f"{ratio:.6f}" if ratio > 0 else "→0"
            rc_str = f"{rate_cycle:.5f}" if rate_cycle != float('inf') else "∞"
            rs_str = f"{rate_step:.5f}" if rate_step != float('inf') else "∞"
            rd_str = f"{rate_delta:.3f}" if rate_delta != float('inf') else "∞"
        else:
            r_str = "---"
            rc_str = "---"
            rs_str = "---"
            rd_str = "---"

        print(f"  {i:6d}  {step:4d}  {gauge:>12.8f}  {r_str:>10s}  {rc_str:>12s}  {rs_str:>12s}  {rd_str:>8s}")


# ═══════════════════════════════════════════════════════════════
#  Summary: average rate per order-cycle for each conjugacy class
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SUMMARY: AVERAGE DECAY RATE PER ORDER-CYCLE")
print("="*80)

print(f"\n  {'element':>8s}  {'order':>5s}  {'irrep':>10s}  {'rate/cycle':>12s}  {'rate/step':>12s}  {'rate/Δ':>8s}")
print("-" * 65)

for name, (order, gauge_at_phase, _) in decay_data.items():
    if order == 1:
        continue

    rates = []
    for i in range(1, len(gauge_at_phase)):
        prev = gauge_at_phase[i-1][1]
        curr = gauge_at_phase[i][1]
        if prev > 1e-10 and curr > 1e-10:
            rates.append(-math.log(curr / prev))

    if rates:
        avg_rate = sum(rates) / len(rates)
        avg_per_step = avg_rate / order
        avg_per_delta = avg_per_step / DELTA
        irrep = "sign" if order == 3 else "standard"
        print(f"  {name:>8s}  {order:>5d}  {irrep:>10s}  {avg_rate:>12.5f}  {avg_per_step:>12.5f}  {avg_per_delta:>8.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test the hypothesis: rate per cycle = constant?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("IS THE RATE PER ORDER-CYCLE UNIVERSAL?")
print("="*80)

all_rates_per_cycle = []
for name, (order, gauge_at_phase, _) in decay_data.items():
    if order == 1:
        continue
    for i in range(1, min(3, len(gauge_at_phase))):  # first 2 cycles only (cleanest)
        prev = gauge_at_phase[i-1][1]
        curr = gauge_at_phase[i][1]
        if prev > 1e-10 and curr > 1e-10:
            rate = -math.log(curr / prev)
            all_rates_per_cycle.append((name, order, i, rate))

if all_rates_per_cycle:
    rates_only = [r for _, _, _, r in all_rates_per_cycle]
    avg = sum(rates_only) / len(rates_only)
    std = (sum((r - avg)**2 for r in rates_only) / len(rates_only)) ** 0.5

    print(f"\n  All elements, first 2 order-cycles:")
    for name, order, cycle, rate in all_rates_per_cycle:
        print(f"    {name:>8s} (order {order}), cycle {cycle}: rate = {rate:.5f}")

    print(f"\n  Average rate per order-cycle: {avg:.5f} ± {std:.5f}")
    print(f"  Candidate: 1.0 → diff = {avg - 1.0:+.5f}")
    print(f"  Candidate: 2Δ = {2*DELTA:.5f} → diff = {avg - 2*DELTA:+.5f}")
    print(f"  Candidate: 3Δ = {3*DELTA:.5f} → diff = {avg - 3*DELTA:+.5f}")
    print(f"  Candidate: Δ×π = {DELTA*math.pi:.5f} → diff = {avg - DELTA*math.pi:+.5f}")
    print(f"  Candidate: ln(H+1) = {math.log(H+1):.5f} → diff = {avg - math.log(H+1):+.5f}")
    print(f"  Candidate: ln(H²+1) = {math.log(H**2+1):.5f} → diff = {avg - math.log(H**2+1):+.5f}")


# ═══════════════════════════════════════════════════════════════
#  Cyclotomic connection
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CYCLOTOMIC POLYNOMIALS IN THE STRUCTURE")
print("="*80)

print(f"""
  K* = Φ₆(H) / (H(H²+1))     where Φ₆(x) = x²-x+1
  K_peak = Φ₃(H) / (2H(H²+1)) where Φ₃(x) = x²+x+1

  At H=3:
    Φ₃(3) = 13  (nucleation numerator)
    Φ₆(3) = 7   (K* numerator)

  Φ₃(x) · Φ₆(x) = x⁴+x²+1
  At H=3: 13 × 7 = 91 = 81+9+1 = 3⁴+3²+1 ✓

  Φ₃ is the minimal polynomial of primitive 3rd roots of unity (ω, ω²)
  Φ₆ is the minimal polynomial of primitive 6th roots of unity (-ω, -ω²)

  The 3rd roots of unity ARE the eigenvalues of the 3-cycle (012).
  The 6th roots of unity = (3rd roots) * (plus/minus 1) = include the sign character.

  Connection: the SIGN eigenvalue of the 3-cycle under composition
  should be related to Φ₆, and the STANDARD eigenvalue to Φ₃.
""")

print(f"  Φ₃(H)/Φ₆(H) = {(H**2+H+1)/(H**2-H+1):.6f} = 13/7 = {13/7:.6f}")
print(f"  K_peak/K* = Φ₃(H)/(2·Φ₆(H)) = 13/14 = {13/14:.6f}")
print(f"  2·K_peak/K* = Φ₃(H)/Φ₆(H) = {(H**2+H+1)/(H**2-H+1):.6f}")

# The ratio Φ₃/Φ₆ should appear somewhere in the decay rates
ratio = (H**2 + H + 1) / (H**2 - H + 1)
print(f"\n  Φ₃/Φ₆ = {ratio:.6f}")
print(f"  ln(Φ₃/Φ₆) = {math.log(ratio):.6f} = {math.log(ratio)/DELTA:.3f}Δ")


print(f"\n\n{'='*80}")
print("WHAT THE UNIVERSAL DECAY REVEALS")
print("="*80)
