"""
Definitive test: is the bijection eigenvalue decay rate universally 2Δ?

Previous finding: inverse gives 2.00Δ exactly, identity averages 1.98Δ
over 20 seeds. But 20 seeds is too few. This tests ALL 6 S₃ elements
over 200 seeds each.

If the mean converges to 2.00Δ for all bijections, that's a theorem.
If different permutation types give different rates, the representation
theory thread opens further.

Also: do the two phase families (0.12π vs 0.85π) correlate with the
sign representation? Transpositions have sign -1, even permutations +1.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, DELTA, schmidt_number
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)

N_SEEDS = 200

# All 6 elements of S₃ as correlation matrices
S3_ELEMENTS = {
    "identity":   torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32),
    "swap(0,1)":  torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "swap(0,2)":  torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),  # = "inverse"
    "swap(1,2)":  torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
    "cycle(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "cycle(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

# Sign representation: even permutations = +1, odd = -1
SIGNS = {
    "identity": +1, "swap(0,1)": -1, "swap(0,2)": -1, "swap(1,2)": -1,
    "cycle(012)": +1, "cycle(021)": +1,
}

# Conjugacy classes
CONJUGACY = {
    "identity": "e", "swap(0,1)": "transposition", "swap(0,2)": "transposition",
    "swap(1,2)": "transposition", "cycle(012)": "3-cycle", "cycle(021)": "3-cycle",
}


def get_eigendata(joint):
    eigenvalues, _ = torch.linalg.eig(joint)
    mags = eigenvalues.abs()
    idx = mags.argsort(descending=True)
    return eigenvalues[idx]


def measure_rates_and_phases(joint, n_steps=15):
    """Measure eigenvalue decay rates and phase steps under self-composition."""
    current = joint.clone()
    records = []
    for step in range(n_steps):
        ev = get_eigendata(current)
        records.append({
            "magnitudes": ev.abs().tolist(),
            "phases": [math.atan2(e.imag.item(), e.real.item()) / math.pi for e in ev],
        })
        current = compose(current, joint)

    # Decay rate for λ₁/λ₀
    rates = {}
    for i in range(1, MASS_DIM):
        rs = []
        for step in range(2, min(10, n_steps)):
            r_prev = records[step-1]["magnitudes"][i] / max(records[step-1]["magnitudes"][0], 1e-15)
            r_curr = records[step]["magnitudes"][i] / max(records[step]["magnitudes"][0], 1e-15)
            if r_prev > 1e-12 and r_curr > 1e-12:
                rs.append(math.log(r_prev / r_curr))
        rates[i] = sum(rs) / len(rs) if rs else 0

    # Phase steps
    phase_steps = {}
    for i in range(MASS_DIM):
        ps = []
        for step in range(1, min(8, n_steps)):
            dp = records[step]["phases"][i] - records[step-1]["phases"][i]
            while dp > 1: dp -= 2
            while dp < -1: dp += 2
            ps.append(dp)
        phase_steps[i] = sum(ps) / len(ps) if ps else 0

    return rates, phase_steps


# ═══════════════════════════════════════════════════════════════
#  Main: 200 seeds × 6 permutations
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print(f"BIJECTION RATE TEST: {N_SEEDS} SEEDS × 6 S₃ ELEMENTS")
print(f"Question: is r₁ = 2Δ universally for all bijections?")
print(f"Δ = {DELTA:.6f},  2Δ = {2*DELTA:.6f}")
print("=" * 80)

all_results = {}

for name, corr in S3_ELEMENTS.items():
    r1_vals = []
    p1_vals = []
    schmidt_vals = []

    for seed in range(N_SEEDS):
        ent = Entangler(corr, seed=seed).build()
        rates, phases = measure_rates_and_phases(ent.joint, n_steps=12)
        r1_vals.append(rates[1] / DELTA)
        p1_vals.append(phases[1])
        schmidt_vals.append(schmidt_number(ent.joint))

    r1_avg = sum(r1_vals) / len(r1_vals)
    r1_std = (sum((v - r1_avg)**2 for v in r1_vals) / len(r1_vals)) ** 0.5
    p1_avg = sum(p1_vals) / len(p1_vals)
    p1_std = (sum((v - p1_avg)**2 for v in p1_vals) / len(p1_vals)) ** 0.5
    sc_avg = sum(schmidt_vals) / len(schmidt_vals)

    all_results[name] = {
        "r1_avg": r1_avg, "r1_std": r1_std, "r1_vals": r1_vals,
        "p1_avg": p1_avg, "p1_std": p1_std, "p1_vals": p1_vals,
        "schmidt_avg": sc_avg,
        "sign": SIGNS[name], "conjugacy": CONJUGACY[name],
    }

# ═══════════════════════════════════════════════════════════════
#  Results table
# ═══════════════════════════════════════════════════════════════

print(f"\n  {'Element':>12s} {'Class':>13s} {'Sign':>5s} {'r₁/Δ':>12s} {'φ₁/π':>12s} {'Schmidt':>8s}")
print("  " + "-" * 75)

for name in S3_ELEMENTS:
    r = all_results[name]
    print(f"  {name:>12s} {r['conjugacy']:>13s} {r['sign']:>+5d} "
          f"{r['r1_avg']:>6.3f}±{r['r1_std']:.3f} "
          f"{r['p1_avg']:>6.4f}±{r['p1_std']:.4f} "
          f"{r['schmidt_avg']:>8.3f}")


# ═══════════════════════════════════════════════════════════════
#  Grouped by conjugacy class
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("GROUPED BY CONJUGACY CLASS")
print("="*80)

for cls in ["e", "transposition", "3-cycle"]:
    members = [n for n in S3_ELEMENTS if CONJUGACY[n] == cls]
    all_r1 = []
    all_p1 = []
    for n in members:
        all_r1.extend(all_results[n]["r1_vals"])
        all_p1.extend(all_results[n]["p1_vals"])

    r1_avg = sum(all_r1) / len(all_r1)
    r1_std = (sum((v - r1_avg)**2 for v in all_r1) / len(all_r1)) ** 0.5
    r1_sem = r1_std / len(all_r1)**0.5
    p1_avg = sum(all_p1) / len(all_p1)
    p1_std = (sum((v - p1_avg)**2 for v in all_p1) / len(all_p1)) ** 0.5

    deviation_from_2 = abs(r1_avg - 2.0) / r1_sem if r1_sem > 0 else 0

    print(f"\n  {cls} (n={len(all_r1)}):")
    print(f"    r₁/Δ = {r1_avg:.4f} ± {r1_std:.4f}  (SEM = {r1_sem:.4f})")
    print(f"    Deviation from 2.0: {deviation_from_2:.1f}σ")
    print(f"    φ₁/π = {p1_avg:.4f} ± {p1_std:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Grouped by sign representation
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("GROUPED BY SIGN REPRESENTATION")
print("Question: do even and odd permutations have different phase families?")
print("="*80)

for sign_val, sign_name in [(+1, "even (sign +1)"), (-1, "odd (sign -1)")]:
    members = [n for n in S3_ELEMENTS if SIGNS[n] == sign_val]
    all_r1 = []
    all_p1 = []
    for n in members:
        all_r1.extend(all_results[n]["r1_vals"])
        all_p1.extend(all_results[n]["p1_vals"])

    r1_avg = sum(all_r1) / len(all_r1)
    p1_avg = sum(all_p1) / len(all_p1)
    p1_std = (sum((v - p1_avg)**2 for v in all_p1) / len(all_p1)) ** 0.5

    print(f"\n  {sign_name} ({', '.join(members)}):")
    print(f"    r₁/Δ = {r1_avg:.4f}")
    print(f"    φ₁/π = {p1_avg:.4f} ± {p1_std:.4f}")


# ═══════════════════════════════════════════════════════════════
#  The definitive question: is 2Δ universal?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("VERDICT: IS r₁ = 2Δ UNIVERSAL FOR ALL BIJECTIONS?")
print("="*80)

all_r1_all = []
for n in S3_ELEMENTS:
    all_r1_all.extend(all_results[n]["r1_vals"])

grand_avg = sum(all_r1_all) / len(all_r1_all)
grand_std = (sum((v - grand_avg)**2 for v in all_r1_all) / len(all_r1_all)) ** 0.5
grand_sem = grand_std / len(all_r1_all)**0.5

print(f"\n  Grand mean r₁/Δ = {grand_avg:.5f} ± {grand_std:.5f}")
print(f"  SEM = {grand_sem:.5f}")
print(f"  Deviation from 2.000: {abs(grand_avg - 2.0)/grand_sem:.1f}σ")
print(f"  95% CI: [{grand_avg - 1.96*grand_sem:.4f}, {grand_avg + 1.96*grand_sem:.4f}]")

if abs(grand_avg - 2.0) < 2 * grand_sem:
    print(f"\n  ★ CONSISTENT with r₁ = 2Δ universal for all bijections")
elif abs(grand_avg - 2.0) < 5 * grand_sem:
    print(f"\n  ◇ Marginal — within 5σ but not 2σ of 2Δ")
else:
    print(f"\n  ✗ NOT consistent with r₁ = 2Δ universal (>{abs(grand_avg - 2.0)/grand_sem:.0f}σ)")
