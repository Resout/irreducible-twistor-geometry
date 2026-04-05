"""
What determines the type-specific eigenvalue decay rate?

Under self-composition, |λ₁|/|λ₀| decays at a constant rate:
  proportional: 0.503 (1.89Δ)
  inverse:      0.531 (2.00Δ — EXACTLY 2Δ!)
  modular:      0.831 (3.13Δ)
  exponential:  0.753 (2.83Δ)
  logarithmic:  0.909 (3.42Δ)
  quadratic:    0.638 (2.40Δ)

What determines these rates? Is it:
  (a) The eigenvalue phases? (phase rotation → different decay)
  (b) The correlation matrix structure?
  (c) The deviation from normality?
  (d) The singular value spectrum?

Also: the PHASES rotate linearly. What determines the phase step?
Is the eigenvalue of composition (magnitude + phase) a complete
invariant of the crystal?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, BORN_FLOOR, K_STAR, DELTA, schmidt_number
from solver.crystals import Entangler, compose, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)


def get_eigendata(joint):
    """Get eigenvalues sorted by magnitude."""
    eigenvalues, _ = torch.linalg.eig(joint)
    mags = eigenvalues.abs()
    idx = mags.argsort(descending=True)
    return eigenvalues[idx]


def composition_eigenvalue_rates(joint, n_steps=15):
    """Track eigenvalue ratios and phases through composition."""
    current = joint.clone()
    records = []

    for step in range(n_steps):
        ev = get_eigendata(current)
        records.append({
            "step": step,
            "eigenvalues": ev.clone(),
            "magnitudes": ev.abs().tolist(),
            "phases": [math.atan2(e.imag.item(), e.real.item()) / math.pi for e in ev],
        })
        current = compose(current, joint)

    # Compute decay rates from eigenvalue ratios
    rates = {}
    for i in range(1, MASS_DIM):
        rs = []
        for step in range(2, min(10, n_steps)):
            r_prev = records[step-1]["magnitudes"][i] / records[step-1]["magnitudes"][0]
            r_curr = records[step]["magnitudes"][i] / records[step]["magnitudes"][0]
            if r_prev > 1e-12 and r_curr > 1e-12:
                rs.append(math.log(r_prev / r_curr))
        rates[i] = sum(rs) / len(rs) if rs else 0

    # Compute phase steps
    phase_steps = {}
    for i in range(MASS_DIM):
        ps = []
        for step in range(1, min(8, n_steps)):
            dp = records[step]["phases"][i] - records[step-1]["phases"][i]
            # Unwrap
            while dp > 1: dp -= 2
            while dp < -1: dp += 2
            ps.append(dp)
        phase_steps[i] = sum(ps) / len(ps) if ps else 0

    return records, rates, phase_steps


# ═══════════════════════════════════════════════════════════════
#  Collect rates and phase steps for all crystal types
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("EIGENVALUE DECAY RATES AND PHASE STEPS BY CRYSTAL TYPE")
print("=" * 70)

results = {}
for name, corr in RELATIONSHIP_SIGNATURES.items():
    ent = Entangler(corr).build()
    records, rates, phase_steps = composition_eigenvalue_rates(ent.joint)

    results[name] = {
        "rates": rates,
        "phase_steps": phase_steps,
        "eigenvalues": records[0]["eigenvalues"],
        "schmidt": schmidt_number(ent.joint),
        "corr": corr,
    }

print(f"\n  {'Crystal':>15s} {'Schmidt':>8s} {'r₁/Δ':>8s} {'r₂/Δ':>8s} {'r₃/Δ':>8s} {'φ₁ step':>8s} {'φ₂ step':>8s} {'φ₃ step':>8s}")
print("  " + "-" * 75)

for name in RELATIONSHIP_SIGNATURES:
    r = results[name]
    rates = r["rates"]
    ps = r["phase_steps"]
    print(f"  {name:>15s} {r['schmidt']:>8.3f} "
          f"{rates[1]/DELTA:>8.3f} {rates[2]/DELTA:>8.3f} {rates[3]/DELTA:>8.3f} "
          f"{ps[1]:>8.4f} {ps[2]:>8.4f} {ps[3]:>8.4f}")


# ═══════════════════════════════════════════════════════════════
#  What correlates with the decay rate?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("WHAT DETERMINES THE DECAY RATE?")
print("=" * 70)

# Collect candidate predictors
data = []
for name in RELATIONSHIP_SIGNATURES:
    r = results[name]
    corr = r["corr"]

    # Correlation matrix properties
    corr_diag_sum = sum(corr[i, i].item() for i in range(H))
    corr_off_diag = corr.sum().item() - corr_diag_sum
    corr_trace = corr_diag_sum / H  # average diagonal
    corr_entropy = 0
    for i in range(H):
        row = corr[i] / corr[i].sum()
        for j in range(H):
            if row[j] > 0:
                corr_entropy -= row[j].item() * math.log(row[j].item())
    corr_entropy /= H  # average row entropy

    # Frobenius distance from identity
    id_corr = torch.eye(H)
    corr_dist_id = (corr - id_corr).pow(2).sum().sqrt().item()

    # Eigenvalues of the correlation matrix itself
    corr_eig = torch.linalg.eigvalsh(corr)
    corr_condition = corr_eig[-1].item() / corr_eig[0].item() if corr_eig[0] > 0 else float('inf')

    data.append({
        "name": name,
        "r1": r["rates"][1],
        "r1_over_delta": r["rates"][1] / DELTA,
        "schmidt": r["schmidt"],
        "corr_trace": corr_trace,
        "corr_entropy": corr_entropy,
        "corr_dist_id": corr_dist_id,
        "corr_condition": corr_condition,
    })

# Print candidate predictors
print(f"\n  {'Crystal':>15s} {'r₁/Δ':>8s} {'Schmidt':>8s} {'diag':>6s} {'entropy':>8s} {'d(id)':>8s} {'cond':>8s}")
print("  " + "-" * 65)
for d in data:
    print(f"  {d['name']:>15s} {d['r1_over_delta']:>8.3f} {d['schmidt']:>8.3f} "
          f"{d['corr_trace']:>6.2f} {d['corr_entropy']:>8.4f} "
          f"{d['corr_dist_id']:>8.4f} {d['corr_condition']:>8.2f}")

# Correlations
def corr_coef(xs, ys):
    n = len(xs)
    mx, my = sum(xs)/n, sum(ys)/n
    c = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx = sum((x-mx)**2 for x in xs)**0.5
    sy = sum((y-my)**2 for y in ys)**0.5
    return c/(sx*sy) if sx>1e-10 and sy>1e-10 else 0

r1_vals = [d["r1_over_delta"] for d in data]
print(f"\n  Correlations with r₁/Δ:")
for pred in ["schmidt", "corr_trace", "corr_entropy", "corr_dist_id", "corr_condition"]:
    pred_vals = [d[pred] for d in data]
    r = corr_coef(r1_vals, pred_vals)
    print(f"    {pred:>15s}: r = {r:+.3f}")


# ═══════════════════════════════════════════════════════════════
#  The eigenvalue as invariant: (|λ₁/λ₀|, φ₁ step)
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("EIGENVALUE SIGNATURE: (decay_rate, phase_step)")
print("Is this pair a complete invariant of the crystal?")
print("=" * 70)

print(f"\n  {'Crystal':>15s} {'(r₁/Δ, φ₁/π)':>20s} {'(r₂/Δ, φ₂/π)':>20s} {'(r₃/Δ, φ₃/π)':>20s}")
print("  " + "-" * 80)

for name in RELATIONSHIP_SIGNATURES:
    r = results[name]
    sigs = []
    for i in range(1, MASS_DIM):
        rate_ratio = r["rates"][i] / DELTA
        phase = r["phase_steps"][i]
        sigs.append(f"({rate_ratio:.2f}, {phase:+.3f})")
    print(f"  {name:>15s} {sigs[0]:>20s} {sigs[1]:>20s} {sigs[2]:>20s}")


# ═══════════════════════════════════════════════════════════════
#  Do different seeds give different rates?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SEED DEPENDENCE OF EIGENVALUE DECAY RATES")
print("=" * 70)

id_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

seed_rates = []
for seed in range(20):
    ent = Entangler(id_corr, seed=seed).build()
    _, rates, phases = composition_eigenvalue_rates(ent.joint, n_steps=12)
    seed_rates.append({
        "seed": seed,
        "r1": rates[1] / DELTA,
        "r2": rates[2] / DELTA,
        "r3": rates[3] / DELTA,
        "p1": phases[1],
        "p2": phases[2],
        "p3": phases[3],
    })

# Statistics
for key in ["r1", "r2", "r3"]:
    vals = [sr[key] for sr in seed_rates]
    avg = sum(vals) / len(vals)
    std = (sum((v - avg)**2 for v in vals) / len(vals)) ** 0.5
    print(f"  {key}/Δ: {avg:.4f} ± {std:.4f} (CV = {std/abs(avg)*100:.1f}%)")

for key in ["p1", "p2", "p3"]:
    vals = [sr[key] for sr in seed_rates]
    avg = sum(vals) / len(vals)
    std = (sum((v - avg)**2 for v in vals) / len(vals)) ** 0.5
    print(f"  φ{key[-1]}/π: {avg:.4f} ± {std:.4f}")


# ═══════════════════════════════════════════════════════════════
#  The magic of the inverse: why EXACTLY 2Δ?
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("WHY DOES THE INVERSE CRYSTAL GIVE EXACTLY 2Δ?")
print("=" * 70)

# The inverse correlation is the anti-diagonal identity.
# It maps low→high, high→low, mid→mid.
# This is the ONLY non-trivial involution in S₃.

# What about the other involutions?
# S₃ has 3 involutions (transpositions): (01), (02), (12)
# Plus the identity (order 1) and two 3-cycles (order 3)

involutions = {
    "swap(0,2)": torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32),  # inverse
    "swap(0,1)": torch.tensor([[0,1,0],[1,0,0],[0,0,1]], dtype=torch.float32),
    "swap(1,2)": torch.tensor([[1,0,0],[0,0,1],[0,1,0]], dtype=torch.float32),
}

three_cycles = {
    "cycle(012)": torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32),
    "cycle(021)": torch.tensor([[0,0,1],[1,0,0],[0,1,0]], dtype=torch.float32),
}

print(f"\n  Involutions (order 2):")
for name, corr in involutions.items():
    ent = Entangler(corr).build()
    _, rates, _ = composition_eigenvalue_rates(ent.joint, n_steps=12)
    print(f"    {name}: r₁/Δ = {rates[1]/DELTA:.4f}")

print(f"\n  3-cycles (order 3):")
for name, corr in three_cycles.items():
    ent = Entangler(corr).build()
    _, rates, _ = composition_eigenvalue_rates(ent.joint, n_steps=12)
    print(f"    {name}: r₁/Δ = {rates[1]/DELTA:.4f}")

print(f"\n  Identity (order 1):")
ent = Entangler(torch.eye(3)).build()
_, rates, _ = composition_eigenvalue_rates(ent.joint, n_steps=12)
print(f"    identity: r₁/Δ = {rates[1]/DELTA:.4f}")


print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
