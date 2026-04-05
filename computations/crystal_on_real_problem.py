"""
The crystal framework applied to a REAL math problem.

We've discovered: the crystal classifies by information structure,
the eigenvalues decay perfectly exponentially, composition is
decorrelation, the Born floor is a projector.

Now: what does all this look like on an actual problem?

Problem: "Find the remainder when 7^2025 is divided by 100."

This is a modular exponentiation problem. The crystal framework
was BUILT for this. But instead of just computing pow(7, 2025, 100),
let's watch the crystal at every level — eigenvalues, θ-fingerprint,
Born probabilities, slot measurement, phase structure — and see
if the abstract findings manifest in a concrete computation.

Then: a problem the crystal CAN'T solve directly.
"How many integers from 1 to 1000 have digit sum equal to 10?"
This has no period-3 structure. What does the crystal SAY about it?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import H, MASS_DIM, BORN_FLOOR, DELTA, K_STAR, schmidt_number, born_probabilities
from solver.crystals import Entangler, compose, classify_relationship, slot_measure
from solver.functional import build_function_crystal, function_to_correlation, _detect_output_period
from solver.resonance import find_resonant_frame
from solver import compute

torch.set_grad_enabled(False)


def theta_fingerprint(joint):
    fp = torch.stack([
        joint[0, 3], joint[1, 3], joint[2, 3],
        joint[3, 0], joint[3, 1], joint[3, 2],
    ])
    norm = fp.abs().pow(2).sum().sqrt()
    if norm > 1e-10:
        fp = fp / norm
    return fp

def theta_asymmetry(fp):
    return (fp[:3].abs() - fp[3:].abs()).abs().sum().item()

def full_crystal_report(fc, name="crystal"):
    """Complete measurement of a crystal."""
    joint = fc.joint
    sn = schmidt_number(joint)
    rel, rel_conf = classify_relationship(joint)
    sm = slot_measure(joint)
    dom_slot = max(sm, key=sm.get)
    fp = theta_fingerprint(joint)
    asymm = theta_asymmetry(fp)

    # Eigenvalues
    ev, _ = torch.linalg.eig(joint)
    ev_mags = ev.abs()
    idx = ev_mags.argsort(descending=True)
    ev = ev[idx]

    # Block decomposition
    total = joint.abs().pow(2).sum().item()
    hh = joint[:H, :H].abs().pow(2).sum().item() / total
    theta_total = 1 - hh

    # Marginal Born(θ)
    marginal = joint.abs().pow(2).sum(dim=1)
    marginal_born = marginal / marginal.sum()
    born_theta = marginal_born[H].item()

    print(f"\n  ═══ {name} ═══")
    print(f"    Schmidt: {sn:.3f}")
    print(f"    Relationship: {rel} (conf={rel_conf:.3f})")
    print(f"    Slot: {dom_slot} ({sm['generator']:.3f}/{sm['spectral']:.3f}/{sm['orbit']:.3f})")
    print(f"    θ-asymmetry: {asymm:.4f}")
    print(f"    h×h: {hh*100:.1f}%, θ-total: {theta_total*100:.1f}%")
    print(f"    Marginal Born(θ): {born_theta:.4f} ({born_theta/BORN_FLOOR:.1f}× floor)")
    print(f"    Eigenvalues |λᵢ|: [{', '.join(f'{e.abs().item():.4f}' for e in ev)}]")
    print(f"    |λ₁/λ₀|: {ev[1].abs().item()/ev[0].abs().item():.4f}")

    return sn, rel, dom_slot


# ═══════════════════════════════════════════════════════════════
#  Problem 1: 7^2025 mod 100
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("PROBLEM 1: Find the remainder when 7^2025 is divided by 100")
print("=" * 70)

# Step 0: exact answer
exact = pow(7, 2025, 100)
print(f"\n  Exact answer: {exact}")

# Step 1: Build the crystal
domain = range(1, 301)
f = lambda x: pow(7, int(x), 100)
fc = build_function_crystal("7^x mod 100", f, domain,
                            var_input="exponent", var_output="remainder")

# Step 2: Detect period
outputs = [f(x) for x in range(1, 101)]
period = _detect_output_period(outputs)
n_distinct = len(set(outputs))
print(f"\n  Output period (≤3): {period}")
print(f"  Distinct values: {n_distinct}")

# Full period detection
for p in range(1, 50):
    if all(outputs[i] == outputs[i % p] for i in range(p, min(len(outputs), p*5))):
        full_period = p
        break
else:
    full_period = None
print(f"  Full period: {full_period}")
print(f"  3 | period? {full_period and full_period % 3 == 0}")

# Step 3: Full crystal report
sn, rel, slot = full_crystal_report(fc, "7^x mod 100 (raw)")

# Step 4: Resonance search
print(f"\n  ── Resonance search ──")
res = find_resonant_frame(f, domain, max_multiplier=12, name="7^x mod 100")
if res.best:
    print(f"    Best frame: k={res.best.multiplier}, Schmidt={res.best.schmidt:.3f}, "
          f"period={res.best.period}")

# Step 5: Can the crystal answer?
result = fc.query(2025.0)
print(f"\n  ── Crystal query: 7^2025 mod 100 ──")
print(f"    Dominant bin: {result['dominant_bin']}")
print(f"    Confidence: {result['confidence']:.3f}")
print(f"    Output bins: {result['output_probs']}")

# Step 6: What the crystal DOES know
print(f"\n  ── What the crystal knows ──")
# The crystal knows the modular structure even if it can't give the exact answer
for x in [1, 2, 3, 4, 5, 10, 100]:
    r = fc.query(float(x))
    exact_val = f(x)
    print(f"    7^{x:3d} mod 100 = {exact_val:2d}, crystal bin={r['dominant_bin']} (conf={r['confidence']:.2f})")


# ═══════════════════════════════════════════════════════════════
#  Problem 2: Count integers 1-1000 with digit sum = 10
# ═══════════════════════════════════════════════════════════════

print("\n\n" + "=" * 70)
print("PROBLEM 2: How many integers from 1 to 1000 have digit sum = 10?")
print("(No periodic structure. What does the crystal see?)")
print("=" * 70)

# Exact answer
exact2 = sum(1 for n in range(1, 1001) if compute.digit_sum(n) == 10)
print(f"\n  Exact answer: {exact2}")

# Build crystal: f(n) = digit_sum(n)
f2 = lambda n: compute.digit_sum(max(1, int(n)))
fc2 = build_function_crystal("digit_sum", f2, range(1, 1001))

sn2, rel2, slot2 = full_crystal_report(fc2, "digit_sum(n)")

# Build crystal: f(n) = 1 if digit_sum(n) == 10 else 0
f3 = lambda n: 1 if compute.digit_sum(max(1, int(n))) == 10 else 0
fc3 = build_function_crystal("digit_sum==10", f3, range(1, 1001))

sn3, rel3, slot3 = full_crystal_report(fc3, "digit_sum(n) == 10 indicator")

# What does the crystal say about THIS problem?
print(f"\n  ── Crystal diagnosis ──")
print(f"    digit_sum crystal: Schmidt={sn2:.3f}, slot={slot2}")
print(f"    indicator crystal: Schmidt={sn3:.3f}, slot={slot3}")
if sn3 < 1.5:
    print(f"    Crystal says: NO STRUCTURE DETECTED (Schmidt below threshold)")
    print(f"    Recommendation: exact enumeration, not crystal computation")
elif sn3 < 2.2:
    print(f"    Crystal says: WEAK structure (Schmidt in noise range)")
    print(f"    Recommendation: use with caution, verify against enumeration")
else:
    print(f"    Crystal says: STRUCTURE DETECTED")


# ═══════════════════════════════════════════════════════════════
#  Problem 3: Fibonacci mod 7
# ═══════════════════════════════════════════════════════════════

print("\n\n" + "=" * 70)
print("PROBLEM 3: Find F_100 mod 7 (Fibonacci number mod 7)")
print("(Pisano period = 16, NOT divisible by 3)")
print("=" * 70)

exact3 = compute.linear_recurrence_mod([1, 1], [0, 1], 100, 7)
print(f"\n  Exact answer: F_100 mod 7 = {exact3}")

# Build crystal from Fibonacci mod 7
fib_vals = [0, 1]
for i in range(2, 301):
    fib_vals.append((fib_vals[-1] + fib_vals[-2]) % 7)
f4 = lambda n: fib_vals[int(n)] if int(n) < len(fib_vals) else 0

fc4 = build_function_crystal("F_n mod 7", f4, range(301))
sn4, rel4, slot4 = full_crystal_report(fc4, "F_n mod 7")

# Resonance search
res4 = find_resonant_frame(f4, range(1, 301), max_multiplier=16, name="F_n mod 7")
if res4.best:
    print(f"\n    Best resonant frame: k={res4.best.multiplier}, "
          f"Schmidt={res4.best.schmidt:.3f}, period={res4.best.period}")

# Pisano period
pisano = 16
print(f"    Pisano period π(7) = {pisano}")
print(f"    3 | 16? {16 % 3 == 0} — NO")
print(f"    2 | 16? {16 % 2 == 0} — YES → partial resonance possible")

# Crystal diagnosis
print(f"\n  ── Crystal diagnosis ──")
if sn4 < 2.2:
    print(f"    Crystal says: below noise (Schmidt={sn4:.3f})")
    print(f"    Recommendation: matrix exponentiation fallback")
else:
    print(f"    Crystal says: structure detected (Schmidt={sn4:.3f})")


# ═══════════════════════════════════════════════════════════════
#  Summary: the crystal as triage system
# ═══════════════════════════════════════════════════════════════

print("\n\n" + "=" * 70)
print("THE CRYSTAL AS TRIAGE SYSTEM")
print("=" * 70)

problems = [
    ("7^2025 mod 100", sn, rel, slot, "period 4, non-resonant"),
    ("digit_sum == 10", sn3, rel3, slot3, "non-periodic"),
    ("F_100 mod 7", sn4, rel4, slot4, "Pisano 16, partial"),
]

print(f"\n  {'Problem':>20s} {'Schmidt':>8s} {'Rel':>13s} {'Slot':>10s} {'Diagnosis':>30s}")
print("  " + "-" * 85)
for name, sn_val, rel_val, slot_val, diag in problems:
    if sn_val > 2.2:
        action = "CRYSTAL CAN HELP"
    elif sn_val > 1.5:
        action = "WEAK — verify"
    else:
        action = "FALLBACK to exact"
    print(f"  {name:>20s} {sn_val:>8.3f} {rel_val:>13s} {slot_val:>10s} {diag:>30s} → {action}")
