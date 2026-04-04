"""
Can twistor portals amplify weak signals above the noise floor?

The crystal's noise floor: Schmidt ≈ 2.0 from random data.
Period-3 functions: Schmidt ≈ 3.27 (well above noise).
Period-5 functions: Schmidt ≈ 1.93 (AT noise floor).

Question: does multi-path propagation differentiate signal from noise
when the SINGLE crystal can't?

Test:
1. Build crystals from a period-5 function (weak signal, Schmidt ≈ 1.9)
2. Build crystals from random data (noise, Schmidt ≈ 2.0)
3. Propagate both through multi-path
4. Does multi-path differentiate signal from noise?

The key insight from Principle 40: multi-path saturates the Born ceiling
for CONSISTENT signals. If 10 seeds of a period-5 crystal all agree on
the same structure (even weakly), their DS combination should reinforce
the signal while noise (inconsistent across seeds) should cancel.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import random
from solver.algebra import (H, MASS_DIM, DELTA, BORN_FLOOR,
                            schmidt_number, born_probabilities,
                            ds_combine, enforce_born_floor, discount_mass)
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)


def correlation_from_function(f, domain_size=300):
    """Build a [3,3] correlation matrix from a function f: int → int."""
    # Bin inputs and outputs into 3 bins
    inputs = list(range(domain_size))
    outputs = [f(x) for x in inputs]

    # Bin by terciles
    sorted_in = sorted(inputs)
    sorted_out = sorted(outputs)
    in_cuts = [sorted_in[domain_size//3], sorted_in[2*domain_size//3]]
    out_cuts = [sorted_out[domain_size//3], sorted_out[2*domain_size//3]]

    def bin_val(v, cuts):
        if v < cuts[0]: return 0
        if v < cuts[1]: return 1
        return 2

    corr = torch.zeros(H, H, dtype=torch.float32)
    for x, y in zip(inputs, outputs):
        bi = bin_val(x, in_cuts)
        bo = bin_val(y, out_cuts)
        corr[bi, bo] += 1

    corr = corr / corr.sum().clamp(min=1e-10)
    return corr


def multi_path_propagate(observation, crystals, discount=0.3):
    """Propagate through multiple crystals via DS combination."""
    results = []
    for c in crystals:
        result = observation @ c
        re_sum = result.real.sum()
        if abs(re_sum) > 1e-8:
            result = result / re_sum
        result = enforce_born_floor(result)
        results.append(result)

    combined = results[0]
    for r in results[1:]:
        discounted = discount_mass(r, discount)
        combined, _K = ds_combine(combined, discounted)
        combined = enforce_born_floor(combined)
    return combined


# ═══════════════════════════════════════════════════════════════
#  Build crystals from functions of different periods
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("TWISTOR SIGNAL AMPLIFICATION")
print("=" * 80)

# Period-3 function (strong signal)
def f_period3(x):
    return pow(2, x, 7)  # period 3

# Period-5 function (weak signal)
def f_period5(x):
    return pow(3, x, 11)  # period 5

# Period-4 function (invisible to H=3)
def f_period4(x):
    return pow(2, x, 5)  # period 4

# Random function (noise)
random.seed(42)
rand_table = {x: random.randint(0, 100) for x in range(300)}
def f_random(x):
    return rand_table.get(x % 300, 0)

functions = {
    "period-3 (2^x mod 7)": f_period3,
    "period-5 (3^x mod 11)": f_period5,
    "period-4 (2^x mod 5)": f_period4,
    "random": f_random,
}


# ═══════════════════════════════════════════════════════════════
#  Single crystal: Schmidt and confidence for each function
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Single crystal (seed 42) ---")
print(f"  {'Function':>25s} {'Schmidt':>8s} {'dom':>4s} {'conf':>6s}")
print("  " + "-" * 50)

for name, f in functions.items():
    corr = correlation_from_function(f)
    ent = Entangler(corr, seed=42).build()
    sc = schmidt_number(ent.joint)

    # Propagate a uniform observation
    obs = torch.tensor([0.25, 0.25, 0.25, 0.25], dtype=torch.cfloat)
    result = obs @ ent.joint
    result = result / result.real.sum()
    result = enforce_born_floor(result)
    bp = born_probabilities(result)
    dom = bp[:H].argmax().item()
    conf = bp[dom].item()

    print(f"  {name:>25s} {sc:>8.3f} h{dom:>3d} {conf:>6.3f}")


# ═══════════════════════════════════════════════════════════════
#  Multi-path: does signal separate from noise?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("MULTI-PATH AMPLIFICATION: SIGNAL vs NOISE")
print("="*80)

obs = torch.tensor([0.25, 0.25, 0.25, 0.25], dtype=torch.cfloat)

for n_paths in [1, 5, 10, 20, 50]:
    print(f"\n  {n_paths} paths:")
    print(f"    {'Function':>25s} {'dom':>4s} {'conf':>6s} {'certainty':>10s}")
    print("    " + "-" * 50)

    for name, f in functions.items():
        # Build n_paths crystals with different seeds
        corr = correlation_from_function(f)
        crystals = []
        for seed in range(n_paths):
            ent = Entangler(corr, seed=seed).build()
            crystals.append(ent.joint)

        # Multi-path propagation
        result = multi_path_propagate(obs, crystals, discount=0.3)
        bp = born_probabilities(result)
        dom = bp[:H].argmax().item()
        conf = bp[dom].item()
        cert = 1 - bp[H].item()

        print(f"    {name:>25s} h{dom:>3d} {conf:>6.3f} {cert:>10.3f}")


# ═══════════════════════════════════════════════════════════════
#  The key test: does multi-path SEPARATE period-5 from random?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("SEPARATION: WEAK SIGNAL vs NOISE")
print("="*80)

# For each, build 50 crystals and measure the SPREAD of Born(dominant)
for name, f in [("period-5 (signal)", f_period5), ("random (noise)", f_random)]:
    corr = correlation_from_function(f)

    # Single-crystal spread
    single_confs = []
    for seed in range(50):
        ent = Entangler(corr, seed=seed).build()
        result = obs @ ent.joint
        result = result / result.real.sum()
        result = enforce_born_floor(result)
        bp = born_probabilities(result)
        single_confs.append(bp[:H].max().item())

    avg_s = sum(single_confs) / len(single_confs)
    std_s = (sum((c - avg_s)**2 for c in single_confs) / len(single_confs)) ** 0.5

    # Multi-path (10 paths) spread: use different sets of 10 seeds
    multi_confs = []
    for trial in range(20):
        crystals = []
        for seed in range(trial*10, trial*10 + 10):
            ent = Entangler(corr, seed=seed).build()
            crystals.append(ent.joint)
        result = multi_path_propagate(obs, crystals, discount=0.3)
        bp = born_probabilities(result)
        multi_confs.append(bp[:H].max().item())

    avg_m = sum(multi_confs) / len(multi_confs)
    std_m = (sum((c - avg_m)**2 for c in multi_confs) / len(multi_confs)) ** 0.5

    print(f"\n  {name}:")
    print(f"    Single: conf = {avg_s:.4f} ± {std_s:.4f}")
    print(f"    Multi(10): conf = {avg_m:.4f} ± {std_m:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Does CONSISTENCY of Born across seeds distinguish signal from noise?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("CONSISTENCY: DO SIGNAL CRYSTALS AGREE ACROSS SEEDS?")
print("="*80)

for name, f in functions.items():
    corr = correlation_from_function(f)

    # For each seed, which h-bin is dominant?
    dominant_counts = [0, 0, 0]
    for seed in range(50):
        ent = Entangler(corr, seed=seed).build()
        result = obs @ ent.joint
        result = result / result.real.sum()
        result = enforce_born_floor(result)
        bp = born_probabilities(result)
        dom = bp[:H].argmax().item()
        dominant_counts[dom] += 1

    max_agreement = max(dominant_counts) / 50.0
    print(f"  {name:>25s}: bin counts = {dominant_counts}, max agreement = {max_agreement:.0%}")


print(f"\n{'='*80}")
print("FINDINGS")
print("="*80)
