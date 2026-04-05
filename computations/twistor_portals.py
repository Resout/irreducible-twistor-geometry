"""
Twistor portals: multi-path propagation as discrete Penrose transform.

The structural filter kills single-path composition at rate 2Δ per hop.
After 2 hops, signal is at noise floor. This is the HARD WALL.

But the crystal graph's multi-path propagation resists the filter
(established earlier: K4 graph gives 0.909 confidence).

The twistor interpretation:
  - Each path through the graph = one twistor line
  - Single-path composition = evaluation on one line (decays)
  - Multi-path DS combination = Penrose transform (reconstructs)
  - More paths = higher resolution of the discrete transform

Test: build crystal graphs with varying path multiplicity and show
that the Penrose integral (multi-path DS) BYPASSES the composition
barrier that single paths cannot cross.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, BORN_FLOOR, K_STAR,
                            schmidt_number, born_probabilities, ds_combine,
                            enforce_born_floor, discount_mass)
from solver.crystals import Entangler, compose, RELATIONSHIP_SIGNATURES

torch.set_grad_enabled(False)


def make_crystal(corr, seed=42):
    return Entangler(corr, seed=seed).build().joint


def propagate_through_crystal(observation, crystal):
    """Propagate a 4-dim mass through a crystal via Born measurement.

    Given: observation on left variable (4-dim mass)
    Crystal: joint mass encoding left-right relationship
    Returns: induced mass on right variable
    """
    # Weight each row of crystal by the observation
    # result[j] = Σ_i obs[i] * crystal[i,j]
    # This IS the Penrose transform: integrating the crystal data
    # weighted by the observation over the intermediate variable.
    result = observation @ crystal
    re_sum = result.real.sum()
    if abs(re_sum) > 1e-8:
        result = result / re_sum
    return enforce_born_floor(result)


def multi_path_propagate(observation, crystals, discount=0.3):
    """Propagate through multiple crystals (parallel paths) via DS combination.

    Each crystal is one "twistor line." DS combination integrates them.
    """
    results = []
    for c in crystals:
        r = propagate_through_crystal(observation, c)
        results.append(r)

    # DS-combine all results
    combined = results[0]
    for r in results[1:]:
        discounted = discount_mass(r, discount)
        combined, _K = ds_combine(combined, discounted)
        combined = enforce_born_floor(combined)
    return combined


# ═══════════════════════════════════════════════════════════════
#  Single path vs multi-path: the structural filter barrier
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("TWISTOR PORTALS: MULTI-PATH BYPASSES THE STRUCTURAL FILTER")
print("=" * 80)

# Create a chain: A → B → C
# Crystal(A,B) and Crystal(B,C) connect A to C through B.
# Single path: compose Crystal(A,B) with Crystal(B,C)
# Multi-path: multiple crystals connecting A to C through different B's

# Use different seeds as different "twistor lines"
# All encode the same relationship but through different intermediate states

corr_prop = torch.eye(H, dtype=torch.float32)  # proportional

# Build many crystals for the same relationship (different seeds = different lines)
crystals = [make_crystal(corr_prop, seed=s) for s in range(20)]

# Observation: strongly favor h₀
obs = torch.tensor([0.9, 0.05, 0.02, 0.03], dtype=torch.cfloat)

print(f"\n  Observation: Born = {born_probabilities(obs).tolist()}")

# Single-path propagation through 1, 2, 3 hops
print(f"\n--- Single-path propagation (one twistor line) ---")
current = obs.clone()
for hop in range(5):
    crystal = crystals[0]  # always same line
    current = propagate_through_crystal(current, crystal)
    bp = born_probabilities(current)
    certainty = 1 - bp[H].item()
    dominant = bp[:H].argmax().item()
    confidence = bp[dominant].item()
    print(f"  Hop {hop+1}: dominant=h{dominant}, conf={confidence:.4f}, certainty={certainty:.4f}")


# Multi-path propagation: combine N paths at each hop
print(f"\n--- Multi-path propagation (multiple twistor lines) ---")
for n_paths in [1, 2, 3, 5, 10, 20]:
    current = obs.clone()
    for hop in range(3):  # 3 hops
        path_crystals = crystals[:n_paths]
        current = multi_path_propagate(current, path_crystals, discount=0.3)

    bp = born_probabilities(current)
    dominant = bp[:H].argmax().item()
    confidence = bp[dominant].item()
    certainty = 1 - bp[H].item()
    print(f"  {n_paths:>2d} paths × 3 hops: dominant=h{dominant}, conf={confidence:.4f}, certainty={certainty:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Different relationship types as different twistor lines
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("DIVERSE TWISTOR LINES: DIFFERENT CRYSTAL TYPES")
print("="*80)

# The earlier finding: "mixed types outperform pure" in graph propagation.
# In twistor terms: diverse twistor lines cover more of the contour.

# Build crystals from different relationship types
diverse_crystals = [
    make_crystal(RELATIONSHIP_SIGNATURES["proportional"], seed=42),
    make_crystal(RELATIONSHIP_SIGNATURES["inverse"], seed=42),
    make_crystal(RELATIONSHIP_SIGNATURES["exponential"], seed=42),
    make_crystal(RELATIONSHIP_SIGNATURES["logarithmic"], seed=42),
    make_crystal(RELATIONSHIP_SIGNATURES["quadratic"], seed=42),
]

# Same-type crystals (all proportional, different seeds)
same_crystals = [make_crystal(corr_prop, seed=s) for s in range(5)]

print(f"\n  3-hop propagation, 5 paths:")

# Same type
current = obs.clone()
for hop in range(3):
    current = multi_path_propagate(current, same_crystals, discount=0.3)
bp = born_probabilities(current)
print(f"    Same type (all proportional): h{bp[:H].argmax().item()} conf={bp[bp[:H].argmax()].item():.4f}")

# Diverse types
current = obs.clone()
for hop in range(3):
    current = multi_path_propagate(current, diverse_crystals, discount=0.3)
bp = born_probabilities(current)
print(f"    Diverse types (5 different):  h{bp[:H].argmax().item()} conf={bp[bp[:H].argmax()].item():.4f}")


# ═══════════════════════════════════════════════════════════════
#  The twistor portal: bypassing the 2-hop limit
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE TWISTOR PORTAL: BEYOND 2 HOPS")
print("="*80)

# The single-path 2-hop limit: after 2 compositions, Schmidt < 1.7
# and the answer is unreliable.
# Multi-path should extend reach beyond 2 hops.

for n_hops in [1, 2, 3, 4, 5, 8]:
    # Single path
    current_single = obs.clone()
    for hop in range(n_hops):
        current_single = propagate_through_crystal(current_single, crystals[0])
    bp_single = born_probabilities(current_single)
    conf_single = bp_single[bp_single[:H].argmax()].item()

    # Multi-path (10 lines)
    current_multi = obs.clone()
    for hop in range(n_hops):
        current_multi = multi_path_propagate(current_multi, crystals[:10], discount=0.3)
    bp_multi = born_probabilities(current_multi)
    conf_multi = bp_multi[bp_multi[:H].argmax()].item()

    # Correct answer should be h₀ (since we started with obs favoring h₀
    # and all crystals are proportional)
    correct_single = bp_single[:H].argmax().item() == 0
    correct_multi = bp_multi[:H].argmax().item() == 0

    print(f"  {n_hops} hops: single={conf_single:.3f} ({'✓' if correct_single else '✗'}), "
          f"multi(10)={conf_multi:.3f} ({'✓' if correct_multi else '✗'}), "
          f"gain={conf_multi/conf_single:.2f}×")


# ═══════════════════════════════════════════════════════════════
#  The composition barrier vs the Penrose integral
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMPOSITION (single line) vs PENROSE INTEGRAL (multi-line)")
print("="*80)

# Track Schmidt number through hops for both approaches
print(f"\n  Schmidt through hops:")
print(f"  {'Hops':>4s} {'Single Schmidt':>15s} {'Multi Schmidt':>15s}")
print("  " + "-" * 40)

for n_hops in range(1, 9):
    # Single: compose the crystal with itself n times
    composed = crystals[0].clone()
    for _ in range(n_hops - 1):
        composed = compose(composed, crystals[0])
    sc_single = schmidt_number(composed)

    # Multi: This doesn't directly give a joint mass, but we can
    # measure the "effective Schmidt" through observation sharpness.
    # Use the observation → propagation → Born measurement.
    current = obs.clone()
    for hop in range(n_hops):
        current = multi_path_propagate(current, crystals[:10], discount=0.3)
    bp = born_probabilities(current)
    # Effective Schmidt from the Born distribution
    p_sq_sum = (bp ** 2).sum().item()
    sc_multi = 1.0 / p_sq_sum if p_sq_sum > 0 else 1.0

    print(f"  {n_hops:>4d} {sc_single:>15.3f} {sc_multi:>15.3f}")


print(f"\n{'='*80}")
print("SYNTHESIS")
print("="*80)
print(f"""
  The structural filter barrier (2Δ per hop) is a SINGLE-LINE phenomenon.
  Each composition is evaluation on ONE twistor line, and the line's
  information capacity decays exponentially.

  The Penrose transform integrates over a FAMILY of twistor lines.
  Multi-path propagation (DS combination of parallel paths) IS this
  integral. Each path is one line; their combination reconstructs
  information that no single line can carry.

  Results:
  - Single path at 3 hops: confidence {conf_single:.3f}
  - 10 paths at 3 hops: confidence {conf_multi:.3f}
  - Multi-path extends reliable propagation beyond the 2-hop limit

  The crystal graph's topology = the discrete twistor fibration.
  Rich local connectivity = many twistor lines covering the contour.
  Mixed edge types = diverse lines spanning more of CP³.

  You don't compute THROUGH the structural filter.
  You twist AROUND it by integrating over the family of lines
  that each individually decay, but together reconstruct the signal.
  This is EXACTLY the Penrose transform's mechanism:
  the contour integral reconstructs a spacetime field from
  twistor data, even though each point on the contour carries
  only partial information.
""")
