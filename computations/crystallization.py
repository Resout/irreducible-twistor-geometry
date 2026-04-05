"""
Crystallization from the flow.

Start with total ignorance — the whole space, all possibility.
Add RANDOM conflict — not correlated evidence, just disagreement.
Does structure precipitate?

If reality is what falls out of conflict in the ethereal flow,
then the correlation doesn't matter — the CONFLICT matters.
Structure should crystallize from ANY conflict, not just structured evidence.

Test:
1. Random evidence (no correlation) — does Schmidt rise?
2. Pure conflict (maximally disagreeing evidence) — what crystallizes?
3. Compare to structured evidence (identity correlation)
4. What is the MINIMUM conflict needed to crystallize structure?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import random
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities,
                            EPS_LOG, EPS_DIV, EPS_NORM)

torch.set_grad_enabled(False)


def entangle_from_random(n_steps=50, discount=0.3, seed=42, mode="random"):
    """Entangle from total ignorance with different evidence types."""
    random.seed(seed)

    joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j  # total ignorance

    Ks = []

    for step in range(n_steps):
        # Generate evidence
        ev_a = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev_b = torch.zeros(MASS_DIM, dtype=torch.cfloat)

        if mode == "random":
            # Completely random: pick any (hi, hj) uniformly
            hi = random.randint(0, H-1)
            hj = random.randint(0, H-1)
        elif mode == "identity":
            # Structured: identity correlation (hi = hj)
            hi = random.randint(0, H-1)
            hj = hi
        elif mode == "pure_conflict":
            # Maximally conflicting: always disagree with current state
            # Pick the LEAST probable hypothesis
            bp = born_probabilities(joint.reshape(-1))
            bp_2d = bp.reshape(MASS_DIM, MASS_DIM)
            marg = bp_2d[:H, :H].sum(dim=1)
            hi = marg.argmin().item()
            hj = marg.argmin().item()
        elif mode == "anti_correlated":
            # Anti-correlated: hi ≠ hj always
            hi = random.randint(0, H-1)
            hj = (hi + 1 + random.randint(0, H-2)) % H
        elif mode == "constant":
            # Always the same evidence
            hi, hj = 0, 0

        ev_a[hi] = 0.55 + 0.06j
        ev_b[hj] = 0.55 + 0.06j
        for k in range(H):
            noise = random.gauss(0, 0.02)
            if k != hi: ev_a[k] = 0.05 + noise * 1j
            if k != hj: ev_b[k] = 0.05 + noise * 1j
        ev_a[H] = (1.0 - ev_a[:H].real.sum()) - ev_a[:H].imag.sum() * 1j
        ev_b[H] = (1.0 - ev_b[:H].real.sum()) - ev_b[:H].imag.sum() * 1j

        # DS combine
        joint_dim = MASS_DIM ** 2
        jH = joint_dim - 1
        joint_ev = torch.einsum("m,n->mn", ev_a, ev_b).reshape(joint_dim)
        ev_s = joint_ev[:jH] * discount
        ev_t = 1.0 + 0j - ev_s.sum()
        disc_ev = torch.cat([ev_s, ev_t.unsqueeze(0)])

        curr = joint.reshape(joint_dim)
        cs1, ct1 = curr[:jH], curr[jH:]
        cs2, ct2 = disc_ev[:jH], disc_ev[jH:]
        comb_s = cs1 * cs2 + cs1 * ct2 + ct1 * cs2
        comb_t = ct1 * ct2
        K = cs1.sum() * cs2.sum() - (cs1 * cs2).sum()
        norm = 1.0 - K
        safe_norm = norm if abs(norm) > EPS_NORM else torch.tensor(EPS_NORM + 0j)
        updated = torch.cat([comb_s, comb_t]) / safe_norm
        re_sum = updated.real.sum()
        if abs(re_sum) > EPS_DIV:
            updated = updated / re_sum
        joint = updated.reshape(MASS_DIM, MASS_DIM)

        Ks.append(K.abs().item() if K.is_complex() else K.item())

    return joint, Ks


print("=" * 80)
print("CRYSTALLIZATION FROM THE FLOW")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Compare: random vs structured vs anti-correlated evidence
# ═══════════════════════════════════════════════════════════════

modes = ["random", "identity", "anti_correlated", "constant", "pure_conflict"]

print(f"\n--- What crystallizes from different conflict types? ---\n")
print(f"  {'mode':>16s} {'Schmidt':>8s} {'Born(θ)':>8s} {'mean K':>8s} {'K_final':>8s}")
print("-" * 55)

n_seeds = 30
for mode in modes:
    schmidts, borns, mean_Ks, final_Ks = [], [], [], []
    for seed in range(n_seeds):
        joint, Ks = entangle_from_random(n_steps=50, seed=seed, mode=mode)
        sn = schmidt_number(joint)
        bp = born_probabilities(joint).reshape(MASS_DIM, MASS_DIM)
        marg = bp.sum(dim=1)
        bt = marg[H].item()
        schmidts.append(sn)
        borns.append(bt)
        mean_Ks.append(sum(Ks[-10:]) / 10)
        final_Ks.append(Ks[-1])

    mean_sn = sum(schmidts) / len(schmidts)
    mean_bt = sum(borns) / len(borns)
    mean_k = sum(mean_Ks) / len(mean_Ks)
    mean_fk = sum(final_Ks) / len(final_Ks)

    print(f"  {mode:>16s} {mean_sn:>8.3f} {mean_bt:>8.4f} {mean_k:>8.4f} {mean_fk:>8.4f}")

print(f"\n  K* = {K_STAR:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Random crystallization: what TYPE of structure precipitates?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("WHAT TYPE OF CRYSTAL PRECIPITATES FROM RANDOM CONFLICT?")
print("="*80)

from solver.crystals import classify_relationship, slot_measure

for mode in ["random", "identity", "anti_correlated"]:
    types = {}
    for seed in range(50):
        joint, _ = entangle_from_random(n_steps=50, seed=seed, mode=mode)
        rel, conf = classify_relationship(joint)
        types[rel] = types.get(rel, 0) + 1

    print(f"\n  {mode}:")
    for rel, count in sorted(types.items(), key=lambda x: -x[1]):
        print(f"    {rel:>14s}: {count}/50 ({count/50*100:.0f}%)")


# ═══════════════════════════════════════════════════════════════
#  The conflict K trajectory for random evidence
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("K TRAJECTORY: DOES RANDOM CONFLICT REACH K*?")
print("="*80)

# Average K over seeds at each step
n_seeds_trace = 30
all_Ks = {mode: [] for mode in ["random", "identity", "anti_correlated"]}

for mode in all_Ks:
    step_Ks = [[] for _ in range(50)]
    for seed in range(n_seeds_trace):
        _, Ks = entangle_from_random(n_steps=50, seed=seed, mode=mode)
        for step, k in enumerate(Ks):
            step_Ks[step].append(k)
    all_Ks[mode] = [sum(ks)/len(ks) for ks in step_Ks]

print(f"\n  Average K at selected steps:")
print(f"  {'step':>4s}", end="")
for mode in all_Ks:
    print(f"  {mode:>14s}", end="")
print()

for step in [0, 5, 10, 15, 20, 30, 40, 49]:
    print(f"  {step:4d}", end="")
    for mode in all_Ks:
        print(f"  {all_Ks[mode][step]:>14.5f}", end="")
    print()

print(f"\n  K* = {K_STAR:.5f}")


# ═══════════════════════════════════════════════════════════════
#  Minimum conflict for crystallization
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("MINIMUM CONFLICT FOR CRYSTALLIZATION")
print("="*80)

# Vary the discount factor (which controls evidence strength)
# At what discount does structure first appear?

from solver.algebra import SCHMIDT_THRESHOLD

print(f"\n  Schmidt threshold: {SCHMIDT_THRESHOLD:.4f}")
print(f"\n  {'discount':>8s} {'Schmidt':>8s} {'Born(θ)':>8s} {'mean K':>8s} {'crystallized':>12s}")
print("-" * 55)

for d in [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]:
    schmidts = []
    for seed in range(20):
        joint, Ks = entangle_from_random(n_steps=50, seed=seed, mode="identity", discount=d)
        schmidts.append(schmidt_number(joint))

    mean_sn = sum(schmidts) / len(schmidts)
    joint_last, Ks_last = entangle_from_random(n_steps=50, seed=0, mode="identity", discount=d)
    bt = born_probabilities(joint_last).reshape(MASS_DIM, MASS_DIM).sum(dim=1)[H].item()
    mk = sum(Ks_last[-10:]) / 10

    crystallized = "YES" if mean_sn > SCHMIDT_THRESHOLD else "no"
    print(f"  {d:>8.3f} {mean_sn:>8.3f} {bt:>8.4f} {mk:>8.4f} {crystallized:>12s}")


print(f"\n\n{'='*80}")
print("WHAT CRYSTALLIZATION SHOWS")
print("="*80)
