"""
When does the gauge direction lock in?

Hypothesis: the first evidence sample determines the root direction.
If the first sample is (h₀,h₀), the crystal aligns with the root
vector associated with h₀. This should happen within the first few
steps, then be irreversible (Principle 14: the entangler must commit).

Track the gauge direction step by step and measure:
1. At what step does the direction stabilize?
2. Does the first evidence sample predict the final direction?
3. What's the angular velocity of the direction during entangling?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import random as _random
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, EPS_LOG, EPS_DIV, EPS_NORM)

torch.set_grad_enabled(False)


def standard_angle(joint):
    """Get the angle of the standard projection."""
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]
    signs_list = [1, -1, -1, -1, 1, 1]

    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM, dtype=M.dtype)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M @ P.T

    triv = sum(permute(joint, p) for p in perms) / 6
    sign = sum(s * permute(joint, p) for p, s in zip(perms, signs_list)) / 6
    std = joint - triv - sign

    U, S, Vh = torch.linalg.svd(std)
    if S[0].abs() < 1e-15:
        return None, 0.0

    v = U[:, 0][:H]
    e1 = torch.tensor([1, -1, 0], dtype=torch.cfloat) / math.sqrt(2)
    e2 = torch.tensor([1, 1, -2], dtype=torch.cfloat) / math.sqrt(6)
    c1 = torch.dot(v.conj(), e1).real.item()
    c2 = torch.dot(v.conj(), e2).real.item()

    angle = math.atan2(c2, c1)
    mag = math.sqrt(c1**2 + c2**2)
    return angle, mag


def entangle_step_by_step(correlation, seed=42, n_steps=50, discount=0.3):
    """Entangle one step at a time, recording gauge direction."""
    _random.seed(seed)
    corr = correlation.float()
    flat = corr.reshape(-1)
    flat = flat / flat.sum().clamp(min=EPS_LOG)

    joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j

    trace = []

    for step in range(n_steps):
        # Measure before evidence
        angle, mag = standard_angle(joint)

        # Sample evidence
        idx = _random.choices(range(H * H), weights=flat.tolist())[0]
        hi = idx // H
        hj = idx % H

        trace.append({
            'step': step,
            'angle': angle,
            'magnitude': mag,
            'evidence_pair': (hi, hj),
        })

        # Apply evidence
        ev_a = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev_b = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev_a[hi] = 0.55 + 0.06j
        ev_b[hj] = 0.55 + 0.06j
        for k in range(H):
            noise = _random.gauss(0, 0.02)
            if k != hi: ev_a[k] = 0.05 + noise * 1j
            if k != hj: ev_b[k] = 0.05 + noise * 1j
        ev_a[H] = (1.0 - ev_a[:H].real.sum()) - ev_a[:H].imag.sum() * 1j
        ev_b[H] = (1.0 - ev_b[:H].real.sum()) - ev_b[:H].imag.sum() * 1j

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

    return trace


identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

print("=" * 80)
print("GAUGE DIRECTION LOCKING")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Step-by-step direction for a few seeds
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Step-by-step gauge direction ---\n")

for seed in [0, 1, 2, 42, 99]:
    trace = entangle_step_by_step(identity_corr, seed=seed, n_steps=50)

    print(f"\n  Seed {seed}:")
    print(f"  {'step':>4s}  {'angle/π':>8s}  {'mag':>8s}  {'evidence':>10s}  {'Δangle/π':>10s}")

    prev_angle = None
    for t in trace[:20]:  # first 20 steps
        angle = t['angle']
        mag = t['magnitude']
        pair = t['evidence_pair']

        if angle is not None and prev_angle is not None:
            da = (angle - prev_angle) / math.pi
            # Wrap to [-1, 1]
            while da > 1: da -= 2
            while da < -1: da += 2
            da_str = f"{da:+.4f}"
        else:
            da_str = "---"

        a_str = f"{angle/math.pi:.4f}" if angle is not None else "---"

        print(f"  {t['step']:4d}  {a_str:>8s}  {mag:>8.5f}  ({pair[0]},{pair[1]}){' ':>5s}  {da_str:>10s}")
        prev_angle = angle


# ═══════════════════════════════════════════════════════════════
#  First evidence → final direction correlation
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("FIRST EVIDENCE SAMPLE → FINAL DIRECTION")
print("="*80)

n_seeds = 500
first_evidence = {0: [], 1: [], 2: []}  # h_i → list of final angles

for seed in range(n_seeds):
    trace = entangle_step_by_step(identity_corr, seed=seed, n_steps=50)

    # First evidence pair
    first_h = trace[0]['evidence_pair'][0]  # which h was first

    # Final angle
    final_angle = trace[-1]['angle']
    if final_angle is not None:
        first_evidence[first_h].append(final_angle)

print(f"\n  First evidence sample → final angle distribution:")
for h, angles in first_evidence.items():
    if angles:
        # Circular mean
        mean_sin = sum(math.sin(a) for a in angles) / len(angles)
        mean_cos = sum(math.cos(a) for a in angles) / len(angles)
        mean_angle = math.atan2(mean_sin, mean_cos)
        # Circular spread
        R = math.sqrt(mean_sin**2 + mean_cos**2)

        print(f"\n  First h={h}: n={len(angles)}, mean angle = {mean_angle/math.pi:.4f}π, concentration R = {R:.4f}")

        # Histogram
        n_bins = 6
        bw = 2 * math.pi / n_bins
        counts = [0] * n_bins
        for a in angles:
            b = int(((a % (2*math.pi)) / bw))
            if b >= n_bins: b = n_bins - 1
            counts[b] += 1
        for i in range(n_bins):
            bar = '█' * (counts[i] * 30 // max(max(counts), 1))
            print(f"    [{i*bw/math.pi:.2f}π,{(i+1)*bw/math.pi:.2f}π): {counts[i]:4d} {bar}")


# ═══════════════════════════════════════════════════════════════
#  Expected root angles
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("EXPECTED A₂ ROOT ANGLES")
print("="*80)

# The A₂ roots in the standard basis (e₁, e₂):
# α₁ corresponds to h₀-h₁ difference → e₁ direction → angle 0
# α₂ corresponds to h₁-h₂ difference → angle 2π/3
# α₃ corresponds to h₂-h₀ difference → angle 4π/3

# But the entangler doesn't just create a difference — it creates
# a complex joint mass. The projection onto the standard rep depends
# on the full [4,4] structure.

# The evidence (h₀,h₀) creates mass at joint[0,0]. This breaks symmetry
# in the h₀ direction. The standard projection of the h₀ direction is:
# h₀ = (1,0,0) → standard coords: (1/√2, 1/√6) → angle atan2(1/√6, 1/√2)

e1_h0 = 1 / math.sqrt(2)
e2_h0 = 1 / math.sqrt(6)
angle_h0 = math.atan2(e2_h0, e1_h0)

e1_h1 = -1 / math.sqrt(2)
e2_h1 = 1 / math.sqrt(6)
angle_h1 = math.atan2(e2_h1, e1_h1)

e1_h2 = 0
e2_h2 = -2 / math.sqrt(6)
angle_h2 = math.atan2(e2_h2, e1_h2)

print(f"\n  Standard rep projection of each h direction:")
print(f"  h₀ → angle = {angle_h0/math.pi:.4f}π (= {angle_h0:.4f} rad)")
print(f"  h₁ → angle = {angle_h1/math.pi:.4f}π (= {angle_h1:.4f} rad)")
print(f"  h₂ → angle = {angle_h2/math.pi:.4f}π (= {angle_h2:.4f} rad)")
print(f"\n  Separation: {(angle_h1-angle_h0)/math.pi:.4f}π, {(angle_h2-angle_h1)/math.pi:.4f}π")
print(f"  Expected: 2/3π = 120° intervals")


# ═══════════════════════════════════════════════════════════════
#  Locking time: when does angle stabilize?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("LOCKING TIME: WHEN DOES DIRECTION STABILIZE?")
print("="*80)

locking_steps = []
for seed in range(200):
    trace = entangle_step_by_step(identity_corr, seed=seed, n_steps=50)

    # Find when angle stays within 0.1 rad of final for 10 consecutive steps
    final_angle = trace[-1]['angle']
    if final_angle is None:
        continue

    locked = None
    for start in range(len(trace) - 10):
        all_close = True
        for s in range(start, start + 10):
            if trace[s]['angle'] is None:
                all_close = False
                break
            da = abs(trace[s]['angle'] - final_angle)
            da = min(da, 2*math.pi - da)
            if da > 0.1:
                all_close = False
                break
        if all_close:
            locked = start
            break

    if locked is not None:
        locking_steps.append(locked)

if locking_steps:
    avg_lock = sum(locking_steps) / len(locking_steps)
    sorted_lock = sorted(locking_steps)
    print(f"\n  {len(locking_steps)}/{200} seeds locked (within 0.1 rad of final for 10 steps)")
    print(f"  Average locking step: {avg_lock:.1f}")
    print(f"  Median: {sorted_lock[len(sorted_lock)//2]}")
    print(f"  10th percentile: {sorted_lock[len(sorted_lock)//10]}")
    print(f"  90th percentile: {sorted_lock[9*len(sorted_lock)//10]}")


print(f"\n\n{'='*80}")
print("WHAT THE LOCKING REVEALS")
print("="*80)
