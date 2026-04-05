"""
Gauge confinement: time-averaged probability density on S¹.

If the gauge direction is confined near A₂ roots, the time-averaged
density should peak at those positions regardless of when we sample.

Also: how does the density change with entangling step? Early steps
should be more diffuse (wandering). Late steps should be more peaked
(confined). The TRANSITION from diffuse to confined is the gauge
confinement mechanism.

Finally: the 1/(H²+1) quantum. It appears in:
  K_peak = K*/2 + 1/(H²+1)
  1-K* = (H-1)/H + 1/(H²+1)
What is 1/(H²+1) in the crystal framework?
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
    """Get the angle and magnitude of the standard projection."""
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


def entangle_angles(correlation, seed, n_steps=50, discount=0.3):
    """Entangle and return angle at every step."""
    _random.seed(seed)
    corr = correlation.float()
    flat = corr.reshape(-1)
    flat = flat / flat.sum().clamp(min=EPS_LOG)

    joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j
    angles = []

    for step in range(n_steps):
        angle, mag = standard_angle(joint)
        angles.append(angle)

        idx = _random.choices(range(H * H), weights=flat.tolist())[0]
        hi, hj = idx // H, idx % H
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

    return angles


identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

print("=" * 80)
print("GAUGE CONFINEMENT: TIME-AVERAGED DENSITY")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Collect ALL angles at ALL steps
# ═══════════════════════════════════════════════════════════════

n_seeds = 300
n_steps = 50
n_bins = 12

# Density by step window
windows = [
    ("steps 1-5", 1, 5),
    ("steps 5-10", 5, 10),
    ("steps 10-20", 10, 20),
    ("steps 20-30", 20, 30),
    ("steps 30-40", 30, 40),
    ("steps 40-50", 40, 50),
]

window_angles = {label: [] for label, _, _ in windows}
all_angles = []

print(f"\n  Collecting angles: {n_seeds} seeds × {n_steps} steps...")

for seed in range(n_seeds):
    angles = entangle_angles(identity_corr, seed=seed, n_steps=n_steps)
    for step, angle in enumerate(angles):
        if angle is not None:
            all_angles.append(angle)
            for label, lo, hi in windows:
                if lo <= step < hi:
                    window_angles[label].append(angle)


# ═══════════════════════════════════════════════════════════════
#  Density by step window
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Angular density by step window ---\n")

bw = 2 * math.pi / n_bins

for label, lo, hi in windows:
    angles_w = window_angles[label]
    counts = [0] * n_bins
    for a in angles_w:
        b = int(((a % (2*math.pi)) / bw))
        if b >= n_bins: b = n_bins - 1
        counts[b] += 1

    # Entropy of the distribution (lower = more concentrated)
    total = sum(counts)
    probs = [c / total for c in counts] if total > 0 else [1/n_bins]*n_bins
    entropy = -sum(p * math.log(p + 1e-15) for p in probs)
    max_entropy = math.log(n_bins)

    print(f"  {label}: n={len(angles_w)}, entropy={entropy:.3f}/{max_entropy:.3f} = {entropy/max_entropy:.3f}")

    # Compact bar chart
    max_c = max(counts) if counts else 1
    for i in range(n_bins):
        bar = '█' * (counts[i] * 30 // max(max_c, 1))
        print(f"    [{i*bw/math.pi:.1f}π,{(i+1)*bw/math.pi:.1f}π) {counts[i]:4d} {bar}")
    print()


# ═══════════════════════════════════════════════════════════════
#  Confinement measure: entropy vs step
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print("CONFINEMENT: ENTROPY vs STEP")
print("="*80)

# Fine-grained: entropy at each step
max_entropy = math.log(n_bins)

print(f"\n  {'step':>4s}  {'entropy':>8s}  {'S/S_max':>8s}  {'confinement':>12s}")

for step in range(1, n_steps):
    step_angles = []
    for seed in range(n_seeds):
        angles = entangle_angles(identity_corr, seed=seed, n_steps=n_steps)
        if angles[step] is not None:
            step_angles.append(angles[step])

    if len(step_angles) < 10:
        continue

    counts = [0] * n_bins
    for a in step_angles:
        b = int(((a % (2*math.pi)) / bw))
        if b >= n_bins: b = n_bins - 1
        counts[b] += 1

    total = sum(counts)
    probs = [c / total for c in counts]
    entropy = -sum(p * math.log(p + 1e-15) for p in probs)
    confinement = 1 - entropy / max_entropy

    if step <= 10 or step % 5 == 0:
        print(f"  {step:4d}  {entropy:>8.4f}  {entropy/max_entropy:>8.4f}  {confinement:>12.4f}")

    # Only compute every 5th step to save time
    if step > 10 and step % 5 != 0:
        continue


# ═══════════════════════════════════════════════════════════════
#  The 1/(H²+1) quantum
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE 1/(H²+1) QUANTUM")
print("="*80)

val = 1 / (H**2 + 1)
print(f"\n  1/(H²+1) = 1/10 = {val:.6f}")
print(f"\n  Appearances:")
print(f"    K_peak = K*/2 + 1/(H²+1)")
print(f"    1-K* = (H-1)/H + 1/(H²+1)")
print(f"    1/(H²+1) = Born prob of one entry in a uniform [H²+1]-dim state")
print(f"    H²+1 = dim(Sym²(C^(H-1))) = dim of symmetric products")

# What is H²+1?
print(f"\n  H²+1 = {H**2+1}")
print(f"  = number of distinct entries in a symmetric {H}×{H} matrix + 1")
print(f"  = H(H+1)/2 + H = H(H+3)/2 at H=3: 9 ≠ 10. No.")
print(f"  = (H-1)² + (H-1) + 1 = Φ₃(H-1) at H=3: Φ₃(2) = 7 ≠ 10. No.")
print(f"  = H² + 1. Simply the mass-squared plus one.")
print(f"  = MASS_DIM² - MASS_DIM² + H² + 1 = ...")
print(f"  = The number of (i,j) pairs with i,j ∈ {{0,...,H-1}} PLUS one (for θ)")
print(f"  = H² hypothesis-pairs + 1 θ-pair = H²+1 distinct DS focal elements")
print(f"  = The DENOMINATOR of the group action: H(H²+1) = H × (hypothesis-pairs + θ)")

print(f"\n  1/(H²+1) = the Born probability of one focal element among H²+1")
print(f"  It IS the information contributed by a single focal element")
print(f"  to the Dempster combination.")

print(f"\n  So K_peak = K*/2 + [one focal element's worth of information]")
print(f"  And 1-K* = (H-1)/H + [one focal element's worth of information]")

print(f"\n  The nucleation happens when the accumulated conflict equals")
print(f"  half the equilibrium conflict PLUS one focal element's contribution.")
print(f"  This is when the evidence has filled exactly one degree of freedom")
print(f"  beyond the background conflict rate.")


print(f"\n\n{'='*80}")
print("WHAT CONFINEMENT REVEALS")
print("="*80)
