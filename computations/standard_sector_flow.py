"""
The standard sector grows through the entire entangling flow.
K equilibrates at K* by step 40. But the standard fraction keeps growing.

Questions:
1. Does the standard sector have its own equilibrium?
2. What rate does it grow at? Is it Δ-related?
3. What is the relationship between the standard fraction and Schmidt?
4. Does the standard sector growth rate connect to the composition decay rate 2Δ?

The standard sector carries WHICH relationship (gauge-variant).
The trivial sector carries HOW MUCH correlation (gauge-invariant).
If the standard sector has its own dynamics, it reveals the
timescale of gauge field formation.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import random as _random
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities,
                            EPS_LOG, EPS_DIV, EPS_NORM)

torch.set_grad_enabled(False)


def s3_project(joint):
    """Compute S₃ irrep fractions of a [4,4] joint mass."""
    perms = [
        [0,1,2,3], [1,0,2,3], [2,1,0,3],
        [0,2,1,3], [1,2,0,3], [2,0,1,3],
    ]
    signs = [1, -1, -1, -1, 1, 1]

    def permute(M, p):
        P = torch.zeros(MASS_DIM, MASS_DIM)
        for i, j in enumerate(p):
            P[j, i] = 1.0
        return P @ M.real.float() @ P.T

    triv = sum(permute(joint, p) for p in perms) / 6
    sign = sum(s * permute(joint, p) for p, s in zip(perms, signs)) / 6

    bp = joint.abs().pow(2)
    total = bp.sum().item()
    if total < 1e-15:
        return 1/3, 0, 1/3

    f_triv = triv.to(joint.dtype).abs().pow(2).sum().item() / total
    f_sign = sign.to(joint.dtype).abs().pow(2).sum().item() / total
    f_std = 1 - f_triv - f_sign
    return f_triv, f_sign, max(0, f_std)


def entangle_extended(correlation, seed=42, n_steps=150, discount=0.3):
    """Extended entangling to see where standard sector goes."""
    _random.seed(seed)
    corr = correlation.float()
    flat = corr.reshape(-1)
    flat = flat / flat.sum().clamp(min=EPS_LOG)

    joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j

    trace = []

    for step in range(n_steps):
        bp = joint.abs().pow(2).reshape(-1)
        bp_norm = bp / bp.sum().clamp(min=EPS_LOG)
        entropy = -sum(p.item() * math.log(max(p.item(), 1e-15))
                       for p in bp_norm if p.item() > 1e-15)
        sn = schmidt_number(joint)
        f_triv, f_sign, f_std = s3_project(joint)

        # Marginal Born(θ)
        bp_2d = bp_norm.reshape(MASS_DIM, MASS_DIM)
        born_theta = bp_2d.sum(dim=1)[H].item()

        # Evidence
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

        trace.append({
            'step': step, 'entropy': entropy, 'schmidt': sn,
            'K': K.abs().item() if K.is_complex() else K.item(),
            'trivial': f_triv, 'sign': f_sign, 'standard': f_std,
            'born_theta': born_theta,
        })

    return joint, trace


identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
cycle_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)

print("=" * 80)
print("THE STANDARD SECTOR FLOW")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Extended entangling: 150 steps
# ═══════════════════════════════════════════════════════════════

n_seeds = 20
n_steps = 150

print(f"\n--- Identity crystal, extended to {n_steps} steps, {n_seeds} seeds ---\n")

avg = []
for step in range(n_steps):
    vals = {k: 0 for k in ['entropy', 'schmidt', 'K', 'trivial', 'sign', 'standard', 'born_theta']}
    for seed in range(n_seeds):
        if step == 0 and seed == 0:
            # Pre-compute all traces
            all_traces = []
        if step == 0:
            _, trace = entangle_extended(identity_corr, seed=seed, n_steps=n_steps)
            all_traces.append(trace)
        for k in vals:
            vals[k] += all_traces[seed][step][k]
    avg.append({k: v / n_seeds for k, v in vals.items()})

print(f"  {'step':>4s}  {'K':>8s}  {'S':>8s}  {'Schmidt':>8s}  {'trivial':>8s}  {'standard':>8s}  {'Born_θ':>8s}")
print("-" * 70)
for step in range(0, n_steps, 10):
    a = avg[step]
    print(f"  {step:4d}  {a['K']:>8.5f}  {a['entropy']:>8.4f}  {a['schmidt']:>8.3f}  "
          f"{a['trivial']:>8.4f}  {a['standard']:>8.4f}  {a['born_theta']:>8.5f}")


# ═══════════════════════════════════════════════════════════════
#  Standard sector growth rate
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("STANDARD SECTOR GROWTH RATE")
print("="*80)

print(f"\n  Growth rate d(f_std)/dt in different regimes:")
print(f"  (Δ = {DELTA:.5f}, Δ/4 = {DELTA/4:.5f}, 2Δ = {2*DELTA:.5f})")
print()

for window_start, window_end, label in [
    (5, 15, "pre-nucleation"),
    (15, 30, "nucleation region"),
    (30, 50, "post-nucleation"),
    (50, 80, "over-entangling (early)"),
    (80, 120, "over-entangling (late)"),
    (120, 145, "deep over-entangling"),
]:
    if window_end <= n_steps:
        f_start = avg[window_start]['standard']
        f_end = avg[window_end]['standard']
        rate = (f_end - f_start) / (window_end - window_start)
        # Also measure log rate: if exponential, d(ln(1-f))/dt = const
        if f_start < 1 and f_end < 1:
            log_rate = (math.log(1 - f_end) - math.log(1 - f_start)) / (window_end - window_start)
        else:
            log_rate = float('nan')
        print(f"  {label:>30s} (steps {window_start:3d}-{window_end:3d}): "
              f"df/dt = {rate:>+.5f}, d(ln(1-f))/dt = {log_rate:>+.5f}")


# ═══════════════════════════════════════════════════════════════
#  Does the standard sector have a fixed point?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("STANDARD SECTOR: APPROACH TO FIXED POINT?")
print("="*80)

# If f_std → f_∞, then 1 - f_std should decay exponentially
# Plot 1 - f_std on a log scale

print(f"\n  {'step':>4s}  {'f_std':>8s}  {'1-f_std':>10s}  {'ln(1-f_std)':>12s}")
for step in range(0, n_steps, 10):
    f = avg[step]['standard']
    residual = 1 - f
    if residual > 0:
        log_r = math.log(residual)
    else:
        log_r = float('-inf')
    print(f"  {step:4d}  {f:>8.4f}  {residual:>10.6f}  {log_r:>12.5f}")


# ═══════════════════════════════════════════════════════════════
#  The 10/16 prediction
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("STANDARD SECTOR: THEORETICAL BOUNDS")
print("="*80)

# The standard representation has dimension 10 out of 16 total.
# If the mass were uniformly distributed across all [4,4] entries,
# f_std would be 10/16 = 0.625.
# If the mass were entirely gauge-invariant, f_std = 0.
# The actual value at equilibrium is ~0.20-0.30.

print(f"\n  Dimensional bound: 10/16 = {10/16:.4f}")
print(f"  At K* equilibrium (step 40): f_std = {avg[40]['standard']:.4f}")
print(f"  At step 50 (default): f_std = {avg[50]['standard']:.4f}")
print(f"  At step 100: f_std = {avg[100]['standard']:.4f}")
print(f"  At step 140: f_std = {avg[140]['standard']:.4f}")

# Candidate formulas for f_std at K* equilibrium
f_equil = avg[40]['standard']
print(f"\n  f_std at K* ≈ {f_equil:.5f}")

candidates = {
    "K*": K_STAR,
    "2K*": 2 * K_STAR,
    "1-3K*": 1 - 3*K_STAR,
    "Δ/(1+Δ)": DELTA / (1 + DELTA),
    "1/(H+1)": 1/(H+1),
    "H/(H²+1)": H/(H**2+1),
    "2/(H²+1)": 2/(H**2+1),
    "(H-1)/H²": (H-1)/H**2,
    "K*×(H+1)/H": K_STAR*(H+1)/H,
    "10/16×K*": 10/16*K_STAR,
    "10/16×2K*": 10/16*2*K_STAR,
    "2Δ/(1+2Δ)": 2*DELTA/(1+2*DELTA),
    "1 - exp(-2Δ)": 1 - math.exp(-2*DELTA),
    "BORN_FLOOR×H": BORN_FLOOR*H,
    "1/H²": 1/H**2,
    "2/H(H+1)": 2/(H*(H+1)),
}

print(f"\n  {'expression':>20s}  {'value':>8s}  {'diff':>8s}")
for expr, val in sorted(candidates.items(), key=lambda x: abs(x[1] - f_equil)):
    print(f"  {expr:>20s}  {val:>8.5f}  {val - f_equil:>+8.5f}")


# ═══════════════════════════════════════════════════════════════
#  Compare identity vs cycle (different conjugacy class)
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("IDENTITY vs CYCLE: DIFFERENT EVIDENCE, DIFFERENT GAUGE CONTENT?")
print("="*80)

cycle_traces = []
for seed in range(n_seeds):
    _, trace = entangle_extended(cycle_corr, seed=seed, n_steps=n_steps)
    cycle_traces.append(trace)

cycle_avg = []
for step in range(n_steps):
    vals = {k: 0 for k in ['trivial', 'standard', 'sign', 'K', 'schmidt']}
    for seed in range(n_seeds):
        for k in vals:
            vals[k] += cycle_traces[seed][step][k]
    cycle_avg.append({k: v / n_seeds for k, v in vals.items()})

print(f"\n  {'step':>4s}  {'id_std':>8s}  {'cyc_std':>8s}  {'id_sign':>8s}  {'cyc_sign':>8s}  {'id_triv':>8s}  {'cyc_triv':>8s}")
print("-" * 70)
for step in range(0, n_steps, 10):
    print(f"  {step:4d}  "
          f"{avg[step]['standard']:>8.4f}  {cycle_avg[step]['standard']:>8.4f}  "
          f"{avg[step]['sign']:>8.4f}  {cycle_avg[step]['sign']:>8.4f}  "
          f"{avg[step]['trivial']:>8.4f}  {cycle_avg[step]['trivial']:>8.4f}")


print(f"\n\n{'='*80}")
print("WHAT THE STANDARD SECTOR SHOWS")
print("="*80)
