"""
Corrected S₃ projector: complex-aware.

The old projector permuted M.real.float(), losing complex phase information.
The correct projector permutes the full complex matrix.

Redo the key measurements from Principles 59-60 with the corrected projector.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import random as _random
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, born_probabilities,
                            EPS_LOG, EPS_DIV, EPS_NORM)
from solver.crystals import compose, Entangler

torch.set_grad_enabled(False)


def s3_fractions(joint):
    """CORRECT S₃ projections on complex masses."""
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

    f_triv = triv.abs().pow(2).sum().item() / total
    f_sign = sign.abs().pow(2).sum().item() / total
    f_std = std.abs().pow(2).sum().item() / total
    return f_triv, f_sign, f_std


def entangle_with_irreps(correlation, seed=42, n_steps=70, discount=0.3):
    """Entangle and track irrep decomposition with correct projector."""
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
        f_triv, f_sign, f_std = s3_fractions(joint)

        # Evidence step
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
        })

    return joint, trace


identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
cycle_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)

print("=" * 80)
print("CORRECTED GAUGE FLOW (complex projector)")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Entangling trajectory: identity
# ═══════════════════════════════════════════════════════════════

n_seeds = 30
n_steps = 70

print(f"\n--- Identity crystal, {n_seeds} seeds, corrected projector ---\n")

all_traces = []
for seed in range(n_seeds):
    _, trace = entangle_with_irreps(identity_corr, seed=seed, n_steps=n_steps)
    all_traces.append(trace)

avg = []
for step in range(n_steps):
    avg.append({k: sum(all_traces[s][step][k] for s in range(n_seeds)) / n_seeds
                for k in all_traces[0][0].keys()})

print(f"  {'step':>4s}  {'K':>8s}  {'S':>8s}  {'Schmidt':>8s}  {'trivial':>8s}  {'sign':>8s}  {'standard':>8s}")
print("-" * 70)
for step in range(0, n_steps, 5):
    a = avg[step]
    print(f"  {step:4d}  {a['K']:>8.5f}  {a['entropy']:>8.4f}  {a['schmidt']:>8.3f}  "
          f"{a['trivial']:>8.4f}  {a['sign']:>8.4f}  {a['standard']:>8.4f}")


# ═══════════════════════════════════════════════════════════════
#  Cycle crystal (should show sign sector)
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Cycle crystal, {n_seeds} seeds, corrected projector ---\n")

cycle_traces = []
for seed in range(n_seeds):
    _, trace = entangle_with_irreps(cycle_corr, seed=seed, n_steps=n_steps)
    cycle_traces.append(trace)

cycle_avg = []
for step in range(n_steps):
    cycle_avg.append({k: sum(cycle_traces[s][step][k] for s in range(n_seeds)) / n_seeds
                      for k in cycle_traces[0][0].keys()})

print(f"  {'step':>4s}  {'K':>8s}  {'trivial':>8s}  {'sign':>8s}  {'standard':>8s}")
print("-" * 55)
for step in range(0, n_steps, 5):
    a = cycle_avg[step]
    print(f"  {step:4d}  {a['K']:>8.5f}  {a['trivial']:>8.4f}  {a['sign']:>8.4f}  {a['standard']:>8.4f}")


# ═══════════════════════════════════════════════════════════════
#  Composition cascade with correct projector
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("COMPOSITION CASCADE (corrected projector)")
print("="*80)

for name, corr in [("identity", identity_corr), ("cycle", cycle_corr)]:
    avg_joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(50):
        ent = Entangler(corr, seed=seed).build()
        avg_joint += ent.joint
    avg_joint /= 50

    print(f"\n  {name} (50-seed average):")
    print(f"  {'comp':>4s}  {'trivial':>8s}  {'sign':>8s}  {'standard':>8s}  {'Schmidt':>8s}")

    current = avg_joint.clone()
    for comp in range(15):
        f_t, f_s, f_st = s3_fractions(current)
        sn = schmidt_number(current)
        print(f"  {comp:4d}  {f_t:>8.5f}  {f_s:>8.5f}  {f_st:>8.5f}  {sn:>8.3f}")
        current = compose(current, avg_joint)


# ═══════════════════════════════════════════════════════════════
#  Hysteresis check with correct projector
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("HYSTERESIS WITH CORRECTED PROJECTOR")
print("="*80)

# Find peak and compare ascent vs descent
peak_steps = []
for s in range(n_seeds):
    peak = max(range(n_steps), key=lambda st: all_traces[s][st]['entropy'])
    peak_steps.append(peak)

ascent_data = {}
descent_data = {}
for s in range(n_seeds):
    peak = peak_steps[s]
    for i, t in enumerate(all_traces[s]):
        k_bin = round(t['K'] / 0.01) * 0.01
        entry = (t['entropy'], t['trivial'], t['standard'], t['schmidt'])
        if i < peak:
            ascent_data.setdefault(k_bin, []).append(entry)
        elif i > peak:
            descent_data.setdefault(k_bin, []).append(entry)

common = sorted(set(ascent_data.keys()) & set(descent_data.keys()))

print(f"\n  {'K':>6s}  {'S_asc':>8s}  {'S_des':>8s}  {'std_a':>8s}  {'std_d':>8s}  {'triv_a':>8s}  {'triv_d':>8s}")
print("-" * 65)
for kb in common[:15]:
    a_vals = ascent_data[kb]
    d_vals = descent_data[kb]
    a_avg = [sum(x[i] for x in a_vals)/len(a_vals) for i in range(4)]
    d_avg = [sum(x[i] for x in d_vals)/len(d_vals) for i in range(4)]
    print(f"  {kb:>6.3f}  {a_avg[0]:>8.4f}  {d_avg[0]:>8.4f}  "
          f"{a_avg[2]:>8.4f}  {d_avg[2]:>8.4f}  "
          f"{a_avg[1]:>8.4f}  {d_avg[1]:>8.4f}")


print(f"\n\n{'='*80}")
print("THE CORRECTED PICTURE")
print("="*80)
