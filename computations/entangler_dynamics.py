"""
The entangler dynamics: what happens at each DS step?

The entangler builds a [4,4] crystal through 50 DS combinations.
Each step samples evidence from the correlation and combines it.
The conflict K at each step measures the disagreement.

We know K converges to K* = 7/30 (the eternal).
We know the final crystal's dominant eigenvector has Born(θ) = 10/33.

Questions:
1. At each step, what is the dominant eigenvector's Born(θ)?
2. Does it converge to 10/33 during entangling?
3. What is the trajectory of K through the 50 steps?
4. What is the trajectory of Schmidt through the 50 steps?
5. When does the crystal "become itself"?
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


def entangle_with_trace(correlation, seed=42, n_steps=50, discount=0.3):
    """Run the entangler and record everything at each step."""
    _random.seed(seed)

    corr = correlation.float()
    flat = corr.reshape(-1)
    flat = flat / flat.sum().clamp(min=EPS_LOG)

    joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j

    trace = []

    for step in range(n_steps):
        # Sample evidence
        idx = _random.choices(range(H * H), weights=flat.tolist())[0]
        hi = idx // H
        hj = idx % H

        ev_a = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev_b = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev_a[hi] = 0.55 + 0.06j
        ev_b[hj] = 0.55 + 0.06j
        for k in range(H):
            noise = _random.gauss(0, 0.02)
            if k != hi:
                ev_a[k] = 0.05 + noise * 1j
            if k != hj:
                ev_b[k] = 0.05 + noise * 1j
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

        # Record trace
        sn = schmidt_number(joint)

        # Dominant eigenvector Born(θ)
        try:
            eigvals, eigvecs = torch.linalg.eig(joint)
            dom_idx = eigvals.abs().argmax()
            v0 = eigvecs[:, dom_idx]
            bp_v0 = v0.abs().pow(2) / v0.abs().pow(2).sum()
            born_theta_eig = bp_v0[H].item()
        except:
            born_theta_eig = 0

        # Marginal Born(θ)
        bp_joint = joint.abs().pow(2)
        bp_joint = bp_joint / bp_joint.sum()
        marg = bp_joint.sum(dim=1)
        born_theta_marg = marg[H].item()

        trace.append({
            'step': step,
            'K': K.abs().item() if K.is_complex() else K.item(),
            'schmidt': sn,
            'born_theta_eig': born_theta_eig,
            'born_theta_marg': born_theta_marg,
            'hi': hi, 'hj': hj,
        })

    return joint, trace


print("=" * 80)
print("ENTANGLER DYNAMICS — WATCHING THE CRYSTAL FORM")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Single trajectory (identity crystal, seed 42)
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Single trajectory (identity correlation, seed 42) ---\n")

identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
joint, trace = entangle_with_trace(identity_corr, seed=42)

print(f"  {'step':>4s} {'K':>8s} {'Schmidt':>8s} {'Born(θ)_eig':>12s} {'Born(θ)_marg':>13s} {'sample':>8s}")
print("-" * 60)
for t in trace:
    if t['step'] < 10 or t['step'] % 5 == 0:
        print(f"  {t['step']:4d} {t['K']:>8.4f} {t['schmidt']:>8.3f} "
              f"{t['born_theta_eig']:>12.5f} {t['born_theta_marg']:>13.5f} "
              f"({t['hi']},{t['hj']})")


# ═══════════════════════════════════════════════════════════════
#  Average over seeds: when does Born(θ)_eig converge?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Average Born(θ)_eig over 30 seeds at each step ---\n")

n_seeds = 30
all_traces = []
for seed in range(n_seeds):
    _, trace = entangle_with_trace(identity_corr, seed=seed)
    all_traces.append(trace)

print(f"  {'step':>4s} {'mean_K':>8s} {'mean_Sn':>8s} {'mean_Born(θ)_eig':>18s} {'std':>8s}")
print("-" * 55)

for step in range(50):
    Ks = [all_traces[s][step]['K'] for s in range(n_seeds)]
    Sns = [all_traces[s][step]['schmidt'] for s in range(n_seeds)]
    BTs = [all_traces[s][step]['born_theta_eig'] for s in range(n_seeds)]

    mean_K = sum(Ks) / n_seeds
    mean_Sn = sum(Sns) / n_seeds
    mean_BT = sum(BTs) / n_seeds
    std_BT = (sum((b - mean_BT)**2 for b in BTs) / n_seeds)**0.5

    if step < 10 or step % 5 == 0:
        print(f"  {step:4d} {mean_K:>8.4f} {mean_Sn:>8.3f} {mean_BT:>18.5f} {std_BT:>8.5f}")


# ═══════════════════════════════════════════════════════════════
#  K trajectory: how it approaches K*
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- K trajectory toward K* = {K_STAR:.6f} ---\n")

# Average K over seeds
print(f"  {'step':>4s} {'mean_K':>10s} {'K - K*':>10s} {'|K-K*|/K*':>10s}")
print("-" * 40)

for step in range(50):
    Ks = [all_traces[s][step]['K'] for s in range(n_seeds)]
    mean_K = sum(Ks) / n_seeds
    diff = mean_K - K_STAR
    rel = abs(diff) / K_STAR

    if step < 10 or step % 5 == 0:
        print(f"  {step:4d} {mean_K:>10.6f} {diff:>+10.6f} {rel:>10.4f}")


# ═══════════════════════════════════════════════════════════════
#  The crystal's Born(θ)_eig at different step counts
# ═══════════════════════════════════════════════════════════════

print(f"\n\n--- Born(θ)_eig of the AVERAGED crystal at different step counts ---\n")

for n_steps in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 80, 100]:
    crystal = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(50):
        c, _ = entangle_with_trace(identity_corr, seed=seed, n_steps=n_steps)
        crystal += c
    crystal /= 50

    eigvals, eigvecs = torch.linalg.eig(crystal)
    dom_idx = eigvals.abs().argmax()
    v0 = eigvecs[:, dom_idx]
    bp = v0.abs().pow(2) / v0.abs().pow(2).sum()
    born_theta = bp[H].item()
    sn = schmidt_number(crystal)

    print(f"  {n_steps:3d} steps: Born(θ)_eig = {born_theta:.6f}, Schmidt = {sn:.3f}")

print(f"\n  10/33 = {10/33:.6f}")


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
