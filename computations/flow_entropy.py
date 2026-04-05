"""
The flow yields structure. What does it cost?

Start from total ignorance (S=0, all mass on θ,θ).
As conflict accumulates, mass spreads to h-entries.
Entropy increases (more states occupied).
Structure crystallizes (Schmidt rises).

Track simultaneously:
  - Entropy of the Born distribution (how spread the mass is)
  - Schmidt (how entangled the crystal is)
  - K (the conflict rate)
  - Born(θ) of the marginal (how much ignorance remains)

The three eternals should appear in their natural habitat:
  Born floor (1/27): where Born(θ) stabilizes
  K* (7/30): where conflict equilibrates
  Δ/4: the rate at which structure forms

The UNIVERSAL features (same for all evidence types) = the flow.
The TYPE-DEPENDENT features = the meaning.
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


def entangle_full_trace(correlation, seed=42, n_steps=80, discount=0.3):
    """Entangle and track everything."""
    _random.seed(seed)
    corr = correlation.float()
    flat = corr.reshape(-1)
    flat = flat / flat.sum().clamp(min=EPS_LOG)

    joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j

    trace = []

    for step in range(n_steps):
        # Measure BEFORE this step's evidence
        bp = joint.abs().pow(2).reshape(-1)
        bp_norm = bp / bp.sum().clamp(min=EPS_LOG)

        # Entropy
        entropy = -sum(p.item() * math.log(max(p.item(), 1e-15))
                       for p in bp_norm if p.item() > 1e-15)

        # Schmidt
        sn = schmidt_number(joint) if step % 2 == 0 else trace[-1]['schmidt'] if trace else 1.0

        # Marginal Born(θ)
        bp_2d = bp_norm.reshape(MASS_DIM, MASS_DIM)
        marg = bp_2d.sum(dim=1)
        born_theta = marg[H].item()

        # Sample and combine evidence
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
            'step': step,
            'entropy': entropy,
            'schmidt': sn,
            'K': K.abs().item() if K.is_complex() else K.item(),
            'born_theta': born_theta,
        })

    return joint, trace


identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
inverse_corr = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)
cycle_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)
random_corr = torch.ones(H, H, dtype=torch.float32) / (H * H)

corrs = {
    "identity": identity_corr,
    "inverse": inverse_corr,
    "cycle": cycle_corr,
    "random": random_corr,
}


print("=" * 80)
print("THE FLOW YIELDS STRUCTURE")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Entropy trajectories for different evidence types
# ═══════════════════════════════════════════════════════════════

n_seeds = 20
n_steps = 70

avg_traces = {}
for name, corr in corrs.items():
    all_t = []
    for seed in range(n_seeds):
        _, trace = entangle_full_trace(corr, seed=seed, n_steps=n_steps)
        all_t.append(trace)

    # Average
    avg = []
    for step in range(n_steps):
        avg.append({
            'entropy': sum(all_t[s][step]['entropy'] for s in range(n_seeds)) / n_seeds,
            'schmidt': sum(all_t[s][step]['schmidt'] for s in range(n_seeds)) / n_seeds,
            'K': sum(all_t[s][step]['K'] for s in range(n_seeds)) / n_seeds,
            'born_theta': sum(all_t[s][step]['born_theta'] for s in range(n_seeds)) / n_seeds,
        })
    avg_traces[name] = avg


print(f"\n--- Entropy S(step) for different evidence types ---\n")
print(f"  {'step':>4s}", end="")
for name in corrs:
    print(f"  {name:>10s}", end="")
print(f"  {'max=ln16':>10s}")
print("-" * 60)

for step in range(0, n_steps, 5):
    print(f"  {step:4d}", end="")
    for name in corrs:
        print(f"  {avg_traces[name][step]['entropy']:>10.4f}", end="")
    print(f"  {math.log(16):>10.4f}")


# ═══════════════════════════════════════════════════════════════
#  What is UNIVERSAL vs TYPE-DEPENDENT?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("UNIVERSAL vs TYPE-DEPENDENT")
print("="*80)

print(f"\n  At step 40 (near K* equilibrium):")
print(f"  {'quantity':>12s}", end="")
for name in corrs:
    print(f"  {name:>10s}", end="")
print(f"  {'spread':>10s}")

for qty in ['entropy', 'schmidt', 'K', 'born_theta']:
    vals = [avg_traces[name][40][qty] for name in corrs]
    spread = max(vals) - min(vals)
    rel_spread = spread / (sum(vals)/len(vals)) * 100 if sum(vals) > 0 else 0

    print(f"  {qty:>12s}", end="")
    for v in vals:
        print(f"  {v:>10.4f}", end="")
    print(f"  {rel_spread:>9.1f}%")


# ═══════════════════════════════════════════════════════════════
#  The entropy cost of crystallization
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE ENTROPY COST OF CRYSTALLIZATION")
print("="*80)

# How much entropy does the flow gain when it yields a crystal?
# Start: S = 0 (total ignorance concentrated on θ,θ)
# After crystallization: S = S_crystal

for name in corrs:
    s_initial = avg_traces[name][0]['entropy']
    s_equil = avg_traces[name][40]['entropy']
    s_late = avg_traces[name][-1]['entropy']

    print(f"\n  {name}:")
    print(f"    S(0)  = {s_initial:.4f} (initial)")
    print(f"    S(40) = {s_equil:.4f} (equilibrium)")
    print(f"    S(69) = {s_late:.4f} (late)")
    print(f"    ΔS(crystallization) = {s_equil - s_initial:.4f}")

# Maximum possible entropy
s_max = math.log(MASS_DIM ** 2)  # ln(16)
print(f"\n  S_max = ln(16) = {s_max:.4f}")
print(f"  S_uniform = ln(16) = {s_max:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Where do the three eternals appear?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE THREE ETERNALS IN THE FLOW")
print("="*80)

# Use identity crystal as reference
t = avg_traces["identity"]

# 1. Born floor: where does Born(θ) of the MARGINAL reach 1/27?
print(f"\n  Born(θ) trajectory toward 1/27 = {BORN_FLOOR:.5f}:")
for step in range(0, n_steps, 5):
    bt = t[step]['born_theta']
    diff = bt - BORN_FLOOR
    print(f"    step {step:2d}: Born(θ) = {bt:.5f} ({'above' if diff > 0 else 'BELOW'} floor by {abs(diff):.5f})")

# 2. K*: already shown to converge by step 40-45
print(f"\n  K → K* = {K_STAR:.5f}:")
for step in [10, 20, 30, 40, 50, 60]:
    if step < n_steps:
        print(f"    step {step:2d}: K = {t[step]['K']:.5f}")

# 3. Δ/4 growth rate: entropy growth rate
print(f"\n  Entropy growth rate (should show Δ/4 = {DELTA/4:.5f}):")
for step in range(5, min(40, n_steps), 5):
    ds = (t[step]['entropy'] - t[step-5]['entropy']) / 5
    print(f"    steps {step-5}-{step}: dS/dt = {ds:.5f}")


# ═══════════════════════════════════════════════════════════════
#  The flow at equilibrium: what does the crystal look like?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE EQUILIBRIUM CRYSTAL FROM EACH EVIDENCE TYPE")
print("="*80)

for name, corr in corrs.items():
    crystal = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    for seed in range(30):
        c, _ = entangle_full_trace(corr, seed=seed, n_steps=50)
        crystal += c
    crystal /= 30

    bp = born_probabilities(crystal).reshape(MASS_DIM, MASS_DIM)
    sn = schmidt_number(crystal)

    # Diagonal dominance
    diag_sum = sum(bp[i,i].item() for i in range(MASS_DIM))
    off_diag = 1 - diag_sum

    print(f"\n  {name}:")
    print(f"    Schmidt = {sn:.3f}")
    print(f"    Diagonal sum = {diag_sum:.4f}, Off-diagonal = {off_diag:.4f}")
    print(f"    Born(θ,θ) = {bp[H,H].item():.5f}")
    print(f"    h-diagonal: [{', '.join(f'{bp[i,i].item():.4f}' for i in range(H))}]")


print(f"\n\n{'='*80}")
print("WHAT THE FLOW SHOWS")
print("="*80)
