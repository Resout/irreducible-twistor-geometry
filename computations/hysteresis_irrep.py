"""
Hysteresis in the S₃ irrep decomposition during entangling.

The entropy ascent (S rising, K < K_peak) and descent (S falling, K > K_peak)
pass through the same K values but at different entropies. What carries the
difference? Track the S₃ irrep fractions along the flow.

If the trivial/standard/sign fractions differ between ascent and descent
at the same K, the hysteresis lives in the gauge structure.
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


def s3_projectors():
    """Build S₃ projectors on [4,4] joint mass space.

    S₃ acts by permuting the first 3 indices (h₀,h₁,h₂), fixing θ (index 3).
    """
    # The 6 elements of S₃ as 4×4 permutation matrices
    perms_cycles = [
        [0,1,2,3],  # identity
        [1,0,2,3],  # (01)
        [2,1,0,3],  # (02)
        [0,2,1,3],  # (12)
        [1,2,0,3],  # (012)
        [2,0,1,3],  # (021)
    ]
    signs = [1, -1, -1, -1, 1, 1]  # sign of each permutation

    def perm_matrix(p):
        M = torch.zeros(MASS_DIM, MASS_DIM)
        for i, j in enumerate(p):
            M[j, i] = 1.0
        return M

    perm_mats = [perm_matrix(p) for p in perms_cycles]

    # Action on [4,4] joint mass: (σ · M)_{ij} = M_{σ⁻¹(i), σ⁻¹(j)}
    def act_on_joint(P, joint):
        return P @ joint @ P.T

    # Trivial projector: average over all elements
    def project_trivial(joint):
        result = torch.zeros_like(joint)
        for P in perm_mats:
            result += act_on_joint(P, joint.real.float()).to(joint.dtype)
        return result / 6

    # Sign projector: weighted average with sign character
    def project_sign(joint):
        result = torch.zeros_like(joint)
        for P, s in zip(perm_mats, signs):
            result += s * act_on_joint(P, joint.real.float()).to(joint.dtype)
        return result / 6

    # Standard projector: I - trivial - sign
    def project_standard(joint):
        return joint - project_trivial(joint) - project_sign(joint)

    return project_trivial, project_sign, project_standard


def irrep_fractions(joint, proj_triv, proj_sign, proj_std):
    """Fraction of Born mass in each irrep sector."""
    bp = joint.abs().pow(2)
    total = bp.sum().item()
    if total < 1e-15:
        return 1/3, 1/3, 1/3

    triv = proj_triv(joint).abs().pow(2).sum().item() / total
    sign = proj_sign(joint).abs().pow(2).sum().item() / total
    std = proj_std(joint).abs().pow(2).sum().item() / total

    return triv, sign, std


def entangle_with_irreps(correlation, seed=42, n_steps=70, discount=0.3):
    """Entangle and track irrep decomposition."""
    _random.seed(seed)
    corr = correlation.float()
    flat = corr.reshape(-1)
    flat = flat / flat.sum().clamp(min=EPS_LOG)

    joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j

    proj_triv, proj_sign, proj_std = s3_projectors()
    trace = []

    for step in range(n_steps):
        bp = joint.abs().pow(2).reshape(-1)
        bp_norm = bp / bp.sum().clamp(min=EPS_LOG)
        entropy = -sum(p.item() * math.log(max(p.item(), 1e-15))
                       for p in bp_norm if p.item() > 1e-15)
        sn = schmidt_number(joint)

        f_triv, f_sign, f_std = irrep_fractions(joint, proj_triv, proj_sign, proj_std)

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


# Run for identity correlation
identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

print("=" * 80)
print("HYSTERESIS IN THE IRREP DECOMPOSITION")
print("=" * 80)

n_seeds = 30
n_steps = 70

# Collect all traces
all_traces = []
for seed in range(n_seeds):
    _, trace = entangle_with_irreps(identity_corr, seed=seed, n_steps=n_steps)
    all_traces.append(trace)

# Average trajectory
avg = []
for step in range(n_steps):
    avg.append({k: sum(all_traces[s][step][k] for s in range(n_seeds)) / n_seeds
                for k in all_traces[0][0].keys()})

print(f"\n--- Averaged trajectory (identity, {n_seeds} seeds) ---\n")
print(f"  {'step':>4s}  {'K':>8s}  {'S':>8s}  {'Schmidt':>8s}  {'trivial':>8s}  {'sign':>8s}  {'standard':>8s}")
print("-" * 65)
for step in range(0, n_steps, 3):
    a = avg[step]
    print(f"  {step:4d}  {a['K']:>8.5f}  {a['entropy']:>8.4f}  {a['schmidt']:>8.3f}  "
          f"{a['trivial']:>8.4f}  {a['sign']:>8.4f}  {a['standard']:>8.4f}")


# Ascent vs descent comparison at the same K
print(f"\n\n{'='*80}")
print("IRREP FRACTIONS: ASCENT vs DESCENT")
print("="*80)

# Find peak step for each seed
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
        entry = (t['entropy'], t['trivial'], t['sign'], t['standard'], t['schmidt'])
        if i < peak:
            ascent_data.setdefault(k_bin, []).append(entry)
        elif i > peak:
            descent_data.setdefault(k_bin, []).append(entry)

common = sorted(set(ascent_data.keys()) & set(descent_data.keys()))

print(f"\n  {'K':>6s}  {'S_asc':>8s}  {'S_des':>8s}  {'triv_a':>8s}  {'triv_d':>8s}  {'std_a':>8s}  {'std_d':>8s}  {'Sc_a':>8s}  {'Sc_d':>8s}")
print("-" * 85)

for kb in common[:20]:
    a_vals = ascent_data[kb]
    d_vals = descent_data[kb]
    a_avg = [sum(x[i] for x in a_vals)/len(a_vals) for i in range(5)]
    d_avg = [sum(x[i] for x in d_vals)/len(d_vals) for i in range(5)]

    print(f"  {kb:>6.3f}  {a_avg[0]:>8.4f}  {d_avg[0]:>8.4f}  "
          f"{a_avg[1]:>8.4f}  {d_avg[1]:>8.4f}  "
          f"{a_avg[3]:>8.4f}  {d_avg[3]:>8.4f}  "
          f"{a_avg[4]:>8.4f}  {d_avg[4]:>8.4f}")


print(f"\n\n{'='*80}")
print("WHAT THE HYSTERESIS SHOWS")
print("="*80)
