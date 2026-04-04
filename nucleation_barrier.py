"""
The nucleation barrier: what K-value produces the entropy peak?

Principle 57 showed entropy peaks at step ~15-20 then falls.
K converges to K* at step ~40. The entropy peak PRECEDES equilibrium.

Questions:
  1. At what K does entropy peak? Is it a universal K-threshold?
  2. Does the entropy-vs-K curve (not entropy-vs-step) collapse
     across evidence types? If so, K parameterizes the entire flow.
  3. What is the entropy at K* equilibrium? Is it a simple function of H?

If K is the only universal (Principle 57), then S(K) should be
universal too — the entropy is a function of the conflict, not of time.
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
    """Entangle and track the flow."""
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

        bp_2d = bp_norm.reshape(MASS_DIM, MASS_DIM)
        marg = bp_2d.sum(dim=1)
        born_theta = marg[H].item()

        # Sample evidence
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


# Evidence types
identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)
inverse_corr = torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)
cycle_corr = torch.tensor([[0,1,0],[0,0,1],[1,0,0]], dtype=torch.float32)
random_corr = torch.ones(H, H, dtype=torch.float32) / (H * H)

# Also test asymmetric evidence
exponential_corr = torch.tensor([[1,0,0],[.6,.4,0],[0,0,1]], dtype=torch.float32)

corrs = {
    "identity": identity_corr,
    "inverse": inverse_corr,
    "cycle": cycle_corr,
    "random": random_corr,
    "exponential": exponential_corr,
}


print("=" * 80)
print("THE NUCLEATION BARRIER: S(K) ACROSS EVIDENCE TYPES")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Collect S(K) curves for each evidence type
# ═══════════════════════════════════════════════════════════════

n_seeds = 30
n_steps = 70

all_curves = {}
for name, corr in corrs.items():
    sk_pairs = []  # (K, entropy) pairs
    for seed in range(n_seeds):
        _, trace = entangle_full_trace(corr, seed=seed, n_steps=n_steps)
        for t in trace:
            sk_pairs.append((t['K'], t['entropy'], t['step'], t['schmidt']))

    all_curves[name] = sk_pairs


# ═══════════════════════════════════════════════════════════════
#  1. At what K does entropy peak?
# ═══════════════════════════════════════════════════════════════

print(f"\n--- K-value at entropy peak ---\n")
for name in corrs:
    # Find the entropy peak for each seed
    peak_Ks = []
    for seed in range(n_seeds):
        _, trace = entangle_full_trace(corrs[name], seed=seed, n_steps=n_steps)
        peak_step = max(range(n_steps), key=lambda s: trace[s]['entropy'])
        peak_Ks.append(trace[peak_step]['K'])

    avg_K = sum(peak_Ks) / len(peak_Ks)
    std_K = (sum((k - avg_K)**2 for k in peak_Ks) / len(peak_Ks)) ** 0.5
    print(f"  {name:>12s}: K_peak = {avg_K:.5f} ± {std_K:.5f}")

print(f"\n  Reference: K* = {K_STAR:.5f}")
print(f"  K*/2 = {K_STAR/2:.5f}")
print(f"  K*/3 = {K_STAR/3:.5f}")


# ═══════════════════════════════════════════════════════════════
#  2. S(K) curve: bin by K, average S
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("S(K) CURVE: ENTROPY AS A FUNCTION OF CONFLICT")
print("="*80)

# K bins
K_bins = [i * 0.02 for i in range(15)]  # 0.00, 0.02, ..., 0.28

print(f"\n  {'K_bin':>8s}", end="")
for name in corrs:
    print(f"  {name:>10s}", end="")
print(f"  {'spread%':>10s}")
print("-" * 80)

for kb in K_bins:
    k_lo, k_hi = kb, kb + 0.02
    row_vals = []
    for name in corrs:
        in_bin = [s for k, s, step, sc in all_curves[name] if k_lo <= k < k_hi]
        if in_bin:
            avg_s = sum(in_bin) / len(in_bin)
        else:
            avg_s = float('nan')
        row_vals.append(avg_s)

    valid = [v for v in row_vals if not math.isnan(v)]
    mean_valid = sum(valid) / len(valid) if valid else 0
    if len(valid) >= 2 and abs(mean_valid) > 1e-10:
        spread = (max(valid) - min(valid)) / mean_valid * 100
    else:
        spread = 0.0

    print(f"  [{kb:.2f},{kb+0.02:.2f})", end="")
    for v in row_vals:
        if math.isnan(v):
            print(f"  {'---':>10s}", end="")
        else:
            print(f"  {v:>10.4f}", end="")
    print(f"  {spread:>9.1f}%")


# ═══════════════════════════════════════════════════════════════
#  3. Entropy at K* equilibrium
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("ENTROPY AT K* EQUILIBRIUM")
print("="*80)

for name in corrs:
    # Get entropy at steps 35-45 (near K* equilibrium)
    equil_entropies = []
    for seed in range(n_seeds):
        _, trace = entangle_full_trace(corrs[name], seed=seed, n_steps=n_steps)
        for t in trace[35:46]:
            equil_entropies.append(t['entropy'])

    avg_s = sum(equil_entropies) / len(equil_entropies)
    print(f"  {name:>12s}: S(K*) = {avg_s:.5f}")

# Candidate expressions
print(f"\n  Candidate expressions for S(K*) of bijection crystals:")
s_bijection = sum(
    sum(t['entropy'] for t in entangle_full_trace(
        identity_corr, seed=s, n_steps=n_steps)[1][35:46])
    for s in range(n_seeds)
) / (n_seeds * 11)

candidates = {
    "ln(H+1) = ln(4)": math.log(H + 1),
    "ln(H²) = ln(9)": math.log(H**2),
    "H × Δ": H * DELTA,
    "2Δ + ln(H)": 2 * DELTA + math.log(H),
    "ln(2π)": math.log(2 * math.pi),
    "(H+1)/H × ln(H)": (H+1)/H * math.log(H),
    "ln(H) + Δ": math.log(H) + DELTA,
    "π/√3": math.pi / math.sqrt(3),
    "2 × ln(H)": 2 * math.log(H),
    "ln(H³) = ln(27)": math.log(H**3),
    "(1-K*) × ln(16)": (1-K_STAR) * math.log(16),
    "S_max × (1-K*)": math.log(16) * (1-K_STAR),
    "S_max × K*": math.log(16) * K_STAR,
    "S_max - 2Δ": math.log(16) - 2*DELTA,
    "Δ × H(H+1)/2": DELTA * H*(H+1)/2,
}

print(f"\n  Measured S(K*) for bijections: {s_bijection:.5f}")
print(f"\n  {'expression':>25s}  {'value':>8s}  {'diff':>8s}  {'sigma':>6s}")
for expr, val in sorted(candidates.items(), key=lambda x: abs(x[1] - s_bijection)):
    diff = val - s_bijection
    print(f"  {expr:>25s}  {val:>8.5f}  {diff:>+8.5f}")


# ═══════════════════════════════════════════════════════════════
#  4. The K-entropy phase portrait
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE K-ENTROPY PHASE PORTRAIT")
print("="*80)

# For each evidence type, show the trajectory in (K, S) space
# This should reveal if the dynamics follows a universal curve

print(f"\n  Identity crystal: (K, S) trajectory averaged over {n_seeds} seeds")
avg_trace = []
for step in range(n_steps):
    all_K = []
    all_S = []
    for seed in range(n_seeds):
        _, trace = entangle_full_trace(identity_corr, seed=seed, n_steps=n_steps)
        all_K.append(trace[step]['K'])
        all_S.append(trace[step]['entropy'])
    avg_trace.append((sum(all_K)/len(all_K), sum(all_S)/len(all_S)))

print(f"\n  {'step':>4s}  {'K':>8s}  {'S':>8s}  {'dS/dK':>10s}")
for step in range(0, n_steps, 3):
    K, S = avg_trace[step]
    if step >= 3:
        K_prev, S_prev = avg_trace[step-3]
        dK = K - K_prev
        dS = S - S_prev
        dSdK = dS / dK if abs(dK) > 1e-8 else float('inf')
        print(f"  {step:4d}  {K:>8.5f}  {S:>8.4f}  {dSdK:>+10.3f}")
    else:
        print(f"  {step:4d}  {K:>8.5f}  {S:>8.4f}")


# ═══════════════════════════════════════════════════════════════
#  5. Entropy-K relationship: is S = f(K) a universal function?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("IS S = f(K) UNIVERSAL?")
print("="*80)

# Compare the S(K) relationship during the ASCENT (before peak)
# vs the DESCENT (after peak) — hysteresis would break universality

print(f"\n  Identity crystal: S(K) on ascent vs descent")
print(f"  {'K':>8s}  {'S_ascent':>10s}  {'S_descent':>10s}  {'ratio':>8s}")

# Find peak step
peak_steps = []
for seed in range(n_seeds):
    _, trace = entangle_full_trace(identity_corr, seed=seed, n_steps=n_steps)
    peak_step = max(range(n_steps), key=lambda s: trace[s]['entropy'])
    peak_steps.append(peak_step)
avg_peak = int(sum(peak_steps) / len(peak_steps))

# Collect ascent and descent data
ascent_data = {}   # K_bin → [S values]
descent_data = {}

for seed in range(n_seeds):
    _, trace = entangle_full_trace(identity_corr, seed=seed, n_steps=n_steps)
    peak = max(range(n_steps), key=lambda s: trace[s]['entropy'])
    for i, t in enumerate(trace):
        k_bin = round(t['K'] / 0.01) * 0.01  # bin to nearest 0.01
        if i < peak:
            ascent_data.setdefault(k_bin, []).append(t['entropy'])
        elif i > peak:
            descent_data.setdefault(k_bin, []).append(t['entropy'])

common_bins = sorted(set(ascent_data.keys()) & set(descent_data.keys()))
for kb in common_bins[:15]:
    s_asc = sum(ascent_data[kb]) / len(ascent_data[kb])
    s_des = sum(descent_data[kb]) / len(descent_data[kb])
    ratio = s_des / s_asc if s_asc > 0.01 else float('nan')
    print(f"  {kb:>8.3f}  {s_asc:>10.4f}  {s_des:>10.4f}  {ratio:>8.4f}")


print(f"\n\n{'='*80}")
print("WHAT THE NUCLEATION BARRIER SHOWS")
print("="*80)
