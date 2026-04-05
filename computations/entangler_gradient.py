"""
Is the entangler a gradient flow?

Born(θ)_eig descends monotonically from 1.0 to 0.0 during entangling.
If this is a gradient flow on a potential V(Born(θ)), then:

  d(Born(θ))/dt = -∂V/∂(Born(θ))

The descent rate should depend on Born(θ) itself — a self-consistent
dynamics. At the K* equilibrium, the rate should vanish (or change
character).

Questions:
1. What is the functional form of Born(θ)(step)?
2. Is it exponential? Power law? Logistic?
3. Does the descent rate depend only on Born(θ), or also on step number?
4. Is there a natural time variable that makes the descent a straight line?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import random as _random
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            schmidt_number, EPS_LOG, EPS_DIV, EPS_NORM)

torch.set_grad_enabled(False)


def entangle_born_trace(correlation, seed=42, n_steps=100, discount=0.3):
    """Run entangler, return Born(θ)_eig at each step."""
    _random.seed(seed)

    corr = correlation.float()
    flat = corr.reshape(-1)
    flat = flat / flat.sum().clamp(min=EPS_LOG)

    joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j

    born_trace = []

    for step in range(n_steps):
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

        try:
            eigvals, eigvecs = torch.linalg.eig(joint)
            dom_idx = eigvals.abs().argmax()
            v0 = eigvecs[:, dom_idx]
            bp = v0.abs().pow(2) / v0.abs().pow(2).sum()
            born_trace.append(bp[H].item())
        except:
            born_trace.append(born_trace[-1] if born_trace else 1.0)

    return born_trace


identity_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)


print("=" * 80)
print("IS THE ENTANGLER A GRADIENT FLOW?")
print("=" * 80)


# ═══════════════════════════════════════════════════════════════
#  Average trajectory over many seeds
# ═══════════════════════════════════════════════════════════════

n_seeds = 40
n_steps = 80
all_traces = []
for seed in range(n_seeds):
    trace = entangle_born_trace(identity_corr, seed=seed, n_steps=n_steps)
    all_traces.append(trace)

# Average at each step
avg_born = []
for step in range(n_steps):
    vals = [all_traces[s][step] for s in range(n_seeds)]
    avg_born.append(sum(vals) / len(vals))


# ═══════════════════════════════════════════════════════════════
#  Test 1: Is it exponential? log(Born(θ)) vs step
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Test 1: Exponential fit — log(Born(θ)) vs step ---\n")

print(f"  {'step':>4s} {'Born(θ)':>10s} {'ln(Born)':>10s} {'Δ_ln':>10s}")
print("-" * 40)
for i in range(0, n_steps, 5):
    b = avg_born[i]
    lb = math.log(max(b, 1e-10))
    d_ln = (math.log(max(avg_born[i], 1e-10)) - math.log(max(avg_born[max(0,i-5)], 1e-10))) / 5 if i > 0 else 0
    print(f"  {i:4d} {b:>10.5f} {lb:>10.4f} {d_ln:>10.5f}")


# ═══════════════════════════════════════════════════════════════
#  Test 2: Is it a power law? log(Born(θ)) vs log(step)
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Test 2: Power law — log(Born(θ)) vs log(step) ---\n")

print(f"  {'step':>4s} {'Born(θ)':>10s} {'ln(step+1)':>10s} {'ln(Born)':>10s} {'slope':>8s}")
print("-" * 50)
for i in [1, 2, 5, 10, 15, 20, 30, 40, 50, 60, 70]:
    if i < n_steps:
        b = avg_born[i]
        ls = math.log(i + 1)
        lb = math.log(max(b, 1e-10))
        if i > 1:
            prev_i = max(1, i - 5)
            slope = (math.log(max(avg_born[i], 1e-10)) - math.log(max(avg_born[prev_i], 1e-10))) / \
                    (math.log(i+1) - math.log(prev_i+1))
        else:
            slope = 0
        print(f"  {i:4d} {b:>10.5f} {ls:>10.4f} {lb:>10.4f} {slope:>8.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 3: Is the descent rate a function of Born(θ)?
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Test 3: Descent rate vs current Born(θ) ---\n")

# Compute -Δ(Born)/Δstep at each step
print(f"  {'step':>4s} {'Born(θ)':>10s} {'-dB/dt':>10s} {'-dB/B':>10s} {'K':>8s}")
print("-" * 50)

# Also get K from a reference trace
_, full_trace = __import__('computations.entangler_dynamics', fromlist=['entangle_with_trace']).entangle_with_trace(identity_corr, seed=42, n_steps=n_steps)

# Hmm, can't import like that. Let me just compute the rates.
for i in range(1, n_steps, 3):
    b = avg_born[i]
    b_prev = avg_born[i-1]
    rate = -(b - b_prev)  # positive when born decreases
    rel_rate = rate / max(b, 1e-10)  # relative rate
    print(f"  {i:4d} {b:>10.5f} {rate:>10.5f} {rel_rate:>10.5f}")


# ═══════════════════════════════════════════════════════════════
#  Test 4: Is -d(ln Born)/dt = constant? (exponential decay)
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Test 4: -d(ln Born)/dt vs step ---\n")

rates = []
for i in range(1, min(60, n_steps)):
    b = avg_born[i]
    b_prev = avg_born[i-1]
    if b > 0.01 and b_prev > 0.01:
        r = -(math.log(b) - math.log(b_prev))
        rates.append((i, b, r))

print(f"  {'step':>4s} {'Born(θ)':>10s} {'-d(lnB)/dt':>12s}")
print("-" * 30)
for i, b, r in rates[::3]:
    print(f"  {i:4d} {b:>10.5f} {r:>12.5f}")

if rates:
    early_rates = [r for _, _, r in rates[:15]]
    late_rates = [r for _, _, r in rates[30:50]]
    print(f"\n  Early mean (steps 1-15): {sum(early_rates)/len(early_rates):.5f}")
    if late_rates:
        print(f"  Late mean (steps 30-50):  {sum(late_rates)/len(late_rates):.5f}")
    print(f"  Δ/4 = {DELTA/4:.5f}")
    print(f"  Δ = {DELTA:.5f}")


# ═══════════════════════════════════════════════════════════════
#  Test 5: What coordinate makes the descent linear?
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Test 5: Transformed coordinates ---\n")

# Try logit(Born(θ)) = ln(B/(1-B))
print(f"  Logit transform: logit(Born) vs step")
print(f"  {'step':>4s} {'Born(θ)':>10s} {'logit':>10s} {'Δ_logit':>10s}")
print("-" * 40)
for i in range(0, min(60, n_steps), 5):
    b = avg_born[i]
    if 0.01 < b < 0.99:
        logit = math.log(b / (1 - b))
        if i > 0 and 0.01 < avg_born[i-5] < 0.99:
            prev_logit = math.log(avg_born[i-5] / (1 - avg_born[i-5]))
            d_logit = (logit - prev_logit) / 5
        else:
            d_logit = 0
        print(f"  {i:4d} {b:>10.5f} {logit:>10.4f} {d_logit:>10.5f}")


# ═══════════════════════════════════════════════════════════════
#  Test 6: Born(θ) ≈ 1/(1 + α·step^β)?
# ═══════════════════════════════════════════════════════════════

print(f"\n--- Test 6: Fit Born(θ) = 1/(1 + α·step^β) ---\n")

# If Born = 1/(1+α·t^β), then 1/Born - 1 = α·t^β
# ln(1/Born - 1) = ln(α) + β·ln(t)
print(f"  {'step':>4s} {'1/Born-1':>10s} {'ln(1/B-1)':>10s} {'ln(step)':>10s}")
print("-" * 45)
for i in [1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 50, 60]:
    if i < n_steps:
        b = avg_born[i]
        if b > 0.01:
            inv_b_minus_1 = 1/b - 1
            print(f"  {i:4d} {inv_b_minus_1:>10.4f} {math.log(max(inv_b_minus_1, 1e-10)):>10.4f} {math.log(i):>10.4f}")

# Linear regression on ln(1/B-1) vs ln(t) for steps 5-50
xs, ys = [], []
for i in range(5, min(50, n_steps)):
    b = avg_born[i]
    if b > 0.01:
        xs.append(math.log(i))
        ys.append(math.log(1/b - 1))

if len(xs) > 2:
    n = len(xs)
    sx = sum(xs); sy = sum(ys); sxx = sum(x**2 for x in xs); sxy = sum(x*y for x,y in zip(xs,ys))
    beta = (n*sxy - sx*sy) / (n*sxx - sx**2)
    ln_alpha = (sy - beta*sx) / n
    alpha = math.exp(ln_alpha)

    print(f"\n  Fit: Born(θ) ≈ 1/(1 + {alpha:.4f} × step^{beta:.4f})")
    print(f"  β = {beta:.4f}")
    print(f"  α = {alpha:.4f}")

    # Check quality of fit
    print(f"\n  {'step':>4s} {'measured':>10s} {'fitted':>10s} {'error':>8s}")
    for i in [5, 10, 20, 30, 40, 50]:
        if i < n_steps:
            b_meas = avg_born[i]
            b_fit = 1 / (1 + alpha * i**beta)
            err = abs(b_meas - b_fit) / b_meas * 100
            print(f"  {i:4d} {b_meas:>10.5f} {b_fit:>10.5f} {err:>7.1f}%")


print(f"\n\n{'='*80}")
print("FINDINGS")
print("="*80)
