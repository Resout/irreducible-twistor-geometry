"""
Does Δ appear in the entangler's dynamics?

The structural filter rate Δ = -ln(1 - K*) ≈ 0.266 governs
composition decay. The entangler operates at K*. Is Δ visible
in the entangler's step-by-step evolution?

Possible manifestations:
1. Schmidt growth rate per step → Δ?
2. θ-fingerprint convergence rate → Δ?
3. Information gain per step (entropy decrease) → Δ?
4. The ratio between early and late K → Δ?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import random as _random
from solver.algebra import (
    H, MASS_DIM, K_STAR, DELTA, BORN_FLOOR,
    schmidt_number, born_probabilities, born_entropy,
    EPS_LOG, EPS_DIV, EPS_NORM,
)
from solver.crystals import Entangler, classify_relationship

torch.set_grad_enabled(False)

id_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)


def theta_fingerprint(joint):
    fp = torch.stack([
        joint[0, 3], joint[1, 3], joint[2, 3],
        joint[3, 0], joint[3, 1], joint[3, 2],
    ])
    norm = fp.abs().pow(2).sum().sqrt()
    if norm > 1e-10:
        fp = fp / norm
    return fp


def run_entangler_detailed(corr, n_steps=50, discount=0.3, seed=42):
    """Run entangler recording detailed measurements at each step."""
    _random.seed(seed)

    joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
    joint[-1, -1] = 1.0 + 0j

    flat = corr.reshape(-1)
    flat = flat / flat.sum().clamp(min=EPS_LOG)

    records = []

    for step in range(n_steps):
        idx = _random.choices(range(H * H), weights=flat.tolist())[0]
        hi, hj = idx // H, idx % H

        ev_a = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev_b = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev_a[hi] = 0.55 + 0.06j; ev_b[hj] = 0.55 + 0.06j
        for k in range(H):
            noise = _random.gauss(0, 0.02)
            if k != hi: ev_a[k] = 0.05 + noise * 1j
            if k != hj: ev_b[k] = 0.05 + noise * 1j
        ev_a[H] = (1.0 - ev_a[:H].real.sum()) - ev_a[:H].imag.sum() * 1j
        ev_b[H] = (1.0 - ev_b[:H].real.sum()) - ev_b[:H].imag.sum() * 1j

        joint_dim = MASS_DIM ** 2; jH = joint_dim - 1
        joint_ev = torch.einsum("m,n->mn", ev_a, ev_b).reshape(joint_dim)
        ev_s = joint_ev[:jH] * discount
        ev_t = 1.0 + 0j - ev_s.sum()
        disc_ev = torch.cat([ev_s, ev_t.unsqueeze(0)])

        curr = joint.reshape(joint_dim)
        cs1, ct1 = curr[:jH], curr[jH:]
        cs2, ct2 = disc_ev[:jH], disc_ev[jH:]
        K = cs1.sum() * cs2.sum() - (cs1 * cs2).sum()
        K_val = abs(K.item()) if isinstance(K, torch.Tensor) else abs(K)

        comb_s = cs1 * cs2 + cs1 * ct2 + ct1 * cs2
        comb_t = ct1 * ct2
        norm = 1.0 - K
        safe_norm = norm if abs(norm) > EPS_NORM else torch.tensor(EPS_NORM + 0j)
        updated = torch.cat([comb_s, comb_t]) / safe_norm
        re_sum = updated.real.sum()
        if abs(re_sum) > EPS_DIV:
            updated = updated / re_sum
        joint = updated.reshape(MASS_DIM, MASS_DIM)

        sn = schmidt_number(joint)
        # θ-mass fraction
        theta_mass = joint[3, :].abs().pow(2).sum().item() / joint.abs().pow(2).sum().item()
        # Singular values
        _, S, _ = torch.linalg.svd(joint)
        sv = S.abs().tolist()

        records.append({
            "step": step, "schmidt": sn, "K": K_val,
            "theta_mass": theta_mass, "sv": sv,
        })

    return records


# ═══════════════════════════════════════════════════════════════
#  Test 1: Schmidt growth rate
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("TEST 1: SCHMIDT GROWTH RATE")
print(f"If Δ = {DELTA:.4f} governs entangling, Schmidt should")
print("grow at a rate related to Δ.")
print("=" * 70)

# Average over many seeds
n_seeds = 50
avg_schmidts = [0.0] * 50

for seed in range(n_seeds):
    records = run_entangler_detailed(id_corr, n_steps=50, seed=seed)
    for r in records:
        avg_schmidts[r["step"]] += r["schmidt"] / n_seeds

print(f"\n  {'step':>5s} {'avg_Schmidt':>12s} {'ΔSchmidt':>10s} {'growth_rate':>12s}")
print("  " + "-" * 42)
for i in range(0, 50, 5):
    ds = avg_schmidts[i] - avg_schmidts[i-1] if i > 0 else 0
    # Growth rate relative to remaining capacity (3.27 - current)
    remaining = 3.272 - avg_schmidts[i]
    growth_rate = ds / remaining if remaining > 0.01 else 0
    print(f"  {i:>5d} {avg_schmidts[i]:>12.4f} {ds:>10.4f} {growth_rate:>12.4f}")

# Fit: Schmidt(n) = S_max * (1 - exp(-rate * n))
# If rate = Δ, then Schmidt(n) should approach 3.27 exponentially
# Let's measure the effective rate

# Schmidt excess above ignorance (1.0): S_excess = Schmidt - 1
# If S_excess ≈ S_max_excess * (1 - exp(-rate * n)):
# ln(1 - S_excess/S_max_excess) = -rate * n

S_max = 3.27
rates = []
for i in range(5, 45, 5):
    s_excess = avg_schmidts[i] - 1.0
    s_max_excess = S_max - 1.0
    ratio = s_excess / s_max_excess
    if 0 < ratio < 1:
        rate = -math.log(1 - ratio) / i
        rates.append(rate)
        if i % 10 == 0:
            print(f"  At step {i}: effective rate = {rate:.4f}")

if rates:
    avg_rate = sum(rates) / len(rates)
    print(f"\n  Average growth rate: {avg_rate:.4f}")
    print(f"  Δ = {DELTA:.4f}")
    print(f"  Ratio rate/Δ = {avg_rate/DELTA:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 2: θ-mass decay rate
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 2: θ-MASS FRACTION DECAY")
print("The θ-row mass fraction should decay as the crystal")
print("builds up structure. Is the decay rate Δ?")
print("=" * 70)

avg_theta = [0.0] * 50
for seed in range(n_seeds):
    records = run_entangler_detailed(id_corr, n_steps=50, seed=seed)
    for r in records:
        avg_theta[r["step"]] += r["theta_mass"] / n_seeds

print(f"\n  {'step':>5s} {'θ-mass':>10s} {'ln(θ-mass)':>12s}")
print("  " + "-" * 30)
for i in range(0, 50, 5):
    ln_theta = math.log(avg_theta[i]) if avg_theta[i] > 1e-10 else -99
    print(f"  {i:>5d} {avg_theta[i]:>10.4f} {ln_theta:>12.4f}")

# Fit exponential decay: θ(n) = θ_0 * exp(-rate * n)
# ln(θ) = ln(θ_0) - rate * n
# Rate from consecutive points
theta_rates = []
for i in range(5, 45, 5):
    if avg_theta[i] > 0 and avg_theta[i-5] > 0:
        rate = -(math.log(avg_theta[i]) - math.log(avg_theta[i-5])) / 5
        theta_rates.append(rate)

if theta_rates:
    avg_theta_rate = sum(theta_rates) / len(theta_rates)
    print(f"\n  θ-decay rate: {avg_theta_rate:.4f}")
    print(f"  Δ = {DELTA:.4f}")
    print(f"  Ratio θ-rate/Δ = {avg_theta_rate/DELTA:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 3: Singular value dynamics
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 3: SINGULAR VALUE EVOLUTION")
print("How do the 4 singular values evolve during entangling?")
print("=" * 70)

# Single trajectory
records = run_entangler_detailed(id_corr, n_steps=50, seed=42)

print(f"\n  {'step':>5s} {'σ₀':>8s} {'σ₁':>8s} {'σ₂':>8s} {'σ₃':>8s} {'σ₀/σ₃':>8s} {'Schmidt':>8s}")
print("  " + "-" * 55)
for r in records:
    if r["step"] % 5 == 0 or r["step"] < 3:
        sv = r["sv"]
        ratio = sv[0] / sv[3] if sv[3] > 1e-10 else float('inf')
        print(f"  {r['step']:>5d} {sv[0]:>8.4f} {sv[1]:>8.4f} "
              f"{sv[2]:>8.4f} {sv[3]:>8.4f} {ratio:>8.2f} {r['schmidt']:>8.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 4: Composition decay rate measurement
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 4: COMPOSITION DECAY — direct measurement of Δ")
print("Build a crystal, compose it with itself N times,")
print("measure Schmidt decay. Is the rate exactly Δ?")
print("=" * 70)

from solver.crystals import compose

# Build identity crystal
ent = Entangler(id_corr).build()
joint_0 = ent.joint.clone()
sn_0 = schmidt_number(joint_0)

print(f"\n  Starting crystal: Schmidt = {sn_0:.4f}")
print(f"\n  {'step':>5s} {'Schmidt':>9s} {'excess':>9s} {'ln(excess)':>12s} {'rate':>8s}")
print("  " + "-" * 48)

current = joint_0.clone()
prev_ln = None
for step in range(8):
    sn = schmidt_number(current)
    excess = sn - 1.0
    ln_excess = math.log(excess) if excess > 1e-6 else -99
    rate = -(ln_excess - prev_ln) if prev_ln is not None else 0
    prev_ln = ln_excess
    print(f"  {step:>5d} {sn:>9.4f} {excess:>9.4f} {ln_excess:>12.4f} {rate:>8.4f}")
    current = compose(current, joint_0)

print(f"\n  Predicted rate Δ = {DELTA:.4f}")
print(f"  Predicted rate 2Δ = {2*DELTA:.4f} (if Schmidt ∝ |m|², doubling)")


# ═══════════════════════════════════════════════════════════════
#  Test 5: The Δ-K*-Schmidt triangle
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 5: THE Δ–K*–SCHMIDT TRIANGLE")
print("Three quantities, all from H=3. How are they connected?")
print("=" * 70)

print(f"\n  H = {H}")
print(f"  K* = (H²-H+1)/(H(H²+1)) = {K_STAR:.6f} = 7/30")
print(f"  Δ  = -ln(1-K*) = {DELTA:.6f}")
print(f"  Schmidt_max ≈ 3.27 (empirical, from entangler)")
print(f"  BORN_FLOOR = 1/H³ = {BORN_FLOOR:.6f}")

# Theoretical Schmidt maximum from K*
# At K*, the next combination brings factor (1-K*) to the
# off-diagonal singletons. The diagonal is preserved.
# Schmidt = 1/Σ(p_i²) where p_i are normalized singular values.

# For a [4,4] matrix with maximally spread singular values:
# 4 equal SVs → Schmidt = 4 (maximum possible)
# The Born floor limits this.

# What is the theoretical Schmidt at K* equilibrium?
# S = 1/(p₀² + p₁² + p₂² + p₃²)
# where Σpᵢ = 1 and p₃ ≥ BORN_FLOOR

# For H=3: maximum Schmidt occurs when p₀=p₁=p₂ are as equal
# as possible while respecting p₃ ≥ 1/27
# If p₃ = 1/27: remaining = 26/27, split 3 ways: p_i = 26/81 each
# S = 1/(3*(26/81)² + (1/27)²) = 1/(3*0.1030 + 0.00137) = 1/0.3103 = 3.223

S_theoretical = 1.0 / (3 * (26/81)**2 + (1/27)**2)
print(f"\n  Theoretical Schmidt at Born floor = 1/H³:")
print(f"    S = 1/(3×(26/81)² + (1/27)²) = {S_theoretical:.4f}")
print(f"    Empirical: ≈ 3.27")
print(f"    Gap: {3.27 - S_theoretical:.4f}")

# What Born(θ) gives Schmidt = 3.27?
# S = 1/(3*p_h² + p_θ²) where 3*p_h + p_θ = 1
# With equal p_h: p_h = (1-p_θ)/3
# S = 1/(3*((1-p_θ)/3)² + p_θ²) = 1/((1-p_θ)²/3 + p_θ²)
# 1/3.27 = (1-p_θ)²/3 + p_θ² = (1-2p_θ+p_θ²)/3 + p_θ²
# = 1/3 - 2p_θ/3 + p_θ²/3 + p_θ² = 1/3 - 2p_θ/3 + 4p_θ²/3

target = 1.0 / 3.27
# 4p²/3 - 2p/3 + 1/3 = target
# 4p² - 2p + 1 = 3*target
# 4p² - 2p + (1 - 3*target) = 0
a, b, c_coeff = 4, -2, 1 - 3*target
disc = b**2 - 4*a*c_coeff
if disc >= 0:
    p_theta = (-b - math.sqrt(disc)) / (2*a)
    print(f"\n  Born(θ) that gives Schmidt 3.27: {p_theta:.6f}")
    print(f"  Born floor 1/H³ = {BORN_FLOOR:.6f}")
    print(f"  Ratio: {p_theta / BORN_FLOOR:.3f}")


print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
