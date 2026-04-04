"""
The mass gap as relaxation rate.

Born(θ) = 1/27 is the equilibrium. Perturb it. Measure the
relaxation rate. That rate is the spectral gap.

The mass gap is the cost of changing HOW MUCH you know.
You can change WHAT you know freely (singletons track evidence).
But changing the LEVEL of ignorance costs energy proportional to Δ.

Also: is the compositional trace limit exactly 1/(H+1)?
And: what does the θ-fingerprint look like at the fixed point?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import random as _random
from solver.algebra import (
    H, MASS_DIM, K_STAR, DELTA, BORN_FLOOR,
    schmidt_number, born_probabilities, born_entropy,
    ds_combine, enforce_born_floor, discount_mass,
    EPS_LOG, EPS_DIV, EPS_NORM,
)
from solver.crystals import Entangler, compose

torch.set_grad_enabled(False)


# ═══════════════════════════════════════════════════════════════
#  Test 1: Relaxation of Born(θ) to 1/27
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("TEST 1: BORN(θ) RELAXATION — the mass gap as a rate")
print("=" * 70)

# Prepare a consistent evidence stream
_random.seed(42)
evidence_stream = []
for _ in range(200):
    hi = _random.randint(0, H - 1)
    ev = torch.zeros(MASS_DIM, dtype=torch.cfloat)
    ev[hi] = 0.6 + 0j
    for k in range(H):
        if k != hi:
            ev[k] = 0.05 + 0j
    ev[H] = 1.0 - ev[:H].sum()
    evidence_stream.append(ev)

# Start with different Born(θ) values
print(f"\n  Born floor = {BORN_FLOOR:.6f} = 1/{H**3}")
print(f"\n  Relaxation from various initial Born(θ):")

for initial_theta_born in [0.50, 0.30, 0.15, 0.08, 0.05, 0.037037, 0.02, 0.01]:
    # Construct mass with this Born(θ)
    # Born(θ) = |m_θ|² / Σ|m_i|²
    # Set singletons equal, adjust θ
    # If Born(θ) = p, and Born(h_i) = (1-p)/3 each:
    # |m_θ|² / (3|m_h|² + |m_θ|²) = p
    # Let m_h = sqrt((1-p)/3), m_θ = sqrt(p) (real masses)
    p = initial_theta_born
    m_h = math.sqrt((1 - p) / 3)
    m_theta = math.sqrt(p)
    mass = torch.tensor([m_h, m_h, m_h, m_theta], dtype=torch.cfloat)

    # Verify
    bp = born_probabilities(mass)
    assert abs(bp[3].item() - p) < 0.001, f"Born(θ) = {bp[3].item()}, expected {p}"

    # Apply evidence and track Born(θ)
    thetas = [bp[3].item()]
    for step in range(150):
        ev = evidence_stream[step]
        disc_ev = torch.zeros_like(ev)
        disc_ev[:H] = ev[:H] * 0.3
        disc_ev[H] = 1.0 - disc_ev[:H].sum()
        mass, K = ds_combine(mass, disc_ev)
        bp = born_probabilities(mass)
        thetas.append(bp[3].item())

    # Report relaxation
    print(f"\n  Initial Born(θ) = {initial_theta_born:.4f}:")
    for step in [0, 1, 2, 3, 5, 10, 20, 50, 100, 150]:
        if step < len(thetas):
            excess = thetas[step] - BORN_FLOOR
            ln_excess = math.log(abs(excess)) if abs(excess) > 1e-12 else -99
            print(f"    step {step:3d}: Born(θ) = {thetas[step]:.6f}, "
                  f"excess = {excess:+.6f}, ln|excess| = {ln_excess:.4f}")

    # Compute relaxation rate from steps 1-20
    rates = []
    for i in range(1, min(20, len(thetas))):
        ex_prev = abs(thetas[i-1] - BORN_FLOOR)
        ex_curr = abs(thetas[i] - BORN_FLOOR)
        if ex_prev > 1e-12 and ex_curr > 1e-12:
            rate = math.log(ex_prev / ex_curr)
            rates.append(rate)
    if rates:
        avg_rate = sum(rates) / len(rates)
        print(f"    Relaxation rate (steps 1-20): {avg_rate:.4f}")
        print(f"    Δ = {DELTA:.4f}, ratio = {avg_rate/DELTA:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 2: Relaxation from BELOW the floor
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 2: RELAXATION FROM BELOW THE BORN FLOOR")
print("What happens if Born(θ) starts BELOW 1/27?")
print("=" * 70)

for initial_theta_born in [0.03, 0.02, 0.01, 0.001]:
    p = initial_theta_born
    m_h = math.sqrt((1 - p) / 3)
    m_theta = math.sqrt(p)
    mass = torch.tensor([m_h, m_h, m_h, m_theta], dtype=torch.cfloat)

    thetas = [born_probabilities(mass)[3].item()]
    for step in range(50):
        ev = evidence_stream[step]
        disc_ev = torch.zeros_like(ev)
        disc_ev[:H] = ev[:H] * 0.3
        disc_ev[H] = 1.0 - disc_ev[:H].sum()
        mass, K = ds_combine(mass, disc_ev)
        thetas.append(born_probabilities(mass)[3].item())

    print(f"\n  Initial Born(θ) = {initial_theta_born:.4f} (below floor {BORN_FLOOR:.4f}):")
    for step in [0, 1, 2, 3, 5, 10, 20, 50]:
        if step < len(thetas):
            print(f"    step {step:3d}: Born(θ) = {thetas[step]:.6f}")


# ═══════════════════════════════════════════════════════════════
#  Test 3: The compositional trace limit
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 3: TRACE UNDER ITERATED COMPOSITION")
print(f"Does trace → 1/(H+1) = {1/(H+1):.6f} exactly?")
print("=" * 70)

id_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

# Test with several different starting crystals
for name, corr in [
    ("identity", torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)),
    ("inverse", torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)),
    ("modular", torch.tensor([[.4,.3,.3],[.3,.4,.3],[.3,.3,.4]], dtype=torch.float32)),
    ("exponential", torch.tensor([[1,0,0],[.6,.4,0],[0,0,1]], dtype=torch.float32)),
]:
    ent = Entangler(corr).build()
    current = ent.joint.clone()

    print(f"\n  {name} crystal:")
    for step in range(15):
        tr = sum(current[i, i] for i in range(MASS_DIM))
        tr_real = tr.real.item()
        diff_from_quarter = tr_real - 1/(H+1)
        print(f"    comp {step:2d}: trace = {tr_real:.6f}, "
              f"diff from 1/(H+1) = {diff_from_quarter:+.6f}")
        current = compose(current, ent.joint)


# ═══════════════════════════════════════════════════════════════
#  Test 4: The fixed-point manifold structure
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 4: STRUCTURE OF THE FIXED-POINT MANIFOLD")
print("Born(θ) = 1/27 is fixed. What's free? What's constrained?")
print("=" * 70)

# Run DS convergence with different evidence BIASES
# (some seeds favor h0, others h1, others h2)
print("\n  Final states from different evidence biases:")
print(f"  {'bias':>10s} {'Born(h0)':>10s} {'Born(h1)':>10s} {'Born(h2)':>10s} {'Born(θ)':>10s} {'entropy':>8s}")

for bias_h in range(H):
    # Evidence biased toward hypothesis bias_h
    _random.seed(42)
    biased_evidence = []
    for _ in range(200):
        ev = torch.zeros(MASS_DIM, dtype=torch.cfloat)
        ev[bias_h] = 0.7 + 0j
        for k in range(H):
            if k != bias_h:
                ev[k] = 0.05 + 0j
        ev[H] = 1.0 - ev[:H].sum()
        biased_evidence.append(ev)

    mass = torch.tensor([0.25, 0.25, 0.25, 0.25], dtype=torch.cfloat)
    for ev in biased_evidence:
        disc_ev = torch.zeros_like(ev)
        disc_ev[:H] = ev[:H] * 0.3
        disc_ev[H] = 1.0 - disc_ev[:H].sum()
        mass, K = ds_combine(mass, disc_ev)

    bp = born_probabilities(mass)
    ent_val = born_entropy(mass).item()
    print(f"  h{bias_h} biased {bp[0].item():>10.6f} {bp[1].item():>10.6f} "
          f"{bp[2].item():>10.6f} {bp[3].item():>10.6f} {ent_val:>8.4f}")

# Mixed evidence (uniform across hypotheses)
mass = torch.tensor([0.25, 0.25, 0.25, 0.25], dtype=torch.cfloat)
_random.seed(42)
for step in range(200):
    hi = _random.randint(0, H - 1)
    ev = torch.zeros(MASS_DIM, dtype=torch.cfloat)
    ev[hi] = 0.6
    for k in range(H):
        if k != hi: ev[k] = 0.05
    ev[H] = 1.0 - ev[:H].sum()
    disc_ev = torch.zeros_like(ev)
    disc_ev[:H] = ev[:H] * 0.3
    disc_ev[H] = 1.0 - disc_ev[:H].sum()
    mass, K = ds_combine(mass, disc_ev)

bp = born_probabilities(mass)
ent_val = born_entropy(mass).item()
print(f"  mixed    {bp[0].item():>10.6f} {bp[1].item():>10.6f} "
      f"{bp[2].item():>10.6f} {bp[3].item():>10.6f} {ent_val:>8.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 5: Born(θ) under the FULL crystal dynamics
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 5: BORN(θ) IN THE JOINT MASS — does it also hit 1/27?")
print("The previous tests used 4-dim masses. The entangler builds")
print("16-dim joint masses. Is Born(θ,θ) also at the floor?")
print("=" * 70)

# Build crystals and check Born(θ,θ) = |joint[3,3]|² / total
for name, corr in [
    ("identity", torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)),
    ("inverse", torch.tensor([[0,0,1],[0,1,0],[1,0,0]], dtype=torch.float32)),
    ("modular", torch.tensor([[.4,.3,.3],[.3,.4,.3],[.3,.3,.4]], dtype=torch.float32)),
]:
    ent = Entangler(corr).build()
    total = ent.joint.abs().pow(2).sum().item()
    theta_theta = ent.joint[3, 3].abs().pow(2).item() / total

    # Marginal Born(θ) — sum over one index
    row_born = ent.joint[3, :].abs().pow(2).sum().item() / total
    col_born = ent.joint[:, 3].abs().pow(2).sum().item() / total

    # What would Born(θ,θ) be if marginals were independent?
    independent_theta_theta = row_born * col_born / (row_born + col_born) if (row_born + col_born) > 0 else 0

    print(f"\n  {name}:")
    print(f"    Born(θ,θ) = {theta_theta:.6f}")
    print(f"    1/H³ = {BORN_FLOOR:.6f}")
    print(f"    1/H⁶ = {1/H**6:.6f}")
    print(f"    Born(θ row) = {row_born:.6f}")
    print(f"    Born(θ col) = {col_born:.6f}")
    print(f"    (1/H³)² = {BORN_FLOOR**2:.6f}")

    # Check: is Born(θ row) = Born(θ)?
    # For the 4-dim marginal: sum over columns
    marginal = ent.joint.abs().pow(2).sum(dim=1)
    marginal_born = marginal / marginal.sum()
    print(f"    4-dim marginal Born = [{marginal_born[0]:.4f}, {marginal_born[1]:.4f}, "
          f"{marginal_born[2]:.4f}, {marginal_born[3]:.4f}]")
    print(f"    Marginal Born(θ) = {marginal_born[3].item():.6f}")


# ═══════════════════════════════════════════════════════════════
#  Test 6: The relaxation rate in the joint space
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 6: RELAXATION IN THE JOINT SPACE")
print("Build an entangler trajectory. Track marginal Born(θ).")
print("Does it converge to 1/27? At what rate?")
print("=" * 70)

_random.seed(42)
joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
joint[-1, -1] = 1.0 + 0j  # Start as pure (θ,θ) — Born(θ) = 1.0

flat = id_corr.reshape(-1)
flat = flat / flat.sum().clamp(min=EPS_LOG)

thetas_joint = []
for step in range(80):
    # Measure marginal Born(θ)
    marginal = joint.abs().pow(2).sum(dim=1)
    marginal_born = marginal / marginal.sum().clamp(min=1e-10)
    thetas_joint.append(marginal_born[3].item())

    # Entangle one step
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
    ev_s = joint_ev[:jH] * 0.3; ev_t = 1.0 + 0j - ev_s.sum()
    disc_ev = torch.cat([ev_s, ev_t.unsqueeze(0)])
    curr = joint.reshape(joint_dim)
    cs1, ct1 = curr[:jH], curr[jH:]
    cs2, ct2 = disc_ev[:jH], disc_ev[jH:]
    K_val = cs1.sum() * cs2.sum() - (cs1 * cs2).sum()
    comb_s = cs1 * cs2 + cs1 * ct2 + ct1 * cs2
    comb_t = ct1 * ct2
    nv = 1.0 - K_val
    snv = nv if abs(nv) > EPS_NORM else torch.tensor(EPS_NORM + 0j)
    upd = torch.cat([comb_s, comb_t]) / snv
    rs = upd.real.sum()
    if abs(rs) > EPS_DIV: upd = upd / rs
    joint = upd.reshape(MASS_DIM, MASS_DIM)

print(f"\n  Marginal Born(θ) during entangling (starting from pure ignorance):")
for step in [0, 1, 2, 3, 5, 10, 20, 30, 40, 50, 60, 70, 79]:
    if step < len(thetas_joint):
        excess = thetas_joint[step] - BORN_FLOOR
        print(f"    step {step:2d}: Born(θ) = {thetas_joint[step]:.6f}, "
              f"excess = {excess:+.6f}")

# Relaxation rate in joint space
rates_j = []
for i in range(5, min(40, len(thetas_joint))):
    ex_prev = abs(thetas_joint[i-1] - BORN_FLOOR)
    ex_curr = abs(thetas_joint[i] - BORN_FLOOR)
    if ex_prev > 1e-6 and ex_curr > 1e-6:
        rate = math.log(ex_prev / ex_curr)
        rates_j.append(rate)

if rates_j:
    avg_rate_j = sum(rates_j) / len(rates_j)
    print(f"\n  Joint-space relaxation rate (steps 5-40): {avg_rate_j:.4f}")
    print(f"  Δ = {DELTA:.4f}, ratio = {avg_rate_j/DELTA:.3f}")
    print(f"  Δ/4 = {DELTA/4:.4f} (growth rate)")


print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
