"""
The entangler's characteristic conflict K ≈ 0.207.

The phase transition analysis showed that the step-by-step conflict K
is CONSTANT across correlation matrices: K ≈ 0.207 regardless of
whether the input is identity, inverse, modular, or mixed.

Is this a fundamental constant? What determines it?

Candidates:
  (H-1)/H²     = 2/9   ≈ 0.2222
  K*            = 7/30  ≈ 0.2333
  2/H³          = 2/27  ≈ 0.0741
  (H-1)/(H²+1) = 2/10  = 0.2000

The answer should come from the entangler's construction: discount=0.3,
evidence strength=0.55 for the dominant + 0.05 for the rest + 0.06j phase.

Let's derive it from first principles.
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
import random as _random
from solver.algebra import H, MASS_DIM, K_STAR, BORN_FLOOR, EPS_LOG, EPS_DIV, EPS_NORM
from solver.crystals import Entangler, schmidt_number

torch.set_grad_enabled(False)


# ═══════════════════════════════════════════════════════════════
#  Test 1: K as a function of discount
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("TEST 1: K vs DISCOUNT PARAMETER")
print("How does the characteristic K depend on the entangler's discount?")
print("=" * 70)

id_corr = torch.tensor([[1,0,0],[0,1,0],[0,0,1]], dtype=torch.float32)

for discount in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    # Run entangler with this discount, collect K values
    all_Ks = []
    for seed in range(20):
        _random.seed(seed)

        joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
        joint[-1, -1] = 1.0 + 0j

        flat = id_corr.reshape(-1)
        flat = flat / flat.sum().clamp(min=EPS_LOG)

        for step in range(50):
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

            K = cs1.sum() * cs2.sum() - (cs1 * cs2).sum()
            all_Ks.append(abs(K.item()) if isinstance(K, torch.Tensor) else abs(K))

            comb_s = cs1 * cs2 + cs1 * ct2 + ct1 * cs2
            comb_t = ct1 * ct2
            norm = 1.0 - K
            safe_norm = norm if abs(norm) > EPS_NORM else torch.tensor(EPS_NORM + 0j)
            updated = torch.cat([comb_s, comb_t]) / safe_norm
            re_sum = updated.real.sum()
            if abs(re_sum) > EPS_DIV:
                updated = updated / re_sum
            joint = updated.reshape(MASS_DIM, MASS_DIM)

    avg_K = sum(all_Ks) / len(all_Ks)
    sn = schmidt_number(joint)
    print(f"  discount={discount:.1f}: avg K = {avg_K:.4f}, final Schmidt = {sn:.3f}")


# ═══════════════════════════════════════════════════════════════
#  Test 2: K in the early vs late steps
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 2: K TRAJECTORY (early vs late steps)")
print("=" * 70)

# Run one entangler, record K at each step
_random.seed(42)
joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
joint[-1, -1] = 1.0 + 0j

flat = id_corr.reshape(-1)
flat = flat / flat.sum().clamp(min=EPS_LOG)

K_trajectory = []
for step in range(200):
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
    ev_s = joint_ev[:jH] * 0.3  # discount
    ev_t = 1.0 + 0j - ev_s.sum()
    disc_ev = torch.cat([ev_s, ev_t.unsqueeze(0)])

    curr = joint.reshape(joint_dim)
    cs1, ct1 = curr[:jH], curr[jH:]
    cs2, ct2 = disc_ev[:jH], disc_ev[jH:]
    K = cs1.sum() * cs2.sum() - (cs1 * cs2).sum()
    K_val = abs(K.item()) if isinstance(K, torch.Tensor) else abs(K)
    K_trajectory.append(K_val)

    comb_s = cs1 * cs2 + cs1 * ct2 + ct1 * cs2
    comb_t = ct1 * ct2
    norm = 1.0 - K
    safe_norm = norm if abs(norm) > EPS_NORM else torch.tensor(EPS_NORM + 0j)
    updated = torch.cat([comb_s, comb_t]) / safe_norm
    re_sum = updated.real.sum()
    if abs(re_sum) > EPS_DIV:
        updated = updated / re_sum
    joint = updated.reshape(MASS_DIM, MASS_DIM)

# Report K in windows
windows = [(0,10), (10,20), (20,50), (50,100), (100,200)]
for lo, hi in windows:
    segment = K_trajectory[lo:hi]
    avg = sum(segment) / len(segment)
    print(f"  Steps {lo:3d}-{hi:3d}: avg K = {avg:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 3: Analytical K from evidence structure
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 3: DERIVE K FROM EVIDENCE STRUCTURE")
print("=" * 70)

# The evidence mass for one step (after discount):
# If hi=0 is selected:
#   ev = [0.55+0.06j, 0.05+noise_j, 0.05+noise_j, θ] * discount
# After discounting by 0.3:
#   singletons: [0.165+0.018j, 0.015+noise_j, 0.015+noise_j]
#   θ = 1 - sum(singletons)

# The JOINT evidence is outer product ev_a ⊗ ev_b, then discounted.
# For identity correlation with (hi,hj) = (i,i):
#   ev_a = ev_b concentrated on hypothesis i

# Let's compute the expected K analytically.
# For the first step (current = ignorance):
#   cs1 = [0,...,0] (all zero singletons)
#   cs2 = discounted evidence singletons
#   K = cs1.sum() * cs2.sum() - (cs1*cs2).sum() = 0 (because cs1 = 0)
# So K=0 at step 1 — confirmed by trajectory.

# For subsequent steps, K depends on the current state.
# In the generic regime (after ~10 steps), the singletons are spread
# across the 15 joint positions. Let's compute K for a "generic" state.

# A key insight: K for DS combination of mass m and evidence e is
#   K = Σ_{i≠j} m(h_i) * e(h_j) = m_total * e_total - Σ_i m(h_i)*e(h_i)
# where m_total = sum of all singletons, e_total = sum of all singletons

# For the discounted evidence (discount=α, dominant strength=s, minor=r):
#   In 16-dim joint space, the evidence has:
#   One dominant at position (hi,hj): α * s² (where s = 0.55+0.06j)
#   Cross terms: various α * s*r products
#   θ-cross terms: α * s * θ_component
#   All minor: α * r² terms

# This is complex. Let me just compute it numerically.

# Average evidence singleton sum and agreement for each evidence type
print("\n  Numerical evidence structure:")

# Generate many evidence samples and compute their properties
n_samples = 10000
_random.seed(42)
sum_s_total = 0
sum_agreement = 0

for _ in range(n_samples):
    hi = _random.randint(0, H-1)
    hj = _random.randint(0, H-1)

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

    joint_ev = torch.einsum("m,n->mn", ev_a, ev_b).reshape(MASS_DIM**2)
    jH = MASS_DIM**2 - 1
    ev_s = joint_ev[:jH] * 0.3  # discount
    ev_t = 1.0 + 0j - ev_s.sum()

    s_total = ev_s.abs().sum().item()
    sum_s_total += s_total

avg_s_total = sum_s_total / n_samples
print(f"  Average evidence singleton total (after discount): {avg_s_total:.6f}")
print(f"  Average evidence θ: {1 - avg_s_total:.6f}")

# For a "typical" mass after ~20 steps, the singleton sum stabilizes.
# Let's measure it directly.
_random.seed(42)
joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
joint[-1, -1] = 1.0 + 0j
flat = id_corr.reshape(-1)
flat = flat / flat.sum().clamp(min=EPS_LOG)

mass_s_totals = []
for step in range(100):
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

    joint_dim = MASS_DIM**2; jH = joint_dim - 1
    joint_ev = torch.einsum("m,n->mn", ev_a, ev_b).reshape(joint_dim)
    ev_s = joint_ev[:jH] * 0.3; ev_t = 1.0 + 0j - ev_s.sum()
    disc_ev = torch.cat([ev_s, ev_t.unsqueeze(0)])
    curr = joint.reshape(joint_dim)
    cs1, ct1 = curr[:jH], curr[jH:]
    cs2, ct2 = disc_ev[:jH], disc_ev[jH:]
    K = cs1.sum() * cs2.sum() - (cs1 * cs2).sum()
    comb_s = cs1 * cs2 + cs1 * ct2 + ct1 * cs2
    comb_t = ct1 * ct2
    norm = 1.0 - K
    safe_norm = norm if abs(norm) > EPS_NORM else torch.tensor(EPS_NORM + 0j)
    updated = torch.cat([comb_s, comb_t]) / safe_norm
    re_sum = updated.real.sum()
    if abs(re_sum) > EPS_DIV: updated = updated / re_sum
    joint = updated.reshape(MASS_DIM, MASS_DIM)

    curr_s_total = joint.reshape(-1)[:jH].abs().sum().item()
    mass_s_totals.append(curr_s_total)

# Report mass singleton total over time
for lo, hi_val in [(0,5), (5,10), (10,20), (20,50), (50,100)]:
    segment = mass_s_totals[lo:hi_val]
    avg = sum(segment) / len(segment)
    print(f"  Steps {lo:3d}-{hi_val:3d}: avg mass singleton total = {avg:.4f}")


# ═══════════════════════════════════════════════════════════════
#  Test 4: K from different evidence strengths
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 4: K vs EVIDENCE STRENGTH")
print("Change the dominant singleton from 0.55 to other values")
print("=" * 70)

for dominant in [0.40, 0.45, 0.50, 0.55, 0.60, 0.70, 0.80, 0.90]:
    minor = (1.0 - dominant) / (H - 1)  # split remaining equally
    all_Ks = []
    for seed in range(10):
        _random.seed(seed)
        joint = torch.zeros(MASS_DIM, MASS_DIM, dtype=torch.cfloat)
        joint[-1, -1] = 1.0 + 0j
        flat = id_corr.reshape(-1); flat = flat / flat.sum().clamp(min=EPS_LOG)

        for step in range(50):
            idx = _random.choices(range(H * H), weights=flat.tolist())[0]
            hi, hj = idx // H, idx % H
            ev_a = torch.zeros(MASS_DIM, dtype=torch.cfloat)
            ev_b = torch.zeros(MASS_DIM, dtype=torch.cfloat)
            ev_a[hi] = dominant + 0.06j; ev_b[hj] = dominant + 0.06j
            for k in range(H):
                noise = _random.gauss(0, 0.02)
                if k != hi: ev_a[k] = minor + noise * 1j
                if k != hj: ev_b[k] = minor + noise * 1j
            ev_a[H] = (1.0 - ev_a[:H].real.sum()) - ev_a[:H].imag.sum() * 1j
            ev_b[H] = (1.0 - ev_b[:H].real.sum()) - ev_b[:H].imag.sum() * 1j

            joint_dim = MASS_DIM**2; jH = joint_dim - 1
            joint_ev = torch.einsum("m,n->mn", ev_a, ev_b).reshape(joint_dim)
            ev_s = joint_ev[:jH] * 0.3; ev_t = 1.0 + 0j - ev_s.sum()
            disc_ev = torch.cat([ev_s, ev_t.unsqueeze(0)])
            curr = joint.reshape(joint_dim)
            cs1, ct1 = curr[:jH], curr[jH:]
            cs2, ct2 = disc_ev[:jH], disc_ev[jH:]
            K = cs1.sum() * cs2.sum() - (cs1 * cs2).sum()
            all_Ks.append(abs(K.item()) if isinstance(K, torch.Tensor) else abs(K))
            comb_s = cs1 * cs2 + cs1 * ct2 + ct1 * cs2
            comb_t = ct1 * ct2
            norm = 1.0 - K
            safe_norm = norm if abs(norm) > EPS_NORM else torch.tensor(EPS_NORM + 0j)
            updated = torch.cat([comb_s, comb_t]) / safe_norm
            re_sum = updated.real.sum()
            if abs(re_sum) > EPS_DIV: updated = updated / re_sum
            joint = updated.reshape(MASS_DIM, MASS_DIM)

    avg_K = sum(all_Ks) / len(all_Ks)
    sn = schmidt_number(joint)
    # Skip first 5 steps (warmup)
    late_Ks = all_Ks[5:]
    avg_K_late = sum(late_Ks) / len(late_Ks)
    print(f"  dominant={dominant:.2f}, minor={minor:.3f}: "
          f"avg K (all)={avg_K:.4f}, avg K (late)={avg_K_late:.4f}, "
          f"final Schmidt={sn:.3f}")

    # Theoretical: for concentrated evidence with singleton total T,
    # and mass with singleton total M, the expected K is:
    # K = M * T - (M*T)/N_eff where N_eff is the effective number of
    # agreeing positions. For uniform mass across 15 singletons:
    # agreement ≈ M*T/15, so K ≈ M*T * (1 - 1/15) = M*T * 14/15


# ═══════════════════════════════════════════════════════════════
#  Test 5: Derive the formula
# ═══════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("TEST 5: THE K FORMULA")
print("=" * 70)

# K = m1_total * m2_total - agreement
# where m1_total = sum of singletons in the current mass
# and m2_total = sum of singletons in the discounted evidence

# For a mass with N=15 singleton positions, if the mass is
# approximately uniform across these positions with total M:
#   m(h_i) ≈ M/15 for each singleton
#   agreement = Σ_i m(h_i)*e(h_i) ≈ M * e_dominant / 15 * (dominant_count)

# For the discounted evidence:
#   e_dominant = 0.3 * (0.55)² ≈ 0.3 * 0.3025 = 0.0908 (for the one position)
#   e_minor ≈ 0.3 * 0.55 * 0.05 = 0.00825 (for cross terms)
#   e_total ≈ 0.3 * [total of outer product singletons]

# The outer product ev_a ⊗ ev_b has:
#   (hi,hj) position: 0.55² + complex ≈ 0.3025
#   (hi, k≠hj): 0.55 * 0.05 ≈ 0.0275 (H-1 of these)
#   (k≠hi, hj): same
#   (k≠hi, k'≠hj): 0.05 * 0.05 ≈ 0.0025 ((H-1)² of these)
#   θ components: various

# Total singleton from outer product:
dom = 0.55
minor_val = 0.05
alpha = 0.3  # discount

# One evidence sample contributes:
s_dom_dom = dom * dom  # 1 position
s_dom_min = dom * minor_val  # 2*(H-1) positions
s_min_min = minor_val * minor_val  # (H-1)² positions
# Plus θ cross terms

n_dom_dom = 1
n_dom_min = 2 * (H - 1)
n_min_min = (H - 1) ** 2

total_singleton_undiscounted = (n_dom_dom * s_dom_dom +
                                 n_dom_min * s_dom_min +
                                 n_min_min * s_min_min)

# Total in the 15 singleton positions (out of 16 total joint positions)
# But we also have θ-row and θ-column positions:
# Positions: (0..2, 0..2) = 9 h-h positions, (0..2, 3) = 3 h-θ, (3, 0..2) = 3 θ-h, (3,3) = θ-θ
# Singletons are ALL except (3,3), so 15 positions.

# h-h positions (9): (hi,hj) is dominant, rest are cross/minor
# h-θ positions (3): ev_a[h] * ev_b[θ]
# θ-h positions (3): ev_a[θ] * ev_b[h]

# ev_b[θ] = 1 - ev_b[:H].real.sum() ≈ 1 - (0.55 + 2*0.05) = 0.35
ev_theta = 1 - dom - (H-1) * minor_val  # ≈ 0.35

s_h_theta = dom * ev_theta  # 1 position (hi, θ)
s_h_theta_minor = minor_val * ev_theta  # (H-1) positions
s_theta_h = ev_theta * dom  # 1 position (θ, hj)
s_theta_h_minor = ev_theta * minor_val  # (H-1) positions

total_15 = (n_dom_dom * s_dom_dom +
            n_dom_min * s_dom_min +
            n_min_min * s_min_min +
            1 * s_h_theta + (H-1) * s_h_theta_minor +
            1 * s_theta_h + (H-1) * s_theta_h_minor)

total_15_discounted = alpha * total_15
ev_theta_total = 1.0 - total_15_discounted

print(f"  Evidence structure (before discount):")
print(f"    Dominant (1 pos):     {s_dom_dom:.4f}")
print(f"    Dom×min ({n_dom_min} pos):  {s_dom_min:.4f} each, total {n_dom_min * s_dom_min:.4f}")
print(f"    Min×min ({n_min_min} pos):  {s_min_min:.4f} each, total {n_min_min * s_min_min:.4f}")
print(f"    H×θ (1+{H-1} pos):   {s_h_theta:.4f} + {(H-1)*s_h_theta_minor:.4f}")
print(f"    θ×H (1+{H-1} pos):   {s_theta_h:.4f} + {(H-1)*s_theta_h_minor:.4f}")
print(f"    Total 15 singletons:  {total_15:.4f}")
print(f"  After discount (α={alpha}):")
print(f"    Singleton total:      {total_15_discounted:.4f}")
print(f"    θ total:              {ev_theta_total:.4f}")

# For K computation:
# K = M * E - A where M = mass singleton total, E = evidence singleton total, A = agreement
# In steady state, M stabilizes (the mass reaches a characteristic distribution).
# Let's see what the steady-state M is:

print(f"\n  Measured mass singleton total in steady state:")
print(f"    Steps 20-50: {sum(mass_s_totals[20:50])/30:.4f}")
print(f"    Steps 50-100: {sum(mass_s_totals[50:100])/50:.4f}")

M_steady = sum(mass_s_totals[20:100]) / 80
E = total_15_discounted

# Agreement depends on how concentrated the mass is.
# For a mass concentrated on the SAME positions as the evidence:
# agreement ≈ M * dominant_fraction_of_evidence
# For a mass spread uniformly: agreement ≈ M * E / 15

# Uniform mass: K = M*E*(1 - 1/15) = M*E*14/15
K_uniform = M_steady * E * (1 - 1/15)
# Concentrated mass: K = M*E*(1 - agreement_fraction)
# where agreement_fraction depends on alignment

print(f"\n  Predicted K:")
print(f"    M (steady) = {M_steady:.4f}")
print(f"    E (per step) = {E:.4f}")
print(f"    K (uniform mass, 14/15) = {K_uniform:.4f}")
print(f"    K (observed) = 0.207")
print(f"    Ratio observed/predicted = {0.207 / K_uniform:.3f}" if K_uniform > 0 else "")

# What agreement fraction gives K = 0.207?
# K = M*E - A = M*E*(1 - A/(M*E))
# 0.207 = M*E*(1 - f) where f = agreement fraction
# f = 1 - 0.207/(M*E)
f_agreement = 1 - 0.207 / (M_steady * E) if M_steady * E > 0 else 0
print(f"    Agreement fraction that gives K=0.207: {f_agreement:.4f}")
print(f"    Equivalent to: 1 agreement position out of {1/f_agreement:.1f} total" if f_agreement > 0 else "")


print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
