"""
Connection between the crystal spectral gap and the paper's mass gap.

Crystal spectral gap: ln(H)/2 = ln(3)/2 ≈ 0.5493
  (eigenvalue ratio |λ₁/λ₀| = 1/√H of seed-averaged crystal)

Paper's mass gap: Δ_gap = -ln(λ₀_Jacobian) = 1.2626
  (second eigenvalue of DS iteration Jacobian at K* fixed point)

Conjecture: Δ_gap = (H/2) × ln(H) × (1-K*)
  = (23/20) × ln(3)
  = crystal_gap × H × (1-K*)
  = 1.2634...

Does it match to high precision?
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import math
from solver.algebra import (H, MASS_DIM, DELTA, K_STAR, BORN_FLOOR,
                            ds_combine, enforce_born_floor, born_probabilities,
                            EPS_LOG, EPS_NORM)

torch.set_grad_enabled(False)


# ═══════════════════════════════════════════════════════════════
#  Step 1: Reproduce the DS Jacobian eigenvalue
# ═══════════════════════════════════════════════════════════════

print("=" * 80)
print("REPRODUCING THE DS JACOBIAN EIGENVALUE AT K* FIXED POINT")
print("=" * 80)

# The DS iteration at the fixed point: m → DS(m, e)
# where e is the evidence mass and m is the current state
# At equilibrium with K*, the Jacobian has eigenvalue λ₀

# Method: find the fixed point mass m*, then compute Jacobian

# The fixed point satisfies: m* = DS(m*, e*) where e* produces K=K*
# From the memory: λ₀ = 0.2829102944

# Let's reproduce this via finite-difference Jacobian
# First, find the fixed point by iteration

def ds_step(m, evidence, floor=BORN_FLOOR):
    """One DS combination step."""
    combined, K = ds_combine(m, evidence, floor=floor)
    return combined, K


# Build evidence that produces K = K*
# For a mass m with Born probs p, and evidence e with probs q:
# K = sum_i(p_i) * sum_j(q_j) - sum_i(p_i * q_i)  [for singletons only]
# K* = 7/30

# Start from ignorance, iterate to find fixed point
m = torch.zeros(1, MASS_DIM, dtype=torch.cfloat)
m[0, -1] = 1.0 + 0j

# Use canonical evidence: concentrated on h₀ with Born floor
# Let's use the actual DS equation
# At fixed point, the mass m* has specific Born probabilities
# p_dom (dominant singleton), p_weak (other singletons), p_theta

# From the paper: at K*=7/30, the fixed point has:
# p₀ = (1-K*) × dominant_share, etc.
# But we need the actual iteration

# Simple approach: start with a specific evidence and iterate
def make_evidence_for_K(K_target, m_state):
    """Build an evidence mass that produces conflict K_target with m_state."""
    p_state = born_probabilities(m_state).squeeze()

    # We want: K = (sum_i p_i)(sum_j q_j) - sum_i(p_i q_j)
    # For singletons only (i, j < H)
    # Try: q concentrated on h₀
    ev = torch.zeros(1, MASS_DIM, dtype=torch.cfloat)
    ev[0, 0] = 0.9 + 0j
    ev[0, 1] = 0.03 + 0j
    ev[0, 2] = 0.03 + 0j
    ev[0, -1] = 1.0 - 0.96 + 0j  # = 0.04

    return ev


# Better: use power iteration on the DS map to find the fixed point
# Start with some initial mass and iterate DS with self-evidence

print(f"\n--- Finding DS fixed point by iteration ---")

# The fixed point from the paper has K=K* in steady state
# Let's build it directly

# At the fixed point with K*=7/30:
# The evidence and state must be self-consistent
# For the SELF-COMBINATION version: m -> DS(m, m)

# Actually, from the paper: the mass gap is measured from EXTERNAL evidence
# Self-combination gives K*=0. The K*=7/30 comes from gauge field evidence.

# Let's compute the Jacobian of m -> DS(m, e) where e is FIXED
# at a canonical value.

# From paper_spectral_gap_findings: "K*=7/30 equilibrium requires EXTERNAL evidence"
# The DS iteration is m_{n+1} = DS(m_n, e) for fixed e

# What evidence e produces K=K* at the fixed point?
# At the fixed point, the combined mass equals the input mass:
# m* = DS(m*, e) = (m* ⊗ e) / (1 - K(m*, e))

# Let's search for the fixed point numerically
def ds_iterate(evidence, n_steps=200):
    """Iterate m -> DS(m, evidence) from ignorance."""
    m = torch.zeros(1, MASS_DIM, dtype=torch.cfloat)
    m[0, -1] = 1.0 + 0j

    for step in range(n_steps):
        m_new, K = ds_combine(m, evidence)
        m = m_new

    return m, K


# Try various evidence shapes and find one giving K* at fixed point
print(f"\n  Searching for evidence that gives K=K* at fixed point...")

best_ev = None
best_diff = 1.0

for p_dom in [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97]:
    ev = torch.zeros(1, MASS_DIM, dtype=torch.cfloat)
    p_rest = (1.0 - p_dom) / H  # split remainder equally
    for i in range(H):
        ev[0, i] = (p_dom if i == 0 else p_rest) + 0j
    ev[0, -1] = 1.0 - ev[0, :H].real.sum() + 0j

    m_fp, K_fp = ds_iterate(ev, n_steps=500)
    K_val = K_fp.abs().item()
    diff = abs(K_val - K_STAR)
    if diff < best_diff:
        best_diff = diff
        best_ev = ev.clone()
    # print(f"    p_dom={p_dom:.2f}: K={K_val:.6f} (diff={diff:.4e})")


# Fine-tune with bisection
lo_p, hi_p = 0.5, 0.99
for _ in range(50):
    mid_p = (lo_p + hi_p) / 2
    ev = torch.zeros(1, MASS_DIM, dtype=torch.cfloat)
    p_rest = (1.0 - mid_p) / H
    for i in range(H):
        ev[0, i] = (mid_p if i == 0 else p_rest) + 0j
    ev[0, -1] = 1.0 - ev[0, :H].real.sum() + 0j

    m_fp, K_fp = ds_iterate(ev, n_steps=500)
    K_val = K_fp.abs().item()

    if K_val < K_STAR:
        lo_p = mid_p
    else:
        hi_p = mid_p

ev_star = torch.zeros(1, MASS_DIM, dtype=torch.cfloat)
p_dom_star = (lo_p + hi_p) / 2
p_rest_star = (1.0 - p_dom_star) / H
for i in range(H):
    ev_star[0, i] = (p_dom_star if i == 0 else p_rest_star) + 0j
ev_star[0, -1] = 1.0 - ev_star[0, :H].real.sum() + 0j

m_star, K_star = ds_iterate(ev_star, n_steps=1000)
K_val = K_star.abs().item()

print(f"\n  Fixed point found:")
print(f"    Evidence p_dom = {p_dom_star:.10f}")
print(f"    K at fixed point = {K_val:.10f}")
print(f"    K* (target)      = {K_STAR:.10f}")
print(f"    Difference       = {abs(K_val - K_STAR):.2e}")

print(f"\n  Fixed point mass m*:")
p_star = born_probabilities(m_star).squeeze()
for i in range(MASS_DIM):
    print(f"    p_{i} = {p_star[i].item():.10f}")


# ═══════════════════════════════════════════════════════════════
#  Step 2: Compute the Jacobian at the fixed point
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("JACOBIAN AT THE FIXED POINT")
print("="*80)

# The DS map F: m -> DS(m, e*) is a map from R^{H} -> R^{H}
# (in the L1=1 tangent space, so effectively H-dimensional)
# Its Jacobian J = dF/dm has eigenvalues that determine stability

# Finite difference: J_{ij} ≈ (F(m* + ε·e_j) - F(m*))_i / ε

eps = 1e-7
dim = H  # working in the H-dimensional subspace (theta determined by L1=1)

J = torch.zeros(dim, dim)

for j in range(dim):
    # Perturb m* in direction j (keep L1=1 by adjusting theta)
    m_plus = m_star.clone()
    m_plus[0, j] += eps
    m_plus[0, -1] -= eps  # maintain L1=1

    m_minus = m_star.clone()
    m_minus[0, j] -= eps
    m_minus[0, -1] += eps

    # Apply DS map
    Fm_plus, _ = ds_combine(m_plus, ev_star)
    Fm_minus, _ = ds_combine(m_minus, ev_star)

    # Jacobian column
    for i in range(dim):
        J[i, j] = (Fm_plus[0, i].real - Fm_minus[0, i].real) / (2 * eps)

print(f"\n  Jacobian matrix (3×3, L1=1 tangent space):")
for i in range(dim):
    row = "  ".join(f"{J[i,j].item():>10.6f}" for j in range(dim))
    print(f"    [{row}]")

eigvals_J = torch.linalg.eigvals(J)
eigvals_J_abs = eigvals_J.abs()
eigvals_J_sorted, _ = eigvals_J_abs.sort(descending=True)

print(f"\n  Jacobian eigenvalues:")
for i in range(dim):
    ev = eigvals_J[i]
    print(f"    λ_{i} = {ev.real.item():.10f} + {ev.imag.item():.10f}i  (|λ| = {eigvals_J_abs[i].item():.10f})")

lambda_0 = eigvals_J_sorted[0].item()
delta_gap = -math.log(lambda_0)

print(f"\n  Dominant eigenvalue: λ₀ = {lambda_0:.10f}")
print(f"  Mass gap: Δ_gap = -ln(λ₀) = {delta_gap:.10f}")


# ═══════════════════════════════════════════════════════════════
#  Step 3: Test the connection formula
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("THE CONNECTION FORMULA")
print("="*80)

crystal_gap = math.log(H) / 2  # ln(3)/2
formula_1 = crystal_gap * H * (1 - K_STAR)  # = (23/20) × ln(3)
formula_2 = (23/20) * math.log(3)  # same thing, explicit

print(f"\n  Crystal spectral gap = ln(H)/2 = {crystal_gap:.10f}")
print(f"  H × (1-K*) = {H * (1-K_STAR):.10f} = 23/10")
print(f"\n  Formula: Δ_gap = crystal_gap × H × (1-K*)")
print(f"           = (ln(H)/2) × H × (1-K*)")
print(f"           = (H/2) × ln(H) × (1-K*)")
print(f"           = (23/20) × ln(3)")
print(f"           = {formula_1:.10f}")
print(f"\n  Measured: Δ_gap = {delta_gap:.10f}")
print(f"  Difference: {formula_1 - delta_gap:.2e}")
print(f"  Relative: {abs(formula_1 - delta_gap)/delta_gap:.2e}")

# Check against the memory's value
memory_lambda_0 = 0.2829102944
memory_delta = -math.log(memory_lambda_0)
print(f"\n  Paper's value (from memory): λ₀ = {memory_lambda_0}")
print(f"  Paper's Δ_gap = -ln({memory_lambda_0}) = {memory_delta:.10f}")
print(f"  Formula value = {formula_1:.10f}")
print(f"  Difference from paper: {formula_1 - memory_delta:.6e}")
print(f"  Relative: {abs(formula_1 - memory_delta)/memory_delta:.6e}")


# ═══════════════════════════════════════════════════════════════
#  Step 4: What would make it exact?
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("EXACTNESS CHECK")
print("="*80)

# If Δ_gap = (23/20) × ln(3) exactly, then:
# λ₀ = exp(-23/20 × ln(3)) = 3^{-23/20} = 3^{-1.15}
exact_lambda = 3 ** (-23/20)
print(f"\n  If exact: λ₀ = 3^{{-23/20}} = {exact_lambda:.10f}")
print(f"  Paper's:  λ₀ = {memory_lambda_0:.10f}")
print(f"  Our computation: λ₀ = {math.exp(-delta_gap):.10f}")
print(f"  Difference (exact vs paper): {exact_lambda - memory_lambda_0:.6e}")
print(f"  Relative: {abs(exact_lambda - memory_lambda_0)/memory_lambda_0:.6e}")

# If not exact, what IS the exact relationship?
# The ratio Δ_gap / (ln(3)/2) should equal H(1-K*) = 23/10 = 2.3
measured_ratio = delta_gap / crystal_gap
print(f"\n  Δ_gap / crystal_gap = {measured_ratio:.10f}")
print(f"  H(1-K*) = 23/10 = {H*(1-K_STAR):.10f}")
print(f"  Difference: {measured_ratio - H*(1-K_STAR):.6e}")

# Also check against other candidates
print(f"\n  Other candidates for the ratio:")
candidates = {
    "23/10": 23/10,
    "13/5 × (23/26)": 13/5 * 23/26,
    "H²/(H-1)×(1-K*)": H**2/(H-1) * (1-K_STAR),
    "2(H²+1)/H²×(1-K*)": 2*(H**2+1)/H**2 * (1-K_STAR),
    "H/(1-BORN)": H / (1-BORN_FLOOR),
    "H×(H-1+K*)/(H-1)": H * (H-1+K_STAR)/(H-1),
}

for name, val in sorted(candidates.items(), key=lambda x: abs(x[1] - measured_ratio)):
    diff = measured_ratio - val
    print(f"    {name:>30s} = {val:.10f} (diff = {diff:+.6e})")


# ═══════════════════════════════════════════════════════════════
#  Step 5: The physical meaning
# ═══════════════════════════════════════════════════════════════

print(f"\n\n{'='*80}")
print("PHYSICAL MEANING")
print("="*80)

print(f"""
  If Δ_gap = (H/2) × ln(H) × (1-K*) is exact (or approximately so):

  The mass gap decomposes as:

    Δ_gap = (crystal spectral gap) × (effective dimension)

  where:
    crystal spectral gap = ln(H)/2 = rate of rank reduction in a single crystal
    effective dimension  = H(1-K*) = # hypotheses × non-conflicting fraction
                         = 23/10 at H=3

  Physical interpretation:
    A single crystal loses structure at rate ln(H)/2 per composition.
    The DS iteration involves H(1-K*) effective channels.
    The TOTAL decay rate = single-channel rate × channel count.

  This is analogous to:
    decay rate of coupled system = single-mode rate × number of modes

  The mass gap is the TOTAL decorrelation rate of the gauge field,
  summing over all non-conflicting channels.
""")


print(f"\n{'='*80}")
print("WHAT THE GAP CONNECTION REVEALS")
print("="*80)
