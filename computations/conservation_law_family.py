"""
Characterise the family of DS fixed points parametrised by evidence.

For each evidence distribution, find the unique fixed point and record:
  - K* (conflict)
  - λ₀ (spectral gap)
  - The relationship between K* and evidence structure

Key question: does the self-consistency condition K*=7/30 have a
geometric interpretation within this family?
"""

import numpy as np
from scipy.optimize import fsolve

H = 3
FLOOR = 1.0 / H**3

def enforce_floor(m):
    s, theta = m[:3], m[3]
    ssq = sum(si**2 for si in s)
    total_sq = ssq + theta**2
    born = theta**2 / total_sq if total_sq > 0 else 1
    if born >= FLOOR:
        return m.copy()
    ss = sum(s)
    if ss < 1e-15:
        return m.copy()
    r = ssq / ss**2
    a_c = 26 - r
    b_c = 2 * r
    c_c = -r
    disc = b_c**2 - 4 * a_c * c_c
    t = (-b_c + np.sqrt(disc)) / (2 * a_c)
    alpha = (1 - t) / ss
    return np.array([s[0]*alpha, s[1]*alpha, s[2]*alpha, t])

def ds_combine_raw(m, e):
    s, theta = m[:3], m[3]
    se, phi = e[:3], e[3]
    s_new = s * se + s * phi + theta * se
    theta_new = theta * phi
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    denom = 1.0 - K
    m_out = np.zeros(4)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom
    total = np.sum(m_out)
    if total > 0:
        m_out = m_out / total
    return m_out, K

def full_step(m, e):
    m_ds, K = ds_combine_raw(m, e)
    return enforce_floor(m_ds), K

def born(m):
    return m[3]**2 / np.sum(m**2) if np.sum(m**2) > 0 else 1

def K_conflict(m, e):
    return sum(m[i] * e[j] for i in range(3) for j in range(3) if i != j)

def compute_jacobian(m, e, eps=1e-8):
    """3x3 projected Jacobian on L₁=1 tangent space."""
    # Tangent vectors: change s₁ (compensate θ), change s₂ (compensate θ)
    # At S₂ symmetric point, project to {δs₁, δs₂, δθ} with δs₁+2δs₂+δθ=0
    m0, _ = full_step(m, e)
    J = np.zeros((3, 3))  # 3x3 in (s₁, s₂, θ) with constraint s₁+2s₂+θ=1

    for j in range(3):
        idx = [0, 1, 3][j]  # s₁, s₂, θ
        m_p = m.copy()
        m_p[idx] += eps
        # Compensate: adjust θ (or s₁) to keep L₁=1
        if idx != 3:
            m_p[3] -= eps
        else:
            m_p[0] -= eps
        m_p_out, _ = full_step(m_p, e)
        for i in range(3):
            idx_i = [0, 1, 3][i]
            J[i, j] = (m_p_out[idx_i] - m0[idx_i]) / eps
    return J

# ============================================================
# Find fixed points across a fine grid of evidence parameters
# ============================================================
# Parametrise evidence by (e₁, e₂) with φ = 1 - e₁ - 2e₂
# S₂ symmetric: e₂ = e₃

print("=" * 80)
print("FAMILY OF FIXED POINTS: K* vs evidence structure")
print("=" * 80)

# Method: for a given evidence e, iterate DS to convergence
def find_fixed_point(e, m_init=None, max_iter=500):
    """Find fixed point by iteration."""
    if m_init is None:
        m = np.array([0.7, 0.05, 0.05, 0.2])
    else:
        m = m_init.copy()
    for i in range(max_iter):
        m_new, K = full_step(m, e)
        if np.max(np.abs(m_new - m)) < 1e-14:
            return m_new, K, True
        m = m_new
    return m, K_conflict(m, e), False

# Sweep: parametrise evidence by dominant fraction p_dom
# e = (p_dom * (1-φ), p_weak * (1-φ), p_weak * (1-φ), φ)
# where p_weak = (1-p_dom)/2

results = []

print(f"\n{'p_dom':>8s} {'phi':>8s} {'K*':>12s} {'K*-7/30':>12s} "
      f"{'delta':>8s} {'lambda0':>10s} {'Delta':>8s} {'theta*':>8s}")
print("-" * 95)

for p_dom in np.linspace(0.40, 0.99, 120):
    for phi in np.linspace(0.02, 0.50, 50):
        p_weak = (1 - p_dom) / 2
        e1 = p_dom * (1 - phi)
        e2 = p_weak * (1 - phi)
        if e1 <= 0 or e2 <= 0 or phi <= 0:
            continue
        e = np.array([e1, e2, e2, phi])
        if abs(np.sum(e) - 1) > 1e-10:
            continue

        m_fp, K_val, converged = find_fixed_point(e)
        if not converged or K_val <= 0 or K_val >= 1:
            continue

        # Check it's actually at Born floor
        b = born(m_fp)
        if abs(b - FLOOR) > 1e-6:
            continue

        # Compute diagonal content delta
        delta = (m_fp[0]*e[0] + m_fp[1]*e[1] + m_fp[2]*e[2]) / \
                (sum(m_fp[i]*e[j] for i in range(3) for j in range(3)))

        # Compute spectral gap
        J = compute_jacobian(m_fp, e)
        eigs = np.linalg.eigvals(J)
        eigs_real = np.sort(np.abs(eigs))[::-1]
        lam0 = eigs_real[0] if len(eigs_real) > 0 else 0

        gap = -np.log(lam0) if lam0 > 0 and lam0 < 1 else 0

        results.append({
            'p_dom': p_dom, 'phi': phi,
            'K': K_val, 'delta': delta,
            'lam0': lam0, 'gap': gap,
            'theta': m_fp[3], 'm': m_fp, 'e': e
        })

# Sort by K* and print a representative subset
results.sort(key=lambda r: r['K'])

# Find the one closest to 7/30
best_idx = min(range(len(results)), key=lambda i: abs(results[i]['K'] - 7/30))

# Print every ~50th result plus the best
printed = set()
step = max(1, len(results) // 40)
for i in list(range(0, len(results), step)) + [best_idx]:
    if i in printed:
        continue
    printed.add(i)
    r = results[i]
    marker = " <--- 7/30" if abs(r['K'] - 7/30) < 0.002 else ""
    print(f"{r['p_dom']:8.4f} {r['phi']:8.4f} {r['K']:12.8f} {r['K']-7/30:12.2e} "
          f"{r['delta']:8.4f} {r['lam0']:10.6f} {r['gap']:8.4f} {r['theta']:8.5f}{marker}")

print(f"\nTotal fixed points found: {len(results)}")

# ============================================================
# Key analysis: what is special about K* = 7/30?
# ============================================================
print("\n" + "=" * 80)
print("ANALYSIS: What is special about K* = 7/30?")
print("=" * 80)

K_vals = [r['K'] for r in results]
gaps = [r['gap'] for r in results]
deltas = [r['delta'] for r in results]
lam0s = [r['lam0'] for r in results]
phis = [r['phi'] for r in results]

print(f"\nK* range: [{min(K_vals):.6f}, {max(K_vals):.6f}]")
print(f"Δ range:  [{min(gaps):.6f}, {max(gaps):.6f}]")
print(f"λ₀ range: [{min(lam0s):.6f}, {max(lam0s):.6f}]")

# Check: at K*=7/30, what's the relationship between phi and K?
# The θ fixed-point equation gives: θ* = floor(θ*φ/(1-K*))
# At the floor: 1-K = φ would hold WITHOUT the floor.
# WITH the floor: the relationship is modified.

r_best = results[best_idx]
print(f"\nAt K* ≈ 7/30:")
print(f"  K*    = {r_best['K']:.10f}")
print(f"  φ*    = {r_best['phi']:.10f}")
print(f"  1-K*  = {1-r_best['K']:.10f}")
print(f"  δ     = {r_best['delta']:.10f}")
print(f"  λ₀    = {r_best['lam0']:.10f}")
print(f"  Δ     = {r_best['gap']:.6f}")
print(f"  θ*    = {r_best['theta']:.10f}")
print(f"  m*    = {r_best['m']}")
print(f"  e*    = {r_best['e']}")

# Check the conservation law at EVERY fixed point
print("\n" + "=" * 80)
print("CONSERVATION LAW TEST: K*(H²+1) - η*H² at each fixed point")
print("=" * 80)

eta = (H-1)**2 / H**3  # = 4/27
print(f"η = {eta:.10f}")
print(f"If conservation law is universal: K*(H²+1) - η*H² = K*·10 - η·9 should = 1")

cons_vals = [r['K'] * (H**2 + 1) - eta * H**2 for r in results]
print(f"\nK*(H²+1) - η·H² across all {len(results)} fixed points:")
print(f"  Range: [{min(cons_vals):.6f}, {max(cons_vals):.6f}]")
print(f"  At K*=7/30: {r_best['K']*(H**2+1) - eta*H**2:.10f}")
print(f"  Mean:  {np.mean(cons_vals):.6f}")
print(f"  Std:   {np.std(cons_vals):.6f}")

# The conservation law value = 1 ONLY at K*=7/30
# What does it equal at other fixed points?
close_to_1 = [r for r in results if abs(r['K']*(H**2+1) - eta*H**2 - 1) < 0.01]
print(f"\n  Fixed points where conservation law ≈ 1: {len(close_to_1)} "
      f"(out of {len(results)})")
if close_to_1:
    K_at_1 = [r['K'] for r in close_to_1]
    print(f"  Their K* values: [{min(K_at_1):.6f}, {max(K_at_1):.6f}]")
    print(f"  All near 7/30 = {7/30:.6f}? {all(abs(k - 7/30) < 0.01 for k in K_at_1)}")

# ============================================================
# THE BIG QUESTION: Is K*=7/30 extremal in some sense?
# ============================================================
print("\n" + "=" * 80)
print("IS K*=7/30 EXTREMAL? (max gap, min λ₀, etc.)")
print("=" * 80)

# Find the fixed point with maximum spectral gap
max_gap_idx = max(range(len(results)), key=lambda i: results[i]['gap'])
r_maxgap = results[max_gap_idx]
print(f"\nMax spectral gap Δ = {r_maxgap['gap']:.6f} at K* = {r_maxgap['K']:.6f}")
print(f"  (7/30 gives Δ = {r_best['gap']:.6f})")

# Find the fixed point with minimum λ₀
min_lam_idx = min(range(len(results)), key=lambda i: results[i]['lam0'])
r_minlam = results[min_lam_idx]
print(f"\nMin λ₀ = {r_minlam['lam0']:.6f} at K* = {r_minlam['K']:.6f}")

# Check: is K*=7/30 where δ = some special value?
print(f"\nδ at K*=7/30: {r_best['delta']:.6f}")
print(f"1/H = {1/H:.6f}")
print(f"δ range: [{min(deltas):.6f}, {max(deltas):.6f}]")

# Check: relationship between K and phi across the family
print("\n" + "=" * 80)
print("K* vs φ relationship across the family")
print("=" * 80)
# Group by similar K values and check what constraint K*=7/30 imposes on φ
for K_target in [0.05, 0.10, 0.15, 7/30, 0.30, 0.35, 0.40]:
    nearby = [r for r in results if abs(r['K'] - K_target) < 0.005]
    if nearby:
        phis_at_K = [r['phi'] for r in nearby]
        gaps_at_K = [r['gap'] for r in nearby]
        print(f"  K*≈{K_target:.4f}: {len(nearby):3d} points, "
              f"φ∈[{min(phis_at_K):.4f},{max(phis_at_K):.4f}], "
              f"Δ∈[{min(gaps_at_K):.4f},{max(gaps_at_K):.4f}]")
