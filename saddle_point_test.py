"""
Test: Is the DS equilibrium the saddle point of the YM action?

If K is the integration constant (action density at the critical point),
then:
  1. The conservation law should be δS=0 (stationarity)
  2. The spectral gap should come from δ²S (Hessian)
  3. The path integral Z = ∫DA e^{-S} ≈ e^{-S*} / √det(δ²S)
     (Gaussian/saddle-point approximation)

Concrete test: does the Hessian of K (as a function of the mass
function m) at the equilibrium give eigenvalues related to λ₀=0.283?

If K plays the role of the action, then:
  - The DS transfer operator T(m) = DS(m,e)/(1-K) is the RENORMALIZED dynamics
  - The Jacobian of T at the fixed point gives λ₀
  - But the Hessian of K at the fixed point should ALSO give information
    about the spectral gap

Let's compute both and see if they're related.
"""
import numpy as np

H = 3
FLOOR = 1.0 / H**3

def enforce_floor(m):
    s = m[:3].copy()
    lo, hi = m[3], 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        s_scale = (1.0 - mid) / np.sum(s) if np.sum(s) > 0 else 0
        s_trial = s * s_scale
        born = mid**2 / (np.sum(s_trial**2) + mid**2)
        if born < FLOOR: lo = mid
        else: hi = mid
    theta_new = (lo + hi) / 2
    s_scale = (1.0 - theta_new) / np.sum(s) if np.sum(s) > 0 else 0
    m_out = np.zeros(4)
    m_out[:3] = s * s_scale
    m_out[3] = theta_new
    return m_out

def ds_combine(m, e, apply_floor=True):
    s, theta = m[:3], m[3]
    se, theta_e = e[:3], e[3]
    s_new = s * se + s * theta_e + theta * se
    theta_new = theta * theta_e
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    denom = 1.0 - K
    m_out = np.zeros(4)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom
    if apply_floor:
        born_theta = m_out[3]**2 / np.sum(m_out**2)
        if born_theta < FLOOR:
            m_out = enforce_floor(m_out)
    return m_out, K

# K*=7/30 equilibrium
p_dom_exact = 0.9322756157
p_weak_exact = (1 - p_dom_exact) / 2
raw = np.array([np.sqrt(p_dom_exact*(1-FLOOR)), np.sqrt(p_weak_exact*(1-FLOOR)),
                np.sqrt(p_weak_exact*(1-FLOOR)), np.sqrt(FLOOR)])
e_exact = raw / np.sum(raw)

m = np.array([0.4, 0.2, 0.2, 0.2])
for _ in range(1000):
    m_new, _ = ds_combine(m, e_exact)
    if np.max(np.abs(m_new - m)) < 1e-15:
        break
    m = m_new
m_star = m_new

_, K_star = ds_combine(m_star, e_exact, apply_floor=False)
print(f"Fixed point: m* = {m_star}")
print(f"K* = {K_star:.10f}")

# ============================================================
# Compute the HESSIAN of K(m) at the fixed point
# K(m) = Σ_{i≠j} s_i * e_j(m) where e depends on m through
# the evidence structure.
#
# For FIXED evidence: K(m) = Σ_{i≠j} m_i * e_j is LINEAR in m.
# So the Hessian of K w.r.t. m is ZERO.
#
# The interesting object is the Hessian of K along the
# co-evolution trajectory: K(m, e(m)).
# ============================================================

print("\n" + "=" * 60)
print("TEST 1: K as function of m (fixed evidence)")
print("=" * 60)

def K_of_m(m, e):
    """K as a function of m with fixed evidence e."""
    s = m[:3]
    return sum(s[i] * e[j] for i in range(3) for j in range(3) if i != j)

# Gradient of K w.r.t. m (fixed evidence)
eps = 1e-8
grad_K = np.zeros(4)
for j in range(4):
    mp = m_star.copy(); mp[j] += eps
    mm = m_star.copy(); mm[j] -= eps
    grad_K[j] = (K_of_m(mp, e_exact[:3]) - K_of_m(mm, e_exact[:3])) / (2*eps)

print(f"∇K at m* = {grad_K}")
print(f"|∇K| = {np.linalg.norm(grad_K):.6f}")

# K is linear in m (with fixed e), so the gradient is constant
# and the Hessian is zero. Not useful.

# ============================================================
# TEST 2: The "action" is not K itself but the LOGARITHMIC action
# In the path integral: Z = ∫ e^{-S} ≈ e^{-nS*} for n steps
# The DS dynamics give: (1-K)^n = e^{-nΔ_filter}
# So the "action per step" is -ln(1-K) = Δ_filter
# ============================================================

print("\n" + "=" * 60)
print("TEST 2: -ln(1-K) as the action per step")
print("=" * 60)

# The transfer operator eigenvalue λ₀ gives:
# perturbation after n steps ~ λ₀^n = e^{-nΔ}
# The total mass after n steps: (1-K)^n = e^{-nΔ_filter}
#
# If K is the "constant dropped" and -ln(1-K) is the "action per step",
# then the spectral gap Δ should be related to the SECOND DERIVATIVE
# of -ln(1-K) at K=K*.

# -ln(1-K) at K=K*=7/30:
S_star = -np.log(1 - K_star)
print(f"S* = -ln(1-K*) = {S_star:.6f}")
print(f"Δ (spectral gap) = {-np.log(0.2829):.6f}")
print(f"Ratio Δ/S* = {-np.log(0.2829)/S_star:.4f}")

# d/dK[-ln(1-K)] = 1/(1-K)
# d²/dK²[-ln(1-K)] = 1/(1-K)²
# At K=7/30: 1/(1-7/30)² = (30/23)² = 900/529 = 1.7013

hessian_S = 1/(1-K_star)**2
print(f"\nHessian of S=-ln(1-K) at K*:")
print(f"  d²S/dK² = 1/(1-K*)² = {hessian_S:.6f}")

# The spectral gap from the Hessian would be:
# Δ_hessian = sqrt(d²S/dK²) / (normalization factor)
# Let's see if there's a simple relationship

print(f"\n  sqrt(d²S/dK²) = {np.sqrt(hessian_S):.6f}")
print(f"  Δ = {-np.log(0.2829):.6f}")
print(f"  Ratio = {-np.log(0.2829)/np.sqrt(hessian_S):.6f}")

# ============================================================
# TEST 3: The FULL test — saddle point of the "DS action"
# ============================================================
print("\n" + "=" * 60)
print("TEST 3: Saddle point structure")
print("=" * 60)

# Define the "DS action" as S(m) = -ln(1-K(m,e(m)))
# This is the per-step action along the co-evolution.
#
# At the fixed point: S(m*) = -ln(1-K*) = -ln(23/30) = 0.266
#
# The stationarity condition: dS/dm = 0 at m*
# This should be equivalent to the fixed-point condition DS(m*,e*)=m*.
#
# The Hessian: d²S/dm² at m* should give the spectral gap.

# Let's define S(m) = -ln(1-K(m)) where K is computed from DS dynamics
# at the co-evolving equilibrium.

def S_action(m_test):
    """Action = -ln(1-K) at mass function m with co-evolving evidence."""
    # Use the K*=7/30 evidence (fixed, not co-evolving, for simplicity)
    _, K = ds_combine(m_test, e_exact, apply_floor=False)
    if K >= 1:
        return 100  # regularize
    return -np.log(1 - K)

# Gradient of S at m*
grad_S = np.zeros(4)
for j in range(4):
    mp = m_star.copy(); mp[j] += eps
    mm = m_star.copy(); mm[j] -= eps
    grad_S[j] = (S_action(mp) - S_action(mm)) / (2*eps)

print(f"∇S at m* = {grad_S}")
print(f"|∇S| = {np.linalg.norm(grad_S):.6f}")
print(f"∇S is {'approximately zero' if np.linalg.norm(grad_S) < 0.1 else 'NOT zero'}")

# Hessian of S at m*
H_S = np.zeros((4,4))
S0 = S_action(m_star)
for i in range(4):
    for j in range(4):
        mpp = m_star.copy(); mpp[i] += eps; mpp[j] += eps
        mpm = m_star.copy(); mpm[i] += eps; mpm[j] -= eps
        mmp = m_star.copy(); mmp[i] -= eps; mmp[j] += eps
        mmm = m_star.copy(); mmm[i] -= eps; mmm[j] -= eps
        H_S[i,j] = (S_action(mpp) - S_action(mpm) - S_action(mmp) + S_action(mmm)) / (4*eps**2)

print(f"\nHessian of S:")
for i in range(4):
    print(f"  [{H_S[i,0]:+.4f}  {H_S[i,1]:+.4f}  {H_S[i,2]:+.4f}  {H_S[i,3]:+.4f}]")

# Project to L1=1 tangent space
V = np.array([[1,0,0,-1], [0,1,0,-1], [0,0,1,-1]], dtype=float).T
VTV_inv = np.linalg.inv(V.T @ V)
H_proj = VTV_inv @ V.T @ H_S @ V

print(f"\nProjected Hessian (3x3):")
for i in range(3):
    print(f"  [{H_proj[i,0]:+.6f}  {H_proj[i,1]:+.6f}  {H_proj[i,2]:+.6f}]")

evals_H = np.linalg.eigvals(H_proj)
evals_H_sorted = np.sort(np.real(evals_H))[::-1]
print(f"\nHessian eigenvalues: {evals_H_sorted}")

# ============================================================
# TEST 4: Compare transfer operator eigenvalue to Hessian
# ============================================================
print("\n" + "=" * 60)
print("TEST 4: Hessian vs transfer operator")
print("=" * 60)

# Transfer operator eigenvalue (from earlier): λ₀ = 0.2829
# Spectral gap: Δ = -ln(λ₀) = 1.263

# If the path integral is Z ≈ e^{-nS*} · (det H)^{-1/2}
# then the mass gap comes from the Hessian eigenvalues.
# In the Gaussian approximation:
# Δ_gaussian = min(eigenvalues of H) / 2  (or some function)

lambda_0 = 0.2829
Delta = -np.log(lambda_0)

print(f"Transfer operator: λ₀ = {lambda_0:.4f}, Δ = {Delta:.4f}")
print(f"Hessian eigenvalues: {evals_H_sorted}")
print()

# Check various relationships
for i, ev in enumerate(evals_H_sorted):
    if ev > 0:
        print(f"  H_eig[{i}] = {ev:.6f}")
        print(f"    sqrt(H_eig) = {np.sqrt(ev):.6f}")
        print(f"    H_eig/Δ = {ev/Delta:.6f}")
        print(f"    exp(-H_eig) = {np.exp(-ev):.6f}")
        print(f"    λ₀ = {lambda_0:.6f}")
        print(f"    exp(-sqrt(H_eig)) = {np.exp(-np.sqrt(ev)):.6f}")
        print()

# ============================================================
# TEST 5: The key relationship
# ============================================================
print("=" * 60)
print("TEST 5: Is there a clean relationship?")
print("=" * 60)

# In the saddle-point approximation of a path integral:
# Z = ∫ e^{-S[φ]} Dφ ≈ e^{-S[φ*]} · (2π)^{n/2} / sqrt(det(S''[φ*]))
#
# The correlation function at distance d:
# <φ(0)φ(d)> ~ e^{-m·d} where m = mass gap
#
# For a quadratic action S = (1/2)φ·H·φ:
# <φ(0)φ(d)> = (H^{-1})_d ~ e^{-√(min eigenvalue of H) · d}
# if H is a lattice Laplacian.
#
# But our H is a 3x3 matrix (finite dimensional), not a lattice operator.
# The eigenvalues of our projected Hessian describe the curvature of
# the action in the 3 independent directions on the L1=1 surface.

# The transfer operator approach:
# <O(0)O(d)> = Σ_i c_i λ_i^d
# mass gap = -ln(λ_0) where λ_0 is the largest eigenvalue < 1

# For the relationship to work, we need:
# λ_0 = e^{-something from Hessian}
# or equivalently: Δ = something from Hessian

# Let's check: is Δ = S* * (something involving Hessian eigenvalues)?
print(f"  S* = {S_star:.6f}")
print(f"  Δ = {Delta:.6f}")
print(f"  Δ/S* = {Delta/S_star:.6f}")
print(f"  This ratio is {Delta/S_star:.4f}")
print()

# Δ/S* ≈ 4.75 (we knew this — it's the ratio of spectral gap to filter rate)
# Is 4.75 related to the Hessian eigenvalues?
if evals_H_sorted[0] > 0:
    print(f"  Hessian max eigenvalue = {evals_H_sorted[0]:.6f}")
    print(f"  S* × Hessian_max = {S_star * evals_H_sorted[0]:.6f}")
    print(f"  Δ = {Delta:.6f}")
    print(f"  Match? {abs(S_star * evals_H_sorted[0] - Delta) < 0.1}")
