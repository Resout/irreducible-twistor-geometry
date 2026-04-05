"""
Spectral gap at the K*=7/30 self-consistent equilibrium.

The conservation law gives K*=7/30 as the equilibrium conflict.
The evidence distribution at this equilibrium is (eq:fpd in the paper):
  p_dom = 9/10, p_weak = 1/20 each (for the two weak hypotheses).

We need:
1. Construct the evidence mass function from these Born probabilities
2. Find the DS fixed point m* under this evidence with Born floor 1/27
3. Verify K(m*, e) = 7/30 at this fixed point
4. Linearize the transfer operator T: m -> DS(m, e) around m*
5. Compute eigenvalues of the Jacobian -> spectral gap
"""
import numpy as np
from scipy.optimize import fsolve

H = 3
FLOOR = 1.0 / H**3  # = 1/27

# ============================================================
# DS combination rule (real-valued masses)
# ============================================================
def ds_combine(m, e, apply_floor=True):
    """
    Combine mass function m = (s1, s2, s3, theta) with evidence e.
    Returns (m_out, K).
    """
    s = m[:3]
    theta = m[3]
    se = e[:3]
    theta_e = e[3]

    # Pre-normalisation
    s_new = s * se + s * theta_e + theta * se
    theta_new = theta * theta_e

    # Conflict
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)

    # Normalise
    denom = 1.0 - K
    m_out = np.zeros(4)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom

    # Born floor enforcement
    if apply_floor:
        born_theta = m_out[3]**2 / np.sum(m_out**2)
        if born_theta < FLOOR:
            m_out = enforce_floor(m_out)

    return m_out, K

def enforce_floor(m):
    """Enforce Born(theta) >= 1/27 by boosting theta, maintaining L1=1."""
    # Binary search for the right theta
    s = m[:3].copy()
    total = np.sum(m)

    # We want |theta|^2 / sum(|m_i|^2) >= 1/27
    # With L1=1: s1+s2+s3+theta = 1
    # Born(theta) = theta^2 / (s1^2 + s2^2 + s3^2 + theta^2)
    # We boost theta and rescale s_i proportionally to maintain L1=1

    lo, hi = m[3], 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        s_scale = (1.0 - mid) / np.sum(s) if np.sum(s) > 0 else 0
        s_trial = s * s_scale
        born = mid**2 / (np.sum(s_trial**2) + mid**2)
        if born < FLOOR:
            lo = mid
        else:
            hi = mid

    theta_new = (lo + hi) / 2
    s_scale = (1.0 - theta_new) / np.sum(s) if np.sum(s) > 0 else 0
    m_out = np.zeros(4)
    m_out[:3] = s * s_scale
    m_out[3] = theta_new
    return m_out

# ============================================================
# Part 1: Construct evidence from the fixed-point distribution
# ============================================================
# Paper eq:fpd: p_dom = H^2/(H^2+1) = 9/10, p_weak = 1/((H-1)(H^2+1)) = 1/20
# These are Born probabilities for the 3 hypotheses.
# But what about theta? The evidence mass function also has an ignorance component.

# At the equilibrium, the evidence should reproduce the state.
# The Born probabilities from eq:fpd sum to 1 over just the 3 singletons:
# 9/10 + 1/20 + 1/20 = 1. So Born(theta_e) contributes additionally.

# Actually, Born probabilities of ALL 4 components sum to 1.
# If p1 + p2 + p3 + p_theta = 1 and p1=9/10, p2=p3=1/20, then
# p_theta = 1 - 9/10 - 1/20 - 1/20 = 0.
# But floor says p_theta >= 1/27 > 0!

# Resolution: the Born floor IS active. The distribution (9/10, 1/20, 1/20)
# is the pre-floor distribution. After floor enforcement:
# p_theta = 1/27, and the singleton probabilities get rescaled.

# Let's find the mass function with these properties.
# With L1=1: s1 + s2 + s3 + theta = 1
# Born(si) = si^2 / (s1^2 + s2^2 + s3^2 + theta^2)
# Born(theta) = theta^2 / (s1^2 + s2^2 + s3^2 + theta^2) = 1/27

# Let the dominant hypothesis be s1, weak ones s2=s3.
# s1 + 2*s2 + theta = 1
# s1^2 / L2^2 : s2^2 / L2^2 = 9/10 : 1/20 => s1^2 / s2^2 = 18
# => s1 = sqrt(18) * s2 = 3*sqrt(2) * s2

# Born(theta) = 1/27:
# theta^2 / (s1^2 + 2*s2^2 + theta^2) = 1/27
# 27*theta^2 = s1^2 + 2*s2^2 + theta^2
# 26*theta^2 = s1^2 + 2*s2^2 = 18*s2^2 + 2*s2^2 = 20*s2^2
# theta^2 = 20*s2^2/26 = 10*s2^2/13
# theta = s2 * sqrt(10/13)

# L1: s1 + 2*s2 + theta = 1
# 3*sqrt(2)*s2 + 2*s2 + sqrt(10/13)*s2 = 1
# s2 * (3*sqrt(2) + 2 + sqrt(10/13)) = 1

coeff = 3*np.sqrt(2) + 2 + np.sqrt(10/13)
s2_eq = 1.0 / coeff
s1_eq = 3*np.sqrt(2) * s2_eq
theta_eq = np.sqrt(10/13) * s2_eq

print("=" * 60)
print("EVIDENCE MASS FUNCTION AT EQUILIBRIUM")
print("=" * 60)
e_eq = np.array([s1_eq, s2_eq, s2_eq, theta_eq])
print(f"  e = ({s1_eq:.6f}, {s2_eq:.6f}, {s2_eq:.6f}, {theta_eq:.6f})")
print(f"  L1 = {np.sum(e_eq):.10f}")
L2sq = np.sum(e_eq**2)
born = e_eq**2 / L2sq
print(f"  Born probs: ({born[0]:.6f}, {born[1]:.6f}, {born[2]:.6f}, {born[3]:.6f})")
print(f"  Born(theta) = {born[3]:.6f} (floor = {FLOOR:.6f})")
print(f"  Singleton ratio: p1/p2 = {born[0]/born[1]:.2f} (expect 18)")

# ============================================================
# Part 2: Find the DS fixed point under this evidence
# ============================================================
print("\n" + "=" * 60)
print("FINDING DS FIXED POINT")
print("=" * 60)

# Iterate DS combination with fixed evidence e_eq
m = np.array([0.4, 0.2, 0.2, 0.2])  # initial guess
for i in range(500):
    m_new, K = ds_combine(m, e_eq)
    diff = np.max(np.abs(m_new - m))
    if i < 5 or i % 50 == 0 or diff < 1e-14:
        print(f"  Step {i}: m = ({m_new[0]:.8f}, {m_new[1]:.8f}, {m_new[2]:.8f}, {m_new[3]:.8f}), K = {K:.8f}, diff = {diff:.2e}")
    if diff < 1e-15:
        print(f"  Converged at step {i}")
        break
    m = m_new

m_star = m_new
_, K_star = ds_combine(m_star, e_eq, apply_floor=False)
print(f"\n  Fixed point: m* = ({m_star[0]:.10f}, {m_star[1]:.10f}, {m_star[2]:.10f}, {m_star[3]:.10f})")
print(f"  K* = {K_star:.10f}")
print(f"  7/30 = {7/30:.10f}")
print(f"  Difference: {abs(K_star - 7/30):.2e}")

born_star = m_star**2 / np.sum(m_star**2)
print(f"  Born(theta) = {born_star[3]:.6f}")

# Also try self-evidence for comparison
print("\n--- Self-evidence comparison ---")
m_self = np.array([0.4, 0.2, 0.2, 0.2])
for i in range(500):
    m_self_new, K_self = ds_combine(m_self, m_self)
    diff = np.max(np.abs(m_self_new - m_self))
    if diff < 1e-15:
        break
    m_self = m_self_new
_, K_self_fp = ds_combine(m_self, m_self, apply_floor=False)
print(f"  Self-evidence FP: K* = {K_self_fp:.6f}")

# ============================================================
# Part 3: Compute the Jacobian of the transfer operator
# ============================================================
print("\n" + "=" * 60)
print("JACOBIAN OF TRANSFER OPERATOR")
print("=" * 60)

def transfer_op(m, e):
    """Full transfer operator: m -> DS(m, e) with floor."""
    return ds_combine(m, e, apply_floor=True)[0]

# Numerical Jacobian via finite differences
eps = 1e-8
J = np.zeros((4, 4))
f0 = transfer_op(m_star, e_eq)
for j in range(4):
    m_pert = m_star.copy()
    m_pert[j] += eps
    # Maintain L1=1 by adjusting another component
    # Actually, the constraint is L1=1. Perturbations must live in the
    # tangent space of the L1=1 constraint surface.
    # The tangent space has directions where sum(dm) = 0.
    # So we perturb m_j by +eps and m[3] (theta) by -eps to stay on L1=1.
    # But this couples the perturbation. Better: use unconstrained Jacobian
    # and then project.

    # Actually, let's compute the full 4x4 Jacobian first, then analyse
    # eigenvalues on the constraint surface.
    f_pert = transfer_op(m_pert, e_eq)
    J[:, j] = (f_pert - f0) / eps

print("  Full 4x4 Jacobian:")
for i in range(4):
    print(f"    [{J[i,0]:+.6f}  {J[i,1]:+.6f}  {J[i,2]:+.6f}  {J[i,3]:+.6f}]")

evals_full = np.linalg.eigvals(J)
print(f"\n  Eigenvalues (full): {np.sort(np.abs(evals_full))[::-1]}")

# The L1=1 constraint means we should look at the Jacobian restricted
# to the tangent space of {sum(m)=1}. This is 3-dimensional.
# Project: basis vectors e1-e4, e2-e4, e3-e4 (all have sum=0)
print("\n  Projected to tangent space of L1=1:")

# Tangent vectors: v1 = e1-e4, v2 = e2-e4, v3 = e3-e4
V = np.array([[1,0,0,-1], [0,1,0,-1], [0,0,1,-1]], dtype=float).T  # 4x3

# Project Jacobian: J_proj = V^T @ J @ V  (but need pseudoinverse for projection)
# Actually, J maps R^4 -> R^4. We want the restriction to the hyperplane sum=0.
# The constraint surface has tangent vectors in V's column space.
# J_proj = (V^T V)^{-1} V^T J V
VTV_inv = np.linalg.inv(V.T @ V)
J_proj = VTV_inv @ V.T @ J @ V

print("  3x3 projected Jacobian:")
for i in range(3):
    print(f"    [{J_proj[i,0]:+.8f}  {J_proj[i,1]:+.8f}  {J_proj[i,2]:+.8f}]")

evals_proj = np.linalg.eigvals(J_proj)
evals_abs = np.sort(np.abs(evals_proj))[::-1]
print(f"\n  Eigenvalues (projected): {evals_abs}")
print(f"  |λ₀| = {evals_abs[0]:.8f}")
print(f"  |λ₁| = {evals_abs[1]:.8f}")
print(f"  |λ₂| = {evals_abs[2]:.8f}")

# ============================================================
# Part 4: Spectral gap
# ============================================================
print("\n" + "=" * 60)
print("SPECTRAL GAP")
print("=" * 60)

lambda_0 = evals_abs[0]
if lambda_0 < 1:
    Delta = -np.log(lambda_0)
    print(f"  λ₀ = {lambda_0:.8f} < 1  =>  GAP EXISTS")
    print(f"  Δ = -ln(λ₀) = {Delta:.8f}")
    print(f"  Structural filter rate = -ln(1-7/30) = {-np.log(23/30):.8f}")
    print(f"  Ratio Δ/(-ln(1-K*)) = {Delta / (-np.log(23/30)):.4f}")
else:
    print(f"  λ₀ = {lambda_0:.8f} >= 1  =>  NO GAP (marginal or unstable)")

# ============================================================
# Part 5: Sensitivity analysis — vary evidence distribution
# ============================================================
print("\n" + "=" * 60)
print("SENSITIVITY: EVIDENCE DISTRIBUTION SCAN")
print("=" * 60)

# What if we use self-evidence at the fixed point?
print("\n--- Fixed point under self-evidence ---")
m_se = m_self.copy()
eps2 = 1e-8
J_se = np.zeros((4, 4))
f0_se = transfer_op(m_se, m_se)
# For self-evidence, the Jacobian includes both m and e=m changing together.
# But the PHYSICAL spectral gap uses FIXED evidence (perturbation of state only).
# So we linearize T(m) = DS(m, e*) with e* = m* fixed.
for j in range(4):
    m_pert = m_se.copy()
    m_pert[j] += eps2
    f_pert = transfer_op(m_pert, m_se)  # evidence stays at m_se
    J_se[:, j] = (f_pert - f0_se) / eps2

J_se_proj = VTV_inv @ V.T @ J_se @ V
evals_se = np.sort(np.abs(np.linalg.eigvals(J_se_proj)))[::-1]
print(f"  Self-evidence eigenvalues: {evals_se}")
if evals_se[0] < 1:
    print(f"  λ₀ = {evals_se[0]:.6f}, Δ = {-np.log(evals_se[0]):.6f}")

# ============================================================
# Part 6: Verify K at various evidence strengths
# ============================================================
print("\n" + "=" * 60)
print("EVIDENCE STRENGTH SCAN")
print("=" * 60)
print(f"  {'p_dom':>8}  {'K*':>10}  {'λ₀':>10}  {'Δ':>10}  {'floor?':>6}")

for p_dom in [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]:
    p_weak = (1 - p_dom) / 2
    # Construct evidence with these Born probs and Born(theta)=1/27
    # si^2/L2^2 = pi (for i=1,2,3), theta^2/L2^2 = 1/27
    # With L1=1, this is a constrained system.
    # s1^2 : s2^2 : s3^2 : theta^2 = p_dom : p_weak : p_weak : 1/27
    # But sum must be 1/(1+1/27*L2^2/theta^2)... this is circular.
    #
    # Simpler: set Born probs for singletons only (3 components),
    # then add theta to satisfy floor.
    # p_dom + 2*p_weak + p_theta = 1 where p_theta = 1/27.
    # But p_dom + 2*p_weak already = 1 in paper's formula.
    # So with floor: scale down singleton probs.

    p_theta = FLOOR
    scale = 1 - p_theta
    p1 = p_dom * scale
    p2 = p_weak * scale

    # m_i = sqrt(p_i) * L2, with L1 = sum(m_i) = 1
    # m_i proportional to sqrt(p_i), then normalise to L1=1
    raw = np.array([np.sqrt(p1), np.sqrt(p2), np.sqrt(p2), np.sqrt(p_theta)])
    e_test = raw / np.sum(raw)

    # Find fixed point
    m_test = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(500):
        m_new, _ = ds_combine(m_test, e_test)
        if np.max(np.abs(m_new - m_test)) < 1e-15:
            break
        m_test = m_new

    _, K_test = ds_combine(m_test, e_test, apply_floor=False)

    # Jacobian
    J_test = np.zeros((4, 4))
    f0_test = transfer_op(m_test, e_test)
    for j in range(4):
        mp = m_test.copy()
        mp[j] += eps
        fp = transfer_op(mp, e_test)
        J_test[:, j] = (fp - f0_test) / eps

    J_test_proj = VTV_inv @ V.T @ J_test @ V
    evals_test = np.sort(np.abs(np.linalg.eigvals(J_test_proj)))[::-1]
    lam0 = evals_test[0]
    delta = -np.log(lam0) if lam0 < 1 else float('nan')

    born_test = m_test**2 / np.sum(m_test**2)
    floor_active = born_test[3] <= FLOOR + 1e-6

    print(f"  {p_dom:8.3f}  {K_test:10.6f}  {lam0:10.6f}  {delta:10.6f}  {'yes' if floor_active else 'no':>6}")

# ============================================================
# Part 7: The physical spectral gap — at the EXACT K*=7/30 point
# ============================================================
print("\n" + "=" * 60)
print("FINDING EVIDENCE THAT GIVES EXACTLY K*=7/30")
print("=" * 60)

# We need to find the evidence e such that the fixed point m* of
# DS(m, e) gives K(m*, e) = 7/30 exactly.
# From the scan above, we can see what p_dom gives K*=7/30.

from scipy.optimize import brentq

def K_at_equilibrium(p_dom_val):
    """Given dominant Born probability, find K* at the DS fixed point."""
    p_weak_val = (1 - p_dom_val) / 2
    p_theta = FLOOR
    scale = 1 - p_theta
    p1 = p_dom_val * scale
    p2 = p_weak_val * scale

    raw = np.array([np.sqrt(p1), np.sqrt(p2), np.sqrt(p2), np.sqrt(p_theta)])
    if np.sum(raw) == 0:
        return 0
    e = raw / np.sum(raw)

    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(1000):
        m_new, _ = ds_combine(m, e)
        if np.max(np.abs(m_new - m)) < 1e-15:
            break
        m = m_new

    _, K = ds_combine(m, e, apply_floor=False)
    return K

# Search for p_dom that gives K*=7/30
target_K = 7.0 / 30.0
print(f"  Target K* = {target_K:.10f}")

# First, scan to find the bracket
p_vals = np.linspace(0.4, 0.99, 50)
K_vals = [K_at_equilibrium(p) for p in p_vals]

for i in range(len(p_vals)):
    if i > 0 and (K_vals[i-1] - target_K) * (K_vals[i] - target_K) < 0:
        print(f"  Bracket found: p_dom in [{p_vals[i-1]:.4f}, {p_vals[i]:.4f}]")
        print(f"  K values: [{K_vals[i-1]:.6f}, {K_vals[i]:.6f}]")

        # Refine with Brent's method
        p_exact = brentq(lambda p: K_at_equilibrium(p) - target_K,
                         p_vals[i-1], p_vals[i], xtol=1e-12)
        K_check = K_at_equilibrium(p_exact)
        print(f"  Exact p_dom = {p_exact:.10f}")
        print(f"  K* at this p_dom = {K_check:.10f}")
        print(f"  |K* - 7/30| = {abs(K_check - target_K):.2e}")

        # Now compute the spectral gap at this exact point
        p_w = (1 - p_exact) / 2
        p_theta = FLOOR
        sc = 1 - p_theta
        raw = np.array([np.sqrt(p_exact*sc), np.sqrt(p_w*sc), np.sqrt(p_w*sc), np.sqrt(p_theta)])
        e_exact = raw / np.sum(raw)

        m = np.array([0.4, 0.2, 0.2, 0.2])
        for _ in range(1000):
            m_new, _ = ds_combine(m, e_exact)
            if np.max(np.abs(m_new - m)) < 1e-15:
                break
            m = m_new
        m_exact = m_new

        print(f"\n  Fixed point: m* = ({m_exact[0]:.10f}, {m_exact[1]:.10f}, {m_exact[2]:.10f}, {m_exact[3]:.10f})")
        print(f"  Evidence:    e  = ({e_exact[0]:.10f}, {e_exact[1]:.10f}, {e_exact[2]:.10f}, {e_exact[3]:.10f})")

        born_exact = m_exact**2 / np.sum(m_exact**2)
        print(f"  Born(m*): ({born_exact[0]:.6f}, {born_exact[1]:.6f}, {born_exact[2]:.6f}, {born_exact[3]:.6f})")

        # Jacobian at exact point
        J_exact = np.zeros((4, 4))
        f0_exact = transfer_op(m_exact, e_exact)
        for j in range(4):
            mp = m_exact.copy()
            mp[j] += eps
            fp = transfer_op(mp, e_exact)
            J_exact[:, j] = (fp - f0_exact) / eps

        J_exact_proj = VTV_inv @ V.T @ J_exact @ V
        evals_exact = np.linalg.eigvals(J_exact_proj)
        evals_exact_abs = np.sort(np.abs(evals_exact))[::-1]

        print(f"\n  Eigenvalues: {evals_exact}")
        print(f"  |eigenvalues|: {evals_exact_abs}")

        lambda_0_exact = evals_exact_abs[0]
        if lambda_0_exact < 1:
            Delta_exact = -np.log(lambda_0_exact)
            print(f"\n  *** SPECTRAL GAP AT K*=7/30 ***")
            print(f"  λ₀ = {lambda_0_exact:.10f}")
            print(f"  Δ = -ln(λ₀) = {Delta_exact:.10f}")
            print(f"  For comparison:")
            print(f"    Structural filter rate = -ln(23/30) = {-np.log(23/30):.10f}")
            print(f"    Self-evidence gap (from memory) = 0.196")
            print(f"    Ratio Δ/filter_rate = {Delta_exact / (-np.log(23/30)):.6f}")
        else:
            print(f"\n  λ₀ = {lambda_0_exact:.10f} >= 1 — NO GAP")

        break
else:
    print("  No bracket found! Dumping K* values:")
    for p, K in zip(p_vals, K_vals):
        print(f"    p_dom={p:.3f}  K*={K:.6f}")
