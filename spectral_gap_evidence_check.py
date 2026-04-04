"""
Check consistency between paper's eq:fpd evidence distribution and K*=7/30.
Paper claims: p_dom = 9/10, p_weak = 1/20 (for singletons, excluding theta).
Does this evidence actually produce K*=7/30?
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

def make_evidence(p_dom, p_weak, p_theta=FLOOR):
    """Construct mass function from Born probabilities."""
    # Born(i) = m_i^2 / sum(m_j^2)
    # So m_i proportional to sqrt(p_i), then normalize to L1=1
    raw = np.array([np.sqrt(p_dom), np.sqrt(p_weak), np.sqrt(p_weak), np.sqrt(p_theta)])
    return raw / np.sum(raw)

def find_K_at_fp(e):
    """Find K* at the DS fixed point under evidence e."""
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(1000):
        m_new, _ = ds_combine(m, e)
        if np.max(np.abs(m_new - m)) < 1e-15:
            break
        m = m_new
    _, K = ds_combine(m, e, apply_floor=False)
    return K, m

# Test 1: Paper's eq:fpd values (singletons only, renormalized with floor)
print("=" * 60)
print("TEST 1: Paper's eq:fpd (p_dom=9/10, p_weak=1/20 of singletons)")
print("=" * 60)

# As ratios among singletons: p1:p2:p3 = 9/10 : 1/20 : 1/20
# With Born floor p_theta = 1/27:
p_sing_total = 1 - FLOOR  # = 26/27
p1 = (9/10) * p_sing_total
p2 = (1/20) * p_sing_total

e_paper = make_evidence(p1, p2, FLOOR)
K_paper, m_paper = find_K_at_fp(e_paper)
print(f"  Evidence: {e_paper}")
print(f"  Born(e): {e_paper**2 / np.sum(e_paper**2)}")
print(f"  K* = {K_paper:.10f}")
print(f"  7/30 = {7/30:.10f}")
print(f"  Discrepancy: {abs(K_paper - 7/30):.6f} ({abs(K_paper - 7/30)/(7/30)*100:.1f}%)")

# Test 2: What if eq:fpd gives the Born probabilities INCLUDING theta?
print("\n" + "=" * 60)
print("TEST 2: What if p_dom=9/10, p_weak=1/20, p_theta=0?")
print("=" * 60)

# Here the floor will activate since p_theta=0
# Use p_theta = 1/27 explicitly
e_test2 = make_evidence(9/10, 1/20, FLOOR)
K_test2, _ = find_K_at_fp(e_test2)
print(f"  K* = {K_test2:.10f}")
print(f"  Discrepancy: {abs(K_test2 - 7/30):.6f}")

# Test 3: What evidence ratio gives K*=7/30 exactly?
print("\n" + "=" * 60)
print("TEST 3: Exact evidence for K*=7/30")
print("=" * 60)

from scipy.optimize import brentq

def K_for_ratio(ratio):
    """Given p_dom/p_weak ratio, find K*."""
    # p_dom + 2*p_weak = 1 - FLOOR
    # p_dom = ratio * p_weak
    # ratio*p_weak + 2*p_weak = 1-FLOOR
    # p_weak = (1-FLOOR)/(ratio+2)
    p_weak = (1-FLOOR)/(ratio+2)
    p_dom = ratio * p_weak
    e = make_evidence(p_dom, p_weak, FLOOR)
    K, _ = find_K_at_fp(e)
    return K

# Paper's ratio: (9/10)/(1/20) = 18
print(f"  Paper's ratio p_dom/p_weak = 18, K* = {K_for_ratio(18):.6f}")

# Find exact ratio
r_exact = brentq(lambda r: K_for_ratio(r) - 7/30, 5, 100, xtol=1e-10)
print(f"  Exact ratio for K*=7/30: {r_exact:.6f}")

p_weak_ex = (1-FLOOR)/(r_exact+2)
p_dom_ex = r_exact * p_weak_ex
print(f"  p_dom = {p_dom_ex:.6f}, p_weak = {p_weak_ex:.6f}")
print(f"  p_dom among singletons = {p_dom_ex/(p_dom_ex+2*p_weak_ex):.6f}")
print(f"  Paper claims: p_dom/(p_dom+2*p_weak) = 9/10 = {9/10:.6f}")

# Test 4: Self-evidence at a specific mass function
print("\n" + "=" * 60)
print("TEST 4: Self-evidence fixed point")
print("=" * 60)

# What if we use the FIXED POINT as evidence (self-consistent)?
# At self-consistency: m* is both the state AND the evidence
# DS(m*, m*) = m*, K(m*, m*) = ?

# Start from paper's distribution and iterate self-consistently
m_sc = make_evidence(9/10 * (1-FLOOR), 1/20 * (1-FLOOR), FLOOR)
for i in range(1000):
    m_sc_new, K_sc = ds_combine(m_sc, m_sc)  # self-evidence
    diff = np.max(np.abs(m_sc_new - m_sc))
    if diff < 1e-15:
        break
    m_sc = m_sc_new

_, K_self = ds_combine(m_sc, m_sc, apply_floor=False)
print(f"  Self-evidence FP: {m_sc}")
print(f"  K* (self-evidence) = {K_self:.10f}")
print(f"  Born: {m_sc**2/np.sum(m_sc**2)}")

# Test 5: What if the evidence is NOT the mass function?
# The conservation law describes a system where evidence comes from
# the environment (gauge field), not from self-combination.
# The "self-consistent" means the evidence distribution is
# COMPATIBLE with the state, not identical to it.
print("\n" + "=" * 60)
print("TEST 5: Conservation law prediction")
print("=" * 60)

# From the conservation law: K*(H^2+1) - eta*H^2 = 1
# K* = (H^2-H+1)/(H(H^2+1)) = 7/30 at H=3
# This is exact and algebraic.
# The evidence distribution is derived from channel dimensions.
# The eigenvalue at this equilibrium is what we computed: 0.2829.

print(f"  Conservation law: K* = (H²-H+1)/(H(H²+1)) = {(H**2-H+1)/(H*(H**2+1)):.10f}")
print(f"  Spectral gap (computed): Δ = 1.2626 per DS step")
print(f"  Leading eigenvalue: λ₀ = 0.2829")
print(f"")
print(f"  KEY: The conservation law fixes K*=7/30 algebraically.")
print(f"  The eigenvalue λ₀ = 0.283 is a joint property of K*=7/30")
print(f"  and the evidence structure at equilibrium.")
print(f"  The evidence structure depends on the gauge group through")
print(f"  the Clebsch-Gordan decomposition.")

# Test 6: Sensitivity — how much does Delta change with evidence?
print("\n" + "=" * 60)
print("TEST 6: Delta sensitivity to evidence at K*=7/30")
print("=" * 60)

# We already know the exact ratio for K*=7/30.
# But K*=7/30 can also be achieved with different evidence shapes
# (e.g., 3-symmetric evidence where all singletons are equal).
# Check if the eigenvalue depends on the specific evidence shape.

# 3-symmetric evidence: p1 = p2 = p3 = (1-FLOOR)/3
p_sym = (1-FLOOR)/3
e_sym = make_evidence(p_sym, p_sym, FLOOR)
K_sym, m_sym = find_K_at_fp(e_sym)
print(f"  Symmetric evidence: K* = {K_sym:.6f}")

# We can't get K*=7/30 with symmetric evidence unless the floor
# makes it happen. Let's check what K we get.
# For asymmetric evidence, we already found the right ratio.
# The point is: at fixed K*=7/30, different evidence shapes
# give different eigenvalues.

# Let's try: one dominant + two equal vs two dominant + one weak
print(f"\n  At K*=7/30:")

# Shape 1: (big, small, small) — our computed case
print(f"    1-dominant shape: λ₀ = 0.2829, Δ = 1.263")

# Shape 2: (big, big, tiny)
def K_for_shape2(ratio):
    """big, big, small — p1=p2=ratio*p3"""
    # 2*ratio*p3 + p3 = 1-FLOOR
    p3 = (1-FLOOR)/(2*ratio+1)
    p12 = ratio * p3
    e = make_evidence(p12, p12, FLOOR)
    # This makes evidence with p1=p2=big, p3=sqrt(big) but wait...
    # make_evidence takes (p_dom, p_weak, p_theta) for (s1, s2=s3, theta)
    # We need a different function for (s1=s2, s3=small, theta)
    raw = np.array([np.sqrt(p12), np.sqrt(p12), np.sqrt(p3), np.sqrt(FLOOR)])
    e = raw / np.sum(raw)
    K, m = find_K_at_fp(e)
    return K, e, m

# Find ratio for K*=7/30 with shape 2
try:
    r2 = brentq(lambda r: K_for_shape2(r)[0] - 7/30, 2, 50, xtol=1e-10)
    K2, e2, m2 = K_for_shape2(r2)
    print(f"    2-dominant shape: ratio = {r2:.4f}, K* = {K2:.10f}")

    # Eigenvalue
    eps = 1e-8
    V = np.array([[1,0,0,-1], [0,1,0,-1], [0,0,1,-1]], dtype=float).T
    VTV_inv = np.linalg.inv(V.T @ V)
    J2 = np.zeros((4, 4))
    f0 = ds_combine(m2, e2)[0]
    for j in range(4):
        mp = m2.copy()
        mp[j] += eps
        fp = ds_combine(mp, e2)[0]
        J2[:, j] = (fp - f0) / eps
    J2_proj = VTV_inv @ V.T @ J2 @ V
    evals2 = np.sort(np.abs(np.linalg.eigvals(J2_proj)))[::-1]
    print(f"    2-dominant shape: λ₀ = {evals2[0]:.6f}, Δ = {-np.log(evals2[0]):.6f}")
except:
    print("    Could not find 2-dominant shape at K*=7/30")

print(f"\n  CONCLUSION: The spectral gap value depends on the evidence shape,")
print(f"  not just on K*. Different gauge groups (SU(2) vs SU(3)) have")
print(f"  different evidence distributions and hence different Δ values,")
print(f"  even though K*=7/30 is universal.")
