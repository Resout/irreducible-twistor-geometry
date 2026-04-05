"""
Verify that Born floor IS the mechanism:
Compare eigenvalues WITH and WITHOUT floor at the SAME fixed point.
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

# K*=7/30 evidence
p_dom_exact = 0.9322756157
p_weak_exact = (1 - p_dom_exact) / 2
raw = np.array([np.sqrt(p_dom_exact*(1-FLOOR)), np.sqrt(p_weak_exact*(1-FLOOR)),
                np.sqrt(p_weak_exact*(1-FLOOR)), np.sqrt(FLOOR)])
e_exact = raw / np.sum(raw)

# Fixed point WITH floor
m = np.array([0.4, 0.2, 0.2, 0.2])
for _ in range(1000):
    m_new, _ = ds_combine(m, e_exact, apply_floor=True)
    if np.max(np.abs(m_new - m)) < 1e-15:
        break
    m = m_new
m_star_floor = m_new

# Is the floor actually active at m*?
born_theta = m_star_floor[3]**2 / np.sum(m_star_floor**2)
print(f"Born(θ) at m* = {born_theta:.8f}")
print(f"Floor = {FLOOR:.8f}")
print(f"Floor active? {born_theta <= FLOOR + 1e-8}")

# Jacobian WITH floor
eps = 1e-8
V = np.array([[1,0,0,-1], [0,1,0,-1], [0,0,1,-1]], dtype=float).T
VTV_inv = np.linalg.inv(V.T @ V)

def get_eigenvalues(m_fp, e, use_floor):
    J = np.zeros((4, 4))
    f0 = ds_combine(m_fp, e, apply_floor=use_floor)[0]
    for j in range(4):
        mp = m_fp.copy()
        mp[j] += eps
        fp = ds_combine(mp, e, apply_floor=use_floor)[0]
        J[:, j] = (fp - f0) / eps
    J_proj = VTV_inv @ V.T @ J @ V
    return np.sort(np.abs(np.linalg.eigvals(J_proj)))[::-1]

evals_with = get_eigenvalues(m_star_floor, e_exact, True)
evals_without = get_eigenvalues(m_star_floor, e_exact, False)

print(f"\nAt the SAME fixed point m*:")
print(f"  WITH floor:    λ = {evals_with}")
print(f"  WITHOUT floor: λ = {evals_without}")
print(f"\n  WITH floor:    λ₀ = {evals_with[0]:.8f}, Δ = {-np.log(evals_with[0]):.8f}")
if evals_without[0] < 1:
    print(f"  WITHOUT floor: λ₀ = {evals_without[0]:.8f}, Δ = {-np.log(evals_without[0]):.8f}")
else:
    print(f"  WITHOUT floor: λ₀ = {evals_without[0]:.8f} ≥ 1 — MARGINAL/UNSTABLE")

# Now find the no-floor fixed point with the same evidence
print("\n" + "=" * 60)
print("No-floor fixed point (same evidence)")
print("=" * 60)

m_nf = np.array([0.6, 0.15, 0.15, 0.1])
for i in range(5000):
    m_nf_new, _ = ds_combine(m_nf, e_exact, apply_floor=False)
    diff = np.max(np.abs(m_nf_new - m_nf))
    if diff < 1e-15:
        print(f"  Converged at step {i}")
        break
    m_nf = m_nf_new
    if i % 500 == 0:
        print(f"  Step {i}: m = {m_nf_new}, diff = {diff:.2e}")

print(f"  No-floor FP: {m_nf}")
_, K_nf = ds_combine(m_nf, e_exact, apply_floor=False)
print(f"  K* (no floor) = {K_nf:.10f}")
born_nf = m_nf**2 / np.sum(m_nf**2) if np.sum(m_nf**2) > 0 else np.zeros(4)
print(f"  Born: {born_nf}")

if np.sum(m_nf**2) > 1e-30:
    evals_nf = get_eigenvalues(m_nf, e_exact, False)
    print(f"  Eigenvalues: {evals_nf}")
    if evals_nf[0] >= 1 - 1e-6:
        print(f"  λ₀ ≈ 1 — MARGINAL (no gap without floor)")
    elif evals_nf[0] < 1:
        print(f"  λ₀ = {evals_nf[0]:.8f}, Δ = {-np.log(evals_nf[0]):.8f}")
else:
    print("  Fixed point degenerate — cannot compute eigenvalues")

# ============================================================
# Key question: at the k-singleton FP, does eigenvalue -> 1 without floor?
# ============================================================
print("\n" + "=" * 60)
print("K-singleton fixed point test")
print("=" * 60)

# The k-singleton FP has mass concentrated on k hypotheses.
# At k=1: m = (1, 0, 0, 0) — all mass on one hypothesis
# This is a fixed point of DS with any evidence (self-sorting).
# The eigenvalue here should be 1 without floor.

m_1sing = np.array([1 - 1e-10, 1e-11, 1e-11, 1e-10 - 2e-11])
# L1 check
m_1sing = m_1sing / np.sum(m_1sing)

# Evidence: use something nontrivial
e_test = np.array([0.5, 0.2, 0.2, 0.1])

# Iterate to approach a singleton FP
m_it = m_1sing.copy()
for i in range(100):
    m_it, _ = ds_combine(m_it, e_test, apply_floor=False)

print(f"  1-singleton FP (no floor): {m_it}")
evals_1s = get_eigenvalues(m_it, e_test, False)
print(f"  Eigenvalues: {evals_1s}")

# 2-singleton FP: start with mass on 2 hypotheses
e_2 = np.array([0.45, 0.45, 0.05, 0.05])
m_2s = np.array([0.45, 0.45, 0.05, 0.05])
for _ in range(500):
    m_new, _ = ds_combine(m_2s, e_2, apply_floor=False)
    if np.max(np.abs(m_new - m_2s)) < 1e-15:
        break
    m_2s = m_new

print(f"\n  2-singleton FP (no floor): {m_2s}")
evals_2s_nf = get_eigenvalues(m_2s, e_2, False)
print(f"  Eigenvalues (no floor): {evals_2s_nf}")

# Same FP with floor
m_2sf = np.array([0.45, 0.45, 0.05, 0.05])
for _ in range(500):
    m_new, _ = ds_combine(m_2sf, e_2, apply_floor=True)
    if np.max(np.abs(m_new - m_2sf)) < 1e-15:
        break
    m_2sf = m_new

print(f"  2-singleton FP (with floor): {m_2sf}")
evals_2s_f = get_eigenvalues(m_2sf, e_2, True)
print(f"  Eigenvalues (with floor): {evals_2s_f}")
