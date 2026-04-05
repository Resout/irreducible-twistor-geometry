"""
Focused check: does the DS Schwinger function decay exponentially?

The connected 2-point function centered at the EQUILIBRIUM value:
  S₂ᶜ(t) = ⟨[O(m₀) - O(m*)] · [O(Φᵗ(m₀)) - O(m*)]⟩

should decay as λ₀^t = e^{-Δt}.
"""

import numpy as np
from scipy.optimize import brentq

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
    se, phi = e[:3], e[3]
    s_new = s * se + s * phi + theta * se
    theta_new = theta * phi
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i != j)
    denom = 1.0 - K
    m_out = np.zeros(4)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom
    if apply_floor:
        born_theta = m_out[3]**2 / np.sum(m_out**2)
        if born_theta < FLOOR:
            m_out = enforce_floor(m_out)
    return m_out, K

def ds_step(m, e):
    return ds_combine(m, e, apply_floor=True)[0]

# K*=7/30 equilibrium
def K_at_equilibrium(p_dom_val):
    p_weak_val = (1 - p_dom_val) / 2
    scale = 1 - FLOOR
    raw = np.array([np.sqrt(p_dom_val*scale), np.sqrt(p_weak_val*scale),
                    np.sqrt(p_weak_val*scale), np.sqrt(FLOOR)])
    e = raw / np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(1000):
        m_new, _ = ds_combine(m, e)
        if np.max(np.abs(m_new - m)) < 1e-15: break
        m = m_new
    _, K = ds_combine(m, e, apply_floor=False)
    return K

p_dom = brentq(lambda p: K_at_equilibrium(p) - 7.0/30, 0.90, 0.95)
scale = 1 - FLOOR
raw = np.array([np.sqrt(p_dom*scale), np.sqrt((1-p_dom)/2*scale),
                np.sqrt((1-p_dom)/2*scale), np.sqrt(FLOOR)])
e_star = raw / np.sum(raw)
m = np.array([0.4, 0.2, 0.2, 0.2])
for _ in range(1000):
    m_new, _ = ds_combine(m, e_star)
    if np.max(np.abs(m_new - m)) < 1e-15: break
    m = m_new
m_star = m_new

print(f"Equilibrium: m* = [{m_star[0]:.6f}, {m_star[1]:.6f}, {m_star[2]:.6f}, {m_star[3]:.6f}]")

# Sample initial conditions
rng = np.random.default_rng(42)
n_samples = 50000
n_steps = 30

samples = []
while len(samples) < n_samples:
    m = rng.dirichlet([1, 1, 1, 1])
    born = m[3]**2 / np.sum(m**2)
    if born >= FLOOR:
        samples.append(m)
initial = np.array(samples)

# Observable
def O(m): return m[0]  # s₁

O_star = O(m_star)
print(f"O(m*) = {O_star:.6f}")

# Evolve and compute CORRECTLY CENTERED connected correlator
print(f"\n{'t':>4}  {'S₂ᶜ(t)':>14}  {'ratio(t/t-1)':>14}  {'λ₀^t':>14}  {'ratio/λ₀':>10}")

lambda_0 = 0.2829103473

current = initial.copy()
O_init_centered = np.array([O(m) - O_star for m in initial])  # O(m₀) - O(m*)

prev_S2c = None
for t in range(n_steps + 1):
    # O(Φᵗ(m₀)) - O(m*)
    O_t_centered = np.array([O(m) - O_star for m in current])

    # S₂ᶜ(t) = ⟨[O(m₀)-O(m*)] · [O(Φᵗ(m₀))-O(m*)]⟩
    S2c = np.mean(O_init_centered * O_t_centered)

    if t > 0 and prev_S2c != 0 and abs(prev_S2c) > 1e-15:
        ratio = S2c / prev_S2c
    else:
        ratio = float('nan')

    predicted = lambda_0**t
    ratio_over_lambda = ratio / lambda_0 if not np.isnan(ratio) and lambda_0 != 0 else float('nan')

    if t <= 15 or t % 5 == 0:
        print(f"{t:>4}  {S2c:>14.6e}  {ratio:>14.6f}  {predicted:>14.6e}  {ratio_over_lambda:>10.4f}")

    prev_S2c = S2c

    # Advance all trajectories one step
    if t < n_steps:
        current = np.array([ds_step(m, e_star) for m in current])

# Fit decay rate from first 10 non-zero points
print("\n--- Decay rate extraction ---")
current = initial.copy()
S2c_values = []
for t in range(n_steps + 1):
    O_t = np.array([O(m) - O_star for m in current])
    S2c = np.mean(O_init_centered * O_t)
    S2c_values.append(S2c)
    if t < n_steps:
        current = np.array([ds_step(m, e_star) for m in current])

valid = [(t, np.log(abs(S2c_values[t]))) for t in range(1, 15) if abs(S2c_values[t]) > 1e-15]
if len(valid) >= 3:
    ts = np.array([v[0] for v in valid])
    log_S = np.array([v[1] for v in valid])
    coeffs = np.polyfit(ts, log_S, 1)
    Delta_measured = -coeffs[0]
    Delta_exact = -np.log(lambda_0)
    print(f"Measured Δ = {Delta_measured:.6f}")
    print(f"Exact Δ   = {Delta_exact:.6f}")
    print(f"Error     = {abs(Delta_measured - Delta_exact)/Delta_exact*100:.2f}%")

    # Per-step ratio
    ratios = [S2c_values[t+1]/S2c_values[t] for t in range(1, 10)
              if abs(S2c_values[t]) > 1e-15]
    print(f"\nPer-step ratios S₂ᶜ(t+1)/S₂ᶜ(t):")
    for i, r in enumerate(ratios):
        print(f"  t={i+1}→{i+2}: {r:.6f}  (λ₀ = {lambda_0:.6f})")
