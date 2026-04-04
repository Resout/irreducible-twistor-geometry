"""
OS4 on R⁴: does the mass gap survive decompactification?

The question: on S⁴ (compact), the spectral gap Δ = 1.263 is proven.
On R⁴ = S⁴ \ {∞}, does the gap survive?

Standard approach: show the gap is uniform in volume.
DS approach: the gap is a fibre property — it doesn't reference the base.

This computation checks: for an N-site coupled system (mimicking
the decompactification limit), does the Koopman spectral gap
remain bounded below by a positive constant as N → ∞?

The multi-site proposition says the rate converges to λ₀ = 0.2817
at all N tested. Here we check the SCHWINGER FUNCTION decay
on the coupled system — does the CORRELATOR decay exponentially
with a rate that doesn't vanish as N grows?
"""

import numpy as np
from scipy.optimize import brentq

H, FLOOR = 3, 1.0/27

def enforce_floor(m):
    s = m[:3].copy()
    lo, hi = m[3], 1.0
    for _ in range(50):
        mid = (lo+hi)/2
        ss = (1-mid)/np.sum(s)*s if np.sum(s)>0 else s
        if mid**2/(np.sum(ss**2)+mid**2) < FLOOR: lo=mid
        else: hi=mid
    tn = (lo+hi)/2
    ss = (1-tn)/np.sum(s)*s if np.sum(s)>0 else s
    return np.array([ss[0],ss[1],ss[2],tn])

def ds_combine(m, e):
    s,th = m[:3],m[3]; se,ph = e[:3],e[3]
    sn = s*se + s*ph + th*se; tn = th*ph
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
    mo = np.append(sn,tn)/(1-K)
    if mo[3]**2/np.sum(mo**2) < FLOOR: mo = enforce_floor(mo)
    return mo

# Find K*=7/30 equilibrium
def K_eq(pd):
    pw=(1-pd)/2; sc=1-FLOOR
    raw=np.array([np.sqrt(pd*sc),np.sqrt(pw*sc),np.sqrt(pw*sc),np.sqrt(FLOOR)])
    e=raw/np.sum(raw); m=np.array([.4,.2,.2,.2])
    for _ in range(500):
        mn = ds_combine(m, e)
        if np.max(np.abs(mn-m))<1e-15: break
        m=mn
    s,th=m[:3],m[3]; se,ph=e[:3],e[3]
    return sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)

pd=brentq(lambda p: K_eq(p)-7/30, .9, .95)
sc=1-FLOOR; pw=(1-pd)/2
raw=np.array([np.sqrt(pd*sc),np.sqrt(pw*sc),np.sqrt(pw*sc),np.sqrt(FLOOR)])
e_star=raw/np.sum(raw)
m=np.array([.4,.2,.2,.2])
for _ in range(500):
    mn=ds_combine(m,e_star)
    if np.max(np.abs(mn-m))<1e-15: break
    m=mn
m_star=mn

print(f"Single-site equilibrium: m* = [{m_star[0]:.6f}, {m_star[1]:.6f}, {m_star[2]:.6f}, {m_star[3]:.6f}]")
print(f"λ₀ = 0.2829, Δ = 1.263\n")

# ============================================================
# Multi-site coupled system: ring of N sites
# ============================================================

def coupled_step(lattice, N):
    """One sweep of coupled DS on a ring of N sites."""
    new_lattice = np.zeros_like(lattice)
    for i in range(N):
        # Evidence = average of neighbors
        left = lattice[(i-1) % N]
        right = lattice[(i+1) % N]
        e = (left + right) / 2
        e = e / np.sum(e)  # L₁ normalize
        new_lattice[i] = ds_combine(lattice[i], e)
    return new_lattice

def find_coupled_equilibrium(N, n_sweeps=200):
    """Find the coupled equilibrium for a ring of N sites."""
    rng = np.random.default_rng(42 + N)
    lattice = np.zeros((N, 4))
    for i in range(N):
        m = rng.dirichlet([1,1,1,1])
        while m[3]**2/np.sum(m**2) < FLOOR:
            m = rng.dirichlet([1,1,1,1])
        lattice[i] = m

    for _ in range(n_sweeps):
        lattice = coupled_step(lattice, N)
    return lattice

print("=" * 70)
print("MULTI-SITE SCHWINGER FUNCTION DECAY")
print("=" * 70)

lambda_0_single = 0.2829103473

for N in [4, 8, 16, 32, 64]:
    print(f"\n--- N = {N} sites (ring) ---")

    # Find equilibrium
    eq = find_coupled_equilibrium(N, n_sweeps=300)
    eq_s1 = eq[:, 0]  # s₁ at each site

    # Check: are all sites at the same equilibrium?
    spread = np.std(eq_s1)
    print(f"  Equilibrium s₁ spread: {spread:.6e}")

    # Perturb site 0 and measure how the perturbation decays
    n_trials = 200
    n_steps = 30
    decay_rates = []

    for trial in range(n_trials):
        rng = np.random.default_rng(1000*N + trial)

        # Start at equilibrium, perturb site 0
        lattice = eq.copy()
        perturbation = rng.normal(0, 0.01, 4)
        perturbation[3] = -np.sum(perturbation[:3])  # keep L₁
        lattice[0] += perturbation
        # Re-enforce L₁ and positivity
        lattice[0] = np.abs(lattice[0])
        lattice[0] /= np.sum(lattice[0])
        if lattice[0][3]**2/np.sum(lattice[0]**2) < FLOOR:
            lattice[0] = enforce_floor(lattice[0])

        # Track deviation at site 0
        devs = []
        for t in range(n_steps):
            dev = np.linalg.norm(lattice[0] - eq[0])
            devs.append(dev)
            lattice = coupled_step(lattice, N)

        # Compute per-step decay ratio after transient settles
        if len(devs) > 5 and devs[3] > 1e-15:
            ratios = [devs[t+1]/devs[t] for t in range(3, min(10, len(devs)-1))
                      if devs[t] > 1e-15]
            if ratios:
                decay_rates.append(np.mean(ratios))

    if decay_rates:
        mean_rate = np.mean(decay_rates)
        std_rate = np.std(decay_rates)
        gap = -np.log(mean_rate) if mean_rate > 0 and mean_rate < 1 else float('nan')
        print(f"  Mean decay rate: {mean_rate:.6f} ± {std_rate:.6f}")
        print(f"  Gap: Δ = {gap:.4f}  (single-site: 1.263)")
        print(f"  Rate / λ₀: {mean_rate / lambda_0_single:.4f}  (should be ~1 if gap is preserved)")
    else:
        print(f"  Could not extract decay rate")


# ============================================================
# Spatial correlator: C(r) = ⟨O(0)O(r)⟩ - ⟨O⟩² at equilibrium
# ============================================================
print("\n" + "=" * 70)
print("SPATIAL CORRELATOR AT EQUILIBRIUM (N=64 ring)")
print("=" * 70)

N = 64
eq = find_coupled_equilibrium(N, n_sweeps=500)

# Perturb site 0, evolve, measure correlation at distance r
n_trials = 500
n_equil_steps = 5  # let system settle after perturbation

spatial_corr = np.zeros(N//2 + 1)
spatial_count = np.zeros(N//2 + 1)

for trial in range(n_trials):
    rng = np.random.default_rng(5000 + trial)

    lattice = eq.copy()
    # Small perturbation at site 0
    pert = rng.normal(0, 0.01, 4)
    pert[3] = -np.sum(pert[:3])
    lattice[0] += pert
    lattice[0] = np.abs(lattice[0])
    lattice[0] /= np.sum(lattice[0])
    if lattice[0][3]**2/np.sum(lattice[0]**2) < FLOOR:
        lattice[0] = enforce_floor(lattice[0])

    # Initial perturbation at each site
    delta_0 = lattice[:, 0] - eq[:, 0]  # s₁ deviation

    # Evolve a few steps
    for _ in range(n_equil_steps):
        lattice = coupled_step(lattice, N)

    # Deviation after evolution
    delta_t = lattice[:, 0] - eq[:, 0]

    # Spatial correlator: ⟨δ(0) · δ(r)⟩
    for r in range(N//2 + 1):
        spatial_corr[r] += delta_0[0] * delta_t[r]
        spatial_count[r] += 1

spatial_corr /= spatial_count

print(f"\n  {'r':>4}  {'C(r)':>14}  {'|C(r)|':>14}  {'ln|C(r)|':>12}")
for r in range(min(20, N//2 + 1)):
    if abs(spatial_corr[r]) > 1e-18:
        print(f"  {r:>4}  {spatial_corr[r]:>14.6e}  {abs(spatial_corr[r]):>14.6e}  {np.log(abs(spatial_corr[r])):>12.4f}")
    else:
        print(f"  {r:>4}  {spatial_corr[r]:>14.6e}  {'~0':>14}  {'---':>12}")

# Fit exponential decay
valid = [(r, np.log(abs(spatial_corr[r]))) for r in range(1, N//4)
         if abs(spatial_corr[r]) > 1e-18]
if len(valid) >= 3:
    rs = np.array([v[0] for v in valid[:10]])
    log_C = np.array([v[1] for v in valid[:10]])
    coeffs = np.polyfit(rs, log_C, 1)
    m_spatial = -coeffs[0]
    print(f"\n  Spatial mass: m = {m_spatial:.4f} per site")
    print(f"  (If this is positive, correlations decay exponentially in space)")
    print(f"  (This is the mass gap on the lattice)")
