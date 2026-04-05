"""
Multi-site spectral gap: geometric approach
============================================
Per-step bounds in the Hilbert metric fail because:
1. The evidence is state-dependent (not a fixed linear map)
2. The Hilbert metric diverges when components approach zero
3. The evidence response at equilibrium is enormous (380x)

New approach: find the RIGHT metric and the RIGHT Lyapunov function.

Key observations:
- The state space B = {Born(Theta) >= 1/27} is COMPACT in CP^3
- The DS+floor dynamics is continuous on this compact set
- Global convergence is proven: all trajectories converge to m*
- The coupled system preserves spatial uniformity (if all sites equal, they stay equal)

Strategy: measure contraction in metrics that respect the compact geometry.
"""

import numpy as np
from numpy.linalg import norm, eigvalsh

# === Core DS machinery ===

def ds_combine(m, e):
    th_m, s1m, s2m, s3m = m
    th_e, s1e, s2e, s3e = e
    K = s1m*s2e + s1m*s3e + s2m*s1e + s2m*s3e + s3m*s1e + s3m*s2e
    if abs(K) >= 1.0:
        return m.copy()
    inv = 1.0 / (1.0 - K)
    return np.array([
        inv * th_m * th_e,
        inv * (s1m*s1e + s1m*th_e + th_m*s1e),
        inv * (s2m*s2e + s2m*th_e + th_m*s2e),
        inv * (s3m*s3e + s3m*th_e + th_m*s3e),
    ])

def born_floor(m, H=3):
    th = m[0]
    s = m[1:]
    s_sq = np.dot(s, s)
    total_sq = th**2 + s_sq
    if total_sq < 1e-30:
        return np.array([1.0, 0.0, 0.0, 0.0])
    born = th**2 / total_sq
    if born < 1.0/H**3:
        th_new = np.sqrt(s_sq / (H**3 - 1))
        m = np.array([th_new, m[1], m[2], m[3]])
    total = np.sum(np.abs(m))
    if total > 0:
        m = m / total
    return m

def ds_step(m, e, H=3):
    return born_floor(ds_combine(m, e), H)

# Fubini-Study distance on CP^3 (well-defined, bounded by pi/2)
def fs_distance(m1, m2):
    n1 = m1 / norm(m1)
    n2 = m2 / norm(m2)
    cos_d = np.clip(abs(np.dot(n1, n2)), 0, 1)
    return np.arccos(cos_d)

# Euclidean distance after L1 normalization
def l2_distance(m1, m2):
    return norm(m1 - m2)

# === Find the true equilibrium with K*=7/30 evidence ===
def find_equilibrium():
    """Find m* and e* at the K*=7/30 self-consistent equilibrium."""
    # Start from the known equilibrium structure:
    # dominant singleton ~0.596, two weak ~0.140, theta ~0.123
    m = np.array([0.123, 0.596, 0.140, 0.140])
    m /= np.sum(m)
    m = born_floor(m)

    # Self-consistent evidence that produces K*=7/30
    # From the paper: singleton concentration 93.2%, Born prob 0.932
    e = np.array([0.05, 0.85, 0.05, 0.05])
    e /= np.sum(e)
    e = born_floor(e)

    # Iterate to find fixed point
    for _ in range(200):
        m = ds_step(m, e)

    # Verify K*
    th_m, s1m, s2m, s3m = m
    th_e, s1e, s2e, s3e = e
    K = s1m*s2e + s1m*s3e + s2m*s1e + s2m*s3e + s3m*s1e + s3m*s2e
    return m, e, K

m_star, e_star, K_star = find_equilibrium()
print(f"Equilibrium: m* = {m_star}")
print(f"Evidence:    e* = {e_star}")
print(f"K* = {K_star:.6f} (target: {7/30:.6f})")

# ============================================================
# PART 1: Which metric contracts on the coupled system?
# ============================================================
print("\n" + "="*70)
print("PART 1: Contraction in different metrics")
print("="*70)

def coupled_step_ring(lattice):
    N = len(lattice)
    new = np.zeros_like(lattice)
    for i in range(N):
        e_avg = 0.5 * (lattice[(i-1) % N] + lattice[(i+1) % N])
        new[i] = ds_step(lattice[i], e_avg)
    return new

np.random.seed(42)
N = 16

# Two lattices starting from different random configs
lat1 = np.array([born_floor(np.random.dirichlet([2,2,2,2])) for _ in range(N)])
lat2 = np.array([born_floor(np.random.dirichlet([2,2,2,2])) for _ in range(N)])

print(f"\n1D ring N={N}: tracking distance between two trajectories")
print(f"{'step':>5} {'max_FS':>10} {'ratio_FS':>10} {'max_L2':>10} {'ratio_L2':>10} {'sum_L2':>10} {'ratio_sum':>10}")

prev_fs = prev_l2 = prev_sum = None
for step in range(60):
    max_fs = max(fs_distance(lat1[i], lat2[i]) for i in range(N))
    max_l2 = max(l2_distance(lat1[i], lat2[i]) for i in range(N))
    sum_l2 = sum(l2_distance(lat1[i], lat2[i]) for i in range(N))

    if step % 5 == 0:
        r_fs = max_fs / prev_fs if prev_fs and prev_fs > 1e-15 else float('nan')
        r_l2 = max_l2 / prev_l2 if prev_l2 and prev_l2 > 1e-15 else float('nan')
        r_sum = sum_l2 / prev_sum if prev_sum and prev_sum > 1e-15 else float('nan')
        print(f"{step:5d} {max_fs:10.6f} {r_fs:10.6f} {max_l2:10.6f} {r_l2:10.6f} {sum_l2:10.6f} {r_sum:10.6f}")

    prev_fs = max_fs
    prev_l2 = max_l2
    prev_sum = sum_l2

    lat1 = coupled_step_ring(lat1)
    lat2 = coupled_step_ring(lat2)

# ============================================================
# PART 2: Distance from equilibrium (not between two trajectories)
# ============================================================
print("\n" + "="*70)
print("PART 2: Convergence to uniform equilibrium")
print("="*70)

# A single lattice converging to uniform m*
np.random.seed(7)
lat = np.array([born_floor(np.random.dirichlet([2,2,2,2])) for _ in range(N)])

print(f"\n1D ring N={N}: distance from uniform m*")
print(f"{'step':>5} {'max_FS_to_m*':>14} {'ratio':>10} {'spatial_var':>14} {'ratio':>10}")

prev_d = prev_v = None
for step in range(80):
    max_d = max(fs_distance(lat[i], m_star) for i in range(N))
    # Spatial variance: how different are sites from each other?
    mean_m = np.mean(lat, axis=0)
    spatial_var = np.sqrt(sum(l2_distance(lat[i], mean_m)**2 for i in range(N)) / N)

    if step % 5 == 0:
        r_d = max_d / prev_d if prev_d and prev_d > 1e-15 else float('nan')
        r_v = spatial_var / prev_v if prev_v and prev_v > 1e-15 else float('nan')
        print(f"{step:5d} {max_d:14.6f} {r_d:10.6f} {spatial_var:14.8f} {r_v:10.6f}")

    prev_d = max_d
    prev_v = spatial_var
    lat = coupled_step_ring(lat)

# ============================================================
# PART 3: Two-step contraction in L2
# ============================================================
print("\n" + "="*70)
print("PART 3: Multi-step contraction")
print("="*70)

# Does the k-step map contract for some k?
np.random.seed(13)

def k_step_ring(lattice, k):
    for _ in range(k):
        lattice = coupled_step_ring(lattice)
    return lattice

for k in [1, 2, 3, 5, 10]:
    ratios = []
    for trial in range(500):
        lat_a = np.array([born_floor(np.random.dirichlet([2,2,2,2])) for _ in range(N)])
        lat_b = np.array([born_floor(np.random.dirichlet([2,2,2,2])) for _ in range(N)])

        d_before = sum(l2_distance(lat_a[i], lat_b[i])**2 for i in range(N))

        lat_a = k_step_ring(lat_a, k)
        lat_b = k_step_ring(lat_b, k)

        d_after = sum(l2_distance(lat_a[i], lat_b[i])**2 for i in range(N))

        if d_before > 1e-20:
            ratios.append(d_after / d_before)

    print(f"k={k:2d}: d²_after/d²_before  mean={np.mean(ratios):.6f}  "
          f"max={np.max(ratios):.6f}  P(ratio<1)={np.mean(np.array(ratios)<1):.3f}")

# ============================================================
# PART 4: Lyapunov function approach
# ============================================================
print("\n" + "="*70)
print("PART 4: Lyapunov function candidates")
print("="*70)

# Candidate 1: Total inter-site L2 variance
# V1 = sum_x ||m_x - m_bar||^2 where m_bar = spatial mean

# Candidate 2: Total FS distance from m*
# V2 = sum_x d_FS(m_x, m*)^2

# Candidate 3: Maximum component ratio spread across sites
# V3 = max_{x,y} d_H(m_x, m_y)

np.random.seed(99)

print(f"\nTracking three Lyapunov candidates on N={N} ring:")
print(f"{'step':>5} {'V1_spatial':>12} {'r1':>8} {'V2_to_m*':>12} {'r2':>8} {'V3_maxpair':>12} {'r3':>8}")

lat = np.array([born_floor(np.random.dirichlet([2,2,2,2])) for _ in range(N)])

prev1 = prev2 = prev3 = None
for step in range(80):
    mean_m = np.mean(lat, axis=0)
    V1 = sum(l2_distance(lat[i], mean_m)**2 for i in range(N))
    V2 = sum(fs_distance(lat[i], m_star)**2 for i in range(N))

    # V3: max pairwise Hilbert distance
    max_pair = 0
    for i in range(N):
        for j in range(i+1, N):
            d = fs_distance(lat[i], lat[j])
            max_pair = max(max_pair, d)
    V3 = max_pair

    if step % 5 == 0:
        r1 = V1 / prev1 if prev1 and prev1 > 1e-20 else float('nan')
        r2 = V2 / prev2 if prev2 and prev2 > 1e-20 else float('nan')
        r3 = V3 / prev3 if prev3 and prev3 > 1e-20 else float('nan')
        print(f"{step:5d} {V1:12.8f} {r1:8.5f} {V2:12.8f} {r2:8.5f} {V3:12.8f} {r3:8.5f}")

    prev1 = V1
    prev2 = V2
    prev3 = V3
    lat = coupled_step_ring(lat)

# ============================================================
# PART 5: Jacobian eigenvalues of the coupled map at equilibrium
# ============================================================
print("\n" + "="*70)
print("PART 5: Coupled Jacobian spectrum at uniform equilibrium")
print("="*70)

# At uniform equilibrium: all sites = m*, evidence = m* (since avg of neighbors = m*)
# Perturbation: delta m_x at each site
# The linearized map has block structure: local + coupling

# Compute single-site Jacobian by finite differences
eps = 1e-7

def compute_jacobian(m_base, e_base):
    """4x4 Jacobian d(ds_step)/dm at (m_base, e_base)"""
    J = np.zeros((4, 4))
    f0 = ds_step(m_base, e_base)
    for j in range(4):
        m_pert = m_base.copy()
        m_pert[j] += eps
        m_pert /= np.sum(m_pert)  # re-normalize
        m_pert = born_floor(m_pert)
        f1 = ds_step(m_pert, e_base)
        J[:, j] = (f1 - f0) / eps
    return J

def compute_jacobian_evidence(m_base, e_base):
    """4x4 Jacobian d(ds_step)/de at (m_base, e_base)"""
    J = np.zeros((4, 4))
    f0 = ds_step(m_base, e_base)
    for j in range(4):
        e_pert = e_base.copy()
        e_pert[j] += eps
        e_pert /= np.sum(e_pert)
        e_pert = born_floor(e_pert)
        f1 = ds_step(m_base, e_pert)
        J[:, j] = (f1 - f0) / eps
    return J

# At uniform equilibrium with self-evidence
m_eq_self = m_star.copy()
e_eq_self = m_star.copy()  # When all neighbors = m*, evidence = m*

J_local = compute_jacobian(m_eq_self, e_eq_self)
J_evidence = compute_jacobian_evidence(m_eq_self, e_eq_self)

print(f"\nSingle-site at self-evidence equilibrium:")
print(f"Local Jacobian eigenvalues: {sorted(abs(np.linalg.eigvals(J_local)), reverse=True)}")
print(f"Evidence Jacobian eigenvalues: {sorted(abs(np.linalg.eigvals(J_evidence)), reverse=True)}")

# The coupled system on a ring of N sites:
# delta m_x^{n+1} = J_local * delta m_x^n + (J_evidence / 2) * (delta m_{x-1}^n + delta m_{x+1}^n)
# This is a circulant system! Fourier transform diagonalizes it.
# In Fourier: delta_hat_k^{n+1} = [J_local + J_evidence * cos(2*pi*k/N)] * delta_hat_k^n

print(f"\nCoupled ring N={N}: Fourier block eigenvalues")
print(f"{'mode k':>8} {'cos(2πk/N)':>12} {'max |λ|':>10}")

max_total = 0
for k in range(N):
    cos_k = np.cos(2 * np.pi * k / N)
    # The coupled block for mode k
    M_k = J_local + J_evidence * cos_k
    evals = np.linalg.eigvals(M_k)
    max_eval = np.max(np.abs(evals))
    max_total = max(max_total, max_eval)
    if k <= N//2:
        print(f"{k:8d} {cos_k:12.6f} {max_eval:10.6f}")

print(f"\nMaximum |eigenvalue| across all modes: {max_total:.6f}")
print(f"All modes decay: {max_total < 1}")

# ============================================================
# PART 6: Same analysis with K*=7/30 evidence
# ============================================================
print("\n" + "="*70)
print("PART 6: Coupled Jacobian at K*=7/30 equilibrium (non-self evidence)")
print("="*70)

# The PHYSICAL equilibrium has m* ≠ e*
# But on a coupled lattice where evidence = neighbor average,
# at uniform equilibrium ALL sites = m* and evidence = m*.
# The K*=7/30 equilibrium requires EXTERNAL evidence from the gauge field.
#
# On the coupled lattice, the evidence at site x is the average of
# neighboring mass functions. At uniform equilibrium, this is m* itself.
# The self-evidence K is NOT 7/30 — it's whatever K(m*, m*) gives.

th, s1, s2, s3 = m_star
K_self = s1*s2 + s1*s3 + s2*s1 + s2*s3 + s3*s1 + s3*s2
print(f"Self-evidence K(m*, m*) = {K_self:.6f}")
print(f"This is NOT K* = 7/30 = {7/30:.6f}")
print(f"The gauge field provides ADDITIONAL evidence beyond neighbor coupling.")

# Actually, the physical picture is:
# At each site, evidence comes from TWO sources:
# 1. The gauge field conditional distribution C_beta (external, gives K*=7/30)
# 2. The spatial coupling (neighbor mass functions)
#
# The paper's Proposition 6.9 considers the FULL coupled system
# where evidence = f(neighbors, gauge field).
# At uniform equilibrium with gauge evidence e*:
# Each site sees e* from the gauge field.
# Each site sees m* from all neighbors.
# The combined evidence is some combination of e* and m*.

# For the pure spatial coupling (no external evidence):
print(f"\nPure spatial coupling (evidence = neighbor average):")
J_spatial = J_local + J_evidence  # cos(0)=1 for k=0 mode
evals_spatial = np.linalg.eigvals(J_spatial)
print(f"k=0 mode eigenvalues: {sorted(abs(evals_spatial), reverse=True)}")
print(f"(k=0 is the uniform perturbation — should have eigenvalue related to self-evidence)")

# ============================================================
# PART 7: The real question — eigenvalues of full coupled+gauge system
# ============================================================
print("\n" + "="*70)
print("PART 7: Full system with gauge evidence")
print("="*70)

# Model: at each step, site x receives evidence that is a mixture of
# (a) gauge-field evidence e* and (b) spatial neighbor average
#
# Simplest model: e_x = alpha * e* + (1-alpha) * avg_neighbors
# At alpha=1: pure gauge, no spatial coupling
# At alpha=0: pure spatial coupling
#
# Physical: the gauge field produces correlations between sites;
# the spatial coupling is the lattice propagation.

# With gauge evidence e*:
J_local_gauge = compute_jacobian(m_star, e_star)
J_ev_gauge = compute_jacobian_evidence(m_star, e_star)

print(f"With gauge evidence e*:")
print(f"Local Jacobian eigenvalues: {sorted(abs(np.linalg.eigvals(J_local_gauge)), reverse=True)}")
print(f"Evidence Jacobian eigenvalues: {sorted(abs(np.linalg.eigvals(J_ev_gauge)), reverse=True)}")

# For mixed evidence: e_x = alpha * e* + (1-alpha) * avg(neighbors)
# Jacobian wrt m_x: J_local_gauge (independent of alpha for the local part)
# Jacobian wrt m_y (neighbor): (1-alpha)/z * J_ev_gauge

for alpha in [0.0, 0.25, 0.5, 0.75, 1.0]:
    print(f"\nalpha={alpha:.2f}: e_x = {alpha}*e* + {1-alpha}*avg_neighbors")
    max_eval_all = 0
    for k in range(N):
        cos_k = np.cos(2 * np.pi * k / N)
        if alpha < 1:
            M_k = J_local_gauge + (1-alpha) * J_ev_gauge * cos_k
        else:
            M_k = J_local_gauge
        evals = np.linalg.eigvals(M_k)
        max_eval_all = max(max_eval_all, np.max(np.abs(evals)))
    print(f"  max |eigenvalue| = {max_eval_all:.6f}, contracts: {max_eval_all < 1}")

# ============================================================
# PART 8: Direct numerical eigenvalues of full 4N x 4N system
# ============================================================
print("\n" + "="*70)
print("PART 8: Full coupled system eigenvalues (direct)")
print("="*70)

for N_test in [4, 8, 16]:
    # Build the full 4N x 4N Jacobian at uniform equilibrium with self-evidence
    full_J = np.zeros((4*N_test, 4*N_test))

    for x in range(N_test):
        # Local contribution
        full_J[4*x:4*(x+1), 4*x:4*(x+1)] = J_local

        # Neighbor contributions (evidence Jacobian * 1/2)
        x_left = (x - 1) % N_test
        x_right = (x + 1) % N_test
        full_J[4*x:4*(x+1), 4*x_left:4*(x_left+1)] += 0.5 * J_evidence
        full_J[4*x:4*(x+1), 4*x_right:4*(x_right+1)] += 0.5 * J_evidence

    evals = np.linalg.eigvals(full_J)
    evals_abs = np.sort(np.abs(evals))[::-1]

    print(f"\nN={N_test}: Full {4*N_test}x{4*N_test} Jacobian")
    print(f"  Top 10 |eigenvalues|: {evals_abs[:10].round(6)}")
    print(f"  Max |eigenvalue|: {evals_abs[0]:.6f}")
    print(f"  All < 1: {evals_abs[0] < 1}")
    print(f"  # eigenvalues > 0.9: {np.sum(evals_abs > 0.9)}")
    print(f"  # eigenvalues > 0.5: {np.sum(evals_abs > 0.5)}")
