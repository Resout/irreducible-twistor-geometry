#!/usr/bin/env python3
"""
GLUEBALL BAND STRUCTURE FROM FOURIER-DIAGONALISED TRANSFER MATRIX
===================================================================

The single-site coupled Jacobian gives the k=0 (spatially uniform) glueball
spectrum. The full spectrum requires the Fourier decomposition:

    M_k = J_m + J_e * cos(k)

where J_m is the mass-Jacobian (how the map depends on the state) and J_e is
the evidence-Jacobian (how the map depends on the neighbor's state used as
evidence). Each k in [0, pi] gives a different effective transfer matrix.

The band structure -- eigenvalues of M_k as a function of k -- gives the
full glueball spectrum. The minimum of each eigenvalue band is a glueball mass.

For SU(3): rank 2, internal Jacobian 8x8, Fourier bands sweep from k=0 to k=pi.

Lattice QCD targets (Morningstar-Peardon 1999, Chen et al. 2006):
  0++:   1.000 (reference)
  2++:   1.40 +/- 0.04
  0-+:   1.50 +/- 0.04
  0++*:  1.56 +/- 0.11
  1+-:   1.74 +/- 0.04
  2-+:   1.82 +/- 0.06
  3++:   1.85 +/- 0.05
"""

import numpy as np
from scipy.optimize import brentq

# ============================================================
# DS framework constants
# ============================================================
H = 3
FLOOR = 1.0 / H**3
K_STAR = 7.0 / 30
g_coupling = K_STAR

# ============================================================
# Core DS dynamics (same as glueball_mass_ratios.py)
# ============================================================
def ds_combine(m, e):
    s, th = m[:3], m[3]
    se, ph = e[:3], e[3]
    s_new = s * se + s * ph + th * se
    th_new = th * ph
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    d = 1.0 - K
    if abs(d) < 1e-15:
        return m.copy()
    out = np.zeros(4)
    out[:3] = s_new / d
    out[3] = th_new / d
    born = out[3]**2 / np.sum(out**2) if np.sum(out**2) > 0 else 1.0
    if born < FLOOR:
        abs_s = np.abs(out[:3])
        ss = np.sum(abs_s)
        sq = np.sum(abs_s**2)
        if ss > 1e-15:
            r = sq / ss**2
            a_c = 26.0 - r
            b_c = 2.0 * r
            c_c = -r
            disc = b_c**2 - 4 * a_c * c_c
            tn = (-b_c + np.sqrt(disc)) / (2 * a_c)
            sc = (1.0 - tn) / ss
            out[:3] = abs_s * sc
            out[3] = tn
    return out

def find_single_eq():
    def K_at(p_dom):
        p_w = (1.0 - p_dom) / 2.0
        sc = 1.0 - FLOOR
        raw = np.array([np.sqrt(p_dom * sc), np.sqrt(p_w * sc),
                        np.sqrt(p_w * sc), np.sqrt(FLOOR)])
        e = raw / np.sum(raw)
        m = np.array([0.4, 0.2, 0.2, 0.2])
        for _ in range(8000):
            m2 = ds_combine(m, e)
            if np.max(np.abs(m2 - m)) < 1e-15:
                break
            m = m2
        s, th = m[:3], m[3]
        se, ph = e[:3], e[3]
        return sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    p_dom = brentq(lambda p: K_at(p) - K_STAR, 0.90, 0.96, xtol=1e-14)
    p_w = (1.0 - p_dom) / 2.0
    sc = 1.0 - FLOOR
    raw = np.array([np.sqrt(p_dom * sc), np.sqrt(p_w * sc),
                    np.sqrt(p_w * sc), np.sqrt(FLOOR)])
    e_star = raw / np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(8000):
        m2 = ds_combine(m, e_star)
        if np.max(np.abs(m2 - m)) < 1e-15:
            break
        m = m2
    return m, e_star

m0, e0 = find_single_eq()

# ============================================================
# Cartan matrices
# ============================================================
def cartan_A(n):
    C = np.zeros((n, n))
    for i in range(n):
        C[i, i] = 2
    for i in range(n - 1):
        C[i, i + 1] = -1
        C[i + 1, i] = -1
    return C

# ============================================================
# SEPARATE J_m and J_e computation
# ============================================================
def coupled_step_full(states, evidences, cartan, g_val):
    """Coupled step with explicit evidence array (not tied to states)."""
    r = len(states)
    new = []
    for i in range(r):
        e_i = evidences[i].copy()
        for j in range(r):
            if i != j and abs(cartan[i, j]) > 0.01:
                e_i = e_i + g_val * cartan[i, j] * (states[j] - 0.25)
        e_i = np.abs(e_i)
        e_i = np.maximum(e_i, 1e-10)
        e_i /= np.sum(e_i)
        new.append(ds_combine(states[i], e_i))
    return new

def find_coupled_eq(cartan, g_val, n_iter=15000):
    r = cartan.shape[0]
    states = [m0.copy() for _ in range(r)]
    for step in range(n_iter):
        new = coupled_step_full(states, [e0.copy() for _ in range(r)], cartan, g_val)
        diff = max(np.max(np.abs(new[i] - states[i])) for i in range(r))
        states = new
        if diff < 1e-15:
            return states, True
    return states, False

def compute_Jm_Je(states_eq, cartan, g_val, eps=1e-7):
    """
    Compute J_m (Jacobian w.r.t. own state) and J_e (Jacobian w.r.t. neighbor states)
    separately at equilibrium.

    For the coupled map F(m_1,...,m_r) with evidence e0:
      J_m = dF/d(own state), holding neighbor states fixed
      J_e = dF/d(neighbor states), holding own state fixed

    For the Fourier decomposition:
      M_k = J_m + J_e * cos(k)

    J_m is block-diagonal (each site's self-Jacobian)
    J_e encodes the inter-site coupling
    """
    r = len(states_eq)
    dim = 4 * r
    x0 = np.concatenate(states_eq)

    def F(x):
        ss = [x[4*i:4*(i+1)] for i in range(r)]
        return np.concatenate(coupled_step_full(ss, [e0.copy() for _ in range(r)], cartan, g_val))

    # Full Jacobian
    J_full = np.zeros((dim, dim))
    for j in range(dim):
        xp = x0.copy(); xp[j] += eps
        xm = x0.copy(); xm[j] -= eps
        J_full[:, j] = (F(xp) - F(xm)) / (2 * eps)

    # J_m: block-diagonal part (site i depends on site i)
    J_m = np.zeros((dim, dim))
    for i in range(r):
        J_m[4*i:4*(i+1), 4*i:4*(i+1)] = J_full[4*i:4*(i+1), 4*i:4*(i+1)]

    # J_e: off-block-diagonal part (site i depends on site j, j != i)
    J_e = J_full - J_m

    return J_full, J_m, J_e

def fourier_spectrum(J_m, J_e, n_k=200):
    """
    Compute eigenvalues of M_k = J_m + J_e * cos(k) for k in [0, pi].
    Returns array of shape (n_k, dim) with eigenvalue magnitudes.
    """
    dim = J_m.shape[0]
    k_vals = np.linspace(0, np.pi, n_k)
    all_evals = np.zeros((n_k, dim))
    all_evals_complex = np.zeros((n_k, dim), dtype=complex)

    for ik, k in enumerate(k_vals):
        M_k = J_m + J_e * np.cos(k)
        evals = np.linalg.eigvals(M_k)
        # Sort by magnitude descending
        idx = np.argsort(np.abs(evals))[::-1]
        all_evals[ik, :] = np.abs(evals[idx])
        all_evals_complex[ik, :] = evals[idx]

    return k_vals, all_evals, all_evals_complex

# ============================================================
# SO(4) decomposition of eigenmodes
# ============================================================
def so4_decompose(eigvec, r):
    """
    Decompose an eigenmode of the 4r-dimensional Jacobian under SO(4).

    Each site has 4 components: (s1, s2, s3, theta).
    Under SO(4) = SU(2)_L x SU(2)_R:
      - (s1, s2, s3) transform as the adjoint of SU(2) = (3,1)+(1,3) (gauge)
      - theta transforms as (1,1) (scalar)

    The (3,3) graviton sector requires CROSS-SITE coupling of the s-components,
    which appears only when the eigenmode has nontrivial spatial structure (k != 0).

    Returns: fraction of eigenmode in scalar (theta) vs gauge (s) sectors.
    """
    # For each site, compute the fraction in s vs theta
    s_weight = 0.0
    th_weight = 0.0
    for i in range(r):
        block = eigvec[4*i:4*(i+1)]
        s_weight += np.sum(np.abs(block[:3])**2)
        th_weight += np.abs(block[3])**2

    total = s_weight + th_weight
    if total < 1e-20:
        return 0.0, 0.0
    return s_weight / total, th_weight / total

# ============================================================
# MAIN COMPUTATION: SU(3) BAND STRUCTURE
# ============================================================
print("=" * 80)
print("GLUEBALL BAND STRUCTURE: SU(3)")
print("Fourier decomposition M_k = J_m + J_e cos(k)")
print("=" * 80)

cartan_su3 = cartan_A(2)
states_eq, conv = find_coupled_eq(cartan_su3, g_coupling)
print(f"\nEquilibrium converged: {conv}")
for i, s in enumerate(states_eq):
    print(f"  Site {i}: m = {s}")

J_full, J_m, J_e = compute_Jm_Je(states_eq, cartan_su3, g_coupling)

print(f"\nJ_full spectral radius: {np.max(np.abs(np.linalg.eigvals(J_full))):.8f}")
print(f"J_m spectral radius:    {np.max(np.abs(np.linalg.eigvals(J_m))):.8f}")
print(f"J_e spectral radius:    {np.max(np.abs(np.linalg.eigvals(J_e))):.8f}")
print(f"J_e / J_m ratio:        {np.max(np.abs(np.linalg.eigvals(J_e))) / np.max(np.abs(np.linalg.eigvals(J_m))):.4f}")

# Compute band structure
k_vals, bands, bands_complex = fourier_spectrum(J_m, J_e, n_k=500)

print(f"\n{'='*80}")
print("BAND STRUCTURE: eigenvalue magnitudes at key k values")
print(f"{'='*80}")

# Print at selected k values
for k_target in [0.0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, 2*np.pi/3, np.pi]:
    ik = np.argmin(np.abs(k_vals - k_target))
    k = k_vals[ik]
    evals_k = bands[ik]

    # Filter physical eigenvalues
    physical = [(mag, bands_complex[ik, j]) for j, mag in enumerate(evals_k)
                if 1e-8 < mag < 1.0 - 1e-8]
    physical.sort(key=lambda x: -x[0])

    print(f"\nk = {k:.4f} (k/pi = {k/np.pi:.3f}):")
    for j, (mag, ev) in enumerate(physical):
        delta = -np.log(mag)
        print(f"  mode {j}: |lambda| = {mag:.8f}  Delta = {delta:.6f}  "
              f"lambda = {ev.real:+.6f}{ev.imag:+.6f}i")

# ============================================================
# BAND MINIMA AND MASS RATIOS
# ============================================================
print(f"\n{'='*80}")
print("BAND MINIMA (mass = min_k Delta_k for each band)")
print(f"{'='*80}")

# For each eigenvalue band, find the minimum Delta (= maximum |lambda|)
n_bands = bands.shape[1]
band_maxlam = np.zeros(n_bands)  # max |lambda| across k for each band
band_k_at_max = np.zeros(n_bands)  # k where max occurs
band_mindelta = np.zeros(n_bands)  # min Delta = -ln(max |lambda|)

for b in range(n_bands):
    valid = bands[:, b] > 1e-8  # exclude near-zero bands
    if np.any(valid):
        max_idx = np.argmax(bands[valid, b])
        k_indices = np.where(valid)[0]
        band_maxlam[b] = bands[k_indices[max_idx], b]
        band_k_at_max[b] = k_vals[k_indices[max_idx]]
    else:
        band_maxlam[b] = 0.0
        band_k_at_max[b] = 0.0

# Filter physical bands
physical_bands = [(b, band_maxlam[b], band_k_at_max[b])
                  for b in range(n_bands)
                  if 1e-6 < band_maxlam[b] < 1.0 - 1e-6]
physical_bands.sort(key=lambda x: -x[1])  # sort by max |lambda| descending (= lightest first)

# Reference mass = lightest band
if physical_bands:
    ref_delta = -np.log(physical_bands[0][1])

    print(f"\n{'Band':>6s} {'max|lambda|':>12s} {'k_max':>8s} {'k_max/pi':>8s} "
          f"{'Delta_min':>10s} {'m/m_0':>8s}")
    print(f"{'─'*60}")

    all_ratios = []
    for rank_idx, (b, maxlam, k_max) in enumerate(physical_bands):
        delta_min = -np.log(maxlam)
        ratio = delta_min / ref_delta
        all_ratios.append((ratio, delta_min, k_max, b))
        print(f"{rank_idx:6d} {maxlam:12.8f} {k_max:8.4f} {k_max/np.pi:8.3f} "
              f"{delta_min:10.6f} {ratio:8.4f}")

# ============================================================
# COMPARISON WITH LATTICE QCD
# ============================================================
print(f"\n{'='*80}")
print("COMPARISON WITH LATTICE QCD (SU(3))")
print(f"{'='*80}")

lattice = [
    ('0++',   1.000, 0.00),
    ('2++',   1.40,  0.04),
    ('0-+',   1.50,  0.04),
    ('0++*',  1.56,  0.11),
    ('1+-',   1.74,  0.04),
    ('2-+',   1.82,  0.06),
    ('3++',   1.85,  0.05),
    ('0++**', 2.12,  0.10),
]

if physical_bands:
    ds_ratios = sorted(all_ratios, key=lambda x: x[0])

    print(f"\n{'Lattice':>10s} {'lat m/m0':>10s} {'DS mode':>8s} {'DS m/m0':>10s} "
          f"{'k_max/pi':>8s} {'dev(sigma)':>10s} {'match':>6s}")
    print(f"{'─'*70}")

    for jpc, lat_r, lat_err in lattice:
        # Find closest DS ratio
        if ds_ratios:
            best_idx = min(range(len(ds_ratios)), key=lambda i: abs(ds_ratios[i][0] - lat_r))
            ds_r, ds_delta, ds_k, ds_band = ds_ratios[best_idx]
            dev = abs(ds_r - lat_r)
            sigma = dev / lat_err if lat_err > 0 else float('inf')
            match = "YES" if sigma < 2.0 else "no"
            print(f"{jpc:>10s} {lat_r:10.3f} {best_idx:8d} {ds_r:10.4f} "
                  f"{ds_k/np.pi:8.3f} {sigma:10.1f} {match:>6s}")
        else:
            print(f"{jpc:>10s} {lat_r:10.3f} {'---':>8s}")

# ============================================================
# BAND STRUCTURE VISUALIZATION (text-based)
# ============================================================
print(f"\n{'='*80}")
print("BAND STRUCTURE: Delta(k) for each physical band")
print("(lower Delta = lighter mass)")
print(f"{'='*80}")

# Print band structure at 20 k points
k_sample = np.linspace(0, np.pi, 20)
print(f"\n{'k/pi':>6s}", end="")
for rank_idx in range(min(8, len(physical_bands))):
    print(f" {'band'+str(rank_idx):>10s}", end="")
print()
print("─" * (6 + 11 * min(8, len(physical_bands))))

for k_target in k_sample:
    ik = np.argmin(np.abs(k_vals - k_target))
    print(f"{k_vals[ik]/np.pi:6.3f}", end="")

    # Get eigenvalues at this k, sorted by magnitude
    evals_k = bands[ik]
    physical_k = sorted([mag for mag in evals_k if 1e-8 < mag < 1.0 - 1e-8], reverse=True)

    for j in range(min(8, len(physical_k))):
        delta = -np.log(physical_k[j])
        print(f" {delta:10.4f}", end="")
    print()

# ============================================================
# EIGENMODE ANALYSIS AT k=pi/2 (WHERE NEW STATES EMERGE)
# ============================================================
print(f"\n{'='*80}")
print("EIGENMODE ANALYSIS AT k = pi/2")
print(f"{'='*80}")

ik_half = np.argmin(np.abs(k_vals - np.pi/2))
M_half = J_m + J_e * np.cos(np.pi/2)  # = J_m (since cos(pi/2) = 0)
M_pi = J_m + J_e * np.cos(np.pi)      # = J_m - J_e (alternating mode)
M_0 = J_m + J_e                        # = J_full (uniform mode)

for label, M in [("k=0 (M_0=J_m+J_e)", M_0),
                  ("k=pi/2 (M=J_m)", M_half),
                  ("k=pi (M=J_m-J_e)", M_pi)]:
    evals = np.linalg.eigvals(M)
    evals_phys = sorted([e for e in evals if 1e-8 < abs(e) < 1-1e-8],
                        key=lambda x: -abs(x))

    # Also compute eigenvectors for SO(4) decomposition
    evals_full, evecs_full = np.linalg.eig(M)

    print(f"\n{label}:")
    for j, e in enumerate(evals_phys[:8]):
        delta = -np.log(abs(e))

        # Find corresponding eigenvector
        idx = np.argmin(np.abs(np.abs(evals_full) - abs(e)))
        evec = evecs_full[:, idx]
        s_frac, th_frac = so4_decompose(evec, cartan_su3.shape[0])

        print(f"  mode {j}: |lambda|={abs(e):.8f} Delta={delta:.6f} "
              f"m/m0={delta/(-np.log(max(abs(ep) for ep in evals_phys))):.4f} "
              f"gauge_frac={s_frac:.3f} scalar_frac={th_frac:.3f} "
              f"lambda={e.real:+.6f}{e.imag:+.6f}i")

# ============================================================
# WORST-CASE k=pi MODE (alternating, highest mass states)
# ============================================================
print(f"\n{'='*80}")
print("k=pi MODE (alternating pattern = heaviest states)")
print(f"{'='*80}")

evals_pi = np.linalg.eigvals(M_pi)
evals_pi_phys = sorted([e for e in evals_pi if 1e-8 < abs(e) < 1-1e-8],
                        key=lambda x: -abs(x))

if evals_pi_phys:
    ref = -np.log(abs(evals_pi_phys[0]))
    for j, e in enumerate(evals_pi_phys[:8]):
        delta = -np.log(abs(e))
        print(f"  mode {j}: |lambda|={abs(e):.8f} Delta={delta:.6f} "
              f"m/m_0(k=pi)={delta/ref:.4f}")

print(f"\n{'='*80}")
print("DONE")
print(f"{'='*80}")
