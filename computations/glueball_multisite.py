#!/usr/bin/env python3
"""
GLUEBALL SPECTRUM FROM MULTI-SITE STRIP TRANSFER MATRIX
=========================================================

The single-plaquette Jacobian (8x8 for SU(3)) gives 4 physical states —
the root algebra spectrum. The full glueball tower requires spatial depth.

A strip of N spatial sites, each carrying a rank-r internal (Cartan-coupled)
DS system, gives a transfer matrix of dimension 4rN x 4rN. The physical
eigenvalues of this matrix give O(rN) glueball masses.

Physical picture:
- Each spatial site carries the SU(3) root algebra (rank 2, 8 internal DOF)
- Neighboring spatial sites couple through DS evidence exchange
- The strip transfer matrix propagates states across one Euclidean time step
- Its eigenvalue spectrum = the glueball mass spectrum

For SU(3) with N spatial sites: Jacobian is 8N x 8N.
N=1: 4 physical states (already computed)
N=2: up to 8 physical states
N=3: up to 12 physical states
N=4: up to 16 physical states
...enough to fill the lattice spectrum.

The spatial coupling between sites uses the SAME DS evidence mechanism
that couples root subalgebras within a site, but now coupling site i
to site i+1 (nearest-neighbor on the strip).
"""

import numpy as np
from scipy.optimize import brentq

# ============================================================
# Framework constants
# ============================================================
H = 3
FLOOR = 1.0 / H**3
K_STAR = 7.0 / 30
g_internal = K_STAR   # Cartan coupling within a site
g_spatial = K_STAR     # spatial coupling between sites

# ============================================================
# Core DS dynamics
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

# ============================================================
# Single-site equilibrium
# ============================================================
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
# Cartan matrix for SU(3)
# ============================================================
def cartan_A(n):
    C = np.zeros((n, n))
    for i in range(n):
        C[i, i] = 2
    for i in range(n - 1):
        C[i, i + 1] = -1
        C[i + 1, i] = -1
    return C

cartan_su3 = cartan_A(2)
rank = 2

# ============================================================
# Multi-site strip: N spatial sites, each with rank-r internal structure
# ============================================================
# State layout: x = [site_0_root_0, site_0_root_1, site_1_root_0, site_1_root_1, ...]
# Each root has 4 components. Total dimension = 4 * rank * N_sites.

def multisite_step(x, N_sites, cartan, g_int, g_spat, e_base):
    """
    One step of the multi-site strip transfer matrix.

    Each site has `rank` root subalgebras, each with 4 DS components.

    Coupling structure:
    1. INTERNAL (Cartan): root i at site s couples to root j at site s
       via cartan[i,j] with strength g_int
    2. SPATIAL: root i at site s couples to root i at site s+1 and s-1
       (same root, neighboring site) with strength g_spat
    """
    r = cartan.shape[0]
    dim_site = 4 * r  # 8 for SU(3)

    # Parse state into sites and roots
    def get_root(site, root):
        return x[site * dim_site + root * 4 : site * dim_site + root * 4 + 4]

    out = np.zeros_like(x)

    for s in range(N_sites):
        for i in range(r):
            # Start with base evidence
            e_i = e_base.copy()

            # Internal coupling: other roots at same site
            for j in range(r):
                if i != j and abs(cartan[i, j]) > 0.01:
                    neighbor = get_root(s, j)
                    e_i = e_i + g_int * cartan[i, j] * (neighbor - 0.25)

            # Spatial coupling: same root at neighboring sites
            for ds in [-1, 1]:
                s_neighbor = s + ds
                if 0 <= s_neighbor < N_sites:
                    neighbor = get_root(s_neighbor, i)
                    # Spatial coupling: same root, neighboring site
                    # Use -1 off-diagonal (like A_n Cartan for the spatial chain)
                    e_i = e_i + g_spat * (-1.0) * (neighbor - 0.25)

            # Enforce positivity and normalization
            e_i = np.abs(e_i)
            e_i = np.maximum(e_i, 1e-10)
            e_i /= np.sum(e_i)

            # DS combine
            m_i = get_root(s, i)
            new_i = ds_combine(m_i, e_i)

            out[s * dim_site + i * 4 : s * dim_site + i * 4 + 4] = new_i

    return out

def find_multisite_eq(N_sites, cartan, g_int, g_spat, n_iter=20000):
    """Find equilibrium of the multi-site strip."""
    r = cartan.shape[0]
    dim = 4 * r * N_sites

    # Initialize all sites at single-site equilibrium
    x = np.zeros(dim)
    for s in range(N_sites):
        for i in range(r):
            x[s * 4 * r + i * 4 : s * 4 * r + i * 4 + 4] = m0.copy()

    for step in range(n_iter):
        x_new = multisite_step(x, N_sites, cartan, g_int, g_spat, e0)
        diff = np.max(np.abs(x_new - x))
        x = x_new
        if diff < 1e-14:
            return x, True, step
    return x, False, n_iter

def multisite_jacobian(x_eq, N_sites, cartan, g_int, g_spat, eps=1e-7):
    """Full Jacobian of the multi-site transfer map."""
    dim = len(x_eq)

    def F(x):
        return multisite_step(x, N_sites, cartan, g_int, g_spat, e0)

    J = np.zeros((dim, dim))
    for j in range(dim):
        xp = x_eq.copy(); xp[j] += eps
        xm = x_eq.copy(); xm[j] -= eps
        J[:, j] = (F(xp) - F(xm)) / (2 * eps)

    return J

def extract_spectrum(J, ref_delta=None):
    """Extract physical eigenvalues and mass ratios from Jacobian."""
    evals = np.linalg.eigvals(J)
    mags = np.abs(evals)

    # Physical: not near-zero and not near-unity
    physical = []
    for ev, mag in zip(evals, mags):
        if 1e-6 < mag < 1.0 - 1e-6:
            delta = -np.log(mag)
            physical.append((delta, mag, ev))

    physical.sort(key=lambda x: x[0])  # sort by mass (ascending)

    if ref_delta is None and physical:
        ref_delta = physical[0][0]

    return physical, ref_delta

# ============================================================
# MAIN COMPUTATION
# ============================================================
print("=" * 90)
print("GLUEBALL SPECTRUM FROM MULTI-SITE STRIP TRANSFER MATRIX")
print(f"Framework: H=3, K*=7/30, Born=1/27, g_int=g_spat={K_STAR:.6f}")
print(f"Gauge group: SU(3), rank={rank}, internal dim per site = {4*rank}")
print("=" * 90)

# Lattice QCD targets
lattice = [
    ('0++',   1.000, 0.00),
    ('2++',   1.40,  0.04),
    ('0-+',   1.50,  0.04),
    ('0++*',  1.56,  0.11),
    ('1+-',   1.74,  0.04),
    ('2-+',   1.82,  0.06),
    ('3++',   1.85,  0.05),
    ('0++**', 2.12,  0.10),
    ('1--',   2.13,  0.09),
    ('2--',   2.12,  0.10),
]

# Reference delta from single-site
ref_delta_global = None

for N in [1, 2, 3, 4, 5]:
    dim = 4 * rank * N
    print(f"\n\n{'#'*90}")
    print(f"# N = {N} spatial sites, Jacobian {dim}x{dim}")
    print(f"{'#'*90}")

    x_eq, conv, steps = find_multisite_eq(N, cartan_su3, g_internal, g_spatial)
    print(f"Converged: {conv} (steps: {steps})")

    # Check uniformity of equilibrium
    for s in range(N):
        site_state = x_eq[s*4*rank:(s+1)*4*rank]
        print(f"  Site {s}: {site_state[:4]}  {site_state[4:8]}")

    J = multisite_jacobian(x_eq, N, cartan_su3, g_internal, g_spatial)
    print(f"Jacobian dim: {J.shape[0]}x{J.shape[1]}")
    print(f"Spectral radius: {np.max(np.abs(np.linalg.eigvals(J))):.8f}")

    physical, ref_d = extract_spectrum(J)

    # Use N=1 lightest as global reference
    if N == 1 and physical:
        ref_delta_global = physical[0][0]

    n_phys = len(physical)
    print(f"Physical eigenvalues: {n_phys} (out of {dim})")

    if physical:
        print(f"\n{'Mode':>6s} {'|lambda|':>12s} {'Delta':>10s} {'m/m_0(local)':>12s} "
              f"{'m/m_0(global)':>13s} {'Re(lambda)':>12s} {'Im(lambda)':>12s}")
        print(f"{'─'*80}")

        local_ref = physical[0][0]
        for k, (delta, mag, ev) in enumerate(physical[:20]):  # cap at 20
            local_ratio = delta / local_ref if local_ref > 0 else 0
            global_ratio = delta / ref_delta_global if ref_delta_global and ref_delta_global > 0 else 0
            print(f"{k:6d} {mag:12.8f} {delta:10.6f} {local_ratio:12.4f} "
                  f"{global_ratio:13.4f} {ev.real:12.8f} {ev.imag:12.8f}")

    # Comparison with lattice
    if physical and ref_delta_global:
        print(f"\n--- Lattice comparison (ratios relative to lightest at N={N}) ---")
        local_ref = physical[0][0]
        ds_ratios = [(p[0] / local_ref, p[0], p[1]) for p in physical]

        print(f"{'Lattice':>10s} {'lat m/m0':>10s} {'Best DS':>10s} {'dev(sigma)':>10s} {'match':>6s}")
        print(f"{'─'*50}")

        for jpc, lat_r, lat_err in lattice:
            if ds_ratios:
                best = min(ds_ratios, key=lambda x: abs(x[0] - lat_r))
                dev = abs(best[0] - lat_r)
                sigma = dev / lat_err if lat_err > 0 else float('inf')
                match = "HIT" if sigma < 1.5 else ("near" if sigma < 3.0 else "no")
                print(f"{jpc:>10s} {lat_r:10.3f} {best[0]:10.4f} {sigma:10.1f} {match:>6s}")

print(f"\n{'='*90}")
print("SUMMARY: How mass ratios evolve with strip width N")
print(f"{'='*90}")

# Collect all spectra for summary
print(f"\n{'N':>3s} {'n_phys':>7s} {'lightest Delta':>14s} {'ratios (first 8)':>50s}")
print(f"{'─'*80}")

for N in [1, 2, 3, 4, 5]:
    dim = 4 * rank * N
    x_eq, conv, steps = find_multisite_eq(N, cartan_su3, g_internal, g_spatial)
    J = multisite_jacobian(x_eq, N, cartan_su3, g_internal, g_spatial)
    physical, ref_d = extract_spectrum(J)

    if physical:
        local_ref = physical[0][0]
        ratios = [p[0] / local_ref for p in physical[:8]]
        ratio_str = ", ".join(f"{r:.3f}" for r in ratios)
        print(f"{N:3d} {len(physical):7d} {local_ref:14.6f} {ratio_str}")
