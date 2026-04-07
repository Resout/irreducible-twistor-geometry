#!/usr/bin/env python3
"""
GLUEBALL MASS RATIOS FROM DS TRANSFER OPERATOR
================================================

Computes the full eigenvalue spectrum of the coupled DS transfer operator
for SU(N) gauge theories at the K*=7/30 equilibrium, and forms mass ratios
Delta_k / Delta_0 for comparison with lattice QCD glueball spectra.

The key test: dimensionless mass ratios require no scale-setting.
If these match lattice data, it's the Mercury perihelion of this framework.

Lattice QCD reference values (Morningstar & Peardon 1999, Chen et al. 2006):
  m(2++)/m(0++) = 1.40 +/- 0.04
  m(0-+)/m(0++) = 1.50 +/- 0.04
  m(0++*)/m(0++) = 1.56 +/- 0.11
  m(1+-)/m(0++) = 1.74 +/- 0.04
  m(2-+)/m(0++) = 1.82 +/- 0.06
  m(3++)/m(0++) = 1.85 +/- 0.05
"""

import numpy as np
from scipy.optimize import brentq

# ============================================================
# DS framework constants
# ============================================================
H = 3
FLOOR = 1.0 / H**3  # 1/27
K_STAR = 7.0 / 30
g_coupling = K_STAR  # coupling = K* = 7/30

# ============================================================
# Core DS dynamics
# ============================================================
def ds_combine(m, e):
    """One DS combination step with Born floor enforcement."""
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
    """Find the single-site DS equilibrium at K*=7/30."""
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
print(f"Single-site equilibrium: m* = {m0}")
print(f"Single-site evidence:    e* = {e0}")
K_check = sum(m0[i] * e0[j] for i in range(3) for j in range(3) if i != j)
print(f"K* check: {K_check:.10f} (target: {K_STAR:.10f})")

# ============================================================
# Cartan matrices
# ============================================================
def cartan_A(n):
    """A_n Cartan matrix (for SU(n+1))."""
    C = np.zeros((n, n))
    for i in range(n):
        C[i, i] = 2
    for i in range(n - 1):
        C[i, i + 1] = -1
        C[i + 1, i] = -1
    return C

def cartan_G2():
    return np.array([[2, -3], [-1, 2]], dtype=float)

def cartan_F4():
    return np.array([[2,-1,0,0],[-1,2,-2,0],[0,-1,2,-1],[0,0,-1,2]], dtype=float)

def cartan_E(n):
    C = np.zeros((n, n))
    for i in range(n):
        C[i, i] = 2
    links = {
        6: [(0,1),(1,2),(2,3),(3,4),(2,5)],
        7: [(0,1),(1,2),(2,3),(3,4),(4,5),(2,6)],
        8: [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(2,7)],
    }
    for i, j in links[n]:
        C[i, j] = -1
        C[j, i] = -1
    return C

# ============================================================
# Coupled dynamics
# ============================================================
def coupled_step(states, e_base, cartan, g_val):
    """One coupled DS step: each site's evidence is modified by neighbours via Cartan matrix."""
    r = len(states)
    new = []
    for i in range(r):
        e_i = e_base.copy()
        for j in range(r):
            if i != j and abs(cartan[i, j]) > 0.01:
                e_i = e_i + g_val * cartan[i, j] * (states[j] - 0.25)
        e_i = np.abs(e_i)
        e_i = np.maximum(e_i, 1e-10)
        e_i /= np.sum(e_i)
        new.append(ds_combine(states[i], e_i))
    return new

def find_coupled_eq(cartan, g_val, n_iter=15000):
    """Find coupled equilibrium by iteration."""
    r = cartan.shape[0]
    states = [m0.copy() for _ in range(r)]
    for step in range(n_iter):
        new = coupled_step(states, e0, cartan, g_val)
        diff = max(np.max(np.abs(new[i] - states[i])) for i in range(r))
        states = new
        if diff < 1e-15:
            return states, True
    return states, False

def full_jacobian(states, cartan, g_val, eps=1e-7):
    """Full coupled Jacobian: 4r x 4r matrix, ALL eigenvalues."""
    r = len(states)
    dim = 4 * r
    x0 = np.concatenate(states)

    def F(x):
        ss = [x[4*i:4*(i+1)] for i in range(r)]
        return np.concatenate(coupled_step(ss, e0, cartan, g_val))

    J = np.zeros((dim, dim))
    for j in range(dim):
        xp = x0.copy()
        xp[j] += eps
        xm = x0.copy()
        xm[j] -= eps
        J[:, j] = (F(xp) - F(xm)) / (2 * eps)

    return J

# ============================================================
# Compute full spectrum for a given group
# ============================================================
def glueball_spectrum(name, cartan, g_val):
    """Compute full eigenvalue spectrum and mass ratios."""
    states, conv = find_coupled_eq(cartan, g_val)
    J = full_jacobian(states, cartan, g_val)

    evals = np.linalg.eigvals(J)

    # Sort by magnitude (descending)
    mags = np.abs(evals)
    idx = np.argsort(mags)[::-1]
    evals_sorted = evals[idx]
    mags_sorted = mags[idx]

    # Mass gaps: Delta_k = -ln|lambda_k| for |lambda_k| > threshold
    threshold = 1e-10
    mass_gaps = []
    for k, (ev, mag) in enumerate(zip(evals_sorted, mags_sorted)):
        if mag > threshold and mag < 1.0 - 1e-10:  # exclude unit eigenvalues and zeros
            delta = -np.log(mag)
            mass_gaps.append((delta, ev, mag))

    # Sort by mass (ascending Delta = ascending mass)
    mass_gaps.sort(key=lambda x: x[0])

    return {
        'name': name,
        'converged': conv,
        'rank': cartan.shape[0],
        'dim_jacobian': 4 * cartan.shape[0],
        'eigenvalues': evals_sorted,
        'mass_gaps': mass_gaps,
    }

def print_spectrum(result):
    """Print the full glueball spectrum with mass ratios."""
    print(f"\n{'='*80}")
    print(f"GLUEBALL SPECTRUM: {result['name']}")
    print(f"Rank = {result['rank']}, Jacobian dim = {result['dim_jacobian']}, "
          f"Converged = {result['converged']}")
    print(f"{'='*80}")

    gaps = result['mass_gaps']
    if len(gaps) == 0:
        print("No non-trivial mass gaps found.")
        return

    delta_0 = gaps[0][0]  # lightest mass = reference

    print(f"\n{'Mode':>6s} {'|lambda|':>12s} {'Re(lambda)':>12s} {'Im(lambda)':>12s} "
          f"{'Delta':>10s} {'m/m_0':>8s}")
    print(f"{'─'*70}")

    for k, (delta, ev, mag) in enumerate(gaps):
        ratio = delta / delta_0 if delta_0 > 0 else 0
        print(f"{k:6d} {mag:12.8f} {ev.real:12.8f} {ev.imag:12.8f} "
              f"{delta:10.6f} {ratio:8.4f}")

    # Print mass ratios summary
    print(f"\n--- MASS RATIOS (relative to lightest) ---")
    print(f"{'Mode':>6s} {'m/m_0':>10s} {'Delta':>12s}")
    print(f"{'─'*35}")
    for k, (delta, ev, mag) in enumerate(gaps[:12]):
        ratio = delta / delta_0
        print(f"{k:6d} {ratio:10.4f} {delta:12.6f}")

    return gaps, delta_0


# ============================================================
# MAIN COMPUTATION
# ============================================================

print("\n" + "="*80)
print("GLUEBALL MASS RATIOS FROM DS TRANSFER OPERATOR")
print("Framework: H=3, K*=7/30, Born floor=1/27, g=K*=7/30")
print("="*80)

# --- SU(2): rank 1, Jacobian 4x4 ---
print("\n\n" + "#"*80)
print("# SU(2)  [rank 1]")
print("#"*80)
result_su2 = glueball_spectrum("SU(2)", cartan_A(1), g_coupling)
gaps_su2, d0_su2 = print_spectrum(result_su2)

# --- SU(3): rank 2, Jacobian 8x8 ---
print("\n\n" + "#"*80)
print("# SU(3)  [rank 2]  <-- THE MAIN PREDICTION")
print("#"*80)
result_su3 = glueball_spectrum("SU(3)", cartan_A(2), g_coupling)
gaps_su3, d0_su3 = print_spectrum(result_su3)

# --- SU(4): rank 3, Jacobian 12x12 ---
print("\n\n" + "#"*80)
print("# SU(4)  [rank 3]")
print("#"*80)
result_su4 = glueball_spectrum("SU(4)", cartan_A(3), g_coupling)
gaps_su4, d0_su4 = print_spectrum(result_su4)

# --- SU(5): rank 4, Jacobian 16x16 ---
print("\n\n" + "#"*80)
print("# SU(5)  [rank 4]")
print("#"*80)
result_su5 = glueball_spectrum("SU(5)", cartan_A(4), g_coupling)
gaps_su5, d0_su5 = print_spectrum(result_su5)

# --- G2: rank 2 ---
print("\n\n" + "#"*80)
print("# G2  [rank 2, exceptional]")
print("#"*80)
result_g2 = glueball_spectrum("G2", cartan_G2(), g_coupling)
gaps_g2, d0_g2 = print_spectrum(result_g2)

# --- F4: rank 4 ---
print("\n\n" + "#"*80)
print("# F4  [rank 4, exceptional, tightest]")
print("#"*80)
result_f4 = glueball_spectrum("F4", cartan_F4(), g_coupling)
gaps_f4, d0_f4 = print_spectrum(result_f4)

# ============================================================
# COMPARISON WITH LATTICE QCD
# ============================================================
print("\n\n" + "="*80)
print("COMPARISON WITH LATTICE QCD (SU(3))")
print("="*80)

lattice_ratios = {
    '0++':    (1.000, 0.00, 'scalar (reference)'),
    '2++':    (1.40,  0.04, 'tensor'),
    '0-+':    (1.50,  0.04, 'pseudoscalar'),
    '0++*':   (1.56,  0.11, 'excited scalar'),
    '1+-':    (1.74,  0.04, 'exotic vector'),
    '2-+':    (1.82,  0.06, ''),
    '3++':    (1.85,  0.05, ''),
    '0++**':  (2.12,  0.10, '2nd excited scalar'),
}

if gaps_su3:
    d0 = gaps_su3[0][0]
    ds_ratios = [g[0] / d0 for g in gaps_su3]

    print(f"\n{'Lattice J^PC':>14s} {'Lattice m/m0':>12s} {'DS mode':>8s} {'DS m/m0':>10s} {'Deviation':>10s}")
    print(f"{'─'*60}")

    for i, (jpc, (lat_ratio, lat_err, desc)) in enumerate(lattice_ratios.items()):
        if i < len(ds_ratios):
            ds_r = ds_ratios[i]
            dev = abs(ds_r - lat_ratio)
            sigma = dev / lat_err if lat_err > 0 else float('inf')
            print(f"{jpc:>14s} {lat_ratio:12.3f} {i:8d} {ds_r:10.4f} {sigma:8.1f}σ")
        else:
            print(f"{jpc:>14s} {lat_ratio:12.3f} {'---':>8s} {'---':>10s}")

# ============================================================
# EIGENVALUE STRUCTURE ANALYSIS
# ============================================================
print("\n\n" + "="*80)
print("EIGENVALUE STRUCTURE ANALYSIS")
print("="*80)

for result in [result_su2, result_su3, result_su4, result_su5, result_g2, result_f4]:
    evals = result['eigenvalues']
    mags = np.abs(evals)
    n_unit = np.sum(mags > 0.999)
    n_zero = np.sum(mags < 1e-10)
    n_physical = len(evals) - n_unit - n_zero

    # Real vs complex eigenvalues
    n_real = np.sum(np.abs(evals.imag) < 1e-10)
    n_complex = len(evals) - n_real

    print(f"\n{result['name']:>8s}: {len(evals)} eigenvalues = "
          f"{n_unit} unit + {n_physical} physical + {n_zero} zero | "
          f"{n_real} real + {n_complex} complex")

    physical = sorted([e for e in evals if 1e-10 < abs(e) < 0.999],
                      key=lambda x: -abs(x))
    if physical:
        print(f"         Physical eigenvalues (by |λ|):")
        for e in physical[:8]:
            print(f"           λ = {e.real:+10.6f} {e.imag:+10.6f}i  |λ| = {abs(e):.8f}  "
                  f"Δ = {-np.log(abs(e)):.6f}")
