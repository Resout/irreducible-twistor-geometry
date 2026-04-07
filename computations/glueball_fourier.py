#!/usr/bin/env python3
"""
============================================================================
GLUEBALL SPECTRUM VIA FOURIER-DECOMPOSED TRANSFER OPERATOR
============================================================================

Theorem 8.3: On a ring of N sites, the coupled Jacobian is block-circulant
and Fourier-decomposes into N independent blocks:

    M_k = J_m + J_e cos(2πk/N),   k = 0, ..., N-1

where:
    J_m = dΦ/dm|_{m*,e*}  — local Jacobian (4×4)
    J_e = dΦ/de|_{m*,e*}  — evidence Jacobian (4×4)

Each M_k is 4×4 with eigenvalues λ_i(k). The mass gap at momentum k is:
    Δ_i(k) = -ln|λ_i(k)|

The glueball mass spectrum consists of ALL (k, i) pairs.
Different k carry spatial momentum; different i carry fibre quantum numbers.

Quantum number assignment via SO(4) = (SU(2)_L × SU(2)_R)/Z₂:
    Eigenvector in θ direction → (1,1) → 0++ scalar
    Isotropic section perturbation → breathing mode → 0++* excited scalar
    Anisotropic section perturbation → (3,3) → 2++ tensor

Lattice QCD reference (Morningstar & Peardon 1999, Chen et al. 2006):
    m(0++)/m(0++) = 1.000  (reference)
    m(2++)/m(0++) = 1.40 ± 0.04
    m(0-+)/m(0++) = 1.50 ± 0.04
    m(0++*)/m(0++) = 1.56 ± 0.11
    m(1+-)/m(0++) = 1.74 ± 0.04
============================================================================
"""

import numpy as np
from scipy.optimize import brentq

# ============================================================
# SECTION 1: DS CORE
# ============================================================

H = 3
FLOOR = 1.0 / H**3  # 1/27
K_STAR = 7.0 / 30


def ds_combine(m, e):
    """DS combination + Born floor enforcement."""
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
            sc_f = (1.0 - tn) / ss
            out[:3] = abs_s * sc_f
            out[3] = tn
    return out


# ============================================================
# SECTION 2: EQUILIBRIUM
# ============================================================

def find_equilibrium():
    """Find (m*, e*) at K* = 7/30."""
    def K_at(p_dom):
        p_w = (1.0 - p_dom) / 2.0
        sc = 1.0 - FLOOR
        raw = np.array([np.sqrt(p_dom * sc), np.sqrt(p_w * sc),
                        np.sqrt(p_w * sc), np.sqrt(FLOOR)])
        e = raw / np.sum(raw)
        m = np.array([0.4, 0.2, 0.2, 0.2])
        for _ in range(10000):
            m2 = ds_combine(m, e)
            if np.max(np.abs(m2 - m)) < 1e-15:
                break
            m = m2
        s = m[:3]; se = e[:3]
        return sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)

    p_dom = brentq(lambda p: K_at(p) - K_STAR, 0.90, 0.96, xtol=1e-14)
    p_w = (1.0 - p_dom) / 2.0
    sc = 1.0 - FLOOR
    raw = np.array([np.sqrt(p_dom * sc), np.sqrt(p_w * sc),
                    np.sqrt(p_w * sc), np.sqrt(FLOOR)])
    e_star = raw / np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(10000):
        m2 = ds_combine(m, e_star)
        if np.max(np.abs(m2 - m)) < 1e-15:
            break
        m = m2
    return m, e_star


# ============================================================
# SECTION 3: JACOBIAN DECOMPOSITION  J_m and J_e
# ============================================================

def compute_jacobian(m_star, e_star, eps=1e-8):
    """Compute J_m = dΦ/dm and J_e = dΦ/de at (m*, e*).

    J_m: how the DS+floor output changes when the INPUT mass is perturbed
         (evidence held fixed at e*)
    J_e: how the DS+floor output changes when the EVIDENCE is perturbed
         (input mass held fixed at m*)

    Both are 4×4 matrices.
    """
    Jm = np.zeros((4, 4))
    Je = np.zeros((4, 4))

    f0 = ds_combine(m_star, e_star)

    # J_m: perturb m, hold e fixed
    for j in range(4):
        mp = m_star.copy(); mp[j] += eps
        mm = m_star.copy(); mm[j] -= eps
        fp = ds_combine(mp, e_star)
        fm = ds_combine(mm, e_star)
        Jm[:, j] = (fp - fm) / (2 * eps)

    # J_e: perturb e, hold m fixed
    for j in range(4):
        ep = e_star.copy(); ep[j] += eps
        em = e_star.copy(); em[j] -= eps
        fp = ds_combine(m_star, ep)
        fm = ds_combine(m_star, em)
        Je[:, j] = (fp - fm) / (2 * eps)

    return Jm, Je


# ============================================================
# SECTION 4: FOURIER TRANSFER OPERATOR
# ============================================================

def transfer_matrix(Jm, Je, k, N):
    """M_k = J_m + J_e · cos(2πk/N)  (Theorem 8.3)."""
    return Jm + Je * np.cos(2 * np.pi * k / N)


def compute_spectrum(Jm, Je, N=512):
    """Full glueball spectrum from Fourier modes k = 0, ..., N/2.

    Returns list of (delta, k, eigenvalue, eigenvector) sorted by mass.
    """
    states = []

    for k in range(N // 2 + 1):  # k and N-k give same |eigenvalues|
        Mk = transfer_matrix(Jm, Je, k, N)
        evals, evecs = np.linalg.eig(Mk)

        for i in range(4):
            mag = abs(evals[i])
            if 1e-10 < mag < 1.0 - 1e-10:
                delta = -np.log(mag)
                states.append({
                    'delta': delta,
                    'k': k,
                    'momentum': 2 * np.pi * k / N,
                    'eigenvalue': evals[i],
                    'mag': mag,
                    'eigvec': evecs[:, i].real,
                })

    states.sort(key=lambda x: x['delta'])
    return states


def classify_so4(eigvec):
    """Classify eigenvector by SO(4) content.

    Under Pauli embedding M = θI + s·σ:
      δm = (δs₁, δs₂, δs₃, δθ)

    (1,1) trace: pure δθ perturbation → 0++ scalar
    (3,3) symmetric traceless: anisotropic δs → 2++ tensor
    Mixed: various
    """
    v = np.abs(eigvec)
    norm2 = np.sum(v**2)
    if norm2 < 1e-30:
        return '?', 0, 0

    theta_frac = v[3]**2 / norm2
    section_frac = np.sum(v[:3]**2) / norm2

    # Section anisotropy
    if np.sum(v[:3]**2) > 1e-20:
        s_norm = v[:3] / np.sqrt(np.sum(v[:3]**2))
        anisotropy = 1.0 - 3 * np.min(s_norm**2)  # 0 = isotropic, 1 = single direction
    else:
        anisotropy = 0

    if theta_frac > 0.6:
        label = '0++'     # θ-dominated → scalar
    elif anisotropy > 0.5:
        label = '2++'     # anisotropic section → tensor
    elif anisotropy < 0.2:
        label = '0++*'    # isotropic section → excited scalar
    else:
        label = '0-+'     # mixed

    return label, theta_frac, section_frac


# ============================================================
# SECTION 5: MAIN COMPUTATION
# ============================================================

print("=" * 80)
print("GLUEBALL SPECTRUM VIA FOURIER TRANSFER OPERATOR")
print("Theorem 8.3: M_k = J_m + J_e cos(2πk/N)")
print("=" * 80)

# Step 1: Equilibrium
m_star, e_star = find_equilibrium()
K_check = sum(m_star[i] * e_star[j] for i in range(3) for j in range(3) if i != j)
print(f"\nm* = [{', '.join(f'{x:.8f}' for x in m_star)}]")
print(f"e* = [{', '.join(f'{x:.8f}' for x in e_star)}]")
print(f"K*  = {K_check:.10f}  (target: {K_STAR:.10f})")

# Step 2: Jacobian decomposition
Jm, Je = compute_jacobian(m_star, e_star)

print(f"\n{'─'*50}")
print(f"J_m (local Jacobian, 4×4):")
for row in Jm:
    print(f"  [{', '.join(f'{x:+.6f}' for x in row)}]")
print(f"  eigenvalues: {np.sort(np.abs(np.linalg.eigvals(Jm)))[::-1]}")

print(f"\nJ_e (evidence Jacobian, 4×4):")
for row in Je:
    print(f"  [{', '.join(f'{x:+.6f}' for x in row)}]")
print(f"  eigenvalues: {np.sort(np.abs(np.linalg.eigvals(Je)))[::-1]}")

# Sanity: J_m + J_e = single-site Jacobian
J_single = Jm + Je
evals_single = np.linalg.eigvals(J_single)
print(f"\nJ_m + J_e (k=0, single-site):")
print(f"  eigenvalues: {np.sort(np.abs(evals_single))[::-1]}")
print(f"  (expect: 0.2829, 0.2813, ~0, ~0)")

print(f"\n||J_m - J_e|| = {np.linalg.norm(Jm - Je):.6f}")
if np.linalg.norm(Jm - Je) > 0.01:
    print(f"  → J_m ≠ J_e: asymmetric coupling (physical equilibrium, e* ≠ m*)")

# Step 3: Dispersion relation
N = 512
print(f"\n{'='*80}")
print(f"DISPERSION RELATION (N = {N})")
print(f"{'='*80}")

print(f"\n{'k':>6s} {'p=2πk/N':>10s} {'cos(p)':>10s} "
      f"{'|λ₀|':>10s} {'Δ₀':>10s} {'|λ₁|':>10s} {'Δ₁':>10s}")
print("─" * 75)

k_samples = list(range(0, 33)) + [64, 128, N//4, N//2]
for k in sorted(set(k_samples)):
    if k > N // 2:
        continue
    Mk = transfer_matrix(Jm, Je, k, N)
    evals = np.linalg.eigvals(Mk)
    mags = np.sort(np.abs(evals))[::-1]

    p = 2 * np.pi * k / N
    cos_p = np.cos(p)

    entries = []
    for mag in mags[:2]:
        if mag > 1e-10 and mag < 2.0:
            entries.append((mag, -np.log(mag) if mag < 1 else -np.log(mag)))
        else:
            entries.append((0, float('inf')))

    print(f"{k:6d} {p:10.4f} {cos_p:10.4f} "
          f"{entries[0][0]:10.6f} {entries[0][1]:10.4f} "
          f"{entries[1][0]:10.6f} {entries[1][1]:10.4f}")

# Step 4: Full glueball spectrum
print(f"\n{'='*80}")
print(f"GLUEBALL MASS SPECTRUM")
print(f"{'='*80}")

states = compute_spectrum(Jm, Je, N)

if not states:
    print("No physical states found!")
else:
    delta_ref = states[0]['delta']
    print(f"\nLightest: Δ₀ = {delta_ref:.6f} at k={states[0]['k']} "
          f"(p = {states[0]['momentum']:.4f})")

    # Deduplicate near-identical mass ratios
    shown_ratios = []
    unique_states = []
    for s in states:
        ratio = s['delta'] / delta_ref
        if not any(abs(ratio - r) < 0.005 for r in shown_ratios):
            shown_ratios.append(ratio)
            unique_states.append(s)

    print(f"\n{'#':>4s} {'k':>5s} {'p':>8s} {'|λ|':>10s} {'Δ':>10s} "
          f"{'m/m₀':>8s} {'θ-frac':>8s} {'s-frac':>8s} {'J^PC':>6s}")
    print("─" * 75)

    for i, s in enumerate(unique_states[:20]):
        ratio = s['delta'] / delta_ref
        label, tf, sf = classify_so4(s['eigvec'])
        print(f"{i:4d} {s['k']:5d} {s['momentum']:8.4f} {s['mag']:10.6f} "
              f"{s['delta']:10.4f} {ratio:8.4f} {tf:8.3f} {sf:8.3f} {label:>6s}")

    # Step 5: Lattice comparison
    print(f"\n{'='*80}")
    print(f"COMPARISON WITH LATTICE QCD")
    print(f"{'='*80}")

    lattice = [
        ('0++',   1.000, 0.00),
        ('2++',   1.400, 0.04),
        ('0-+',   1.500, 0.04),
        ('0++*',  1.560, 0.11),
        ('1+-',   1.740, 0.04),
        ('2-+',   1.820, 0.06),
        ('3++',   1.850, 0.05),
    ]

    ds_mass_ratios = [s['delta'] / delta_ref for s in unique_states]

    print(f"\n{'State':>8s} {'Lattice':>8s} {'±':>6s} "
          f"{'Best DS':>10s} {'DS k':>6s} {'Dev':>8s}")
    print("─" * 55)

    for jpc, lat_r, lat_err in lattice:
        # Find closest DS ratio
        best_idx = min(range(len(ds_mass_ratios)),
                       key=lambda i: abs(ds_mass_ratios[i] - lat_r))
        best_r = ds_mass_ratios[best_idx]
        best_s = unique_states[best_idx]

        if lat_err > 0:
            sigma = abs(best_r - lat_r) / lat_err
            dev_str = f"{sigma:.1f}σ"
        else:
            dev_str = "ref"

        print(f"{jpc:>8s} {lat_r:8.3f} {lat_err:6.3f} "
              f"{best_r:10.4f} {best_s['k']:6d} {dev_str:>8s}")

print(f"\n{'='*80}")
print("DONE")
print(f"{'='*80}")
