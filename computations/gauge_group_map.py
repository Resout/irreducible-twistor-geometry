#!/usr/bin/env python3
"""
GAUGE GROUP -> MODULI CURVE MAP
================================
Derive the analytical relationship between gauge group G and
mass gap Delta(G) in the DS framework.

Key discovery: DS spectral radius is invariant under Dynkin diagram folding.
This reduces the full classification to ADE (simply-laced) algebras.
"""

import numpy as np
from scipy.optimize import brentq

# ============================================================
# DS CORE
# ============================================================
H = 3
FLOOR = 1.0 / H**3
g = 7.0 / 30

def ds_combine(m, e):
    s, th = m[:3], m[3]
    se, ph = e[:3], e[3]
    s_new = s * se + s * ph + th * se
    th_new = th * ph
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    d = 1.0 - K
    if abs(d) < 1e-15: return m.copy()
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
            disc = b_c**2 - 4*a_c*c_c
            tn = (-b_c + np.sqrt(disc)) / (2*a_c)
            sc = (1.0 - tn) / ss
            out[:3] = abs_s * sc
            out[3] = tn
    return out

def find_equilibrium():
    def K_at(p_dom):
        p_w = (1.0 - p_dom) / 2.0
        sc = 1.0 - FLOOR
        raw = np.array([np.sqrt(p_dom*sc), np.sqrt(p_w*sc),
                        np.sqrt(p_w*sc), np.sqrt(FLOOR)])
        e = raw / np.sum(raw)
        m = np.array([0.4, 0.2, 0.2, 0.2])
        for _ in range(8000):
            m2 = ds_combine(m, e)
            if np.max(np.abs(m2 - m)) < 1e-15: break
            m = m2
        s, th = m[:3], m[3]
        se, ph = e[:3], e[3]
        return sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
    p_dom = brentq(lambda p: K_at(p) - 7.0/30, 0.90, 0.96, xtol=1e-14)
    p_w = (1.0 - p_dom) / 2.0
    sc = 1.0 - FLOOR
    raw = np.array([np.sqrt(p_dom*sc), np.sqrt(p_w*sc),
                    np.sqrt(p_w*sc), np.sqrt(FLOOR)])
    e_star = raw / np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(8000):
        m2 = ds_combine(m, e_star)
        if np.max(np.abs(m2 - m)) < 1e-15: break
        m = m2
    return m, e_star

m0, e0 = find_equilibrium()
print(f"Equilibrium found: m* = {m0}")
print(f"                   e* = {e0}")

# ============================================================
# SINGLE-SITE AND CROSS JACOBIANS
# ============================================================
eps = 1e-8

# Single-site Jacobian (no coupling)
J_uncoupled = np.zeros((4, 4))
for j in range(4):
    mp = m0.copy(); mp[j] += eps
    mm = m0.copy(); mm[j] -= eps
    J_uncoupled[:, j] = (ds_combine(mp, e0) - ds_combine(mm, e0)) / (2*eps)

rho_single = np.max(np.abs(np.linalg.eigvals(J_uncoupled)))
print(f"\nSingle-site rho = {rho_single:.15f}")
print(f"Single-site Delta = {-np.log(rho_single):.15f}")

# Coupled single step: site i affected by site j through Cartan entry a_ij
def coupled_output(m_i, m_j, e_base, g_val, a_ij):
    e_i = e_base.copy() + g_val * a_ij * (m_j - 0.25)
    e_i = np.abs(e_i)
    e_i = np.maximum(e_i, 1e-10)
    e_i /= np.sum(e_i)
    return ds_combine(m_i, e_i)

# Self-Jacobian in coupled context (with one neighbor at a_ij = -1)
J_self = np.zeros((4, 4))
a_ij = -1.0
for j in range(4):
    mp = m0.copy(); mp[j] += eps
    mm = m0.copy(); mm[j] -= eps
    J_self[:, j] = (coupled_output(mp, m0, e0, g, a_ij) -
                     coupled_output(mm, m0, e0, g, a_ij)) / (2*eps)

# Cross-Jacobian: d(output_i)/d(m_j) with a_ij = -1
J_cross = np.zeros((4, 4))
for j in range(4):
    mp = m0.copy(); mp[j] += eps
    mm = m0.copy(); mm[j] -= eps
    J_cross[:, j] = (coupled_output(m0, mp, e0, g, a_ij) -
                      coupled_output(m0, mm, e0, g, a_ij)) / (2*eps)

print(f"\nJ_self eigenvalues:  {sorted(np.abs(np.linalg.eigvals(J_self)), reverse=True)}")
print(f"J_cross eigenvalues: {sorted(np.abs(np.linalg.eigvals(J_cross)), reverse=True)}")

# ============================================================
# FOURIER ANALYSIS: M(k) = J_self + 2*cos(k) * J_cross
# ============================================================
print("\n" + "="*70)
print("FOURIER MODE ANALYSIS: rho(k) = spectral_radius(J_self + 2cos(k) J_cross)")
print("="*70)

k_vals = np.linspace(0, np.pi, 201)
rhos = []
for k in k_vals:
    M = J_self + 2*np.cos(k) * J_cross
    rhos.append(np.max(np.abs(np.linalg.eigvals(M))))
rhos = np.array(rhos)

k_worst = k_vals[np.argmax(rhos)]
rho_max = np.max(rhos)
print(f"Worst mode: k = {k_worst:.6f} (pi = {np.pi:.6f})")
print(f"rho_max (= rho_A_inf) = {rho_max:.15f}")
print(f"Delta_A_inf = {-np.log(rho_max):.15f}")
print(f"|Delta_A_inf - 1/4| = {abs(-np.log(rho_max) - 0.25):.6e}")
print()

# Eigenvalues of M at worst mode
M_worst = J_self + 2*np.cos(k_worst) * J_cross
evals_worst = np.linalg.eigvals(M_worst)
print(f"M(k=pi) eigenvalues: {sorted(np.abs(evals_worst), reverse=True)}")
print()

# ============================================================
# FINITE A_n: exact Fourier modes
# ============================================================
print("="*70)
print("FINITE A_n: rho vs exact Fourier prediction")
print("="*70)

def coupled_step(states, e_base, cartan, g_val):
    r = len(states)
    new = []
    for i in range(r):
        e_i = e_base.copy()
        for j in range(r):
            if i != j and abs(cartan[i,j]) > 0.01:
                e_i = e_i + g_val * cartan[i,j] * (states[j] - 0.25)
        e_i = np.abs(e_i)
        e_i = np.maximum(e_i, 1e-10)
        e_i /= np.sum(e_i)
        new.append(ds_combine(states[i], e_i))
    return new

def find_coupled_eq(cartan, g_val, n_iter=8000):
    r = cartan.shape[0]
    states = [m0.copy() for _ in range(r)]
    for step in range(n_iter):
        new = coupled_step(states, e0, cartan, g_val)
        diff = max(np.max(np.abs(new[i] - states[i])) for i in range(r))
        states = new
        if diff < 1e-14: return states, True
    return states, False

def coupled_rho(states, cartan, g_val, eps=1e-7):
    r = len(states)
    dim = 4 * r
    x0 = np.concatenate(states)
    def F(x):
        ss = [x[4*i:4*(i+1)] for i in range(r)]
        return np.concatenate(coupled_step(ss, e0, cartan, g_val))
    J = np.zeros((dim, dim))
    for j in range(dim):
        xp = x0.copy(); xp[j] += eps
        xm = x0.copy(); xm[j] -= eps
        J[:, j] = (F(xp) - F(xm)) / (2*eps)
    return np.max(np.abs(np.linalg.eigvals(J)))

def cartan_A(n):
    C = np.zeros((n,n))
    for i in range(n): C[i,i] = 2
    for i in range(n-1): C[i,i+1]=-1; C[i+1,i]=-1
    return C

print(f"\n{'n':>4s} {'rank':>4s} {'rho_exact':>14s} {'rho_fourier':>14s} {'diff':>10s} {'Delta':>10s}")
for N in [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50]:
    rank = N - 1
    cartan = cartan_A(rank) if rank > 0 else np.array([[2.0]])
    states, conv = find_coupled_eq(cartan, g)
    rho_exact = coupled_rho(states, cartan, g)

    # Fourier prediction: modes k_m = m*pi/N, m=1,...,N-1
    rho_fourier = 0
    for m in range(1, N):
        k = m * np.pi / N
        M = J_self + 2*np.cos(k) * J_cross
        rho_k = np.max(np.abs(np.linalg.eigvals(M)))
        rho_fourier = max(rho_fourier, rho_k)

    diff = abs(rho_exact - rho_fourier)
    delta = -np.log(rho_exact)
    print(f"{N:4d} {rank:4d} {rho_exact:14.10f} {rho_fourier:14.10f} {diff:10.2e} {delta:10.6f}")

# ============================================================
# D-SERIES: fork analysis
# ============================================================
print("\n" + "="*70)
print("D-SERIES ANALYSIS")
print("="*70)

def cartan_D(n):
    C = np.zeros((n,n))
    for i in range(n): C[i,i]=2
    for i in range(n-2): C[i,i+1]=-1; C[i+1,i]=-1
    C[n-3,n-1]=-1; C[n-1,n-3]=-1
    return C

for n in [4, 5, 6, 8, 10, 15, 20]:
    cartan = cartan_D(n)
    states, conv = find_coupled_eq(cartan, g)
    rho = coupled_rho(states, cartan, g)
    delta = -np.log(rho)
    print(f"D_{n:2d}: rho = {rho:.10f}, Delta = {delta:.6f}")

# ============================================================
# E-SERIES
# ============================================================
print("\n" + "="*70)
print("E-SERIES ANALYSIS")
print("="*70)

def cartan_E(n):
    C = np.zeros((n,n))
    for i in range(n): C[i,i]=2
    links = {
        6: [(0,1),(1,2),(2,3),(3,4),(2,5)],
        7: [(0,1),(1,2),(2,3),(3,4),(4,5),(2,6)],
        8: [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(2,7)],
    }
    for i,j in links[n]: C[i,j]=-1; C[j,i]=-1
    return C

for n in [6, 7, 8]:
    cartan = cartan_E(n)
    states, conv = find_coupled_eq(cartan, g)
    rho = coupled_rho(states, cartan, g)
    delta = -np.log(rho)

    # Graph spectrum
    evals_C = sorted(np.linalg.eigvals(cartan).real)
    adj = np.eye(n)*2 - cartan  # adjacency matrix
    evals_adj = sorted(np.linalg.eigvals(adj).real, reverse=True)

    print(f"E_{n}: rho = {rho:.10f}, Delta = {delta:.6f}")
    print(f"      Adjacency spectrum: {[f'{e:.4f}' for e in evals_adj]}")

# ============================================================
# SATURATION VALUES
# ============================================================
print("\n" + "="*70)
print("SATURATION VALUES AND ANALYTICAL CANDIDATES")
print("="*70)

# A-series limit
M_pi = J_self - 2 * J_cross  # k = pi
rho_A = np.max(np.abs(np.linalg.eigvals(M_pi)))
delta_A = -np.log(rho_A)
print(f"\nA_inf: rho = {rho_A:.15f}")
print(f"       Delta = {delta_A:.15f}")
print(f"       1/4 = {0.25:.15f}")
print(f"       |Delta - 1/4| = {abs(delta_A - 0.25):.2e}")

# D-series limit (compute D_30)
cartan_d30 = cartan_D(30)
states_d30, _ = find_coupled_eq(cartan_d30, g)
rho_D = coupled_rho(states_d30, cartan_d30, g)
delta_D = -np.log(rho_D)
print(f"\nD_30:  rho = {rho_D:.15f}")
print(f"       Delta = {delta_D:.15f}")

# Analytical candidates for D_inf Delta
candidates = {
    "K*/2 = 7/60": 7/60,
    "1/9": 1/9,
    "1/8": 1/8,
    "2/17": 2/17,
    "ln(9/8)": np.log(9/8),
    "7/60": 7/60,
    "1/(2*H+1+H)": 1/7,
}
print(f"\n  Analytical candidates for D_inf Delta = {delta_D:.8f}:")
for name, val in candidates.items():
    print(f"    {name:20s} = {val:.8f}, diff = {abs(delta_D - val):.2e}")
