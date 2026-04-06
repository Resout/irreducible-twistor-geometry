#!/usr/bin/env python3
"""
GAUGE GROUP -> MASS GAP: Analytical Fourier decomposition (v2)

Fix: compute J_self and J_cross at the ACTUAL coupled equilibrium
of a long chain, not at the single-site equilibrium.
"""

import numpy as np
from scipy.optimize import brentq

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

def find_single_eq():
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

m0, e0 = find_single_eq()

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

def coupled_rho_and_jacobian(states, cartan, g_val, eps=1e-7):
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
    evals = np.abs(np.linalg.eigvals(J))
    return np.max(evals), J

def cartan_A(n):
    C = np.zeros((n,n))
    for i in range(n): C[i,i] = 2
    for i in range(n-1): C[i,i+1]=-1; C[i+1,i]=-1
    return C

# ============================================================
# EXTRACT J_self AND J_cross FROM LONG CHAIN
# ============================================================
N_chain = 30  # long enough to extract bulk behavior
cartan = cartan_A(N_chain)
states, conv = find_coupled_eq(cartan, g)
rho_full, J_full = coupled_rho_and_jacobian(states, cartan, g)

# Extract blocks from the middle of the chain
mid = N_chain // 2
J_self_chain = J_full[4*mid:4*(mid+1), 4*mid:4*(mid+1)]
J_left = J_full[4*mid:4*(mid+1), 4*(mid-1):4*mid]
J_right = J_full[4*mid:4*(mid+1), 4*(mid+1):4*(mid+2)]

print("=== A-CHAIN BULK JACOBIAN BLOCKS (extracted from A_30 middle) ===")
print(f"J_self eigenvalues:  {sorted(np.abs(np.linalg.eigvals(J_self_chain)), reverse=True)}")
print(f"J_left eigenvalues:  {sorted(np.abs(np.linalg.eigvals(J_left)), reverse=True)}")
print(f"J_right eigenvalues: {sorted(np.abs(np.linalg.eigvals(J_right)), reverse=True)}")
print(f"J_left ~= J_right? {np.max(np.abs(J_left - J_right)):.2e}")
J_cross_chain = (J_left + J_right) / 2  # should be identical by symmetry

print()

# Fourier analysis with correct blocks
print("=== FOURIER ANALYSIS WITH CORRECT BULK BLOCKS ===")
k_vals = np.linspace(0, np.pi, 201)
rhos = []
for k in k_vals:
    M = J_self_chain + 2*np.cos(k) * J_cross_chain
    rhos.append(np.max(np.abs(np.linalg.eigvals(M))))
rhos = np.array(rhos)

k_worst = k_vals[np.argmax(rhos)]
rho_fourier_inf = np.max(rhos)
delta_fourier_inf = -np.log(rho_fourier_inf)

print(f"Worst mode: k = {k_worst:.6f}")
print(f"rho_A_inf (Fourier) = {rho_fourier_inf:.15f}")
print(f"Delta_A_inf = {delta_fourier_inf:.15f}")
print(f"Actual A_30 rho     = {rho_full:.15f}")
print(f"Match: {abs(rho_fourier_inf - rho_full):.2e}")
print()

# Check Delta_A_inf = 1/4?
print(f"|Delta_A_inf - 1/4| = {abs(delta_fourier_inf - 0.25):.6e}")

# Finite-n Fourier prediction
print(f"\n{'n':>4s} {'rank':>4s} {'rho_exact':>14s} {'rho_fourier':>14s} {'diff':>10s} {'Delta':>10s}")
for N in [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50]:
    rank = N - 1
    if rank == 0:
        rho_exact = np.max(np.abs(np.linalg.eigvals(J_self_chain)))
        rho_fourier = rho_exact
    else:
        c = cartan_A(rank)
        st, cv = find_coupled_eq(c, g)
        rho_exact = coupled_rho_and_jacobian(st, c, g)[0]
        # Fourier modes: k_m = m*pi/N, m=1,...,N-1
        rho_fourier = 0
        for m in range(1, N):
            k = m * np.pi / N
            M = J_self_chain + 2*np.cos(k) * J_cross_chain
            rho_k = np.max(np.abs(np.linalg.eigvals(M)))
            rho_fourier = max(rho_fourier, rho_k)

    diff = abs(rho_exact - rho_fourier)
    delta = -np.log(rho_exact)
    print(f"{N:4d} {rank:4d} {rho_exact:14.10f} {rho_fourier:14.10f} {diff:10.2e} {delta:10.6f}")

# ============================================================
# D-SERIES: extract fork Jacobian from D_20
# ============================================================
print("\n" + "="*70)
print("D-SERIES: FORK NODE ANALYSIS")
print("="*70)

def cartan_D(n):
    C = np.zeros((n,n))
    for i in range(n): C[i,i]=2
    for i in range(n-2): C[i,i+1]=-1; C[i+1,i]=-1
    C[n-3,n-1]=-1; C[n-1,n-3]=-1
    return C

cartan_d = cartan_D(20)
states_d, _ = find_coupled_eq(cartan_d, g)
rho_d, J_d = coupled_rho_and_jacobian(states_d, cartan_d, g)

# Fork node is n-3 = 17, connected to 16, 18, 19
fork = 17
J_self_fork = J_d[4*fork:4*(fork+1), 4*fork:4*(fork+1)]
J_cross_fork_16 = J_d[4*fork:4*(fork+1), 4*16:4*17]  # chain neighbor
J_cross_fork_18 = J_d[4*fork:4*(fork+1), 4*18:4*19]  # branch 1
J_cross_fork_19 = J_d[4*fork:4*(fork+1), 4*19:4*20]  # branch 2

print(f"D_20 rho = {rho_d:.10f}, Delta = {-np.log(rho_d):.6f}")
print(f"Fork node {fork}: connected to {fork-1}, {fork+1}, {19}")
print(f"J_self_fork eigenvalues: {sorted(np.abs(np.linalg.eigvals(J_self_fork)), reverse=True)}")
print(f"J_cross_16 eigenvalues:  {sorted(np.abs(np.linalg.eigvals(J_cross_fork_16)), reverse=True)}")
print(f"J_cross_18 eigenvalues:  {sorted(np.abs(np.linalg.eigvals(J_cross_fork_18)), reverse=True)}")
print(f"J_cross_19 eigenvalues:  {sorted(np.abs(np.linalg.eigvals(J_cross_fork_19)), reverse=True)}")

# Check if all cross blocks are the same
print(f"J_cross_18 ~= J_cross_19? {np.max(np.abs(J_cross_fork_18 - J_cross_fork_19)):.2e}")
print(f"J_cross_16 ~= J_cross_18? {np.max(np.abs(J_cross_fork_16 - J_cross_fork_18)):.2e}")

# Also extract a bulk chain node from D_20
chain_mid = 8  # middle of the chain part
J_self_d_chain = J_d[4*chain_mid:4*(chain_mid+1), 4*chain_mid:4*(chain_mid+1)]
J_cross_d_left = J_d[4*chain_mid:4*(chain_mid+1), 4*(chain_mid-1):4*chain_mid]
print(f"\nBulk chain node {chain_mid}:")
print(f"J_self eigenvalues: {sorted(np.abs(np.linalg.eigvals(J_self_d_chain)), reverse=True)}")
print(f"Same as A-chain J_self? {np.max(np.abs(J_self_d_chain - J_self_chain)):.2e}")

# ============================================================
# KEY QUESTION: Does the equilibrium state differ from m0?
# ============================================================
print("\n" + "="*70)
print("EQUILIBRIUM STATE COMPARISON")
print("="*70)

# A-chain
for i in [0, mid, N_chain-1]:
    diff = np.max(np.abs(states[i] - m0))
    print(f"A_30 node {i:2d}: |m - m0| = {diff:.2e}")

# D-chain
for i in [0, chain_mid, fork, 18, 19]:
    diff = np.max(np.abs(states_d[i] - m0))
    print(f"D_20 node {i:2d}: |m - m0| = {diff:.2e}")

# ============================================================
# FULL EIGENVALUE SPECTRUM OF COUPLED SYSTEM
# ============================================================
print("\n" + "="*70)
print("EIGENVALUE SPECTRUM: A_30 (top 20)")
print("="*70)
evals_A30 = sorted(np.abs(np.linalg.eigvals(J_full)), reverse=True)
for i, ev in enumerate(evals_A30[:20]):
    print(f"  lambda_{i+1:2d} = {ev:.10f}  (Delta = {-np.log(ev):.6f})")

print(f"\n  rho = {evals_A30[0]:.10f}")
print(f"  Gap ratio lambda_2/lambda_1 = {evals_A30[1]/evals_A30[0]:.6f}")
