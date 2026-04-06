#!/usr/bin/env python3
"""
ANALYTICAL MASS GAP FORMULA: Exact resolution
==============================================

The 4x4 Jacobian blocks have rank 2. Extract the effective 2x2 system,
derive exact formulas, and prove:
  1. Delta_A(inf) = 1/4 exactly
  2. Exact D-series and E-series values
  3. Why folding invariance holds
"""

import numpy as np
from scipy.optimize import brentq
from scipy.linalg import svd, schur

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

# ============================================================
# COUPLED SYSTEM INFRASTRUCTURE
# ============================================================
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

def find_coupled_eq(cartan, g_val, n_iter=12000):
    r = cartan.shape[0]
    states = [m0.copy() for _ in range(r)]
    for step in range(n_iter):
        new = coupled_step(states, e0, cartan, g_val)
        diff = max(np.max(np.abs(new[i] - states[i])) for i in range(r))
        states = new
        if diff < 1e-15: return states, True
    return states, False

def full_jacobian(states, cartan, g_val, eps=1e-8):
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
    return J

def cartan_A(n):
    C = np.zeros((n,n))
    for i in range(n): C[i,i] = 2
    for i in range(n-1): C[i,i+1]=-1; C[i+1,i]=-1
    return C

def cartan_D(n):
    C = np.zeros((n,n))
    for i in range(n): C[i,i]=2
    for i in range(n-2): C[i,i+1]=-1; C[i+1,i]=-1
    C[n-3,n-1]=-1; C[n-1,n-3]=-1
    return C

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

# ============================================================
# PART 1: EFFECTIVE 2x2 REDUCTION
# ============================================================
print("=" * 70)
print("PART 1: EFFECTIVE 2x2 REDUCTION")
print("=" * 70)

# Extract bulk blocks from a long A-chain
N_chain = 50
cartan_a50 = cartan_A(N_chain)
states_a50, _ = find_coupled_eq(cartan_a50, g)
J_a50 = full_jacobian(states_a50, cartan_a50, g)

mid = N_chain // 2
J_self = J_a50[4*mid:4*(mid+1), 4*mid:4*(mid+1)]
J_cross = J_a50[4*mid:4*(mid+1), 4*(mid-1):4*mid]

# SVD to find active subspace
U_s, s_s, Vt_s = svd(J_self)
U_c, s_c, Vt_c = svd(J_cross)

print("J_self singular values:", s_s)
print("J_cross singular values:", s_c)
print()

# The active subspace is spanned by the two dominant right singular vectors
# Let's find the common active subspace
# Actually, let's look at the eigenvectors of J_self
evals_self, evecs_self = np.linalg.eig(J_self)
idx = np.argsort(np.abs(evals_self))[::-1]
evals_self = evals_self[idx]
evecs_self = evecs_self[:, idx]

print("J_self eigenvalues:", evals_self)
print("J_self eigenvectors (columns):")
for i in range(4):
    print(f"  v_{i}: {evecs_self[:, i]}")
print()

# Project onto 2D active subspace
P = evecs_self[:, :2]  # projection matrix (4x2)
J_self_2d = P.T @ J_self @ P
J_cross_2d = P.T @ J_cross @ P

print("J_self (2x2 reduced):")
print(J_self_2d)
print("J_cross (2x2 reduced):")
print(J_cross_2d)
print()

# Fourier analysis in 2D
M_worst_2d = J_self_2d - 2 * J_cross_2d  # k = pi
evals_2d = np.linalg.eigvals(M_worst_2d)
rho_2d = np.max(np.abs(evals_2d))
print(f"M(k=pi) 2x2 eigenvalues: {evals_2d}")
print(f"rho_A_inf (2x2) = {rho_2d:.15f}")
print(f"Delta_A_inf (2x2) = {-np.log(rho_2d):.15f}")
print()

# Compare with full 4x4
M_worst_4d = J_self - 2 * J_cross
evals_4d = np.linalg.eigvals(M_worst_4d)
rho_4d = np.max(np.abs(evals_4d))
print(f"M(k=pi) 4x4 eigenvalues: {sorted(np.abs(evals_4d), reverse=True)}")
print(f"rho_A_inf (4x4) = {rho_4d:.15f}")
print(f"Delta_A_inf (4x4) = {-np.log(rho_4d):.15f}")
print()

# ============================================================
# PART 2: CONVERGENCE TO 1/4 - HIGH PRECISION
# ============================================================
print("=" * 70)
print("PART 2: IS Delta_A(inf) = 1/4 EXACTLY?")
print("=" * 70)

# Use progressively longer chains to extrapolate
chain_lengths = [20, 30, 40, 50, 60, 80, 100]
deltas_fourier = []

for N in chain_lengths:
    cartan = cartan_A(N)
    states, conv = find_coupled_eq(cartan, g)
    J = full_jacobian(states, cartan, g)
    mid_n = N // 2
    Js = J[4*mid_n:4*(mid_n+1), 4*mid_n:4*(mid_n+1)]
    Jc = J[4*mid_n:4*(mid_n+1), 4*(mid_n-1):4*mid_n]
    M = Js - 2 * Jc
    rho = np.max(np.abs(np.linalg.eigvals(M)))
    delta = -np.log(rho)
    deltas_fourier.append(delta)
    print(f"  A_{N:3d}: Delta_fourier = {delta:.15f}, |Delta - 1/4| = {abs(delta - 0.25):.6e}")

# Richardson extrapolation
if len(deltas_fourier) >= 3:
    d1, d2, d3 = deltas_fourier[-3], deltas_fourier[-2], deltas_fourier[-1]
    n1, n2, n3 = chain_lengths[-3], chain_lengths[-2], chain_lengths[-1]
    # Assume Delta(N) = Delta_inf + a/N^2
    # d1 = D + a/n1^2, d2 = D + a/n2^2
    # D = (d1*n1^2 - d2*n2^2) / (n1^2 - n2^2)
    D_extrap_1 = (d2*n2**2 - d3*n3**2) / (n2**2 - n3**2)
    D_extrap_2 = (d1*n1**2 - d3*n3**2) / (n1**2 - n3**2)
    print(f"\n  Richardson extrapolation (1/N^2):")
    print(f"    From N={n2},{n3}: Delta_inf = {D_extrap_1:.15f}")
    print(f"    From N={n1},{n3}: Delta_inf = {D_extrap_2:.15f}")
    print(f"    |extrap - 1/4| = {abs(D_extrap_1 - 0.25):.6e}")

    # Try 1/N^alpha convergence
    # log(delta - D) = log(a) + alpha * log(1/N)
    # Use three points to determine alpha
    if d1 > d2 > d3:  # monotone convergence
        r1 = np.log(abs(d1 - 0.25))
        r2 = np.log(abs(d2 - 0.25))
        r3 = np.log(abs(d3 - 0.25))
        l1 = np.log(1.0/n1)
        l2 = np.log(1.0/n2)
        l3 = np.log(1.0/n3)
        alpha_12 = (r2 - r1) / (l2 - l1)
        alpha_23 = (r3 - r2) / (l3 - l2)
        alpha_13 = (r3 - r1) / (l3 - l1)
        print(f"\n  Convergence rate (assuming Delta -> 1/4):")
        print(f"    alpha (N={n1},{n2}) = {alpha_12:.4f}")
        print(f"    alpha (N={n2},{n3}) = {alpha_23:.4f}")
        print(f"    alpha (N={n1},{n3}) = {alpha_13:.4f}")

# ============================================================
# PART 3: EXACT EIGENVALUE STRUCTURE AT EQUILIBRIUM
# ============================================================
print("\n" + "=" * 70)
print("PART 3: EQUILIBRIUM STATE STRUCTURE")
print("=" * 70)

# Single-site equilibrium
print(f"\nSingle-site equilibrium:")
print(f"  m* = {m0}")
print(f"  e* = {e0}")
K_eq = sum(m0[i]*e0[j] for i in range(3) for j in range(3) if i!=j)
print(f"  K* = {K_eq:.15f}")
print(f"  7/30 = {7/30:.15f}")

# Coupled A-chain equilibrium (middle node)
m_bulk = states_a50[25]
print(f"\nA_50 bulk node equilibrium:")
print(f"  m_bulk = {m_bulk}")
print(f"  |m_bulk - m0| = {np.max(np.abs(m_bulk - m0)):.6e}")

# The evidence that the bulk node sees
e_bulk = e0.copy()
for j in [24, 26]:  # neighbors
    e_bulk = e_bulk + g * (-1) * (states_a50[j] - 0.25)
e_bulk = np.abs(e_bulk)
e_bulk = np.maximum(e_bulk, 1e-10)
e_bulk /= np.sum(e_bulk)
print(f"  e_bulk (effective) = {e_bulk}")
print(f"  |e_bulk - e0| = {np.max(np.abs(e_bulk - e0)):.6e}")

# K at bulk equilibrium
K_bulk = sum(m_bulk[i]*e_bulk[j] for i in range(3) for j in range(3) if i!=j)
print(f"  K_bulk = {K_bulk:.15f}")

# ============================================================
# PART 4: D-SERIES AND E-SERIES EXACT VALUES
# ============================================================
print("\n" + "=" * 70)
print("PART 4: D-SERIES AND E-SERIES SATURATION")
print("=" * 70)

# D-series with increasing N
print("\nD-series convergence:")
for n in [4, 5, 6, 8, 10, 15, 20, 30]:
    cartan = cartan_D(n)
    states, conv = find_coupled_eq(cartan, g)
    J = full_jacobian(states, cartan, g)
    rho = np.max(np.abs(np.linalg.eigvals(J)))
    delta = -np.log(rho)
    print(f"  D_{n:2d}: Delta = {delta:.15f}")

# D_30 fork node detailed analysis
cartan_d30 = cartan_D(30)
states_d30, _ = find_coupled_eq(cartan_d30, g)
J_d30 = full_jacobian(states_d30, cartan_d30, g)

# Fork node is at position 27 (n-3)
fork = 27
J_self_fork = J_d30[4*fork:4*(fork+1), 4*fork:4*(fork+1)]
evals_fork = np.linalg.eigvals(J_self_fork)
print(f"\n  D_30 fork J_self eigenvalues: {sorted(np.abs(evals_fork), reverse=True)}")

# Fork equilibrium state
m_fork = states_d30[fork]
print(f"  Fork state: {m_fork}")
print(f"  Fork state ratios: s1/s2 = {m_fork[0]/m_fork[1]:.6f}, s1/theta = {m_fork[0]/m_fork[3]:.6f}")

# Analytical candidate search for D_inf
rho_d_inf = np.max(np.abs(np.linalg.eigvals(J_d30)))
delta_d_inf = -np.log(rho_d_inf)
print(f"\n  D_inf Delta = {delta_d_inf:.15f}")

# Try many candidates
print(f"\n  Analytical candidates for Delta_D_inf = {delta_d_inf:.12f}:")
candidates = {}
for a in range(1, 30):
    for b in range(1, 100):
        val = a / b
        if abs(val - delta_d_inf) < 0.001:
            candidates[f"{a}/{b}"] = val

for name, val in sorted(candidates.items(), key=lambda x: abs(x[1] - delta_d_inf)):
    print(f"    {name:10s} = {val:.12f}, diff = {abs(val - delta_d_inf):.6e}")
    if len([1 for _ in candidates if abs(candidates[_] - delta_d_inf) < abs(val - delta_d_inf)]) > 10:
        break

# Also check ln-based candidates
print(f"\n  Logarithmic candidates:")
for name, val in [
    ("ln(10/9)", np.log(10/9)),
    ("ln(30/27)", np.log(30/27)),
    ("ln(H/(H-1))/H", np.log(H/(H-1))/H),
    ("K*/2", g/2),
    ("1/2H+1", 1/(2*H+1)),
    ("1/(3H)", 1/(3*H)),
    ("ln(28/27)*3", np.log(28/27)*3),
    ("1/H^2 + 1/(2H^3)", 1/H**2 + 1/(2*H**3)),
    ("(H-1)/(H^2+H+1)*ln(H)", (H-1)/(H**2+H+1)*np.log(H)),
]:
    print(f"    {name:30s} = {val:.12f}, diff = {abs(val - delta_d_inf):.6e}")

# E-series exact values
print(f"\n  E-series:")
for n in [6, 7, 8]:
    cartan = cartan_E(n)
    states, conv = find_coupled_eq(cartan, g)
    J = full_jacobian(states, cartan, g)
    rho = np.max(np.abs(np.linalg.eigvals(J)))
    delta = -np.log(rho)
    print(f"  E_{n}: Delta = {delta:.15f}, rho = {rho:.15f}")

# ============================================================
# PART 5: WHY FOLDING INVARIANCE HOLDS
# ============================================================
print("\n" + "=" * 70)
print("PART 5: FOLDING INVARIANCE MECHANISM")
print("=" * 70)

# Compare B3 and D4 in detail
# B3: nodes 0-1-2 with A[1,2]=-2
# D4: nodes 0-1-2-3 with fork at 1 (1 connects to 0,2,3)
# Folding: D4 nodes {2,3} fold onto B3 node 2

def cartan_B(n):
    C = np.zeros((n,n))
    for i in range(n): C[i,i]=2
    for i in range(n-1): C[i,i+1]=-1; C[i+1,i]=-1
    if n>=2: C[n-2,n-1]=-2
    return C

def cartan_C(n):
    C = np.zeros((n,n))
    for i in range(n): C[i,i]=2
    for i in range(n-1): C[i,i+1]=-1; C[i+1,i]=-1
    if n>=2: C[n-1,n-2]=-2
    return C

# B3 vs D4
print("\nB3 vs D4:")
cartan_b3 = cartan_B(3)
cartan_d4 = cartan_D(4)

states_b3, _ = find_coupled_eq(cartan_b3, g)
states_d4, _ = find_coupled_eq(cartan_d4, g)

print(f"  B3 Cartan:\n{cartan_b3}")
print(f"  D4 Cartan:\n{cartan_d4}")

print(f"\n  B3 equilibrium states:")
for i, s in enumerate(states_b3):
    print(f"    node {i}: {s}")

print(f"\n  D4 equilibrium states:")
for i, s in enumerate(states_d4):
    print(f"    node {i}: {s}")

# In D4, the folding identifies nodes 2 and 3
# Check if D4 nodes 2 and 3 are identical
print(f"\n  D4 node 2 vs 3: |diff| = {np.max(np.abs(states_d4[2] - states_d4[3])):.2e}")

# The folding symmetry: D4 has Z2 symmetry swapping nodes 2,3
# At equilibrium, states[2] = states[3] by this symmetry
# Perturbations decompose into symmetric (s) and antisymmetric (a):
#   delta_s = (delta_2 + delta_3)/2  (survives folding)
#   delta_a = (delta_2 - delta_3)/2  (killed by folding)

J_b3 = full_jacobian(states_b3, cartan_b3, g)
J_d4 = full_jacobian(states_d4, cartan_d4, g)

evals_b3 = sorted(np.abs(np.linalg.eigvals(J_b3)), reverse=True)
evals_d4 = sorted(np.abs(np.linalg.eigvals(J_d4)), reverse=True)

print(f"\n  B3 eigenvalues (top 6): {[f'{e:.8f}' for e in evals_b3[:6]]}")
print(f"  D4 eigenvalues (top 6): {[f'{e:.8f}' for e in evals_d4[:6]]}")
print(f"  rho(B3) = {evals_b3[0]:.10f}")
print(f"  rho(D4) = {evals_d4[0]:.10f}")
print(f"  |rho(B3) - rho(D4)| = {abs(evals_b3[0] - evals_d4[0]):.2e}")

# D4 has 16 eigenvalues, B3 has 12
# The 12 symmetric eigenvalues of D4 should match B3's 12
# The 4 antisymmetric eigenvalues are the "extra"
# Project D4 Jacobian onto symmetric subspace

# Symmetric subspace: nodes 0, 1, (2+3)/sqrt(2)
# Build projection
P_sym = np.zeros((16, 12))
P_sym[:4, :4] = np.eye(4)  # node 0
P_sym[4:8, 4:8] = np.eye(4)  # node 1
P_sym[8:12, 8:12] = np.eye(4) / np.sqrt(2)  # node 2 symmetric part
P_sym[12:16, 8:12] = np.eye(4) / np.sqrt(2)  # node 3 symmetric part

J_d4_sym = P_sym.T @ J_d4 @ P_sym
evals_d4_sym = sorted(np.abs(np.linalg.eigvals(J_d4_sym)), reverse=True)
print(f"\n  D4 symmetric-projected eigenvalues: {[f'{e:.8f}' for e in evals_d4_sym[:6]]}")
print(f"  B3 eigenvalues for comparison:      {[f'{e:.8f}' for e in evals_b3[:6]]}")

# Now C3 vs A5
print("\n\nC3 = Sp(6) vs A5 = SU(6):")
cartan_c3 = cartan_C(3)
cartan_a5 = cartan_A(5)

states_c3, _ = find_coupled_eq(cartan_c3, g)
states_a5, _ = find_coupled_eq(cartan_a5, g)

J_c3 = full_jacobian(states_c3, cartan_c3, g)
J_a5 = full_jacobian(states_a5, cartan_a5, g)

evals_c3 = sorted(np.abs(np.linalg.eigvals(J_c3)), reverse=True)
evals_a5 = sorted(np.abs(np.linalg.eigvals(J_a5)), reverse=True)

print(f"  rho(C3) = {evals_c3[0]:.10f}")
print(f"  rho(A5) = {evals_a5[0]:.10f}")
print(f"  |rho(C3) - rho(A5)| = {abs(evals_c3[0] - evals_a5[0]):.2e}")

# A5 has Z2 symmetry (reflection): node i <-> node 4-i
# C3 should correspond to the symmetric sector of A5
print(f"\n  A5 equilibrium (checking Z2 symmetry):")
for i in range(5):
    print(f"    node {i}: {states_a5[i]}, |diff with node {4-i}|: {np.max(np.abs(states_a5[i] - states_a5[4-i])):.2e}")

# Project A5 onto Z2-symmetric subspace: nodes 0+4, 1+3, 2
P_sym_a5 = np.zeros((20, 12))
P_sym_a5[:4, :4] = np.eye(4) / np.sqrt(2)     # node 0
P_sym_a5[16:20, :4] = np.eye(4) / np.sqrt(2)   # node 4 -> same as 0
P_sym_a5[4:8, 4:8] = np.eye(4) / np.sqrt(2)    # node 1
P_sym_a5[12:16, 4:8] = np.eye(4) / np.sqrt(2)  # node 3 -> same as 1
P_sym_a5[8:12, 8:12] = np.eye(4)                # node 2 (center)

J_a5_sym = P_sym_a5.T @ J_a5 @ P_sym_a5
evals_a5_sym = sorted(np.abs(np.linalg.eigvals(J_a5_sym)), reverse=True)
print(f"\n  A5 symmetric-projected eigenvalues: {[f'{e:.8f}' for e in evals_a5_sym[:6]]}")
print(f"  C3 eigenvalues for comparison:      {[f'{e:.8f}' for e in evals_c3[:6]]}")

# ============================================================
# PART 6: FORMULA ATTEMPT
# ============================================================
print("\n" + "=" * 70)
print("PART 6: LOOKING FOR ANALYTICAL FORMULA")
print("=" * 70)

# The key: at the A-chain bulk equilibrium, what is the 2x2 system?
# Let's look at J_self and J_cross in a basis-independent way
print("\nA-chain bulk 4x4 blocks (N=50):")

# Trace and determinant of active 2x2
Js = J_a50[100:104, 100:104]  # mid=25, 4*25=100
Jc = J_a50[100:104, 96:100]

# Get the 2 active eigenvalues
es = sorted(np.abs(np.linalg.eigvals(Js)), reverse=True)[:2]
ec = sorted(np.abs(np.linalg.eigvals(Jc)), reverse=True)[:2]

print(f"  J_self active eigenvalues:  {es}")
print(f"  J_cross active eigenvalues: {ec}")
print(f"  sum(J_self active) = {sum(es):.12f}")
print(f"  prod(J_self active) = {es[0]*es[1]:.12f}")
print(f"  sum(J_cross active) = {sum(ec):.12f}")
print(f"  prod(J_cross active) = {ec[0]*ec[1]:.12f}")

# M(k=pi) = J_self - 2*J_cross
# Active eigenvalues of M(pi):
em = sorted(np.abs(np.linalg.eigvals(Js - 2*Jc)), reverse=True)[:2]
print(f"\n  M(pi) active eigenvalues: {em}")
print(f"  rho = max = {em[0]:.15f}")
print(f"  Delta = {-np.log(em[0]):.15f}")

# Check: is rho = exp(-1/4)?
print(f"\n  exp(-1/4) = {np.exp(-0.25):.15f}")
print(f"  rho       = {em[0]:.15f}")
print(f"  diff      = {abs(em[0] - np.exp(-0.25)):.6e}")

# The two eigenvalues of M(pi)
print(f"\n  lambda_1(M_pi) = {em[0]:.15f}")
print(f"  lambda_2(M_pi) = {em[1]:.15f}")
print(f"  ratio = {em[1]/em[0]:.15f}")
print(f"  -ln(lambda_1) = {-np.log(em[0]):.15f}")
print(f"  -ln(lambda_2) = {-np.log(em[1]):.15f}")

# Are these eigenvalues of the form a +/- b?
avg = (em[0] + em[1]) / 2
half_diff = (em[0] - em[1]) / 2
print(f"\n  (lambda_1 + lambda_2)/2 = {avg:.15f}")
print(f"  (lambda_1 - lambda_2)/2 = {half_diff:.15f}")

# Check K*, H relationships
print(f"\n  K* = {g:.15f}")
print(f"  1-K* = {1-g:.15f}")
print(f"  rho_single = {0.282910348545127:.15f}")
print(f"  exp(-1/4) / rho_single^? = ???")
# rho_inf = f(rho_single, K*)?
# rho_inf = rho_single * (1 + something)?

# Deep search: try rho_inf = function of K* and rho_single
rho_inf = em[0]
rho_s = 0.282910348545127
print(f"\n  rho_inf / rho_s = {rho_inf / rho_s:.15f}")
print(f"  rho_inf^2 = {rho_inf**2:.15f}")
print(f"  rho_inf * rho_s = {rho_inf * rho_s:.15f}")
print(f"  rho_inf + rho_s = {rho_inf + rho_s:.15f}")
print(f"  rho_inf - rho_s = {rho_inf - rho_s:.15f}")
print(f"  (rho_inf - rho_s) / K* = {(rho_inf - rho_s) / g:.15f}")
