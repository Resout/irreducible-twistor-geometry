"""
Popov Second Variation in DS Variables

Tests whether the DS transfer operator at K*=7/30 equilibrium
is the exponential of the Popov Hessian restricted to the DS fibre.

CLAIM: exp(-H_Popov) = J_DS
where H_Popov = K + ad(a*) is the J-holomorphic Chern-Simons Hessian.

If true: the DS QFT IS quantised Yang-Mills (the identification theorem).

9 Stages:
  1. Equilibrium m*, e* at K*=7/30
  2. Pauli embedding M*, M*^{-1}
  3. 3x3 DS Jacobian (analytical), eigendecomposition
  4. Wirtinger anti-holomorphic Jacobian J_anti
  5. Equilibrium (0,1)-form a* = M*^{-1} · Pauli(J_anti)
  6. su(2) decomposition of a*, adjoint commutator ad(a*)
  7. Popov Hessian decomposition: H = -ln(J) vs K + ad(a*)
  8. Eigenvalue splitting test
  9. Full reconstruction: exp(-(K + ad(a*))) vs J
"""

import numpy as np
from scipy.linalg import expm, logm
np.set_printoptions(precision=10, linewidth=100)

# === Pauli matrices ===
sig1 = np.array([[0, 1], [1, 0]], dtype=complex)
sig2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
sig3 = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
pauli = [sig1, sig2, sig3]
basis = [sig1, sig2, sig3, I2]

def to_M(m):
    return (m[3]*I2 + m[0]*sig1 + m[1]*sig2 + m[2]*sig3) / np.sqrt(2)

def ds_combine(m, e):
    s, th = m[:3], m[3]
    es, ph = e[:3], e[3]
    sn = np.array([s[i]*es[i] + s[i]*ph + th*es[i] for i in range(3)])
    tn = th * ph
    K = sum(s[i]*es[j] for i in range(3) for j in range(3) if i != j)
    d = 1.0 - K
    return np.array(list(sn/d) + [tn/d]), K

def enforce_floor(m):
    s, th = m[:3], m[3]
    ssq = sum(si**2 for si in s)
    born = th**2 / (ssq + th**2) if (ssq + th**2) > 0 else 1
    if born >= 1.0/27:
        return m.copy()
    ss = sum(s)
    r = ssq / ss**2
    disc = (2*r)**2 + 4*(26-r)*r
    t = (-2*r + np.sqrt(disc)) / (2*(26-r))
    alpha = (1 - t) / ss
    return np.array([s[0]*alpha, s[1]*alpha, s[2]*alpha, t])

def full_step(m, e):
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds), K

def enforce_floor_complex(m):
    """Floor for complex masses (small imaginary perturbations)."""
    s, th = m[:3], m[3]
    sq_abs = sum(abs(si)**2 for si in s)
    abs_th_sq = abs(th)**2
    total = sq_abs + abs_th_sq
    born = abs_th_sq / total if total > 1e-300 else 1
    if born >= 1.0/27 - 1e-14:
        return m.copy()
    ss = sum(s)
    r_abs = sq_abs / abs(ss)**2 if abs(ss) > 1e-300 else 26.0
    t_abs = 1.0 / (np.sqrt(26.0/r_abs) + 1)
    th_new = (th / abs(th)) * t_abs if abs(th) > 1e-300 else t_abs
    alpha = (1 - th_new) / ss if abs(ss) > 1e-300 else 0
    return np.array([s[0]*alpha, s[1]*alpha, s[2]*alpha, th_new], dtype=complex)

def ds_combine_complex(m, e):
    s, th = m[:3], m[3]
    es, ph = e[:3], e[3]
    sn = np.array([s[i]*es[i] + s[i]*ph + th*es[i] for i in range(3)])
    tn = th * ph
    K = sum(s[i]*es[j] for i in range(3) for j in range(3) if i != j)
    d = 1.0 - K
    return np.array(list(sn/d) + [tn/d], dtype=complex), K

def full_step_complex(m, e):
    m_ds, K = ds_combine_complex(m, e)
    return enforce_floor_complex(m_ds), K

def make_evidence(p_dom):
    floor = 1.0/27
    p_weak = (1 - p_dom)/2
    sc = 1 - floor
    raw = np.array([np.sqrt(p_dom*sc), np.sqrt(p_weak*sc),
                     np.sqrt(p_weak*sc), np.sqrt(floor)])
    return raw / raw.sum()


# =====================================================================
print("=" * 70)
print("POPOV SECOND VARIATION IN DS VARIABLES")
print("=" * 70)

# === Stage 1: Equilibrium ===
print("\nStage 1: Finding K*=7/30 equilibrium (secant method)...")
TARGET_K = 7.0/30

def K_at_pdom(p):
    e = make_evidence(p)
    m = np.array([0.4, 0.15, 0.15, 0.3])
    for _ in range(10000):
        m, _ = full_step(m, e)
    _, K = ds_combine(m, e)
    return K

# Secant method starting near known solution (p_dom ~ 0.932)
p0, p1 = 0.930, 0.935
K0, K1 = K_at_pdom(p0), K_at_pdom(p1)
for step in range(50):
    f0, f1 = K0 - TARGET_K, K1 - TARGET_K
    if abs(f1 - f0) < 1e-30:
        break
    p2 = p1 - f1 * (p1 - p0) / (f1 - f0)
    p2 = max(0.5, min(0.999, p2))  # safety clamp
    p0, K0 = p1, K1
    p1 = p2
    K1 = K_at_pdom(p1)
    if abs(p1 - p0) < 1e-14:
        break

p_dom = p1
e_star = make_evidence(p_dom)
m = np.array([0.4, 0.15, 0.15, 0.3])
for _ in range(10000):
    m, _ = full_step(m, e_star)
m_star = m.copy()
_, K_star = ds_combine(m_star, e_star)

print(f"  m* = [{', '.join(f'{x:.10f}' for x in m_star)}]")
print(f"  e* = [{', '.join(f'{x:.10f}' for x in e_star)}]")
print(f"  K* = {K_star:.15f}  (target: {TARGET_K:.15f})")
print(f"  |K* - 7/30| = {abs(K_star - TARGET_K):.2e}")
print(f"  p_dom = {p_dom:.15f}")


# === Stage 2: Pauli embedding ===
print("\nStage 2: Pauli embedding...")
M_star = to_M(m_star)
M_star_inv = np.linalg.inv(M_star)
det_M = np.linalg.det(M_star)
print(f"  det(M*) = {det_M:.10f}")

# DS output (before floor)
m_ds, _ = ds_combine(m_star, e_star)
M_ds = to_M(m_ds)
M_ds_inv = np.linalg.inv(M_ds)


# === Stage 3: 3x3 DS Jacobian (analytical) ===
print("\nStage 3: DS Jacobian (analytical, projected to L₁=1)...")

s1v, s2v, s3v, thv = m_star
e1v, e2v, e3v, phv = e_star
oneMinusK = 1.0 - K_star

S_vals = [s1v*(e1v+phv)+thv*e1v, s2v*(e2v+phv)+thv*e2v,
          s3v*(e3v+phv)+thv*e3v, thv*phv]
N_vals = [S_vals[i] / oneMinusK for i in range(4)]

# DS Jacobian (4x4)
dSdm = np.zeros((4, 4))
for i in range(3):
    dSdm[i, i] = e_star[i] + phv
    dSdm[i, 3] = e_star[i]
dSdm[3, 3] = phv
dKdm = np.array([e2v+e3v, e1v+e3v, e1v+e2v, 0.0])

J_DS = np.zeros((4, 4))
for i in range(4):
    for j in range(4):
        J_DS[i, j] = (dSdm[i, j] * oneMinusK + S_vals[i] * dKdm[j]) / oneMinusK**2

# Floor Jacobian
Ss = N_vals[0] + N_vals[1] + N_vals[2]
Sq = N_vals[0]**2 + N_vals[1]**2 + N_vals[2]**2
R = Sq / Ss**2
u_fl = np.sqrt(26*R)
t_val = (u_fl - R) / (26 - R)
g = (1 - t_val) / Ss

dtdR = (338 + 13*R - 26*u_fl) / (u_fl * (26 - R)**2)
dRdN = [2*(N_vals[i]*Ss - Sq)/Ss**3 for i in range(3)]
dtdN = [dtdR * dRdN[i] for i in range(3)]
dgdN = [(-dtdN[i]*Ss - (1-t_val))/Ss**2 for i in range(3)]

J_floor = np.zeros((4, 4))
for i in range(3):
    for j in range(3):
        J_floor[i, j] = (g if i == j else 0) + N_vals[i]*dgdN[j]
J_floor[3, 0] = dtdN[0]; J_floor[3, 1] = dtdN[1]; J_floor[3, 2] = dtdN[2]

# Total Jacobian
J_total = J_floor @ J_DS

# Floor correction Jacobian (for Ward connection basis)
J_corr = (J_floor - np.eye(4)) @ J_DS

# Project to L₁=1 tangent space: v_i = e_i - e_4
P_proj = np.zeros((4, 3))
P_proj[0, 0] = P_proj[1, 1] = P_proj[2, 2] = 1.0
P_proj[3, :] = -1.0
J3 = P_proj.T @ J_total @ P_proj

evals_raw, evecs_raw = np.linalg.eig(J3)
idx = np.argsort(-np.abs(evals_raw))
evals = evals_raw[idx].real
evecs = evecs_raw[:, idx].real

print(f"  λ₀ = {evals[0]:.10f}")
print(f"  λ₁ = {evals[1]:.10f}")
print(f"  λ₂ = {evals[2]:.10f}")
print(f"  Δ₀ = -ln(λ₀) = {-np.log(abs(evals[0])):.10f}")
print(f"  Δ₁ = -ln(λ₁) = {-np.log(abs(evals[1])):.10f}")
print(f"  Splitting: (Δ₁-Δ₀)/Δ₀ = {(-np.log(abs(evals[1]))+np.log(abs(evals[0])))/-np.log(abs(evals[0]))*100:.4f}%")

# Numerical cross-check of Jacobian
eps_num = 1e-7
J4_num = np.zeros((4, 4))
for j in range(4):
    mp = m_star.copy(); mp[j] += eps_num
    mm = m_star.copy(); mm[j] -= eps_num
    fp, _ = full_step(mp, e_star)
    fm, _ = full_step(mm, e_star)
    J4_num[:, j] = (fp - fm) / (2*eps_num)

J3_num = P_proj.T @ J4_num @ P_proj
evals_num = sorted(np.linalg.eigvals(J3_num).real, key=lambda x: -abs(x))
print(f"\n  Numerical cross-check:")
print(f"    λ₀(num) = {evals_num[0]:.10f}")
print(f"    λ₁(num) = {evals_num[1]:.10f}")
print(f"    λ₂(num) = {evals_num[2]:.10f}")
print(f"    ||J_analytical - J_numerical|| = {np.linalg.norm(J_total - J4_num):.2e}")


# === Stage 4: Anti-holomorphic Jacobian J_anti (ANALYTICAL) ===
print("\nStage 4: Wirtinger anti-holomorphic Jacobian J_anti (analytical)...")

# The anti-holomorphic floor Jacobian at w = DS(m*, e*) has RANK 1.
# All rows are determined by the theta-row:
#   dtheta_out/dw_bar_j for j = s1,s2,s3,theta
#   ds_i_out/dw_bar_j = -N[i]/Ss * dtheta_out/dw_bar_j
#
# From Wirtinger calculus on floor(w):
#   theta_out = (w4/|w4|) * c,  c = sqrt(sum|w_i|^2/26)
#   dtheta_out/dw_bar_si = w4 * w_i / (52 * |w4| * c)  [paper eq 968]
#   dtheta_out/dw_bar_th = -c * w4^2 / (2|w4|^3)       [paper eq 964]
# At real w (DS output is real): w4=theta_ds>0, simplifies.

# Quantities at DS output point
N_s = N_vals[:3]  # singleton components of DS output
N_th = N_vals[3]   # theta component of DS output
Ss_fl = sum(N_s)
c_fl = np.sqrt(sum(n**2 for n in N_s) / 26)

# Theta-row of anti-holomorphic floor Jacobian
dth_dw_bar = np.zeros(4)
for i in range(3):
    dth_dw_bar[i] = N_s[i] / (52 * c_fl)    # d(theta_out)/d(s_bar_i)
dth_dw_bar[3] = -c_fl / (2 * N_th)           # d(theta_out)/d(theta_bar)

# Full anti-holomorphic floor Jacobian (rank 1)
J_floor_anti = np.zeros((4, 4))
for i in range(3):
    J_floor_anti[i, :] = -N_s[i] / Ss_fl * dth_dw_bar  # singleton rows
J_floor_anti[3, :] = dth_dw_bar  # theta row

# J_anti = J_floor_anti . J_DS  (chain rule: floor's anti-hol part × DS holomorphic)
J_anti = J_floor_anti @ J_DS

sv_anti = np.linalg.svd(J_anti, compute_uv=False)
print(f"  ||J_anti|| = {np.linalg.norm(J_anti):.6f}")
print(f"  Singular values: [{', '.join(f'{v:.6f}' for v in sv_anti)}]")
rank_ratio = sv_anti[1]/sv_anti[0] if sv_anti[0] > 1e-15 else 0
print(f"  σ₂/σ₁ = {rank_ratio:.2e}  (rank {'1' if rank_ratio < 0.01 else '≥2'})")
print(f"  Floor anti-hol Jacobian rank: {np.linalg.matrix_rank(J_floor_anti, tol=1e-10)}")
print(f"  dθ/dw̄ = [{', '.join(f'{v:.8f}' for v in dth_dw_bar)}]")


# === Stage 5: Equilibrium (0,1)-form a* ===
print("\nStage 5: Equilibrium (0,1)-form a*...")

# The (0,1)-form at equilibrium: a* = M_ds^{-1} · sum_j basis[j] · J_anti[j,:] / sqrt(2)
# But J_anti has rank 1, so pick the dominant direction.
# More precisely: for each "direction" mu, form the connection matrix B[mu]

# Using the floor correction approach (as in ward_highprecision):
# B[mu] = sum_j M_ds_inv · basis[j] · J_corr[j,mu] / sqrt(2)
# This gives the connection per spacetime direction mu

# For the Wirtinger version: same structure but with J_anti
# a*[mu] = sum_j M_ds_inv · basis[j] · J_anti[j,mu] / sqrt(2)

a_star_per_dir = []
for mu in range(4):
    a_mu = np.zeros((2, 2), dtype=complex)
    for j in range(4):
        a_mu += M_ds_inv @ basis[j] * (J_anti[j, mu] / np.sqrt(2))
    a_star_per_dir.append(a_mu)

# Combined a* (Frobenius sum over all directions)
a_star_total_sq = sum(np.trace(a.conj().T @ a).real for a in a_star_per_dir)
print(f"  ||a*||² (sum over directions) = {a_star_total_sq:.6f}")

# Since J_anti has rank 1, pick the dominant singular vector as THE direction
U_anti, S_anti, Vh_anti = np.linalg.svd(J_anti)
u_dominant = U_anti[:, 0] * S_anti[0]  # dominant image in C^4

# The connection form from the dominant direction
delta_M = sum(u_dominant[j] * basis[j] for j in range(4)) / np.sqrt(2)
a_star = M_ds_inv @ delta_M

print(f"  a* (from dominant SVD direction):")
print(f"    [{a_star[0,0]:.8f}  {a_star[0,1]:.8f}]")
print(f"    [{a_star[1,0]:.8f}  {a_star[1,1]:.8f}]")
print(f"  ||a*|| = {np.linalg.norm(a_star):.10f}")
print(f"  tr(a*) = {np.trace(a_star):.10f}")


# === Stage 6: su(2) decomposition and adjoint commutator ===
print("\nStage 6: su(2) decomposition of a* and adjoint commutator...")

# Decompose: a* = c₀I + c₁σ₁ + c₂σ₂ + c₃σ₃
c0 = np.trace(a_star) / 2
c = np.array([np.trace(a_star @ pauli[i]) / 2 for i in range(3)])

print(f"  c₀ (scalar)  = {c0:.10f}")
print(f"  c₁ (σ₁)     = {c[0]:.10f}")
print(f"  c₂ (σ₂)     = {c[1]:.10f}")
print(f"  c₃ (σ₃)     = {c[2]:.10f}")
print(f"  ||c||        = {np.linalg.norm(c):.10f}")

scalar_frac = abs(c0)**2 / (abs(c0)**2 + np.sum(np.abs(c)**2))
print(f"  Scalar fraction:  {scalar_frac:.4f}")
print(f"  su(2) fraction:   {1-scalar_frac:.4f}")

# Adjoint commutator: [a*_su2, sigma_j] = 2i sum_k (c x e_j)_k sigma_k
# Matrix ad_{kj} = 2i sum_i eps_{ijk} c_i
eps_tensor = np.zeros((3, 3, 3))
eps_tensor[0, 1, 2] = eps_tensor[1, 2, 0] = eps_tensor[2, 0, 1] = 1
eps_tensor[0, 2, 1] = eps_tensor[2, 1, 0] = eps_tensor[1, 0, 2] = -1

ad_a = np.zeros((3, 3), dtype=complex)
for k in range(3):
    for j in range(3):
        ad_a[k, j] = 2j * sum(c[i] * eps_tensor[i, j, k] for i in range(3))

print(f"\n  ad(a*) = ")
for i in range(3):
    print(f"    [{', '.join(f'{ad_a[i,j]:12.8f}' for j in range(3))}]")

ad_evals = np.linalg.eigvals(ad_a)
ad_evals_sorted = sorted(ad_evals, key=lambda x: -abs(x))
print(f"  Eigenvalues of ad(a*): [{', '.join(f'{v:.8f}' for v in ad_evals_sorted)}]")
print(f"  |max eigenvalue| = {max(abs(ad_evals)):.10f}")


# === Stage 7: Popov Hessian decomposition ===
print("\n" + "=" * 70)
print("Stage 7: Popov Hessian decomposition")
print("=" * 70)

# The DS Hessian (matrix logarithm): H = -ln(J3)
# J3 has eigenvalues {lambda_0, lambda_1, ~0}
# Work in eigenbasis for the 2D active subspace

Delta_0 = -np.log(abs(evals[0]))
Delta_1 = -np.log(abs(evals[1]))

print(f"\n  DS Hessian eigenvalues:")
print(f"    E₀ = -ln(λ₀) = {Delta_0:.10f}")
print(f"    E₁ = -ln(λ₁) = {Delta_1:.10f}")
print(f"    E₂ = -ln(λ₂) → ∞  (conservation law)")

# Full 3x3 Hessian in original basis (regularise λ₂)
lam2_reg = max(abs(evals[2]), 1e-15)
H3 = evecs @ np.diag([-np.log(abs(evals[0])),
                        -np.log(abs(evals[1])),
                        -np.log(lam2_reg)]) @ np.linalg.inv(evecs)

print(f"\n  H₃ (DS Hessian in original basis, regularised):")
for i in range(3):
    print(f"    [{', '.join(f'{H3[i,j]:12.6f}' for j in range(3))}]")

# Project ad(a*) to real part for comparison with real H3
ad_real = ad_a.real
ad_imag = ad_a.imag

print(f"\n  ad(a*) real part norm:  {np.linalg.norm(ad_real):.6f}")
print(f"  ad(a*) imag part norm:  {np.linalg.norm(ad_imag):.6f}")

# Required kinetic: K = H - ad(a*)
# Since H is real, and ad(a*) may be complex:
K3 = H3 - ad_a
print(f"\n  Required K = H - ad(a*):")
for i in range(3):
    print(f"    [{', '.join(f'{K3[i,j]:12.6f}' for j in range(3))}]")

K_evals = np.linalg.eigvals(K3)
print(f"  K eigenvalues: [{', '.join(f'{v:.6f}' for v in sorted(K_evals, key=lambda x: -abs(x)))}]")

# Check: is K Hermitian? symmetric?
K_herm = (K3 + K3.conj().T) / 2
K_aherm = (K3 - K3.conj().T) / 2
print(f"  ||K_hermitian|| = {np.linalg.norm(K_herm):.6f}")
print(f"  ||K_anti-hermitian|| = {np.linalg.norm(K_aherm):.6f}")


# === Stage 8: Eigenvalue splitting test ===
print("\n" + "=" * 70)
print("Stage 8: Eigenvalue splitting test")
print("=" * 70)

splitting = Delta_1 - Delta_0
mean_Delta = (Delta_0 + Delta_1) / 2

print(f"\n  Mean energy:     E_mean = {mean_Delta:.10f}")
print(f"  Splitting:       ΔE = E₁ - E₀ = {splitting:.10f}")
print(f"  Relative split:  ΔE/E_mean = {splitting/mean_Delta:.6f}")
print(f"  2||c||:          {2*np.linalg.norm(c):.10f}")
print(f"  |max ad eval|:   {max(abs(ad_evals)):.10f}")

# Check if splitting ≈ 2||c|| (simple commutator model)
if np.linalg.norm(c) > 1e-10:
    print(f"\n  ΔE / (2||c||) = {splitting / (2*np.linalg.norm(c)):.6f}")
    print(f"  ΔE / |max ad eval| = {splitting / max(abs(ad_evals)):.6f}")


# === Stage 9: Full reconstruction test ===
print("\n" + "=" * 70)
print("Stage 9: Full reconstruction — exp(-H_Popov) vs J_DS")
print("=" * 70)

# Test 1: Direct reconstruction from K + ad(a*)
# K was defined as K = H - ad(a*), so by construction exp(-(K+ad)) = exp(-H) = J
# This is tautological. The real test is whether K has Popov structure.

# Test 2: Does K decompose into a recognizable kinetic operator?
# A Popov kinetic operator for constant-on-fibre perturbations is N·δa
# where N is the Nijenhuis tensor. This should be rank-1 (since dbar Phi is rank 1).

print(f"\n  Test A: Is K rank-1 (Nijenhuis structure)?")
K_sv = np.linalg.svd(K3, compute_uv=False)
print(f"    Singular values of K: [{', '.join(f'{v:.6f}' for v in K_sv)}]")
if abs(K_sv[0]) > 1e-10:
    print(f"    σ₂/σ₁ = {abs(K_sv[1]/K_sv[0]):.6f}")
    print(f"    σ₃/σ₁ = {abs(K_sv[2]/K_sv[0]):.6f}")
print(f"    K is rank-1: {'NO' if abs(K_sv[1]/K_sv[0]) > 0.1 else 'YES'}")

# Test 3: Scalar kinetic model — K = alpha * I on active subspace
# If K ≈ αI, then ad(a*) alone determines the splitting
print(f"\n  Test B: Is K ≈ αI (scalar kinetic)?")
K_trace = np.trace(K3) / 3
K_traceless = K3 - K_trace * np.eye(3)
print(f"    tr(K)/3 = {K_trace:.6f}")
print(f"    ||K - (tr/3)I|| / ||K|| = {np.linalg.norm(K_traceless)/np.linalg.norm(K3):.6f}")
print(f"    K is scalar: {'YES' if np.linalg.norm(K_traceless)/np.linalg.norm(K3) < 0.1 else 'NO'}")

# Test 4: Does the kinetic part relate to ||dbar Phi||?
dbar_norm = np.sqrt(a_star_total_sq)
print(f"\n  Test C: Kinetic eigenvalues vs ||∂̄Φ||")
print(f"    ||∂̄Φ|| (from a*) = {dbar_norm:.6f}")
print(f"    ||J_anti||        = {np.linalg.norm(J_anti):.6f}")
print(f"    K eigenvalues:     {[f'{v:.6f}' for v in sorted(K_evals, key=lambda x: abs(x))]}")

# Test 5: Alternative decomposition — J = J_DS_hol * J_floor_anti
# The DS Jacobian is the holomorphic part; the floor is the anti-holomorphic part.
# In the Popov picture, the holomorphic part is trivial (flat connection),
# and the floor part is the non-integrable deformation.
print(f"\n  Test D: Factorisation J_total = J_floor · J_DS")
print(f"    J_DS  eigenvalues (4x4): {sorted(np.abs(np.linalg.eigvals(J_DS)))[::-1]}")
print(f"    J_floor eigenvalues:     {sorted(np.abs(np.linalg.eigvals(J_floor)))[::-1]}")

# Project J_floor and J_DS separately to L₁=1
J3_DS = P_proj.T @ J_DS @ P_proj
J3_floor = P_proj.T @ J_floor @ P_proj
evals_DS = np.linalg.eigvals(J3_DS)
evals_floor = np.linalg.eigvals(J3_floor)
print(f"    J_DS(3x3) eigenvalues:   {sorted(np.abs(evals_DS))[::-1]}")
print(f"    J_floor(3x3) eigenvalues:{sorted(np.abs(evals_floor))[::-1]}")

# The Popov Hessian should be: H = -ln(J_floor) - ln(J_DS) (if they commuted)
# In practice: H = -ln(J_floor · J_DS)
# Decomposition: -ln(J_DS) is the "holomorphic kinetic" part
#                -ln(J_floor) is the "Nijenhuis correction" part
# Check if -ln(J_floor) ≈ ad(a*) (the commutator/interaction)
# and -ln(J_DS) ≈ K (the kinetic/mass term)

print(f"\n  Test E: DS as kinetic, floor as interaction?")
try:
    H_DS_3 = -logm(J3_DS.astype(complex))
    H_floor_3 = -logm(J3_floor.astype(complex))
    print(f"    -ln(J_DS) eigenvalues:   {sorted(np.abs(np.linalg.eigvals(H_DS_3)))[::-1]}")
    print(f"    -ln(J_floor) eigenvalues:{sorted(np.abs(np.linalg.eigvals(H_floor_3)))[::-1]}")

    # Check: H_DS ≈ K?
    print(f"    ||H_DS - K|| / ||K|| = {np.linalg.norm(H_DS_3 - K3)/np.linalg.norm(K3):.6f}")
    # Check: H_floor ≈ ad(a*)?
    print(f"    ||H_floor - ad(a*)|| / ||ad|| = {np.linalg.norm(H_floor_3 - ad_a)/np.linalg.norm(ad_a):.6f}")
except Exception as ex:
    print(f"    Matrix logarithm failed: {ex}")


# === Summary ===
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
DS transfer operator:
  λ₀ = {evals[0]:.10f}   →  E₀ = {Delta_0:.10f}
  λ₁ = {evals[1]:.10f}   →  E₁ = {Delta_1:.10f}
  λ₂ = {evals[2]:.10f}   →  E₂ = ∞ (conservation law)

Equilibrium (0,1)-form a*:
  ||a*|| = {np.linalg.norm(a_star):.6f}
  su(2) content ||c|| = {np.linalg.norm(c):.6f}
  Scalar fraction: {scalar_frac:.4f}

Adjoint commutator ad(a*):
  Max |eigenvalue| = {max(abs(ad_evals)):.6f}

Popov Hessian H = K + ad(a*):
  K eigenvalues: {[f'{v:.4f}' for v in sorted(K_evals, key=lambda x: -abs(x))]}

IDENTIFICATION TEST:
  If K has Nijenhuis structure (rank-1, from Born floor)
  and ad(a*) provides the eigenvalue splitting,
  then exp(-(K + ad(a*))) = J_DS = exp(-H_Popov),
  and the DS QFT IS quantised Yang-Mills.
""")

print("=" * 70)
print("DONE")
print("=" * 70)
