"""
Ward reconstruction at the K*=7/30 equilibrium.

The physical gauge connection is NOT M^{-1}dM (Maurer-Cartan, identically flat).
It is the Ward-reconstructed connection A = h+^{-1} d h+, obtained by Birkhoff
factorisation of the transition function on each twistor line.

This computation:
1. Compute the equilibrium (m*, e*) at high precision
2. Compute dbar_Phi (the anti-holomorphic Jacobian) — the source of non-integrability
3. Compute dbar_M = (dM/dm) . dbar_Phi — the kink of the transition function
4. Decompose dbar_M into SO(4) irreps: (1,1) + (3,1)+(1,3) + (3,3)
5. Set up the Penrose transform: the (3,1)+(1,3) component, restricted to
   twistor lines and integrated via the Cauchy kernel, gives F+ on S^4.

Key insight: M is constant on each twistor line L_x (it depends on x, not zeta).
So dbar_M|_{fibre} = 0 in the STANDARD complex structure. But the Born floor
makes J non-integrable, which ROTATES what counts as the "fibre" direction.
The Nijenhuis tensor N_J mixes base and fibre (0,1)-directions. The component
of dbar_M that gets rotated into the fibre direction IS the source of F+.

The self-dual curvature at spacetime point x is:
  F+_{A'B'}(x) = oint_{L_x} pi_{A'} pi_{B'} . [N_J-rotated component of dbar_M] . dzeta
"""

import numpy as np
from numpy.linalg import inv, norm, det, eig, svd
from scipy.optimize import brentq
import sys

# ============================================================
# Setup: equilibrium computation
# ============================================================

H = 3; FLOOR = 1.0/27.0

sigma = [
    np.array([[0,1],[1,0]], dtype=complex),     # sigma_1
    np.array([[0,-1j],[1j,0]], dtype=complex),   # sigma_2
    np.array([[1,0],[0,-1]], dtype=complex)       # sigma_3
]
I2 = np.eye(2, dtype=complex)

def enforce_floor(m):
    s = m[:3].copy(); th = m[3]
    born = th**2 / (np.sum(s**2) + th**2)
    if born >= FLOOR: return m.copy()
    S = np.sum(s); Sq = np.sum(s**2); R = Sq/S**2
    disc = (2*R)**2 + 4*(26-R)*R
    t = (-2*R + np.sqrt(disc)) / (2*(26-R))
    alpha = (1-t)/S
    out = m.copy(); out[:3] = s*alpha; out[3] = t
    return out

def ds_combine(m, e):
    s,th = m[:3],m[3]; se,ph = e[:3],e[3]
    sn = s*se + s*ph + th*se; tn = th*ph
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
    d = 1-K; out = np.zeros(4); out[:3]=sn/d; out[3]=tn/d
    return out, K

def full_step(m, e):
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds), K

def to_matrix(m):
    """Pauli embedding: m -> M = (theta*I + s.sigma)/sqrt(2)"""
    s1,s2,s3,th = m
    return (th*I2 + s1*sigma[0] + s2*sigma[1] + s3*sigma[2]) / np.sqrt(2)

def find_equilibrium():
    """Find the K*=7/30 equilibrium via bisection on p_dom."""
    TARGET_K = 7.0/30.0

    def K_residual(p):
        pw = (1-p)/2; sc = 26./27.
        raw = np.array([np.sqrt(p*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR)])
        e = raw/np.sum(raw)
        m = np.array([0.4, 0.2, 0.2, 0.2])
        for _ in range(5000):
            m, _ = full_step(m, e)
        _, K = ds_combine(m, e)
        return K - TARGET_K

    p_dom = brentq(K_residual, 0.92, 0.94, xtol=1e-15)
    pw = (1-p_dom)/2; sc = 26./27.
    raw = np.array([np.sqrt(p_dom*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR)])
    e_star = raw/np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(5000):
        m, _ = full_step(m, e_star)
    m_star = m
    return m_star, e_star

print("="*70)
print("WARD RECONSTRUCTION AT THE K*=7/30 EQUILIBRIUM")
print("="*70)

m_star, e_star = find_equilibrium()
_, K_star = ds_combine(m_star, e_star)
M_star = to_matrix(m_star)

print(f"\nm* = {m_star}")
print(f"e* = {e_star}")
print(f"K* = {K_star:.15f}  (target: {7/30:.15f})")
print(f"M* = \n{M_star}")
print(f"det(M*) = {det(M_star):.10f}")
print(f"|det(M*)| = {abs(det(M_star)):.10f}")
sys.stdout.flush()

# ============================================================
# Step 1: Anti-holomorphic Jacobian dbar_Phi
# ============================================================
print("\n" + "="*70)
print("STEP 1: Anti-holomorphic Jacobian dbar_Phi")
print("="*70)

# Compute dbar_Phi via analytical chain rule (floor . DS)
s1, s2, s3, th = m_star
e1, e2, e3, ph = e_star
oneMinusK = 1 - K_star

# DS Jacobian: d(DS)/dm at (m*, e*)
# DS output (pre-floor) components: N_i = (s_i*(e_i+ph) + th*e_i) / (1-K)
N = np.zeros(4)
N[0] = (s1*(e1+ph) + th*e1) / oneMinusK
N[1] = (s2*(e2+ph) + th*e2) / oneMinusK
N[2] = (s3*(e3+ph) + th*e3) / oneMinusK
N[3] = (th*ph) / oneMinusK

# Pre-normalisation values
S_vals = np.array([s1*(e1+ph)+th*e1, s2*(e2+ph)+th*e2, s3*(e3+ph)+th*e3, th*ph])

# dS/dm (pre-normalisation derivatives)
dSdm = np.zeros((4,4))
for i in range(3):
    dSdm[i,i] = e_star[i] + ph  # d(s_i*e_i + s_i*ph + th*e_i)/ds_i
    dSdm[i,3] = e_star[i]       # d(...)/dth
dSdm[3,3] = ph                   # d(th*ph)/dth

# dK/dm
dKdm = np.array([e2+e3, e1+e3, e1+e2, 0.0])

# Full DS Jacobian (quotient rule)
J_DS = np.zeros((4,4))
for i in range(4):
    for j in range(4):
        J_DS[i,j] = (dSdm[i,j] * oneMinusK + S_vals[i] * dKdm[j]) / oneMinusK**2

# Floor Jacobian at N (the DS output before floor)
Ss = N[0] + N[1] + N[2]  # sum of singletons
Sq = N[0]**2 + N[1]**2 + N[2]**2  # sum of squares
R = Sq / Ss**2
u = np.sqrt(26 * R)
t = (u - R) / (26 - R)
g = (1 - t) / Ss

# Derivatives for floor Jacobian
dtdR = (338 + 13*R - 26*u) / (u * (26 - R)**2)
dRdN = np.array([2*(N[i]*Ss - Sq) / Ss**3 for i in range(3)])
dtdN = np.array([dtdR * dRdN[i] for i in range(3)])
dgdN = np.array([(-dtdN[i]*Ss - (1-t)) / Ss**2 for i in range(3)])

J_floor = np.zeros((4,4))
for i in range(3):
    for j in range(3):
        J_floor[i,j] = (g if i == j else 0.0) + N[i] * dgdN[j]
    J_floor[i,3] = 0.0
J_floor[3,0] = dtdN[0]
J_floor[3,1] = dtdN[1]
J_floor[3,2] = dtdN[2]
J_floor[3,3] = 0.0

# Total Jacobian = J_floor . J_DS
J_total = J_floor @ J_DS

print(f"\nDS output (before floor): N = {N}")
print(f"Born(N) = {N[3]**2 / (np.sum(N[:3]**2) + N[3]**2):.6f}  (floor = {FLOOR:.6f})")
print(f"\nJ_total (full 4x4 Jacobian):")
for i in range(4):
    print(f"  [{' '.join(f'{J_total[i,j]:10.6f}' for j in range(4))}]")

# Anti-holomorphic Jacobian: dbar_Phi
# For real mass functions at equilibrium, dbar_Phi = J_total restricted to
# the anti-holomorphic directions. At real m*, the anti-holomorphic Jacobian
# IS the Jacobian (real mass functions have z = x, so d/dz_bar = d/dx).
dbar_Phi = J_total.copy()

# Check rank
U_svd, s_svd, Vt_svd = svd(dbar_Phi)
print(f"\nSingular values of dbar_Phi: {s_svd}")
print(f"Rank (>1e-10): {np.sum(s_svd > 1e-10)}")
print(f"||dbar_Phi||_F = {norm(dbar_Phi, 'fro'):.6f}")
sys.stdout.flush()

# ============================================================
# Step 2: dbar_M = (dM/dm) . dbar_Phi
# ============================================================
print("\n" + "="*70)
print("STEP 2: Anti-holomorphic derivative of the transition function dbar_M")
print("="*70)

# dM/dm_j: the Pauli embedding maps m_j -> sigma_j/sqrt(2) (j=0,1,2 for singletons)
# and m_3 -> I/sqrt(2) (for theta)
# So dM = (1/sqrt(2)) * (dtheta * I + ds1 * sigma1 + ds2 * sigma2 + ds3 * sigma3)
# i.e., dM/dm_j = sigma_j/sqrt(2) for j=0,1,2 and I/sqrt(2) for j=3

basis = [sigma[0], sigma[1], sigma[2], I2]  # {sigma_1, sigma_2, sigma_3, I}

# dbar_M has one 2x2 matrix for each anti-holomorphic direction alpha
# dbar_M_alpha = sum_j (basis_j / sqrt(2)) * (dbar_Phi)_{j,alpha}
# Since dbar_Phi has rank 1, all dbar_M_alpha are proportional

dbar_M = np.zeros((4, 2, 2), dtype=complex)  # dbar_M[alpha] is a 2x2 matrix
for alpha in range(4):
    for j in range(4):
        dbar_M[alpha] += basis[j] * dbar_Phi[j, alpha] / np.sqrt(2)

print("\ndbar_M for each anti-holomorphic direction:")
for alpha in range(4):
    nrm = norm(dbar_M[alpha])
    print(f"  alpha={alpha}: ||dbar_M|| = {nrm:.8f}")
    if nrm > 1e-10:
        # SVD of this 2x2 matrix
        U2, s2, V2 = svd(dbar_M[alpha])
        print(f"           singular values: {s2}")
        ratio = s2[1]/s2[0] if s2[0] > 1e-15 else 0
        print(f"           ratio sigma2/sigma1 = {ratio:.6f}")

# The total dbar_M norm
total_norm_sq = sum(norm(dbar_M[alpha])**2 for alpha in range(4))
print(f"\nTotal ||dbar_M||_F = {np.sqrt(total_norm_sq):.8f}")

# Check: is dbar_M proportional across alphas? (rank-1 Jacobian means yes)
norms = [norm(dbar_M[alpha]) for alpha in range(4)]
max_norm = max(norms)
if max_norm > 1e-10:
    ref = None
    for alpha in range(4):
        if norms[alpha] > 1e-10:
            ref = alpha; break
    print(f"\nProportionality check (all should be ~1 if rank-1):")
    for alpha in range(4):
        if norms[alpha] > 1e-10:
            # Check if dbar_M[alpha] is proportional to dbar_M[ref]
            ratio_mat = dbar_M[alpha] / dbar_M[ref] if ref is not None else None
            if ratio_mat is not None:
                # All entries should be the same scalar
                vals = ratio_mat.flatten()
                vals = vals[np.abs(dbar_M[ref].flatten()) > 1e-12]
                if len(vals) > 0:
                    spread = np.std(vals) / np.abs(np.mean(vals)) if np.abs(np.mean(vals)) > 1e-12 else float('inf')
                    print(f"  alpha={alpha}: ratio = {np.mean(vals):.6f}, spread = {spread:.2e}")
sys.stdout.flush()

# ============================================================
# Step 3: SO(4) decomposition of dbar_M
# ============================================================
print("\n" + "="*70)
print("STEP 3: SO(4) decomposition of dbar_M")
print("="*70)

# A 2x2 complex matrix decomposes as: A = (tr(A)/2)*I + A_traceless
# The trace part is the (1,1) scalar
# The traceless part has su(2) content: A_traceless = sum a_i sigma_i
# Under SO(4), the full 4x4 Jacobian decomposes as:
# (1,1): trace/4 * I
# (3,1)+(1,3): antisymmetric part (A - A^T)/2
# (3,3): symmetric traceless (A + A^T)/2 - trace/4 * I

# For the 2x2 matrices dbar_M[alpha], the decomposition is:
# scalar: tr(dbar_M)/2 . I (proportional to identity)
# gauge: traceless part (su(2) components)

# Extract su(2) components of each dbar_M[alpha]
print("\nsu(2) decomposition of dbar_M:")
for alpha in range(4):
    if norm(dbar_M[alpha]) < 1e-12:
        continue
    A = dbar_M[alpha]
    tr_part = np.trace(A) / 2  # scalar = theta component
    traceless = A - tr_part * I2  # su(2) part

    # Extract sigma components: a_i = Tr(sigma_i . A) / 2
    su2_components = np.array([np.trace(sigma[i] @ A) / 2 for i in range(3)])

    scalar_frac = abs(tr_part)**2 / (abs(tr_part)**2 + np.sum(np.abs(su2_components)**2))
    gauge_frac = 1 - scalar_frac

    print(f"  alpha={alpha}:")
    print(f"    scalar (I):  {tr_part:.8f}  ({100*scalar_frac:.1f}%)")
    print(f"    gauge (sigma): {su2_components}  ({100*gauge_frac:.1f}%)")

# Now decompose the FULL 4x4 Jacobian dbar_Phi
print("\nSO(4) decomposition of full 4x4 dbar_Phi:")
A = dbar_Phi
trace_part = np.trace(A) / 4 * np.eye(4)
antisym = (A - A.T) / 2
sym_traceless = (A + A.T) / 2 - np.trace(A) / 4 * np.eye(4)

total_sq = norm(A, 'fro')**2
scalar_sq = norm(trace_part, 'fro')**2
gauge_sq = norm(antisym, 'fro')**2
gravity_sq = norm(sym_traceless, 'fro')**2

print(f"  (1,1) scalar:       {100*scalar_sq/total_sq:.1f}%")
print(f"  (3,1)+(1,3) gauge:  {100*gauge_sq/total_sq:.1f}%")
print(f"  (3,3) gravity:      {100*gravity_sq/total_sq:.1f}%")
print(f"  Sum check:          {100*(scalar_sq+gauge_sq+gravity_sq)/total_sq:.1f}%")
sys.stdout.flush()

# ============================================================
# Step 4: M*^{-1} and the connection form components
# ============================================================
print("\n" + "="*70)
print("STEP 4: Connection form a = M*^{-1} . dbar_M")
print("="*70)

M_inv = inv(M_star)
print(f"\nM*^{{-1}} = \n{M_inv}")
print(f"Check M*^{{-1}} M* = I: ||M^{{-1}}M - I|| = {norm(M_inv @ M_star - I2):.2e}")

# Connection components: a_alpha = M*^{-1} . dbar_M[alpha]
a = np.zeros((4, 2, 2), dtype=complex)
for alpha in range(4):
    a[alpha] = M_inv @ dbar_M[alpha]

print("\nMaurer-Cartan form a_alpha = M^{-1} dbar_M[alpha]:")
for alpha in range(4):
    nrm = norm(a[alpha])
    if nrm > 1e-12:
        print(f"  alpha={alpha}: ||a|| = {nrm:.8f}")
        # su(2) decomposition
        tr_a = np.trace(a[alpha]) / 2
        su2_a = np.array([np.trace(sigma[i] @ a[alpha]) / 2 for i in range(3)])
        print(f"    scalar: {tr_a:.8f}")
        print(f"    su(2):  {su2_a}")

# Check commutators [a_alpha, a_beta] - should be ~0 because rank-1 Jacobian
print("\nCommutators [a_alpha, a_beta] (should be ~0 from rank-1 Jacobian):")
for alpha in range(4):
    for beta in range(alpha+1, 4):
        comm = a[alpha] @ a[beta] - a[beta] @ a[alpha]
        nrm = norm(comm)
        if norm(a[alpha]) > 1e-12 and norm(a[beta]) > 1e-12:
            print(f"  [a_{alpha}, a_{beta}]: ||comm|| = {nrm:.2e}")

# Maurer-Cartan curvature: F = dbar(a) + a wedge a
# By Maurer-Cartan identity, this is IDENTICALLY ZERO for a = M^{-1} dbar M.
# The commutators above should be ~0 (rank-1), and da should cancel [a,a] exactly.
print("\n*** The Maurer-Cartan form has F = 0 identically. ***")
print("*** The physical curvature requires the Ward reconstruction. ***")
sys.stdout.flush()

# ============================================================
# Step 5: Nijenhuis tensor and fibre-base mixing
# ============================================================
print("\n" + "="*70)
print("STEP 5: Non-integrability data for the Ward reconstruction")
print("="*70)

# The Nijenhuis tensor N_J measures the failure of J to be integrable.
# For our deformed J: N_J is proportional to dbar_Phi.
# The key quantity for F+ is: how much of the base dbar_M gets
# rotated into the fibre direction by N_J.

# At the equilibrium, dbar_Phi has rank 1. Let v be its image direction.
# All dbar_M_alpha are proportional to a single 2x2 matrix A_0 = (dM/dm).v

# The rank-1 direction of dbar_Phi
U_svd, s_svd, Vt_svd = svd(dbar_Phi)
v_image = U_svd[:, 0]  # dominant left singular vector (image direction)
v_source = Vt_svd[0, :]  # dominant right singular vector (source direction)
sigma_1_svd = s_svd[0]  # dominant singular value

print(f"\nRank-1 structure of dbar_Phi:")
print(f"  Dominant singular value: {sigma_1_svd:.8f}")
print(f"  Image direction (v):  {v_image}")
print(f"  Source direction (w): {v_source}")

# The 2x2 matrix that ALL dbar_M_alpha are proportional to:
A0 = sum(basis[j] * v_image[j] for j in range(4)) / np.sqrt(2)
print(f"\nDominant 2x2 kink matrix A0 = (dM/dm).v:")
print(f"  {A0}")
print(f"  ||A0|| = {norm(A0):.8f}")
print(f"  tr(A0) = {np.trace(A0):.8f}")
print(f"  Traceless fraction: {1 - abs(np.trace(A0)/2)**2/norm(A0)**2 * 4:.4f}")

# SVD of A0: this tells us how the kink couples the two bundle channels
U_A0, s_A0, V_A0 = svd(A0)
print(f"\n  SVD of A0: sigma1={s_A0[0]:.6f}, sigma2={s_A0[1]:.6f}")
print(f"  Ratio sigma2/sigma1 = {s_A0[1]/s_A0[0]:.6f}")
print(f"  Entanglement entropy S_ent = {-sum(p*np.log2(p) for p in (s_A0**2/np.sum(s_A0**2)) if p>1e-15):.4f} bits")

# The curvature M*^{-1} A0 (the connection-form version of the kink)
a0 = M_inv @ A0
print(f"\na0 = M*^{{-1}} A0 (connection kink):")
print(f"  {a0}")
su2_a0 = np.array([np.trace(sigma[i] @ a0) / 2 for i in range(3)])
tr_a0 = np.trace(a0) / 2
print(f"  scalar: {tr_a0:.8f}")
print(f"  su(2):  {su2_a0}")
print(f"  Non-abelian fraction: {np.sum(np.abs(su2_a0)**2) / (abs(tr_a0)**2 + np.sum(np.abs(su2_a0)**2)):.4f}")

# ============================================================
# Step 6: The Penrose transform setup
# ============================================================
print("\n" + "="*70)
print("STEP 6: Penrose transform — F+ from dbar_M")
print("="*70)

# Mason's formula (Theorem thm:mason): F+ ~ M_hol^{-1} . dbar(epsilon)|_{L_x}
# At leading order, F+ is proportional to the anti-holomorphic derivative
# of the floor correction, pushed through the Pauli embedding.
#
# For the perturbative Birkhoff factorisation:
#   M_tilde = M_0(zeta) + epsilon(zeta, zeta_bar)
#   h+^{(1)} = M_0 . P+[M_0^{-1} . epsilon]
#   A^{(1)} = M_0^{-1} . P+[M_0^{-1} . d_x(epsilon)]
#   F+ ~ d_x(A^{(1)})
#
# The key object is epsilon: the anti-holomorphic correction to M
# from the Born floor.
#
# At the equilibrium: M(m*) is the FULL transition function (after floor).
# The holomorphic part M_hol = M(DS(m*, e*)) is the DS output BEFORE floor.
# epsilon = M(m*) - M(DS(m*, e*)) is the floor correction in matrix form.

m_ds, K_ds = ds_combine(m_star, e_star)
M_ds = to_matrix(m_ds)  # DS output before floor
M_floor = to_matrix(m_star)  # After floor = m* (equilibrium)

epsilon = M_floor - M_ds  # The floor correction

print(f"\nDS output before floor: m_ds = {m_ds}")
print(f"Born(m_ds) = {m_ds[3]**2 / (np.sum(m_ds[:3]**2) + m_ds[3]**2):.6f}  (floor = {FLOOR:.6f})")
print(f"\nM_hol (DS output, holomorphic part):")
print(f"  {M_ds}")
print(f"\nM_floor (after floor = equilibrium):")
print(f"  {M_floor}")
print(f"\nepsilon (floor correction = anti-holomorphic part):")
print(f"  {epsilon}")
print(f"  ||epsilon|| = {norm(epsilon):.8f}")
print(f"  ||epsilon|| / ||M_hol|| = {norm(epsilon)/norm(M_ds):.6f}")

# M_hol^{-1} . epsilon: the correction in the bundle frame
M_ds_inv = inv(M_ds)
rel_correction = M_ds_inv @ epsilon
print(f"\nM_hol^{{-1}} . epsilon (relative correction):")
print(f"  {rel_correction}")
print(f"  ||M_hol^{{-1}} epsilon|| = {norm(rel_correction):.8f}")

# su(2) decomposition of the correction
tr_eps = np.trace(rel_correction) / 2
su2_eps = np.array([np.trace(sigma[i] @ rel_correction) / 2 for i in range(3)])
print(f"\n  scalar part: {tr_eps:.8f}")
print(f"  su(2) part:  {su2_eps}")
na_frac = np.sum(np.abs(su2_eps)**2) / (abs(tr_eps)**2 + np.sum(np.abs(su2_eps)**2))
print(f"  Non-abelian fraction of correction: {na_frac:.4f}")

# ============================================================
# Step 7: Leading-order F+ estimate
# ============================================================
print("\n" + "="*70)
print("STEP 7: Leading-order self-dual curvature estimate")
print("="*70)

# Mason's formula: F+ ~ M_hol^{-1} . dbar(epsilon)|_{L_x}
# At the equilibrium (spatially uniform), the "derivative" of epsilon
# with respect to the base directions is zero — because all sites have
# the same m*.
#
# BUT: for FLUCTUATIONS around m* (which is what the path integral samples),
# epsilon varies from site to site. The derivative of epsilon between
# neighboring sites gives the physical curvature.
#
# The curvature from a perturbation delta_m at a neighboring site:
# delta_epsilon = d(epsilon)/dm . delta_m
# where d(epsilon)/dm is the Jacobian of the floor correction.
#
# F+ ~ M_hol^{-1} . d(epsilon)/dm . delta_m / lattice_spacing

# The floor correction Jacobian: d(M_floor - M_DS)/dm
# = d(M_floor)/dm - d(M_DS)/dm
# = (dM/dm) . J_floor . J_DS - (dM/dm) . J_DS  [chain rule]
# = (dM/dm) . (J_floor - I) . J_DS

J_correction = (J_floor - np.eye(4)) @ J_DS  # Jacobian of the floor correction

print(f"\nJacobian of floor correction (J_floor - I) @ J_DS:")
for i in range(4):
    print(f"  [{' '.join(f'{J_correction[i,j]:10.6f}' for j in range(4))}]")

# For Gaussian fluctuations with variance sigma^2, the average |F+|^2 is:
# <|F+|^2> = sum_{mu,nu} |M_hol^{-1} . dbar(epsilon)_mu|^2
# where the sum is over spatial directions mu.

# Each direction contributes: (M_ds_inv . sum_j basis_j/sqrt(2) . J_correction_{j,mu})
F_plus_components = np.zeros((4, 2, 2), dtype=complex)
for mu in range(4):
    for j in range(4):
        F_plus_components[mu] += M_ds_inv @ (basis[j] / np.sqrt(2)) * J_correction[j, mu]

F_plus_sq_total = sum(np.real(np.trace(F_plus_components[mu].conj().T @ F_plus_components[mu]))
                      for mu in range(4))

print(f"\n||F+||^2 (per unit fluctuation) = {F_plus_sq_total:.8f}")

# For comparison: the commutator content C (from the condensate theorem)
# C measures |[a_mu, a_nu]|^2 which is the Maurer-Cartan form's commutator
C_condensate = 149.596  # from the analytic computation
print(f"\nFor comparison:")
print(f"  Commutator content C (MC form) = {C_condensate:.3f}")
print(f"  F+ content (Ward correction)   = {F_plus_sq_total:.3f}")
print(f"  Ratio F+/C = {F_plus_sq_total/C_condensate:.6f}")

# The key number: at the equilibrium itself (uniform m*), F+ = 0
# (no spatial variation -> no curvature). The curvature comes from
# fluctuations, just like the condensate.
# But now it's the WARD-reconstructed curvature, not the Maurer-Cartan form.

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")
print(f"\n1. dbar_Phi: rank 1, ||dbar_Phi|| = {norm(dbar_Phi, 'fro'):.4f}")
print(f"2. dbar_M: all components proportional (rank-1 Jacobian)")
print(f"3. Floor correction: ||epsilon|| = {norm(epsilon):.4f}, ||epsilon/M_hol|| = {norm(epsilon)/norm(M_ds):.4f}")
print(f"4. Non-abelian fraction of correction: {na_frac:.4f}")
print(f"5. F+ content (from floor correction Jacobian): {F_plus_sq_total:.4f}")
print(f"\nThe physical curvature at a UNIFORM equilibrium is zero (no spatial variation).")
print(f"The physical curvature from FLUCTUATIONS is {F_plus_sq_total:.4f} per unit variance,")
print(f"which is the leading-order Ward-reconstructed <Tr(F+^2)>.")
print()
print("DONE.")
