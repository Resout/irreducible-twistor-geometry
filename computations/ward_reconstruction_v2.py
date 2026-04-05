"""
Ward reconstruction v2: algebraic structure of the F+ coefficient.

From v1: the Ward-reconstructed F+ content is 12.163 per unit fluctuation variance.
The Maurer-Cartan commutator content is 149.596 (= C from thm:condensate).
Ratio: 12.163 / 149.596 = 0.0813.

Question: does 12.163 or the ratio 0.0813 have clean algebraic form at H=3?

The condensate coefficient was C * det(M*)^2 = 8332/625 = 4(3n^2+2n+3)/(n-1)^2
with n = H^3-1 = 26.

Is there a similar formula for the Ward coefficient?
"""

import numpy as np
from numpy.linalg import inv, norm, det, svd
from scipy.optimize import brentq
from fractions import Fraction
import sys

H = 3; FLOOR = 1.0/27.0

sigma = [
    np.array([[0,1],[1,0]], dtype=complex),
    np.array([[0,-1j],[1j,0]], dtype=complex),
    np.array([[1,0],[0,-1]], dtype=complex)
]
I2 = np.eye(2, dtype=complex)
basis = [sigma[0], sigma[1], sigma[2], I2]

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
    s1,s2,s3,th = m
    return (th*I2 + s1*sigma[0] + s2*sigma[1] + s3*sigma[2]) / np.sqrt(2)

def find_equilibrium():
    TARGET_K = 7.0/30.0
    def K_residual(p):
        pw = (1-p)/2; sc = 26./27.
        raw = np.array([np.sqrt(p*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR)])
        e = raw/np.sum(raw)
        m = np.array([0.4, 0.2, 0.2, 0.2])
        for _ in range(5000): m, _ = full_step(m, e)
        _, K = ds_combine(m, e)
        return K - TARGET_K
    p_dom = brentq(K_residual, 0.92, 0.94, xtol=1e-15)
    pw = (1-p_dom)/2; sc = 26./27.
    raw = np.array([np.sqrt(p_dom*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR)])
    e_star = raw/np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(5000): m, _ = full_step(m, e_star)
    return m, e_star

m_star, e_star = find_equilibrium()
_, K_star = ds_combine(m_star, e_star)

print("="*70)
print("WARD F+ COEFFICIENT: ALGEBRAIC STRUCTURE")
print("="*70)
print(f"m* = {m_star}")
print(f"K* = {K_star:.15f}")

# ============================================================
# Compute the floor correction epsilon = M(m*) - M(DS(m*,e*))
# ============================================================

m_ds, _ = ds_combine(m_star, e_star)
M_star = to_matrix(m_star)
M_ds = to_matrix(m_ds)
epsilon = M_star - M_ds

# M_ds inverse
M_ds_inv = inv(M_ds)
det_M_ds = det(M_ds)
det_M_star = det(M_star)

print(f"\ndet(M*) = {det_M_star:.10f}")
print(f"det(M_ds) = {det_M_ds:.10f}")
print(f"||epsilon|| = {norm(epsilon):.10f}")
print(f"||epsilon/M_ds|| = {norm(M_ds_inv @ epsilon):.10f}")

# ============================================================
# Analytical Jacobians
# ============================================================

s1, s2, s3, th = m_star
e1, e2, e3, ph = e_star
oneMinusK = 1 - K_star

# DS pre-normalisation values
S_vals = np.array([s1*(e1+ph)+th*e1, s2*(e2+ph)+th*e2, s3*(e3+ph)+th*e3, th*ph])
N = S_vals / oneMinusK  # DS output (before floor)

# DS Jacobian
dSdm = np.zeros((4,4))
for i in range(3):
    dSdm[i,i] = e_star[i] + ph
    dSdm[i,3] = e_star[i]
dSdm[3,3] = ph
dKdm = np.array([e2+e3, e1+e3, e1+e2, 0.0])

J_DS = np.zeros((4,4))
for i in range(4):
    for j in range(4):
        J_DS[i,j] = (dSdm[i,j] * oneMinusK + S_vals[i] * dKdm[j]) / oneMinusK**2

# Floor Jacobian
Ss = N[0] + N[1] + N[2]
Sq = N[0]**2 + N[1]**2 + N[2]**2
R = Sq / Ss**2
u = np.sqrt(26 * R)
t = (u - R) / (26 - R)
g = (1 - t) / Ss

dtdR = (338 + 13*R - 26*u) / (u * (26 - R)**2)
dRdN = np.array([2*(N[i]*Ss - Sq) / Ss**3 for i in range(3)])
dtdN = np.array([dtdR * dRdN[i] for i in range(3)])
dgdN = np.array([(-dtdN[i]*Ss - (1-t)) / Ss**2 for i in range(3)])

J_floor = np.zeros((4,4))
for i in range(3):
    for j in range(3):
        J_floor[i,j] = (g if i == j else 0.0) + N[i] * dgdN[j]
J_floor[3,:3] = dtdN
J_floor[3,3] = 0.0

# Floor correction Jacobian: d(epsilon)/dm = d(M_floor - M_DS)/dm
# = (dM/dm)(J_floor - I)(J_DS)... wait, epsilon = M(floor(DS(m))) - M(DS(m))
# So d(epsilon)/dm = (dM/dm)|_{m*} . J_total - (dM/dm)|_{DS(m*)} . J_DS
# where J_total = J_floor . J_DS

# Actually: epsilon = M(m*) - M(m_ds) where m* = floor(m_ds)
# d(epsilon)/dm evaluated at m* gives the SPATIAL derivative of epsilon
# when the input m varies spatially.
# epsilon(m) = M(floor(DS(m, e*))) - M(DS(m, e*))
# d(epsilon)/dm = (dM/dm)|_{m*} . J_total - (dM/dm)|_{m_ds} . J_DS

# But (dM/dm)_j = basis_j / sqrt(2) is the SAME linear map regardless of the point.
# So d(epsilon)/dm = (1/sqrt(2)) sum_j basis_j . [J_total - J_DS]_{j,.}
# = (1/sqrt(2)) sum_j basis_j . [(J_floor - I) . J_DS]_{j,.}

J_correction = (J_floor - np.eye(4)) @ J_DS

# ============================================================
# F+ coefficient: Tr(F+† F+) per unit variance
# ============================================================

# For a Gaussian perturbation delta_m with <delta_m_i delta_m_j> = sigma^2 delta_ij,
# the Ward-reconstructed F+ in direction mu is:
#   F+_mu = M_ds^{-1} . sum_j (basis_j/sqrt(2)) . J_correction_{j,mu} . delta_m_mu
#
# The average |F+|^2 uses the Wick theorem for <delta_m_mu delta_m_nu>:
#   <Tr(F+† F+)> = sigma^4 . sum_{mu,nu} Tr(B_mu† B_nu) delta_{mu,nu}
# where B_mu = M_ds^{-1} . sum_j (basis_j/sqrt(2)) . J_correction_{j,mu}

# Actually, for the curvature from spatial variation, we need F_{mu,nu} = [A_mu, A_nu]
# at leading order (the da term and the [a,a] term both contribute for the Ward connection).
# But at leading order in epsilon, the Ward connection is:
#   A_mu = M_ds^{-1} . d(epsilon)_mu = M_ds^{-1} . sum_j (basis_j/sqrt(2)) . J_correction_{j,mu} . delta_m_mu
#
# For INDEPENDENT Gaussian perturbations in each spatial direction:
#   F_{mu,nu} = d_mu A_nu - d_nu A_mu + [A_mu, A_nu]
# At leading order in fluctuations, the [A,A] term is O(delta_m^2) and the dA terms are O(delta_m).
# But dA requires second derivatives of the perturbation... hmm.
#
# Actually, for the WARD-reconstructed connection, the curvature is NOT computed
# from A = h+^{-1} dh+ by the standard F = dA + A^2 formula. The Penrose transform
# gives F DIRECTLY from the cohomology class, without going through A.
#
# The leading-order F+ from Mason's formula:
#   F+ proportional to M_hol^{-1} . dbar(epsilon)
# where dbar(epsilon) is the anti-holomorphic derivative of the floor correction.
#
# For our case: at a single spacetime point, dbar(epsilon) involves the spatial
# variation of epsilon (how the floor correction changes between neighboring sites).
# This is exactly J_correction applied to the spatial perturbation.

# The F+ COEFFICIENT (per unit fluctuation variance):
# C_Ward = sum_{i,k} Tr(B_i† B_k . B_k† B_i) -- no, this is commutator-like
#
# Actually, following the same logic as the condensate theorem:
# For independent Gaussian perturbations v, w in R^4 with <v_i v_j> = sigma^2 delta_ij,
# <|[B_v, B_w]|^2> = C_Ward . sigma^4
# where B_v = M_ds^{-1} . sum_j (basis_j/sqrt(2)) . sum_mu J_correction_{j,mu} v_mu
#
# But this is the commutator content of the WARD-corrected connection, not of M^{-1}dM.

# Let me compute both:
# (a) The norm-squared of the Ward connection components (analogous to |A|^2)
# (b) The commutator content (analogous to |[A,A]|^2 = |F|^2 for Ward)

# Ward connection basis vectors: B_mu = M_ds^{-1} . sum_j (basis_j/sqrt(2)) . J_correction_{j,mu}
B = np.zeros((4, 2, 2), dtype=complex)
for mu in range(4):
    for j in range(4):
        B[mu] += M_ds_inv @ (basis[j] / np.sqrt(2)) * J_correction[j, mu]

# (a) Connection norm: sum_mu Tr(B_mu† B_mu)
A_norm_sq = sum(np.real(np.trace(B[mu].conj().T @ B[mu])) for mu in range(4))
print(f"\n||A_Ward||^2 (connection norm) = {A_norm_sq:.10f}")

# (b) Commutator content: sum_{mu<nu} Tr([B_mu, B_nu]† [B_mu, B_nu])
comm_content = 0.0
for mu in range(4):
    for nu in range(mu+1, 4):
        comm = B[mu] @ B[nu] - B[nu] @ B[mu]
        comm_content += np.real(np.trace(comm.conj().T @ comm))

print(f"Ward commutator content = {comm_content:.10f}")

# (c) Full Wick contraction for <|F|^2> with Gaussian fluctuations
# Following the condensate theorem:
# <|F_{mu,nu}|^2> = <|[A_mu, A_nu]|^2> for independent perturbations
# = sum_{i,k} Tr([B_i, B_k]† [B_i, B_k]) by Wick
# This IS comm_content above (summed over mu<nu pairs).

# But the FULL F = dA + [A,A]. For the Ward connection, at leading order in epsilon:
# dA is O(d^2(epsilon)) and [A,A] is O(epsilon^2). Both are second order.
# The leading-order curvature from Mason's formula is F+ ~ epsilon itself (not its derivative).
# So the direct norm of B gives F+, not the commutator of B.

# Let me compute the direct F+ norm (what Mason's formula gives):
# F+_mu = B_mu (the Ward connection component in direction mu)
# <|F+|^2> = sum_mu Tr(B_mu† B_mu) * sigma^4 = A_norm_sq * sigma^4

print(f"\n--- COMPARISON ---")
print(f"Ward connection norm ||A||^2:     {A_norm_sq:.8f}")
print(f"Ward commutator |[A,A]|^2:       {comm_content:.8f}")
print(f"MC commutator content C:          149.596")
print(f"Ratio Ward_conn/MC_comm:          {A_norm_sq/149.596:.8f}")
print(f"Ratio Ward_comm/MC_comm:          {comm_content/149.596:.8f}")

# ============================================================
# Check algebraic structure of the Ward coefficient
# ============================================================
print(f"\n{'='*70}")
print("ALGEBRAIC STRUCTURE CHECK")
print(f"{'='*70}")

# The MC condensate: C * det(M*)^2 = 8332/625
# Check: Ward_coefficient * det(M_ds)^2 = ?
det_M_ds_sq = abs(det_M_ds)**2

ward_times_det_sq = A_norm_sq * det_M_ds_sq
comm_times_det_sq = comm_content * det_M_ds_sq

print(f"\ndet(M_ds)^2 = {det_M_ds_sq:.10f}")
print(f"det(M*)^2   = {abs(det_M_star)**2:.10f}")

print(f"\nWard_conn * det(M_ds)^2  = {ward_times_det_sq:.10f}")
print(f"Ward_comm * det(M_ds)^2  = {comm_times_det_sq:.10f}")
print(f"MC_comm * det(M*)^2      = {149.596 * abs(det_M_star)**2:.10f}")

# Check simple fractions near these values
for name, val in [("Ward_conn * det(M_ds)^2", ward_times_det_sq),
                  ("Ward_comm * det(M_ds)^2", comm_times_det_sq)]:
    print(f"\n{name} = {val:.10f}")
    # Try to identify as p/q for small p,q
    best_frac = None; best_err = 1.0
    for q in range(1, 1000):
        p = round(val * q)
        err = abs(val - p/q)
        if err < best_err:
            best_err = err
            best_frac = (p, q)
        if err < 1e-6:
            print(f"  ≈ {p}/{q} = {p/q:.10f}  (error = {err:.2e})")
            break
    else:
        if best_frac:
            p, q = best_frac
            print(f"  Best: {p}/{q} = {p/q:.10f}  (error = {abs(val-p/q):.2e})")

# Also check ratios to known constants
print(f"\nRatios to framework constants:")
print(f"  Ward_conn / (7/30) = {A_norm_sq / (7/30):.8f}")
print(f"  Ward_conn / (23/30) = {A_norm_sq / (23/30):.8f}")
print(f"  Ward_conn / (4/27) = {A_norm_sq / (4/27):.8f}")
print(f"  Ward_conn * K* = {A_norm_sq * 7/30:.8f}")
print(f"  Ward_conn * (1-K*) = {A_norm_sq * 23/30:.8f}")
print(f"  Ward_comm / K*^2 = {comm_content / (7/30)**2:.8f}")
print(f"  Ward_comm / (K*(1-K*)) = {comm_content / (7/30 * 23/30):.8f}")

# The Frobenius identity: K = 1/2 ||M-E||^2 + 1/2 K_self(m) + 1/2 K_self(e) - (th-ph)^2
# Check at equilibrium:
K_self_m = (1 - m_star[3])**2 - np.sum(m_star[:3]**2)
K_self_e = (1 - e_star[3])**2 - np.sum(e_star[:3]**2)
M_E_dist_sq = norm(M_star - to_matrix(e_star))**2
th_diff_sq = (m_star[3] - e_star[3])**2

K_from_identity = 0.5 * M_E_dist_sq + 0.5 * K_self_m + 0.5 * K_self_e - th_diff_sq
print(f"\nFrobenius identity check:")
print(f"  K_self(m*) = {K_self_m:.10f}")
print(f"  K_self(e*) = {K_self_e:.10f}")
print(f"  ||M*-E*||^2 = {M_E_dist_sq:.10f}")
print(f"  K from identity = {K_from_identity:.10f}")
print(f"  K* = {K_star:.10f}")
print(f"  Match: {abs(K_from_identity - K_star) < 1e-10}")

# The floor correction in terms of known quantities
print(f"\nFloor correction structure:")
print(f"  ||epsilon||^2 = {norm(epsilon)**2:.10f}")
print(f"  ||epsilon||^2 / det(M*)^2 = {norm(epsilon)**2 / abs(det_M_star)**2:.10f}")
print(f"  ||epsilon||^2 / K* = {norm(epsilon)**2 / (7/30):.10f}")
print(f"  ||epsilon||^2 / (1/27) = {norm(epsilon)**2 / (1/27):.10f}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
