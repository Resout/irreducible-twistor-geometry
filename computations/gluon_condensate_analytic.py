"""
Gluon condensate: ANALYTIC computation of C = <|[a,a]|^2> / sigma^4.

For small perturbations delta_m around m*:
  delta_M = (delta_th * I + delta_s . sigma) / sqrt(2)
  a = M*^{-1} delta_M

The commutator [a_mu, a_nu] for two independent perturbations:
  [a_mu, a_nu] = M*^{-1} [delta_M_mu, M*^{-1} delta_M_nu]
               + terms that cancel in the commutator

Actually, more carefully:
  a_mu = M*^{-1} delta_M_mu
  a_nu = M*^{-1} delta_M_nu
  [a_mu, a_nu] = M*^{-1} delta_M_mu M*^{-1} delta_M_nu - M*^{-1} delta_M_nu M*^{-1} delta_M_mu
               = M*^{-1} (delta_M_mu M*^{-1} delta_M_nu - delta_M_nu M*^{-1} delta_M_mu)

For the Pauli embedding: delta_M = (v_0 I + v_1 sigma_1 + v_2 sigma_2 + v_3 sigma_3)/sqrt(2)
where v = (delta_th, delta_s1, delta_s2, delta_s3).

The commutator [delta_M_mu, A] for any 2x2 matrix A:
  [v.sigma, A] involves the su(2) structure constants.

C = <Tr([a_mu, a_nu]^dag [a_mu, a_nu])> / sigma^4
where the average is over independent Gaussian v_mu, v_nu with variance sigma^2.

This is a quartic moment of a Gaussian — computable exactly via Isserlis/Wick theorem.
"""

import numpy as np
from numpy.linalg import inv, norm
from scipy.optimize import brentq
import sys

H = 3; FLOOR = 1.0/27.0

sigma_mat = [
    np.array([[0,1],[1,0]], dtype=complex),
    np.array([[0,-1j],[1j,0]], dtype=complex),
    np.array([[1,0],[0,-1]], dtype=complex)
]
I2 = np.eye(2, dtype=complex)

def enforce_floor(m):
    s = m[:3].copy(); th = m[3]
    born = th**2 / (np.sum(s**2) + th**2)
    if born >= FLOOR: return m.copy()
    S = np.sum(s); Sq = np.sum(s**2); R = Sq/S**2
    t = (np.sqrt(26*R) - R)/(26-R); alpha = (1-t)/S
    out = m.copy(); out[:3] = s*alpha; out[3] = t
    return out

def to_matrix(m):
    s1,s2,s3,th = m
    return (th*I2 + s1*sigma_mat[0] + s2*sigma_mat[1] + s3*sigma_mat[2]) / np.sqrt(2)

def find_eq():
    def ds_step(m, e):
        s,th = m[:3],m[3]; se,ph = e[:3],e[3]
        sn = s*se + s*ph + th*se; tn = th*ph
        K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
        d = 1-K; out = np.zeros(4); out[:3]=sn/d; out[3]=tn/d
        return enforce_floor(out)
    def K_res(p):
        pw=(1-p)/2; sc=26./27.
        raw = np.array([np.sqrt(p*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR)])
        e = raw/np.sum(raw); m = np.array([0.4,0.2,0.2,0.2])
        for _ in range(3000): m = ds_step(m,e)
        K = sum(m[i]*e[j] for i in range(3) for j in range(3) if i!=j)
        return K - 7./30.
    p = brentq(K_res, 0.92, 0.94, xtol=1e-15)
    pw=(1-p)/2; sc=26./27.
    raw = np.array([np.sqrt(p*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR)])
    e = raw/np.sum(raw); m = np.array([0.4,0.2,0.2,0.2])
    for _ in range(3000): m = ds_step(m,e)
    return m, e

m_star, e_star = find_eq()
M_star = to_matrix(m_star)
M_inv = inv(M_star)

print(f"m* = {m_star}")
print(f"M* = \n{M_star}")
print(f"M*^{{-1}} = \n{M_inv}")
print(f"det(M*) = {np.linalg.det(M_star):.8f}")
print(); sys.stdout.flush()

# ============================================================
# ANALYTIC STRUCTURE
# ============================================================
#
# Let v = (v0, v1, v2, v3) be a random perturbation with <v_i v_j> = sigma^2 delta_{ij}.
# delta_M(v) = (v0*I + v1*s1 + v2*s2 + v3*s3) / sqrt(2)
# a(v) = M*^{-1} delta_M(v)
#
# For two independent perturbations v (direction mu) and w (direction nu):
# [a(v), a(w)] = M*^{-1} delta_M(v) M*^{-1} delta_M(w) - M*^{-1} delta_M(w) M*^{-1} delta_M(v)
#
# |[a,a]|^2 = Tr([a(v),a(w)]^dag [a(v),a(w)])
#
# This is QUARTIC in v,w. The Gaussian average factorizes by Wick's theorem:
# <v_i v_j w_k w_l> = <v_i v_j><w_k w_l> + <v_i w_k><v_j w_l> + <v_i w_l><v_j w_k>
# Since v and w are independent: <v_i w_k> = 0
# So: <v_i v_j w_k w_l> = sigma^4 * delta_{ij} * delta_{kl}
#
# Therefore:
# <|[a,a]|^2> = sigma^4 * sum_{i,j,k,l} delta_{ij} delta_{kl} * [coefficient of v_i v_j w_k w_l in |[a,a]|^2]
#             = sigma^4 * sum_{i,k} [coefficient of v_i^2 w_k^2 in |[a,a]|^2]
#
# So C = sum_{i=0}^{3} sum_{k=0}^{3} Tr([A_i, A_k]^dag [A_i, A_k])
# where A_i = M*^{-1} * (basis_matrix_i) / sqrt(2)
# and basis_matrix_0 = I, basis_matrix_{1,2,3} = sigma_{1,2,3}

print("="*60)
print("ANALYTIC COMPUTATION OF C")
print("="*60)
print(); sys.stdout.flush()

# Basis matrices for the Pauli embedding
basis = [I2, sigma_mat[0], sigma_mat[1], sigma_mat[2]]

# A_i = M*^{-1} basis_i / sqrt(2)
A = [M_inv @ basis[i] / np.sqrt(2) for i in range(4)]

# C = sum_{i,k} Tr([A_i, A_k]^dag [A_i, A_k])
C_exact = 0.0
C_components = np.zeros((4, 4))

for i in range(4):
    for k in range(4):
        comm = A[i] @ A[k] - A[k] @ A[i]
        val = np.real(np.trace(comm.conj().T @ comm))
        C_components[i, k] = val
        C_exact += val

print(f"C = sum_{{i,k}} Tr([A_i, A_k]^dag [A_i, A_k])")
print()
print(f"Component matrix C_{{ik}} = Tr([A_i, A_k]^dag [A_i, A_k]):")
for i in range(4):
    print(f"  [{C_components[i,0]:.8f} {C_components[i,1]:.8f} {C_components[i,2]:.8f} {C_components[i,3]:.8f}]")
print()
print(f"C = {C_exact:.10f}")
print()

# Verify: the diagonal is zero ([A_i, A_i] = 0) and it's symmetric
print(f"Diagonal (should be 0): {[f'{C_components[i,i]:.2e}' for i in range(4)]}")
print(f"Symmetric: max|C_ik - C_ki| = {max(abs(C_components[i,k] - C_components[k,i]) for i in range(4) for k in range(4)):.2e}")
print()

# Cross-check with Monte Carlo
print("Cross-checking against Monte Carlo (sigma=0.01)...")
rng = np.random.default_rng(42)
n_mc = 200000
F2_mc = []
for _ in range(n_mc):
    v = rng.standard_normal(4)  # perturbation 1
    w = rng.standard_normal(4)  # perturbation 2
    dM_v = sum(v[i] * basis[i] for i in range(4)) / np.sqrt(2)
    dM_w = sum(w[i] * basis[i] for i in range(4)) / np.sqrt(2)
    a_v = M_inv @ dM_v
    a_w = M_inv @ dM_w
    comm = a_v @ a_w - a_w @ a_v
    F2_mc.append(np.real(np.trace(comm.conj().T @ comm)))

C_mc = np.mean(F2_mc)  # This should equal C (since v,w are unit variance)
print(f"  C (analytic) = {C_exact:.10f}")
print(f"  C (MC, {n_mc} samples) = {C_mc:.10f}")
print(f"  Agreement: {abs(C_exact - C_mc)/C_exact*100:.2f}%")
print()

# ============================================================
# Express C in terms of M* components
# ============================================================

print("="*60)
print("EXPRESSING C IN TERMS OF EQUILIBRIUM")
print("="*60)
print()

# M*^{-1} = adj(M*) / det(M*)
# For M* = (th*I + s.sigma)/sqrt(2):
# det(M*) = (th^2 - s1^2 - s2^2 - s3^2) / 2
# M*^{-1} = (th*I - s.sigma) / (sqrt(2) * det(M*))

th, s1, s2, s3 = m_star[3], m_star[0], m_star[1], m_star[2]
det_M = (th**2 - s1**2 - s2**2 - s3**2) / 2
print(f"det(M*) = (th^2 - |s|^2)/2 = {det_M:.10f}")
print(f"th = {th:.10f}")
print(f"|s|^2 = {s1**2 + s2**2 + s3**2:.10f}")
print(f"th^2 = {th**2:.10f}")
print()

# C depends on M*^{-1}, which is proportional to 1/det(M*).
# C scales as |M*^{-1}|^4 ~ 1/det(M*)^4.
# But det(M*) is determined by the equilibrium, so C is a function of m*.

# Let's compute C / |det(M*)|^{-4} to see if there's a clean coefficient
C_normalized = C_exact * det_M**2  # C ~ |M^{-1}|^4 ~ det^{-4}, so C * det^2 ~ det^{-2}...
# Actually: A_i ~ M^{-1} ~ 1/det, so [A_i, A_k] ~ 1/det^2, |[A,A]|^2 ~ 1/det^4
# So C * det^4 should be a pure geometric number

C_pure = C_exact * abs(det_M)**4
print(f"C * |det(M*)|^4 = {C_pure:.10f}")
print()

# Also try: C * det(M*)^2 (since |[A,A]|^2 involves Tr, which adds a factor)
C_pure2 = C_exact * det_M**2
print(f"C * det(M*)^2 = {C_pure2:.10f}")

# The pure geometric part comes from the Pauli commutators:
# [sigma_i, sigma_j] = 2i eps_{ijk} sigma_k
# [I, anything] = 0
# So only the (i,k) pairs with i,k in {1,2,3} contribute (the A_0 = M^{-1}I/sqrt(2)
# commutes with everything up to the M^{-1} factor)

print()
print("Singleton-only block (i,k in {1,2,3}):")
C_singleton = sum(C_components[i,k] for i in range(1,4) for k in range(1,4))
print(f"  C_singleton = {C_singleton:.10f}")
print(f"  C_full = {C_exact:.10f}")
print(f"  Fraction: {C_singleton/C_exact*100:.2f}%")
print()

# The identity-mixed terms: [A_0, A_k] where A_0 = M^{-1}*I/sqrt(2)
C_identity = sum(C_components[0,k] + C_components[k,0] for k in range(1,4))
print(f"  C_identity_mixed = {C_identity:.10f}")
print(f"  Fraction: {C_identity/C_exact*100:.2f}%")
print()

# ============================================================
# Can we get C in terms of s1, s2, s3, theta analytically?
# ============================================================

print("="*60)
print("SYMBOLIC COMPUTATION")
print("="*60)
print()

# The key: A_i = M^{-1} sigma_i / sqrt(2) for i=1,2,3 and A_0 = M^{-1} I / sqrt(2)
# M^{-1} = (th*I - s.sigma) / (th^2 - |s|^2)  (times sqrt(2))
# Actually: M = (th*I + s.sigma)/sqrt(2), so M^{-1} = sqrt(2) * (th*I - s.sigma) / (th^2 - |s|^2)

# A_i = M^{-1} sigma_i / sqrt(2) = (th*I - s.sigma) sigma_i / (th^2 - |s|^2)
# Using sigma_j sigma_i = delta_{ji} I + i eps_{jik} sigma_k:
# (th*I - s_j sigma_j) sigma_i = th sigma_i - s_j (delta_{ji} I + i eps_{jik} sigma_k)
#                                = -s_i I + th sigma_i - i s_j eps_{jik} sigma_k
#                                = -s_i I + th sigma_i + i (s x e_i) . sigma
# where e_i is the unit vector in direction i.

# So: A_i = [-s_i I + th sigma_i + i (s x e_i).sigma] / D
# where D = th^2 - |s|^2.

# [A_i, A_k] = [th sigma_i + i(s x e_i).sigma, th sigma_k + i(s x e_k).sigma] / D^2
# (the -s_i I and -s_k I terms commute with everything)

# [sigma_i, sigma_k] = 2i eps_{ikl} sigma_l
# [sigma_i, (s x e_k).sigma] = 2i (s x e_k)_l [by the same formula, contracted]
# Actually: [(s x e_k).sigma, sigma_i] = ... this gets messy.
# Let's just verify numerically that C has a clean form.

# Try: C = 8 * |s|^2 / det(M*)^4 ? or something similar
s_sq = s1**2 + s2**2 + s3**2
print(f"|s|^2 = {s_sq:.10f}")
print(f"det(M*) = {det_M:.10f}")
print(f"det(M*)^2 = {det_M**2:.10f}")
print()

# Try various rational combinations
candidates = {
    "8|s|^2 / D^4": 8 * s_sq / det_M**4,
    "4|s|^2 / D^4": 4 * s_sq / det_M**4,
    "16|s|^2 / D^4": 16 * s_sq / det_M**4,
    "4(3|s|^2 + th^2) / D^4": 4 * (3*s_sq + th**2) / det_M**4,
    "8(|s|^2 + th^2) / D^4": 8 * (s_sq + th**2) / det_M**4,
    "4/D^4": 4 / det_M**4,
    "8/D^4": 8 / det_M**4,
    "4|s|^2/(D^2*(th^2+|s|^2))": 4*s_sq / (det_M**2 * (th**2 + s_sq)),
    "8*th^2*|s|^2/D^4": 8 * th**2 * s_sq / det_M**4,
    "4*(th^2+|s|^2)*|s|^2/D^4": 4 * (th**2+s_sq) * s_sq / det_M**4,
}

print(f"C = {C_exact:.10f}")
print()
print("Testing rational forms:")
for name, val in candidates.items():
    ratio = C_exact / val if abs(val) > 1e-30 else float('inf')
    match = abs(ratio - 1.0) < 0.001
    print(f"  {name:40s} = {val:.10f}  ratio = {ratio:.8f}  {'*** MATCH ***' if match else ''}")

print()

# Let me try to compute C symbolically for the SYMMETRIC case s2=s3, general s1, th
# Using the structure of M^{-1}
print("="*60)
print("DIRECT SYMBOLIC EVALUATION")
print("="*60)
print()

# M^{-1} in the Pauli basis:
# M^{-1} = (th*I - s1*sigma_1 - s2*sigma_2 - s3*sigma_3) * sqrt(2) / (th^2 - |s|^2)
# Let D = th^2 - |s|^2 (= 2*det(M))

D = th**2 - s_sq  # = 2 * det_M

# A_i = M^{-1} * basis_i / sqrt(2)
# For basis_0 = I: A_0 = (th*I - s.sigma) / D
# For basis_i (i=1,2,3): A_i = (th*I - s.sigma) sigma_i / D

# [A_i, A_k] for i,k in {1,2,3}:
# = [(th I - s.sigma) sigma_i, (th I - s.sigma) sigma_k] / D^2
# = [(th sigma_i - s.sigma sigma_i), (th sigma_k - s.sigma sigma_k)] / D^2
# Using sigma_j sigma_i = delta_{ji} I + i eps_{jil} sigma_l:
# s.sigma sigma_i = s_i I + i (s x e_i).sigma
# So: (th I - s.sigma) sigma_i = th sigma_i - s_i I - i(s x e_i).sigma
# Let P_i = th sigma_i - s_i I - i(s x e_i).sigma

# Then [A_i, A_k] = [P_i, P_k] / D^2
# [P_i, P_k] = [th sigma_i, th sigma_k] + [th sigma_i, -i(s x e_k).sigma]
#            + [-i(s x e_i).sigma, th sigma_k] + [-i(s x e_i).sigma, -i(s x e_k).sigma]
# Note: [-s_i I, anything] = 0

# [sigma_i, sigma_k] = 2i eps_{ikl} sigma_l
# So [th sigma_i, th sigma_k] = 2i th^2 eps_{ikl} sigma_l

# [sigma_i, (s x e_k).sigma] = sum_m (s x e_k)_m [sigma_i, sigma_m]
#                              = sum_m (s x e_k)_m 2i eps_{iml} sigma_l
#                              = 2i (e_i x (s x e_k)).sigma
# Using BAC-CAB: e_i x (s x e_k) = s(e_i.e_k) - e_k(e_i.s) = s*delta_{ik} - e_k*s_i

# [th sigma_i, -i(s x e_k).sigma] = -i * th * 2i (s delta_{ik} - e_k s_i).sigma
#                                  = 2th (s delta_{ik} - e_k s_i).sigma

# [-i(s x e_i).sigma, th sigma_k] = -[-i(s x e_i).sigma, th sigma_k]... wait, same structure
# = i * th * 2i (s delta_{ik} - e_i s_k).sigma  ... let me just compute numerically.

# Compute [P_i, P_k] directly for all i,k
P = []
for i in range(3):
    ei = np.zeros(3); ei[i] = 1.0
    s_vec = np.array([s1, s2, s3])
    cross_s_ei = np.cross(s_vec, ei)  # s x e_i
    Pi = th * sigma_mat[i] - m_star[i] * I2 - 1j * sum(cross_s_ei[l]*sigma_mat[l] for l in range(3))
    P.append(Pi)

print("Commutator matrix [P_i, P_k]:")
comm_PP = np.zeros((3,3))
for i in range(3):
    for k in range(3):
        comm_ik = P[i] @ P[k] - P[k] @ P[i]
        comm_PP[i,k] = np.real(np.trace(comm_ik.conj().T @ comm_ik))

C_from_P = np.sum(comm_PP) / D**4
print(f"sum Tr([P_i,P_k]^dag [P_i,P_k]) = {np.sum(comm_PP):.10f}")
print(f"D^4 = {D**4:.10f}")
print(f"C_singleton = sum/D^4 = {C_from_P:.10f}")
print(f"C_singleton (direct) = {C_singleton:.10f}")
print(f"Agreement: {abs(C_from_P - C_singleton)/C_singleton*100:.4f}%")
print()

# Now compute the identity-mixed terms
# [A_0, A_k] where A_0 = (th I - s.sigma)/D, A_k = P_k/D
# [A_0, A_k] = [(th I - s.sigma), P_k] / D^2
# = [th I, P_k] - [s.sigma, P_k] / D^2
# = -[s.sigma, P_k] / D^2  (since [I, anything] = 0)

S_mat = sum(m_star[i]*sigma_mat[i] for i in range(3))  # s.sigma

comm_SP = np.zeros(3)
for k in range(3):
    comm_sk = S_mat @ P[k] - P[k] @ S_mat
    comm_SP[k] = np.real(np.trace(comm_sk.conj().T @ comm_sk))

C_identity_from_P = 2 * np.sum(comm_SP) / D**4  # factor 2 because [A_0,A_k] + [A_k,A_0]
print(f"Identity-mixed: 2*sum Tr([S,P_k]^dag [S,P_k]) / D^4 = {C_identity_from_P:.10f}")
print(f"C_identity (direct) = {C_identity:.10f}")
print(f"Agreement: {abs(C_identity_from_P - C_identity)/max(C_identity,1e-30)*100:.4f}%")
print()

C_total_from_P = C_from_P + C_identity_from_P
print(f"C_total (from P decomposition) = {C_total_from_P:.10f}")
print(f"C_total (direct) = {C_exact:.10f}")
print(f"Agreement: {abs(C_total_from_P - C_exact)/C_exact*100:.4f}%")
print()

# ============================================================
# Now try to simplify the numerator
# ============================================================

print("="*60)
print("SIMPLIFICATION")
print("="*60)
print()

numerator = C_exact * D**4
print(f"C * D^4 = {numerator:.10f}")
print(f"C * (2*det(M*))^4 = {C_exact * (2*det_M)**4:.10f}")
print()

# Try: is numerator = f(th, |s|^2, s1) for known framework quantities?
# At the symmetric equilibrium s2=s3, the relevant quantities are:
# th, s1, s2 (=s3), D = th^2 - s1^2 - 2*s2^2

print(f"th = {th:.10f}")
print(f"s1 = {s1:.10f}")
print(f"s2 = s3 = {s2:.10f}")
print(f"D = {D:.10f}")
print(f"|s|^2 = {s_sq:.10f}")
print()

# The Pauli commutator structure:
# [sigma_i, sigma_j] = 2i eps_{ijk} sigma_k
# Tr(sigma_k^dag sigma_l) = 2 delta_{kl}
# So Tr([sigma_i, sigma_j]^dag [sigma_i, sigma_j]) = |2i eps_{ijl}|^2 * 2 = 8 for i!=j, 0 for i=j
# And sum_{i,j} Tr([sigma_i, sigma_j]^dag [sigma_i, sigma_j]) = 8 * 6 = 48

# For the full A_i = M^{-1} basis_i / sqrt(2):
# If M^{-1} = c*I (scalar), then [A_i, A_j] = c^2/2 * [basis_i, basis_j]
# and C would be c^4/4 * 48 = 12 c^4 (only sigma-sigma commutators)

# But M^{-1} is NOT proportional to I in general.
# M^{-1} = sqrt(2) * (th I - s.sigma) / D
# Its "scalar part" is th*sqrt(2)/D and its "vector part" is -s*sqrt(2)/D

# Try a decomposition: M^{-1} = alpha*I + beta*n.sigma
alpha_inv = th * np.sqrt(2) / D
beta_inv = -np.sqrt(2) * np.sqrt(s_sq) / D
print(f"M^{{-1}} decomposition: alpha = {alpha_inv:.10f}, |beta| = {abs(beta_inv):.10f}")
print(f"alpha^2 = {alpha_inv**2:.10f}")
print(f"beta^2 = {beta_inv**2:.10f}")
print(f"alpha^2 + beta^2 = {alpha_inv**2 + beta_inv**2:.10f}")
print(f"2/D^2 = {2/D**2:.10f}")
# alpha^2 + beta^2 = 2*(th^2 + |s|^2)/D^2 = 2*(sum|m|^2_Born_denominator)/D^2

# Actually 2*(th^2+|s|^2)/D^2 = 2*sum_i |m_i|^2 / D^2 where sum|m|^2 is the L2 norm
L2_sq = th**2 + s_sq
print(f"L2^2 = th^2 + |s|^2 = {L2_sq:.10f}")
print(f"2*L2^2/D^2 = {2*L2_sq/D**2:.10f}")
print()

# For M^{-1} = alpha*I + beta*n.sigma, with n = -s/|s|:
# The commutator [M^{-1} sigma_i, M^{-1} sigma_j] involves cross terms
# This is getting algebraically heavy. Let me just check:
# is C * D^4 a polynomial in th, s1, s2?

# Evaluate at several equilibria by varying the evidence slightly
print("Testing whether C*D^4 is universal (independent of equilibrium):")
print()

# At the actual equilibrium:
print(f"  Actual equilibrium: C*D^4 = {numerator:.10f}")

# At the symmetric point (all s_i equal):
m_sym = np.array([0.25, 0.25, 0.25, 0.25])
m_sym = enforce_floor(m_sym)
M_sym = to_matrix(m_sym)
Mi_sym = inv(M_sym)
D_sym = m_sym[3]**2 - np.sum(m_sym[:3]**2)
A_sym = [Mi_sym @ basis[i] / np.sqrt(2) for i in range(4)]
C_sym = sum(np.real(np.trace((A_sym[i]@A_sym[k]-A_sym[k]@A_sym[i]).conj().T @ (A_sym[i]@A_sym[k]-A_sym[k]@A_sym[i]))) for i in range(4) for k in range(4))
print(f"  Symmetric (0.25,0.25,0.25,0.25 floored): C*D^4 = {C_sym * D_sym**4:.10f}")

# At a different asymmetric point
m_test = np.array([0.5, 0.15, 0.15, 0.2])
m_test = enforce_floor(m_test)
M_test = to_matrix(m_test)
Mi_test = inv(M_test)
D_test = m_test[3]**2 - np.sum(m_test[:3]**2)
A_test = [Mi_test @ basis[i] / np.sqrt(2) for i in range(4)]
C_test = sum(np.real(np.trace((A_test[i]@A_test[k]-A_test[k]@A_test[i]).conj().T @ (A_test[i]@A_test[k]-A_test[k]@A_test[i]))) for i in range(4) for k in range(4))
print(f"  Asymmetric (0.5,0.15,0.15,0.2 floored): C*D^4 = {C_test * D_test**4:.10f}")

# Another test
m_test2 = np.array([0.6, 0.1, 0.1, 0.2])
m_test2 = enforce_floor(m_test2)
M_test2 = to_matrix(m_test2)
Mi_test2 = inv(M_test2)
D_test2 = m_test2[3]**2 - np.sum(m_test2[:3]**2)
A_test2 = [Mi_test2 @ basis[i] / np.sqrt(2) for i in range(4)]
C_test2 = sum(np.real(np.trace((A_test2[i]@A_test2[k]-A_test2[k]@A_test2[i]).conj().T @ (A_test2[i]@A_test2[k]-A_test2[k]@A_test2[i]))) for i in range(4) for k in range(4))
print(f"  Another (0.6,0.1,0.1,0.2 floored): C*D^4 = {C_test2 * D_test2**4:.10f}")

print()

# Check if C*D^4 = 8*|s|^2 or some other simple form at each point
for label, m_pt in [("equilibrium", m_star), ("symmetric", m_sym), ("asymmetric", m_test), ("another", m_test2)]:
    th_pt = m_pt[3]; s_sq_pt = np.sum(m_pt[:3]**2)
    D_pt = th_pt**2 - s_sq_pt
    M_pt = to_matrix(m_pt); Mi_pt = inv(M_pt)
    A_pt = [Mi_pt @ basis[i] / np.sqrt(2) for i in range(4)]
    C_pt = sum(np.real(np.trace((A_pt[i]@A_pt[k]-A_pt[k]@A_pt[i]).conj().T @ (A_pt[i]@A_pt[k]-A_pt[k]@A_pt[i]))) for i in range(4) for k in range(4))
    CD4 = C_pt * D_pt**4
    print(f"  {label:15s}: C*D^4 = {CD4:.6f}, 8|s|^2 = {8*s_sq_pt:.6f}, ratio = {CD4/(8*s_sq_pt) if s_sq_pt > 1e-10 else 0:.6f}")

print()
print("DONE.")
