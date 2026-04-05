"""
Ward reconstruction: the Penrose contour integral.

The physical curvature F+ on spacetime is obtained by integrating the
non-integrability data over each twistor line L_x ~ CP^1.

At UNIFORM equilibrium: F = 0 (no spatial variation).
For FLUCTUATIONS: F != 0. The curvature comes from how the floor correction
varies spatially, translated to zeta-dependence through the incidence relation.

The incidence relation:  omega^A = i x^{AA'} pi_{A'}
In affine coords (zeta = pi_1'/pi_0'):
  z1_bar = -i/sqrt(2) [(x0+x3) + (x1-ix2) zeta_bar]
  z2_bar = -i/sqrt(2) [(x1+ix2) + (x0-x3) zeta_bar]

On the equator of CP^1 (zeta = e^{ialpha}), zeta_bar = 1/zeta.

The floor correction's spatial derivative, projected through the incidence
relation, becomes a Laurent polynomial in zeta on the equator. The Penrose
integral picks up specific Laurent coefficients = F+ components.
"""

import numpy as np
from numpy.linalg import inv, norm, det
from scipy.optimize import brentq
import sys

# ============================================================
# Setup (same as before)
# ============================================================
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
M_star = to_matrix(m_star)
m_ds, _ = ds_combine(m_star, e_star)
M_ds = to_matrix(m_ds)
M_ds_inv = inv(M_ds)
epsilon = M_star - M_ds

print("="*70)
print("PENROSE CONTOUR INTEGRAL FOR F+")
print("="*70)
print(f"K* = {K_star:.10f}")
print(f"||epsilon|| = {norm(epsilon):.6f}")
sys.stdout.flush()

# ============================================================
# Compute the floor correction Jacobian (same as v1/v2)
# ============================================================
s1, s2, s3, th = m_star
e1, e2, e3, ph = e_star
oneMinusK = 1 - K_star

S_vals = np.array([s1*(e1+ph)+th*e1, s2*(e2+ph)+th*e2, s3*(e3+ph)+th*e3, th*ph])
N = S_vals / oneMinusK

dSdm = np.zeros((4,4))
for i in range(3):
    dSdm[i,i] = e_star[i] + ph; dSdm[i,3] = e_star[i]
dSdm[3,3] = ph
dKdm = np.array([e2+e3, e1+e3, e1+e2, 0.0])

J_DS = np.zeros((4,4))
for i in range(4):
    for j in range(4):
        J_DS[i,j] = (dSdm[i,j] * oneMinusK + S_vals[i] * dKdm[j]) / oneMinusK**2

Ss = N[0]+N[1]+N[2]; Sq = N[0]**2+N[1]**2+N[2]**2; R = Sq/Ss**2
u = np.sqrt(26*R); t = (u-R)/(26-R); g = (1-t)/Ss
dtdR = (338+13*R-26*u)/(u*(26-R)**2)
dRdN = np.array([2*(N[i]*Ss-Sq)/Ss**3 for i in range(3)])
dtdN = np.array([dtdR*dRdN[i] for i in range(3)])
dgdN = np.array([(-dtdN[i]*Ss-(1-t))/Ss**2 for i in range(3)])

J_floor = np.zeros((4,4))
for i in range(3):
    for j in range(3): J_floor[i,j] = (g if i==j else 0.0) + N[i]*dgdN[j]
J_floor[3,:3] = dtdN; J_floor[3,3] = 0.0

J_correction = (J_floor - np.eye(4)) @ J_DS

# ============================================================
# The Ward connection basis vectors B_mu
# B_mu = M_ds^{-1} . sum_j (basis_j/sqrt(2)) . J_correction_{j,mu}
# ============================================================
B = np.zeros((4, 2, 2), dtype=complex)
for mu in range(4):
    for j in range(4):
        B[mu] += M_ds_inv @ (basis[j] / np.sqrt(2)) * J_correction[j, mu]

# ============================================================
# INCIDENCE RELATION
# ============================================================
# At x0 = 0, the twistor line L_0 is parameterized by zeta.
# The anti-holomorphic derivatives dz_bar/dx^mu give the
# zeta-dependent coefficients:
#
# dz1_bar/dx^mu = -i/sqrt(2) * [delta_{mu,0} + delta_{mu,3}
#                                + (delta_{mu,1} - i*delta_{mu,2})*zeta_bar]
#
# dz2_bar/dx^mu = -i/sqrt(2) * [(delta_{mu,1} + i*delta_{mu,2})
#                                + (delta_{mu,0} - delta_{mu,3})*zeta_bar]
#
# On the equator: zeta = e^{ialpha}, zeta_bar = e^{-ialpha} = 1/zeta.
#
# So dz_bar/dx^mu = A_mu + C_mu / zeta  (Laurent polynomial in zeta)

# Coefficients A_mu (zeta^0 term) and C_mu (zeta^{-1} term)
# for dz1_bar/dx^mu:
A1 = np.zeros(4, dtype=complex)  # constant part of dz1_bar/dx^mu
C1 = np.zeros(4, dtype=complex)  # 1/zeta part of dz1_bar/dx^mu

A1[0] = -1j/np.sqrt(2)  # mu=0: -i/sqrt(2) * 1
A1[1] = 0                # mu=1: -i/sqrt(2) * 0 (zeta_bar term)
A1[2] = 0                # mu=2: similar
A1[3] = -1j/np.sqrt(2)  # mu=3: -i/sqrt(2) * 1

C1[0] = 0                # mu=0: no zeta_bar term
C1[1] = -1j/np.sqrt(2)  # mu=1: -i/sqrt(2) * 1 * zeta_bar = ... * 1/zeta
C1[2] = 1/np.sqrt(2)     # mu=2: -i/sqrt(2) * (-i) * zeta_bar
C1[3] = 0                # mu=3: no zeta_bar term

# for dz2_bar/dx^mu:
A2 = np.zeros(4, dtype=complex)
C2 = np.zeros(4, dtype=complex)

A2[0] = 0
A2[1] = -1j/np.sqrt(2)  # mu=1: -i/sqrt(2) * 1
A2[2] = 1/np.sqrt(2)     # mu=2: -i/sqrt(2) * (i)
A2[3] = 0

C2[0] = -1j/np.sqrt(2)  # mu=0: -i/sqrt(2) * 1 * zeta_bar
C2[1] = 0
C2[2] = 0
C2[3] = 1j/np.sqrt(2)   # mu=3: -i/sqrt(2) * (-1) * zeta_bar

print("\n" + "="*70)
print("INCIDENCE RELATION COEFFICIENTS")
print("="*70)
for mu in range(4):
    print(f"  mu={mu}: dz1_bar/dx = {A1[mu]:.4f} + {C1[mu]:.4f}/zeta")
    print(f"         dz2_bar/dx = {A2[mu]:.4f} + {C2[mu]:.4f}/zeta")

# ============================================================
# The integrand on L_0
# ============================================================
# The floor correction derivative in spacetime direction mu,
# expressed on L_0, is:
#
#   d(epsilon)/dz_bar_i = sum_mu (d(epsilon)/dx^mu) * (dx^mu/dz_bar_i)
#
# where d(epsilon)/dx^mu = sum_j (basis_j/sqrt(2)) * J_correction_{j,k} * (dm_k/dx^mu)
#
# For a Gaussian perturbation in the mu-th direction:
#   dm_k/dx^mu = delta_{k,mu_mapped} (perturbation in one mass component)
#
# But actually, for a general spatial variation of m:
#   d(epsilon)/dx^mu involves J_correction applied to dm/dx^mu.
#
# The F+ Penrose integral involves:
#   F+_{A'B'} = oint pi_{A'} pi_{B'} * M_ds^{-1} * N_J * (d_epsilon/dx) * dz/dz_bar * dzeta
#
# At leading order, the integrand on L_0 is a Laurent polynomial in zeta.
# The Penrose integral picks up the residue at zeta = 0.

# For each spacetime direction mu, the floor correction variation is B[mu].
# Through the incidence relation, B[mu] gets ζ-dependent coefficients:
#
# Contribution to integrand from direction mu:
#   B[mu] * (A1[mu] + C1[mu]/zeta) for the z1_bar channel
#   B[mu] * (A2[mu] + C2[mu]/zeta) for the z2_bar channel
#
# The total integrand for the spin-1 Penrose transform:
#   rho(zeta) = sum_mu B[mu] * [(A1[mu] + C1[mu]/zeta) + (A2[mu] + C2[mu]/zeta)]
#             = sum_mu B[mu] * [(A1[mu]+A2[mu]) + (C1[mu]+C2[mu])/zeta]

# Laurent expansion: rho(zeta) = rho_0 + rho_{-1}/zeta
rho_0 = np.zeros((2,2), dtype=complex)  # coefficient of zeta^0
rho_m1 = np.zeros((2,2), dtype=complex)  # coefficient of zeta^{-1}

for mu in range(4):
    rho_0 += B[mu] * (A1[mu] + A2[mu])
    rho_m1 += B[mu] * (C1[mu] + C2[mu])

print("\n" + "="*70)
print("LAURENT EXPANSION OF INTEGRAND ON L_0")
print("="*70)
print(f"\nrho_0 (constant term):")
print(f"  {rho_0}")
print(f"  ||rho_0|| = {norm(rho_0):.8f}")

print(f"\nrho_{{-1}} (1/zeta term):")
print(f"  {rho_m1}")
print(f"  ||rho_{{-1}}|| = {norm(rho_m1):.8f}")

# ============================================================
# PENROSE INTEGRAL: F+ components
# ============================================================
# F+_{A'B'}(0) = oint pi_{A'} pi_{B'} * rho(zeta) * dzeta / (2pi i)
#
# pi_{0'} = 1, pi_{1'} = zeta (in affine coords)
#
# F+_{0'0'} = oint 1 * rho(zeta) * dzeta/(2pi i)
#           = Res_{zeta=0} [rho_0 + rho_{-1}/zeta]
#           = rho_{-1}  (the 1/zeta coefficient)
#
# F+_{0'1'} = oint zeta * rho(zeta) * dzeta/(2pi i)
#           = Res_{zeta=0} [zeta * rho_0 + rho_{-1}]
#           = rho_{-1} ... wait, Res of zeta*rho_0 is zero (no 1/zeta)
#           Actually: zeta * [rho_0 + rho_{-1}/zeta] = rho_0*zeta + rho_{-1}
#           Res at zeta=0: 0 (no 1/zeta term)
#           Hmm, so F+_{0'1'} = 0?
#
# Wait, I need to be more careful. The contour integral on CP^1:
# oint f(zeta) dzeta/(2pi i) = sum of residues inside the contour.
# On the equator |zeta|=1, the contour encloses zeta=0.
# Residue of f at zeta=0 = coefficient of 1/zeta in Laurent expansion.
#
# For f = rho_0 + rho_{-1}/zeta:
#   Res_{zeta=0} f = rho_{-1}
#
# For f = zeta * (rho_0 + rho_{-1}/zeta) = rho_0*zeta + rho_{-1}:
#   Res_{zeta=0} (rho_0*zeta + rho_{-1}) = 0 (no 1/zeta term)
#
# For f = zeta^2 * (rho_0 + rho_{-1}/zeta) = rho_0*zeta^2 + rho_{-1}*zeta:
#   Res_{zeta=0} = 0

print("\n" + "="*70)
print("PENROSE INTEGRAL: F+ COMPONENTS")
print("="*70)

# F+_{A'B'} = residue of pi_{A'} pi_{B'} * rho(zeta) at zeta=0
F_plus = np.zeros((2, 2, 2, 2), dtype=complex)  # F+_{A'B'} is a 2x2 matrix for each (A',B')

# (0',0'): pi^2 = 1 -> residue of rho = rho_{-1}
F_plus_00 = rho_m1

# (0',1'): pi_{0'} pi_{1'} = zeta -> residue of zeta * rho = 0
F_plus_01 = np.zeros((2,2), dtype=complex)

# (1',1'): pi_{1'}^2 = zeta^2 -> residue of zeta^2 * rho = 0
F_plus_11 = np.zeros((2,2), dtype=complex)

print(f"\nF+_{{0'0'}} = rho_{{-1}}:")
print(f"  {F_plus_00}")
print(f"  ||F+_{{0'0'}}|| = {norm(F_plus_00):.8f}")

print(f"\nF+_{{0'1'}} = 0 (by residue calculation)")
print(f"F+_{{1'1'}} = 0 (by residue calculation)")

# Total |F+|^2 = |F+_{0'0'}|^2 + 2|F+_{0'1'}|^2 + |F+_{1'1'}|^2
F_plus_sq = np.real(np.trace(F_plus_00.conj().T @ F_plus_00))
print(f"\n|F+|^2 = Tr(F+_{{0'0'}}^dag F+_{{0'0'}}) = {F_plus_sq:.10f}")

# su(2) decomposition of F+
tr_F = np.trace(F_plus_00) / 2
su2_F = np.array([np.trace(sigma[i] @ F_plus_00) / 2 for i in range(3)])
print(f"\nsu(2) decomposition of F+:")
print(f"  scalar (tr/2): {tr_F:.8f}")
print(f"  su(2) components: {su2_F}")
na_frac = np.sum(np.abs(su2_F)**2) / (abs(tr_F)**2 + np.sum(np.abs(su2_F)**2)) if abs(tr_F)**2 + np.sum(np.abs(su2_F)**2) > 1e-20 else 0
print(f"  Non-abelian fraction: {na_frac:.4f}")

# ============================================================
# NUMERICAL VERIFICATION via direct quadrature
# ============================================================
print("\n" + "="*70)
print("NUMERICAL VERIFICATION: direct quadrature on CP^1")
print("="*70)

N_quad = 1000
F_quad_00 = np.zeros((2,2), dtype=complex)
F_quad_01 = np.zeros((2,2), dtype=complex)
F_quad_11 = np.zeros((2,2), dtype=complex)

for k in range(N_quad):
    alpha = 2 * np.pi * k / N_quad
    zeta = np.exp(1j * alpha)
    zeta_bar = 1.0 / zeta  # on the equator

    # Integrand: rho(zeta) = rho_0 + rho_{-1}/zeta
    rho = rho_0 + rho_m1 / zeta

    # pi_{0'} = 1, pi_{1'} = zeta
    dzeta = 1j * zeta * (2*np.pi/N_quad)  # dzeta = i*zeta*dalpha

    F_quad_00 += 1.0 * 1.0 * rho * dzeta / (2*np.pi*1j)
    F_quad_01 += 1.0 * zeta * rho * dzeta / (2*np.pi*1j)
    F_quad_11 += zeta * zeta * rho * dzeta / (2*np.pi*1j)

print(f"\nQuadrature F+_{{0'0'}}:")
print(f"  {F_quad_00}")
print(f"  ||F+_{{0'0'}}|| = {norm(F_quad_00):.8f}")
print(f"  Agrees with analytic: {norm(F_quad_00 - F_plus_00) < 1e-6}")

print(f"\nQuadrature F+_{{0'1'}}:")
print(f"  {F_quad_01}")
print(f"  ||F+_{{0'1'}}|| = {norm(F_quad_01):.8f}")

print(f"\nQuadrature F+_{{1'1'}}:")
print(f"  {F_quad_11}")
print(f"  ||F+_{{1'1'}}|| = {norm(F_quad_11):.8f}")

# ============================================================
# The F+ coefficient per unit fluctuation
# ============================================================
print("\n" + "="*70)
print("F+ COEFFICIENT (Ward-reconstructed, per unit fluctuation)")
print("="*70)

# |F+|^2 per unit (dm/dx)^2 fluctuation:
# The B[mu] are response of the Ward connection to dm in direction mu.
# Through the incidence relation, only the C (1/zeta) part contributes.
# F+_{0'0'} = sum_mu B[mu] * (C1[mu] + C2[mu])

F_plus_ward = F_plus_sq
print(f"\n|F+_Ward|^2 (per unit spatial derivative variance) = {F_plus_ward:.10f}")

# Compare with the Maurer-Cartan commutator content
C_mc = 149.596
print(f"|[a,a]|^2 (MC commutator content) = {C_mc:.3f}")
print(f"Ratio F+_Ward / MC_comm = {F_plus_ward / C_mc:.8f}")

# Check algebraic form
print(f"\nAlgebraic checks:")
print(f"  F+_Ward * det(M_ds)^2 = {F_plus_ward * abs(det(M_ds))**2:.10f}")
print(f"  F+_Ward * det(M*)^2 = {F_plus_ward * abs(det(M_star))**2:.10f}")
print(f"  F+_Ward * K* = {F_plus_ward * 7/30:.10f}")
print(f"  F+_Ward / K* = {F_plus_ward / (7/30):.10f}")
print(f"  F+_Ward / (1-K*) = {F_plus_ward / (23/30):.10f}")

# Per unit sigma^4 (fluctuation variance from transfer operator):
sigma_eff_sq = 1.0 / (1 - 0.2829**2)
sigma_eff_4 = sigma_eff_sq**2
print(f"\nsigma_eff^2 = 1/(1-lambda_0^2) = {sigma_eff_sq:.6f}")
print(f"<|F+|^2>_Ward = F+_coeff * sigma_eff^4 = {F_plus_ward * sigma_eff_4:.6f}")
print(f"<|[a,a]|^2>_MC = C * sigma_eff^4 = {C_mc * sigma_eff_4:.6f}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
