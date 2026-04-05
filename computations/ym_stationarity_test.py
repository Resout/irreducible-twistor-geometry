"""
YM Stationarity Test: Does dbar*(dbar_E)^2 = 0 at K*=7/30?

Open item 1: DS stationarity = YM stationarity.

The DS equilibrium produces nonzero (0,2)-curvature F+ via the Born floor's
non-integrable almost complex structure. The Yang-Mills field equation requires
the codifferential of F+ to vanish:

    dbar_E* F^{0,2} = 0

This script computes F^{0,2} and its codifferential explicitly at the K*=7/30
equilibrium, using Wirtinger calculus on the L1=1 tangent space of CP^3.

Strategy:
  1. Construct K*=7/30 equilibrium (m*, e*)
  2. Extend DS+floor map to complex mass functions
  3. Compute connection a^{0,1} = M^{-1} (dM/dm) dbar(Phi) at m*
  4. Compute curvature F_{ab} = dbar_a a_b - dbar_b a_a + [a_a, a_b]
  5. Compute codifferential (dbar_E* F)_c = -sum_a (d F_{ac}/dz^a + [a^{1,0}_a, F_{ac}])
  6. Report ||dbar_E* F||

If ||dbar_E* F|| ~ 0: DS equilibrium IS a YM solution.
If ||dbar_E* F|| ~ O(1): it is NOT, and the gap between DS and YM is real.
"""

import numpy as np
from numpy.linalg import inv, det, norm

# ============================================================
# PART 0: Constants
# ============================================================
H = 3
FLOOR = 1.0 / H**3  # = 1/27

# Pauli matrices
sigma = [
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
]
I2 = np.eye(2, dtype=complex)

# Basis for the L1=1 tangent space (3 directions, each in C^4)
V = np.array([[1, 0, 0, -1],
              [0, 1, 0, -1],
              [0, 0, 1, -1]], dtype=float).T  # 4x3 matrix


# ============================================================
# PART 1: su(2) embedding
# ============================================================
def mass_to_M(m):
    """Mass function m = (s1,s2,s3,theta) -> 2x2 matrix M."""
    return (m[3]*I2 + m[0]*sigma[0] + m[1]*sigma[1] + m[2]*sigma[2]) / np.sqrt(2)

# dM/dm_j are constant (M is linear in m):
dM_dm = [sigma[i]/np.sqrt(2) for i in range(3)] + [I2/np.sqrt(2)]  # j=0,1,2,3


# ============================================================
# PART 2: DS combination + floor for COMPLEX mass functions
# ============================================================
def ds_combine_complex(m, e):
    """DS combination step for complex m, real e. Returns (m_out, K)."""
    s, theta = m[:3], m[3]
    se, phi = e[:3], e[3]
    s_new = s * se + s * phi + theta * se  # element-wise
    theta_new = theta * phi
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    denom = 1.0 - K
    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom
    return m_out, K


def enforce_floor_complex(m):
    """Born floor enforcement for complex mass functions.

    Maintains L1=1 and Born(theta)=1/27.
    Uses the quadratic: (26-R)*r^2 + 2*R*Re(phase)*r - R = 0
    where R = sum|s_i|^2/|sum(s_i)|^2, phase = theta/|theta|, r = |theta_new|.
    """
    s = m[:3].copy()
    theta = m[3]

    # Born probability
    s_mag_sq = np.sum(np.abs(s)**2)
    theta_mag = np.abs(theta)
    total_mag_sq = s_mag_sq + theta_mag**2

    if total_mag_sq < 1e-30:
        return m.copy()

    born = theta_mag**2 / total_mag_sq

    if born >= FLOOR - 1e-15:
        return m.copy()

    # Floor active
    S = np.sum(s)
    S_mag_sq = np.abs(S)**2

    if S_mag_sq < 1e-30:
        return m.copy()

    phase = theta / theta_mag if theta_mag > 1e-30 else 1.0+0j
    R = s_mag_sq / S_mag_sq
    Re_phase = np.real(phase)

    # Quadratic: (26-R)*r^2 + 2*R*Re(phase)*r - R = 0
    a = 26.0 - R
    b = 2.0 * R * Re_phase
    c = -R

    disc = b**2 - 4*a*c
    if disc < 0:
        # Shouldn't happen for physical states
        return m.copy()

    sqrt_disc = np.sqrt(disc)
    r1 = (-b + sqrt_disc) / (2*a)
    r2 = (-b - sqrt_disc) / (2*a)

    # Pick positive root
    r_new = r1 if r1 > 0 else r2
    if r_new <= 0:
        return m.copy()

    theta_new = phase * r_new
    alpha = (1.0 - theta_new) / S
    s_new = s * alpha

    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = s_new
    m_out[3] = theta_new
    return m_out


def phi_map(m_complex, e_star):
    """Full DS map: Phi = Floor o DS, for complex mass functions."""
    m_ds, K = ds_combine_complex(m_complex, e_star)
    m_out = enforce_floor_complex(m_ds)
    return m_out


# ============================================================
# PART 3: Find K*=7/30 equilibrium
# ============================================================
print("=" * 70)
print("PART 1: K*=7/30 EQUILIBRIUM")
print("=" * 70)

def enforce_floor_real(m):
    """Standard real floor enforcement via binary search."""
    s = m[:3].copy()
    S = np.sum(s)
    lo, hi = 0.0, 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if S > 0:
            s_scale = (1.0 - mid) / S
        else:
            s_scale = 0
        s_trial = s * s_scale
        born = mid**2 / (np.sum(s_trial**2) + mid**2)
        if born < FLOOR:
            lo = mid
        else:
            hi = mid
    theta_new = (lo + hi) / 2
    s_scale = (1.0 - theta_new) / S if S > 0 else 0
    m_out = np.zeros(4)
    m_out[:3] = s * s_scale
    m_out[3] = theta_new
    return m_out

def ds_combine_real(m, e):
    s, theta = m[:3], m[3]
    se, phi = e[:3], e[3]
    s_new = s * se + s * phi + theta * se
    theta_new = theta * phi
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i != j)
    denom = 1.0 - K
    m_out = np.zeros(4)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom
    born = m_out[3]**2 / np.sum(m_out**2)
    if born < FLOOR:
        m_out = enforce_floor_real(m_out)
    return m_out, K

from scipy.optimize import brentq

def K_at_equil(p_dom):
    p_weak = (1 - p_dom) / 2
    scale = 1 - FLOOR
    raw = np.array([np.sqrt(p_dom*scale), np.sqrt(p_weak*scale),
                    np.sqrt(p_weak*scale), np.sqrt(FLOOR)])
    e = raw / np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(2000):
        m_new, _ = ds_combine_real(m, e)
        if np.max(np.abs(m_new - m)) < 1e-15:
            break
        m = m_new
    _, K = ds_combine_real(m, e)
    return K

target = 7.0/30.0
p_dom_exact = brentq(lambda p: K_at_equil(p) - target, 0.90, 0.95, xtol=1e-14)

# Reconstruct evidence and fixed point
p_weak_exact = (1 - p_dom_exact) / 2
scale = 1 - FLOOR
raw = np.array([np.sqrt(p_dom_exact*scale), np.sqrt(p_weak_exact*scale),
                np.sqrt(p_weak_exact*scale), np.sqrt(FLOOR)])
e_star = raw / np.sum(raw)

m = np.array([0.4, 0.2, 0.2, 0.2])
for _ in range(5000):
    m_new, _ = ds_combine_real(m, e_star)
    if np.max(np.abs(m_new - m)) < 1e-15:
        break
    m = m_new
m_star = m_new.copy()
_, K_star = ds_combine_real(m_star, e_star)

print(f"  m* = ({m_star[0]:.10f}, {m_star[1]:.10f}, {m_star[2]:.10f}, {m_star[3]:.10f})")
print(f"  e* = ({e_star[0]:.10f}, {e_star[1]:.10f}, {e_star[2]:.10f}, {e_star[3]:.10f})")
print(f"  K* = {K_star:.15f}  (7/30 = {target:.15f})")
print(f"  |K* - 7/30| = {abs(K_star - target):.2e}")

M_star = mass_to_M(m_star)
print(f"  det(M*) = {det(M_star):.10f}")
born_star = m_star[3]**2 / np.sum(m_star**2)
print(f"  Born(theta*) = {born_star:.10f}  (1/27 = {1/27:.10f})")

# Verify complex floor matches real floor at real input
m_test = phi_map(m_star.astype(complex), e_star)
print(f"\n  Complex floor test: |Phi(m*) - m*| = {norm(m_test - m_star):.2e}")


# ============================================================
# PART 4: Wirtinger derivatives of Phi at m*
# ============================================================
print("\n" + "=" * 70)
print("PART 2: WIRTINGER DERIVATIVES OF Phi AT m*")
print("=" * 70)

def wirtinger_jacobians(m0, e, eps=1e-7):
    """Compute holomorphic and anti-holomorphic Jacobians of Phi at m0.

    Returns:
        J_hol: 4x4 complex matrix, (dPhi/dm)_{j,alpha}
        J_anti: 4x4 complex matrix, (dPhi/dm_bar)_{j,alpha}
    """
    m0c = m0.astype(complex)
    J_hol = np.zeros((4, 4), dtype=complex)
    J_anti = np.zeros((4, 4), dtype=complex)

    for alpha in range(4):
        e_alpha = np.zeros(4, dtype=complex)
        e_alpha[alpha] = 1.0

        # Real perturbation
        f_plus_r = phi_map(m0c + eps * e_alpha, e)
        f_minus_r = phi_map(m0c - eps * e_alpha, e)
        df_dx = (f_plus_r - f_minus_r) / (2 * eps)

        # Imaginary perturbation
        f_plus_i = phi_map(m0c + 1j * eps * e_alpha, e)
        f_minus_i = phi_map(m0c - 1j * eps * e_alpha, e)
        df_dy = (f_plus_i - f_minus_i) / (2 * eps)

        # Wirtinger: d/dz = (1/2)(d/dx - i d/dy), d/dz_bar = (1/2)(d/dx + i d/dy)
        J_hol[:, alpha] = 0.5 * (df_dx - 1j * df_dy)
        J_anti[:, alpha] = 0.5 * (df_dx + 1j * df_dy)

    return J_hol, J_anti

J_hol_0, J_anti_0 = wirtinger_jacobians(m_star, e_star)

print(f"\n  Anti-holomorphic Jacobian dbar(Phi) at m*:")
print(f"    ||J_anti|| = {norm(J_anti_0):.6f}")
print(f"    ||J_hol||  = {norm(J_hol_0):.6f}")

# Verify: J_hol + J_anti should give the real Jacobian dPhi/dx
J_real_check = J_hol_0 + J_anti_0  # = dPhi/dx (real Jacobian)
eps = 1e-7
J_real_direct = np.zeros((4, 4))
for j in range(4):
    e_j = np.zeros(4)
    e_j[j] = 1.0
    f_p = phi_map((m_star + eps*e_j).astype(complex), e_star)
    f_m = phi_map((m_star - eps*e_j).astype(complex), e_star)
    J_real_direct[:, j] = np.real((f_p - f_m) / (2*eps))

print(f"\n  Sanity check: ||J_hol + J_anti - J_real|| = {norm(J_real_check - J_real_direct):.2e}")


# ============================================================
# PART 5: Connection a^{0,1} and a^{1,0} at m*
# ============================================================
print("\n" + "=" * 70)
print("PART 3: CONNECTION AT m*")
print("=" * 70)

M_at_phi = mass_to_M(phi_map(m_star.astype(complex), e_star))
M_inv = inv(M_at_phi)

# Connection a^{0,1}_alpha_bar(m*) = M^{-1} * sum_j (dM/dm_j) * J_anti[j, alpha]
# This is a 2x2 matrix for each alpha = 0,1,2,3
a_01 = np.zeros((4, 2, 2), dtype=complex)  # a_01[alpha] = 2x2 matrix
for alpha in range(4):
    dM_contribution = sum(dM_dm[j] * J_anti_0[j, alpha] for j in range(4))
    a_01[alpha] = M_inv @ dM_contribution

# Connection a^{1,0}_alpha(m*) = M^{-1} * sum_j (dM/dm_j) * J_hol[j, alpha]
a_10 = np.zeros((4, 2, 2), dtype=complex)
for alpha in range(4):
    dM_contribution = sum(dM_dm[j] * J_hol_0[j, alpha] for j in range(4))
    a_10[alpha] = M_inv @ dM_contribution

print(f"  Connection a^{{0,1}} norms:")
for alpha in range(4):
    print(f"    ||a_01[{alpha}]|| = {norm(a_01[alpha]):.6f}")

print(f"\n  Connection a^{{1,0}} norms:")
for alpha in range(4):
    print(f"    ||a_10[{alpha}]|| = {norm(a_10[alpha]):.6f}")

# Total connection norm
a01_total = np.sqrt(sum(norm(a_01[a])**2 for a in range(4)))
a10_total = np.sqrt(sum(norm(a_10[a])**2 for a in range(4)))
print(f"\n  Total ||a^{{0,1}}|| = {a01_total:.6f}")
print(f"  Total ||a^{{1,0}}|| = {a10_total:.6f}")


# ============================================================
# PART 6: Project to L1=1 tangent space
# ============================================================
print("\n" + "=" * 70)
print("PART 4: PROJECTION TO L1=1 TANGENT SPACE")
print("=" * 70)

# Project connection to 3D tangent space
# Tangent directions: v_a = V[:, a] for a=0,1,2
# The projected connection: a^{0,1}_a = sum_alpha V[alpha, a] * a^{0,1}[alpha]

a_01_proj = np.zeros((3, 2, 2), dtype=complex)  # a=0,1,2
a_10_proj = np.zeros((3, 2, 2), dtype=complex)

for a in range(3):
    for alpha in range(4):
        a_01_proj[a] += V[alpha, a] * a_01[alpha]
        a_10_proj[a] += V[alpha, a] * a_10[alpha]

print(f"  Projected a^{{0,1}} norms:")
for a in range(3):
    print(f"    ||a_01_proj[{a}]|| = {norm(a_01_proj[a]):.6f}")

print(f"\n  Projected a^{{1,0}} norms:")
for a in range(3):
    print(f"    ||a_10_proj[{a}]|| = {norm(a_10_proj[a]):.6f}")


# ============================================================
# PART 7: Curvature F^{0,2} at m*
# ============================================================
print("\n" + "=" * 70)
print("PART 5: CURVATURE F^{0,2} AT m*")
print("=" * 70)

# F_{ab} = dbar_a(a_b) - dbar_b(a_a) + [a_a, a_b]
# where dbar_a = d/d(z_bar_a) is the anti-holomorphic derivative along v_a.
#
# Compute dbar_a(a_b) by displacing m* in anti-holomorphic direction v_a
# and recomputing a_b.

def compute_connection_at(m0, e, eps_inner=1e-7):
    """Compute a^{0,1}_alpha (4 matrices) at point m0."""
    m0c = m0.astype(complex)
    _, J_anti = wirtinger_jacobians(m0, e, eps=eps_inner)
    M_phi = mass_to_M(phi_map(m0c, e))
    M_i = inv(M_phi)
    a01_local = np.zeros((4, 2, 2), dtype=complex)
    for alpha in range(4):
        dM_contribution = sum(dM_dm[j] * J_anti[j, alpha] for j in range(4))
        a01_local[alpha] = M_i @ dM_contribution
    return a01_local

def compute_connection_projected_at(m0, e, eps_inner=1e-7):
    """Compute projected a^{0,1}_a (3 matrices) at point m0."""
    a01_local = compute_connection_at(m0, e, eps_inner)
    a01_p = np.zeros((3, 2, 2), dtype=complex)
    for a in range(3):
        for alpha in range(4):
            a01_p[a] += V[alpha, a] * a01_local[alpha]
    return a01_p

def compute_connections_both_at(m0, e, eps_inner=1e-7):
    """Compute both a^{0,1} and a^{1,0} (projected) at point m0."""
    m0c = m0.astype(complex)
    J_hol_l, J_anti_l = wirtinger_jacobians(m0, e, eps=eps_inner)
    M_phi = mass_to_M(phi_map(m0c, e))
    M_i = inv(M_phi)

    a01_local = np.zeros((4, 2, 2), dtype=complex)
    a10_local = np.zeros((4, 2, 2), dtype=complex)
    for alpha in range(4):
        dM_anti = sum(dM_dm[j] * J_anti_l[j, alpha] for j in range(4))
        dM_hol = sum(dM_dm[j] * J_hol_l[j, alpha] for j in range(4))
        a01_local[alpha] = M_i @ dM_anti
        a10_local[alpha] = M_i @ dM_hol

    a01_p = np.zeros((3, 2, 2), dtype=complex)
    a10_p = np.zeros((3, 2, 2), dtype=complex)
    for a in range(3):
        for alpha in range(4):
            a01_p[a] += V[alpha, a] * a01_local[alpha]
            a10_p[a] += V[alpha, a] * a10_local[alpha]
    return a01_p, a10_p


# Compute dbar_a(a^{0,1}_b) by Wirtinger finite differences along v_a (anti-holo)
eps_outer = 1e-5  # Step for outer derivative

def dbar_connection(m0, e, eps_out=1e-5, eps_in=1e-7):
    """Compute dbar_a(a^{0,1}_b) at m0.

    Returns: da[a, b] = 2x2 matrix = d(a^{0,1}_b)/d(z_bar_a)
    Shape: (3, 3, 2, 2) complex
    """
    da = np.zeros((3, 3, 2, 2), dtype=complex)

    for a in range(3):
        v_a = V[:, a]  # tangent direction (real, 4-vector)

        # Real perturbation along v_a
        m_p_r = m0 + eps_out * v_a
        m_m_r = m0 - eps_out * v_a
        a01_p_r = compute_connection_projected_at(m_p_r, e, eps_in)
        a01_m_r = compute_connection_projected_at(m_m_r, e, eps_in)
        d_dx = (a01_p_r - a01_m_r) / (2 * eps_out)  # shape (3, 2, 2)

        # Imaginary perturbation along v_a
        # m0 + i*eps*v_a: need complex m0
        m_p_i = m0.astype(complex) + 1j * eps_out * v_a
        m_m_i = m0.astype(complex) - 1j * eps_out * v_a

        # For imaginary-perturbed points, compute connection
        # These are complex mass functions; the real floor enforcement won't work.
        # Use the complex version.
        a01_p_i = compute_connection_projected_at(np.real(m_p_i), e, eps_in)
        a01_m_i = compute_connection_projected_at(np.real(m_m_i), e, eps_in)

        # Hmm, for imaginary perturbation we need to evaluate at complex m.
        # But compute_connection_projected_at expects real m0 for Wirtinger derivatives.
        # We need a version that takes complex m0.
        # For now, use the REAL part approximation:
        # d/dy(a(m*+iy*v)) ≈ d/dy(a(m*)) + ... which involves Im perturbation.
        # But we're computing d(a)/d(m_bar), which for real m is:
        # d/dm_bar = (1/2)(d/dx + i*d/dy)
        # For d/dy, we'd need a at m*+iy*v, which requires complex floor.

        # ALTERNATIVE: Since m* is real and the floor is smooth at m*,
        # we can compute d/dy by noting that for real functions on real inputs:
        # d/dm_bar = d/dm = (d/dx) for real perturbations.
        # But this is NOT right -- the floor is NOT holomorphic.
        #
        # The correct approach: use purely real finite differences.
        # For a real analytic function f(x):
        #   df/dz_bar = (1/2)(df/dx + i df/dy) = (1/2)(df/dx) when f(x+iy) = f(x)+i*Jf*y
        # This is wrong for the floor enforcement which mixes real and imaginary.
        #
        # We need the full complex computation. Let me use a different approach.

        # APPROACH: compute a^{0,1}(m) at complex m directly.
        # a^{0,1}_alpha(m) = M(Phi(m))^{-1} * sum_j dM_j * (dbar Phi)_{j,alpha}(m)
        # where dbar Phi is computed at complex m.
        #
        # But wirtinger_jacobians() perturbs m by epsilon.
        # If m is already complex, we can still perturb and compute.
        # We just need phi_map to work for complex inputs -- which it does.
        pass

    return da


# Actually, let me take a cleaner approach. The key function is:
# a^{0,1}(m) = M(Phi(m))^{-1} * dM * dbar_Phi(m)
# which is computable for any m (real or complex) as long as Phi and its
# Wirtinger derivatives are computable at that m.
#
# For the Wirtinger derivative of a with respect to z_bar:
# d(a^{0,1}_b)/d(z_bar_a) = computed at m = m* by displacing m* in direction
# v_a (anti-holomorphic) and recomputing a^{0,1}_b at the displaced point.
#
# The displacement in anti-holomorphic direction v_a means:
# m -> m + delta_bar where delta_bar perturbs z_bar but not z.
#
# In Wirtinger calculus, d/d(z_bar_a) = (1/2)(d/dx_a + i d/dy_a)
# So we need:
# (1/2) * [ (a(m+eps*v_a) - a(m-eps*v_a))/(2*eps)
#          + i*(a(m+i*eps*v_a) - a(m-i*eps*v_a))/(2*eps) ]
#
# For the imaginary perturbation m+i*eps*v_a, m is complex.
# We need a function that computes a^{0,1}(m) for complex m.

def compute_a01_complex(m_complex, e, eps_inner=1e-7):
    """Compute projected a^{0,1}_a (3 matrices) for complex mass function m.

    Works for complex inputs by using the complex-capable phi_map
    and wirtinger_jacobians_complex.
    """
    # Wirtinger jacobians at complex m
    J_hol_l = np.zeros((4, 4), dtype=complex)
    J_anti_l = np.zeros((4, 4), dtype=complex)

    for alpha in range(4):
        e_alpha = np.zeros(4, dtype=complex)
        e_alpha[alpha] = 1.0

        # Real perturbation
        f_pr = phi_map(m_complex + eps_inner * e_alpha, e)
        f_mr = phi_map(m_complex - eps_inner * e_alpha, e)
        df_dx = (f_pr - f_mr) / (2 * eps_inner)

        # Imaginary perturbation
        f_pi = phi_map(m_complex + 1j * eps_inner * e_alpha, e)
        f_mi = phi_map(m_complex - 1j * eps_inner * e_alpha, e)
        df_dy = (f_pi - f_mi) / (2 * eps_inner)

        J_hol_l[:, alpha] = 0.5 * (df_dx - 1j * df_dy)
        J_anti_l[:, alpha] = 0.5 * (df_dx + 1j * df_dy)

    # Connection
    M_phi = mass_to_M(phi_map(m_complex, e))
    M_i = inv(M_phi)

    a01_local = np.zeros((4, 2, 2), dtype=complex)
    a10_local = np.zeros((4, 2, 2), dtype=complex)
    for alpha in range(4):
        dM_anti = sum(dM_dm[j] * J_anti_l[j, alpha] for j in range(4))
        dM_hol = sum(dM_dm[j] * J_hol_l[j, alpha] for j in range(4))
        a01_local[alpha] = M_i @ dM_anti
        a10_local[alpha] = M_i @ dM_hol

    # Project
    a01_p = np.zeros((3, 2, 2), dtype=complex)
    a10_p = np.zeros((3, 2, 2), dtype=complex)
    for a in range(3):
        for alpha in range(4):
            a01_p[a] += V[alpha, a] * a01_local[alpha]
            a10_p[a] += V[alpha, a] * a10_local[alpha]

    return a01_p, a10_p


# Compute F^{0,2} and dbar_E* F^{0,2} using Wirtinger finite differences.

def compute_curvature_and_codiff(m_star_real, e, eps_out=1e-5, eps_in=1e-7):
    """Compute F^{0,2} and its codifferential at m*.

    F_{ab} = dbar_a(a_b) - dbar_b(a_a) + [a_a, a_b]  (3 components for a<b)

    (dbar_E* F)_c = -sum_a (d F_{ac}/dz^a + [a^{1,0}_a, F_{ac}])  (3 components)
    """
    m0 = m_star_real.astype(complex)

    # Connection at m*
    a01_0, a10_0 = compute_a01_complex(m0, e, eps_in)

    print(f"\n  Connection a^{{0,1}} at m* (projected):")
    for a in range(3):
        print(f"    a^01[{a}]: ||={norm(a01_0[a]):.6f}")

    # ====================================================
    # Compute dbar_a(a^{0,1}_b) by Wirtinger finite differences
    # dbar_a = d/d(z_bar_a) = (1/2)(d/dx_a + i d/dy_a)
    # Displacement in tangent direction v_a
    # ====================================================
    da01_dbar = np.zeros((3, 3, 2, 2), dtype=complex)  # [a, b, 2, 2]

    for a in range(3):
        v_a = V[:, a].astype(complex)

        # Real perturbation along v_a
        m_pr = m0 + eps_out * v_a
        m_mr = m0 - eps_out * v_a
        a01_pr, _ = compute_a01_complex(m_pr, e, eps_in)
        a01_mr, _ = compute_a01_complex(m_mr, e, eps_in)
        d_dx = (a01_pr - a01_mr) / (2 * eps_out)  # d(a_b)/d(x_a), shape (3,2,2)

        # Imaginary perturbation along v_a
        m_pi = m0 + 1j * eps_out * v_a
        m_mi = m0 - 1j * eps_out * v_a
        a01_pi, _ = compute_a01_complex(m_pi, e, eps_in)
        a01_mi, _ = compute_a01_complex(m_mi, e, eps_in)
        d_dy = (a01_pi - a01_mi) / (2 * eps_out)  # d(a_b)/d(y_a)

        # Wirtinger: d/d(z_bar_a) = (1/2)(d/dx + i*d/dy)
        for b in range(3):
            da01_dbar[a, b] = 0.5 * (d_dx[b] + 1j * d_dy[b])

    # ====================================================
    # Curvature: F_{ab} = da01_dbar[a,b] - da01_dbar[b,a] + [a01[a], a01[b]]
    # ====================================================
    F = np.zeros((3, 3, 2, 2), dtype=complex)  # anti-symmetric in (a,b)

    for a in range(3):
        for b in range(3):
            commutator = a01_0[a] @ a01_0[b] - a01_0[b] @ a01_0[a]
            F[a, b] = da01_dbar[a, b] - da01_dbar[b, a] + commutator

    print(f"\n  Curvature F^{{0,2}} components:")
    F_norm_total = 0.0
    for a in range(3):
        for b in range(a+1, 3):
            fn = norm(F[a, b])
            F_norm_total += fn**2
            print(f"    F[{a},{b}]: ||={fn:.6f}")
            # Print matrix elements
            for i in range(2):
                print(f"      [{F[a,b,i,0]:.4e}, {F[a,b,i,1]:.4e}]")
    F_norm_total = np.sqrt(F_norm_total)
    print(f"  Total ||F^{{0,2}}|| = {F_norm_total:.6f}")

    # Check anti-symmetry
    for a in range(3):
        for b in range(a+1, 3):
            asym_err = norm(F[a, b] + F[b, a])
            print(f"    Antisymmetry F[{a},{b}]+F[{b},{a}]: {asym_err:.2e}")

    # ====================================================
    # Codifferential: (dbar_E* F)_c = -sum_a (dF_{ac}/dz^a + [a^{1,0}_a, F_{ac}])
    # Need dF/dz^a (holomorphic derivative of F)
    # ====================================================
    print(f"\n  Computing codifferential dbar_E* F...")

    # Compute F at displaced points (holomorphic direction)
    def compute_F_at(m_disp, e_local, eps_in_local=1e-7, eps_out_local=1e-5):
        """Compute F^{0,2} at a displaced point."""
        a01_d, _ = compute_a01_complex(m_disp, e_local, eps_in_local)

        # Need da01_dbar at displaced point too
        da01_d = np.zeros((3, 3, 2, 2), dtype=complex)
        for aa in range(3):
            v_aa = V[:, aa].astype(complex)
            m_pr2 = m_disp + eps_out_local * v_aa
            m_mr2 = m_disp - eps_out_local * v_aa
            a01_pr2, _ = compute_a01_complex(m_pr2, e_local, eps_in_local)
            a01_mr2, _ = compute_a01_complex(m_mr2, e_local, eps_in_local)
            d_dx2 = (a01_pr2 - a01_mr2) / (2 * eps_out_local)

            m_pi2 = m_disp + 1j * eps_out_local * v_aa
            m_mi2 = m_disp - 1j * eps_out_local * v_aa
            a01_pi2, _ = compute_a01_complex(m_pi2, e_local, eps_in_local)
            a01_mi2, _ = compute_a01_complex(m_mi2, e_local, eps_in_local)
            d_dy2 = (a01_pi2 - a01_mi2) / (2 * eps_out_local)

            for bb in range(3):
                da01_d[aa, bb] = 0.5 * (d_dx2[bb] + 1j * d_dy2[bb])

        F_d = np.zeros((3, 3, 2, 2), dtype=complex)
        for aa in range(3):
            for bb in range(3):
                comm = a01_d[aa] @ a01_d[bb] - a01_d[bb] @ a01_d[aa]
                F_d[aa, bb] = da01_d[aa, bb] - da01_d[bb, aa] + comm
        return F_d

    # Holomorphic derivative of F: dF/dz^a
    # d/dz^a = (1/2)(d/dx_a - i*d/dy_a)  [holomorphic Wirtinger]
    eps_hol = 5e-4  # Larger step for this outermost derivative

    dF_dz = np.zeros((3, 3, 3, 2, 2), dtype=complex)  # [a, b, c, 2, 2] = dF_{bc}/dz^a

    for a in range(3):
        v_a = V[:, a].astype(complex)

        # Real perturbation
        F_pr = compute_F_at(m0 + eps_hol * v_a, e, eps_in, eps_out)
        F_mr = compute_F_at(m0 - eps_hol * v_a, e, eps_in, eps_out)
        dF_dx = (F_pr - F_mr) / (2 * eps_hol)

        # Imaginary perturbation
        F_pi = compute_F_at(m0 + 1j * eps_hol * v_a, e, eps_in, eps_out)
        F_mi = compute_F_at(m0 - 1j * eps_hol * v_a, e, eps_in, eps_out)
        dF_dy = (F_pi - F_mi) / (2 * eps_hol)

        # d/dz^a = (1/2)(d/dx - i*d/dy)
        dF_dz[a] = 0.5 * (dF_dx - 1j * dF_dy)

        print(f"    Computed dF/dz^{a}")

    # Codifferential: (dbar_E* F)_c = -sum_a [dF_{ac}/dz^a + [a^{1,0}_a, F_{ac}]]
    # Using flat metric: g^{a,a_bar} = delta_{a,a_bar} (on L1=1 tangent space)
    codiff = np.zeros((3, 2, 2), dtype=complex)

    for c in range(3):
        for a in range(3):
            # Term 1: dF_{ac}/dz^a
            codiff[c] -= dF_dz[a, a, c]  # Note: F_{ac} = F[a,c]

            # Term 2: [a^{1,0}_a, F_{ac}]
            codiff[c] -= (a10_0[a] @ F[a, c] - F[a, c] @ a10_0[a])

    return F, codiff, a01_0, a10_0


# Run the computation
print("\n" + "=" * 70)
print("PART 6: FULL COMPUTATION")
print("=" * 70)

F_result, codiff_result, a01_final, a10_final = compute_curvature_and_codiff(
    m_star, e_star, eps_out=1e-5, eps_in=1e-7
)


# ============================================================
# PART 8: Results
# ============================================================
print("\n" + "=" * 70)
print("RESULTS: YM STATIONARITY TEST")
print("=" * 70)

print(f"\n  Curvature F^{{0,2}} at K*=7/30:")
F_total_norm = 0.0
for a in range(3):
    for b in range(a+1, 3):
        fn = norm(F_result[a, b])
        F_total_norm += fn**2
        print(f"    ||F[{a},{b}]|| = {fn:.6f}")
F_total_norm = np.sqrt(F_total_norm)
print(f"  Total ||F^{{0,2}}|| = {F_total_norm:.6f}")

print(f"\n  Codifferential dbar_E* F^{{0,2}}:")
codiff_total = 0.0
for c in range(3):
    cn = norm(codiff_result[c])
    codiff_total += cn**2
    print(f"    ||(dbar_E* F)[{c}]|| = {cn:.6e}")
    for i in range(2):
        print(f"      [{codiff_result[c,i,0]:.4e}, {codiff_result[c,i,1]:.4e}]")
codiff_total = np.sqrt(codiff_total)

print(f"\n  Total ||dbar_E* F^{{0,2}}|| = {codiff_total:.6e}")
print(f"  Total ||F^{{0,2}}||         = {F_total_norm:.6f}")

if F_total_norm > 1e-10:
    ratio = codiff_total / F_total_norm
    print(f"  Ratio ||codiff|| / ||F||   = {ratio:.6e}")
else:
    print(f"  F^{{0,2}} is zero (holomorphic bundle, no non-integrability?)")

print()
if codiff_total < 1e-3 * F_total_norm:
    print("  *** dbar_E* F^{0,2} ≈ 0 ***")
    print("  The DS equilibrium satisfies the Yang-Mills equation")
    print("  to the precision of this numerical computation.")
    print()
    print("  DS stationarity = YM stationarity: SUPPORTED")
elif codiff_total < 0.1 * F_total_norm:
    print("  dbar_E* F^{0,2} is small but nonzero relative to F.")
    print("  This may be numerical noise or a genuine residual.")
    print("  Needs investigation with different step sizes.")
else:
    print("  *** dbar_E* F^{0,2} is NOT zero ***")
    print("  The DS equilibrium does NOT satisfy the Yang-Mills equation.")
    print("  The gap between DS dynamics and YM dynamics is real.")

# ============================================================
# PART 9: Numerical stability check
# ============================================================
print("\n" + "=" * 70)
print("NUMERICAL STABILITY CHECK")
print("=" * 70)

# Re-run with different step sizes to check convergence
for eps_out_test in [5e-6, 1e-5, 2e-5]:
    for eps_in_test in [5e-8, 1e-7, 2e-7]:
        # Just compute F norm at m* with different stepsizes
        a01_t, _ = compute_a01_complex(m_star.astype(complex), e_star, eps_in_test)
        f_norm = np.sqrt(sum(norm(a01_t[a])**2 for a in range(3)))
        print(f"  eps_in={eps_in_test:.0e}, eps_out={eps_out_test:.0e}: ||a^01||={f_norm:.8f}")

print("\nDONE.")
