"""
Diagnostic: separate the two contributions to F^{0,2} at K*=7/30.

F_{ab} = [dbar_a(a_b) - dbar_b(a_a)]  +  [a_a, a_b]
       = (derivative part)              +  (commutator part)

If both are O(1) and nearly cancel to give ||F|| ~ 10^{-5},
that's a genuine near-cancellation with physical meaning.

If one is O(10^{-5}) and the other is O(10^{-5}), the small F
is simply because both terms are small.

Also: compute F at displaced points to check if F is small
only at m* or everywhere nearby.
"""
import numpy as np
from numpy.linalg import inv, det, norm

H = 3
FLOOR = 1.0 / H**3

sigma = [
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
]
I2 = np.eye(2, dtype=complex)
V = np.array([[1, 0, 0, -1], [0, 1, 0, -1], [0, 0, 1, -1]], dtype=float).T
dM_dm = [sigma[i]/np.sqrt(2) for i in range(3)] + [I2/np.sqrt(2)]

def mass_to_M(m):
    return (m[3]*I2 + m[0]*sigma[0] + m[1]*sigma[1] + m[2]*sigma[2]) / np.sqrt(2)

def ds_combine_complex(m, e):
    s, theta = m[:3], m[3]
    se, phi = e[:3], e[3]
    s_new = s * se + s * phi + theta * se
    theta_new = theta * phi
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    denom = 1.0 - K
    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom
    return m_out, K

def enforce_floor_complex(m):
    s = m[:3].copy()
    theta = m[3]
    s_mag_sq = np.sum(np.abs(s)**2)
    theta_mag = np.abs(theta)
    total_mag_sq = s_mag_sq + theta_mag**2
    if total_mag_sq < 1e-30:
        return m.copy()
    born = theta_mag**2 / total_mag_sq
    if born >= FLOOR - 1e-15:
        return m.copy()
    S = np.sum(s)
    S_mag_sq = np.abs(S)**2
    if S_mag_sq < 1e-30:
        return m.copy()
    phase = theta / theta_mag if theta_mag > 1e-30 else 1.0+0j
    R = s_mag_sq / S_mag_sq
    Re_phase = np.real(phase)
    a_c = 26.0 - R
    b_c = 2.0 * R * Re_phase
    c_c = -R
    disc = b_c**2 - 4*a_c*c_c
    if disc < 0:
        return m.copy()
    r1 = (-b_c + np.sqrt(disc)) / (2*a_c)
    r2 = (-b_c - np.sqrt(disc)) / (2*a_c)
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

def phi_map(m_complex, e):
    m_ds, K = ds_combine_complex(m_complex, e)
    return enforce_floor_complex(m_ds)

# Equilibrium
def enforce_floor_real(m):
    s = m[:3].copy()
    S = np.sum(s)
    lo, hi = 0.0, 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        s_scale = (1.0 - mid) / S if S > 0 else 0
        s_trial = s * s_scale
        born = mid**2 / (np.sum(s_trial**2) + mid**2)
        if born < FLOOR: lo = mid
        else: hi = mid
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
        if np.max(np.abs(m_new - m)) < 1e-15: break
        m = m_new
    _, K = ds_combine_real(m, e)
    return K

p_dom_exact = brentq(lambda p: K_at_equil(p) - 7.0/30.0, 0.90, 0.95, xtol=1e-14)
p_weak_exact = (1 - p_dom_exact) / 2
scale = 1 - FLOOR
raw = np.array([np.sqrt(p_dom_exact*scale), np.sqrt(p_weak_exact*scale),
                np.sqrt(p_weak_exact*scale), np.sqrt(FLOOR)])
e_star = raw / np.sum(raw)

m = np.array([0.4, 0.2, 0.2, 0.2])
for _ in range(5000):
    m_new, _ = ds_combine_real(m, e_star)
    if np.max(np.abs(m_new - m)) < 1e-15: break
    m = m_new
m_star = m_new.copy()

print(f"m* = ({m_star[0]:.10f}, {m_star[1]:.10f}, {m_star[2]:.10f}, {m_star[3]:.10f})")

# ============================================================
# Compute connection at a point
# ============================================================
def wirtinger_jacobians(m0, e, eps=1e-7):
    m0c = m0.astype(complex) if np.isrealobj(m0) else m0
    J_hol = np.zeros((4, 4), dtype=complex)
    J_anti = np.zeros((4, 4), dtype=complex)
    for alpha in range(4):
        e_a = np.zeros(4, dtype=complex)
        e_a[alpha] = 1.0
        f_pr = phi_map(m0c + eps * e_a, e)
        f_mr = phi_map(m0c - eps * e_a, e)
        df_dx = (f_pr - f_mr) / (2 * eps)
        f_pi = phi_map(m0c + 1j * eps * e_a, e)
        f_mi = phi_map(m0c - 1j * eps * e_a, e)
        df_dy = (f_pi - f_mi) / (2 * eps)
        J_hol[:, alpha] = 0.5 * (df_dx - 1j * df_dy)
        J_anti[:, alpha] = 0.5 * (df_dx + 1j * df_dy)
    return J_hol, J_anti

def get_projected_connection(m0, e, eps=1e-7):
    """Returns (a01_proj, a10_proj) each shape (3,2,2) complex."""
    m0c = m0.astype(complex) if np.isrealobj(m0) else m0
    J_hol, J_anti = wirtinger_jacobians(m0, e, eps)
    M_phi = mass_to_M(phi_map(m0c, e))
    M_i = inv(M_phi)
    a01 = np.zeros((4, 2, 2), dtype=complex)
    a10 = np.zeros((4, 2, 2), dtype=complex)
    for alpha in range(4):
        dM_anti = sum(dM_dm[j] * J_anti[j, alpha] for j in range(4))
        dM_hol = sum(dM_dm[j] * J_hol[j, alpha] for j in range(4))
        a01[alpha] = M_i @ dM_anti
        a10[alpha] = M_i @ dM_hol
    a01_p = np.zeros((3, 2, 2), dtype=complex)
    a10_p = np.zeros((3, 2, 2), dtype=complex)
    for a in range(3):
        for alpha in range(4):
            a01_p[a] += V[alpha, a] * a01[alpha]
            a10_p[a] += V[alpha, a] * a10[alpha]
    return a01_p, a10_p


# ============================================================
# DIAGNOSTIC 1: Separate commutator and derivative contributions
# ============================================================
print("\n" + "=" * 70)
print("DIAGNOSTIC 1: Commutator vs derivative in F = dbar(a) + [a,a]")
print("=" * 70)

a01_0, a10_0 = get_projected_connection(m_star, e_star)

# Commutator [a_a, a_b]
print("\nCommutator norms ||[a_a, a_b]||:")
for a in range(3):
    for b in range(a+1, 3):
        comm = a01_0[a] @ a01_0[b] - a01_0[b] @ a01_0[a]
        print(f"  [{a},{b}]: {norm(comm):.6e}")
        print(f"    {comm[0,0]:.4e}, {comm[0,1]:.4e}")
        print(f"    {comm[1,0]:.4e}, {comm[1,1]:.4e}")

# Derivative dbar_a(a_b) - dbar_b(a_a)
eps_out = 1e-5
eps_in = 1e-7

print(f"\nDerivative norms ||dbar_a(a_b) - dbar_b(a_a)||:")
print(f"  (eps_out={eps_out:.0e}, eps_in={eps_in:.0e})")

# Compute dbar_a(a_b) at m*
da_dbar = np.zeros((3, 3, 2, 2), dtype=complex)

for a in range(3):
    v_a = V[:, a].astype(complex)
    m0c = m_star.astype(complex)

    # Real perturbation
    a01_pr, _ = get_projected_connection(np.real(m0c + eps_out * v_a), e_star, eps_in)
    a01_mr, _ = get_projected_connection(np.real(m0c - eps_out * v_a), e_star, eps_in)
    d_dx = (a01_pr - a01_mr) / (2 * eps_out)

    # Imaginary perturbation
    a01_pi, _ = get_projected_connection(m0c + 1j * eps_out * v_a, e_star, eps_in)
    a01_mi, _ = get_projected_connection(m0c - 1j * eps_out * v_a, e_star, eps_in)
    d_dy = (a01_pi - a01_mi) / (2 * eps_out)

    for b in range(3):
        da_dbar[a, b] = 0.5 * (d_dx[b] + 1j * d_dy[b])

for a in range(3):
    for b in range(a+1, 3):
        deriv_part = da_dbar[a, b] - da_dbar[b, a]
        comm_part = a01_0[a] @ a01_0[b] - a01_0[b] @ a01_0[a]
        F_ab = deriv_part + comm_part

        print(f"\n  F[{a},{b}] decomposition:")
        print(f"    ||derivative part|| = {norm(deriv_part):.6e}")
        print(f"    ||commutator part|| = {norm(comm_part):.6e}")
        print(f"    ||F = deriv + comm|| = {norm(F_ab):.6e}")
        if max(norm(deriv_part), norm(comm_part)) > 1e-10:
            print(f"    cancellation ratio = {norm(F_ab)/max(norm(deriv_part), norm(comm_part)):.6e}")


# ============================================================
# DIAGNOSTIC 2: F at displaced points
# ============================================================
print("\n" + "=" * 70)
print("DIAGNOSTIC 2: ||F|| at displaced points")
print("=" * 70)

def compute_F_at_point(m0, e, eps_out=1e-5, eps_in=1e-7):
    """Compute F^{0,2} at a point."""
    a01_loc, _ = get_projected_connection(m0, e, eps_in)

    da_loc = np.zeros((3, 3, 2, 2), dtype=complex)
    m0c = m0.astype(complex) if np.isrealobj(m0) else m0
    for a in range(3):
        v_a = V[:, a].astype(complex)
        a01_pr, _ = get_projected_connection(np.real(m0c + eps_out*v_a), e, eps_in)
        a01_mr, _ = get_projected_connection(np.real(m0c - eps_out*v_a), e, eps_in)
        d_dx = (a01_pr - a01_mr) / (2*eps_out)

        a01_pi, _ = get_projected_connection(m0c + 1j*eps_out*v_a, e, eps_in)
        a01_mi, _ = get_projected_connection(m0c - 1j*eps_out*v_a, e, eps_in)
        d_dy = (a01_pi - a01_mi) / (2*eps_out)

        for b in range(3):
            da_loc[a, b] = 0.5 * (d_dx[b] + 1j * d_dy[b])

    F_loc = np.zeros((3, 3, 2, 2), dtype=complex)
    for a in range(3):
        for b in range(3):
            comm = a01_loc[a] @ a01_loc[b] - a01_loc[b] @ a01_loc[a]
            F_loc[a, b] = da_loc[a, b] - da_loc[b, a] + comm

    total = 0
    for a in range(3):
        for b in range(a+1, 3):
            total += norm(F_loc[a, b])**2
    return np.sqrt(total), F_loc

# F at m*
F_norm_star, _ = compute_F_at_point(m_star, e_star)
print(f"\n  At m*: ||F|| = {F_norm_star:.6e}")

# F at displaced points along each tangent direction
for disp in [1e-4, 1e-3, 1e-2, 5e-2]:
    print(f"\n  Displacement = {disp:.0e}:")
    for a in range(3):
        m_disp = m_star + disp * V[:, a]
        # Check L1 norm
        if abs(np.sum(m_disp) - 1.0) > 1e-10:
            print(f"    WARNING: L1 = {np.sum(m_disp):.10f}")
        fn, _ = compute_F_at_point(m_disp, e_star)
        print(f"    direction v_{a}: ||F|| = {fn:.6e}")


# ============================================================
# DIAGNOSTIC 3: Step size convergence for F
# ============================================================
print("\n" + "=" * 70)
print("DIAGNOSTIC 3: Step size convergence")
print("=" * 70)

for eps_o in [5e-6, 1e-5, 2e-5, 5e-5]:
    for eps_i in [5e-8, 1e-7, 2e-7]:
        fn, _ = compute_F_at_point(m_star, e_star, eps_out=eps_o, eps_in=eps_i)
        print(f"  eps_out={eps_o:.0e}, eps_in={eps_i:.0e}: ||F|| = {fn:.6e}")

print("\nDONE.")
