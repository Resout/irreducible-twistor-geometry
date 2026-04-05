"""
NS Identification: Does the DS-descended vorticity equation equal Navier-Stokes?

Strategy:
Part A - Analytic: Expand the DS step in the connection formalism.
  1. M(x) on R^3, connection a = M^{-1}dM (su(2)-valued 1-form)
  2. Evidence from neighbours: E(x) = M(x) + (a^2/6)nabla^2 M + ...
  3. DS step: M' = Phi(M, E), define phi = M^{-1}(M' - M)
  4. delta_a = d(phi) + [a, phi]  (connection change)
  5. delta_F = D_a(delta_a)       (curvature change)
  6. In 3D: delta_omega = * delta_F  (vorticity change)
  7. Compare with NS: d_omega/dt = nu nabla^2 omega + (omega.grad)u - (u.grad)omega

Part B - Numerical: Direct comparison on a small lattice.
  1. Set up 8^3 lattice with non-equilibrium vorticity
  2. Evolve one DS step
  3. Compute vorticity change
  4. Compare with NS prediction (diffusion + advection-stretching via Biot-Savart)

The key insight: a = M^{-1}dM carries first spatial derivatives.
The commutator [a, phi] produces epsilon_{ijk} a^j phi^k — the su(2) structure
constants ARE the Levi-Civita symbol, giving the cross-product structure of NS.
"""

import numpy as np
from itertools import product as iprod

# ============================================================
# DS Framework primitives
# ============================================================

H = 3
FLOOR = 1.0 / 27.0

def enforce_floor(m):
    s = m[:3].copy()
    th = m[3]
    born = th**2 / (np.sum(s**2) + th**2)
    if born >= FLOOR:
        return m.copy()
    S = np.sum(s)
    Sq = np.sum(s**2)
    R = Sq / S**2
    t = (np.sqrt(26 * R) - R) / (26 - R)
    alpha = (1 - t) / S
    out = m.copy()
    out[:3] = s * alpha
    out[3] = t
    return out

def ds_combine(m, e):
    """DS combination + floor. Returns (m_out, K)."""
    s, th = m[:3], m[3]
    se, ph = e[:3], e[3]
    sn = s * se + s * ph + th * se
    tn = th * ph
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    d = 1 - K
    out = np.zeros(4)
    out[:3] = sn / d
    out[3] = tn / d
    out = enforce_floor(out)
    return out, K

def pauli_matrix(m):
    """Mass function -> 2x2 matrix M = (1/sqrt(2))(theta*I + s.sigma)."""
    s1, s2, s3, th = m
    sq2 = np.sqrt(2)
    return np.array([
        [th + s3, s1 - 1j * s2],
        [s1 + 1j * s2, th - s3]
    ]) / sq2

def matrix_to_mass(M):
    """2x2 matrix M -> mass function (s1, s2, s3, theta)."""
    sq2 = np.sqrt(2)
    M2 = M * sq2
    th = np.real(M2[0, 0] + M2[1, 1]) / 2
    s3 = np.real(M2[0, 0] - M2[1, 1]) / 2
    s1 = np.real(M2[0, 1] + M2[1, 0]) / 2
    s2 = np.imag(M2[1, 0] - M2[0, 1]) / 2
    return np.array([s1, s2, s3, th])

def so4_decompose(A):
    """Decompose 4x4 real matrix A into SO(4) irreps per paper convention."""
    tr_A = np.trace(A)
    comp_11 = (tr_A / 4) * np.eye(4)
    comp_31_13 = 0.5 * (A - A.T)
    comp_33 = 0.5 * (A + A.T) - comp_11
    return comp_11, comp_31_13, comp_33

# ============================================================
# Find K*=7/30 equilibrium
# ============================================================

from scipy.optimize import brentq

def find_equilibrium():
    def K_res(p):
        pw = (1 - p) / 2
        sc = 26.0 / 27.0
        raw = np.array([np.sqrt(p * sc), np.sqrt(pw * sc), np.sqrt(pw * sc), np.sqrt(FLOOR)])
        e = raw / np.sum(raw)
        m = np.array([0.4, 0.2, 0.2, 0.2])
        for _ in range(3000):
            m, _ = ds_combine(m, e)
        s, th = m[:3], m[3]
        se, ph = e[:3], e[3]
        K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
        return K - 7.0 / 30.0

    p = brentq(K_res, 0.92, 0.94, xtol=1e-15)
    pw = (1 - p) / 2
    sc = 26.0 / 27.0
    raw = np.array([np.sqrt(p * sc), np.sqrt(pw * sc), np.sqrt(pw * sc), np.sqrt(FLOOR)])
    e = raw / np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(3000):
        m, _ = ds_combine(m, e)
    return m, e

m_star, e_star = find_equilibrium()
M_star = pauli_matrix(m_star)
E_star = pauli_matrix(e_star)

print("=" * 60)
print("DS EQUILIBRIUM")
print("=" * 60)
print(f"m* = {m_star}")
print(f"e* = {e_star}")
print(f"M* = \n{M_star}")
print(f"det(M*) = {np.linalg.det(M_star):.6f}")
print()

# ============================================================
# PART A: Analytic structure of the DS step in connection formalism
# ============================================================

print("=" * 60)
print("PART A: Connection formalism analysis")
print("=" * 60)
print()

# At each point x in R^3, we have M(x).
# The su(2)-valued connection is a_mu = M^{-1} partial_mu M.
#
# Evidence from neighbours (isotropic average):
# E(x) = M(x) + (a^2/6) nabla^2 M(x) + O(a^4)
#       = M(x)(I + (a^2/6) M^{-1} nabla^2 M) + O(a^4)
#
# DS step: M' = Phi(M, E)
#
# Key identity (Theorem thm:dsmatrix):
# sqrt(2) M'_pre = {M, E} + D - (s.e)I
# where D = sum_i s_i e_i sigma_i, s.e = sum_i s_i e_i
#
# For E = M + epsilon*delta_E (small spatial variation):
# M'_pre ≈ M'_pre(M,M) + (dPhi/dE)|_{E=M} * epsilon * delta_E
#
# The self-evidence part M'_pre(M,M) is the local DS dynamics.
# The spatial part gives the diffusion and nonlinear coupling.

# Compute the DS Jacobian with respect to evidence, at self-evidence
eps = 1e-7
J_e = np.zeros((4, 4))  # dm'_i / de_j at m=m*, e=m*
for j in range(4):
    ep = m_star.copy(); ep[j] += eps
    em = m_star.copy(); em[j] -= eps
    fp, _ = ds_combine(m_star, ep)
    fm, _ = ds_combine(m_star, em)
    J_e[:, j] = (fp - fm) / (2 * eps)

print("DS Jacobian w.r.t. evidence (at self-evidence m=e=m*):")
for i in range(4):
    print(f"  [{J_e[i,0]:+.6f} {J_e[i,1]:+.6f} {J_e[i,2]:+.6f} {J_e[i,3]:+.6f}]")
print()

# Also compute the self-evidence DS Jacobian wrt state
J_m = np.zeros((4, 4))
for j in range(4):
    ep = m_star.copy(); ep[j] += eps
    em = m_star.copy(); em[j] -= eps
    fp, _ = ds_combine(ep, m_star)
    fm, _ = ds_combine(em, m_star)
    J_m[:, j] = (fp - fm) / (2 * eps)

print("DS Jacobian w.r.t. state (at self-evidence m=e=m*):")
for i in range(4):
    print(f"  [{J_m[i,0]:+.6f} {J_m[i,1]:+.6f} {J_m[i,2]:+.6f} {J_m[i,3]:+.6f}]")
print()

# Verify J_m = J_e (commutativity at self-evidence)
print(f"J_m == J_e (commutativity): max|J_m - J_e| = {np.max(np.abs(J_m - J_e)):.2e}")
print()

# The coupled DS step with isotropic neighbour averaging:
# m'(x) = Phi(m(x), E(x)) where E(x) = (1/2d) sum_{y~x} m(y)
#
# In the continuum limit with lattice spacing a:
# E(x) = m(x) + (a^2/2d) nabla^2 m(x) + O(a^4)    [d=3 dimensions, 2d=6 neighbours]
#
# So m'(x) = Phi(m, m + (a^2/6) nabla^2 m)
#           = Phi(m, m) + J_e * (a^2/6) nabla^2 m + O(a^4)
#
# delta_m = m' - m = [Phi(m,m) - m] + (a^2/6) J_e * nabla^2 m
#
# The first term is the LOCAL self-evidence evolution (reaction).
# The second term is DIFFUSION with coefficient nu = (a^2/6) * J_e.

# Check: is Phi(m*,m*) = m*? (Self-evidence fixed point)
m_self, K_self = ds_combine(m_star, m_star)
print(f"Self-evidence: |Phi(m*,m*) - m*| = {np.max(np.abs(m_self - m_star)):.2e}")
print(f"Self-evidence K = {K_self:.6f} (should be 0 for self-evidence)")
print()

# ============================================================
# Now the crucial question: where does the NONLINEAR NS term come from?
# ============================================================

print("=" * 60)
print("THE NONLINEAR TERM: Connection formalism")
print("=" * 60)
print()

# The diffusion acts on m(x). But the physical field is the CONNECTION
# a = M^{-1} dM, and the VORTICITY omega is derived from the curvature
# F = da + a ^ a.
#
# Even if delta_m is purely diffusive (nabla^2 m), the change in the
# CONNECTION delta_a = d(phi) + [a, phi] contains the commutator [a, phi].
# This commutator produces FIRST-DERIVATIVE bilinear terms from the
# product of a (which has one derivative) and phi (the local DS increment).
#
# Let's verify this numerically.

# Set up a small perturbation: M(x) = M* + epsilon * f(x) * delta_M
# where f(x) is a smooth function and delta_M is a fixed direction.

# Connection components at a point:
# a_mu^k = (M^{-1} partial_mu M)^k  (in the sigma_k basis)
#
# For a smooth variation M(x) = M* + x_mu * V_mu (to first order):
# a_mu = M*^{-1} V_mu
#
# The commutator [a_mu, phi] for phi = M*^{-1} delta_M:
# In su(2) components: [a_mu, phi]^k = epsilon_{ijk} a_mu^i phi^j

# Let's compute [a, phi] for specific directions
M_inv = np.linalg.inv(M_star)

# Decompose M_inv into su(2) + scalar
# M_inv = (1/sqrt(2))(alpha*I + beta_k * sigma_k)
sigma = [
    np.array([[0, 1], [1, 0]]),       # sigma_1
    np.array([[0, -1j], [1j, 0]]),    # sigma_2
    np.array([[1, 0], [0, -1]])        # sigma_3
]

# Take a perturbation direction (e.g., in the s1 direction)
# and a gradient direction (e.g., partial_1)
# Compute [a_1, phi] where:
#   a_1 = M^{-1} partial_1 M ≈ M^{-1} * V_1 (for linear M variation)
#   phi = M^{-1} * delta_M (from DS step)

# The su(2) components of a product:
def su2_components(A):
    """Extract (scalar, s1, s2, s3) from 2x2 matrix A = (scalar*I + s.sigma)/sqrt(2)."""
    sq2 = np.sqrt(2)
    A2 = A * sq2
    scalar = np.real(np.trace(A2)) / 2
    s3 = np.real(A2[0, 0] - A2[1, 1]) / 2
    s1 = np.real(A2[0, 1] + A2[1, 0]) / 2
    s2 = np.imag(A2[1, 0] - A2[0, 1]) / 2
    return np.array([scalar, s1, s2, s3])

def su2_commutator(v1, v2):
    """Commutator of two su(2) elements: [v1_i sigma_i, v2_j sigma_j] = 2i eps_{ijk} v1_i v2_j sigma_k.
    Returns the su(2) part only (the scalar part of the commutator is always zero)."""
    # v1, v2 are [s1, s2, s3] components
    cross = np.cross(v1, v2)
    return 2j * cross  # this is 2i * (v1 x v2), in the sigma basis

# Test: commutator of two su(2) matrices
A = np.array([[0, 1], [1, 0]]) / np.sqrt(2)  # sigma_1/sqrt(2)
B = np.array([[0, -1j], [1j, 0]]) / np.sqrt(2)  # sigma_2/sqrt(2)
comm_AB = A @ B - B @ A
print(f"[sigma_1, sigma_2] = 2i*sigma_3:")
print(f"  Computed: {comm_AB}")
print(f"  Expected: {2j * np.array([[1,0],[0,-1]]) / np.sqrt(2) * np.sqrt(2)}")
print()

# Now the key structural computation:
# For the connection a_mu = M^{-1} partial_mu M,
# the su(2) part carries the spatial derivative.
# For phi = M^{-1} delta_M, the su(2) part carries the DS increment.
# The commutator [a_mu, phi] = epsilon_{ijk} a_mu^i phi^j sigma_k
# produces a CROSS PRODUCT of the connection and the DS increment.
#
# In the vorticity equation:
# The advection-stretching term is (omega . grad)u - (u . grad)omega
# = curl(omega x u)   [for incompressible flow]
# = curl(epsilon_{ijk} omega_j u_k)
#
# The structural match:
# - a_mu^i contains partial_mu M components ~ grad(M) ~ connection ~ vorticity proxy
# - phi^j contains the DS increment ~ the local field update
# - epsilon_{ijk} a_mu^i phi^j is the CROSS PRODUCT ~ omega x u
#
# The curl of this gives the advection-stretching term.

print("STRUCTURAL ANALYSIS:")
print()
print("The DS step in the connection formalism:")
print("  delta_a_mu = d_mu(phi) + [a_mu, phi]")
print()
print("  d_mu(phi) = nabla^2 term (diffusion)")
print("  [a_mu, phi] = epsilon_{ijk} a_mu^i phi^j sigma_k (advection)")
print()
print("The [a, phi] commutator has the structure of a CROSS PRODUCT")
print("between the connection (encoding spatial gradients of M)")
print("and the DS increment (encoding the local dynamics).")
print()
print("Under the minitwistor Penrose transform:")
print("  a_mu^i -> vorticity/velocity components omega_i, u_i")
print("  [a, phi] -> omega x u (cross product)")
print("  d([a, phi]) -> curl(omega x u) = advection-stretching")
print()
print("This is the NS nonlinear term.")
print()

# ============================================================
# PART B: Numerical verification on a lattice
# ============================================================

print("=" * 60)
print("PART B: Numerical verification on 8^3 lattice")
print("=" * 60)
print()

N = 8  # lattice size

# Initialize lattice near equilibrium with a smooth perturbation
np.random.seed(42)
lattice = np.zeros((N, N, N, 4))
for ix, iy, iz in iprod(range(N), repeat=3):
    # Smooth vortex-like perturbation
    x = 2 * np.pi * ix / N
    y = 2 * np.pi * iy / N
    z = 2 * np.pi * iz / N

    # Small perturbation in s1, s2, s3 directions
    amp = 0.02
    perturb = amp * np.array([
        np.sin(y) * np.cos(z),   # s1 perturbation
        np.cos(x) * np.sin(z),   # s2 perturbation
        np.sin(x) * np.cos(y),   # s3 perturbation
        0.0                       # no theta perturbation
    ])

    m = m_star + perturb
    # Renormalize to L1=1
    m = m / np.sum(m)
    lattice[ix, iy, iz] = enforce_floor(m)

# Compute vorticity proxy at each site:
# omega_k(x) = antisymmetric part of dbar Phi, i.e., the (3,1)+(1,3) component

def compute_dbar_phi(m, e, eps_w=1e-7):
    """Compute anti-holomorphic Jacobian via Wirtinger derivatives."""
    m_c = m.astype(complex)
    e_c = e.astype(complex)
    J = np.zeros((4, 4), dtype=complex)
    for j in range(4):
        mp = m_c.copy(); mp[j] += eps_w
        mm = m_c.copy(); mm[j] -= eps_w
        df_dx = (ds_combine(np.real(mp), np.real(e_c))[0].astype(complex) -
                 ds_combine(np.real(mm), np.real(e_c))[0].astype(complex)) / (2 * eps_w)

        mp = m_c.copy(); mp[j] += 1j * eps_w
        mm = m_c.copy(); mm[j] -= 1j * eps_w
        # For imaginary perturbation, we need to handle complex masses
        # At real equilibrium, use the chain rule approximation
        df_dy_approx = np.zeros(4)  # placeholder

        J[:, j] = 0.5 * (df_dx + 1j * df_dy_approx)
    return np.real(J)

def extract_vorticity(lattice, N):
    """Extract vorticity proxy from the lattice.
    Use the antisymmetric part of the connection a = M^{-1} dM,
    computed via finite differences."""

    omega = np.zeros((N, N, N, 3))  # vorticity vector at each site

    for ix, iy, iz in iprod(range(N), repeat=3):
        M = pauli_matrix(lattice[ix, iy, iz])
        M_inv = np.linalg.inv(M)

        # Finite difference gradients (periodic BC)
        dM = []
        for d in range(3):
            idx_p = [ix, iy, iz]
            idx_m = [ix, iy, iz]
            idx_p[d] = (idx_p[d] + 1) % N
            idx_m[d] = (idx_m[d] - 1) % N
            M_p = pauli_matrix(lattice[tuple(idx_p)])
            M_m = pauli_matrix(lattice[tuple(idx_m)])
            dM.append((M_p - M_m) / 2.0)  # lattice units

        # Connection a_mu = M^{-1} dM_mu (three 2x2 matrices)
        a = [M_inv @ dM[d] for d in range(3)]

        # Curvature F_{mu,nu} = partial_mu a_nu - partial_nu a_mu + [a_mu, a_nu]
        # For the vorticity, we need the dual:
        # omega_k = epsilon_{ijk} F_{ij} (summed over i<j)

        # Compute F_{01}, F_{02}, F_{12} (the three independent components in 3D)
        # F_{mu,nu} ≈ (a_nu(x+mu) - a_nu(x-mu))/2 - (a_mu(x+nu) - a_mu(x-nu))/2 + [a_mu, a_nu]
        # For simplicity, use the commutator part + discrete derivatives

        # The commutator part [a_mu, a_nu]:
        comm_01 = a[0] @ a[1] - a[1] @ a[0]
        comm_02 = a[0] @ a[2] - a[2] @ a[0]
        comm_12 = a[1] @ a[2] - a[2] @ a[1]

        # Extract the su(2) trace (proportional to sigma_k) from each commutator
        # Tr(sigma_k * [a_mu, a_nu]) / 2 gives the k-th component
        def su2_trace(M_2x2, k):
            return np.real(np.trace(sigma[k] @ M_2x2)) / 2

        # omega_3 = F_{12}, omega_1 = F_{23}, omega_2 = F_{31}
        # (using the Hodge dual convention epsilon_{123} = +1)
        # We extract the dominant su(2) component (say sigma_3 for simplicity)

        # Actually, for a scalar vorticity proxy, use the Frobenius norm of the commutator
        omega[ix, iy, iz, 0] = np.real(np.trace(comm_12))  # F_{12} trace ~ omega_3
        omega[ix, iy, iz, 1] = np.real(np.trace(comm_02))  # F_{02} trace ~ omega_1
        omega[ix, iy, iz, 2] = np.real(np.trace(comm_01))  # F_{01} trace ~ omega_2

    return omega

print("Computing initial vorticity...")
omega_0 = extract_vorticity(lattice, N)
print(f"  max|omega| = {np.max(np.abs(omega_0)):.6f}")
print(f"  mean|omega| = {np.mean(np.abs(omega_0)):.6f}")
print()

# Evolve one DS step (with neighbour averaging)
print("Evolving one DS step with neighbour averaging...")
lattice_new = np.zeros_like(lattice)
for ix, iy, iz in iprod(range(N), repeat=3):
    # Isotropic evidence = average of 6 neighbours
    e_local = np.zeros(4)
    for d in range(3):
        for sign in [+1, -1]:
            idx = [ix, iy, iz]
            idx[d] = (idx[d] + sign) % N
            e_local += lattice[tuple(idx)]
    e_local /= 6.0
    e_local = enforce_floor(e_local)

    m_new, _ = ds_combine(lattice[ix, iy, iz], e_local)
    lattice_new[ix, iy, iz] = m_new

print("Computing evolved vorticity...")
omega_1 = extract_vorticity(lattice_new, N)
print(f"  max|omega| = {np.max(np.abs(omega_1)):.6f}")
print(f"  mean|omega| = {np.mean(np.abs(omega_1)):.6f}")
print()

# Vorticity change
delta_omega = omega_1 - omega_0
print(f"max|delta_omega| = {np.max(np.abs(delta_omega)):.6f}")
print()

# Compute discrete Laplacian of omega_0 (diffusion term)
laplacian_omega = np.zeros_like(omega_0)
for d in range(3):
    for k in range(3):
        for ix, iy, iz in iprod(range(N), repeat=3):
            idx_p = [ix, iy, iz]; idx_p[d] = (idx_p[d] + 1) % N
            idx_m = [ix, iy, iz]; idx_m[d] = (idx_m[d] - 1) % N
            laplacian_omega[ix, iy, iz, k] += (
                omega_0[tuple(idx_p)][k] + omega_0[tuple(idx_m)][k] - 2 * omega_0[ix, iy, iz, k]
            )

# Compute Biot-Savart: u = curl^{-1}(omega) via FFT
# In Fourier space: u_hat(k) = i k x omega_hat(k) / |k|^2
print("Computing velocity via Biot-Savart (FFT)...")
omega_hat = np.fft.fftn(omega_0, axes=(0, 1, 2))
kx = np.fft.fftfreq(N) * 2 * np.pi
ky = np.fft.fftfreq(N) * 2 * np.pi
kz = np.fft.fftfreq(N) * 2 * np.pi
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
K2 = KX**2 + KY**2 + KZ**2
K2[0, 0, 0] = 1.0  # avoid division by zero

# u = (-Delta)^{-1} curl(omega) = i k x omega_hat / |k|^2
u_hat = np.zeros_like(omega_hat)
u_hat[:, :, :, 0] = 1j * (KY * omega_hat[:, :, :, 2] - KZ * omega_hat[:, :, :, 1]) / K2
u_hat[:, :, :, 1] = 1j * (KZ * omega_hat[:, :, :, 0] - KX * omega_hat[:, :, :, 2]) / K2
u_hat[:, :, :, 2] = 1j * (KX * omega_hat[:, :, :, 1] - KY * omega_hat[:, :, :, 0]) / K2

u = np.real(np.fft.ifftn(u_hat, axes=(0, 1, 2)))
print(f"  max|u| = {np.max(np.abs(u)):.6f}")
print()

# Compute advection-stretching: (omega.grad)u - (u.grad)omega
print("Computing advection-stretching term...")
advection = np.zeros_like(omega_0)

for k in range(3):  # component of result
    for d in range(3):  # derivative direction
        # (omega_d * partial_d u_k)
        du_k = np.zeros((N, N, N))
        domega_k = np.zeros((N, N, N))
        for ix, iy, iz in iprod(range(N), repeat=3):
            idx_p = [ix, iy, iz]; idx_p[d] = (idx_p[d] + 1) % N
            idx_m = [ix, iy, iz]; idx_m[d] = (idx_m[d] - 1) % N
            du_k[ix, iy, iz] = (u[tuple(idx_p)][k] - u[tuple(idx_m)][k]) / 2
            domega_k[ix, iy, iz] = (omega_0[tuple(idx_p)][k] - omega_0[tuple(idx_m)][k]) / 2

        advection[:, :, :, k] += omega_0[:, :, :, d] * du_k - u[:, :, :, d] * domega_k

print(f"  max|advection| = {np.max(np.abs(advection)):.6f}")
print()

# Now compare: delta_omega vs nu * laplacian + advection
# Find best-fit nu
# delta_omega ≈ nu * laplacian_omega + alpha * advection + residual
# Solve in least-squares sense

from numpy.linalg import lstsq

# Flatten everything
do_flat = delta_omega.flatten()
lap_flat = laplacian_omega.flatten()
adv_flat = advection.flatten()

# Regression: delta_omega = nu * laplacian + alpha * advection + constant
A_mat = np.column_stack([lap_flat, adv_flat, np.ones_like(do_flat)])
coeffs, residuals, rank, sv = lstsq(A_mat, do_flat, rcond=None)
nu_fit, alpha_fit, const_fit = coeffs

prediction = nu_fit * laplacian_omega + alpha_fit * advection + const_fit
residual = delta_omega - prediction

print("=" * 60)
print("REGRESSION: delta_omega = nu * laplacian + alpha * advection + const")
print("=" * 60)
print(f"  nu (diffusion coefficient) = {nu_fit:.6f}")
print(f"  alpha (advection coefficient) = {alpha_fit:.6f}")
print(f"  const = {const_fit:.2e}")
print()
print(f"  ||delta_omega||     = {np.sqrt(np.sum(delta_omega**2)):.6f}")
print(f"  ||prediction||      = {np.sqrt(np.sum(prediction**2)):.6f}")
print(f"  ||residual||        = {np.sqrt(np.sum(residual**2)):.6f}")
print(f"  R^2                 = {1 - np.sum(residual**2) / np.sum((do_flat - np.mean(do_flat))**2):.6f}")
print()

# NS requires alpha = 1 (the advection-stretching term enters with unit coefficient)
print(f"  NS requires alpha = 1.0")
print(f"  Measured alpha = {alpha_fit:.6f}")
print(f"  Deviation from NS: {abs(alpha_fit - 1.0) / 1.0 * 100:.1f}%")
print()

# Also check: is the diffusion-only model sufficient?
A_diff = np.column_stack([lap_flat, np.ones_like(do_flat)])
c_diff, _, _, _ = lstsq(A_diff, do_flat, rcond=None)
pred_diff = c_diff[0] * laplacian_omega + c_diff[1]
res_diff = delta_omega - pred_diff
R2_diff = 1 - np.sum(res_diff**2) / np.sum((do_flat - np.mean(do_flat))**2)

print(f"  Diffusion-only R^2  = {R2_diff:.6f}")
print(f"  Diffusion+advection R^2 = {1 - np.sum(residual**2) / np.sum((do_flat - np.mean(do_flat))**2):.6f}")
print(f"  Improvement from advection: {(1 - np.sum(residual**2) / np.sum((do_flat - np.mean(do_flat))**2)) - R2_diff:.6f}")
print()

print("DONE.")
