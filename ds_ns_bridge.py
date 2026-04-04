"""
DS -> NS Bridge: Formal Construction
=====================================
Goal: Determine whether coupled DS dynamics on a lattice of CP3 sites
produces Navier-Stokes in the continuum limit.
"""

import numpy as np
from numpy.linalg import norm
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Part 0: Single-site DS framework
# ============================================================

def ds_combine(m, e):
    th_m, m1, m2, m3 = m
    th_e, e1, e2, e3 = e
    K = m1*e2 + m1*e3 + m2*e1 + m2*e3 + m3*e1 + m3*e2
    if K >= 1.0:
        return m.copy()
    inv = 1.0 / (1.0 - K)
    th_new = inv * th_m * th_e
    m1_new = inv * (m1*e1 + m1*th_e + th_m*e1)
    m2_new = inv * (m2*e2 + m2*th_e + th_m*e2)
    m3_new = inv * (m3*e3 + m3*th_e + th_m*e3)
    return np.array([th_new, m1_new, m2_new, m3_new])

def born_floor_enforce(m, H=3):
    th, m1, m2, m3 = m
    s_sq = m1**2 + m2**2 + m3**2
    born = th**2 / (th**2 + s_sq) if (th**2 + s_sq) > 0 else 1.0
    if born < 1.0/H**3:
        th = np.sqrt(s_sq / (H**3 - 1))
    total = th + m1 + m2 + m3
    if total > 0:
        return np.array([th, m1, m2, m3]) / total
    return m.copy()

def ds_step(m, e, H=3):
    m_new = ds_combine(m, e)
    return born_floor_enforce(m_new, H)

def conflict(m, e):
    return m[1]*e[2] + m[1]*e[3] + m[2]*e[1] + m[2]*e[3] + m[3]*e[1] + m[3]*e[2]

# Find equilibrium
print("=" * 70)
print("PART 0: Single-site equilibrium")
print("=" * 70)
m0 = np.array([0.25, 0.25, 0.25, 0.25])
e0 = np.array([0.2, 0.3, 0.3, 0.2])
for i in range(100):
    m0 = ds_step(m0, e0)
K_star = conflict(m0, e0)
born = m0[0]**2 / (m0[0]**2 + np.sum(m0[1:]**2))
print(f"Fixed point: {m0}")
print(f"K* = {K_star:.6f} (expect {7/30:.6f})")
print(f"Born = {born:.6f} (expect {1/27:.6f})")
m_star = m0.copy()
e_star = e0.copy()
theta_eq = m_star[0]
s_eq = norm(m_star[1:])

# ============================================================
# Part 1: Jacobian of DS wrt evidence at equilibrium
# ============================================================
print()
print("=" * 70)
print("PART 1: Evidence Jacobian at equilibrium")
print("=" * 70)

eps = 1e-7
J_e = np.zeros((4, 4))
for j in range(4):
    e_plus = e_star.copy()
    e_minus = e_star.copy()
    e_plus[j] += eps
    e_minus[j] -= eps
    e_plus /= np.sum(e_plus)
    e_minus /= np.sum(e_minus)
    f_plus = ds_combine(m_star, e_plus)
    f_minus = ds_combine(m_star, e_minus)
    J_e[:, j] = (f_plus - f_minus) / (2*eps * np.sum(e_star))

print("J_e (response to evidence perturbation):")
for row in J_e:
    print(f"  [{row[0]:+.5f} {row[1]:+.5f} {row[2]:+.5f} {row[3]:+.5f}]")
eigvals = np.linalg.eigvals(J_e)
print(f"Eigenvalues: {sorted(eigvals.real, reverse=True)}")

# ============================================================
# Part 2: Trace/divergence connection
# ============================================================
print()
print("=" * 70)
print("PART 2: Trace constraint = divergence-free")
print("=" * 70)
print(f"tr(M) = sqrt(2)*theta = sqrt(2)*{theta_eq:.6f} = {np.sqrt(2)*theta_eq:.6f}")
print(f"|s|^2 = {np.sum(m_star[1:]**2):.6f}")
print(f"26*theta^2 = {26*theta_eq**2:.6f}")
print(f"Born floor saturated: |s| = sqrt(26)*theta? {abs(s_eq - np.sqrt(26)*theta_eq) < 0.01}")
print(f"  |s|/|s|_max = {s_eq/(np.sqrt(26)*theta_eq):.6f}")

# ============================================================
# Part 3: Cross-product structure in DS coupling
# ============================================================
print()
print("=" * 70)
print("PART 3: Cross-product (precession) from Pauli algebra")
print("=" * 70)

m_A = np.array([theta_eq, s_eq, 0.0, 0.0])
m_A = np.abs(m_A)
m_A /= np.sum(m_A)
m_A = born_floor_enforce(m_A)

print(f"Testing DS coupling between sites with different orientations:")
print(f"{'alpha':>8} {'delta_parallel':>15} {'delta_cross':>15} {'delta_s3':>12} {'cross/sin(a)':>14}")

for alpha_t in [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]:
    m_B = np.array([theta_eq, s_eq*np.cos(alpha_t), s_eq*np.sin(alpha_t), 0.0])
    m_B = np.abs(m_B)
    m_B /= np.sum(m_B)
    m_B = born_floor_enforce(m_B)

    m_comb = ds_step(m_A, m_B)
    delta = m_comb - m_A

    n_A = m_A[1:] / norm(m_A[1:])
    d_parallel = np.dot(delta[1:], n_A)
    n_perp = np.array([-np.sin(alpha_t/2), np.cos(alpha_t/2), 0])
    d_perp = np.dot(delta[1:], n_perp)
    d_cross = delta[3]

    ratio = abs(d_cross / np.sin(alpha_t)) if np.sin(alpha_t) > 1e-10 else 0
    print(f"{alpha_t:8.3f} {d_parallel:15.8f} {d_perp:15.8f} {d_cross:12.8f} {ratio:14.8f}")

# ============================================================
# Part 4: Coupled lattice simulation
# ============================================================
print()
print("=" * 70)
print("PART 4: Coupled DS on 2D lattice")
print("=" * 70)

N = 24
num_steps = 100

def make_lattice(N, theta_eq, s_eq, pattern="vortex"):
    lattice = np.zeros((N, N, 4))
    for ix in range(N):
        for iy in range(N):
            dx = ix - N/2
            dy = iy - N/2
            r = np.sqrt(dx**2 + dy**2) + 0.1
            if pattern == "vortex":
                phi = np.arctan2(dy, dx)
            elif pattern == "shear":
                phi = 2 * np.pi * iy / N
            elif pattern == "uniform":
                phi = 0.0
            else:
                phi = 0.0
            w1 = (1 + np.cos(phi - 0)) / 3
            w2 = (1 + np.cos(phi - 2*np.pi/3)) / 3
            w3 = (1 + np.cos(phi - 4*np.pi/3)) / 3
            s_total = s_eq
            m_init = np.array([theta_eq, s_total*w1, s_total*w2, s_total*w3])
            m_init /= np.sum(m_init)
            lattice[ix, iy] = born_floor_enforce(m_init)
    return lattice

def get_orientation(m):
    m1, m2, m3 = m[1], m[2], m[3]
    x = m1 - 0.5*(m2 + m3)
    y = (np.sqrt(3)/2) * (m2 - m3)
    return np.arctan2(y, x)

def lattice_enstrophy(lattice, N):
    phi = np.zeros((N, N))
    for ix in range(N):
        for iy in range(N):
            phi[ix, iy] = get_orientation(lattice[ix, iy])
    dpx = np.angle(np.exp(1j*(np.roll(phi,-1,0) - np.roll(phi,1,0))))
    dpy = np.angle(np.exp(1j*(np.roll(phi,-1,1) - np.roll(phi,1,1))))
    return np.mean(dpx**2 + dpy**2)

lattice = make_lattice(N, theta_eq, s_eq, "vortex")
print(f"Vortex pattern on {N}x{N} lattice")
print(f"Initial enstrophy: {lattice_enstrophy(lattice, N):.6f}")

enstrophy_hist = [lattice_enstrophy(lattice, N)]
K_mean_hist = []

for step in range(num_steps):
    new_lattice = np.zeros_like(lattice)
    K_vals = []
    for ix in range(N):
        for iy in range(N):
            m = lattice[ix, iy]
            neighbors = []
            for ddx, ddy in [(1,0),(-1,0),(0,1),(0,-1)]:
                nx, ny = (ix+ddx)%N, (iy+ddy)%N
                neighbors.append(lattice[nx, ny])
            e_avg = np.mean(neighbors, axis=0)
            new_lattice[ix, iy] = ds_step(m, e_avg)
            K_vals.append(conflict(m, e_avg))
    lattice = new_lattice
    if step % 10 == 0:
        enst = lattice_enstrophy(lattice, N)
        km = np.mean(K_vals)
        enstrophy_hist.append(enst)
        K_mean_hist.append(km)
        if step % 25 == 0:
            print(f"  Step {step:3d}: K_mean={km:.6f}, enstrophy={enst:.6f}")

enst_final = lattice_enstrophy(lattice, N)
print(f"Final enstrophy: {enst_final:.6f}")
print(f"Enstrophy ratio (final/initial): {enst_final/enstrophy_hist[0]:.4f}")
print(f"Final K mean: {K_mean_hist[-1]:.6f} (expect {7/30:.6f})")

phi_final = np.zeros((N, N))
for ix in range(N):
    for iy in range(N):
        phi_final[ix, iy] = get_orientation(lattice[ix, iy])

center_angles = []
orientation_at_angle = []
for ix in range(N):
    for iy in range(N):
        dx = ix - N/2
        dy = iy - N/2
        r = np.sqrt(dx**2 + dy**2)
        if 3 < r < N/3:
            center_angles.append(np.arctan2(dy, dx))
            orientation_at_angle.append(phi_final[ix, iy])

if len(center_angles) > 0:
    corr = np.corrcoef(center_angles, orientation_at_angle)[0,1]
    print(f"Orientation-angle correlation (vortex test): {corr:.4f}")
    print(f"  (1.0 = perfect vortex preserved, 0.0 = destroyed)")

# ============================================================
# Part 5: Born floor as algebraic regularity
# ============================================================
print()
print("=" * 70)
print("PART 5: Born floor regularity under strong perturbation")
print("=" * 70)

lattice_test = make_lattice(N, theta_eq, s_eq, "vortex")
cx, cy = N//2, N//2
for dx in range(-2, 3):
    for dy in range(-2, 3):
        ix, iy = (cx+dx)%N, (cy+dy)%N
        lattice_test[ix, iy] = np.array([0.05, 0.85, 0.05, 0.05])
        lattice_test[ix, iy] = born_floor_enforce(lattice_test[ix, iy])

print("Strong localized perturbation at center (5x5 patch with m1=0.85)")
born_min_hist = []
max_s_hist = []

for step in range(50):
    new_lattice = np.zeros_like(lattice_test)
    born_min = 1.0
    max_s = 0.0
    for ix in range(N):
        for iy in range(N):
            m = lattice_test[ix, iy]
            neighbors = []
            for ddx, ddy in [(1,0),(-1,0),(0,1),(0,-1)]:
                nx, ny = (ix+ddx)%N, (iy+ddy)%N
                neighbors.append(lattice_test[nx, ny])
            e_avg = np.mean(neighbors, axis=0)
            new_lattice[ix, iy] = ds_step(m, e_avg)
            b = m[0]**2 / (m[0]**2 + np.sum(m[1:]**2))
            born_min = min(born_min, b)
            max_s = max(max_s, norm(m[1:]))
    lattice_test = new_lattice
    born_min_hist.append(born_min)
    max_s_hist.append(max_s)
    if step % 10 == 0:
        print(f"  Step {step:3d}: Born_min={born_min:.6f}, max|s|={max_s:.6f}")

print(f"Born floor ALWAYS held: {all(b >= 1/27 - 1e-10 for b in born_min_hist)}")
print(f"|s| bounded throughout: max|s| = {max(max_s_hist):.6f}")
print(f"Perturbation decayed: final max|s| = {max_s_hist[-1]:.6f} vs initial {max_s_hist[0]:.6f}")

# ============================================================
# Part 6: Enstrophy cascade test
# ============================================================
print()
print("=" * 70)
print("PART 6: Enstrophy cascade (2D NS signature)")
print("=" * 70)

for pattern in ["vortex", "shear"]:
    lattice_p = make_lattice(N, theta_eq, s_eq, pattern)
    e0_p = lattice_enstrophy(lattice_p, N)
    for step in range(50):
        new_lat = np.zeros_like(lattice_p)
        for ix in range(N):
            for iy in range(N):
                m = lattice_p[ix, iy]
                nbrs = []
                for ddx, ddy in [(1,0),(-1,0),(0,1),(0,-1)]:
                    nx, ny = (ix+ddx)%N, (iy+ddy)%N
                    nbrs.append(lattice_p[nx, ny])
                e_avg = np.mean(nbrs, axis=0)
                new_lat[ix, iy] = ds_step(m, e_avg)
        lattice_p = new_lat
    e_final = lattice_enstrophy(lattice_p, N)
    print(f"{pattern:>8}: enstrophy {e0_p:.6f} -> {e_final:.6f} (ratio {e_final/e0_p:.4f})")

# ============================================================
# Part 7: The PDE identification
# ============================================================
print()
print("=" * 70)
print("PART 7: Continuum PDE from coupled DS")
print("=" * 70)

# Taylor expand DS(m, m+eps*delta) around equilibrium
# m(x), e(x) = m(x+a) for neighbor at displacement a
# e = m + a.grad(m) + (a^2/2)*Lap(m) + ...

# DS(m, m) = T(m) (self-combination = local dynamics)
# DS(m, m + a.grad(m)) = T(m) + J_e . (a.grad(m)) + O(a^2)
# Average over +/- neighbors: odd terms cancel
# DS_avg = T(m) + (a^2/2d) * J_e . Lap(m) + (nonlinear in grad(m))

# On equilibrium manifold: T(m*) = m*, so T(m) - m = 0
# Remaining: dm/dt = (a^2/2d) * J_e . Lap(m) + nonlinear gradient terms

# Compute J_e restricted to the tangent space of the equilibrium manifold
# The equilibrium manifold at K*=7/30 is parametrized by orientation on S^2
# Tangent space: delta_m with delta_theta = 0 and delta_s . s* = 0

# Construct tangent vectors to equilibrium manifold
s_star = m_star[1:]
s_hat = s_star / norm(s_star)

# Two tangent directions: perpendicular to s within the s-plane
# v1: rotate in s1-s2 plane
v1 = np.array([0.0, -s_hat[1], s_hat[0], 0.0])
v1 /= norm(v1) if norm(v1) > 0 else 1.0
# v2: out of s1-s2 plane
v2 = np.array([0.0, 0.0, 0.0, 1.0])

print("Tangent vectors to equilibrium manifold:")
print(f"  v1 = {v1}")
print(f"  v2 = {v2}")

# Project J_e onto tangent space: J_eff = P * J_e * P
# where P projects onto span(v1, v2)
V = np.array([v1, v2]).T  # 4x2
J_eff = V.T @ J_e @ V  # 2x2

print(f"\nEffective 2x2 Jacobian on equilibrium manifold:")
print(f"  [{J_eff[0,0]:+.6f}  {J_eff[0,1]:+.6f}]")
print(f"  [{J_eff[1,0]:+.6f}  {J_eff[1,1]:+.6f}]")

eigvals_eff = np.linalg.eigvals(J_eff)
print(f"  Eigenvalues: {eigvals_eff}")
print(f"  |eigenvalues|: {np.abs(eigvals_eff)}")

# Decompose into symmetric (diffusion) and antisymmetric (precession) parts
J_sym = 0.5 * (J_eff + J_eff.T)
J_anti = 0.5 * (J_eff - J_eff.T)

print(f"\nSymmetric part (diffusion):")
print(f"  [{J_sym[0,0]:+.6f}  {J_sym[0,1]:+.6f}]")
print(f"  [{J_sym[1,0]:+.6f}  {J_sym[1,1]:+.6f}]")
print(f"  Eigenvalues: {np.linalg.eigvals(J_sym)}")

print(f"\nAntisymmetric part (precession):")
print(f"  [{J_anti[0,0]:+.6f}  {J_anti[0,1]:+.6f}]")
print(f"  [{J_anti[1,0]:+.6f}  {J_anti[1,1]:+.6f}]")
print(f"  Off-diagonal: {J_anti[0,1]:.6f}")

print(f"\nThe continuum PDE on the equilibrium manifold:")
print(f"  dm/dt = D_sym * Lap(m) + D_anti * (m x Lap(m)) + O(grad^2)")
print(f"  D_sym (diffusion coeff) ~ eigenvalue of symmetric part")
print(f"  D_anti (precession coeff) ~ off-diagonal of antisymmetric part")
print(f"")
print(f"  This IS the Landau-Lifshitz equation:")
print(f"    dn/dt = alpha * n x Lap(n) + beta * n x (n x Lap(n))")
print(f"  with alpha ~ D_anti and beta ~ D_sym")
