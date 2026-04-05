"""
Multi-site spectral gap via holonomy of the eigenvector frame
=============================================================

The claim: transient amplification on a ring of N sites is bounded
by the holonomy of the eigenvector frame, which is topological
(depends on winding number, not on N).

Strategy:
1. At each site i, compute the Jacobian M_i and its eigenvector frame.
2. Compute the rotation angle between consecutive frames.
3. Sum the rotations around the loop — this is the holonomy.
4. Show the holonomy is independent of N.
5. Show the transient amplification is bounded by the holonomy.
"""

import numpy as np
from numpy.linalg import norm, eig, det, svd

# === Core DS machinery ===

def ds_combine(m, e):
    th_m, s1m, s2m, s3m = m
    th_e, s1e, s2e, s3e = e
    K = s1m*s2e + s1m*s3e + s2m*s1e + s2m*s3e + s3m*s1e + s3m*s2e
    if abs(K) >= 1.0:
        return m.copy()
    inv = 1.0 / (1.0 - K)
    return np.array([
        inv * th_m * th_e,
        inv * (s1m*s1e + s1m*th_e + th_m*s1e),
        inv * (s2m*s2e + s2m*th_e + th_m*s2e),
        inv * (s3m*s3e + s3m*th_e + th_m*s3e),
    ])

def born_floor(m, H=3):
    th = m[0]
    s = m[1:]
    s_sq = np.dot(s, s)
    total_sq = th**2 + s_sq
    if total_sq < 1e-30:
        return np.array([1.0, 0.0, 0.0, 0.0])
    born = th**2 / total_sq
    if born < 1.0/H**3:
        th_new = np.sqrt(s_sq / (H**3 - 1))
        m = np.array([th_new, m[1], m[2], m[3]])
    total = np.sum(np.abs(m))
    if total > 0:
        m = m / total
    return m

def ds_step(m, e, H=3):
    return born_floor(ds_combine(m, e), H)

eps_fd = 1e-7

def compute_jacobian_full(m_base, e_base):
    """4x4 Jacobian d(ds_step)/dm at (m_base, e_base), finite diff."""
    J = np.zeros((4, 4))
    f0 = ds_step(m_base, e_base)
    for j in range(4):
        m_p = m_base.copy()
        m_p[j] += eps_fd
        m_p = np.abs(m_p)
        m_p /= np.sum(m_p)
        f1 = ds_step(m_p, e_base)
        J[:, j] = (f1 - f0) / eps_fd
    return J

def coupled_step_ring(lattice):
    N = len(lattice)
    new = np.zeros_like(lattice)
    for i in range(N):
        e_avg = 0.5 * (lattice[(i-1) % N] + lattice[(i+1) % N])
        new[i] = ds_step(lattice[i], e_avg)
    return new

# ============================================================
# PART 1: Frame rotation along a non-uniform ring
# ============================================================
print("=" * 70)
print("PART 1: Eigenvector frame rotation around a ring")
print("=" * 70)

def frame_angle(V1, V2):
    """Angle between two eigenvector frames (4x4 matrices).
    Measures the rotation needed to align one to the other.
    Uses the principal angles via SVD of V1^T V2."""
    # Align signs: eigenvectors have sign ambiguity
    for j in range(V1.shape[1]):
        if np.dot(V1[:, j], V2[:, j]) < 0:
            V2[:, j] *= -1
    C = V1.T @ V2
    # Principal angles
    s = svd(C, compute_uv=False)
    s = np.clip(s, -1, 1)
    angles = np.arccos(s)
    return np.sum(angles)  # Total rotation

def holonomy_angle(frames):
    """Total frame rotation around a closed loop."""
    N = len(frames)
    total = 0
    for i in range(N):
        total += frame_angle(frames[i].copy(), frames[(i+1) % N].copy())
    return total

np.random.seed(42)

for N in [4, 8, 16, 32, 64]:
    # Create a random non-uniform ring
    lattice = np.array([born_floor(np.random.dirichlet([2, 2, 2, 2])) for _ in range(N)])

    # Compute Jacobian and eigenvector frame at each site
    frames = []
    eigenvalue_products = []
    log_det_sum = 0

    for i in range(N):
        e_avg = 0.5 * (lattice[(i-1) % N] + lattice[(i+1) % N])
        J_i = compute_jacobian_full(lattice[i], e_avg)
        vals, vecs = eig(J_i)

        # Sort by eigenvalue magnitude
        order = np.argsort(-np.abs(vals))
        vals = vals[order]
        vecs = vecs[:, order]

        frames.append(np.real(vecs))
        eigenvalue_products.append(np.prod(np.abs(vals)))
        log_det_sum += np.log(np.abs(det(J_i)) + 1e-30)

    hol = holonomy_angle(frames)
    print(f"N={N:3d}:  holonomy = {hol:.4f}  "
          f"holonomy/N = {hol/N:.4f}  "
          f"log|det product| = {log_det_sum:.4f}  "
          f"log|det|/N = {log_det_sum/N:.4f}")

# ============================================================
# PART 2: Holonomy during transient evolution
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Holonomy during transient approach to equilibrium")
print("=" * 70)

N = 16
np.random.seed(7)
lattice = np.array([born_floor(np.random.dirichlet([2, 2, 2, 2])) for _ in range(N)])

print(f"{'step':>5} {'holonomy':>10} {'hol/N':>8} {'max_amplification':>18} {'spatial_var':>12}")

for step in range(50):
    # Compute frames and measure holonomy
    frames = []
    amplifications = []
    for i in range(N):
        e_avg = 0.5 * (lattice[(i-1) % N] + lattice[(i+1) % N])
        J_i = compute_jacobian_full(lattice[i], e_avg)
        vals, vecs = eig(J_i)
        order = np.argsort(-np.abs(vals))
        frames.append(np.real(vecs[:, order]))
        amplifications.append(np.max(np.abs(vals)))

    hol = holonomy_angle(frames)
    spatial_var = np.std([norm(lattice[i] - np.mean(lattice, axis=0)) for i in range(N)])

    # Measure transient amplification: apply perturbation, propagate full loop
    delta = np.random.randn(4) * 0.001
    delta -= np.mean(delta)
    total_amp = 1.0
    d = delta.copy()
    for i in range(N):
        e_avg = 0.5 * (lattice[(i-1) % N] + lattice[(i+1) % N])
        J_i = compute_jacobian_full(lattice[i], e_avg)
        d = J_i @ d
    amp = norm(d) / (norm(delta) + 1e-30)

    if step % 3 == 0:
        print(f"{step:5d} {hol:10.4f} {hol/N:8.4f} {amp:18.4f} {spatial_var:12.6f}")

    lattice = coupled_step_ring(lattice)

# ============================================================
# PART 3: The spectral radius of the holonomy matrix
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Spectral radius of the loop transfer matrix M_N...M_2 M_1")
print("=" * 70)

# The actual quantity that controls perturbation decay:
# Propagate a perturbation all the way around the ring.
# The matrix is M_1 * M_2 * ... * M_N.
# Its spectral radius must be < 1 for loop contraction.

for N in [4, 8, 16, 32]:
    np.random.seed(42)
    lattice = np.array([born_floor(np.random.dirichlet([2, 2, 2, 2])) for _ in range(N)])

    # Evolve for a few steps so it's not totally random
    for _ in range(5):
        lattice = coupled_step_ring(lattice)

    # Compute product of Jacobians around the ring
    M_product = np.eye(4)
    for i in range(N):
        e_avg = 0.5 * (lattice[(i-1) % N] + lattice[(i+1) % N])
        J_i = compute_jacobian_full(lattice[i], e_avg)
        M_product = J_i @ M_product

    vals_product = np.abs(eig(M_product)[0])
    rho = np.max(vals_product)
    rho_per_site = rho ** (1.0/N)

    print(f"N={N:3d}:  ρ(M_1...M_N) = {rho:.6e}  "
          f"ρ^(1/N) = {rho_per_site:.6f}  "
          f"|det|^(1/N) = {np.abs(det(M_product))**(1.0/N):.6e}")

# At equilibrium
print("\nAt uniform equilibrium:")
for N in [4, 8, 16, 32, 64, 128]:
    # All sites equal
    m_eq = born_floor(np.array([0.25, 0.25, 0.25, 0.25]))
    for _ in range(100):
        m_eq = ds_step(m_eq, m_eq)

    J_eq = compute_jacobian_full(m_eq, m_eq)
    vals_eq = np.abs(eig(J_eq)[0])
    rho_single = np.max(vals_eq)

    # Product = J^N
    rho_N = rho_single ** N
    rho_per = rho_N ** (1.0/N)

    print(f"N={N:4d}:  ρ(J) = {rho_single:.6f}  "
          f"ρ(J^N) = {rho_N:.6e}  "
          f"ρ(J^N)^(1/N) = {rho_per:.6f}")

# ============================================================
# PART 4: Connection curvature along the ring
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Connection curvature (frame rotation rate)")
print("=" * 70)

# The curvature of the eigenvector frame connection:
# At each link (i, i+1), measure the rotation of the eigenvector frame.
# The integral around the loop = 2πn (winding number) by Gauss-Bonnet.
# This constrains how much surfing a perturbation can do.

for N in [4, 8, 16, 32, 64]:
    np.random.seed(42)
    lattice = np.array([born_floor(np.random.dirichlet([2, 2, 2, 2])) for _ in range(N)])

    # Evolve 5 steps
    for _ in range(5):
        lattice = coupled_step_ring(lattice)

    frames = []
    for i in range(N):
        e_avg = 0.5 * (lattice[(i-1) % N] + lattice[(i+1) % N])
        J_i = compute_jacobian_full(lattice[i], e_avg)
        vals, vecs = eig(J_i)
        order = np.argsort(-np.abs(vals))
        frames.append(np.real(vecs[:, order]))

    # Link rotations
    link_angles = []
    for i in range(N):
        angle = frame_angle(frames[i].copy(), frames[(i+1) % N].copy())
        link_angles.append(angle)

    total = sum(link_angles)
    print(f"N={N:3d}:  total_rotation = {total:.4f}  "
          f"per_link = {total/N:.4f}  "
          f"total/(2π) = {total/(2*np.pi):.4f}  "
          f"max_link = {max(link_angles):.4f}  "
          f"min_link = {min(link_angles):.4f}")

# ============================================================
# PART 5: The key test — transient amplification vs N
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Maximum transient amplification as a function of N")
print("=" * 70)

# If the holonomy argument works, max amplification should be bounded
# independently of N.

for N in [4, 8, 16, 32, 64]:
    max_amp_over_trials = 0
    for trial in range(50):
        np.random.seed(trial * 1000 + N)
        lat_a = np.array([born_floor(np.random.dirichlet([2, 2, 2, 2])) for _ in range(N)])
        lat_b = np.array([born_floor(np.random.dirichlet([2, 2, 2, 2])) for _ in range(N)])

        d0 = np.sqrt(sum(norm(lat_a[i] - lat_b[i])**2 for i in range(N)))
        max_d = d0

        for step in range(100):
            lat_a = coupled_step_ring(lat_a)
            lat_b = coupled_step_ring(lat_b)
            d = np.sqrt(sum(norm(lat_a[i] - lat_b[i])**2 for i in range(N)))
            max_d = max(max_d, d)

        amp = max_d / (d0 + 1e-30)
        max_amp_over_trials = max(max_amp_over_trials, amp)

    print(f"N={N:3d}:  max transient amplification = {max_amp_over_trials:.2f}")

# ============================================================
# PART 6: Convergence time vs N
# ============================================================
print("\n" + "=" * 70)
print("PART 6: Convergence time (steps to within 1% of equilibrium) vs N")
print("=" * 70)

for N in [4, 8, 16, 32, 64]:
    conv_times = []
    for trial in range(30):
        np.random.seed(trial * 100 + N)
        lat = np.array([born_floor(np.random.dirichlet([2, 2, 2, 2])) for _ in range(N)])

        # Find this ring's equilibrium by running a long time
        lat_eq = lat.copy()
        for _ in range(500):
            lat_eq = coupled_step_ring(lat_eq)
        m_eq_local = lat_eq[0]

        # Now find how long original lattice takes to get close
        lat2 = np.array([born_floor(np.random.dirichlet([2, 2, 2, 2])) for _ in range(N)])
        converged_at = 500
        for step in range(500):
            lat2 = coupled_step_ring(lat2)
            max_d = max(norm(lat2[i] - m_eq_local) for i in range(N))
            if max_d < 0.01:
                converged_at = step
                break

        conv_times.append(converged_at)

    print(f"N={N:3d}:  mean convergence time = {np.mean(conv_times):.1f}  "
          f"max = {np.max(conv_times)}  "
          f"min = {np.min(conv_times)}")
