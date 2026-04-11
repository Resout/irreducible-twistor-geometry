#!/usr/bin/env python3
"""
Born Floor Preservation Under Navier-Stokes Stretching
=======================================================

THE OPEN QUESTION (paper line 1903):
  Does diffusion dominate stretching at |omega| = sqrt(26)?

The DS descended equation preserves Born(m) >= 1/27 (proved).
The NS stretching term (omega . nabla)u is ABSENT from the DS equation
(0.0% Levi-Civita content). When we embed the DS mass function into a
fluid via the minitwistor construction, the stretching term re-enters.

This script investigates whether the Born floor survives.

Key identity (eq:born_omega in the paper):
    |omega|^2 = |s|^2 / theta^2 = (1 - Born) / Born

So Born >= 1/27  <=>  |omega|^2 <= 26  <=>  |omega| <= sqrt(26).

Framework constants:
    H = 3
    K* = 7/30
    lambda_0 = 0.2829  (leading eigenvalue of single-site transfer operator)
    Delta = -ln(lambda_0) = 1.2626  (spectral gap / decay rate)
    C(H) = 0.344  (curvature bound: max ||dbar Phi|| on compact B_Z)
    sqrt(26) ~ 5.099  (max vorticity from Born floor)

NOTE: C(H) = 0.344 is the curvature bound on the anti-holomorphic
Jacobian dbar(Phi), NOT the vorticity bound. The vorticity bound is
sqrt(26) ~ 5.099. The paper is explicit about this distinction (line 1911).
"""

import numpy as np
from scipy.linalg import expm
import warnings
warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════
# PART 0: Framework Constants
# ═══════════════════════════════════════════════════════════════════

H = 3
Kstar = 7 / 30
Born_floor = 1 / H**3  # = 1/27
lambda_0 = 0.2829
lambda_1 = 0.2813
Delta = -np.log(lambda_0)  # spectral gap
C_H = 0.344  # curvature bound (gauge-sector norm of dbar Phi)
omega_max = np.sqrt(H**3 - 1)  # = sqrt(26) ~ 5.099, vorticity bound

print("=" * 72)
print("BORN FLOOR PRESERVATION UNDER NS STRETCHING")
print("=" * 72)
print(f"\nFramework constants:")
print(f"  H             = {H}")
print(f"  K*            = {Kstar:.10f} = 7/30")
print(f"  Born floor    = 1/{H**3} = {Born_floor:.10f}")
print(f"  lambda_0      = {lambda_0}")
print(f"  Delta         = -ln(lambda_0) = {Delta:.4f}")
print(f"  C(H)          = {C_H} (curvature bound on dbar Phi)")
print(f"  omega_max     = sqrt({H**3 - 1}) = {omega_max:.6f}")

# ═══════════════════════════════════════════════════════════════════
# PART 1: Rate Comparison (Enforcement vs Stretching)
# ═══════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("PART 1: RATE COMPARISON")
print("=" * 72)

# The DS floor enforcement rate: how fast does the DS dynamics push
# a state back to Born >= 1/27 when it's at the boundary?
#
# The single-site contraction rate is lambda_0 = 0.2829 per step.
# Perturbations decay as delta(t) ~ delta(0) * lambda_0^t.
# The enforcement rate per step is 1 - lambda_0 = 0.717.
# As a continuous rate: Delta = -ln(lambda_0) = 1.2626.

enforcement_rate_discrete = 1 - lambda_0
enforcement_rate_continuous = Delta

print(f"\nDS floor enforcement rate:")
print(f"  Discrete: 1 - lambda_0 = {enforcement_rate_discrete:.4f} per step")
print(f"  Continuous: Delta = {enforcement_rate_continuous:.4f} per unit time")

# The NS stretching rate: maximum rate at which (omega . nabla)u can
# amplify vorticity.
#
# The stretching term (omega . nabla)u has magnitude bounded by
# |omega| * |nabla u|.
#
# From Biot-Savart: u = K * omega where K is the BS kernel.
# |nabla u| is bounded by |omega| (Calderon-Zygmund, up to constants).
#
# So the maximum stretching rate is proportional to |omega|^2.
# At the Born floor: |omega| = sqrt(26), so |omega|^2 = 26.
#
# But this is the DIMENSIONAL stretching rate. We need to compare
# it to the diffusion rate nu * |nabla^2 omega| / |omega|.

# For the pointwise maximum principle argument:
# d/dt |omega|^2 <= 2 * |omega| * |S| * |omega| - 2*nu*|nabla omega|^2
# where S is the strain rate tensor.
#
# At the Born boundary (|omega|^2 = 26):
# The stretching term: 2 * omega_i * S_ij * omega_j <= 2 * |omega|^2 * |S|
# The strain bound: |S| <= C_d * |omega| (from Calderon-Zygmund)
# So stretching <= 2 * C_d * |omega|^3

# However, the KEY insight from the framework is:
# The Born floor is enforced ALGEBRAICALLY by the fibre structure.
# It's not a dynamical balance -- it's a constraint.
# The question is whether NS dynamics can violate the constraint.

# Let's compute the stretching-to-enforcement ratio at the boundary.
# Using the paper's C(H) for the Jacobian norm:

stretching_rate_CH = C_H**2
ratio_CH = enforcement_rate_continuous / stretching_rate_CH

print(f"\nStretching rate using C(H) = {C_H} (Jacobian curvature bound):")
print(f"  C(H)^2           = {stretching_rate_CH:.6f}")
print(f"  Enforcement/Stretch = {ratio_CH:.2f}")
print(f"  -> Floor enforcement is {ratio_CH:.1f}x faster")

# Using the vorticity bound sqrt(26):
# This is the relevant bound for the NS problem
stretching_rate_omega = omega_max**2  # = 26
ratio_omega = enforcement_rate_continuous / stretching_rate_omega

print(f"\nStretching rate using omega_max = sqrt(26) (vorticity bound):")
print(f"  omega_max^2       = {stretching_rate_omega:.1f}")
print(f"  Enforcement/Stretch = {ratio_omega:.4f}")
print(f"  -> Stretching is {1/ratio_omega:.1f}x faster than enforcement!")

print(f"\nCRITICAL OBSERVATION:")
print(f"  The C(H) = 0.344 bound applies to the GAUGE SECTOR (dbar Phi),")
print(f"  not to the vorticity. The vorticity bound is sqrt(26) ~ 5.099.")
print(f"  At the Born floor, stretching rate (26) >> enforcement rate (1.26).")
print(f"  This means the simple rate comparison does NOT guarantee preservation.")
print(f"  The question is genuinely open.")

# ═══════════════════════════════════════════════════════════════════
# PART 2: Maximum Principle Analysis
# ═══════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("PART 2: MAXIMUM PRINCIPLE ANALYSIS")
print("=" * 72)

# Define Psi(m) = Born(m) - 1/27
# Born(m) = theta^2 / (s1^2 + s2^2 + s3^2 + theta^2)
# Let r^2 = s1^2 + s2^2 + s3^2.  Then Born = theta^2 / (r^2 + theta^2).
# Born = 1/27 when 27*theta^2 = r^2 + theta^2, i.e. r^2 = 26*theta^2.
#
# The Born floor boundary is the sphere r = sqrt(26) * theta in mass space.
# This is a CONE in (s1, s2, s3, theta) space.
#
# d(Born)/dt = d/dt [theta^2 / (r^2 + theta^2)]
#            = [2*theta*dtheta/dt * (r^2 + theta^2) - theta^2 * (2r*dr/dt + 2*theta*dtheta/dt)]
#              / (r^2 + theta^2)^2
#            = [2*theta*dtheta/dt * r^2 - 2*theta^2 * r * dr/dt] / (r^2 + theta^2)^2
#
# On the boundary Born = 1/27 (r^2 = 26*theta^2):
# d(Born)/dt = 2*theta / (27*theta^2)^2 * [26*theta^2 * dtheta/dt - theta^2 * sqrt(26)*theta * dr/dt]
#            (simplified for radial case, r = |s|)

# The NS stretching term (omega . nabla)u affects the vorticity omega = s/theta.
# Under stretching:
#   d(omega)/dt |_stretch = (omega . nabla)u
# This changes s and theta indirectly through the Biot-Savart reconstruction.

# Key structural point: the stretching term preserves the MAGNITUDE of
# vorticity only when omega is aligned with the eigenvector of the strain
# tensor. In general:
#   d|omega|^2/dt |_stretch = 2 * omega_i * S_ij * omega_j
# where S_ij = (d_i u_j + d_j u_i)/2 is the strain rate tensor.

# This can be positive (amplification) or negative (compression).
# The Born floor is violated when |omega|^2 exceeds 26.

# For a maximum principle: we need d|omega|^2/dt <= 0 whenever |omega|^2 = 26.
# The full evolution:
#   d|omega|^2/dt = 2*omega_i*S_ij*omega_j + nu * (nabla^2 |omega|^2 - 2|nabla omega|^2)

# On the boundary |omega|^2 = 26:
# The diffusion contribution to the pointwise maximum of |omega|^2 is <= 0
# (by the maximum principle for the heat equation: nabla^2|omega|^2 <= 0 at a max).
# But omega does NOT satisfy the heat equation -- it satisfies NS.

print("\nMaximum principle on the Born boundary |omega|^2 = 26:")
print()
print("  d|omega|^2/dt = 2*omega_i*S_ij*omega_j  [stretching]")
print("               + nu*(nabla^2|omega|^2 - 2|nabla omega|^2)  [diffusion]")
print()
print("At a spatial maximum of |omega|^2:")
print("  nabla(|omega|^2) = 0  (critical point)")
print("  nabla^2(|omega|^2) <= 0  (maximum => negative semidefinite Hessian)")
print("  So the diffusion term <= -2*nu*|nabla omega|^2 <= 0")
print()
print("The diffusion term ALWAYS opposes growth at a maximum.")
print("The stretching term 2*omega_i*S_ij*omega_j can be positive.")
print()
print("Question: can stretching exceed diffusion at |omega|^2 = 26?")

# ═══════════════════════════════════════════════════════════════════
# PART 3: Numerical Simulation -- DS Dynamics with Stretching
# ═══════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("PART 3: NUMERICAL SIMULATION")
print("=" * 72)

# We simulate on a 1D periodic grid.
# Mass function m = (s1, s2, s3, theta) at each grid point.
# DS dynamics + an added stretching perturbation proportional to
# the Levi-Civita structure.

def born_prob(m):
    """Born probability: theta^2 / |m|^2."""
    s = m[:3]
    theta = m[3]
    return theta**2 / (np.sum(s**2) + theta**2)

def born_prob_field(field):
    """Born probability at each grid point. field shape: (N, 4)."""
    s_sq = np.sum(field[:, :3]**2, axis=1)
    theta_sq = field[:, 3]**2
    return theta_sq / (s_sq + theta_sq)

def ds_combine(m1, m2, K):
    """Dempster-Shafer combination of two mass functions.

    DS rule: m_combined = (m1 * m2) / (1 - K)
    where K is the conflict and * is the pointwise product on
    the mass components, properly normalized.

    For simplicity, we use the linearized version around equilibrium.
    """
    # Product in mass space (component-wise for singletons, special for theta)
    s1, s2, s3, theta1 = m1[0], m1[1], m1[2], m1[3]
    e1, e2, e3, theta2 = m2[0], m2[1], m2[2], m2[3]

    # DS combination (unnormalized)
    out_s = np.array([
        s1*e1 + s1*theta2 + theta1*e1,
        s2*e2 + s2*theta2 + theta1*e2,
        s3*e3 + s3*theta2 + theta1*e3,
    ])
    out_theta = theta1 * theta2

    # Conflict
    conflict = 0.0
    for i in range(3):
        for j in range(3):
            if i != j:
                conflict += m1[i] * m2[j]

    # Normalize
    total = np.sum(out_s) + out_theta
    if total < 1e-15:
        return m1.copy()

    out = np.zeros(4)
    out[:3] = out_s / total
    out[3] = out_theta / total

    # Enforce Born floor
    out = enforce_born_floor(out)
    return out

def enforce_born_floor(m):
    """Project back to Born >= 1/27 if violated."""
    b = born_prob(m)
    if b < Born_floor:
        # Scale s down to satisfy 26*theta^2 = |s|^2
        s = m[:3]
        theta = m[3]
        s_sq = np.sum(s**2)
        if s_sq > 0 and theta > 0:
            max_s_sq = 26 * theta**2
            scale = np.sqrt(max_s_sq / s_sq)
            m = m.copy()
            m[:3] = s * scale
    return m

def equilibrium_mass():
    """Return the K* = 7/30 equilibrium mass function."""
    # At equilibrium: Born = 1/27, all s_i equal
    # theta^2 / (3*s^2 + theta^2) = 1/27
    # 27*theta^2 = 3*s^2 + theta^2
    # 26*theta^2 = 3*s^2
    # s^2 = 26/3 * theta^2
    # Normalize: 3*s^2 + theta^2 = 1
    # 3*(26/3)*theta^2 + theta^2 = 27*theta^2 = 1
    # theta^2 = 1/27, theta = 1/sqrt(27)
    theta = 1.0 / np.sqrt(27.0)
    s = np.sqrt(26.0 / 3.0) * theta
    return np.array([s, s, s, theta])

def apply_ds_step(field, K):
    """Apply one DS combination step: each site combines with its neighbors."""
    N = field.shape[0]
    new_field = np.zeros_like(field)
    for i in range(N):
        # Average of left and right neighbors as evidence
        left = field[(i - 1) % N]
        right = field[(i + 1) % N]
        evidence = 0.5 * (left + right)
        # Normalize evidence
        evidence = evidence / np.sum(np.abs(evidence))
        new_field[i] = ds_combine(field[i], evidence, K)
    return new_field

def apply_ns_stretching(field, dx, epsilon):
    """Apply a mock NS stretching perturbation.

    The NS stretching term (omega . nabla)u involves the Levi-Civita
    structure. We model this as:

    delta_s_i += epsilon * sum_jk eps_ijk * omega_j * (nabla s_k)

    where omega = s/theta and nabla is the discrete gradient.
    This adds the antisymmetric cross-product structure that is
    ABSENT from the DS equation.
    """
    N = field.shape[0]
    new_field = field.copy()
    eps_tensor = np.zeros((3, 3, 3))
    eps_tensor[0, 1, 2] = 1
    eps_tensor[1, 2, 0] = 1
    eps_tensor[2, 0, 1] = 1
    eps_tensor[0, 2, 1] = -1
    eps_tensor[2, 1, 0] = -1
    eps_tensor[1, 0, 2] = -1

    for i in range(N):
        s = field[i, :3]
        theta = field[i, 3]
        if theta < 1e-15:
            continue
        omega = s / theta  # vorticity

        # Discrete gradient of s
        grad_s = (field[(i + 1) % N, :3] - field[(i - 1) % N, :3]) / (2 * dx)

        # Stretching: delta_s += epsilon * (omega x grad_s) in each component
        # Using Levi-Civita: delta_s_a = epsilon * sum_bc eps_abc * omega_b * grad_s_c
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    new_field[i, a] += epsilon * eps_tensor[a, b, c] * omega[b] * grad_s[c]

    return new_field

def apply_diffusion(field, dx, nu):
    """Apply discrete Laplacian diffusion."""
    N = field.shape[0]
    new_field = field.copy()
    for i in range(N):
        laplacian = (field[(i + 1) % N] - 2 * field[i] + field[(i - 1) % N]) / dx**2
        new_field[i] += nu * laplacian
    return new_field

# --- Simulation Parameters ---
N_grid = 64
dx = 1.0 / N_grid
dt = 0.001
nu = 0.01  # viscosity
n_steps = 2000

# --- Initialize near the Born floor ---
m_eq = equilibrium_mass()
print(f"\nEquilibrium mass: s = {m_eq[0]:.6f}, theta = {m_eq[3]:.6f}")
print(f"Born at equilibrium: {born_prob(m_eq):.10f} (should be {Born_floor:.10f})")

# Create a field with spatial variation
np.random.seed(42)
field = np.zeros((N_grid, 4))
for i in range(N_grid):
    # Perturb around equilibrium; push some regions near the Born floor
    pert = 0.15 * np.random.randn(4)
    # Make the perturbation mostly in s-direction (push toward floor)
    pert[3] *= 0.3  # smaller theta perturbation
    field[i] = m_eq + pert
    # Ensure positivity
    field[i] = np.abs(field[i])
    # Normalize
    field[i] = field[i] / np.sum(field[i])

# Enforce Born floor on initial data
for i in range(N_grid):
    field[i] = enforce_born_floor(field[i])

born_init = born_prob_field(field)
print(f"\nInitial Born probabilities:")
print(f"  min = {born_init.min():.6f}, max = {born_init.max():.6f}")
print(f"  mean = {born_init.mean():.6f}")

# --- Run 4 scenarios ---
scenarios = {
    "DS only (no stretching)": {"ds": True, "stretch": False, "diffuse": True, "eps": 0.0},
    "DS + weak stretching (eps=0.01)": {"ds": True, "stretch": True, "diffuse": True, "eps": 0.01},
    "DS + moderate stretching (eps=0.1)": {"ds": True, "stretch": True, "diffuse": True, "eps": 0.1},
    "DS + strong stretching (eps=1.0)": {"ds": True, "stretch": True, "diffuse": True, "eps": 1.0},
}

results = {}

for name, params in scenarios.items():
    print(f"\n--- {name} ---")
    f = field.copy()
    min_borns = []
    mean_borns = []
    max_omegas = []
    floor_violations = []

    for step in range(n_steps):
        # Record diagnostics
        borns = born_prob_field(f)
        min_borns.append(borns.min())
        mean_borns.append(borns.mean())

        s_sq = np.sum(f[:, :3]**2, axis=1)
        theta_sq = f[:, 3]**2
        omega_sq = s_sq / np.maximum(theta_sq, 1e-30)
        max_omegas.append(np.sqrt(omega_sq.max()))

        violations = np.sum(borns < Born_floor - 1e-10)
        floor_violations.append(violations)

        # Apply dynamics
        if params["diffuse"]:
            f = apply_diffusion(f, dx, nu * dt)

        if params["stretch"]:
            f = apply_ns_stretching(f, dx, params["eps"] * dt)

        if params["ds"]:
            # Enforce Born floor (this is the DS floor enforcement)
            for i in range(N_grid):
                f[i] = enforce_born_floor(f[i])
            # Normalize
            for i in range(N_grid):
                total = np.sum(np.abs(f[i]))
                if total > 0:
                    f[i] = np.abs(f[i]) / total

    borns_final = born_prob_field(f)
    results[name] = {
        "min_borns": np.array(min_borns),
        "mean_borns": np.array(mean_borns),
        "max_omegas": np.array(max_omegas),
        "floor_violations": np.array(floor_violations),
    }

    print(f"  Final Born: min={borns_final.min():.6f}, mean={borns_final.mean():.6f}")
    print(f"  Max omega over run: {np.max(max_omegas):.4f} (bound: {omega_max:.4f})")
    print(f"  Floor violations at any step: {np.max(floor_violations)}")
    print(f"  Min Born ever: {np.min(min_borns):.6f} (floor: {Born_floor:.6f})")

# ═══════════════════════════════════════════════════════════════════
# PART 4: Stretching Without Floor Enforcement
# ═══════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("PART 4: STRETCHING WITHOUT BORN FLOOR ENFORCEMENT")
print("(What happens if the floor is NOT enforced?)")
print("=" * 72)

# Run the same simulation but WITHOUT enforcing the Born floor.
# This shows what stretching alone does.

eps_values = [0.01, 0.1, 0.5, 1.0]
for eps in eps_values:
    f = field.copy()
    min_born_no_floor = 1.0
    violated_steps = 0
    max_omega_no_floor = 0.0

    for step in range(n_steps):
        borns = born_prob_field(f)
        min_born_no_floor = min(min_born_no_floor, borns.min())
        if borns.min() < Born_floor - 1e-10:
            violated_steps += 1

        s_sq = np.sum(f[:, :3]**2, axis=1)
        theta_sq = f[:, 3]**2
        omega_sq = s_sq / np.maximum(theta_sq, 1e-30)
        max_omega_no_floor = max(max_omega_no_floor, np.sqrt(omega_sq.max()))

        # Diffusion + stretching, NO floor enforcement
        f = apply_diffusion(f, dx, nu * dt)
        f = apply_ns_stretching(f, dx, eps * dt)

        # Only normalize, do NOT enforce Born floor
        for i in range(N_grid):
            f[i] = np.abs(f[i])
            total = np.sum(f[i])
            if total > 0:
                f[i] = f[i] / total

    print(f"\n  eps = {eps}:")
    print(f"    Min Born: {min_born_no_floor:.6f} (floor: {Born_floor:.6f})")
    print(f"    Max omega: {max_omega_no_floor:.4f} (bound: {omega_max:.4f})")
    print(f"    Steps with violation: {violated_steps}/{n_steps}")

# ═══════════════════════════════════════════════════════════════════
# PART 5: The BKM Bound with Born Floor
# ═══════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("PART 5: BKM BOUND ANALYSIS")
print("=" * 72)

# The Beale-Kato-Majda criterion: NS blows up at time T iff
#   int_0^T ||omega||_infty dt = infinity
#
# If Born >= 1/27 is preserved, then ||omega||_infty <= sqrt(26) for all t.
# Therefore: int_0^T ||omega||_infty dt <= sqrt(26) * T < infinity for all T.
# => NO blowup. QED.
#
# The entire NS regularity question reduces to: is Born >= 1/27 preserved?

print(f"\nBKM criterion: blowup iff int_0^T ||omega||_infty dt = infinity")
print(f"\nIF Born >= 1/27 is preserved:")
print(f"  ||omega||_infty <= sqrt(26) = {omega_max:.6f} for all t")
print(f"  int_0^T ||omega||_infty dt <= sqrt(26) * T < infinity")
print(f"  => NO BLOWUP for any finite T")
print(f"  => GLOBAL REGULARITY")

# Exponential growth scenario (worst case without Born floor):
print(f"\nWithout Born floor (worst case exponential growth):")
print(f"  d|omega|/dt <= |S| * |omega| <= C * |omega|^2")
print(f"  This gives |omega|(t) ~ 1/(T* - t)  (finite-time blowup)")
print(f"  BKM integral diverges logarithmically at T*")

# With Born floor:
print(f"\nWith Born floor (framework prediction):")
print(f"  The algebraic identity |omega|^2 = (1 - Born)/Born")
print(f"  constrains |omega|^2 <= 26 whenever Born >= 1/27.")
print(f"  The question is whether NS dynamics can push Born below 1/27.")

# ═══════════════════════════════════════════════════════════════════
# PART 6: Energy Estimate on the Born Boundary
# ═══════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("PART 6: ENERGY ESTIMATE ON THE BORN BOUNDARY")
print("=" * 72)

# On the Born boundary, Born = 1/27, i.e. |omega|^2 = 26.
# The evolution of |omega|^2 under NS:
#
#   d/dt |omega|^2 = 2 * omega_i * S_ij * omega_j   [stretching]
#                  + nu * Delta(|omega|^2)            [diffusion of |omega|^2]
#                  - 2*nu * |nabla omega|^2           [enstrophy dissipation]
#
# At a spatial maximum of |omega|^2:
#   Delta(|omega|^2) <= 0  (it's a maximum!)
#   |nabla omega|^2 >= 0
#
# So at the maximum:
#   d/dt |omega|^2_max <= 2 * omega_i * S_ij * omega_j
#
# The strain-vorticity alignment angle phi determines the sign:
#   omega_i * S_ij * omega_j = |omega|^2 * lambda_max(S) * cos^2(phi)
#                              ... (sum over all eigenvalues)
#
# For the stretching to push |omega|^2 above 26, we need:
#   omega_i * S_ij * omega_j > 0 at |omega|^2 = 26
#
# This IS possible in general NS flows. The strain eigenvector can
# align with the vorticity to produce amplification.

# Compute the maximum stretching rate for various vorticity levels
omega_values = np.linspace(0, omega_max, 100)
born_values = 1.0 / (1.0 + omega_values**2 / 1.0)  # wrong
born_values = 1.0 / (1.0 + omega_values**2)  # Born = theta^2/(s^2 + theta^2) with |s| = |omega|*theta
# Actually: Born = 1/(1 + |omega|^2) when we set theta = 1 for normalization
# Wait -- |omega|^2 = (1 - Born)/Born, so Born = 1/(1 + |omega|^2)

born_from_omega = 1.0 / (1.0 + omega_values**2)

print(f"\nBorn probability as a function of vorticity magnitude:")
print(f"  |omega| = 0     -> Born = {1.0/(1+0):.6f}")
print(f"  |omega| = 1     -> Born = {1.0/(1+1):.6f}")
print(f"  |omega| = 3     -> Born = {1.0/(1+9):.6f}")
print(f"  |omega| = sqrt(26) -> Born = {1.0/(1+26):.6f} = 1/27")

# The stretching term magnitude at the boundary:
# |2 * omega . S . omega| <= 2 * |omega|^2 * |S|_max
# where |S|_max is the largest eigenvalue of the strain tensor.
#
# From Calderon-Zygmund theory: |S|_max ~ |omega|_max
# (the strain is comparable to the vorticity in L^p norms)
#
# So maximum stretching rate at Born boundary:
# 2 * 26 * sqrt(26) ~ 265

max_stretch_at_boundary = 2 * 26 * np.sqrt(26)
print(f"\nMaximum stretching rate at Born boundary:")
print(f"  2 * |omega|^2 * |S|_max ~ 2 * 26 * sqrt(26) = {max_stretch_at_boundary:.1f}")
print(f"  (dimensionless; depends on spatial gradients)")

# The diffusion term at the boundary:
# -2*nu*|nabla omega|^2 is always <= 0 (dissipative)
# Its magnitude depends on the spatial structure.
#
# For a vortex tube of radius r with |omega| = sqrt(26):
# |nabla omega| ~ sqrt(26)/r
# Diffusion: -2*nu*26/r^2
# Stretching: 2*26*sqrt(26)/r (if strain ~ omega/r)
# Ratio: stretching/diffusion = sqrt(26)*r / nu
# Critical radius: r_c = nu/sqrt(26)
#
# Below r_c, diffusion dominates. Above r_c, stretching dominates.

nu_phys = 1.0  # kinematic viscosity (dimensionless)
r_critical = nu_phys / omega_max
print(f"\nCritical vortex tube radius (nu = {nu_phys}):")
print(f"  r_c = nu / sqrt(26) = {r_critical:.6f}")
print(f"  Below r_c: diffusion dominates stretching")
print(f"  Above r_c: stretching can potentially dominate")
print(f"  This is the standard NS regularity barrier.")

# ═══════════════════════════════════════════════════════════════════
# PART 7: The Sufficient Condition
# ═══════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("PART 7: SUFFICIENT CONDITION FOR BORN FLOOR PRESERVATION")
print("=" * 72)

# The question: does the O(2) bundle structure constrain the strain tensor?
#
# In the minitwistor construction:
#   omega lives in H^1(Z, O(-3))  (Hitchin correspondence)
#   u is reconstructed from omega via Biot-Savart
#   S = (nabla u + nabla u^T)/2 is the strain tensor
#
# Key observation: u and omega are not independent. They are related by
# the Biot-Savart law: u = curl^{-1}(omega) = K * omega.
#
# The alignment cos^2(phi) between omega and the principal strain direction
# is constrained by this relationship.
#
# Tsinober-Galanti observation (DNS): in turbulent flows,
# cos^2(phi) ~ 0.5-0.7 (vorticity preferentially aligns with
# intermediate strain eigenvector, not the maximum).
#
# But this is a statistical observation, not a pointwise bound.

# Framework-specific constraint:
# If omega has a minitwistor representative in B_Z, then the
# Penrose transform constrains the spatial derivatives.
# Specifically: nabla omega = dbar(something on Z), and the
# Born floor constrains the Z-data.

# The critical computation: what is the maximum of
#   omega_i * S_ij * omega_j / |omega|^2
# subject to:
#   div(u) = 0
#   curl(u) = omega
#   Born(omega) = 1/27  (on the boundary)

# This is a constrained optimization problem.
# S_ij = (d_i u_j + d_j u_i)/2  where  u = K * omega.

# For an axisymmetric vortex with |omega| = omega_max:
# S has eigenvalues that depend on the profile.

# Burgers vortex (exact NS solution):
# omega = (0, 0, omega_z(r)) with omega_z = Gamma/(4*pi*nu) * exp(-alpha*r^2/(4*nu))
# Strain: S_rr = -alpha/2, S_zz = alpha (externally imposed)
# Stretching: omega_z * S_zz = alpha * omega_z (amplifies vorticity)
# Balance: alpha * omega = nu * nabla^2 omega => steady state
#
# At steady state, the vorticity is bounded. BUT the external strain alpha
# is an imposed parameter, not self-generated.

# For self-generated strain (no external forcing), the Biot-Savart
# constraint is crucial.

print(f"\nSelf-generated strain constraint (Biot-Savart):")
print(f"  S = symmetric part of nabla(K * omega)")
print(f"  For omega with ||omega||_infty = sqrt(26):")
print()

# Biot-Savart in 3D: u(x) = -1/(4*pi) * integral omega(y) x (x-y)/|x-y|^3 dy
# nabla u is a singular integral (Calderon-Zygmund operator)
# ||S||_infty <= C_CZ * ||omega||_infty  where C_CZ is universal

# The CZ constant for the Biot-Savart gradient:
# In 3D, C_CZ is related to the Riesz transform bound.
# The sharp L^infty -> L^infty bound is logarithmically divergent
# (it doesn't hold in L^infty!). But in BMO:
# ||S||_BMO <= C * ||omega||_infty

# This is WHY NS regularity is hard: you can't bound ||nabla u||_infty
# directly from ||omega||_infty. You need ||omega||_infty control
# in time, which is exactly what BKM gives IF you have it.

print(f"  KEY DIFFICULTY: The CZ bound ||S||_infty <= C * ||omega||_infty")
print(f"  does NOT hold in L^infty (logarithmic divergence).")
print(f"  This is the standard NS regularity barrier.")
print(f"  The Born floor gives ||omega||_infty <= sqrt(26), but this")
print(f"  does NOT directly bound ||S||_infty.")
print()
print(f"  HOWEVER: if Born >= 1/27 is preserved (the assumption),")
print(f"  then ||omega||_infty stays bounded, and BKM gives regularity.")
print(f"  The question is circular: preservation implies regularity")
print(f"  implies preservation. Breaking the circle requires a new idea.")

# ═══════════════════════════════════════════════════════════════════
# PART 8: The Framework's Structural Advantage
# ═══════════════════════════════════════════════════════════════════

print("\n" + "=" * 72)
print("PART 8: WHAT THE FRAMEWORK ACTUALLY PROVIDES")
print("=" * 72)

print(f"""
The framework does NOT solve the NS Millennium Problem. Here is what it
provides and what it does not:

PROVED (unconditional):
  1. The DS descended equation has global regularity (Thm regularity).
     This is a genuine result for a specific PDE, not NS.
  2. The Born floor Born >= 1/27 is preserved under DS dynamics.
  3. The descended equation is NOT NS (0% Levi-Civita content).
  4. IF Born >= 1/27 is preserved under NS, THEN NS is regular.
     (BKM criterion + algebraic identity |omega|^2 = (1-B)/B.)

CONDITIONAL:
  5. Born floor preservation under NS evolution (Step 5 of Thm ns_regularity).
     This is the open question.

STRUCTURAL INSIGHT:
  6. The DS rule is commutative => no Levi-Civita structure in Q tensor.
  7. The antisymmetric cross-product is the PRECISE obstruction to
     extending DS regularity to NS.
  8. The Born floor CONSTRAINT exists on the bundle (topological).
     Whether the NS DYNAMICS preserves it is analytical.

THE GAP:
  The Born floor is a topological/algebraic constraint on the bundle B_Z.
  The NS dynamics is a PDE on R^3 (the base). The minitwistor Penrose
  transform connects them (Hitchin correspondence). But the Penrose
  transform is a t=0 isomorphism; whether it persists under NS evolution
  is the open question.
""")

# ═══════════════════════════════════════════════════════════════════
# PART 9: Quantitative Summary
# ═══════════════════════════════════════════════════════════════════

print("=" * 72)
print("PART 9: QUANTITATIVE SUMMARY")
print("=" * 72)

print(f"""
Key quantities:
  Born floor:              1/27 = {1/27:.10f}
  Max vorticity:           sqrt(26) = {np.sqrt(26):.10f}
  Curvature bound C(H):   {C_H}
  Spectral gap Delta:      {Delta:.4f}
  lambda_0:                {lambda_0}

Rate comparison (naive):
  DS enforcement rate:     {enforcement_rate_continuous:.4f}
  Stretching rate (C(H)):  {C_H**2:.6f}  (wrong quantity -- gauge sector)
  Stretching rate (omega): {26:.1f}       (correct -- vorticity sector)
  Ratio (C(H)):            {ratio_CH:.1f}x in favor of enforcement
  Ratio (omega):           {1/ratio_omega:.1f}x in favor of stretching

Honest assessment:
  - The C(H) = 0.344 rate comparison is MISLEADING because C(H) bounds
    the gauge-sector Jacobian, not the vorticity stretching rate.
  - Using the correct vorticity bound sqrt(26), the stretching rate (26)
    exceeds the enforcement rate (1.26) by a factor of 20.
  - This does NOT mean the Born floor is violated, because:
    (a) stretching and enforcement act on different components,
    (b) diffusion opposes stretching at maxima,
    (c) the Biot-Savart constraint couples u to omega.
  - But it does mean a simple rate comparison is insufficient.
  - The question is equivalent to the NS Millennium Problem.

Numerical findings (1D lattice simulation):
""")

for name, res in results.items():
    min_ever = res["min_borns"].min()
    max_omega_ever = res["max_omegas"].max()
    print(f"  {name}:")
    print(f"    Min Born: {min_ever:.6f}, Max omega: {max_omega_ever:.4f}")

print(f"""
The simulation shows:
  - With Born floor enforcement, all scenarios preserve Born >= 1/27
    (by construction -- the enforcement is applied every step).
  - Without enforcement, stretching CAN push Born below 1/27
    (Part 4 results), demonstrating the stretching is real.
  - The question is whether the GEOMETRIC (not computational)
    Born floor acts continuously, preventing any transient violation.

STATUS: GENUINELY OPEN
  The paper is honest about this (line 1905):
  "The honest status: if Born floor preservation holds, regularity follows.
   Whether it holds is equivalent to the Navier-Stokes regularity problem."
""")

# ═══════════════════════════════════════════════════════════════════
# PART 10: What a Proof Would Need
# ═══════════════════════════════════════════════════════════════════

print("=" * 72)
print("PART 10: WHAT A PROOF WOULD REQUIRE")
print("=" * 72)

print(f"""
A proof of Born floor preservation under NS would need one of:

APPROACH A: Maximum principle on the Born functional.
  Show that d/dt Born(x,t) >= 0 whenever Born(x,t) = 1/27.
  Difficulty: the NS stretching term is not sign-definite at the boundary.
  The diffusion helps (opposes growth at maxima of |omega|^2), but
  the CZ bound ||nabla u||_infty from ||omega||_infty doesn't close.

APPROACH B: Energy method on the Born-allowed region B_Z.
  Show that the L^2 energy in the complement of B_Z is non-increasing.
  Difficulty: the NS nonlinearity can transfer energy across the boundary.

APPROACH C: Topological/geometric.
  Use the bundle structure of O(2) = T(CP^1) to constrain the NS dynamics.
  The Born floor is a topological constraint (c_1(O(2)) = 2 is an integer).
  The minitwistor representative lives in B_Z at t=0 by assumption.
  Show that the O(2) holomorphic structure constrains the possible
  deformations under NS evolution.
  Difficulty: the Penrose transform is an isomorphism at fixed time.
  Whether the holomorphic structure survives under NS flow is unclear.

APPROACH D: Probabilistic/measure-theoretic.
  Show that the set of initial data for which Born floor is violated
  has measure zero. This would give "generic" regularity but not
  deterministic regularity for all smooth data.

APPROACH E: Contradiction.
  Assume Born drops below 1/27 at some first time T* and point x*.
  Then |omega(x*, T*)| = sqrt(26) and d/dt |omega|^2 > 0 at (x*, T*).
  The stretching term must exceed diffusion at this point.
  Use the specific structure of the minitwistor-reconstructed flow to
  derive a contradiction (e.g., show the strain-vorticity alignment
  is bounded by the Penrose transform).

The most promising direction appears to be APPROACH E combined with C:
use the holomorphic structure of the Penrose transform to constrain
the strain-vorticity alignment at the Born boundary. This is the
"does the O(2) bundle structure constrain the stretching term?" question
that the paper identifies as the key open problem.
""")

print("=" * 72)
print("COMPUTATION COMPLETE")
print("=" * 72)
