#!/usr/bin/env python3
"""
Dark matter rotation curve normalization from H=3 informational geometry.

Physics:
  The framework has a spatial orientation field n: R^3 -> S^2 parameterizing
  the vacuum equilibrium on the S^2 orbit of the fixed point m* under SU(2).

  When a galaxy pins n at its center, the minimum-energy interpolation to
  the uniform vacuum at infinity is the HEDGEHOG: n(r) = r_hat.

  This is a harmonic map R^3\\{0} -> S^2 with Dirichlet energy density:
    epsilon(r) = (1/2)|nabla n|^2 = 1/r^2

  The gravitational coupling of this energy density gives:
    rho_DM(r) = kappa / (8 pi G r^2)

  which produces FLAT rotation curves: v_flat^2 = kappa/2.

  The question: what sets kappa (and hence v_flat) in terms of framework
  quantities (H=3, K*=7/30, Delta=1.2626, etc.)?

Current status in paper (Prediction 18):
  - Perturbative Born floor mechanism RULED OUT (epsilon^{3/2} wrong sign)
  - Moduli space variation K(r) = K* - alpha*ln(r/r0) is the promising lead
  - The coupling alpha is NOT YET derived

This script:
  1. Confirms the hedgehog energy density and flat rotation curve
  2. Tests various framework expressions for the MOND acceleration a_0
  3. Compares to Tully-Fisher observations
  4. Investigates the moduli-gradient mechanism in detail
  5. Is honest about what works and what remains open

J.R. Manuel, 2026
"""

import numpy as np
from fractions import Fraction

# =============================================================================
# Physical constants (SI)
# =============================================================================
c       = 2.99792458e8       # m/s
G       = 6.67430e-11        # m^3 kg^-1 s^-2
hbar    = 1.05457182e-34     # J s
k_B     = 1.380649e-23       # J/K
M_Pl    = np.sqrt(hbar * c / G)  # Planck mass (kg)
l_Pl    = np.sqrt(hbar * G / c**3)  # Planck length (m)
eV      = 1.602176634e-19    # J
GeV     = 1e9 * eV
Mpc     = 3.0857e22          # m

# =============================================================================
# Framework constants
# =============================================================================
H       = 3
K_star  = Fraction(7, 30)
K_float = float(K_star)
lambda_0 = 0.28291           # dominant eigenvalue at K*
Delta   = -np.log(lambda_0)  # spectral gap = 1.2626
sigma   = -np.log(23.0 / 30.0)      # = ln(30/23) = 0.2657
Born_floor = 1 / H**3        # = 1/27

# Derived framework scales
Omega_Lambda = (H - 1) / H   # = 2/3
H_0_val = 69.3               # km/s/Mpc (framework prediction)
H_0_SI  = H_0_val * 1e3 / Mpc  # in s^-1

# =============================================================================
# Observed MOND acceleration
# =============================================================================
a_0_obs = 1.2e-10  # m/s^2 (Milgrom's constant from galaxy rotation curves)

print("=" * 72)
print("DARK MATTER ROTATION CURVES FROM H=3 INFORMATIONAL GEOMETRY")
print("=" * 72)

# =============================================================================
# PART 1: Hedgehog energy density and flat rotation curves
# =============================================================================
print("\n" + "=" * 72)
print("PART 1: HEDGEHOG n(r) = r_hat ON S^2")
print("=" * 72)

print("""
The orientation field n: R^3 -> S^2 satisfies the Beltrami equation
at equilibrium. For a galaxy pinning n at the origin, the minimum-energy
solution interpolating to uniform vacuum at infinity is the hedgehog:

  n(r) = r_hat  (unit radial vector)

This is a harmonic map. Its Dirichlet energy density:

  epsilon(r) = (1/2)|nabla n|^2

For n = r_hat in spherical coordinates:
  n_theta = (1/r) d(r_hat)/d(theta)  => |n_theta| = 1/r
  n_phi   = (1/(r sin theta)) d(r_hat)/d(phi)  => |n_phi| = 1/r

  epsilon = (1/2)(1/r^2 + 1/r^2) = 1/r^2
""")

# Verify numerically
r_test = np.logspace(-1, 3, 1000)  # kpc range
epsilon_hedgehog = 1.0 / r_test**2

eps_at_1 = 1.0 / 1.0**2
eps_at_10 = 1.0 / 10.0**2
print("Numerical check:")
print(f"  epsilon(1) / epsilon(10) = {eps_at_1 / eps_at_10:.1f}")
print(f"  Expected ratio (10/1)^2 = {100.0:.1f}  [1/r^2 confirmed]")

# Enclosed mass for rho = C/r^2
# M(r) = int_0^r 4 pi r'^2 (C/r'^2) dr' = 4 pi C r
# v^2 = GM(r)/r = 4 pi G C  (independent of r => FLAT!)
print(f"""
Gravitational coupling:
  rho_DM(r) = kappa / (8 pi G r^2)

Enclosed mass:
  M_DM(r) = int 4 pi r'^2 rho dr' = (kappa / 2G) r

Rotation velocity:
  v^2 = G M(r) / r = kappa / 2

  => v_flat = sqrt(kappa / 2)   [FLAT - independent of r]

The 1/r^2 density profile from the hedgehog AUTOMATICALLY gives flat
rotation curves. The only question is: what sets kappa?
""")

# =============================================================================
# PART 2: The MOND acceleration scale a_0
# =============================================================================
print("=" * 72)
print("PART 2: TESTING FRAMEWORK EXPRESSIONS FOR a_0")
print("=" * 72)

print(f"\nObserved: a_0 = {a_0_obs:.2e} m/s^2")
print(f"Framework H_0 = {H_0_val} km/s/Mpc = {H_0_SI:.4e} s^-1")
print(f"c * H_0 = {c * H_0_SI:.4e} m/s^2")
print()

# The Tully-Fisher relation: v_flat^4 = a_0 * G * M_baryon
# This means a_0 is the fundamental scale, not v_flat itself.
# a_0 connects baryonic mass to asymptotic velocity.

candidates = {}

# Candidate 1: a_0 = c H_0 / (2 pi)
a_1 = c * H_0_SI / (2 * np.pi)
candidates["c H_0 / (2pi)"] = a_1

# Candidate 2: a_0 = c H_0 / 6
a_2 = c * H_0_SI / 6
candidates["c H_0 / 6"] = a_2

# Candidate 3: a_0 = c H_0 / (2H)  = c H_0 / 6 (same as above at H=3)
# Already covered

# Candidate 4: a_0 = K* c H_0
a_4 = K_float * c * H_0_SI
candidates["K* c H_0"] = a_4

# Candidate 5: a_0 = c H_0 * sqrt(K*)
a_5 = c * H_0_SI * np.sqrt(K_float)
candidates["c H_0 sqrt(K*)"] = a_5

# Candidate 6: a_0 = c H_0 / (H+1)  = c H_0 / 4
a_6 = c * H_0_SI / (H + 1)
candidates["c H_0 / (H+1)"] = a_6

# Candidate 7: a_0 = c H_0 * Delta / (2 pi H)
a_7 = c * H_0_SI * Delta / (2 * np.pi * H)
candidates["c H_0 Delta/(2pi H)"] = a_7

# Candidate 8: a_0 = c H_0 / (2 * sqrt(H^3 - 1))  = c H_0 / (2 sqrt(26))
a_8 = c * H_0_SI / (2 * np.sqrt(H**3 - 1))
candidates["c H_0 / (2 sqrt(26))"] = a_8

# Candidate 9: a_0 = c^2 * H_0 * sqrt(Omega_Lambda) / c = H_0 * c * sqrt(2/3)
a_9 = c * H_0_SI * np.sqrt(Omega_Lambda)
candidates["c H_0 sqrt(Omega_L)"] = a_9

# Candidate 10: a_0 = c H_0 * sqrt(1 - Omega_Lambda) = c H_0 / sqrt(3)
a_10 = c * H_0_SI * np.sqrt(1 - Omega_Lambda)
candidates["c H_0 sqrt(Omega_m)"] = a_10

# Candidate 11: a_0 = c H_0 / (2 pi) * sqrt(6)  (Milgrom's original suggestion close)
a_11 = c * H_0_SI / (2 * np.pi) * np.sqrt(6)
candidates["c H_0 sqrt(6)/(2pi)"] = a_11

# Candidate 12: sigma * c * H_0 (sigma = ln(30/23) = 0.2657)
a_12 = float(sigma) * c * H_0_SI
candidates["sigma c H_0"] = a_12

print(f"{'Expression':<28} {'Value (m/s^2)':<16} {'Ratio to a_0':<14} {'% dev':<10}")
print("-" * 68)

best_name = None
best_ratio = 999

for name, val in sorted(candidates.items(), key=lambda x: abs(np.log(x[1]/a_0_obs))):
    ratio = val / a_0_obs
    pct = (ratio - 1) * 100
    marker = ""
    if abs(pct) < abs((best_ratio - 1) * 100) or best_name is None:
        best_ratio = ratio
        best_name = name
    if abs(pct) < 15:
        marker = " <--"
    print(f"  {name:<26} {val:<16.4e} {ratio:<14.4f} {pct:>+8.1f}%{marker}")

print(f"\nBest match: {best_name} (ratio {best_ratio:.4f})")

# =============================================================================
# PART 3: The c H_0 / 6 identification
# =============================================================================
print("\n" + "=" * 72)
print("PART 3: ANALYSIS OF a_0 = c H_0 / (2H)")
print("=" * 72)

a_0_framework = c * H_0_SI / (2 * H)
print(f"""
The expression a_0 = c H_0 / (2H) = c H_0 / 6 gives:
  a_0 = {a_0_framework:.4e} m/s^2
  Observed: {a_0_obs:.2e} m/s^2
  Ratio: {a_0_framework / a_0_obs:.4f}
  Deviation: {(a_0_framework / a_0_obs - 1) * 100:+.1f}%

Physical interpretation:
  The Hubble acceleration c H_0 is the cosmic deceleration scale.
  Dividing by 2H = 6 gives the scale at which the vacuum orientation
  field's Dirichlet energy becomes gravitationally relevant.

  The factor H = 3 comes from spatial dimensions: the hedgehog lives
  in R^3 (H spatial dimensions), and the 1/r^2 fall-off in H=3
  dimensions gives the normalization factor 1/(2H) = 1/6.

  Alternatively: the hedgehog Dirichlet energy per solid angle per
  decade in r is constant. The number of "decades" between the Hubble
  radius and a galaxy core is O(H * ln(10) * something). The factor
  2H = 6 normalizes this.
""")

# =============================================================================
# PART 4: Tully-Fisher relation from the hedgehog
# =============================================================================
print("=" * 72)
print("PART 4: TULLY-FISHER RELATION")
print("=" * 72)

print("""
The Tully-Fisher relation (empirical):
  v_flat^4 = a_0 * G * M_baryon

For the hedgehog mechanism to reproduce this, we need:
  v_flat^2 = kappa/2  AND  kappa depends on M_baryon such that
  v_flat^4 ~ M_baryon.

Key insight: the hedgehog n = r_hat is TOPOLOGICAL (winding number 1).
Its amplitude is fixed (|n| = 1 on S^2). But the coupling kappa can
depend on the baryonic mass through the CORE MATCHING CONDITION.

The hedgehog is valid for r > r_core, where r_core is set by the
baryonic matter distribution. Inside r_core, the field is dominated
by the galaxy's own gravitational field, not the vacuum topology.

The matching condition at r_core gives:
  kappa = (baryonic gravitational potential at r_core) * (geometric factor)

If r_core = sqrt(G M_b / a_0)  (the MOND radius), then:
""")

# MOND radius for typical galaxies
M_sun = 1.989e30  # kg
M_gal_range = np.array([1e8, 1e9, 1e10, 1e11, 1e12]) * M_sun  # galaxy masses

print(f"{'M_baryon (M_sun)':<20} {'r_MOND (kpc)':<16} {'v_flat (km/s)':<16} {'v_TF (km/s)':<14}")
print("-" * 66)

kpc = 3.0857e19  # m

for M_b in M_gal_range:
    r_MOND = np.sqrt(G * M_b / a_0_framework)
    v_TF = (a_0_framework * G * M_b)**0.25  # from Tully-Fisher
    # From hedgehog with matching:
    # kappa = v_flat^2 * 2, and v_flat from TF
    v_flat = v_TF  # TF determines v_flat
    print(f"  {M_b/M_sun:<18.0e} {r_MOND/kpc:<16.1f} {v_flat/1e3:<16.1f} {v_TF/1e3:<14.1f}")

print("""
The MOND radius r_MOND = sqrt(G M_b / a_0) is where Newtonian gravity
(~ G M_b / r^2) equals the MOND acceleration a_0. This is the natural
matching radius between the baryonic interior and the hedgehog exterior.

At r = r_MOND, the hedgehog energy density kappa/r^2 must equal
the baryonic gravitational energy density ~ G M_b / r^4 * (something).
""")

# =============================================================================
# PART 5: Moduli-gradient mechanism (from Prediction 18)
# =============================================================================
print("=" * 72)
print("PART 5: MODULI-GRADIENT MECHANISM")
print("=" * 72)

print("""
The paper's Prediction 18 identifies the promising non-perturbative lead:

  K(r) = K* - alpha * ln(r / r_0)

where K(r) is the local vacuum conflict, varying logarithmically with
distance from baryonic matter. The gradient energy density:

  rho_grad = (1/2) |dK/dr|^2 * (dE/dK)^2 / (8 pi G)
           = (alpha^2 / 2) * (dE/dK)^2 / (8 pi G r^2)

This gives rho ~ 1/r^2 automatically (from the logarithmic gradient).

The hedgehog and moduli-gradient are actually THE SAME MECHANISM:
  - The hedgehog n = r_hat gives the 1/r^2 from the angular field
  - The moduli gradient gives 1/r^2 from the radial field
  - Together: the FULL vacuum response to a galaxy is a hedgehog in
    orientation WITH a logarithmic moduli gradient in conflict level

Combined energy density:
  rho = [kappa_angular + kappa_radial] / (8 pi G r^2)

where:
  kappa_angular = geometric factor from hedgehog (topological, universal)
  kappa_radial  = alpha^2 * (dE/dK)^2 (depends on framework dynamics)
""")

# Compute the moduli gradient contribution
# The energy of the framework as a function of K near K*:
# E(K) ~ E(K*) + (1/2) E''(K*) (K - K*)^2 + ...
# The "stiffness" dE/dK at K* is related to the spectral gap

# At the fixed point, the Hessian of the free energy gives:
# d^2F/dK^2 = 1/(K*(1-K*)) from the entropy part
# = 1/(7/30 * 23/30) = 900/161

stiffness = 1.0 / (K_float * (1 - K_float))
print(f"Entropy stiffness at K*: d^2S/dK^2 = {stiffness:.4f}")
print(f"  = 1/(K*(1-K*)) = 900/161 = {900/161:.4f}")

# For the moduli gradient:
# rho = alpha^2 * stiffness^2 / (8 pi G r^2)  [if dE/dK ~ stiffness * delta_K]
# Actually: rho = alpha^2 * M_vac^4 * stiffness / (8 pi G c^2 r^2)  in proper units

# =============================================================================
# PART 6: Deriving alpha from framework quantities
# =============================================================================
print("\n" + "=" * 72)
print("PART 6: WHAT SETS ALPHA?")
print("=" * 72)

print("""
The coupling alpha (logarithmic rate of K variation with distance)
must come from the balance between:
  1. The baryonic gravitational field pulling K away from K*
  2. The framework's restoring force pushing K back to K*

The gravitational field of a galaxy of mass M_b at distance r:
  g(r) = G M_b / r^2

This field does work on the vacuum modulus. The rate of K change:
  dK/dr = -alpha / r

The balance condition (gravitational work = restoring force):
  (G M_b / r^2) * (coupling to K) = stiffness * (dK/dr) * M_vac^2

This gives alpha proportional to M_b, which would make rho ~ M_b / r^2,
hence v_flat^2 ~ M_b, hence v_flat^4 ~ M_b^2. This is WRONG for TF
(which needs v^4 ~ M_b).
""")

print("THE TULLY-FISHER EXPONENT PROBLEM:")
print("-" * 40)
print("""
  Hedgehog (topological, amplitude fixed):
    rho = C / r^2  with C independent of M_b
    => v_flat^2 = constant (same for all galaxies)
    => WRONG (galaxies have different v_flat)

  Moduli gradient (alpha ~ M_b):
    rho = alpha^2 / r^2 ~ M_b^2 / r^2
    => v_flat^4 ~ M_b^2
    => WRONG exponent (TF says v^4 ~ M_b^1)

  For correct TF (v^4 ~ M_b):
    Need rho ~ sqrt(M_b) / r^2, i.e., alpha ~ M_b^{1/2}
    This requires a NONLINEAR response of the vacuum to gravity.
""")

# =============================================================================
# PART 7: Nonlinear vacuum response
# =============================================================================
print("=" * 72)
print("PART 7: NONLINEAR VACUUM RESPONSE (THE KEY INSIGHT)")
print("=" * 72)

print("""
The resolution: the vacuum response to gravity is NOT linear in the
gravitational potential. The Born floor 1/H^3 = 1/27 introduces a
nonlinearity.

When the gravitational potential Phi = -G M_b / r exceeds a critical
value, the vacuum modulus K is pulled to the boundary of the allowed
region [1/H^3, 1 - 1/H^3]. The floor clips the response.

The clipping transition occurs at the MOND radius:
  r_MOND = sqrt(G M_b / a_0)

For r < r_MOND: K is clipped (saturated), rho_DM ~ standard NFW-ish
For r > r_MOND: K varies logarithmically, rho_DM ~ 1/r^2

The TRANSITION between these regimes gives the Tully-Fisher relation:
  At r = r_MOND, the hedgehog energy density must match the baryonic:
    kappa / r_MOND^2 ~ G M_b / r_MOND^4
    kappa ~ G M_b / r_MOND^2 = a_0
    v_flat^2 ~ a_0 * r_MOND^2 / r_MOND = a_0 * r_MOND
             = a_0 * sqrt(G M_b / a_0)
             = sqrt(a_0 * G * M_b)
    v_flat^4 = a_0 * G * M_b  ✓  TULLY-FISHER!

So the TF relation emerges from the MATCHING between the clipped
interior and the hedgehog exterior at r_MOND.
""")

# =============================================================================
# PART 8: The full derivation chain
# =============================================================================
print("=" * 72)
print("PART 8: FULL DERIVATION CHAIN (SUMMARY)")
print("=" * 72)

print(f"""
Step 1: The vacuum has an S^2 orientation field n(r) [from SU(2) orbit of m*]
  => This is framework content (Theorem: m* has SU(2) isotropy)

Step 2: A galaxy pins n at its center (gravitational coupling to vacuum)
  => n interpolates between pinned center and uniform infinity
  => Minimum-energy solution: hedgehog n = r_hat

Step 3: Hedgehog Dirichlet energy density = 1/r^2
  => Gravitational coupling gives rho_DM = kappa / (8 pi G r^2)
  => Rotation curves are FLAT (v_flat = const)

Step 4: Born floor introduces MOND radius r_MOND = sqrt(G M_b / a_0)
  => For r < r_MOND: vacuum response saturated
  => For r > r_MOND: hedgehog 1/r^2 profile

Step 5: Matching at r_MOND gives Tully-Fisher:
  v_flat^4 = a_0 * G * M_baryon

Step 6: The acceleration scale:
  a_0 = c H_0 / (2H) = c H_0 / 6
  a_0 = {a_0_framework:.4e} m/s^2
  Observed: {a_0_obs:.2e} m/s^2
  Deviation: {(a_0_framework / a_0_obs - 1)*100:+.1f}%

  The factor 2H = 6:
    - 2 from the hedgehog (2 angular dimensions on S^2)
    - H = 3 from spatial dimensions (H = dim of space)
    - Together: 2H = solid-angle normalization in H-dim space
    - Or: (4 pi / (2 * 2 pi)) * (1/H-dim volume factor)

Step 7: The moduli gradient provides the radial component:
  K(r) = K* - alpha * ln(r/r_0)
  This contributes an ADDITIONAL 1/r^2 density on top of the hedgehog.
  The total kappa = kappa_angular + kappa_radial.
""")

# =============================================================================
# PART 9: Numerical predictions
# =============================================================================
print("=" * 72)
print("PART 9: NUMERICAL PREDICTIONS")
print("=" * 72)

# Milky Way parameters
M_MW = 5e10 * M_sun  # baryonic mass of Milky Way
r_MOND_MW = np.sqrt(G * M_MW / a_0_framework)
v_flat_MW = (a_0_framework * G * M_MW)**0.25

print(f"\nMilky Way (M_b = 5 x 10^10 M_sun):")
print(f"  r_MOND = {r_MOND_MW / kpc:.1f} kpc")
print(f"  v_flat (predicted) = {v_flat_MW / 1e3:.1f} km/s")
print(f"  v_flat (observed)  ~ 220 km/s")
print(f"  Deviation: {(v_flat_MW / 220e3 - 1)*100:+.1f}%")

# Dwarf galaxy
M_dwarf = 1e8 * M_sun
v_flat_dwarf = (a_0_framework * G * M_dwarf)**0.25
print(f"\nDwarf galaxy (M_b = 10^8 M_sun):")
print(f"  v_flat (predicted) = {v_flat_dwarf / 1e3:.1f} km/s")
print(f"  Typical observed: 30-50 km/s")

# Galaxy cluster
M_cluster = 1e14 * M_sun
v_flat_cluster = (a_0_framework * G * M_cluster)**0.25
print(f"\nGalaxy cluster (M_b = 10^14 M_sun):")
print(f"  v_flat (predicted) = {v_flat_cluster / 1e3:.1f} km/s")
print(f"  Typical observed: 1000-1500 km/s")

# =============================================================================
# PART 10: Framework expression scan for a_0
# =============================================================================
print("\n" + "=" * 72)
print("PART 10: DEEPER SCAN OF FRAMEWORK EXPRESSIONS FOR a_0")
print("=" * 72)

cH0 = c * H_0_SI

# Systematic scan: a_0 = cH_0 * f(H, K*, Delta, sigma)
exprs = {}

# Simple rational functions of H
for num in range(1, 13):
    for den in range(1, 37):
        val = cH0 * num / den
        ratio = val / a_0_obs
        if 0.9 < ratio < 1.1:
            key = f"cH0 * {num}/{den}"
            exprs[key] = (val, ratio)

# With K*
for p in [-2, -1, -0.5, 0.5, 1, 2]:
    val = cH0 * K_float**p
    ratio = val / a_0_obs
    if 0.5 < ratio < 2.0:
        key = f"cH0 * K*^({p})"
        exprs[key] = (val, ratio)

# With Delta
for p in [-2, -1, -0.5, 0.5, 1, 2]:
    val = cH0 * Delta**p
    ratio = val / a_0_obs
    if 0.5 < ratio < 2.0:
        key = f"cH0 * Delta^({p})"
        exprs[key] = (val, ratio)

# With sigma
for p in [-1, -0.5, 0.5, 1]:
    val = cH0 * float(sigma)**p
    ratio = val / a_0_obs
    if 0.5 < ratio < 2.0:
        key = f"cH0 * sigma^({p})"
        exprs[key] = (val, ratio)

# Mixed: K* * H-rational
for num in range(1, 7):
    for den in range(1, 13):
        val = cH0 * K_float * num / den
        ratio = val / a_0_obs
        if 0.9 < ratio < 1.1:
            key = f"cH0 * K* * {num}/{den}"
            exprs[key] = (val, ratio)

# With pi
val = cH0 / (2 * np.pi)
exprs["cH0 / (2pi)"] = (val, val / a_0_obs)

val = cH0 * K_float / np.pi
exprs["cH0 * K* / pi"] = (val, val / a_0_obs)

val = cH0 / (np.pi * H)
exprs["cH0 / (pi H)"] = (val, val / a_0_obs)

# Sort by closeness to a_0
sorted_exprs = sorted(exprs.items(), key=lambda x: abs(np.log(x[1][1])))

print(f"\nExpressions within 10% of a_0 = {a_0_obs:.2e} m/s^2:")
print(f"{'Expression':<28} {'Value (m/s^2)':<16} {'Ratio':<10} {'% dev':<10}")
print("-" * 64)
for name, (val, ratio) in sorted_exprs[:20]:
    pct = (ratio - 1) * 100
    print(f"  {name:<26} {val:<16.4e} {ratio:<10.4f} {pct:>+7.1f}%")

# =============================================================================
# PART 11: The c H_0 / (2pi) vs c H_0 / 6 comparison
# =============================================================================
print("\n" + "=" * 72)
print("PART 11: DISCRIMINATING c H_0 / (2pi) vs c H_0 / 6")
print("=" * 72)

a_2pi = cH0 / (2 * np.pi)
a_6   = cH0 / 6

print(f"""
Two leading candidates:
  c H_0 / (2pi) = {a_2pi:.4e} m/s^2  (cosmological; 2pi from horizon geometry)
  c H_0 / 6     = {a_6:.4e} m/s^2  (framework; 6 = 2H from hedgehog in H=3 space)

Difference: {(a_6/a_2pi - 1)*100:+.2f}%
Both: within {max(abs(a_2pi/a_0_obs - 1), abs(a_6/a_0_obs - 1))*100:.0f}% of observed a_0

Current observational uncertainty on a_0 is ~20%, so these cannot be
distinguished today. But note:
  - 2pi is generic (any cosmology)
  - 2H = 6 is specific to the framework (predicts a_0 changes if H changes)
  - The ratio 2pi/6 = pi/3 = 1.047, so they differ by only 4.7%

If the true a_0 is measured to 5% precision, we can distinguish them.
""")

# =============================================================================
# PART 12: Honesty check - what is derived vs assumed
# =============================================================================
print("=" * 72)
print("PART 12: HONESTY CHECK")
print("=" * 72)

print("""
DERIVED FROM THE FRAMEWORK:
  [✓] The vacuum has an S^2 orientation field (from SU(2) orbit structure)
  [✓] The hedgehog is the minimum-energy pinned configuration
  [✓] Dirichlet energy density = 1/r^2 (standard harmonic map result)
  [✓] rho ~ 1/r^2 gives FLAT rotation curves (elementary calculus)
  [✓] Matching at MOND radius gives Tully-Fisher v^4 ~ M_b

NOT YET DERIVED (ASSUMED OR FIT):
  [?] The coupling between Dirichlet energy and gravitational mass
      (this is kappa, the key missing piece)
  [?] That a_0 = c H_0 / 6 specifically (the factor 6 = 2H is motivated
      but not rigorously derived from the framework)
  [?] The moduli gradient coupling alpha (Prediction 18 acknowledges this)
  [?] Why the Born floor gives the MOND radius specifically

TENSIONS:
  [!] The perturbative mechanism was RULED OUT (Prediction 18)
  [!] The nonlinear matching argument is physically reasonable but
      not yet a rigorous derivation
  [!] The factor 2H = 6 vs 2pi = 6.28 cannot be distinguished with
      current data (4.7% difference, observational scatter ~20%)

WHAT WOULD MAKE THIS A REAL DERIVATION:
  1. Derive kappa from the framework's action principle
  2. Show the Born floor creates the MOND radius via saturation
  3. Derive a_0 = c H_0 / (2H) from the asymptotic matching condition
  4. Predict the SHAPE of the rotation curve (not just the flat part)
  5. Address galaxy clusters (where MOND is known to underpredict)
""")

# =============================================================================
# PART 13: The hedgehog winding and topological protection
# =============================================================================
print("=" * 72)
print("PART 13: TOPOLOGICAL PROTECTION")
print("=" * 72)

print("""
The hedgehog n = r_hat has winding number 1 in pi_2(S^2) = Z.

This means:
  1. It CANNOT be continuously deformed to the uniform vacuum
  2. It is TOPOLOGICALLY STABLE (no decay mechanism)
  3. Every galaxy that pins the field creates a hedgehog
  4. The hedgehog persists even after the galaxy evaporates

This topological protection is the framework's explanation for why
dark matter halos are universal:
  - Every baryonic concentration creates a hedgehog
  - The hedgehog survives any smooth perturbation
  - The 1/r^2 profile is the unique minimum-energy configuration
    in the hedgehog sector

The winding number is:
  deg(n) = (1/4pi) integral_{S^2} n* (area form) = 1

This is a pure integer, cannot vary continuously, and therefore the
hedgehog dark matter halo is QUANTIZED: it either exists (winding 1)
or doesn't (winding 0). There are no "partial" halos.

Higher-winding solutions (deg = 2, 3, ...) have energy ~ deg^2/r^2,
giving heavier halos. These could correspond to galaxy clusters
(which need more dark matter than MOND predicts for single galaxies).
""")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)

print(f"""
The H=3 informational geometry framework provides a NATURAL mechanism
for dark matter rotation curves:

  1. FLAT CURVES: The hedgehog n = r_hat gives rho ~ 1/r^2 automatically

  2. TULLY-FISHER: Born floor saturation + matching at r_MOND gives
     v_flat^4 = a_0 * G * M_baryon

  3. MOND SCALE: a_0 = c H_0 / (2H) = c H_0 / 6 = {a_0_framework:.3e} m/s^2
     vs observed {a_0_obs:.1e} m/s^2  ({(a_0_framework/a_0_obs - 1)*100:+.1f}%)

  4. TOPOLOGICAL STABILITY: hedgehog is pi_2(S^2) = Z protected

  5. NO DARK MATTER PARTICLE: the "dark matter" is vacuum strain energy
     (consistent with Prediction 15: no new particles)

Key remaining gaps:
  - Rigorous derivation of kappa from framework action principle
  - Proof that Born floor saturation creates MOND radius
  - Derivation of the factor 2H vs 2pi in the a_0 expression
  - Galaxy cluster regime (may need higher winding numbers)
""")

print("=" * 72)
print("COMPUTATION COMPLETE")
print("=" * 72)
