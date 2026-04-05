"""
Penrose transform integral: Step 2.

Find the normalisation: ||∂̄Φ||² → Λ.

For UNIFORM data on CP³, the Penrose transform gives constant
curvature on S⁴. The cosmological constant Λ is:

  Λ = (normalisation) × ||∂̄Φ||²

The normalisation comes from the Penrose transform kernel
integrated over CP¹.

Key reference: the Penrose transform for the graviton maps
H¹(CP³, O(-6)) → massless spin-2 on S⁴.
For CONFORMAL gravity (our case via Mason), the relevant
cohomology is different: it's the deformation of J itself.
"""
import numpy as np

# ============================================================
# THE NORMALISATION FROM FIRST PRINCIPLES
# ============================================================
#
# Standard CP³ with Fubini-Study metric gives S⁴ with the
# round metric. The scalar curvature of S⁴ with radius r is:
#   R = n(n-1)/r² = 4·3/r² = 12/r²
# so Λ = R/4 = 3/r².
#
# The Fubini-Study metric on CP³ has:
#   Scalar curvature R_{FS} = n(n+1)/r² = 3·4/r² = 12/r²
# Wait, for CPⁿ with Fubini-Study of holomorphic sectional
# curvature 4/r²:
#   R_{FS} = 2n(n+1)/r²
# For CP³ (n=3): R_{FS} = 24/r²
#
# The relationship between CP³ and S⁴ curvatures:
# R_{S⁴} = 12/r², R_{CP³} = 24/r²
# So R_{S⁴} = R_{CP³}/2.
# And Λ = R_{S⁴}/4 = R_{CP³}/8 = 3/r².
#
# For unit Fubini-Study radius (r=1):
# R_{CP³} = 24, R_{S⁴} = 12, Λ = 3.

print("=" * 60)
print("STANDARD CURVATURES (unit radius)")
print("=" * 60)
print(f"CP³ Fubini-Study: R = 24")
print(f"S⁴ round:         R = 12")
print(f"Λ = 3")
print()

# ============================================================
# THE DS DEFORMATION
# ============================================================
#
# The standard complex structure J₀ on CP³ gives the round S⁴.
# The DS deformation J' = J₀ + δJ changes the curvature.
#
# For a UNIFORM deformation (constant across CP³):
# The deformed CP³ still gives constant curvature on S⁴.
# The Weyl tensor remains zero (conformally flat).
# Only Λ changes.
#
# The change in Λ comes from the (1,1) component of ∂̄Φ
# (the trace/conformal part). The (3,3) component (graviton)
# gives W ≠ 0 for NON-UNIFORM deformations. But for uniform
# deformations, the (3,3) part integrates to zero over CP¹
# (by symmetry — constant symmetric traceless tensor integrated
# over CP¹ with its standard measure vanishes).
#
# Wait. Let me reconsider.
#
# If the deformation is uniform (constant ∂̄Φ across CP³),
# then restricted to any twistor line L_x, it's also constant.
# The Penrose transform of a constant function on CP¹ is:
#
# For spin-0: ∮ f dζ∧dζ̄ = f · Vol(CP¹) = f · π (area of S²)
# For spin-2: the kernel has specific ζ-dependence.
#
# The linearised Penrose transform for spin-2:
#   h_μν(x) = ∮_{L_x} f(Z) · π^{(A'} π^{B')} dζ
# where π = (ζ, 1) is the spinor on CP¹.
#
# For CONSTANT f:
#   h_{A'B'}(x) = f · ∮ π_{A'} π_{B'} dζ
#
# The integral ∮ π_{A'} π_{B'} dζ over CP¹:
# In components with π = (ζ, 1):
#   ∮ ζ² dζ = 0 (contour integral of holomorphic function)
#   ∮ ζ dζ = 0
#   ∮ 1 dζ = 0
#
# ALL ZERO by Cauchy's theorem (no poles)!
#
# This confirms: constant spin-2 data integrates to ZERO on S⁴.
# The graviton field h_μν = 0 at equilibrium.
# No gravitational waves, no Weyl curvature. Just Λ.
#
# The cosmological constant comes from the SCALAR (spin-0) part,
# which is the (1,1) trace component of ∂̄Φ.
# For spin-0: ∮ f dζ = 2πi × Res(f) = 0 for constant f too...
#
# Hmm. Actually for the Penrose transform of spin-0 (scalar field):
# φ(x) = (1/2πi) ∮ f(Z(ζ)) / (π·ô)^{2s+2} dζ
# For s=0: φ = ∮ f / (something)² dζ
# This involves the twistor function divided by spinor contractions.
#
# For CONSTANT f, the integral gives a spacetime-dependent result
# from the ζ-dependence of the denominator.
#
# Let me think about this differently.

print("=" * 60)
print("RECONSIDERING: HOW Λ EMERGES")
print("=" * 60)
print()

# The Penrose transform maps twistor DATA to spacetime FIELDS.
# The "data" is not just ∂̄Φ but the full structure of the
# deformed almost complex structure.
#
# For the NON-LINEAR Penrose transform (nonlinear graviton):
# A deformed CP³ (with deformed J) maps to a curved S⁴.
# The AMOUNT of deformation determines the curvature.
#
# For a UNIFORM deformation:
# The deformed CP³ is still a homogeneous space.
# Homogeneous deformations of CP³ → homogeneous metrics on S⁴.
# Homogeneous metrics on S⁴ = round metric with some radius.
# The radius (hence Λ) is determined by the deformation strength.
#
# So Λ doesn't come from the Penrose transform INTEGRAL per se.
# It comes from the GLOBAL GEOMETRY of the deformed CP³.
#
# The deformed CP³ has a modified Fubini-Study metric.
# The modification is determined by the Born floor enforcement.
# The scalar curvature of the modified CP³ determines Λ.
#
# This is much simpler than doing a contour integral!

print("The cosmological constant comes from the GLOBAL GEOMETRY")
print("of the deformed CP³, not from a contour integral.")
print()
print("At equilibrium:")
print("  Standard CP³: R_FS = 24 → Λ = 3 (unit radius)")
print("  Deformed CP³: R_FS' = 24 + δR → Λ' = 3 + δΛ")
print()
print("The deformation δR comes from the Born floor's effect on")
print("the Fubini-Study metric.")
print()

# ============================================================
# THE BORN FLOOR'S EFFECT ON THE FUBINI-STUDY METRIC
# ============================================================
#
# The Fubini-Study metric in affine coordinates w_i = Z_i/Z_0:
#   g_{ij̄} = ((1+|w|²)δ_{ij} - w̄_i w_j) / (1+|w|²)²
#
# The Born floor constrains which points of CP³ are "allowed"
# (those with Born(θ) ≥ 1/27). The allowed region is a subset
# of CP³ — the complement of the "forbidden zone" where
# Born(θ) < 1/27.
#
# In the DS dynamics, the mass function is FORCED to stay in
# the allowed region. The floor enforcement acts as a POTENTIAL
# WALL at the boundary Born = 1/27.
#
# The effective geometry "as seen by DS dynamics" is CP³ with
# the Fubini-Study metric PLUS a potential that diverges at
# the floor boundary. This is like a particle in a box:
# the box geometry determines the effective curvature.
#
# For the equilibrium at Born = 1/27:
# The mass function sits at the floor boundary.
# The "box" has a specific shape determined by the floor.
# The effective curvature of this box determines Λ.
#
# The allowed region {Born(θ) ≥ 1/27} in CP³ coordinates:
# Born = |Z⁰|²/|Z|² ≥ 1/27
# In affine coords (w_i = Z_i/Z⁰):
# 1/(1+|w|²) ≥ 1/27
# |w|² ≤ 26
# This is a BALL of radius √26 in the affine chart!

print("=" * 60)
print("THE ALLOWED REGION IN CP³")
print("=" * 60)
print()
print("Born(θ) = |Z⁰|²/|Z|² ≥ 1/27")
print("In affine coordinates w_i = Z_i/Z⁰:")
print("  1/(1+|w|²) ≥ 1/27")
print("  |w|² ≤ 26")
print()
print("The allowed region is a ball B(0, √26) ⊂ C³ ≅ R⁶.")
print(f"Radius: √26 = {np.sqrt(26):.4f}")
print()
print("The Fubini-Study metric on this ball:")
print("At the boundary |w|² = 26:")
print(f"  1+|w|² = 27")
print(f"  g_ij̄ ~ δ_ij/27² (highly contracted)")
print()
print("At the centre |w| = 0:")
print(f"  g_ij̄ = δ_ij (flat)")
print()

# The equilibrium mass function in affine coords:
# m* = (s₁, s₂, s₃, θ) → w_i = s_i/θ
# At K*=7/30: m* = (0.7869, 0.0293, 0.0293, 0.1545)
# w = (0.7869/0.1545, 0.0293/0.1545, 0.0293/0.1545)
#   = (5.093, 0.190, 0.190)
# |w|² = 25.94 + 0.036 + 0.036 = 26.01 ≈ 26

m_star = np.array([0.7869, 0.0293, 0.0293, 0.1545])
w_star = m_star[:3] / m_star[3]
w_sq = np.sum(w_star**2)

print(f"Equilibrium in affine coords:")
print(f"  w* = ({w_star[0]:.4f}, {w_star[1]:.4f}, {w_star[2]:.4f})")
print(f"  |w*|² = {w_sq:.4f}")
print(f"  26 = {26}")
print(f"  |w*|² ≈ 26: {abs(w_sq - 26) < 0.1}")
print()
print("THE EQUILIBRIUM SITS ON THE BOUNDARY OF THE ALLOWED REGION.")
print("Born = 1/27 means |w|² = 26 exactly.")
print()
print("This is BEAUTIFUL: the DS fixed point is on the BOUNDARY")
print("of the ball B(0,√26) ⊂ CP³. The floor holds it there.")
print("The effective geometry is a ball, not all of CP³.")
print()

# ============================================================
# CURVATURE OF THE BALL
# ============================================================
print("=" * 60)
print("CURVATURE OF THE EFFECTIVE GEOMETRY")
print("=" * 60)
print()

# The Fubini-Study metric on the ball B(0,√26):
# This is the Bergman metric of the ball, which is the
# restriction of the Fubini-Study metric.
#
# The scalar curvature of CP³ restricted to B(0,R) with R²=26:
# The Fubini-Study scalar curvature is constant = 24.
# It doesn't change on subsets (it's intrinsic, not extrinsic).
#
# BUT: the dynamics are CONFINED to the ball by the floor.
# The floor enforcement changes the EFFECTIVE metric:
# the dynamics see a metric that includes the floor potential.
#
# This is like quantum mechanics in a box:
# the box doesn't change the metric, but it changes the
# effective Hamiltonian (adds boundary conditions).
#
# For the Penrose transform: what matters is the GLOBAL
# geometry of the deformed CP³. The deformation is:
# J' = J₀ on B(0,√26)
# J' = J₀ + δJ on the boundary ∂B (where floor activates)
#
# The boundary condition is: the mass function is reflected
# at Born = 1/27. This is a NEUMANN-type boundary condition
# on the ball.
#
# Maldacena's result: conformal gravity + Neumann BC = Einstein + Λ.
# Our BC: reflection at Born = 1/27 on B(0,√26) ⊂ CP³.
# The value of Λ is determined by the radius of the ball.

# For a ball of radius R in CP³:
# The "effective Λ" scales as 1/R².
# R² = 26 → Λ_eff ∝ 1/26.

# Standard CP³ (R → ∞, full space): Λ = 3.
# Ball of radius √26: Λ_ball = 3/f(26) for some function f.

# Actually, the simplest relationship:
# Λ_S⁴ = 3/r² where r is the S⁴ radius.
# The S⁴ radius is related to the CP³ "size" by the
# Hopf fibration: Vol(CP³) = π³r⁶/6, Vol(S⁴) = 8π²r⁴/3.
# The ratio Vol(CP³)/Vol(S⁴) = π r²/16.
#
# For the ball B(0,√26), the effective "size" is R = √26.
# In units where the full CP³ has unit radius:
# the ball has radius √26 in the AFFINE chart,
# which corresponds to the entire CP³ minus a small cap.
# (|w|² = 26 means we're using almost all of CP³.)

# Let's compute what fraction of CP³ is excluded:
# The excluded region is {|w|² > 26} in the affine chart.
# In the Fubini-Study metric, the volume element is:
# dV = (1+|w|²)^{-(n+1)} dV_Eucl  for CPⁿ (n=3)
# Vol(CP³) = ∫₀^∞ r⁵/(1+r²)⁴ dr × Vol(S⁵)

# Fraction inside ball:
from scipy import integrate

def integrand_full(r):
    """Volume element of CP³ in radial coordinate."""
    return r**5 / (1 + r**2)**4

vol_total, _ = integrate.quad(integrand_full, 0, np.inf)
vol_ball, _ = integrate.quad(integrand_full, 0, np.sqrt(26))

frac = vol_ball / vol_total
print(f"Volume fraction of CP³ inside B(0,√26):")
print(f"  Vol(ball)/Vol(CP³) = {frac:.6f}")
print(f"  Excluded fraction = {1-frac:.6f}")
print()

# The effective cosmological constant:
# If the excluded fraction is small, Λ barely changes from 3.
# If significant, Λ changes proportionally.

print(f"The Born floor excludes {(1-frac)*100:.2f}% of CP³.")
print(f"This is a {(1-frac)*100:.2f}% correction to the geometry.")
print()

# What is 1/27 in these terms?
# Born = 1/27 corresponds to |w|² = 26.
# Born = 1 (pure θ) corresponds to w = 0 (centre).
# Born = 0 (no θ) corresponds to |w|² → ∞ (boundary of CP³).
# The floor cuts off the high-|w| region.

print("=" * 60)
print("RESULT")
print("=" * 60)
print()
print("At equilibrium, the DS framework produces:")
print("  • Einstein gravity (Theorem thm:einstein)")
print("  • On S⁴ (from CP³ via twistor correspondence)")
print("  • With Λ > 0 (from the Fubini-Study curvature)")
print("  • The floor modifies Λ by excluding a fraction of CP³")
print(f"  • Excluded fraction: {(1-frac)*100:.2f}%")
print()
print("The cosmological constant at the DS equilibrium is:")
print(f"  Λ_DS = 3 × (correction from {(1-frac)*100:.2f}% exclusion)")
print()
print("For the standard Penrose correspondence (unit FS radius):")
print(f"  Λ_standard = 3")
print(f"  Λ_DS ≈ 3 × {frac:.4f} = {3*frac:.4f}")
print(f"  (This is an estimate; the exact relationship requires")
print(f"   the Penrose transform kernel normalisation.)")
