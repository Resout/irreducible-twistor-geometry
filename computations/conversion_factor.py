"""
Derive the DS-to-lattice conversion factor from geometry.

The conversion factor c = a_DS / a_lattice = 1.78 at β=2.3.

a_DS is the physical length of one DS combination step.
a_lattice is the lattice spacing.

The twistor fibration π: CP³ → S⁴ has fibre CP¹ ≅ S².
The DS step acts on the fibre. The Fubini-Study metric on CP³
induces a metric on the fibre. The DS step length in this metric
determines a_DS.

The Fubini-Study metric on CP^n:
  ds² = (|dz|²|z|² - |z̄·dz|²) / |z|⁴

On CP¹ (the fibre): ds² = |dw|² / (1+|w|²)²
where w is an affine coordinate. The diameter is π/2.

The DS step contracts the fibre state by factor λ₀ = 0.2829.
The Fubini-Study distance contracted per step is:
  Δd_FS = d · (1 - λ₀) ≈ 0.717 × d

For a characteristic perturbation at the scale of the equilibrium,
the FS distance is d ~ arctan(perturbation scale).

But the PHYSICAL length is not the FS distance — it's the FS distance
times a dimensional conversion factor that relates the CP³ metric to
the R⁴ metric through the twistor fibration.

The twistor fibration π: CP³ → S⁴ with fibre CP¹.
The metric on S⁴ (radius R) is:
  ds²_{S⁴} = R² dΩ²₄
The metric on CP³ (induced from C⁴):
  ds²_{CP³} = ds²_{FS}

The fibration preserves the metric in the horizontal (base) direction
and scales the vertical (fibre) direction by a factor related to the
twistor radius.

For the standard twistor fibration:
  Horizontal: ds²_base = R² dΩ²₄ (matches S⁴)
  Vertical: ds²_fibre = r² dΩ²₂ (S² fibre with radius r)

The ratio r/R determines the conversion between fibre steps and base steps.

For CP³ with Fubini-Study: the fibre CP¹ has FS radius 1 (in units
where the FS metric has sectional curvature in [1,4]).
The base S⁴ has radius R.

The physical mass gap is m = Δ_fibre / a_base where:
  Δ_fibre = gap in fibre units = 1.263
  a_base = base step length = a_lattice

So c = a_DS/a_lattice = a_fibre/a_base.

The fibre-to-base ratio for the twistor fibration:
the twistor fibration of CP³ → S⁴ is a Riemannian submersion
with totally geodesic fibres. The O'Neill tensor relates the
fibre and base metrics. For CP³ with FS metric normalised to
have sectional curvatures in [1,4]:
  fibre radius = 1 (CP¹ with curvature 4)
  base radius = 2 (S⁴ with curvature 1)

So the fibre-to-base ratio is 1/2.

At β=2.3: a_lattice ≈ 0.17 fm. One DS fibre step corresponds to
a_fibre = ... we need the absolute scale.
"""
import numpy as np

# The CP³ Fubini-Study metric with standard normalization:
# ds² = 4(|dz|²|z|² - |<z̄,dz>|²) / |z|⁴
# (factor of 4 for holomorphic sectional curvature = 1)
#
# The fibration CP³ → S⁴:
# - Base S⁴ has sectional curvature 1 (radius 1)
# - Fibre CP¹ has curvature 4 (radius 1/2)
# - The fibration is a Riemannian submersion
#
# The horizontal lift of a unit tangent vector on S⁴ has
# length 1 in CP³. The vertical fibre vector has length
# proportional to the fibre radius.

# In the physical theory:
# One lattice spacing a on S⁴ corresponds to one horizontal step.
# One DS step corresponds to one fibre step.
# The conversion factor c = fibre step / horizontal step.

# For a Riemannian submersion with totally geodesic fibres:
# The O'Neill A-tensor relates horizontal and vertical.
# For the Hopf fibration S³ → S² (a simpler case):
#   fibre S¹ circumference = 2π
#   base S² area = 4π
#   The fibration has A-tensor related to the curvature.

# For CP³ → S⁴ via the Penrose fibration:
# The fibre is CP¹ ≅ S² with the round metric scaled by factor 1/2.
# One "unit" fibre step is a great circle arc of length π on S² (diameter).
# One "unit" base step is an arc of length π on S⁴.
# The ratio: fibre_step/base_step = (π/2)/(π) = 1/2
# (if we use the full diameter as the step)

# But the DS step is NOT the full diameter. It's a contraction by λ₀.
# The DS step moves the state by a distance proportional to 1-λ₀ = 0.717
# in units of the current perturbation.

# At the equilibrium, the characteristic perturbation scale is set by
# the eigenvalue structure. The two eigenvalues λ₀=0.2829 and λ₁=0.2813
# are nearly degenerate, suggesting the equilibrium is isotropic.

# A different approach: use the ENERGY-LENGTH correspondence.
# In QFT: E = ℏc/λ, where λ is the Compton wavelength.
# In DS: Δ = 1.263 is the gap in "inverse fibre steps."
# In lattice: am = 0.71 is the gap in "inverse lattice spacings."
# The ratio: c = Δ/(am) = 1.263/0.71 = 1.78.
# This IS the conversion factor, by definition.
# It equals fibre_step/lattice_spacing.

# Can we get 1.78 from pure geometry?

# Attempt: the CP³ FS volume is π³/6.
# The S⁴ volume is 8π²/3.
# The fibre CP¹ volume is π (area of S² with curvature 4).
# By the submersion formula: Vol(CP³) = Vol(S⁴) × Vol(CP¹) / (gauge volume)
# π³/6 = (8π²/3) × π / V_gauge
# V_gauge = (8π²/3) × π / (π³/6) = (8π³/3) × (6/π³) = 16
# Hmm, that doesn't simplify nicely.

# Let's try dimensional analysis:
# c = a_DS/a_latt = (fibre scale) / (base scale)
# For CP³ with FS curvature K=4 on fibres, K=1 on base:
# fibre radius = 1/√K_fibre = 1/2
# base radius = 1/√K_base = 1
# Ratio = 1/2.
# But c = 1.78, not 0.5.

# What if the "step" isn't one radius but one geodesic unit?
# One geodesic step in the fibre: the transfer operator eigenvalue λ₀
# gives the step length as -ln(λ₀) = 1.263 in FS units.
# One geodesic step in the base: the lattice spacing a in physical units.
# The conversion: 1 FS fibre unit = c × a base units.

# From the FS metric on the fibre CP¹:
# The geodesic distance between the north pole and a point at
# latitude θ is θ (in FS units with curvature 4).
# The maximum distance is π/2 (to the south pole).

# The DS step of Δ=1.263 in FS units covers what fraction of the fibre?
# 1.263 / (π/2) = 1.263/1.571 = 0.804 = 80.4% of the fibre diameter.
# That's almost the entire fibre per step!

print("FIBRE COVERAGE PER DS STEP:")
print(f"  Δ = 1.263 FS units")
print(f"  CP¹ diameter = π/2 = {np.pi/2:.4f} FS units")
print(f"  Δ / (π/2) = {1.263/(np.pi/2):.4f} = {1.263/(np.pi/2)*100:.1f}% of fibre")
print()

# This means one DS step covers ~80% of the fibre.
# The lattice spacing covers some distance in the base.
# The conversion c = 1.78 means:
# 1 DS fibre step (1.263 FS) = 1.78 lattice base steps (1.78 × a)

# In physical terms:
# 1 DS step = 1.263/c = 0.71 lattice units = am (lattice mass!)
# This is just the definition: am = Δ/c.

# The geometric meaning of c:
# c = Δ/(am) = (FS fibre gap)/(lattice base gap)
# = (fibre step/fibre gap)/(base step/base gap)
# = (1/λ₀)/(1/(am)) ... this is circular.

# I think the honest conclusion is:
# c = 1.78 relates the DS fibre scale to the lattice base scale.
# It can be MEASURED (from lattice data) but its geometric DERIVATION
# requires knowing the absolute physical scale of the twistor fibration,
# which is set by the coupling constant g (or equivalently β).

print("HONEST CONCLUSION:")
print("  c = a_DS/a_lattice = Δ/(am) = 1.263/0.71 = 1.78")
print("  This is the fibre-to-base scale ratio.")
print("  ")
print("  It CANNOT be derived from pure CP³ geometry because")
print("  it depends on the COUPLING CONSTANT β (which sets a_lattice).")
print("  Under asymptotic scaling: a(β) ~ exp(-β/(4b₀))")
print("  So c(β) = Δ/am(β) grows with β.")
print("  ")
print("  At β=2.3: c=1.78")
print("  At β→∞ (continuum): c→∞")
print("  ")
print("  The conversion factor is NOT a property of the DS framework.")
print("  It's a property of the LATTICE REGULARISATION.")
print("  Under Path (B), the DS framework IS the continuum theory;")
print("  the lattice is a computational tool, not the definition.")
print("  c is the ratio of two DIFFERENT things:")
print("    - Δ from the DS algebra (intrinsic, g-independent)")
print("    - am from lattice Monte Carlo (β-dependent)")
print("  ")
print("  The physical mass m = Δ × (energy scale of one DS step)")
print("  The energy scale is set by the twistor fibration's")
print("  embedding into physical spacetime. This is analogous to")
print("  the lattice spacing being set by β — it's an external scale,")
print("  not derivable from the theory alone.")
print("  ")
print("  Bottom line: c is PHYSICAL INPUT, not a derivable quantity.")
print("  The DS framework predicts Δ=1.263 in natural units.")
print("  Converting to GeV requires an energy scale (like ΛQCD).")
print("  This is the SAME situation as in lattice QCD:")
print("  the theory predicts dimensionless ratios; one physical")
print("  input (ΛQCD or a string tension) sets the scale.")
