"""
Penrose transform: Step 3.

Compute the volume fraction ANALYTICALLY.
"""
import numpy as np
from scipy import integrate

# The volume element of CP³ in radial affine coordinates:
# dV_FS = r⁵/(1+r²)⁴ dr × Vol(S⁵)
# (radial part only, angular part is constant)

# Total volume (radial integral):
# ∫₀^∞ r⁵/(1+r²)⁴ dr
# Substitution u = 1+r², du = 2r dr, r² = u-1:
# = ∫₁^∞ (u-1)² / (2u⁴) du
# = (1/2) ∫₁^∞ (u⁻² - 2u⁻³ + u⁻⁴) du
# = (1/2) [-u⁻¹ + u⁻² - u⁻³/3]₁^∞
# = (1/2) [0 - (-1 + 1 - 1/3)]
# = (1/2)(1/3) = 1/6

# Ball volume (R² = H³-1, so 1+R² = H³):
# = (1/2) [-1/(1+R²) + 1/(1+R²)² - 1/(3(1+R²)³) + 1/3]
# = (1/2) [-1/H³ + 1/H⁶ - 1/(3H⁹) + 1/3]

# Volume fraction:
# f(H) = Vol_ball / Vol_total
#       = (1/6)⁻¹ × (1/2)[-1/H³ + 1/H⁶ - 1/(3H⁹) + 1/3]
#       = 3[1/3 - 1/H³ + 1/H⁶ - 1/(3H⁹)]
#       = 1 - 3/H³ + 3/H⁶ - 1/H⁹

# Factor: this is the BINOMIAL EXPANSION of (1 - 1/H³)³ !
# f(H) = (1 - 1/H³)³ = ((H³-1)/H³)³

print("=" * 60)
print("EXACT VOLUME FRACTION")
print("=" * 60)
print()
print("The Fubini-Study volume fraction of CP³ inside the")
print("Born floor ball B(0, √(H³-1)) is:")
print()
print("  f(H) = (1 - 1/H³)³ = ((H³-1)/H³)³")
print()

for H in [2, 3, 4, 5]:
    f_exact = (1 - 1/H**3)**3
    # Numerical check
    def integrand(r):
        return r**5 / (1 + r**2)**4
    R = np.sqrt(H**3 - 1)
    vol_ball, _ = integrate.quad(integrand, 0, R)
    vol_total, _ = integrate.quad(integrand, 0, np.inf)
    f_numerical = vol_ball / vol_total
    print(f"  H={H}: exact = {f_exact:.10f}, numerical = {f_numerical:.10f}, diff = {abs(f_exact-f_numerical):.2e}")

print()

H = 3
f = (1 - 1/H**3)**3
print(f"At H=3:")
print(f"  f(3) = (26/27)³ = {f:.10f}")
print(f"  = (1 - Born_floor)³")
print(f"  Born_floor = 1/27 = {1/27:.10f}")
print(f"  Excluded fraction = 1 - f = {1-f:.10f}")
print(f"  = 1 - (26/27)³ = {1-(26/27)**3:.10f}")
print()

# Simplify: 1 - (26/27)³ = 1 - 17576/19683 = 2107/19683
excluded_num = 27**3 - 26**3
excluded_den = 27**3
print(f"  Excluded = {excluded_num}/{excluded_den}")
print(f"  = {excluded_num/excluded_den:.10f}")
print()

# Factor 27³ - 26³ = (27-26)(27²+27·26+26²) = 1·(729+702+676) = 2107
print(f"  27³ - 26³ = (27-26)(27² + 27·26 + 26²)")
print(f"            = 1 × ({27**2} + {27*26} + {26**2})")
print(f"            = {27**2 + 27*26 + 26**2}")
print(f"  Check: {27**3 - 26**3} = {27**2 + 27*26 + 26**2}")
print()

print("=" * 60)
print("THE RESULT IS EXACT AND ALGEBRAIC")
print("=" * 60)
print()
print("The Born floor at 1/H³ excludes EXACTLY")
print()
print("  (H³ - (H³-1)³/H⁶) of the Fubini-Study volume of CP³")
print()
print("At H=3: excluded fraction = 2107/19683 ≈ 10.70%")
print()
print("This is a theorem, not a numerical observation.")
print("It follows from the integral ∫r⁵/(1+r²)⁴ dr")
print("evaluated at r = √(H³-1), which is elementary calculus.")
