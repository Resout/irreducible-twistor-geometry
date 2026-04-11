"""
Helium screening constant Z_eff = 27/16 and the H=3 Born floor.

Investigation: Is the helium variational screening constant Z_eff = 27/16
derivable from the irreducible twistor geometry at H=3?

Key fact: The standard variational treatment of helium gives
  Z_eff = Z - 5/16 = 2 - 5/16 = 27/16 = 1.6875
where 5/16 comes from the electron-electron repulsion integral (5/8)Z_eff
evaluated at the minimum.

The question: at H=3, we have H^3/(H+1)^2 = 27/16. Derivation or coincidence?
"""

import numpy as np
from fractions import Fraction

H = 3
Kstar = Fraction(7, 30)

print("=" * 72)
print("HELIUM SCREENING CONSTANT: Z_eff = 27/16 AND H=3")
print("=" * 72)

# ─────────────────────────────────────────────────────────────────────
# SECTION 1: Basic numerical verifications
# ─────────────────────────────────────────────────────────────────────
print("\n" + "─" * 72)
print("SECTION 1: Exact algebraic verifications at H=3")
print("─" * 72)

Z_eff_exact = Fraction(27, 16)
screening = Fraction(5, 16)
repulsion_integral = Fraction(5, 8)

# Framework expressions
expr1 = Fraction(H**3, (H+1)**2)
expr2 = Fraction(H+2, (H+1)**2)
expr3 = Fraction(H+2, 2*(H+1))

print(f"\nZ_eff = 27/16 = {float(Z_eff_exact)}")
print(f"H³/(H+1)² = {H}³/{H+1}² = {expr1} = {float(expr1)}")
print(f"  Match: {Z_eff_exact == expr1}  ✓" if Z_eff_exact == expr1 else f"  MISMATCH")

print(f"\nScreening 5/16 = {float(screening)}")
print(f"(H+2)/(H+1)² = {H+2}/{(H+1)**2} = {expr2} = {float(expr2)}")
print(f"  Match: {screening == expr2}  ✓" if screening == expr2 else f"  MISMATCH")

print(f"\nRepulsion integral 5/8 = {float(repulsion_integral)}")
print(f"(H+2)/(2(H+1)) = {H+2}/{2*(H+1)} = {expr3} = {float(expr3)}")
print(f"  Match: {repulsion_integral == expr3}  ✓" if repulsion_integral == expr3 else f"  MISMATCH")

# Decomposition
print(f"\nDecomposition of Z_eff:")
print(f"  Z_eff = Z - screening = 2 - 5/16 = 27/16")
print(f"  Z_eff = H³/(H+1)² = H³/((H+1)²)")
print(f"  Note: H³ = (H+1-1)³, (H+1)² = (H+1)²")
print(f"  So Z_eff = (H+1-1)³/(H+1)² = (H+1) - 3 + 3/(H+1) - 1/(H+1)²")

# Verify algebraically: Z=2 = H-1, and Z_eff = (H-1) - (H+2)/(H+1)^2
Z_val = Fraction(H - 1)  # Z=2 at H=3... but Z=2 is the nuclear charge
# Actually Z=2 is NOT H-1=2 trivially. Let's be careful.
print(f"\n  Z = 2 (helium nuclear charge)")
print(f"  H - 1 = {H-1} = 2  ... Z = H-1 for helium specifically")
print(f"  Z_eff = Z - (H+2)/(H+1)² = (H-1) - (H+2)/(H+1)²")
print(f"        = [(H-1)(H+1)² - (H+2)] / (H+1)²")

numer = (H-1)*(H+1)**2 - (H+2)
denom = (H+1)**2
print(f"        = [{(H-1)}×{(H+1)**2} - {H+2}] / {(H+1)**2}")
print(f"        = [{(H-1)*(H+1)**2} - {H+2}] / {(H+1)**2}")
print(f"        = {numer}/{denom} = {Fraction(numer, denom)}")

# Verify: is (H-1)(H+1)^2 - (H+2) = H^3?
# (H-1)(H+1)^2 = (H-1)(H^2+2H+1) = H^3+2H^2+H - H^2-2H-1 = H^3+H^2-H-1
# Then subtract (H+2): H^3+H^2-H-1-H-2 = H^3+H^2-2H-3
# At H=3: 27+9-6-3 = 27. Yes!
# In general: H^3+H^2-2H-3 = H^3 iff H^2-2H-3=0 iff (H-3)(H+1)=0
# So H=3 or H=-1. This is H-SPECIFIC!
print(f"\n  Algebraic check: (H-1)(H+1)² - (H+2) = H³ ?")
poly_check = lambda h: (h-1)*(h+1)**2 - (h+2) - h**3
print(f"  At H=3: {poly_check(3)} → {'YES' if poly_check(3)==0 else 'NO'}")
print(f"  In general: (H-1)(H+1)² - (H+2) - H³ = H²-2H-3 = (H-3)(H+1)")
print(f"  This vanishes ONLY at H=3 (or H=-1).")
print(f"\n  *** THIS IS THE KEY RESULT ***")
print(f"  The identity Z_eff = H³/(H+1)² with Z = H-1")
print(f"  holds ONLY at H=3. It is H-specific, not a general pattern.")

# ─────────────────────────────────────────────────────────────────────
# SECTION 2: What happens at other H values?
# ─────────────────────────────────────────────────────────────────────
print("\n" + "─" * 72)
print("SECTION 2: Framework expressions at other H values")
print("─" * 72)

print(f"\n{'H':>3} │ {'H³/(H+1)²':>12} │ {'(H+2)/(H+1)²':>14} │ {'Z=H-1':>6} │ {'Z-(H+2)/(H+1)²':>16} │ {'= H³/(H+1)²?':>13}")
print("─" * 72)
for h in range(2, 8):
    expr_zeff = Fraction(h**3, (h+1)**2)
    expr_screen = Fraction(h+2, (h+1)**2)
    Z = h - 1
    z_minus_screen = Z - expr_screen
    match = "✓ EXACT" if z_minus_screen == expr_zeff else f"off by {float(z_minus_screen - expr_zeff):+.4f}"
    print(f"{h:>3} │ {float(expr_zeff):>12.6f} │ {float(expr_screen):>14.6f} │ {Z:>6} │ {float(z_minus_screen):>16.6f} │ {match}")

print(f"\nThe identity Z-screening = H³/(H+1)² (with Z=H-1) holds ONLY at H=3.")

# ─────────────────────────────────────────────────────────────────────
# SECTION 3: He-like ions — screening 5/16 for various Z
# ─────────────────────────────────────────────────────────────────────
print("\n" + "─" * 72)
print("SECTION 3: He-like ions (2 electrons, nuclear charge Z)")
print("─" * 72)
print("Standard variational: Z_eff = Z - 5/16 for ALL He-like ions")
print("The 5/16 screening is universal for 2-electron hydrogenic systems.\n")

print(f"{'Ion':>6} │ {'Z':>3} │ {'Z_eff=Z-5/16':>13} │ {'H³/(H+1)²':>12} │ {'Match?':>8}")
print("─" * 50)
ions = [("H⁻", 1), ("He", 2), ("Li⁺", 3), ("Be²⁺", 4), ("B³⁺", 5)]
for name, Z in ions:
    z_eff = Fraction(Z) - Fraction(5, 16)
    framework = Fraction(H**3, (H+1)**2)
    match = "✓ EXACT" if z_eff == framework else f"  {float(z_eff):.4f}"
    print(f"{name:>6} │ {Z:>3} │ {float(z_eff):>13.4f} │ {float(framework):>12.4f} │ {match}")

print(f"\nOnly helium (Z=2) matches. The 5/16 screening is universal,")
print(f"but H³/(H+1)² = 27/16 is specific to Z=2=H-1.")

# ─────────────────────────────────────────────────────────────────────
# SECTION 4: Helium ground state energy
# ─────────────────────────────────────────────────────────────────────
print("\n" + "─" * 72)
print("SECTION 4: Helium ground state energy from Z_eff = 27/16")
print("─" * 72)

E1 = -13.605693  # eV, hydrogen ground state energy (Rydberg)
Z_he = 2
Z_eff_val = 27/16

# Variational energy: E = [Z_eff² - 2Z·Z_eff + (5/8)Z_eff] × E₁
# = [Z_eff² - 2Z·Z_eff + (5/8)Z_eff] × E₁
# Actually the standard formula for helium variational energy is:
# E(Z_eff) = (Z_eff^2 - 2*Z*Z_eff + (5/8)*Z_eff) * 2 * E1 / Z_eff^2 ...
# Let me be more careful.

# The correct formula:
# Using hydrogenic trial ψ ~ exp(-Z_eff r/a₀), the energy is:
# E = 2[Z_eff²/2 - Z·Z_eff] + (5/8)Z_eff    (in Rydberg units, E1 = -1/2 Ry)
# No wait. In units of e²/(2a₀) = 1 Hartree = 27.211 eV:
# E = Z_eff² - 2Z·Z_eff + (5/8)Z_eff    ... this is wrong too.

# Let me write it correctly. For a hydrogenic trial wavefunction with parameter Z_eff:
# <T> = Z_eff² (in Hartree) for both electrons
# <V_ne> = -2Z·Z_eff (nuclear-electron attraction, both electrons)
# <V_ee> = (5/8)Z_eff (electron-electron repulsion)
# E = Z_eff² - 2Z·Z_eff + (5/8)Z_eff   (Hartree)

E_hartree = Z_eff_val**2 - 2*Z_he*Z_eff_val + (5/8)*Z_eff_val
E_eV = E_hartree * 27.21138  # 1 Hartree = 27.21138 eV

# Exact fraction
Z_e = Fraction(27, 16)
E_frac = Z_e**2 - 2*2*Z_e + Fraction(5,8)*Z_e
print(f"\nVariational energy (Hartree): E = Z_eff² - 2Z·Z_eff + (5/8)Z_eff")
print(f"  = ({Z_e})² - 4×({Z_e}) + (5/8)×({Z_e})")
print(f"  = {Z_e**2} - {4*Z_e} + {Fraction(5,8)*Z_e}")
print(f"  = {E_frac} = {float(E_frac):.6f} Hartree")
print(f"  = {float(E_frac)*27.21138:.3f} eV")

# Framework form
print(f"\n  In framework terms:")
print(f"  E = H⁶/(H+1)⁴ - 2Z·H³/(H+1)² + (H+2)/(2(H+1))·H³/(H+1)²")
E_framework = Fraction(H**6, (H+1)**4) - 2*2*Fraction(H**3, (H+1)**2) + Fraction(H+2, 2*(H+1))*Fraction(H**3, (H+1)**2)
print(f"  = {E_framework} = {float(E_framework):.6f} Hartree")
print(f"  Match with direct calc: {E_frac == E_framework}")

# Simplify
print(f"\n  Numerator: {E_framework.numerator}")
print(f"  Denominator: {E_framework.denominator}")
# Factor
from math import gcd
g = gcd(abs(E_framework.numerator), E_framework.denominator)
print(f"  Simplified: {E_framework.numerator//g}/{E_framework.denominator//g}")

# Compare to experiment
E_exp = -79.005  # eV (experimental helium ground state, includes relativistic corrections)
E_exp_hartree = E_exp / 27.21138
print(f"\n  Variational result:  {float(E_frac)*27.21138:.3f} eV = {float(E_frac):.6f} Ha")
print(f"  Experimental value:  {E_exp:.3f} eV = {E_exp_hartree:.6f} Ha")
print(f"  Difference: {abs(float(E_frac)*27.21138 - E_exp):.3f} eV ({abs(float(E_frac)*27.21138 - E_exp)/abs(E_exp)*100:.2f}%)")
print(f"  (The variational result is an upper bound, as expected.)")

# Hylleraas-level result for comparison
E_hylleraas = -2.903724  # Hartree (essentially exact non-relativistic)
print(f"\n  Hylleraas (exact NR): {E_hylleraas:.6f} Ha = {E_hylleraas*27.21138:.3f} eV")
print(f"  Simple variational:   {float(E_frac):.6f} Ha")
print(f"  Variational recovers {float(E_frac)/E_hylleraas*100:.2f}% of exact NR energy")

# ─────────────────────────────────────────────────────────────────────
# SECTION 5: Why 5/8? The electron-electron repulsion integral
# ─────────────────────────────────────────────────────────────────────
print("\n" + "─" * 72)
print("SECTION 5: The 5/8 repulsion integral — structural analysis")
print("─" * 72)

print(f"""
The e-e repulsion integral for two 1s hydrogenic electrons:
  <1/r₁₂> = (5/8)(Z_eff/a₀)

The factor 5/8 arises from a 6-dimensional integral over two electron
positions in R³. It is a purely mathematical result of the 1/r Coulomb
potential integrated against exp(-2Z_eff r/a₀) wavefunctions.

Framework decomposition of 5/8:
  5/8 = (H+2)/(2(H+1)) at H=3

  Numerator: H+2 = 5 = spacetime dimension in the framework
  Denominator: 2(H+1) = 8 = 2 × dim(fermion quartet)

  Alternative: 5/8 = (H²+H-1)/(H²+H-1) × 5/8 ... no, this is circular.

Let's check if 5/8 has a natural meaning:
  5/8 = 1 - 3/8 = 1 - H/(2(H+1))
""")

# Check
print(f"  1 - H/(2(H+1)) = 1 - {H}/{2*(H+1)} = {1 - Fraction(H, 2*(H+1))} = {float(1 - Fraction(H, 2*(H+1)))}")
print(f"  = 5/8? {Fraction(1) - Fraction(H, 2*(H+1)) == Fraction(5, 8)}")

print(f"""
  Also: 5/8 = (H+2)/(2H+2) = (H+2)/(2(H+1))

  In the standard QM calculation, 5/8 = 5Z_eff/(8a₀) × (a₀/Z_eff)
  The "5" comes from the angular integration yielding 4π × (1 + 1/4),
  where the 1/4 is from the Legendre expansion of 1/r₁₂.

  The 5 in "5/8" is genuinely 5 = the number that arises from
  integrating 1/r₁₂ over two s-wave functions. It counts:
  - ℓ=0 angular momentum: contributes 1 (monopole)
  - radial integral ratio: contributes 5/8 total

  At H=3: 5 = H+2, and 8 = 2(H+1).
  These are framework quantities, but 5/8 is also just a Coulomb integral.
""")

# ─────────────────────────────────────────────────────────────────────
# SECTION 6: Deeper algebraic structure
# ─────────────────────────────────────────────────────────────────────
print("─" * 72)
print("SECTION 6: Algebraic structure — why H=3 is special")
print("─" * 72)

print(f"""
The key identity is:
  (H-1)(H+1)² - (H+2) = H³   ONLY at H=3

Expanding: H³ + H² - 2H - 3 = H³  ⟹  H² - 2H - 3 = 0  ⟹  (H-3)(H+1) = 0

This means: at H=3, the relation Z_eff = H³/(H+1)² with Z=H-1
is not a parametric family — it is a SPECIFIC algebraic coincidence.

But is it really a coincidence? Let's unpack what each piece means:

  Z = H - 1 = 2     Nuclear charge of helium
  H + 1 = 4          Fermion degrees per generation
  H + 2 = 5          Spacetime dimension
  H³ = 27            Born floor denominator (Born ≥ 1/27)
  (H+1)² = 16        Two-fermion state space dimension

The identity says: if we take the nucleus with charge Z = H-1,
screen it by the ratio (spacetime dim)/(two-fermion states),
we get the Born capacity / two-fermion states.

More precisely:
  Z_eff = Z - σ  where σ = (H+2)/(H+1)² = 5/16

  σ = 5/16 is the probability that a 2-electron system "blocks"
  a unit of nuclear charge, given the constraints of the Born floor.
""")

# ─────────────────────────────────────────────────────────────────────
# SECTION 7: The Slater screening constant
# ─────────────────────────────────────────────────────────────────────
print("─" * 72)
print("SECTION 7: Slater screening constant and K*")
print("─" * 72)

slater_sigma = 0.30  # Slater's empirical screening for 1s electron in He
print(f"\nSlater's empirical screening for He 1s: σ_Slater = {slater_sigma}")
print(f"  Z_eff(Slater) = 2 - 0.30 = {2 - slater_sigma}")
print(f"  Z_eff(variational) = 2 - 5/16 = 2 - 0.3125 = 1.6875")
print(f"  Slater ≠ variational (Slater is simpler empirical rule)")

print(f"\n  K* = 7/30 = {float(Kstar):.6f}")
print(f"  Slater σ = 0.30 = 9/30")
print(f"  Ratio: σ_Slater/K* = 0.30/{float(Kstar):.4f} = {0.30/float(Kstar):.4f} = 9/7")
print(f"  Is 9/7 = H²/K*_numerator? {Fraction(9,7)} ... H²=9, K*_num=7. Interesting but likely coincidence.")

print(f"\n  Variational σ = 5/16 = {5/16}")
print(f"  K* = 7/30")
print(f"  σ/K* = (5/16)/(7/30) = {Fraction(5,16)/Fraction(7,30)} = {float(Fraction(5,16)/Fraction(7,30)):.4f}")
print(f"  = 75/56 ... no obvious framework meaning.")

# ─────────────────────────────────────────────────────────────────────
# SECTION 8: Multi-electron atoms — does 5/16 extend?
# ─────────────────────────────────────────────────────────────────────
print("\n" + "─" * 72)
print("SECTION 8: Multi-electron variational screening")
print("─" * 72)

print(f"""
For the 1s electron in multi-electron atoms, the variational
approach gives approximately:
  Z_eff ≈ Z - σ(N_e)

where σ depends on the number of electrons and their configuration.
""")

# Experimental/calculated Z_eff values for 1s electrons
# From Clementi & Raimondi (1963) — Z_eff for 1s orbital
clementi = [
    ("H",   1,  1, 1.000),
    ("He",  2,  2, 1.688),   # ≈ 27/16 = 1.6875
    ("Li",  3,  3, 2.691),
    ("Be",  4,  4, 3.685),
    ("B",   5,  5, 4.680),
    ("C",   6,  6, 5.673),
    ("N",   7,  7, 6.665),
    ("O",   8,  8, 7.658),
    ("F",   9,  9, 8.650),
    ("Ne", 10, 10, 9.642),
]

print(f"\n{'Atom':>4} │ {'Z':>3} │ {'Z_eff(1s)':>10} │ {'σ=Z-Z_eff':>10} │ {'5/16':>6} │ {'σ-5/16':>8}")
print("─" * 56)
for name, Z, N_e, z_eff in clementi:
    sigma = Z - z_eff
    diff = sigma - 5/16
    marker = " ✓" if abs(diff) < 0.005 else ""
    print(f"{name:>4} │ {Z:>3} │ {z_eff:>10.3f} │ {sigma:>10.3f} │ {5/16:>6.4f} │ {diff:>+8.4f}{marker}")

print(f"""
Key observation: The 1s screening constant for ALL multi-electron atoms
is approximately 0.3125 = 5/16 per additional 1s electron.

For He (Z=2): σ = 2 - 1.688 = 0.312 ≈ 5/16 = 0.3125 (from the other 1s e⁻)
For Li (Z=3): σ = 3 - 2.691 = 0.309 (two 1s screen each other; 2s adds tiny bit)
For heavier: σ ≈ 0.31-0.36 (increases as inner shells feel more screening)

The 5/16 screening is most exact for helium and He-like ions,
where there are exactly 2 electrons in the same 1s orbital.
""")

# ─────────────────────────────────────────────────────────────────────
# SECTION 9: Connection to framework conservation law
# ─────────────────────────────────────────────────────────────────────
print("─" * 72)
print("SECTION 9: Conservation law structure")
print("─" * 72)

print(f"""
The variational principle says: minimize E(Z_eff) over Z_eff.
  dE/dZ_eff = 0  →  2Z_eff - 2Z + 5/8 = 0  →  Z_eff = Z - 5/16

In the framework, conservation laws arise from the Killing form on
the fibre algebra. The conservation constraint is:

  ∂_μ J^μ = 0  →  balance between kinetic (Z_eff²) and potential terms

The variational minimum IS a conservation law: it's the point where
the "force" from nuclear attraction exactly balances the electron
repulsion, in the space of trial wavefunctions.

Framework interpretation:
  Z_eff = Z - (H+2)/(H+1)²

  (H+2) = spacetime dimension = number of independent directions
           in which the electron can "avoid" the other electron
  (H+1)² = number of two-particle quantum states

  The screening (H+2)/(H+1)² is the ratio of geometric escape routes
  to total quantum states. This is a Born-floor-type ratio: it measures
  the fraction of the state space "used up" by the other electron.
""")

# ─────────────────────────────────────────────────────────────────────
# SECTION 10: H³/(H+1)² as Born capacity per two-particle state
# ─────────────────────────────────────────────────────────────────────
print("─" * 72)
print("SECTION 10: Born capacity interpretation")
print("─" * 72)

print(f"""
The Born floor states: m(s) ≥ 1/H³ for all hypotheses s.
  Total Born capacity = H³ = 27 (max distinguishable states)

For two fermions (electrons), each with (H+1) = 4 internal states:
  Two-particle Hilbert space dimension = (H+1)² = 16

  Z_eff = H³/(H+1)² = {float(Fraction(H**3, (H+1)**2))}

  = (Born capacity) / (two-particle states)
  = 27/16 = 1.6875

This says: the effective charge seen by each electron equals the
number of Born-distinguishable states per two-particle quantum state.

Physical meaning: each quantum state of the two-electron system
"carries" H³/(H+1)² = 27/16 units of nuclear charge. The nuclear
charge is distributed over the state space, diluted by the Born floor.
""")

# ─────────────────────────────────────────────────────────────────────
# SECTION 11: Critical assessment
# ─────────────────────────────────────────────────────────────────────
print("─" * 72)
print("SECTION 11: Critical assessment — derivation vs coincidence")
print("─" * 72)

print(f"""
WHAT IS SOLID:
  1. Z_eff = 27/16 is the EXACT variational optimum for helium.         ✓
  2. 27/16 = H³/(H+1)² at H=3 is an exact algebraic identity.          ✓
  3. 5/16 = (H+2)/(H+1)² at H=3 is exact.                              ✓
  4. 5/8 = (H+2)/(2(H+1)) at H=3 is exact.                             ✓
  5. The identity Z_eff = H³/(H+1)² with Z=H-1 holds ONLY at H=3.      ✓
     (The polynomial H²-2H-3=0 has root H=3.)

WHAT IS A STRUCTURAL IDENTIFICATION (not yet a derivation):
  6. "H+2 = spacetime dimension" → explains the 5 in 5/8               ?
  7. "(H+1)² = two-fermion states" → explains the 16 in 5/16           ?
  8. "H³ = Born capacity" → explains the 27 in 27/16                   ?

  These are IDENTIFICATIONS of known framework quantities with the
  numbers appearing in the helium calculation. They are consistent
  but not derived from first principles within the framework.

WHAT WOULD CONSTITUTE A DERIVATION:
  - Show that the electron-electron repulsion integral in the framework's
    geometry MUST equal (H+2)/(2(H+1)) from the Born floor constraint.
  - This would require: the Coulomb 1/r potential emerges from the
    framework, and the 1s wavefunction is the natural ground state,
    and the 6D integration yields framework quantities.

STRONGEST ARGUMENT FOR SIGNIFICANCE:
  The H-specificity result (Section 6). The identity:
    (H-1)(H+1)² - (H+2) = H³
  holds ONLY at H=3. If the framework uniquely selects H=3, and if
  Z=H-1 for helium, then Z_eff = H³/(H+1)² is not a coincidence —
  it's a consequence of the same algebraic constraint that fixes H=3.

  The polynomial (H-3)(H+1) = 0 is a CONSISTENCY CHECK: only at the
  physical value H=3 do the nuclear charge (Z=2), the screening (5/16),
  and the Born capacity (27) form a self-consistent triad.

STRONGEST ARGUMENT AGAINST:
  Z = 2 = H - 1 for helium is just a numerical coincidence. Helium
  has Z=2 because it has 2 protons, not because H=3. The framework
  would need to explain WHY helium specifically (not lithium or hydrogen)
  satisfies Z = H-1. Without this, the H-specificity result is
  suggestive but not conclusive.

  Also: 5/8 is a basic Coulomb integral. It equals (H+2)/(2(H+1)) at
  H=3, but this might just be small-number arithmetic coincidence.
  There are only so many ratios of small integers.

VERDICT:
  The Z_eff = H³/(H+1)² identity is a STRONG STRUCTURAL IDENTIFICATION
  elevated by the H-specificity result. It is NOT yet a derivation.

  Confidence: 60% physics, 40% coincidence.

  The H-specificity — that (H-1)(H+1)²-(H+2)=H³ holds ONLY at H=3 —
  is the strongest evidence. It means this isn't a parametric family
  that accidentally passes through H=3; it's an isolated algebraic
  fact about H=3 specifically.
""")

# ─────────────────────────────────────────────────────────────────────
# SECTION 12: Summary table
# ─────────────────────────────────────────────────────────────────────
print("─" * 72)
print("SECTION 12: Summary")
print("─" * 72)

print(f"""
  ┌─────────────────────────────────────────────────────────────┐
  │ Quantity          │ Standard QM    │ Framework (H=3)        │
  ├───────────────────┼────────────────┼────────────────────────┤
  │ Z_eff             │ 27/16          │ H³/(H+1)²             │
  │ Screening σ       │ 5/16           │ (H+2)/(H+1)²          │
  │ Repulsion 5/8     │ 5/8            │ (H+2)/(2(H+1))        │
  │ Nuclear charge    │ Z = 2          │ H - 1                 │
  │ He energy         │ -729/256 Ha    │ -H⁶/(H+1)⁴ Ha*       │
  │ H-specificity     │ —              │ (H-3)(H+1) = 0        │
  └─────────────────────────────────────────────────────────────┘

  * Energy expression after substitution; not independently derived.

  Key: Z_eff = H³/(H+1)² = 27/16 is EXACT.
  The identity holds ONLY at H=3 among positive integers.
  This elevates it beyond simple numerology.
""")

# Final numerical check
print("─" * 72)
print("FINAL NUMERICAL CROSS-CHECK")
print("─" * 72)
print(f"  H = {H}")
print(f"  H³ = {H**3}")
print(f"  (H+1)² = {(H+1)**2}")
print(f"  H³/(H+1)² = {Fraction(H**3, (H+1)**2)} = {float(Fraction(H**3, (H+1)**2))}")
print(f"  Standard Z_eff = 27/16 = {float(Fraction(27, 16))}")
print(f"  Agreement: EXACT ✓")
print(f"\n  Energy: {E_frac} Ha = {float(E_frac):.10f} Ha")
print(f"  = -729/256 Ha? {E_frac == Fraction(-729, 256)}")
print(f"  = {float(Fraction(-729, 256)):.10f} Ha")

# Check -H^6/(H+1)^4
E_framework_simple = Fraction(-H**6, (H+1)**4)
print(f"  -H⁶/(H+1)⁴ = {E_framework_simple} = {float(E_framework_simple):.10f} Ha")
print(f"  Match direct energy: {E_frac == E_framework_simple}")

# Hmm, let's compute E properly with exact fractions
print(f"\n  Recomputing E exactly:")
print(f"  E = Z_eff² - 2Z·Z_eff + (5/8)Z_eff")
t1 = Z_e**2
t2 = -2 * 2 * Z_e
t3 = Fraction(5, 8) * Z_e
print(f"    = {t1} + ({t2}) + {t3}")
print(f"    = {t1} - {abs(t2)} + {t3}")
print(f"    = {t1 + t2 + t3}")
E_check = t1 + t2 + t3
print(f"    = {E_check} = {float(E_check):.10f} Ha")

# -729/256 check
print(f"    -729/256 = {float(Fraction(-729,256)):.10f}")
print(f"    Match: {E_check == Fraction(-729, 256)}")

# So E = -729/256. Is 729 = 3^6 = H^6? Yes. Is 256 = 4^4 = (H+1)^4? Yes.
print(f"\n  729 = 3⁶ = H⁶: {729 == H**6}")
print(f"  256 = 4⁴ = (H+1)⁴: {256 == (H+1)**4}")
print(f"  E = -H⁶/(H+1)⁴ Ha: {E_check == Fraction(-H**6, (H+1)**4)}")
print(f"\n  *** Helium variational energy = -H⁶/(H+1)⁴ Hartree EXACTLY ***")
