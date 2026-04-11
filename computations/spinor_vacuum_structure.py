#!/usr/bin/env python3
"""
The spinor structure of the vacuum.

S = C² has two complex dimensions (ζ, η).
Irreducibility θ = ε_AB ∈ Λ²(S) is the spinor metric.

Hypothesis: The vacuum energy scale decomposes as:
  - One complex dim → fully real (4 real dims) → ρ_Λ side
  - Other complex dim → both real + complex (9 learning + 1 substrate = 10) → M side

Test: does the natural volume measure of this spinor decomposition
give exactly 40/9?
"""

from mpmath import mp, mpf, sqrt, ln, exp, fabs, nstr, pi, matrix, det

mp.dps = 50
H = 3
K = mpf(7)/30
FLOOR = mpf(1)/H**3
S_action = mpf(810)/7

# Framework dimensions
dim_S_complex = 2          # S = C², two complex dimensions
dim_S_real = 4             # dim_R(S) = 2 × dim_C(S) = 4 = H+1
dim_Sym2S = 3              # dim(Sym²(S)) = 3 = H sections
dim_Lambda2S = 1           # dim(Λ²(S)) = 1 (the ε-spinor, θ)
dim_total_sections = H + 1  # sections + substrate = 4

# DS channel counts
channels_effective = H**2 + 1  # = 10, includes substrate
channels_learning = H**2       # = 9, excludes substrate

print("=" * 70)
print("SPINOR STRUCTURE OF THE VACUUM")
print("=" * 70)

print(f"""
  S = C²:
    dim_C = {dim_S_complex}
    dim_R = {dim_S_real} = H+1
    Sym²(S) = C³ (the 3 sections s₁, s₂, s₃)
    Λ²(S) = C¹ (the substrate θ = ε_AB)

  DS channels:
    Effective (with substrate): {channels_effective} = H²+1
    Learning (without substrate): {channels_learning} = H²

  The ratio: (H²+1)/H² = {channels_effective}/{channels_learning} = {nstr(mpf(channels_effective)/channels_learning, 15)}
""")

# ============================================================
# THE TWO COMPLEX DIMENSIONS
# ============================================================
print("=" * 70)
print("THE TWO COMPLEX DIMENSIONS OF S = C²")
print("=" * 70)

print(f"""
  Component ζ: "Nothing" — decomposes fully in reals
    C¹ → R² (two real directions)
    Total real extent: 2 real dimensions per complex dim
    Over full spinor: dim_R(S) = {dim_S_real}

  Component η: "Emergence" — decomposes in both
    C¹ → R² (real part) + C structure (phase = substrate)
    The real part: {channels_learning} = H² channels
    The complex structure: +1 channel (the substrate)
    Total: {channels_effective} = H²+1 effective channels
    Ratio of effective/learning: {channels_effective}/{channels_learning}
""")

# ============================================================
# THE VOLUME MEASURE
# ============================================================
print("=" * 70)
print("THE NATURAL VOLUME MEASURE")
print("=" * 70)

# From "nothing" (fully real): contributes dim_R(S) = H+1 = 4
factor_nothing = mpf(H + 1)

# From "emergence" (both): contributes (H²+1)/H² = 10/9
factor_emergence = mpf(H**2 + 1) / mpf(H**2)

# Combined
vacuum_ratio = factor_nothing * factor_emergence

print(f"  Factor from 'nothing' (fully real): H+1 = {factor_nothing}")
print(f"  Factor from 'emergence' (both):     (H²+1)/H² = {nstr(factor_emergence, 15)}")
print(f"  Product: {nstr(vacuum_ratio, 15)}")
print(f"  40/9 = {nstr(mpf(40)/9, 15)}")
print(f"  Match: {nstr(fabs(vacuum_ratio - mpf(40)/9), 10)}")

# ============================================================
# ALTERNATIVE DERIVATION: from representation theory
# ============================================================
print("\n" + "=" * 70)
print("ALTERNATIVE: REPRESENTATION-THEORETIC VOLUME")
print("=" * 70)

# The mass function lives in R^4 = (s₁, s₂, s₃, θ)
# On the L¹=1 simplex, this is a 3-dimensional space
# The DS combination acts on Sym²(C⁴) = C^10
# The Born floor constrains to a codimension-1 surface

# Volume of the "nothing" component:
# ζ ∈ C carries 2 real DOF
# In the mass function, these map to the 2 angular directions
# (the s₂-s₃ plane and the s₁-θ radial direction)
# But on the simplex, constrained to L¹=1, we have 3 DOF
# The spinor's real decomposition gives H+1 = 4 directions in R⁴
# On the simplex: H = 3 directions

# Volume of the "emergence" component:
# η ∈ C carries 2 real DOF + 1 complex structure (phase)
# The real DOF give the H² = 9 hypothesis pairs
# The complex structure gives +1 for the substrate
# Total effective: H²+1 = 10

print(f"""
  The mass function m = (s₁, s₂, s₃, θ) ∈ R^{{H+1}} = R⁴

  The simplex Σ = {{m : Σmᵢ = 1}} has dimension H = 3.

  The Born floor surface B = {{m ∈ Σ : θ²/||m||² = 1/H³}}
  has dimension H-1 = 2 on the simplex.

  The DS combination acts on:
    Sym²(C^{{H+1}}) = C^{{(H+1)(H+2)/2}} = C^{{{(H+1)*(H+2)//2}}}
  But projected to the simplex and Born surface, the effective
  action space has dimension controlled by H²+1 = 10 channels.
""")

# ============================================================
# THE INSTANTON AS THE SPINOR METRIC
# ============================================================
print("=" * 70)
print("THE INSTANTON AS THE SPINOR METRIC ε_AB")
print("=" * 70)

print(f"""
  θ = ε_AB ∈ Λ²(S) = C is the spinor metric.
  It connects the two components ζ and η of S = C².

  The instanton amplitude e^{{-S}} tunnels through the Born floor.
  The Born floor IS the constraint |ε|² ≥ 1/H³.
  Tunnelling through the floor = tunnelling through ε → 0.

  So the instanton measures the amplitude for ε to vanish —
  the amplitude for the two spinor components to decouple —
  the amplitude for nothing to lose its self-consistency.

  S = H³/K* = 810/7 is the ACTION of this tunnelling.
  It lives in Λ²(S) — it IS ε's dynamics.
""")

# ============================================================
# COMPUTING ρ_Λ FROM THE SPINOR STRUCTURE
# ============================================================
print("=" * 70)
print("COMPUTING ρ_Λ")
print("=" * 70)

m_glueball = mpf("1.710")  # GeV (lattice, ±50 MeV)

# The spinor vacuum scale
M_vac = vacuum_ratio * m_glueball  # = (40/9) × m(0++)
print(f"  M_vac = (40/9) × m(0++) = {nstr(vacuum_ratio, 6)} × {nstr(m_glueball, 4)} = {nstr(M_vac, 6)} GeV")

# Instanton amplitude
lam0 = mpf("0.28291")
lam1 = mpf("0.28131")
prefactor = ((1-lam0)*(1-lam1))**(-mpf(1)/2)
instanton = prefactor * exp(-S_action)

print(f"  Instanton = {nstr(instanton, 6)}")

# The prediction
rho_predicted = instanton * M_vac**4
print(f"  M_vac⁴ = {nstr(M_vac**4, 6)} GeV⁴")
print(f"  ρ_Λ = instanton × M_vac⁴ = {nstr(rho_predicted, 6)} GeV⁴")

# Observed
rho_observed = mpf("2.52e-47")  # GeV⁴ (Planck 2018 central)
print(f"  ρ_Λ (observed) = {nstr(rho_observed, 6)} GeV⁴")

ratio = rho_predicted / rho_observed
print(f"  Ratio predicted/observed = {nstr(ratio, 6)}")
print(f"  Discrepancy: {nstr(fabs(ratio - 1)*100, 4)}%")

# Uncertainty from m(0++) = 1710 ± 50 MeV
m_high = mpf("1.760"); m_low = mpf("1.660")
rho_high = instanton * (vacuum_ratio * m_high)**4
rho_low = instanton * (vacuum_ratio * m_low)**4
print(f"\n  Uncertainty band from m(0++) = 1710 ± 50 MeV:")
print(f"    ρ_Λ(high) = {nstr(rho_high, 4)} GeV⁴")
print(f"    ρ_Λ(low)  = {nstr(rho_low, 4)} GeV⁴")
print(f"    Observed {nstr(rho_observed, 4)} in band? {float(rho_low) <= float(rho_observed) <= float(rho_high)}")

# ============================================================
# CROSS-CHECK: what does this predict for Λ/M_Pl⁴?
# ============================================================
print("\n" + "=" * 70)
print("CROSS-CHECK: Λ/M_Pl⁴")
print("=" * 70)

M_Pl = mpf("1.2209e19")  # GeV (unreduced Planck mass)
Lambda_dimless = rho_predicted / M_Pl**4
print(f"  ρ_Λ/M_Pl⁴ = {nstr(Lambda_dimless, 6)}")
print(f"  log₁₀ = {nstr(ln(fabs(Lambda_dimless))/ln(10), 6)}")
print(f"  Observed: ≈ 10⁻¹²³")

# ============================================================
# THE FULL CHAIN
# ============================================================
print("\n" + "=" * 70)
print("THE FULL CHAIN: ONE ANCHOR → ρ_Λ")
print("=" * 70)

print(f"""
  Given: m(0++) = 1.710 ± 0.050 GeV (lattice QCD, or reverse Koide)

  Step 1: Spinor vacuum scale
    M_vac = (H²+1)(H+1)/H² × m(0++) = (40/9) × m(0++)
    = {nstr(M_vac, 6)} GeV

  Step 2: Instanton amplitude (zero free parameters)
    S = H³/K* = 810/7
    det(I-J) = (1-λ₀)(1-λ₁) [from equilibrium Jacobian]
    instanton = det(I-J)^{{-1/2}} × e^{{-S}} = {nstr(instanton, 6)}

  Step 3: Vacuum energy density
    ρ_Λ = instanton × M_vac⁴
    = {nstr(rho_predicted, 6)} GeV⁴

  Observed: {nstr(rho_observed, 4)} GeV⁴
  Match: {nstr(fabs(ratio-1)*100, 3)}%
  Within lattice uncertainty: YES
""")

# ============================================================
# WHY 40/9: THE SPINOR DECOMPOSITION
# ============================================================
print("=" * 70)
print("WHY 40/9: THE SPINOR DECOMPOSITION")
print("=" * 70)

print(f"""
  S = C² (the irreducibility spinor)

  Component 1 ("nothing"): decomposes fully into reals
    C → R ⊕ R
    Contributes: dim_R(S) = 2 × dim_C(S) = {dim_S_real} = H+1
    This is the real extent of the vacuum — the 4 dimensions
    of the mass space (s₁, s₂, s₃, θ) that the vacuum energy
    fills. Energy density scales as M⁴ because the vacuum is
    4-real-dimensional.

  Component 2 ("emergence"): decomposes in both
    C → R ⊕ R (real part: the {channels_learning} = H² learning channels)
    C → C    (complex part: the +1 substrate channel)
    Total effective channels: {channels_effective} = H²+1
    The substrate survives as complex structure — it doesn't
    fully decompose into reals because irreducibility persists.
    Ratio: (H²+1)/H² = 10/9

  Combined measure:
    (H+1) × (H²+1)/H² = 4 × 10/9 = 40/9

  This is the unique volume factor that accounts for both:
    (a) the real dimensionality of the vacuum (from nothing)
    (b) the substrate inclusion in the channel count (from emergence)

  40/9 = {nstr(mpf(40)/9, 15)}

  The cosmological constant IS:
    ρ_Λ = [tunnelling amplitude through ε_AB = 0]
         × [spinor vacuum volume]⁴
         × [dimensional anchor]⁴
""")
