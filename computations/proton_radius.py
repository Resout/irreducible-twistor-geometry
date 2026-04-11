"""
Proton charge radius from the H=3 framework
============================================

Observation: r_p = (H+1) × ℏ/(m_p c) = 4 × Compton wavelength of proton
Compare to CODATA 2018 (muonic hydrogen): 0.8409 ± 0.0004 fm

This script investigates:
1. Numerical verification
2. Dipole form factor hypothesis
3. H-dependence (is H=3 special?)
4. Alternative scaling formulas
5. Honest assessment of what is derived vs observed
"""

import numpy as np

H = 3

# =============================================================================
# Physical constants
# =============================================================================
hbar_c = 197.3269804  # MeV·fm (CODATA 2018)
m_p = 938.272046      # MeV/c² (proton mass)
m_p_kg = 1.67262192e-27  # kg
hbar = 1.054571817e-34    # J·s
c = 2.99792458e8          # m/s

# Measured proton charge radius
r_p_exp = 0.8409      # fm (muonic hydrogen, CODATA 2018)
r_p_err = 0.0004      # fm

# Also: PRad (2019) electron scattering
r_p_prad = 0.831      # fm
r_p_prad_err = 0.012  # fm

# e-p scattering (pre-2010, larger value)
r_p_old = 0.8751      # fm
r_p_old_err = 0.0061  # fm

print("=" * 72)
print("PROTON CHARGE RADIUS FROM H=3 FRAMEWORK")
print("=" * 72)

# =============================================================================
# 1. Basic numerical verification
# =============================================================================
print("\n" + "=" * 72)
print("1. NUMERICAL VERIFICATION")
print("=" * 72)

lambda_C = hbar_c / m_p  # proton Compton wavelength in fm
print(f"\nProton Compton wavelength: ℏ/(m_p c) = {lambda_C:.6f} fm")
print(f"  = {hbar_c:.4f} / {m_p:.3f} MeV·fm / MeV")

r_p_pred = (H + 1) * lambda_C
print(f"\nPrediction: r_p = (H+1) × ℏ/(m_p c) = {H+1} × {lambda_C:.6f}")
print(f"  r_p(pred) = {r_p_pred:.6f} fm")
print(f"  r_p(exp)  = {r_p_exp:.4f} ± {r_p_err:.4f} fm  [muonic H]")

dev = r_p_pred - r_p_exp
dev_pct = 100 * dev / r_p_exp
dev_sigma = dev / r_p_err
print(f"\n  Deviation: {dev:+.6f} fm  ({dev_pct:+.4f}%)")
print(f"  In units of experimental error: {dev_sigma:+.2f} σ")

print(f"\n  vs PRad (2019):  {r_p_prad:.3f} ± {r_p_prad_err:.3f} fm")
dev_prad = r_p_pred - r_p_prad
print(f"    Deviation: {dev_prad:+.4f} fm  ({100*dev_prad/r_p_prad:+.2f}%)")
print(f"    {dev_prad/r_p_prad_err:+.2f} σ")

print(f"\n  vs old e-p:      {r_p_old:.4f} ± {r_p_old_err:.4f} fm")
dev_old = r_p_pred - r_p_old
print(f"    Deviation: {dev_old:+.4f} fm  ({100*dev_old/r_p_old:+.2f}%)")
print(f"    {dev_old/r_p_old_err:+.2f} σ")

# =============================================================================
# 2. The ratio r_p / lambda_C
# =============================================================================
print("\n" + "=" * 72)
print("2. THE RATIO r_p / λ_C")
print("=" * 72)

ratio_exp = r_p_exp / lambda_C
print(f"\n  r_p(exp) / λ_C = {r_p_exp} / {lambda_C:.6f} = {ratio_exp:.4f}")
print(f"  H + 1 = {H + 1}")
print(f"  Difference from integer: {ratio_exp - (H+1):.4f}")
print(f"  Fractional deviation from {H+1}: {(ratio_exp - (H+1))/(H+1):.6f}")

# =============================================================================
# 3. Dipole form factor analysis
# =============================================================================
print("\n" + "=" * 72)
print("3. DIPOLE FORM FACTOR ANALYSIS")
print("=" * 72)

print("""
The standard dipole form factor for the proton is:
  G_E(q²) = 1 / (1 + q²/Λ²)²

The power 2 and the scale Λ² = 0.71 GeV² are empirical fits.

Framework hypothesis: the form factor arises from the (H+1)-dimensional
mass space structure. Consider:
  G_E(q²) = [1 + q² r₀² / n]^(-n)    with r₀ = ℏ/(m_p c)

For n = (H+1)/2 = 2, this gives the standard dipole. The charge radius is:
  r² = -6 dG_E/dq²|₀ = 6 n r₀² / n = 6 r₀²  ... no, let's be careful.
""")

# For G_E = (1 + q²a²/n)^(-n), we have:
# dG_E/dq² = -a² (1 + q²a²/n)^(-n-1)
# r² = -6 × dG_E/dq²|₀ = 6a²
# So r = √6 × a, and we need a = r/√6

# Standard dipole: G_E = 1/(1 + q²/Λ²)²
# r² = 12/Λ²
# Λ² = 12/r² = 12/(0.8409 fm)² = 12/(0.8409 × 1/197.327 fm×MeV/MeV)²

# Convert r_p to natural units (1/GeV)
r_p_inv_GeV = r_p_exp / (hbar_c * 1e-3)  # fm / (GeV·fm) = 1/GeV
print(f"  r_p in natural units: {r_p_inv_GeV:.4f} GeV⁻¹")

Lambda2_from_r = 12.0 / r_p_inv_GeV**2
print(f"\n  Standard dipole: Λ² = 12/r² = {Lambda2_from_r:.4f} GeV²")
print(f"  Empirical value: Λ² = 0.71 GeV²")
print(f"  Agreement: {100*(Lambda2_from_r - 0.71)/0.71:+.1f}%")

# Framework prediction for Λ²
Lambda2_fw = 12.0 / (r_p_pred / (hbar_c * 1e-3))**2
print(f"\n  Framework dipole: Λ²(fw) = 12/r(fw)² = {Lambda2_fw:.4f} GeV²")

# Now: can we derive Λ² from framework quantities?
# Λ² = 12/(r_p²) = 12/((H+1)² λ_C²) = 12 m_p² / (H+1)²
# In GeV: Λ² = 12 × (0.938272)² / 16
Lambda2_derived = 12 * (m_p * 1e-3)**2 / (H + 1)**2
print(f"\n  Derived: Λ² = 12 m_p² / (H+1)² = {Lambda2_derived:.4f} GeV²")
print(f"  = 3 m_p² / (H+1)² × 4 = 3 m_p² / 4")
print(f"  = {3 * (m_p*1e-3)**2 / 4:.4f} GeV²")

# The 3 in 12 = 4×3 could be H itself
print(f"\n  Note: 12 = 4 × 3 = (H+1) × H")
print(f"  So Λ² = H(H+1) m_p² / (H+1)² = H m_p² / (H+1)")
print(f"  = {H * (m_p*1e-3)**2 / (H+1):.4f} GeV²")
Lambda2_clean = H * (m_p * 1e-3)**2 / (H + 1)
print(f"  vs empirical 0.71 → deviation: {100*(Lambda2_clean - 0.71)/0.71:+.1f}%")

# The dipole power
print(f"\n  Dipole power = 2 = (H+1)/2 = {(H+1)/2}")
print(f"  This is exact: the standard dipole IS the (H+1)/2 power.")

# =============================================================================
# 4. H-scan: is H=3 special?
# =============================================================================
print("\n" + "=" * 72)
print("4. H-DEPENDENCE: IS H=3 SPECIAL?")
print("=" * 72)

print(f"\n  Formula: r_p = (H+1) × ℏ/(m_p c) = (H+1) × {lambda_C:.6f} fm")
print(f"  Measured: {r_p_exp:.4f} fm\n")
print(f"  {'H':>4}  {'H+1':>5}  {'r_pred (fm)':>12}  {'dev (%)':>10}  {'dev (σ)':>10}")
print(f"  {'-'*4}  {'-'*5}  {'-'*12}  {'-'*10}  {'-'*10}")

for h in range(1, 9):
    r_h = (h + 1) * lambda_C
    d_pct = 100 * (r_h - r_p_exp) / r_p_exp
    d_sig = (r_h - r_p_exp) / r_p_err
    marker = " <-- H=3" if h == 3 else ""
    print(f"  {h:4d}  {h+1:5d}  {r_h:12.6f}  {d_pct:+10.4f}  {d_sig:+10.2f}{marker}")

# =============================================================================
# 5. Alternative scaling formulas
# =============================================================================
print("\n" + "=" * 72)
print("5. ALTERNATIVE SCALING FORMULAS")
print("=" * 72)

alternatives = [
    ("H+1",                  H + 1),
    ("H",                    H),
    ("H+2",                  H + 2),
    ("√(H+1)",               np.sqrt(H + 1)),
    ("√(H²+1)",              np.sqrt(H**2 + 1)),
    ("(H+1)/√H",             (H + 1) / np.sqrt(H)),
    ("2H/(H-1)",             2 * H / (H - 1)),
    ("(H²-1)/H",             (H**2 - 1) / H),
    ("H+1/H",                H + 1.0/H),
    ("π",                    np.pi),
    ("4",                    4),
    ("2√(H+1/H)",            2 * np.sqrt((H + 1.0) / H)),
    ("(H+1)²/H²",            (H + 1)**2 / H**2),
    ("H·(H-1)+2",            H * (H - 1) + 2),
    ("√(H(H+1)(H+2)/2)",     np.sqrt(H * (H + 1) * (H + 2) / 2)),
]

print(f"\n  {'Formula':>25}  {'Value':>8}  {'r_pred (fm)':>12}  {'dev (%)':>10}  {'dev (σ)':>10}")
print(f"  {'-'*25}  {'-'*8}  {'-'*12}  {'-'*10}  {'-'*10}")

for name, val in sorted(alternatives, key=lambda x: abs(x[1] * lambda_C - r_p_exp)):
    r_alt = val * lambda_C
    d_pct = 100 * (r_alt - r_p_exp) / r_p_exp
    d_sig = (r_alt - r_p_exp) / r_p_err
    print(f"  {name:>25}  {val:8.4f}  {r_alt:12.6f}  {d_pct:+10.4f}  {d_sig:+10.2f}")

# =============================================================================
# 6. Physical interpretation attempts
# =============================================================================
print("\n" + "=" * 72)
print("6. PHYSICAL INTERPRETATION")
print("=" * 72)

print("""
Attempt A: Constituent counting
  The proton has 3 valence quarks + 1 gluonic substrate = H+1 = 4 constituents.
  Each contributes one Compton wavelength to the spatial extent.
  r_p = (H+1) × λ_C.

  STATUS: Numerically correct. But WHY does each add exactly one λ_C?
  This is a restatement, not a derivation.

Attempt B: Dipole form factor from mass-space dimension
  The form factor G_E = (1 + q²/Λ²)^(-2) has power 2 = (H+1)/2.
  This suggests a (H+1)-dimensional isotropic Gaussian in mass space,
  projected to physical space (dividing by 2 for complex→real).

  If each complex dimension contributes a factor (1+q²a²)^(-1) to G_E,
  then (H+1)/2 = 2 complex dimensions ≡ C² = S (the spinor space).

  The charge radius then becomes:
    r² = 6a² × (H+1)/2 = 3a²(H+1)
  With a² = λ_C² (each dimension's scale is the Compton wavelength):
    r² = 3(H+1)λ_C² = 12λ_C²
    r = 2√3 × λ_C = 0.3642 fm  ... NO, this doesn't work.

  Actually for G_E = (1+q²/Λ²)^(-n), we get r² = 6n/Λ².
  With n=2 and Λ² = m_p²(H+1)/H:
    r² = 12H / [m_p²(H+1)]
       = 12 × 3 / [m_p² × 4] × ℏ²c²
       = 9ℏ²/(m_p²c²)
    r  = 3 × ℏ/(m_p c) = 3λ_C = 0.6309 fm ... NO.

  STATUS: The dipole form factor power 2=(H+1)/2 is suggestive but
  does not independently produce the factor H+1 in the radius.

Attempt C: From the transfer matrix spectrum
  The proton mass satisfies m_p = const × Λ_QCD.
  The charge radius satisfies r_p = const / Λ_QCD.
  So r_p × m_p = dimensionless constant × ℏc.

  In the framework: r_p × m_p c / ℏ = H + 1 = 4.
  This is a DIMENSIONLESS prediction: the product r_p × m_p (in natural
  units) equals H+1.

  The proton is a bound state of the H=3 system. The number H+1 = 4
  counts the real dimension of the base space B = S³ on which the
  transfer operator acts. The proton "fills" the full base space.

  STATUS: This is the most promising route. dim_R(S³) = H+1 = 4 is
  a geometric fact about the framework. The claim is:
    r_p × m_p = dim_R(B) × ℏ/c = (H+1) × ℏ/c
""")

# =============================================================================
# 7. The dimensionless product r_p × m_p
# =============================================================================
print("=" * 72)
print("7. THE DIMENSIONLESS PRODUCT r_p × m_p c / ℏ")
print("=" * 72)

product_exp = r_p_exp * m_p / hbar_c
product_pred = H + 1
print(f"\n  r_p × m_p c / ℏ (measured) = {product_exp:.4f}")
print(f"  H + 1                      = {product_pred}")
print(f"  Deviation: {product_exp - product_pred:+.4f}  ({100*(product_exp - product_pred)/product_pred:+.4f}%)")

print(f"\n  For comparison, other particle radii × mass products:")
# neutron
m_n = 939.565         # MeV
r_n2 = -0.1161        # fm² (negative! = charge on outside)
print(f"  Neutron:  r²_n = {r_n2} fm² (negative, different structure)")

# pion
m_pi = 139.571        # MeV
r_pi = 0.659          # fm ± 0.004
prod_pi = r_pi * m_pi / hbar_c
print(f"  Pion:     r_π × m_π c/ℏ = {prod_pi:.4f}  (= {prod_pi:.3f}, cf. H-1={H-1})")

# kaon
m_K = 493.677         # MeV
r_K = 0.560           # fm ± 0.031
prod_K = r_K * m_K / hbar_c
print(f"  Kaon:     r_K × m_K c/ℏ = {prod_K:.4f}  (= {prod_K:.3f})")

# =============================================================================
# 8. Connection to S³ = SU(2) geometry
# =============================================================================
print("\n" + "=" * 72)
print("8. CONNECTION TO S³ = SU(2) BASE SPACE")
print("=" * 72)

print(f"""
  The base space B = S³ has:
    - Real dimension: dim_R(S³) = H+1 = 4  ... wait, dim(S³)=3, not 4!

  CORRECTION: S³ is a 3-dimensional manifold embedded in R⁴.
    - dim_R(S³) = 3 = H
    - dim_R(R⁴) = 4 = H+1  (the embedding space)

  So the factor H+1 corresponds to the EMBEDDING dimension of B=S³,
  i.e., the ambient space R^(H+1) in which S^H lives.

  Alternative: H+1 = dim_R(C²) = dim_R(S) where S is the spinor space.
  The proton is a spin-1/2 particle. Its charge distribution extends
  over the full real extent of S = C², which has real dimension H+1 = 4.

  This is actually clean:
    r_p = dim_R(S) × ℏ/(m_p c) = (H+1) × λ_C

  where S = C^2 is the fundamental spinor space of the framework.
""")

# =============================================================================
# 9. Sensitivity and robustness
# =============================================================================
print("=" * 72)
print("9. SENSITIVITY ANALYSIS")
print("=" * 72)

print(f"\n  Using different experimental values of r_p:\n")
measurements = [
    ("Muonic H (2010)",        0.84184, 0.00067),
    ("Muonic H (2013)",        0.84087, 0.00039),
    ("CODATA 2018",            0.8409,  0.0004),
    ("Bezginov et al (2019)",  0.833,   0.010),
    ("PRad (2019)",            0.831,   0.012),
    ("Xiong et al (2019)",     0.831,   0.014),
    ("Muonic D → r_p (2016)",  0.8356,  0.0020),
]

print(f"  {'Measurement':>30}  {'r_p (fm)':>10}  {'err':>8}  {'ratio':>8}  {'dev σ':>8}")
print(f"  {'-'*30}  {'-'*10}  {'-'*8}  {'-'*8}  {'-'*8}")
for name, r, err in measurements:
    ratio = r / lambda_C
    dev_s = (r_p_pred - r) / err
    print(f"  {name:>30}  {r:10.5f}  {err:8.5f}  {ratio:8.4f}  {dev_s:+8.2f}")

print(f"\n  Framework prediction: r_p = {r_p_pred:.6f} fm")
print(f"  All muonic measurements cluster near ratio = {r_p_exp/lambda_C:.3f} ≈ H+1 = 4")

# =============================================================================
# 10. What about the neutron?
# =============================================================================
print("\n" + "=" * 72)
print("10. NEUTRON CHARGE RADIUS")
print("=" * 72)

r_n_sq = -0.1161  # fm² (mean square charge radius)
r_n_sq_err = 0.0022
print(f"\n  Neutron ⟨r²⟩_n = {r_n_sq} ± {r_n_sq_err} fm²")
print(f"  (Negative because the negative charge is farther out than positive)")

# If proton r_p² = (H+1)² λ_C², what about neutron?
r_p_sq_pred = r_p_pred**2
lambda_C_n = hbar_c / m_n
print(f"\n  Proton ⟨r²⟩_p (pred) = {r_p_sq_pred:.4f} fm²")
print(f"  Proton ⟨r²⟩_p (exp)  = {r_p_exp**2:.4f} fm²")
print(f"  Neutron λ_C = {lambda_C_n:.6f} fm")

# Foldy term for neutron
mu_n = -1.91304  # nuclear magnetons
foldy = -3 * mu_n / (2 * m_n**2) * hbar_c**2  # in fm²
# Actually Foldy term = -3μ_n/(2m_n²) in natural units
# = -3 × (-1.913) × ℏ²/(2 m_n² c²) ... let me be careful
# Foldy = -3κ_n/(4m_n²) where κ_n = μ_n - 0 = μ_n (anomalous moment)
kappa_n = -1.91304  # anomalous magnetic moment in nuclear magnetons
foldy_term = -kappa_n / (4 * m_n**2) * hbar_c**2 * 6  # Foldy: 3κ/(2M²)
# Actually r²_Foldy = -3κ/(2M²) in units where ℏ=c=1
# With dimensions: -3κ ℏ²/(2 M² c²)
foldy_proper = -3 * kappa_n * hbar_c**2 / (2 * m_n**2)
print(f"\n  Foldy term (neutron): {foldy_proper:.4f} fm²")
print(f"  This accounts for about {100*foldy_proper/r_n_sq:.0f}% of ⟨r²⟩_n")

# =============================================================================
# 11. Summary assessment
# =============================================================================
print("\n" + "=" * 72)
print("11. SUMMARY AND HONEST ASSESSMENT")
print("=" * 72)

print(f"""
  OBSERVATION (robust):
    r_p = (H+1) × ℏ/(m_p c) = 4 × λ_C

    Numerically: {r_p_pred:.6f} fm vs {r_p_exp:.4f} ± {r_p_err:.4f} fm
    Deviation: {dev_pct:+.4f}% = {dev_sigma:+.2f}σ

    The dimensionless product r_p × m_p c/ℏ = {product_exp:.4f} ≈ H+1 = 4

  STRONGEST INTERPRETATION:
    H+1 = 4 = dim_R(C²) = dim_R(S), the real dimension of spinor space.
    The proton's charge distribution extends over the full real extent
    of the fundamental representation space S = C².

    This is a geometric statement: the charge radius measures how many
    real dimensions the proton wave function explores, times the
    natural length scale λ_C = ℏ/(m_p c).

  WHAT IS DERIVED vs WHAT IS IDENTIFIED:
    DERIVED:   S = C² is the fundamental space (from H=3 first principles)
               dim_R(S) = 2(H-1) ... wait, dim_R(C²) = 4 = 2×2
               Actually dim_R(C^n) = 2n, and S = C² → dim_R = 4 = H+1.
               So dim_R(S) = H+1 is a THEOREM (since dim_C(S)=2, dim_R=4=H+1).

    IDENTIFIED: The charge radius equals dim_R(S) × λ_C.
                WHY the charge radius should equal this particular product
                is NOT derived from first principles. It is a numerical
                observation that r_p/λ_C = 3.998 ≈ 4 = H+1 = dim_R(S).

    PARTIAL:   The dipole form factor power 2 = dim_C(S) = (H+1)/2 is
               suggestive but does not independently fix the radius.

  SIGNIFICANCE:
    - The relation r_p × m_p = (H+1)ℏ/c is accurate to 0.04%
    - It involves only the framework's fundamental integer H=3
    - H+1 has clear geometric meaning (real dim of spinor space)
    - But the derivation gap remains: WHY does the charge radius
      equal dim_R(S) Compton wavelengths?

  COMPARISON TO OTHER RELATIONS:
    - Pion:  r_π × m_π c/ℏ = {prod_pi:.3f} (not a small integer)
    - Kaon:  r_K × m_K c/ℏ = {prod_K:.3f} (not a small integer)
    - Only the proton gives a ratio close to an integer, and that
      integer is H+1. This makes it less likely to be coincidence.
""")

# =============================================================================
# 12. Quick check: could K* = 7/30 appear as a correction?
# =============================================================================
print("=" * 72)
print("12. K* = 7/30 AS A CORRECTION?")
print("=" * 72)

K_star = 7.0 / 30.0

# r_p = (H+1)(1 - δ) λ_C, what is δ?
delta = 1 - r_p_exp / r_p_pred
print(f"\n  If r_p = (H+1)(1-δ)λ_C, then δ = {delta:.6f}")
print(f"  K*/H² = {K_star/H**2:.6f}")
print(f"  K*²   = {K_star**2:.6f}")
print(f"  K*/2π = {K_star/(2*np.pi):.6f}")
print(f"  1/(H+1)³ = {1/(H+1)**3:.6f}")

# Check if deviation could be K*-related
print(f"\n  The deviation is only {abs(dev_sigma):.1f}σ, consistent with zero.")
print(f"  Any K*-correction would be speculation on top of an already")
print(f"  marginal deviation. The formula r_p = (H+1)λ_C is sufficient.")

print("\n" + "=" * 72)
print("DONE")
print("=" * 72)
