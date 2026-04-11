#!/usr/bin/env python3
"""
Test whether the vacuum ladder's 3+1 structure satisfies the
self-consistency equation (H-1)² = H+1.

The four rungs: floor, depth, instanton, Λ
The three steps: α₁ = 2.16, β = 16.2, α₂ = 2.45

Hypothesis: the steps are a "meta mass function" with the same
structure as (s₁, s₂, s₃, θ).
"""

from mpmath import mp, mpf, sqrt, ln, exp, fabs, nstr, pi, log

mp.dps = 50
H = 3
K_star = mpf(7)/30
FLOOR = mpf(1)/H**3
S = mpf(810)/7

# The four values
Born_DS = mpf("0.00080948751214898223625")
lambda_0 = mpf("0.28291034660315")
lambda_1 = mpf("0.28131300001289")
instanton = (1-lambda_0)**(-mpf(1)/2) * (1-lambda_1)**(-mpf(1)/2) * exp(-S)
Lambda_obs = mpf("1.17e-123")

# Floor exponents
ln_floor = ln(FLOOR)
e_floor = mpf(1)
e_depth = ln(Born_DS) / ln_floor
e_inst = ln(instanton) / ln_floor
e_lambda = ln(Lambda_obs) / ln_floor

print("=" * 70)
print("THE LADDER AS A META MASS FUNCTION")
print("=" * 70)

print(f"\nFloor exponents (rungs): {nstr(e_floor,6)}, {nstr(e_depth,6)}, {nstr(e_inst,6)}, {nstr(e_lambda,6)}")

# Steps between rungs
s1 = e_depth / e_floor    # = e_depth
s2 = e_inst / e_depth     # the big step
s3 = e_lambda / e_inst    # the second small step

print(f"Steps: α₁={nstr(s1,8)}, β={nstr(s2,8)}, α₂={nstr(s3,8)}")
print(f"Sum of steps: {nstr(s1+s2+s3, 10)}")
print(f"Product of steps: {nstr(s1*s2*s3, 10)} = Lambda exponent")

# ============================================================
# THE 3+1 STRUCTURE
# ============================================================
print("\n" + "=" * 70)
print("3+1 STRUCTURE: Steps as (s₂, s₁, s₃, θ)")
print("=" * 70)

# In the mass function: s₁ >> s₂ ≈ s₃, θ = substrate
# In the ladder: β >> α₁ ≈ α₂, floor = substrate
# Map: β ↔ s₁, α₁ ↔ s₂, α₂ ↔ s₃, 1 ↔ θ

# Normalize like a mass function: total = α₁ + β + α₂ + 1
total = s1 + s2 + s3 + 1
print(f"\n  'Mass function' = (α₁, β, α₂, 1) normalized:")
print(f"  Total = {nstr(total, 10)}")
m_s1 = s2 / total  # β is the dominant section
m_s2 = s1 / total  # α₁ is the small section
m_s3 = s3 / total  # α₂ is the other small section
m_theta = 1 / total  # substrate

print(f"  s₁ (β/total) = {nstr(m_s1, 10)}")
print(f"  s₂ (α₁/total) = {nstr(m_s2, 10)}")
print(f"  s₃ (α₂/total) = {nstr(m_s3, 10)}")
print(f"  θ  (1/total) = {nstr(m_theta, 10)}")

# Check Born probability of this meta mass function
L2 = m_s1**2 + m_s2**2 + m_s3**2 + m_theta**2
meta_born = m_theta**2 / L2
print(f"\n  Born(meta) = θ²/||m||² = {nstr(meta_born, 15)}")
print(f"  1/H³ = {nstr(FLOOR, 15)}")
print(f"  Ratio = {nstr(meta_born/FLOOR, 10)}")

# ============================================================
# SELF-CONSISTENCY: (H-1)² = H+1 APPLIED TO STEPS
# ============================================================
print("\n" + "=" * 70)
print("SELF-CONSISTENCY EQUATION ON STEPS")
print("=" * 70)

# (H-1)² = H+1 for the mass function means the structure is self-referential.
# What's the analogous equation for the steps?

# The step ratio: β/(α geometric mean) should relate to H
alpha_geo = sqrt(s1 * s3)
print(f"\n  Geometric mean of small steps: √(α₁α₂) = {nstr(alpha_geo, 15)}")
print(f"  ln(H²+1) = ln(10) = {nstr(ln(mpf(H**2+1)), 15)}")
print(f"  ratio = {nstr(alpha_geo/ln(mpf(10)), 10)}")

# The big step
print(f"\n  β = {nstr(s2, 15)}")
print(f"  (H+1)² = {(H+1)**2}")
print(f"  ratio = {nstr(s2/(H+1)**2, 10)}")

# KEY TEST: does β = H × α₁ × α₂?
print(f"\n  H × α₁ × α₂ = {nstr(H * s1 * s3, 10)}")
print(f"  β = {nstr(s2, 10)}")
print(f"  ratio = {nstr(s2/(H*s1*s3), 10)}")

# Self-consistency: (something - 1)² = something + 1
# Applied to each step?
for name, val in [("α₁", s1), ("β", s2), ("α₂", s3)]:
    lhs = (val - 1)**2
    rhs = val + 1
    print(f"\n  ({name} - 1)² = {nstr(lhs, 8)},  {name} + 1 = {nstr(rhs, 8)},  ratio = {nstr(lhs/rhs, 8)}")

# What about: (α-1)² = α+1 where α is the geometric mean?
val = alpha_geo
lhs = (val - 1)**2
rhs = val + 1
print(f"\n  (α_geo - 1)² = {nstr(lhs, 10)}")
print(f"  α_geo + 1 = {nstr(rhs, 10)}")
print(f"  ratio = {nstr(lhs/rhs, 10)}")

# ============================================================
# THE 1 = [0, 1/27) + REST DECOMPOSITION
# ============================================================
print("\n" + "=" * 70)
print("UNIT INTERVAL DECOMPOSITION")
print("=" * 70)

# The Born probability ranges [0, 1].
# Floor divides at 1/27.
# Within [0, 1/27), three more values create 4 zones.

zones = [
    ("Λ zone [0, Λ)", Lambda_obs),
    ("Instanton zone [Λ, inst)", instanton - Lambda_obs),
    ("Depth zone [inst, depth)", Born_DS - instanton),
    ("Floor zone [depth, 1/27)", FLOOR - Born_DS),
]

print(f"\n  The nothing region [0, 1/27) = {nstr(FLOOR, 10)} decomposed:")
for name, width in zones:
    frac = width / FLOOR
    print(f"    {name:35s}: width ≈ {nstr(width, 5):>12s}  fraction = {nstr(frac, 5)}")

# The four zone FRACTIONS (of the nothing region):
f1 = (FLOOR - Born_DS) / FLOOR  # largest
f2 = Born_DS / FLOOR  # depth
f3 = instanton / FLOOR  # instanton
f4 = Lambda_obs / FLOOR  # Λ

print(f"\n  As fractions of 1/27:")
print(f"    f₁ (floor zone)    = {nstr(f1, 15)} ≈ 1 - Born_DS/floor")
print(f"    f₂ (depth zone)    = {nstr(f2, 15)} ≈ Born_DS/floor")
print(f"    f₃ (instanton)     = {nstr(f3, 10)}")
print(f"    f₄ (Λ)             = {nstr(f4, 10)}")
print(f"    sum = {nstr(f1+f2+f3+f4, 10)} (should = 1)")

# The Born floor of THESE fractions
L2_f = f1**2 + f2**2 + f3**2 + f4**2
born_f = f4**2 / L2_f
print(f"\n  'Born' of zone fractions = f₄²/||f||² = {nstr(born_f, 10)}")
print(f"  This is essentially (Λ/floor)² / (1-Born_DS/floor)² ≈ Λ²/floor² ≈ {nstr(Lambda_obs**2/FLOOR**2, 5)}")

# ============================================================
# THE FOUR-FOLD SELF-SIMILARITY
# ============================================================
print("\n" + "=" * 70)
print("SELF-SIMILAR STRUCTURE")
print("=" * 70)

# The mass function has (s₁, s₂, s₃, θ) with s₁ >> s₂ = s₃ >> θ
# and Born(m) = 1/27.
#
# Now consider: does the LADDER satisfy an analogous fixed-point equation?
# The mass function is a fixed point of Φ = floor ∘ DS.
# Is the ladder a fixed point of some "meta-Φ"?

# The instanton action S = H³/K* connects the floor to the instanton.
# If the same S connects the instanton to Λ:
# Λ = instanton × f(Born_DS)
# where f encodes the vacuum depth.

# The formula Λ = instanton × Born_DS^(47/2) says:
# The step from instanton to Λ (in log) = (47/2) × ln(Born_DS)
# = (47/2) × (-7.12) = -167.3

# And the step from floor to instanton (in log) = ln(instanton/floor)
# = ln(instanton) - ln(floor)
# = -115.38 - (-3.30) = -112.1

# So the two "big" log-steps are:
# floor → instanton: -112.1
# instanton → Λ: -167.3

# Ratio: 167.3/112.1 = 1.49 ≈ 3/2?
step_fi = ln(instanton) - ln(FLOOR)
step_il = ln(Lambda_obs) - ln(instanton)
print(f"\n  ln-step floor→instanton: {nstr(step_fi, 10)}")
print(f"  ln-step instanton→Λ:     {nstr(step_il, 10)}")
print(f"  ratio = {nstr(step_il/step_fi, 10)}")
print(f"  3/2 = {nstr(mpf(3)/2, 10)}")
print(f"  (H+1)/H = {nstr(mpf(H+1)/H, 10)}")

# ============================================================
# DEFINITIVE TEST: Λ from n = 47/2
# ============================================================
print("\n" + "=" * 70)
print("DEFINITIVE: Λ = instanton × Born_DS^(47/2)")
print("=" * 70)

n_W = mpf(47)/2
Lambda_predicted = instanton * Born_DS**n_W
log10_pred = ln(fabs(Lambda_predicted)) / ln(10)
log10_obs = ln(fabs(Lambda_obs)) / ln(10)

print(f"\n  n = 47/2 = {nstr(n_W, 6)}")
print(f"  Born_DS^(47/2) = {nstr(Born_DS**n_W, 10)}")
print(f"  Λ_predicted = {nstr(Lambda_predicted, 10)}")
print(f"  log10(Λ_pred) = {nstr(log10_pred, 10)}")
print(f"  log10(Λ_obs)  = {nstr(log10_obs, 10)}")
print(f"  Discrepancy: {nstr(fabs(log10_pred - log10_obs), 5)} dex")

# The 47 decomposition
print(f"\n  47 = H(H+1)² - 1 = 3×16 - 1 = {H*(H+1)**2 - 1}")
print(f"  The W multiplier: m_W/m(0++) = 47")
print(f"  47/2 = 23.5 = half the fermion-state count minus 1")

# ============================================================
# WHAT ABOUT Λ = floor^(product of all steps)?
# And the product = ?
# ============================================================
print("\n" + "=" * 70)
print("PRODUCT STRUCTURE")
print("=" * 70)

product = s1 * s2 * s3
print(f"\n  Product of steps = α₁ × β × α₂ = {nstr(product, 15)}")
print(f"  This IS the Lambda exponent (by construction)")

# But can we decompose it as H × α⁴ where α = ln(10)?
print(f"\n  H × ln(10)⁴ = {nstr(H * ln(mpf(10))**4, 10)}")
print(f"  H × α_geo⁴ = {nstr(H * alpha_geo**4, 10)}")
print(f"  product = {nstr(product, 10)}")
print(f"  ratio = {nstr(product/(H * alpha_geo**4), 10)}")

# If β = H × α₁ × α₂ (tested above), then:
# product = α₁ × (H × α₁ × α₂) × α₂ = H × α₁² × α₂²
# = H × (α₁α₂)²
# = H × (geo mean)⁴... wait
# product = H × (α₁α₂)² = H × (√(α₁α₂))⁴ = H × α_geo⁴? No:
# H × α₁² × α₂² ≠ H × (α₁α₂)². Wait:
# H × α₁² × α₂² = H × (α₁ × α₂)²  ← NO, that's H × α₁²α₂²
# and (α₁α₂)² = α₁²α₂². So yes: product = H × (α₁α₂)².
# But α₁α₂ = 2.16 × 2.45 = 5.30, so H × 5.30² = 3 × 28.07 = 84.2.
# Actual product = 85.88. Ratio: 85.88/84.2 = 1.02. Off by 2%.

print(f"\n  H × (α₁α₂)² = {nstr(H * (s1*s3)**2, 10)}")
print(f"  If β = H × α₁α₂ exactly, product = H × (α₁α₂)²")
print(f"  Discrepancy: {nstr(fabs(product - H*(s1*s3)**2)/product*100, 5)}%")

# ============================================================
# RATIO OF BIG STEPS = (H+1)/H = 4/3?
# ============================================================
print("\n" + "=" * 70)
print("RATIO TEST: step_il/step_fi = (H+1)/H?")
print("=" * 70)

# ln steps in natural log
# floor→instanton: ln(instanton) - ln(1/27) = -115.38 + 3.30 = -112.08
# instanton→Λ: ln(Λ) - ln(instanton) = -283.2 + 115.38 = -167.82

ratio_big_steps = step_il / step_fi
target_43 = mpf(H+1)/H

print(f"\n  step_il/step_fi = {nstr(ratio_big_steps, 15)}")
print(f"  (H+1)/H = 4/3 = {nstr(target_43, 15)}")
print(f"  err = {nstr(fabs(ratio_big_steps - target_43)/target_43*100, 5)}%")

# If exactly 4/3: instanton → Λ step is (4/3) × floor → instanton step
# ln(Λ) = ln(floor) + (1 + 4/3) × (ln(instanton) - ln(floor))
# = ln(floor) + (7/3) × (ln(inst) - ln(floor))
# = (1 - 7/3) ln(floor) + (7/3) ln(inst)
# = (-4/3) ln(floor) + (7/3) ln(inst)
# = (7 ln(inst) - 4 ln(floor)) / 3

Lambda_from_ratio = exp((7*ln(instanton) - 4*ln(FLOOR))/3)
log10_ratio_pred = ln(fabs(Lambda_from_ratio))/ln(10)
print(f"\n  If ratio = 4/3 exactly:")
print(f"  Λ = exp[(7·ln(inst) - 4·ln(floor))/3]")
print(f"  log10(Λ) = {nstr(log10_ratio_pred, 10)}")
print(f"  vs observed: {nstr(log10_obs, 10)}")
print(f"  discrepancy: {nstr(fabs(log10_ratio_pred - log10_obs), 5)} dex")

# ============================================================
# COMBINE ALL: best prediction
# ============================================================
print("\n" + "=" * 70)
print("BEST PREDICTIONS COMPARED")
print("=" * 70)

predictions = {
    "n=47/2": instanton * Born_DS**(mpf(47)/2),
    "n=24 (H(H²-1))": instanton * Born_DS**24,
    "n=23 (30(1-K*))": instanton * Born_DS**23,
    "ratio=4/3": Lambda_from_ratio,
    "α²×16": FLOOR**(alpha_geo**2 * 16),
    "H×αgeo⁴": FLOOR**(H * alpha_geo**4),
}

print(f"\n  {'Method':25s} {'log10(Λ)':>12s} {'err (dex)':>10s}")
print(f"  {'-'*25} {'-'*12} {'-'*10}")
for name, val in predictions.items():
    l10 = float(ln(fabs(val))/ln(10))
    err = abs(l10 - float(log10_obs))
    print(f"  {name:25s} {l10:12.2f} {err:10.2f}")

print(f"\n  {'Observed':25s} {float(log10_obs):12.2f}")


if __name__ == "__main__":
    pass  # All code runs at module level
