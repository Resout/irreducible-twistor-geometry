#!/usr/bin/env python3
"""
Spinor projection residual investigation.

The charged lepton lives on C^1 (the η component of the spinor (ζ,η) ∈ S = C²).
The Koide formula operates on sections (degree 2): s_i = ζ², ζη, η².
The mass formula m_k ∝ (1 + √2 cos(θ + 2πk/3))² squares the spinor amplitude.

Question: does projecting from the spinor (degree 1, C^1) to the section
(degree 2, Sym²(C²)) introduce the observed 9.83 × 10⁻⁶ residual?

Key facts:
- On C^1 (H_eff = 1): Born rule is trivial (p = |η|²/|η|² = 1)
- On C² (H_eff = 2): Gleason's theorem does NOT apply (H ≥ 3 required)
- On C⁴ (H = 3): Gleason applies, Born rule is unique
- The Born floor operates on C⁴ but is transparent to the angular mode (s₂-s₃)
- The angular mode IS the direction that separates ζ from η

Investigation: what happens when we work at the spinor level (degree 1)
instead of the section level (degree 2)?
"""

from mpmath import mp, mpf, cos, sin, pi, sqrt, log, fabs, nstr

mp.dps = 50

H = mpf(3)
theta_K = mpf(2) / mpf(9)  # Koide angle, proved exact

print("=" * 70)
print("SPINOR vs SECTION: DEGREE-1 vs DEGREE-2 MASS FORMULA")
print("=" * 70)
print()

# ============================================================
# SECTION-LEVEL (degree 2): standard Koide
# ============================================================

print("SECTION LEVEL (degree 2): √m_k = M × f_k")
print()

def f_section(k):
    """Section-level Koide factor (degree 1 in √m)"""
    return 1 + sqrt(2) * cos(theta_K + 2*pi*k/3)

def x_section(k):
    """Section-level mass factor (degree 2 in m)"""
    return f_section(k)**2

# k=0: tau, k=1: electron, k=2: muon
f_tau = f_section(0)
f_e = f_section(1)
f_mu = f_section(2)

x_tau = x_section(0)
x_e = x_section(1)
x_mu = x_section(2)

print(f"  f_e  = {nstr(f_e, 20)}   (spinor amplitude)")
print(f"  f_μ  = {nstr(f_mu, 20)}")
print(f"  f_τ  = {nstr(f_tau, 20)}")
print(f"  Σf_k = {nstr(f_e + f_mu + f_tau, 20)} = H = {H}")
print()
print(f"  x_e  = f_e²  = {nstr(x_e, 20)}")
print(f"  x_μ  = f_μ²  = {nstr(x_mu, 20)}")
print(f"  x_τ  = f_τ²  = {nstr(x_tau, 20)}")
print(f"  Σx_k = {nstr(x_e + x_mu + x_tau, 20)} = 2H = {2*H}")
print()

ratio_section = x_e / x_mu
print(f"  (m_e/m_μ)_section = x_e/x_μ = {nstr(ratio_section, 20)}")

# Experimental
m_e_exp = mpf('0.51099895069')
m_mu_exp = mpf('105.6583755')
m_tau_exp = mpf('1776.86')
ratio_exp = m_e_exp / m_mu_exp
print(f"  (m_e/m_μ)_exp     =          {nstr(ratio_exp, 20)}")

residual_section = ratio_exp / ratio_section - 1
print(f"  Residual: {nstr(residual_section, 10)} = {float(residual_section):.4e}")
print()

# ============================================================
# SPINOR LEVEL (degree 1): what if mass ~ |amplitude|, not |amplitude|²?
# ============================================================

print("=" * 70)
print("SPINOR LEVEL (degree 1): m_k ~ |f_k|, not f_k²")
print("=" * 70)
print()

# If the "true" mass at the spinor level is proportional to |f_k| (degree 1)
# rather than f_k² (degree 2), then:
# m_e/m_μ = |f_e|/|f_μ| = f_e/f_μ (all positive)

ratio_spinor = fabs(f_e) / fabs(f_mu)
print(f"  (m_e/m_μ)_spinor = |f_e|/|f_μ| = {nstr(ratio_spinor, 20)}")
print(f"  (m_e/m_μ)_exp    =               {nstr(ratio_exp, 20)}")

residual_spinor = ratio_exp / ratio_spinor - 1
print(f"  Residual: {nstr(residual_spinor, 10)} = {float(residual_spinor):.4e}")
print()

# Neither works directly — the spinor ratio is way off because
# m ~ |f| gives m_e/m_μ ~ 0.04/0.58 = 0.069, not 0.00484.
# The physical mass IS m ~ f², not m ~ |f|.

# ============================================================
# MIXED: mass = section, but RATIO correction from spinor geometry
# ============================================================

print("=" * 70)
print("INVESTIGATION: BORN RULE ON C^1 vs C^4")
print("=" * 70)
print()

# On C⁴ (H=3): Born probability p_i = |m_i|²/Σ|m_j|²
# This is the UNIQUE probability rule (Gleason, H≥3)
#
# On C² (spinor space): Gleason does NOT apply.
# Non-Born measures exist: p = |ψ|^α / Σ|ψ_j|^α for any α > 0
#
# On C¹: trivial (one state, p = 1)
#
# The Koide formula uses the section-level Born rule: m_k = M² × f_k²
# where f_k² is the Born weight of generation k.
#
# But if the lepton lives on C¹ ⊂ C² where Born isn't unique,
# the "effective exponent" might differ from 2.

print("If the mass formula has exponent α instead of 2:")
print("  m_k = M^α × |f_k|^α")
print("  m_e/m_μ = |f_e/f_μ|^α")
print()

# What α matches experiment?
# (m_e/m_μ)_exp = |f_e/f_μ|^α
# log(m_e/m_μ) = α × log(|f_e/f_μ|)

ln_ratio_exp = log(ratio_exp)
ln_f_ratio = log(fabs(f_e) / fabs(f_mu))

alpha_eff = ln_ratio_exp / ln_f_ratio
print(f"  ln(m_e/m_μ) = {nstr(ln_ratio_exp, 15)}")
print(f"  ln(|f_e/f_μ|) = {nstr(ln_f_ratio, 15)}")
print(f"  α_eff = {nstr(alpha_eff, 20)}")
print(f"  α_eff - 2 = {nstr(alpha_eff - 2, 10)}")
print(f"  (α_eff - 2)/2 = {float((alpha_eff - 2)/2):.6e}")
print()

# Check: does α_eff give a better m_e/m_μ?
ratio_alpha = (fabs(f_e) / fabs(f_mu))**alpha_eff
print(f"  (m_e/m_μ)_α = {nstr(ratio_alpha, 20)}")
print(f"  (m_e/m_μ)_exp = {nstr(ratio_exp, 20)}")
print(f"  Match: {float(fabs(ratio_alpha/ratio_exp - 1)):.2e}")
print()

# ============================================================
# KEY TEST: does α_eff have a framework interpretation?
# ============================================================

print("=" * 70)
print("FRAMEWORK INTERPRETATION OF α_eff")
print("=" * 70)
print()

print(f"  α_eff = {nstr(alpha_eff, 30)}")
print(f"  2     = {nstr(mpf(2), 30)}")
print(f"  δ = α_eff - 2 = {nstr(alpha_eff - 2, 20)}")
print()

delta = alpha_eff - 2

# Test against framework quantities
K_star = mpf(7) / mpf(30)
sigma = -log(mpf(23)/mpf(30))

tests = [
    ("1/H³", 1/H**3),
    ("K*²", K_star**2),
    ("K*/H³", K_star / H**3),
    ("σ/H³", sigma / H**3),
    ("1/(H³·Φ₄(H))", 1/(H**3 * (H**2+1))),
    ("K*·σ", K_star * sigma),
    ("σ²", sigma**2),
    ("1/h(E₈)", mpf(1)/30),
    ("K*/h(E₈)", K_star/30),
    ("1/(H²(H²+1))", 1/(H**2 * (H**2+1))),
    ("(H-1)/(H³(H²+1))", (H-1)/(H**3*(H**2+1))),
]

print(f"  Testing δ = {float(delta):.6e} against framework quantities:")
print()
for name, val in tests:
    ratio = delta / val if val != 0 else mpf(0)
    print(f"    δ / ({name:20s}) = {float(ratio):12.6f}   ({name} = {float(val):.6e})")

print()

# ============================================================
# DEEPER: the spinor-section map and its Jacobian
# ============================================================

print("=" * 70)
print("SPINOR → SECTION MAP")
print("=" * 70)
print()

# The map π: C² → Sym²(C²) = C³ sends (ζ,η) → (ζ²,ζη,η²)
# This is the Veronese embedding.
#
# At a point (ζ₀, η₀), the Jacobian is:
# dπ = [[2ζ, 0], [η, ζ], [0, 2η]]
#
# The "mass" of the η component (charged lepton) is |η|².
# But measured through the section map, it appears as s₃ = η².
# The section-level Born probability of s₃ is:
#   p₃ = |s₃|²/(|s₁|²+|s₂|²+|s₃|²) = |η|⁴/(|ζ|⁴+|ζ|²|η|²+|η|⁴)
#
# While the spinor-level "weight" of η is:
#   w_η = |η|²/(|ζ|²+|η|²)
#
# These are NOT the same! The Veronese embedding changes the Born weights.

print("Veronese embedding changes Born weights:")
print()
print("  Spinor level: w_η = |η|²/(|ζ|²+|η|²)")
print("  Section level: p₃ = |η|⁴/(|ζ|⁴+|ζη|²+|η|⁴)")
print()

# At the fixed point, the Koide x-values give the section weights.
# What are the corresponding spinor weights?
#
# If x_k = f_k², and the section is s_k ∝ f_k (in some sense),
# then the spinor amplitude is √f_k (degree 1/2??)
#
# Actually: the Koide factor f_k = 1 + √2 cos(θ+2πk/3) is the
# spinor-level amplitude. The mass is f_k² (section level).
# The question is whether the physical measurement projects through
# the Veronese map or directly accesses the spinor.

# At the vacuum: the three generations have spinor amplitudes f_e, f_μ, f_τ
# with Σf_k = 3 = H.
# The section amplitudes are f_k² with Σf_k² = 6 = 2H.

# Ratio at spinor level: f_e/f_μ
# Ratio at section level: f_e²/f_μ² = (f_e/f_μ)²

# The square root of the section ratio IS the spinor ratio:
# √(x_e/x_μ) = f_e/f_μ

# So: (m_e/m_μ) = (f_e/f_μ)² at section level (standard Koide)
# But: √(m_e/m_μ) = f_e/f_μ at spinor level

print("Spinor ratio: f_e/f_μ = √(m_e/m_μ)_Koide")
print(f"  = {nstr(f_e/f_mu, 20)}")
print(f"  √(m_e/m_μ)_exp = {nstr(sqrt(ratio_exp), 20)}")
spinor_residual = sqrt(ratio_exp) / (f_e/f_mu) - 1
print(f"  Spinor-level residual: {float(spinor_residual):.6e}")
print(f"  Section-level residual: {float(residual_section):.6e}")
print(f"  Ratio (section/spinor): {float(residual_section/spinor_residual):.4f}")
print()
print("  If residual is in the SQUARING (spinor→section),")
print("  the section residual should be ~2× the spinor residual.")
print(f"  Observed ratio: {float(residual_section/spinor_residual):.4f}")
print(f"  Expected if squaring: 2.0000")
print()

# ============================================================
# THE VERONESE DISTORTION
# ============================================================

print("=" * 70)
print("VERONESE DISTORTION: Born weights spinor vs section")
print("=" * 70)
print()

# At the Koide fixed point, the three generations have amplitudes:
# Spinor: (f_e, f_μ, f_τ) with Σ = 3
# Section: (f_e², f_μ², f_τ²) with Σ = 6
#
# Born probability at spinor level (on C³ of generations):
# p_k^(spinor) = f_k² / Σf_j² = f_k² / 6 = x_k / (2H)
#
# Born probability at section level (on C³ of sections):
# p_k^(section) = f_k⁴ / Σf_j⁴
#
# These are DIFFERENT. The Veronese embedding changes the Born weights.

sum_f4 = f_e**4 + f_mu**4 + f_tau**4
print(f"  Σf_k⁴ = {nstr(sum_f4, 15)}")
print()

p_e_spinor = f_e**2 / 6
p_mu_spinor = f_mu**2 / 6
p_tau_spinor = f_tau**2 / 6

p_e_section = f_e**4 / sum_f4
p_mu_section = f_mu**4 / sum_f4
p_tau_section = f_tau**4 / sum_f4

print(f"  Born weights (spinor level, p=f²/6):")
print(f"    p_e  = {float(p_e_spinor):.8f}")
print(f"    p_μ  = {float(p_mu_spinor):.8f}")
print(f"    p_τ  = {float(p_tau_spinor):.8f}")
print(f"    sum  = {float(p_e_spinor+p_mu_spinor+p_tau_spinor):.8f}")
print()

print(f"  Born weights (section level, p=f⁴/Σf⁴):")
print(f"    p_e  = {float(p_e_section):.8f}")
print(f"    p_μ  = {float(p_mu_section):.8f}")
print(f"    p_τ  = {float(p_tau_section):.8f}")
print(f"    sum  = {float(p_e_section+p_mu_section+p_tau_section):.8f}")
print()

# The distortion: ratio of section to spinor Born weights
print(f"  Veronese distortion (section/spinor Born weight):")
d_e = p_e_section / p_e_spinor
d_mu = p_mu_section / p_mu_spinor
d_tau = p_tau_section / p_tau_spinor
print(f"    d_e  = {float(d_e):.8f}")
print(f"    d_μ  = {float(d_mu):.8f}")
print(f"    d_τ  = {float(d_tau):.8f}")
print()

# The mass ratio changes by the distortion ratio:
# (m_e/m_μ)_section = (m_e/m_μ)_spinor × (d_e/d_μ)
veronese_correction = d_e / d_mu
print(f"  Veronese correction to m_e/m_μ: d_e/d_μ = {float(veronese_correction):.10f}")
print(f"  This is (m_e/m_μ)_section / (m_e/m_μ)_spinor")
print()

# Does the Veronese correction account for the residual?
# If the TRUE ratio is at the spinor level, and measurement
# projects through Veronese to the section level:
# (m_e/m_μ)_measured = (m_e/m_μ)_Koide × Veronese_correction
# residual = Veronese_correction - 1

print(f"  Veronese residual: {float(veronese_correction - 1):.6e}")
print(f"  Observed residual: {float(residual_section):.6e}")
print(f"  Ratio: {float((veronese_correction-1)/float(residual_section)):.2f}")

print()
print("=" * 70)
