#!/usr/bin/env python3
"""
Neutrino Koide angle investigation.

The charged leptons see Z₃ as rotation: θ_ℓ = 2/9, giving steep hierarchy.
The substrate sees Z₃ as contraction: θ_ν = ???, giving flat hierarchy.

Key question: what angle gives the observed neutrino mass ratios?

Observed:
  Δm²₂₁ = 7.53e-5 eV² (solar)
  Δm²₃₂ = 2.453e-3 eV² (atmospheric)
  m₃ ≈ 0.050 eV (from atmospheric)
  m₂ ≈ 0.0087 eV (from solar)
  m₁ unknown (could be ~0 or ~m₂)

Normal hierarchy: m₁ << m₂ << m₃
"""

from mpmath import mp, mpf, sqrt, cos, sin, pi, nstr, fabs, acos, ln

mp.dps = 30
H = 3
K = mpf(7)/30

# ============================================================
# WHAT ANGLE FITS THE NEUTRINO DATA?
# ============================================================
print("=" * 70)
print("REVERSE-ENGINEERING THE NEUTRINO KOIDE ANGLE")
print("=" * 70)

# The Koide formula: m_k = M² (1 + √2 cos(θ + 2πk/3))²
# The x-factors: x_k = (1 + √2 cos(θ + 2πk/3))²
# Sum: Σx_k = 6 (always, for any θ)
# The hierarchy depends on θ.

# For charged leptons: θ = 2/9 ≈ 0.222
# x_e = 0.00163, x_μ = 0.337, x_τ = 5.662
# Ratio x_τ/x_μ = 16.8, x_μ/x_e = 207

# For neutrinos, we need a MUCH flatter hierarchy.
# If normal ordering: m₃ >> m₂ >> m₁
# m₃ ≈ 0.050 eV, m₂ ≈ 0.0087 eV
# m₃/m₂ ≈ 5.7

# The x-factors are proportional to masses:
# x₃/x₂ = m₃/m₂ ≈ 5.7

# What θ gives x₃/x₂ ≈ 5.7?
# Let's scan θ from 0 to π/3

print("\nScanning θ for neutrino-like hierarchy (x_max/x_mid ≈ 5.7):")
print(f"{'θ':>10s} {'θ (frac)':>12s} {'x₁':>10s} {'x₂':>10s} {'x₃':>10s} {'x₃/x₂':>8s} {'x₂/x₁':>8s}")

best_theta = None
best_err = 100

for i in range(1, 300):
    th = mpf(i) / 300 * pi / 3  # scan 0 to π/3
    x = sorted([(1 + sqrt(2)*cos(th + 2*pi*k/3))**2 for k in range(3)])
    if float(x[0]) < 1e-6:
        continue
    ratio_32 = x[2]/x[1]
    ratio_21 = x[1]/x[0]

    # Target: x₃/x₂ ≈ 5.7
    err = abs(float(ratio_32) - 5.7)
    if err < best_err:
        best_err = err
        best_theta = th
        best_x = x
        best_ratios = (ratio_32, ratio_21)

print(f"\nBest match for x₃/x₂ ≈ 5.7:")
print(f"  θ = {nstr(best_theta, 10)} rad")
print(f"  x = ({nstr(best_x[0],6)}, {nstr(best_x[1],6)}, {nstr(best_x[2],6)})")
print(f"  x₃/x₂ = {nstr(best_ratios[0], 6)}")
print(f"  x₂/x₁ = {nstr(best_ratios[1], 6)}")

# ============================================================
# TEST SPECIFIC FRAMEWORK ANGLES
# ============================================================
print("\n" + "=" * 70)
print("FRAMEWORK ANGLE CANDIDATES")
print("=" * 70)

candidates = {
    "2/9 (charged lepton)": mpf(2)/9,
    "7/60 = K*/(H-1)": mpf(7)/60,
    "π/4 - 2/9": pi/4 - mpf(2)/9,
    "2/(H³) = 2/27": mpf(2)/27,
    "1/(H(H+1)) = 1/12": mpf(1)/12,
    "K*/H = 7/90": K/H,
    "2/(H(H+1)) = 1/6": mpf(1)/6,
    "(H-1)/H² = 2/9 (same as lepton!)": mpf(2)/9,
    "π/H² = π/9": pi/9,
    "1/H = 1/3": mpf(1)/3,
    "π/(H(H+1)) = π/12": pi/12,
    "2/H³ = 2/27": mpf(2)/27,
    "K*²/(H-1) = 49/1800": K**2/(H-1),
    "α_floor × 2/9": mpf("0.868") * mpf(2)/9,
    "ln(1/6) (compression log)": -ln(mpf(1)/6),
    "(1/6)×(2/9) = 1/27": mpf(1)/27,
    "θ_DS/θ* × π = π/6": pi/6,
}

print(f"\n{'Angle':40s} {'θ':>8s} {'x₃/x₂':>8s} {'x₂/x₁':>8s} {'m₃(eV)':>8s}")
print("-" * 70)

# Use m₃ ≈ 0.050 eV to set scale
m3_target = mpf("0.050")

for name, th in candidates.items():
    x = sorted([(1 + sqrt(2)*cos(th + 2*pi*k/3))**2 for k in range(3)])
    if float(x[0]) < 1e-10:
        continue
    r32 = x[2]/x[1]
    r21 = x[1]/x[0]

    # M² = m₃/x₃ (using largest mass = largest x)
    M2 = m3_target / x[2]
    m2 = M2 * x[1]
    m1 = M2 * x[0]

    # Check against observed Δm² values
    dm21_sq = float(m2**2 - m1**2)  # should be ~7.53e-5 eV²
    dm32_sq = float(m3_target**2 - m2**2)  # should be ~2.453e-3 eV²

    flag = ""
    if abs(dm21_sq - 7.53e-5)/7.53e-5 < 0.3 and abs(dm32_sq - 2.453e-3)/2.453e-3 < 0.3:
        flag = " *** MATCH ***"

    print(f"  {name:38s} {float(th):8.4f} {float(r32):8.3f} {float(r21):8.3f} {float(m3_target):8.4f}{flag}")

# ============================================================
# FINE SCAN NEAR THE BEST MATCH
# ============================================================
print("\n" + "=" * 70)
print("FINE SCAN: MATCHING Δm² VALUES")
print("=" * 70)

# Target: Δm²₂₁ = 7.53e-5 eV², Δm²₃₂ = 2.453e-3 eV²
# Ratio: Δm²₂₁/Δm²₃₂ = 0.0307 (observed)
# With m₃ = 0.050:
# m₂² = m₃² - Δm²₃₂ = 0.0025 - 0.002453 = 0.000047 → m₂ = 0.00686 eV...
# wait: m₂ = sqrt(m₃² - Δm²₃₂) = sqrt(0.0025 - 0.002453) = sqrt(0.000047) = 0.00686
# Then m₁² = m₂² - Δm²₂₁ = 0.000047 - 0.0000753 < 0 !!

# Hmm, that means m₁² < 0 which is impossible. Let me reconsider.
# Normal ordering: m₁ < m₂ < m₃
# Δm²₃₂ = m₃² - m₂² = 2.453e-3 → m₂ = sqrt(m₃² - 2.453e-3)
# Δm²₂₁ = m₂² - m₁² = 7.53e-5 → m₁ = sqrt(m₂² - 7.53e-5)

m3 = mpf("0.050")
dm32 = mpf("2.453e-3")
dm21 = mpf("7.53e-5")

m2_sq = m3**2 - dm32
m2 = sqrt(m2_sq) if float(m2_sq) > 0 else mpf(0)
m1_sq = m2_sq - dm21
m1 = sqrt(m1_sq) if float(m1_sq) > 0 else mpf(0)

print(f"\nFrom observed splittings (normal ordering):")
print(f"  m₃ = {nstr(m3, 6)} eV")
print(f"  m₂ = {nstr(m2, 6)} eV")
print(f"  m₁ = {nstr(m1, 6)} eV ({'real' if float(m1_sq)>0 else 'IMAGINARY - m₃ too small!'})")

if float(m1_sq) <= 0:
    # Try larger m₃
    print("\n  m₃ = 0.050 is too small for normal ordering with these splittings.")
    print("  Trying m₃ = 0.051 eV:")
    m3 = mpf("0.051")
    m2_sq = m3**2 - dm32
    m2 = sqrt(m2_sq)
    m1_sq = m2_sq - dm21
    m1 = sqrt(m1_sq) if float(m1_sq) > 0 else mpf(0)
    print(f"  m₃ = {nstr(m3, 6)}, m₂ = {nstr(m2, 6)}, m₁ = {nstr(m1, 6)}")

print(f"\n  Mass ratios:")
print(f"    m₃/m₂ = {nstr(m3/m2, 6)}")
if float(m1) > 0:
    print(f"    m₂/m₁ = {nstr(m2/m1, 6)}")
    print(f"    m₃/m₁ = {nstr(m3/m1, 6)}")

# Now find the Koide angle that reproduces these ratios
print(f"\n  Target x-factor ratios:")
print(f"    x₃/x₂ = m₃/m₂ = {nstr(m3/m2, 6)}")
if float(m1) > 0:
    print(f"    x₂/x₁ = m₂/m₁ = {nstr(m2/m1, 6)}")

# Fine scan
target_ratio = float(m3/m2)
print(f"\n  Fine scanning for x₃/x₂ = {target_ratio:.4f}:")

best_th = None
best_e = 100

for i in range(1, 10000):
    th = mpf(i)/10000 * pi/3
    x = sorted([(1 + sqrt(2)*cos(th + 2*pi*k/3))**2 for k in range(3)])
    r = float(x[2]/x[1])
    e = abs(r - target_ratio)
    if e < best_e:
        best_e = e
        best_th = th
        best_x = x

print(f"  Best θ = {nstr(best_th, 15)}")
print(f"  x = ({nstr(best_x[0],8)}, {nstr(best_x[1],8)}, {nstr(best_x[2],8)})")
print(f"  x₃/x₂ = {nstr(best_x[2]/best_x[1], 8)} (target {target_ratio:.4f})")
if float(best_x[0]) > 1e-10:
    print(f"  x₂/x₁ = {nstr(best_x[1]/best_x[0], 8)}")

# What framework fraction is this angle?
print(f"\n  θ_ν = {nstr(best_th, 15)}")
print(f"  Testing against framework fractions:")
frac_tests = {
    "2/9": mpf(2)/9,
    "1/6": mpf(1)/6,
    "7/60": mpf(7)/60,
    "π/4-2/9": pi/4-mpf(2)/9,
    "1/H²": mpf(1)/H**2,
    "K*/H": K/H,
    "2/H³": mpf(2)/H**3,
    "K*/(H+1)": K/(H+1),
    "2/(H²+1)": mpf(2)/(H**2+1),
    "1/(H+1)": mpf(1)/(H+1),
    "K*²": K**2,
    "2K*/H": 2*K/H,
    "(π-2)/H²": (pi-2)/H**2,
    "θ_ℓ/(H-1) = 1/9": mpf(1)/9,
    "θ_ℓ/H = 2/27": mpf(2)/27,
    "π/(2H²) = π/18": pi/18,
    "K*/(H(H-1)) = 7/180": K/(H*(H-1)),
}

results = []
for name, val in frac_tests.items():
    err = float(fabs(best_th - val)/best_th * 100)
    results.append((err, name, val))

results.sort()
for err, name, val in results[:8]:
    print(f"    {name:25s} = {nstr(val, 10):>12s}  err={err:.2f}%")
