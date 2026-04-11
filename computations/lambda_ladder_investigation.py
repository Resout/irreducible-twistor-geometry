#!/usr/bin/env python3
"""
The four-rung vacuum ladder: floor → depth → instanton → Λ

Test whether the structure is α × (H+1)² × α with computable α.
"""

from mpmath import mp, mpf, sqrt, ln, exp, log, fabs, nstr, pi, power

mp.dps = 50
H = 3
K_star = mpf(7)/30
FLOOR = mpf(1)/H**3

# The four rungs (as floor exponents)
# Floor = (1/27)^1
floor_exp = mpf(1)

# Depth: Born_DS = 0.000809487...
Born_DS = mpf("0.00080948751214898223625")
depth_exp = ln(Born_DS) / ln(FLOOR)

# Instanton: det(I-J)^{-1/2} × e^{-810/7}
lambda_0_ss = mpf("0.28291034660315")
lambda_1_ss = mpf("0.28131300001289")
S = mpf(810) / 7
prefactor = (1 - lambda_0_ss)**(-mpf(1)/2) * (1 - lambda_1_ss)**(-mpf(1)/2)
instanton = prefactor * exp(-S)
instanton_exp = ln(instanton) / ln(FLOOR)

# Observed Λ: ρ_Λ/M_Pl⁴ ≈ 1.15 × 10⁻¹²³
# More precisely: ρ_Λ = 5.35e-10 J/m³, M_Pl = 1.22e19 GeV
# In natural units: ρ_Λ ≈ 2.58e-47 GeV⁴, M_Pl⁴ = 2.21e76 GeV⁴
# ρ_Λ/M_Pl⁴ = 1.17e-123
Lambda_obs = mpf("1.17e-123")
Lambda_exp = ln(Lambda_obs) / ln(FLOOR)

print("=" * 70)
print("THE FOUR-RUNG VACUUM LADDER")
print("=" * 70)

print(f"\n  Floor:     (1/27)^{nstr(floor_exp, 10)}")
print(f"  Depth:     (1/27)^{nstr(depth_exp, 15)}")
print(f"  Instanton: (1/27)^{nstr(instanton_exp, 15)}")
print(f"  Λ_obs:     (1/27)^{nstr(Lambda_exp, 15)}")

steps = [depth_exp/floor_exp, instanton_exp/depth_exp, Lambda_exp/instanton_exp]
print(f"\n  Step ratios:")
print(f"    Floor → Depth:     ×{nstr(steps[0], 15)}")
print(f"    Depth → Instanton: ×{nstr(steps[1], 15)}")
print(f"    Instanton → Λ:    ×{nstr(steps[2], 15)}")

# ============================================================
# TEST: α × (H+1)² × α structure
# ============================================================
print("\n" + "=" * 70)
print("TEST: α × (H+1)² × α STRUCTURE")
print("=" * 70)

fermion_dim = mpf((H+1)**2)  # = 16
alpha_sq = Lambda_exp / (floor_exp * fermion_dim)
alpha = sqrt(alpha_sq)

print(f"\n  If Λ = floor^(α² × (H+1)²):")
print(f"  α² = {Lambda_exp}/(1 × 16) = {nstr(alpha_sq, 15)}")
print(f"  α = {nstr(alpha, 20)}")

# What IS α?
print(f"\n  Testing α against framework numbers:")

tests = {
    "ln(H²+1) = ln(10)": ln(mpf(H**2+1)),
    "ln(H³) = ln(27)": ln(mpf(H**3)),
    "H/√(H-1) = 3/√2": mpf(H)/sqrt(mpf(H-1)),
    "√(H²-H+1) = √7": sqrt(mpf(H**2-H+1)),
    "√(H(H+1)/2) = √6": sqrt(mpf(H*(H+1))/2),
    "depth_exp (= 2.16)": depth_exp,
    "1/K* = 30/7": mpf(1)/K_star,
    "H-K* = 2.767": mpf(H) - K_star,
    "√(H²+K*) = √(9.233)": sqrt(mpf(H**2) + K_star),
    "ln(H²+1)/ln(H²+1)... trivial": ln(mpf(10)),
    "S/5! × H = 810/7/120 × 3": S/120*H,
    "(H³-1)/H² × K* × H": mpf(H**3-1)/H**2 * K_star * H,
    "1 + 1/K* × (H-1)/H²": 1 + 1/K_star * mpf(H-1)/H**2,
    "√(1/K*) = √(30/7)": sqrt(1/K_star),
}

results = []
for name, val in tests.items():
    err = float(fabs(alpha - val) / alpha * 100)
    results.append((err, name, val))

results.sort()
for err, name, val in results[:10]:
    print(f"    α ≈ {name:40s} = {nstr(val, 15):>20s}  err={err:.4f}%")

# ============================================================
# TEST: other structural decompositions of the exponent
# ============================================================
print("\n" + "=" * 70)
print("ALTERNATIVE DECOMPOSITIONS OF floor^{85.9}")
print("=" * 70)

target_exp = Lambda_exp
print(f"\n  Target exponent: {nstr(target_exp, 15)}")

decomps = {
    "1 × 16 × α²": (1, 16, alpha_sq),
    "2.16 × 16.2 × 2.45": (steps[0], steps[1], steps[2]),
    "H × (H³-1) + H²": (None, None, mpf(H*(H**3-1) + H**2)),
    "S × K* × H = 810/7 × 7/30 × 3": (None, None, S * K_star * H),
    "H⁴ + H²/K* = 81 + 9×30/7": (None, None, mpf(H**4) + mpf(H**2)/K_star),
    "S/K* / (H+1)": (None, None, S/K_star/(H+1)),
    "(H²-1) × (H²+1) + H²-1": (None, None, mpf((H**2-1)*(H**2+1) + H**2-1)),
    "h(E₈) × (H-1) + H³": (None, None, mpf(30*(H-1) + H**3)),
    "5! × K* + S/K*/(H+1)": (None, None, 120*K_star + S/K_star/(H+1)),
    "H⁴-1+H²/K*": (None, None, mpf(H**4-1) + mpf(H**2)/K_star),
    "S²/5!": (None, None, S**2/120),
    "S × (1-K*) = 810/7 × 23/30": (None, None, S*(1-K_star)),
    "S + S²/5!": (None, None, S + S**2/120),
    "d_bos × H + H": (None, None, mpf(26*H + H)),
    "d_bos × (H+1/H)": (None, None, mpf(26)*(H + mpf(1)/H)),
}

results2 = []
for name, val_tuple in decomps.items():
    if val_tuple[0] is None:
        val = val_tuple[2]
    else:
        val = val_tuple[0] * val_tuple[1] * val_tuple[2]
    err = float(fabs(target_exp - val) / target_exp * 100)
    results2.append((err, name, float(val)))

results2.sort()
print("\n  Closest structural expressions:")
for err, name, val in results2[:12]:
    print(f"    {name:50s} = {val:>12.4f}  err={err:.4f}%")

# ============================================================
# THE SELF-SIMILAR HYPOTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SELF-SIMILAR HYPOTHESIS")
print("=" * 70)

print(f"""
  The three step ratios are: {nstr(steps[0],6)}, {nstr(steps[1],6)}, {nstr(steps[2],6)}

  If steps 1 and 3 are "the same process" (call it α):
    α_from_step1 = {nstr(steps[0], 10)}
    α_from_step3 = {nstr(steps[2], 10)}
    geometric mean = {nstr(sqrt(steps[0]*steps[2]), 10)}

  If step 2 is (H+1)² = 16:
    actual step 2 = {nstr(steps[1], 10)}
    (H+1)² = 16
    ratio = {nstr(steps[1]/16, 10)}
""")

# What if the "small step" α is exactly the depth exponent p = 2.16?
# Then Λ = floor^(p × 16 × p) = floor^(p² × 16)
p = depth_exp
predicted_exp = p**2 * 16
predicted_Lambda = FLOOR**predicted_exp

print(f"  If α = depth exponent p = {nstr(p, 10)}:")
print(f"    p² × 16 = {nstr(predicted_exp, 10)}")
print(f"    Λ_predicted = (1/27)^{nstr(predicted_exp, 8)} = {nstr(predicted_Lambda, 5)}")
print(f"    log10 = {nstr(ln(fabs(predicted_Lambda))/ln(10), 10)}")
print(f"    Λ_observed ≈ 10^-123")
print(f"    Discrepancy: {nstr(fabs(ln(fabs(predicted_Lambda))/ln(10) - (-123)), 5)} dex")

# What if the big step is NOT (H+1)² but something else?
# If the structure is α, β, α with α = depth_exp:
beta = instanton_exp / depth_exp
print(f"\n  If structure is α, β, α with α = {nstr(depth_exp, 6)}:")
print(f"    β = {nstr(beta, 10)}")
print(f"    α²β = {nstr(depth_exp**2 * beta, 10)}")
print(f"    floor^(α²β) = floor^{nstr(depth_exp**2 * beta, 6)} = 10^{nstr(depth_exp**2 * beta * ln(FLOOR)/ln(10), 6)}")

# β in framework terms
beta_tests = {
    "(H+1)² = 16": mpf(16),
    "H(H+2) = 15": mpf(15),
    "(H²+H+2)/K*×something": mpf(16),
    "S/ln(1/Born_DS) = S/7.12": S/ln(1/Born_DS),
    "S/H² = 810/63": S/H**2,
    "(H+1)²+K*²×H²": mpf(16) + K_star**2*H**2,
}

print(f"\n  β = {nstr(beta, 10)} vs framework:")
for name, val in beta_tests.items():
    err = float(fabs(beta - val)/beta*100)
    print(f"    β ≈ {name:35s} = {nstr(val, 10):>12s}  err={err:.3f}%")

# ============================================================
# EXACT PREDICTION TEST
# ============================================================
print("\n" + "=" * 70)
print("PREDICTION: Λ FROM THE LADDER")
print("=" * 70)

# If Λ = floor^(depth_exp × instanton_exp/depth_exp × depth_exp)
# = floor^(depth_exp² × β)
# And β = S / ln(1/Born_DS):

beta_structural = S / ln(1/Born_DS)
exp_predicted = depth_exp**2 * beta_structural
Lambda_predicted = FLOOR**exp_predicted

print(f"\n  β = S / ln(1/Born_DS) = {nstr(beta_structural, 15)}")
print(f"  Predicted exponent = p² × β = {nstr(exp_predicted, 10)}")
print(f"  Predicted Λ = (1/27)^{nstr(exp_predicted, 6)}")
print(f"  log10(Λ_predicted) = {nstr(exp_predicted * ln(FLOOR)/ln(10), 10)}")
print(f"  log10(Λ_observed) ≈ -123")

# Also test: S itself as the bridge
# S = H³/K* = 810/7 = 115.7
# S / ln(floor) = 810/7 / 3.296 = 35.1 ≈ instanton_exp
S_over_logfloor = S / (-ln(FLOOR))
print(f"\n  S / |ln(floor)| = {nstr(S_over_logfloor, 15)}")
print(f"  instanton_exp = {nstr(instanton_exp, 15)}")
print(f"  ratio = {nstr(instanton_exp/S_over_logfloor, 10)}")

# The instanton exponent ≈ S / |ln(floor)| = H³/K* / ln(H³)
# = H³ / (K* × ln(H³))
# = H³ / (K* × 3ln(H))

print(f"\n  Instanton exp ≈ H³/(K* × 3ln(H)) = 27/(7/30 × 3×ln(3))")
print(f"    = {nstr(mpf(27)/(K_star * 3 * ln(mpf(3))), 15)}")
print(f"  vs actual = {nstr(instanton_exp, 15)}")


if __name__ == "__main__":
    main()
