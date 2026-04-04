"""
Investigation: what is the exact relationship between K*=7/30 and γ_DLM?

Three possibilities:
1. K* IS γ (exact), and γ_DLM has an error in the counting
2. K* and γ are related by an exact transformation we haven't found
3. The match is a coincidence

Let's compute γ_DLM from scratch and look for structure.
"""
import numpy as np
from fractions import Fraction

# ============================================================
# PART 1: Compute γ_DLM from the LQG state-counting equation
# ============================================================
# The Domagala-Lewandowski-Meissner value is the unique positive
# root of the equation relating black hole entropy to area.
#
# The standard form: γ is determined by
#   Σ_{j=1/2,1,3/2,...} (2j+1) exp(-2γ₀ √(j(j+1))) = 1
# where γ₀ = π γ.
#
# Actually, the precise equation from Meissner (2004) is:
# The number of states for a horizon of area A is dominated by
# configurations where each puncture carries j=1/2.
# For the full counting (all j):
#   γ_DLM is the unique positive solution of
#   Σ_j (2j+1) x^{√(j(j+1))} = 1  where x = exp(-2πγ)
#
# But different sources give slightly different forms.
# Let me use the definition that gives 0.23753295796592...

def gamma_dlm_equation(gamma, n_terms=200):
    """The DLM equation: Σ_j (2j+1) exp(-2πγ √(j(j+1))) = 1
    Sum over j = 1/2, 1, 3/2, 2, ... (half-integers)
    """
    total = 0.0
    for k in range(1, n_terms + 1):
        j = k / 2.0  # j = 1/2, 1, 3/2, 2, ...
        total += (2*j + 1) * np.exp(-2 * np.pi * gamma * np.sqrt(j*(j+1)))
    return total - 1.0

# Find γ_DLM by bisection
lo, hi = 0.1, 0.5
for _ in range(200):
    mid = (lo + hi) / 2
    if gamma_dlm_equation(mid) > 0:
        lo = mid
    else:
        hi = mid

gamma_dlm = (lo + hi) / 2
print(f"γ_DLM = {gamma_dlm:.15f}")
print(f"K*    = {7/30:.15f}")
print(f"Ratio γ/K* = {gamma_dlm / (7/30):.15f}")
print(f"Gap = {gamma_dlm - 7/30:.15f}")
print(f"Gap/K* = {(gamma_dlm - 7/30)/(7/30):.6f} = {(gamma_dlm - 7/30)/(7/30)*100:.4f}%")
print()

# ============================================================
# PART 2: Look at the ratio γ/K* for structure
# ============================================================
r = gamma_dlm / (7/30)
print(f"γ_DLM / K* = {r:.15f}")
print()

# Is r related to known constants?
candidates = {
    "1 + 1/(8π²)": 1 + 1/(8*np.pi**2),
    "1 + 1/(4π²)": 1 + 1/(4*np.pi**2),
    "1 + 1/H⁴": 1 + 1/81,
    "1 + K*²": 1 + (7/30)**2,
    "1 + K*/H²": 1 + (7/30)/9,
    "1 + 1/(2H³)": 1 + 1/54,
    "1 + η/H": 1 + (4/27)/3,
    "1 + 1/(H²(H+1))": 1 + 1/36,
    "1 + 1/(H(H²+1))": 1 + 1/30,
    "exp(K*²)": np.exp((7/30)**2),
    "1 + K*/(H²+1)": 1 + (7/30)/10,
    "(H²+1)/(H²+1-K*)": 10/(10-7/30),
    "30/23 × 23/30 × γ/K*": r,  # tautology check
    "1 + 1/(H⁵)": 1 + 1/243,
    "1 + K*³": 1 + (7/30)**3,
    "(1-K*²)^(-1/2)": 1/np.sqrt(1-(7/30)**2),
    "1/(1-K*²/2)": 1/(1-(7/30)**2/2),
}

print("Candidate expressions for γ_DLM/K*:")
for name, val in sorted(candidates.items(), key=lambda x: abs(x[1]-r)):
    diff = abs(val - r) / r * 1e6  # ppm
    if diff < 10000:  # within 1%
        print(f"  {name:30s} = {val:.10f}  ({diff:.0f} ppm)")

print()

# ============================================================
# PART 3: What if the relationship is NOT multiplicative?
# ============================================================
# Maybe γ = f(K*) where f is not K* × constant.
# Try: γ = K_cons(H_eff) for some H_eff

def K_cons(H):
    """Conservation law: K = (H²-H+1)/(H(H²+1))"""
    return (H**2 - H + 1) / (H * (H**2 + 1))

# At what H_eff does K_cons(H_eff) = γ_DLM?
from scipy.optimize import brentq
H_eff_for_gamma = brentq(lambda h: K_cons(h) - gamma_dlm, 2.5, 3.5)
print(f"K_cons(H_eff) = γ_DLM at H_eff = {H_eff_for_gamma:.10f}")
print(f"H - H_eff = {3 - H_eff_for_gamma:.10f}")
print(f"(H - H_eff)/K* = {(3 - H_eff_for_gamma)/(7/30):.10f}")
print(f"(H - H_eff) × H = {(3 - H_eff_for_gamma) * 3:.10f}")
print()

# The paper's dressing: H_eff = H - γ/H
# Check: 3 - γ_DLM/3 = 3 - 0.07918 = 2.92082
H_eff_paper = 3 - gamma_dlm/3
print(f"Paper's H_eff = H - γ/H = {H_eff_paper:.10f}")
print(f"K_cons(paper's H_eff) = {K_cons(H_eff_paper):.15f}")
print(f"γ_DLM               = {gamma_dlm:.15f}")
print(f"Difference = {K_cons(H_eff_paper) - gamma_dlm:.6e}")
print()

# The gap in H_eff:
print(f"Exact H_eff for γ_DLM:     {H_eff_for_gamma:.10f}")
print(f"Paper's H_eff (α=1/H):     {H_eff_paper:.10f}")
print(f"Difference in H_eff:        {H_eff_for_gamma - H_eff_paper:.6e}")
print()

# ============================================================
# PART 4: What α gives exact match?
# ============================================================
# H_eff = H - α·γ
# We need K_cons(3 - α·γ_DLM) = γ_DLM
# Solve for α

alpha_exact = brentq(lambda a: K_cons(3 - a*gamma_dlm) - gamma_dlm, 0.1, 1.0)
print(f"Exact α for γ_DLM: {alpha_exact:.15f}")
print(f"1/H = 1/3 =        {1/3:.15f}")
print(f"Difference: {alpha_exact - 1/3:.6e} ({(alpha_exact - 1/3)/(1/3)*100:.4f}%)")
print()

# Is α_exact a nice number?
print("Is α_exact recognizable?")
candidates_alpha = {
    "1/3": 1/3,
    "1/3 + 1/H⁴": 1/3 + 1/81,
    "1/3 + K*/H³": 1/3 + (7/30)/27,
    "1/3 + 1/(H²(H²+1))": 1/3 + 1/90,
    "K*/(1-K*)": (7/30)/(23/30),
    "7/23": 7/23,
    "1/(3-K*)": 1/(3-7/30),
    "10/30": 10/30,
    "H/(H²-1)": 3/8,
    "K*/ln(1-K*)": -(7/30)/np.log(23/30),
    "1/e": 1/np.e,
    "ln(3)/π": np.log(3)/np.pi,
}

for name, val in sorted(candidates_alpha.items(), key=lambda x: abs(x[1]-alpha_exact)):
    diff_ppm = abs(val - alpha_exact) / alpha_exact * 1e6
    if diff_ppm < 50000:
        print(f"  {name:25s} = {val:.10f}  ({diff_ppm:.0f} ppm)")

print()

# ============================================================
# PART 5: Is γ_DLM = 7/30 × (something from H=3)?
# ============================================================
print("=" * 60)
print("PART 5: Exact expressions")
print("=" * 60)

# γ_DLM = 0.23753295796592
# 7/30 = 0.23333333333333

# Difference = 0.00419962463259
# Is this a nice number?
diff = gamma_dlm - 7/30
print(f"γ_DLM - K* = {diff:.15f}")
print()

# Check: is the difference related to the Born floor?
print(f"1/H³ = {1/27:.15f}")
print(f"diff/floor = {diff/(1/27):.10f}")
print(f"diff × H³ = {diff * 27:.10f}")
print(f"diff × H² = {diff * 9:.10f}")
print(f"diff × 30 = {diff * 30:.10f}")
print(f"diff × 30 × H = {diff * 90:.10f}")
print()

# diff × 30 ≈ 0.126 ≈ ?
# diff × 90 ≈ 0.378 ≈ ?

# What if γ = K* + K*²?
print(f"K* + K*² = {7/30 + (7/30)**2:.15f}")
print(f"γ_DLM   = {gamma_dlm:.15f}")
print(f"Match? {abs(7/30 + (7/30)**2 - gamma_dlm) < 0.001}")
print()

# K* + K*² = 0.2333 + 0.0544 = 0.2878... too high.

# What if γ = K*/(1-K*) × something?
# K*/(1-K*) = 7/23 = 0.30435
# γ/[K*/(1-K*)] = 0.23753/0.30435 = 0.78054

# What about the self-consistent equation directly?
# γ = K_cons(H - γ/H)
# This is γ = f(γ) where f(x) = K_cons(3 - x/3)
# The fixed point is γ_sc.

print("=" * 60)
print("PART 6: The self-consistency iteration")
print("=" * 60)

gamma = 7/30
for i in range(20):
    H_eff = 3 - gamma/3
    gamma_new = K_cons(H_eff)
    print(f"  iter {i:2d}: γ = {gamma:.15f}, H_eff = {H_eff:.10f}")
    if abs(gamma_new - gamma) < 1e-16:
        break
    gamma = gamma_new

gamma_sc = gamma
print(f"\nγ_sc = {gamma_sc:.15f}")
print(f"γ_DLM = {gamma_dlm:.15f}")
print(f"Gap = {gamma_dlm - gamma_sc:.6e}")
print(f"Gap/γ = {(gamma_dlm - gamma_sc)/gamma_dlm * 1e6:.1f} ppm")
print()

# ============================================================
# PART 7: What if the dressing is H_eff = H - γ²/H?
# ============================================================
print("=" * 60)
print("PART 7: Alternative dressings")
print("=" * 60)

for name, dress_fn in [
    ("α=1/H (paper)", lambda g: 3 - g/3),
    ("α=K*", lambda g: 3 - (7/30)*g),
    ("α=1/(H-1)", lambda g: 3 - g/2),
    ("α=1/H²", lambda g: 3 - g/9),
    ("quadratic: H-γ²", lambda g: 3 - g**2),
    ("quadratic: H-γ²/H", lambda g: 3 - g**2/3),
    ("log: H-ln(1+γ)", lambda g: 3 - np.log(1+g)),
    ("sqrt: H-√γ/H", lambda g: 3 - np.sqrt(g)/3),
    ("exact α*", lambda g: 3 - 0.33987*g),
]:
    gamma = 7/30
    for _ in range(100):
        h = dress_fn(gamma)
        if h < 1 or h > 10:
            gamma = float('nan')
            break
        gamma_new = K_cons(h)
        if abs(gamma_new - gamma) < 1e-16:
            break
        gamma = gamma_new

    if not np.isnan(gamma):
        gap_ppm = abs(gamma - gamma_dlm)/gamma_dlm * 1e6
        print(f"  {name:25s}: γ = {gamma:.12f}  gap = {gap_ppm:.1f} ppm")
    else:
        print(f"  {name:25s}: diverged")
