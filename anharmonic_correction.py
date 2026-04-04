"""
Compute the anharmonic correction to the Gaussian mass gap.

Gaussian: Δ_G = 1/(1-K*) = 30/23 ≈ 1.3043
Exact:    Δ   = 1.2626
Ratio:    Δ/Δ_G = 0.9680

The 3.2% correction should come from the cubic and higher terms
in the expansion of S(K) = -ln(1-K) around K*.

In the saddle-point approximation with anharmonic corrections:
the effective mass gap receives corrections from the skewness
of the action around the saddle.
"""
import numpy as np
from fractions import Fraction

K_star = Fraction(7, 30)
one_minus_K = 1 - K_star  # = 23/30

print("S(K) = -ln(1-K) expanded around K*=7/30:")
print(f"  S' = 1/(1-K*) = {1/one_minus_K} = {float(1/one_minus_K):.6f}")
print(f"  S'' = 1/(1-K*)² = {1/one_minus_K**2} = {float(1/one_minus_K**2):.6f}")
print(f"  S''' = 2/(1-K*)³ = {2/one_minus_K**3} = {float(2/one_minus_K**3):.6f}")
print(f"  S'''' = 6/(1-K*)⁴ = {6/one_minus_K**4} = {float(6/one_minus_K**4):.6f}")
print()

# Gaussian mass gap
Delta_G = float(1/one_minus_K)
Delta_exact = 1.2626

print(f"Gaussian Δ_G = 1/(1-K*) = {Delta_G:.6f}")
print(f"Exact Δ = {Delta_exact:.6f}")
print(f"Ratio = {Delta_exact/Delta_G:.6f}")
print(f"Correction = {1 - Delta_exact/Delta_G:.6f} = {(1-Delta_exact/Delta_G)*100:.2f}%")
print()

# The anharmonic correction to the mass gap in the steepest descent:
# For S(K) = -ln(1-K), all derivatives are known:
# S^(n)(K*) = (n-1)! / (1-K*)^n
#
# The leading correction to the Gaussian eigenvalue comes from
# the quartic coupling (S'''' / S''^2):
#
# In a 1D integral Z = ∫ e^{-n·S(K)} dK ≈ e^{-n·S*} √(2π/nS'')
# the correlation function <δK²> = 1/(nS'')
# The mass gap is related to 1/<δK²> = nS''
# The leading correction: 1/(nS'') → 1/(nS'') · (1 - correction)
#
# But our system is not a simple 1D integral. The transfer operator
# eigenvalue λ₀ = 0.2829 comes from the FULL nonlinear dynamics,
# not just from the action curvature.

# Let me check: is there a CLEAN relationship between Δ and 1/(1-K*)?

# Exact: λ₀ = 0.2829102944
# Note: (1-K*)^{something} = λ₀?
# (23/30)^x = 0.2829
# x = ln(0.2829)/ln(23/30) = 1.2626/0.2657 = 4.752

# So λ₀ = (1-K*)^{4.752}

# Is 4.752 related to anything?
ratio = Delta_exact / float(-np.log(float(one_minus_K)))
print(f"λ₀ = (1-K*)^r where r = {ratio:.6f}")
print(f"  = (23/30)^{ratio:.6f}")
print()

# Is r close to a nice number?
# r ≈ 4.752
# H² + 1 - H = 9 + 1 - 3 = 7? No.
# (H²+1)/2 = 5? No.
# H! - 1 = 5? No.
# dim(Sym²) - dim(C^{H+1}) = 10 - 4 = 6? No.

# Let me check: r = Δ/Δ_filter = Δ/(-ln(1-K*))
# This is just the ratio we already know: 4.752

# What if r = dim(Sym²)/2 - 1/(1-K*)?
# 10/2 - 30/23 = 5 - 1.304 = 3.696? No.

# What if r relates to the Hessian eigenvalues?
# Hessian of S on the L1=1 surface had eigenvalues 0.675, ~0, -0.328
# These don't obviously give 4.752.

# Let me try: is Δ = -ln(K*) ?
print(f"-ln(K*) = -ln(7/30) = {-np.log(7/30):.6f}")
print(f"Δ = {Delta_exact:.6f}")
print(f"Match? No ({-np.log(7/30):.3f} vs {Delta_exact:.3f})")
print()

# Is Δ = -ln(K*/(1-K*))?
print(f"-ln(K*/(1-K*)) = -ln(7/23) = {-np.log(7/23):.6f}")
print(f"Δ = {Delta_exact:.6f}")
print(f"Match? No ({-np.log(7/23):.3f} vs {Delta_exact:.3f})")
print()

# Is Δ = ln(H) + ln(1/(1-K*))?
print(f"ln(H) + ln(1/(1-K*)) = ln(3) + ln(30/23) = {np.log(3) + np.log(30/23):.6f}")
print(f"Δ = {Delta_exact:.6f}")
print(f"Match? {abs(np.log(3) + np.log(30/23) - Delta_exact) < 0.01}")
print()

# ln(3) + ln(30/23) = ln(90/23) = ln(3.913) = 1.364. Not 1.263.

# What about Δ = ln(H²+1) - ln(H)?
print(f"ln(H²+1) - ln(H) = ln(10) - ln(3) = {np.log(10) - np.log(3):.6f}")
print(f"Δ = {Delta_exact:.6f}")
print(f"Match? {abs(np.log(10) - np.log(3) - Delta_exact) < 0.02}")
# ln(10/3) = 1.204. Close but not quite.
print()

# ln(10/3) = 1.204, Δ = 1.263. Ratio = 1.049.

# What about Δ = ln(dim(Sym²)/H)?
print(f"ln(dim(Sym²)/H) = ln(10/3) = {np.log(10/3):.6f}")
print(f"Δ = {Delta_exact:.6f}")
print(f"Difference: {Delta_exact - np.log(10/3):.4f} ({(Delta_exact/np.log(10/3)-1)*100:.1f}%)")
print()

# Δ ≈ ln(10/3) to 4.9%. Hmm, not as clean as 1/(1-K*) which is 3.2%.

# The cleanest relationship remains:
# Δ_exact ≈ 1/(1-K*) × 0.968
# or equivalently: λ₀ ≈ e^{-1/(1-K*)} × e^{0.041}

# The 0.968 factor: is this (1 - 1/dim(Sym²))?
print(f"1 - 1/dim(Sym²) = 1 - 1/10 = 0.9")
print(f"Δ/Δ_G = {Delta_exact/Delta_G:.4f}")
print(f"Match? No (0.9 vs 0.968)")
print()

# Is it (1 - K*)?
print(f"1 - K* = {float(one_minus_K):.6f}")
print(f"Δ/Δ_G = {Delta_exact/Delta_G:.6f}")
print(f"(1-K*) × Δ_G = {float(one_minus_K) * Delta_G:.6f}")
print(f"This equals 1. Trivially: (1-K*) × 1/(1-K*) = 1.")
print()

# Summary
print("=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"  K* = 7/30 (algebraic, from Sym²(C⁴))")
print(f"  S* = -ln(23/30) = 0.2657 (action at saddle)")
print(f"  Δ_G = 1/(1-K*) = 30/23 = 1.3043 (Gaussian gap)")
print(f"  Δ   = 1.2626 (exact, from eigenvalue)")
print(f"  Correction: 3.2% (anharmonic)")
print(f"  λ₀ = (1-K*)^{{4.752}} = (23/30)^{{4.752}} (exact)")
print()
print("  The 3.2% anharmonic correction does not appear to have")
print("  a clean closed form. It's a numerical property of the")
print("  specific DS dynamics at the K*=7/30 equilibrium.")
print()
print("  The Gaussian approximation Δ ≈ 1/(1-K*) is the")
print("  leading-order saddle-point result. The exact gap")
print("  requires the full nonlinear eigenvalue computation.")
