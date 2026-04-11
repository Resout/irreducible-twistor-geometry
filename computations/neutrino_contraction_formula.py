#!/usr/bin/env python3
"""
The neutrino mass formula: contraction instead of rotation.

Charged leptons (sections, Z₃ rotation):
  m_k = M² (1 + √2 cos(θ + 2πk/3))²

The substrate sees Z₃ as a contraction (shrink), not a rotation.
What's the contraction analog?

Approach: replace cos with cosh, or replace the circular Z₃ phases
with contractive ones. The substrate projects the rotation into
amplitude modulation.

Key constraint: Δm²₂₁/Δm²₃₂ = 1/33 (Prediction 20 of the paper)
"""

from mpmath import mp, mpf, sqrt, cos, sin, cosh, sinh, exp, pi, nstr, fabs, ln

mp.dps = 30
H = 3
K = mpf(7)/30

# ============================================================
# THE CONTRACTION FORMULA
# ============================================================
print("=" * 70)
print("CONTRACTION FORMULAS FOR NEUTRINO MASSES")
print("=" * 70)

# Approach 1: Replace cos → cosh
# m_k = M² (1 + √2 cosh(θ + 2πk/3))²
# Problem: cosh is always > 1, so (1 + √2 cosh) > 1 + √2 > 2.
# All three masses would be large and similar. This DOES give near-degeneracy.

print("\n--- Approach 1: cos → cosh ---")
for th_name, th in [("2/9", mpf(2)/9), ("7/60", mpf(7)/60), ("1/6", mpf(1)/6)]:
    x = [(1 + sqrt(2)*cosh(th + 2*pi*k/3))**2 for k in range(3)]
    x_sorted = sorted(x)
    r32 = x_sorted[2]/x_sorted[1]
    r21 = x_sorted[1]/x_sorted[0]
    print(f"  θ = {th_name}: x = ({nstr(x_sorted[0],4)}, {nstr(x_sorted[1],4)}, {nstr(x_sorted[2],4)}) r32={nstr(r32,4)} r21={nstr(r21,4)}")

# Approach 2: Replace circular phase with exponential decay
# The substrate doesn't rotate; it contracts.
# Instead of e^{2πik/3} (rotation), use e^{-2πk/3} (decay).
# m_k = M² (1 + √2 × e^{-(θ + 2πk/3)})²
# This gives exponentially decreasing masses.

print("\n--- Approach 2: Exponential decay ---")
for th_name, th in [("2/9", mpf(2)/9), ("7/60", mpf(7)/60)]:
    x = [(1 + sqrt(2)*exp(-(th + 2*pi*k/3)))**2 for k in range(3)]
    x_sorted = sorted(x, reverse=True)  # decreasing
    r = x_sorted[0]/x_sorted[1]
    r2 = x_sorted[1]/x_sorted[2]
    print(f"  θ = {th_name}: x = ({nstr(x_sorted[0],6)}, {nstr(x_sorted[1],6)}, {nstr(x_sorted[2],6)}) r12={nstr(r,4)} r23={nstr(r2,4)}")

# Approach 3: The substrate WITNESSES the rotation via projection
# The Z₃ rotation has eigenvalues 1, ω, ω² where ω = e^{2πi/3}.
# The substrate (Z₃ singlet) couples to the sections via the seesaw.
# The coupling projects the Z₃ rotation onto the trivial representation.
# The projected amplitude for generation k is |1 + ω^k|² / normalization
# or Re(coupling_k).

print("\n--- Approach 3: Z₃ projection onto singlet ---")
# If the sections have phases ω^k, the substrate sees |⟨θ|s_k⟩|²
# The overlap ⟨θ|s_k⟩ depends on the Koide angle θ.
# For the charged leptons: ⟨s_j|s_k⟩ ∝ cos(θ + 2πk/3)
# For the substrate: ⟨θ|s_k⟩ ∝ ? (the seesaw projection)

# The seesaw: m_ν = v² Y^T M_R^{-1} Y where Y is the Yukawa matrix.
# In the Z₃ basis, Y_ik couples generation i of the left-handed ν
# to generation k of the right-handed ν.
# If the right-handed ν is a Z₃ singlet, Y has the form Y_k = y × δ_{k,0}
# (only couples to the invariant combination of sections).

# But that gives only ONE massive neutrino! To get three, we need
# three right-handed neutrinos, each coupling differently.

# In the framework: there are THREE right-handed neutrinos
# (one per generation), each a Z₃ singlet but with different
# seesaw couplings. The couplings are determined by the Koide
# structure of the section space.

# Approach 3a: The seesaw matrix M_ν = v² Y^T M_R^{-1} Y
# where Y is the Koide-structured Yukawa matrix.
# Y_ij = δ_ij × y_i where y_i ∝ √(x_i^{lepton})
# (Yukawa couplings track the charged lepton masses)

print("  Using seesaw with Koide-structured Yukawa couplings:")
theta_l = mpf(2)/9
x_l = [(1 + sqrt(2)*cos(theta_l + 2*pi*k/3))**2 for k in range(3)]
x_l_sorted = sorted(x_l)

# Yukawa couplings proportional to sqrt of charged lepton mass factors
y = [sqrt(x) for x in x_l_sorted]

print(f"  Charged lepton x-factors: {[nstr(x,6) for x in x_l_sorted]}")
print(f"  Yukawa couplings ∝ √x: {[nstr(yi,6) for yi in y]}")

# Diagonal seesaw: m_ν_k ∝ y_k² / M_R = x_k^{lepton} / M_R
# This gives m_ν_k ∝ x_k^{lepton} → same hierarchy as charged leptons!
# That's too steep.

# But if M_R is generation-dependent: M_R_k
# The right-handed neutrino in generation k has mass M_R_k
# Then m_ν_k ∝ y_k² / M_R_k = x_k / M_R_k
# For flat hierarchy: need M_R_k ∝ x_k (larger x → larger M_R → less seesaw)

# Approach 4: The contraction formula
# If the substrate sees a CONTRACTION instead of rotation,
# the mass formula might be:
# m_k = M² (1 + r × f(θ_ν + 2πk/3))²
# where f is a contractive function.

# What if f(x) = cos(x) but with a SMALLER r?
# For charged leptons: r = √2
# For neutrinos: r < √2 (contraction reduces the amplitude)

print("\n--- Approach 4: Reduced amplitude r < √2 ---")
print("  m_k = M² (1 + r cos(θ + 2πk/3))² with variable r")
print()

# Constraint: Δm²₂₁/Δm²₃₂ = 1/33
target_ratio = mpf(1)/33

# For given θ and r, compute the splitting ratio
def splitting_ratio(theta, r):
    x = sorted([(1 + r*cos(theta + 2*pi*k/3))**2 for k in range(3)])
    # x[0] < x[1] < x[2]
    # Δm²₂₁ ∝ x[1] - x[0] (if m² ∝ x)...
    # Actually m_k = M² × x_k, so m_k² = M⁴ × x_k²
    # Δm²₂₁ = M⁴(x₂² - x₁²), Δm²₃₂ = M⁴(x₃² - x₂²)
    dm21 = x[1]**2 - x[0]**2
    dm32 = x[2]**2 - x[1]**2
    if float(dm32) < 1e-20:
        return mpf(0)
    return dm21/dm32

# Scan (θ, r) space
print(f"  Scanning for Δm²₂₁/Δm²₃₂ = 1/33 = {nstr(target_ratio, 6)}:")
print(f"  {'θ':>8s} {'r':>8s} {'ratio':>10s} {'x₃/x₂':>8s} {'x₂/x₁':>8s}")

best = (100, 0, 0)
for i_th in range(1, 200):
    for i_r in range(1, 200):
        th = mpf(i_th)/200 * pi/3
        r = mpf(i_r)/200 * sqrt(2)  # r from 0 to √2
        ratio = splitting_ratio(th, r)
        if float(ratio) > 0:
            err = float(fabs(ratio - target_ratio))
            if err < best[0]:
                x = sorted([(1 + r*cos(th + 2*pi*k/3))**2 for k in range(3)])
                best = (err, th, r, x, ratio)

err, th, r, x, ratio = best
print(f"  {nstr(th,6):>8s} {nstr(r,6):>8s} {nstr(ratio,6):>10s} {nstr(x[2]/x[1],4):>8s} {nstr(x[1]/x[0],4):>8s}")

print(f"\n  Best fit: θ = {nstr(th, 10)}, r = {nstr(r, 10)}")
print(f"  Δm²₂₁/Δm²₃₂ = {nstr(ratio, 8)} (target: {nstr(target_ratio, 8)})")

# What are θ and r in framework terms?
print(f"\n  θ = {nstr(th, 10)} vs framework:")
for name, val in [("2/9", mpf(2)/9), ("π/4-2/9", pi/4-mpf(2)/9), ("7/60", mpf(7)/60), ("1/3", mpf(1)/3), ("π/9", pi/9)]:
    err_pct = float(fabs(th-val)/th*100)
    if err_pct < 20:
        print(f"    {name:15s} = {nstr(val,8):>10s}  err={err_pct:.1f}%")

print(f"\n  r = {nstr(r, 10)} vs framework:")
for name, val in [("√2", sqrt(2)), ("1", mpf(1)), ("√(2/3)", sqrt(mpf(2)/3)),
                   ("1/√2", 1/sqrt(2)), ("K*√2", K*sqrt(2)), ("(H-1)/H×√2", mpf(H-1)/H*sqrt(2)),
                   ("√(2(H-1)/H)", sqrt(2*mpf(H-1)/H)), ("√(2/H)", sqrt(mpf(2)/H))]:
    err_pct = float(fabs(r-val)/r*100)
    if err_pct < 20:
        print(f"    {name:20s} = {nstr(val,8):>10s}  err={err_pct:.1f}%")

# ============================================================
# THE KEY INSIGHT: Q changes, not just θ
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 5: DIFFERENT Q FOR NEUTRINOS")
print("=" * 70)

# For charged leptons: Q = (H-1)/H = 2/3 (proved)
# The Koide formula: m_k = M²(1 + r cos(θ + 2πk/3))²
# with r = √(6Q-2) = √2 when Q = 2/3.
#
# What if neutrinos have a DIFFERENT Q?
# The substrate is a Z₃ singlet. Its Koide ratio might differ
# from the sections' 2/3.
#
# For quark sectors: Q depends on colour charge C and isospin T₃:
# Q = (H-1)/H + C/(H²-1) + T₃×K*/2
#
# For neutrinos: C=0 (singlet), T₃ = +1/2 (weak isospin up)
# Q_ν = (H-1)/H + 0 + (1/2)×K*/2 = 2/3 + K*/4 = 2/3 + 7/120 = 87/120 = 29/40

Q_nu = mpf(H-1)/H + K/4
r_nu = sqrt(6*Q_nu - 2)

print(f"  Neutrino Q = (H-1)/H + K*/4 = 2/3 + 7/120 = {nstr(Q_nu, 10)}")
print(f"  r = √(6Q-2) = √(6×{nstr(Q_nu,6)}-2) = {nstr(r_nu, 10)}")
print(f"  r/√2 = {nstr(r_nu/sqrt(2), 10)}")

# Now scan θ with this r
print(f"\n  With Q_ν = {nstr(Q_nu, 6)}, r = {nstr(r_nu, 6)}:")
print(f"  Scanning θ for Δm²₂₁/Δm²₃₂ = 1/33:")

best2 = (100, 0)
for i in range(1, 10000):
    th = mpf(i)/10000 * pi/2
    ratio = splitting_ratio(th, r_nu)
    if float(ratio) > 0:
        err = float(fabs(ratio - target_ratio))
        if err < best2[0]:
            best2 = (err, th)

err2, th2 = best2
x2 = sorted([(1 + r_nu*cos(th2 + 2*pi*k/3))**2 for k in range(3)])
ratio2 = splitting_ratio(th2, r_nu)

print(f"  Best θ = {nstr(th2, 15)}")
print(f"  Δm²₂₁/Δm²₃₂ = {nstr(ratio2, 8)} (target: 1/33 = {nstr(target_ratio, 8)})")
print(f"  x = ({nstr(x2[0],6)}, {nstr(x2[1],6)}, {nstr(x2[2],6)})")
print(f"  x₃/x₂ = {nstr(x2[2]/x2[1], 6)}, x₂/x₁ = {nstr(x2[1]/x2[0], 6)}")

print(f"\n  θ vs framework:")
for name, val in [("7/60", mpf(7)/60), ("2/9", mpf(2)/9), ("K*/(H-1)", K/(H-1)),
                   ("π/4-2/9", pi/4-mpf(2)/9), ("K*/H", K/H), ("1/9", mpf(1)/9),
                   ("K*²", K**2), ("2/27", mpf(2)/27)]:
    err_pct = float(fabs(th2-val)/th2*100)
    if err_pct < 15:
        print(f"    {name:15s} = {nstr(val,10):>12s}  err={err_pct:.2f}%")
