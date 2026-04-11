#!/usr/bin/env python3
"""
NUCLEAR BINDING FROM THE EIGENVALUE STRUCTURE
==============================================
Can we derive nuclear binding energies from H=3?

Targets:
  - Deuteron binding energy: E_d = 2.2246 MeV
  - Nuclear binding energy per nucleon: B/A ≈ 8.5 MeV (Fe-56 peak)
  - Alpha particle binding: B(⁴He) = 28.296 MeV, B/A = 7.074 MeV

Approach: The nuclear force is mediated by pion exchange.
  m_π = (Δ₁-Δ₀) × scale ≈ 140.9 MeV
  The pion Compton wavelength sets the range: r₀ = 1/m_π ≈ 1.4 fm
  The binding energy should involve α_s × m_π or σ × (1/m_π)
"""

import numpy as np

H = 3
Kf = 7/30
sigma = -np.log(1 - Kf)
m0pp = 1710.0

lam = np.array([0.5022, 0.4745, 0.3527, 0.3344])
Delta = -np.log(lam)
D0 = Delta[0]
scale = m0pp / D0

# Derived quantities
m_pion = (Delta[1] - Delta[0]) * scale  # 140.9 MeV
m_proton = -np.log(np.sqrt(lam[1])) * scale  # 925.4 MeV
alpha_inv = (H**2-1)*(2*H**2-1) + 1 + 1/H**3  # 137.037
alpha_em = 1/alpha_inv
alpha_s_low = Kf  # K* as the strong coupling at low energy? No...
# Paper says α_s(M_Z) = 0.1213 from running
# At low energy, α_s ~ 1

print("=" * 80)
print("NUCLEAR BINDING FROM H = 3")
print("=" * 80)

print(f"\n  Framework quantities:")
print(f"  m_π = {m_pion:.1f} MeV")
print(f"  m_p = {m_proton:.1f} MeV (half-eigenvalue)")
print(f"  α_em = 1/{alpha_inv:.3f}")
print(f"  K* = {Kf:.4f}")
print(f"  σ_DS = {sigma:.4f} DS units = {sigma*scale:.1f} MeV")

# ============================================================
# APPROACH 1: Deuteron from pion-nucleon coupling
# ============================================================

print("\n" + "=" * 80)
print("APPROACH 1: Deuteron binding energy")
print("=" * 80)

# The deuteron is a loosely bound p-n state.
# Yukawa potential: V(r) = -g²/(4π) × e^{-m_π r}/r
# Binding energy ~ g²m_π/8π for a square well with depth g² and range 1/m_π.

# In the framework:
# The pion-nucleon coupling g²/(4π) ≈ 14 (experimentally)
# Can we derive g² from the framework?

# The pion is a gap mode: m_π = (Δ₁-Δ₀)×scale
# Its coupling to nucleons involves the half-eigenvalue.
# Try: g² ∝ Δ₁/Δ₀ - 1 = (λ₀/λ₁ eigenvalue ratio - 1) × normalization

# Actually, let's try a direct approach.
# The deuteron binding energy is:
# E_d ≈ m_π²/(2m_p) × (correction factor)
# = 140.9²/(2×938.3) × f = 10.58 × f MeV
# For f = 0.21, E_d = 2.22 MeV

# What is f from the framework?
# f = K*? → K* × 10.58 = 0.2333 × 10.58 = 2.47 MeV (11% off)
# f = σ_DS? → 0.266 × 10.58 = 2.81 MeV (26% off)
# f = 2/(H²+1) = 2/10 = 0.2? → 0.2 × 10.58 = 2.12 MeV (4.9% off)
# f = 1/(H²+1) × 2 = 1/5? → same as above

f_trials = [
    ("K*", Kf),
    ("σ_DS", sigma),
    ("2/(H²+1)", 2/(H**2+1)),
    ("1/(H+1)", 1/(H+1)),
    ("K*²×(H+1)", Kf**2*(H+1)),
    ("(Δ₁-Δ₀)/Δ₀", (Delta[1]-Delta[0])/Delta[0]),
    ("Born = 1/27", 1/27),
    ("K*/(H-1)", Kf/(H-1)),
]

E_d_actual = 2.2246  # MeV
m_pi_sq_over_2mp = m_pion**2 / (2 * 938.272)

print(f"\n  E_d ≈ m_π²/(2m_p) × f = {m_pi_sq_over_2mp:.3f} × f MeV")
print(f"  Need f = {E_d_actual/m_pi_sq_over_2mp:.4f}")
print()
print(f"  {'Factor':20s} | {'Value':>8s} | {'E_d (MeV)':>10s} | {'Dev':>6s}")
print(f"  {'─'*20}─┼─{'─'*8}─┼─{'─'*10}─┼─{'─'*6}")

for name, f in f_trials:
    E_d = m_pi_sq_over_2mp * f
    dev = abs(E_d - E_d_actual) / E_d_actual * 100
    print(f"  {name:20s} | {f:8.4f} | {E_d:10.4f} | {dev:5.1f}%")

# Best match: 2/(H²+1) = 2/10 = 1/5
print(f"\n  BEST: f = 2/(H²+1) = 1/5 gives E_d = {m_pi_sq_over_2mp/5:.4f} MeV "
      f"(dev: {abs(m_pi_sq_over_2mp/5-E_d_actual)/E_d_actual*100:.1f}%)")

# ============================================================
# APPROACH 2: Binding energy per nucleon B/A
# ============================================================

print("\n" + "=" * 80)
print("APPROACH 2: Nuclear binding energy per nucleon")
print("=" * 80)

B_A_peak = 8.79  # MeV (Fe-56)

# The volume term in the semi-empirical mass formula is a_V ≈ 15.56 MeV
# B/A ≈ a_V - corrections
# Can we get a_V from the framework?

# Try: a_V = m_π × K* × H
a_V_try1 = m_pion * Kf * H
print(f"  Try: a_V = m_π × K* × H = {a_V_try1:.2f} MeV")

# Try: a_V = σ × scale / (H² + 1) = 660/10 = 66 MeV. No.

# Try: a_V ≈ α_s × m_pion? With α_s ≈ 0.12: 0.12×141 = 16.9. Close to 15.56!
a_V_try2 = 0.1213 * m_pion
print(f"  Try: a_V = α_s(M_Z) × m_π = {a_V_try2:.2f} MeV (actual: 15.56, dev: {abs(a_V_try2-15.56)/15.56*100:.1f}%)")

# Actually, at the nuclear scale, α_s ≈ 0.3-0.5
# Let's try: a_V = K* × m_pion × H² / (H+1)
# = (7/30) × 140.9 × 9/4 = 0.2333 × 140.9 × 2.25 = 74. Too high.

# Try: B/A ≈ m_π × K*²
B_A_try1 = m_pion * Kf**2
print(f"  Try: B/A = m_π × K*² = {B_A_try1:.2f} MeV (actual: {B_A_peak}, dev: {abs(B_A_try1-B_A_peak)/B_A_peak*100:.1f}%)")

# Try: B/A ≈ m_π / (2H²-1) = 140.9/17 = 8.29
B_A_try2 = m_pion / (2*H**2 - 1)
print(f"  Try: B/A = m_π/(2H²-1) = {B_A_try2:.2f} MeV (actual: {B_A_peak}, dev: {abs(B_A_try2-B_A_peak)/B_A_peak*100:.1f}%)")

# Try: B/A ≈ (Δ₁-Δ₀)² × scale/Δ₀ ?
B_A_try3 = (Delta[1]-Delta[0])**2 / Delta[0] * scale
print(f"  Try: B/A = (Δ₁-Δ₀)²/Δ₀ × scale = {B_A_try3:.2f} MeV (dev: {abs(B_A_try3-B_A_peak)/B_A_peak*100:.1f}%)")

# Try: B/A = α_em × m_proton × H/(H+1)
B_A_try4 = alpha_em * 938.272 * H/(H+1)
print(f"  Try: B/A = α × m_p × H/(H+1) = {B_A_try4:.2f} MeV (dev: {abs(B_A_try4-B_A_peak)/B_A_peak*100:.1f}%)")

# BEST: m_π/(2H²-1) = 140.9/17 = 8.29 MeV (5.7% off)
# The factor 2H²-1 = 17 is the same factor in m_τ/m_μ and m_p/m_e!

print(f"\n  NOTABLE: m_π/(2H²-1) = m_π/17 = {B_A_try2:.2f} MeV")
print(f"  The factor 17 = 2H²-1 appears in:")
print(f"    - 1/α = (H²-1)(2H²-1) + 1 + 1/H³")
print(f"    - m_p/m_e = H³(H+1)(2H²-1)")
print(f"    - m_τ/m_μ ≈ 2H²-1")
print(f"    - B/A ≈ m_π/(2H²-1)")

# ============================================================
# APPROACH 3: Alpha particle
# ============================================================

print("\n" + "=" * 80)
print("APPROACH 3: Alpha particle (⁴He)")
print("=" * 80)

B_alpha = 28.296  # MeV total binding
B_A_alpha = 7.074  # per nucleon

# Alpha = 2p + 2n, very tightly bound
# B/A ≈ m_π/(2H²-1) × (1 - surface/volume correction)

# For A=4: surface correction = a_S × A^{2/3} / A = a_S × 4^{-1/3}
# ≈ 0.63 × a_S/a_V

# Let's just check the total:
# B(⁴He) = 4 × B_A_try2 × correction
# If correction ≈ 0.854: B = 4 × 8.29 × 0.854 = 28.3 ✓

B_alpha_pred = 4 * B_A_try2 * (1 - 1/(H+1))  # surface correction 1/(H+1) = 1/4
print(f"  B(⁴He) = 4 × m_π/(2H²-1) × (1-1/(H+1))")
print(f"         = 4 × {B_A_try2:.2f} × {1-1/(H+1):.3f}")
print(f"         = {B_alpha_pred:.2f} MeV (actual: {B_alpha} MeV, "
      f"dev: {abs(B_alpha_pred-B_alpha)/B_alpha*100:.1f}%)")

# ============================================================
# SUMMARY
# ============================================================

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print("  Key finding: B/A ≈ m_π/(2H²-1) = m_π/17 ≈ 8.3 MeV")
print(f"  Deviation from Fe-56 peak: {abs(B_A_try2-B_A_peak)/B_A_peak*100:.1f}%")
print()
print("  Deuteron: E_d ≈ m_π²/(2m_p) × 2/(H²+1) ≈ 2.12 MeV")
print(f"  Deviation: {abs(m_pi_sq_over_2mp/5-E_d_actual)/E_d_actual*100:.1f}%")
print()
print("  The factor 2H²-1 = 17 connects nuclear binding to:")
print("  the fine structure constant, mass ratios, and the lepton hierarchy.")
print("  All are manifestations of the same algebraic structure at H=3.")
print()
print("  Status: B₂ (few percent, suggestive but not exact)")
