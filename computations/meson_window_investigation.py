#!/usr/bin/env python3
"""
MESON SPECTRUM IN THE 400-1020 MeV WINDOW
==========================================
The paper identifies 5 unmatched hadrons:
  K± (493.7), K⁰ (497.6), ρ (775.3), ω (782.7), φ (1019.5)

These sit where strange-quark content and chiral dynamics
are expected to shift energies from the leading-order pattern.

Can we find eigenvalue combinations that work?
Or identify the structural reason for the mismatch?
"""

import numpy as np
from itertools import product

H = 3
Kf = 7/30
sigma = -np.log(1 - Kf)  # 0.2657
m0pp = 1710.0

# A2 eigenvalues
lam = np.array([0.5022, 0.4745, 0.3527, 0.3344])
Delta = -np.log(lam)
D0 = Delta[0]
scale = m0pp / D0

print("=" * 80)
print("MESON WINDOW: 400-1020 MeV")
print("=" * 80)

# Framework building blocks
print("\n  Building blocks:")
print(f"  Eigenvalue gaps: Δ₀={Delta[0]:.4f}, Δ₁={Delta[1]:.4f}, "
      f"Δ₂={Delta[2]:.4f}, Δ₃={Delta[3]:.4f}")
print(f"  String tension:  σ = {sigma:.4f}")
print(f"  Scale:           {scale:.1f} MeV per DS unit")
print(f"  Pion mass:       (Δ₁-Δ₀)×scale = {(Delta[1]-Delta[0])*scale:.1f} MeV")

# Target states
targets = [
    ("K±", 493.7, 0.1),
    ("K⁰", 497.6, 0.1),
    ("ρ(770)", 775.3, 0.3),
    ("ω(782)", 782.7, 0.2),
    ("φ(1020)", 1019.5, 0.5),
    ("η (comparison)", 547.9, 1.0),
]

print(f"\n  Target masses:")
for name, mass, unc in targets:
    ds_units = mass / scale
    print(f"    {name:15s}: {mass:8.1f} MeV = {ds_units:.4f} DS units")

# ============================================================
# SEARCH: All combinations λᵢ^(p/q) + nσ
# ============================================================

print("\n" + "=" * 80)
print("SYSTEMATIC SEARCH: λᵢ^(p/q) + nσ")
print("=" * 80)

best_matches = {}

for name, target_mass, _ in targets:
    target_ds = target_mass / scale
    best = None
    best_dev = float('inf')

    for i in range(4):  # eigenvalue index
        for p in range(1, 13):  # numerator
            for q in [1, 2, 3, 4, 6, 12]:  # denominator | 12
                if p > q * 4:  # reasonable range
                    continue
                power = p / q
                for n in range(-2, 5):  # string tension units
                    mass_ds = -np.log(lam[i]**power) + n * sigma
                    if mass_ds <= 0:
                        continue
                    mass_mev = mass_ds * scale
                    dev = abs(mass_mev - target_mass) / target_mass * 100

                    if dev < best_dev:
                        best_dev = dev
                        best = (i, p, q, n, mass_mev, dev)

    # Also try gap differences
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            for n in range(-2, 5):
                mass_ds = abs(Delta[i] - Delta[j]) + n * sigma
                mass_mev = mass_ds * scale
                dev = abs(mass_mev - target_mass) / target_mass * 100
                if dev < best_dev:
                    best_dev = dev
                    best = (-1, i, j, n, mass_mev, dev)  # -1 = gap mode

    # Try Koide-type combinations
    # m_K ~ M²×x (Koide mass parameter for strange sector)
    theta_down = 1/9  # conjectured down-type Koide angle
    x_down = [(1 + np.sqrt(2)*np.cos(theta_down + 2*np.pi*k/3))**2 for k in range(3)]
    M2_down = m0pp * 3/8  # = sin²θ_W × m(0++)
    for k, xk in enumerate(x_down):
        mass_k = M2_down * xk
        dev = abs(mass_k - target_mass) / target_mass * 100
        if dev < best_dev:
            best_dev = dev
            best = (-2, k, 0, 0, mass_k, dev)  # -2 = Koide mode

    best_matches[name] = best

print(f"\n  {'State':15s} | {'Target':>8s} | {'Best match':>10s} | {'Dev':>6s} | Formula")
print(f"  {'─'*15}─┼─{'─'*8}─┼─{'─'*10}─┼─{'─'*6}─┼─{'─'*40}")

for name, target_mass, _ in targets:
    b = best_matches[name]
    if b is None:
        print(f"  {name:15s} | {target_mass:8.1f} | {'---':>10s} | {'---':>6s} |")
        continue

    if b[0] == -1:
        # Gap mode
        _, i, j, n, mass, dev = b
        sigma_str = f"+{n}σ" if n > 0 else (f"{n}σ" if n < 0 else "")
        formula = f"|Δ_{i}-Δ_{j}|{sigma_str}"
    elif b[0] == -2:
        # Koide mode
        _, k, _, _, mass, dev = b
        formula = f"Koide down k={k}: M²_d×x_{k}"
    else:
        i, p, q, n, mass, dev = b
        sigma_str = f"+{n}σ" if n > 0 else (f"{n}σ" if n < 0 else "")
        if q == 1:
            formula = f"λ_{i}^{p}{sigma_str}"
        else:
            formula = f"λ_{i}^({p}/{q}){sigma_str}"

    print(f"  {name:15s} | {target_mass:8.1f} | {mass:10.1f} | {dev:5.1f}% | {formula}")

# ============================================================
# PHYSICAL ANALYSIS
# ============================================================

print("\n" + "=" * 80)
print("PHYSICAL ANALYSIS")
print("=" * 80)

# The kaon
m_pion = (Delta[1] - Delta[0]) * scale
print(f"\n  KAON ANALYSIS:")
print(f"  m_π = {m_pion:.1f} MeV (gap mode: Δ₁-Δ₀)")
print(f"  m_K = {493.7:.1f} MeV")
print(f"  m_K/m_π = {493.7/m_pion:.3f}")
print(f"  m_K²/m_π² = {493.7**2/m_pion**2:.3f}")

# Gell-Mann-Okubo: m_K² ≈ ½(m_π² + m_η²) for pseudoscalar octet
# m_η_pred = 547.9
# m_K² ≈ ½(139.6² + 547.9²) ≈ ½(19488 + 300194) = 159841 → m_K ≈ 400. Too low.
# The GMO formula isn't perfect either.

# The kaon mass in the eigenvalue framework:
# K contains strange quark → the Koide down-type sector
# m_s ≈ 93.5 MeV (from Koide with θ=1/9)
# m_K² ≈ m_s × Λ_QCD (PCAC relation)
# Λ_QCD ≈ m0pp/10 ≈ 171 MeV from the framework?

# Actually, let's try: m_K ≈ m_π + m_s
# = 140 + 93.5 = 233.5. No, too low.

# m_K² ≈ m_π² + 2m_s × Λ (chiral perturbation theory)
# This requires Λ, which we don't have cleanly.

# Let's try a different approach: m_K from eigenvalue fractions
# m_K ≈ λ₀^(1/6) × scale ?
m_K_try = -np.log(lam[0]**(1/6)) * scale
print(f"\n  Try: m_K = -ln(λ₀^(1/6))×scale = {m_K_try:.1f} MeV (dev: {abs(m_K_try-493.7)/493.7*100:.1f}%)")

# m_K ≈ Δ₁/3 × scale ?
m_K_try2 = Delta[1]/3 * scale
print(f"  Try: m_K = Δ₁/3 × scale = {m_K_try2:.1f} MeV (dev: {abs(m_K_try2-493.7)/493.7*100:.1f}%)")

# m_K ≈ (Δ₁-Δ₀+σ) × scale?
m_K_try3 = (Delta[1] - Delta[0] + sigma) * scale
print(f"  Try: m_K = (Δ₁-Δ₀+σ)×scale = {m_K_try3:.1f} MeV (dev: {abs(m_K_try3-493.7)/493.7*100:.1f}%)")

# m_K ≈ σ×2 × scale (2 string tension units)?
m_K_try4 = 2 * sigma * scale
print(f"  Try: m_K = 2σ×scale = {m_K_try4:.1f} MeV (dev: {abs(m_K_try4-493.7)/493.7*100:.1f}%)")

# The ρ meson
print(f"\n  ρ(770) ANALYSIS:")
print(f"  m_ρ = 775.3 MeV")
# ρ is the lightest vector meson (J^PC = 1--)
# m_ρ ≈ 2 × m_constituent_quark ≈ 2 × 330 = 660 MeV? No, too low.
# m_ρ ≈ σ × 3?
m_rho_try = 3 * sigma * scale
print(f"  Try: m_ρ = 3σ×scale = {m_rho_try:.1f} MeV (dev: {abs(m_rho_try-775.3)/775.3*100:.1f}%)")

# m_ρ ≈ ½Δ₁ × scale ?
m_rho_try2 = Delta[1]/2 * scale
print(f"  Try: m_ρ = Δ₁/2 × scale = {m_rho_try2:.1f} MeV (dev: {abs(m_rho_try2-775.3)/775.3*100:.1f}%)")

# m_ρ ≈ λ₁^(1/3) ?
m_rho_try3 = -np.log(lam[1]**(1/3)) * scale
print(f"  Try: m_ρ = -ln(λ₁^(1/3))×scale = {m_rho_try3:.1f} MeV (dev: {abs(m_rho_try3-775.3)/775.3*100:.1f}%)")

# The φ meson
print(f"\n  φ(1020) ANALYSIS:")
print(f"  m_φ = 1019.5 MeV (ss̄ state)")
# φ is almost pure ss̄
# m_φ ≈ 2m_s(constituent) ≈ 2 × 500 = 1000 MeV (close!)
# In the framework: m_φ ≈ ?

# m_φ ≈ Δ₀ × 2/3 × scale ?
m_phi_try = Delta[0] * 2/3 * scale
print(f"  Try: m_φ = 2Δ₀/3 × scale = {m_phi_try:.1f} MeV (dev: {abs(m_phi_try-1019.5)/1019.5*100:.1f}%)")

# m_φ ≈ λ₀^(1/3) + σ ?
m_phi_try2 = (-np.log(lam[0]**(1/3)) + sigma) * scale
print(f"  Try: m_φ = -ln(λ₀^(1/3))+σ = {m_phi_try2:.1f} MeV (dev: {abs(m_phi_try2-1019.5)/1019.5*100:.1f}%)")

# KSRF relation: m_ρ² = 2f_π²g_ρππ² where f_π ≈ 93 MeV
# From the framework: f_π ≈ m_π/K* ?
f_pi_try = m_pion / (np.sqrt(2) * Kf)
print(f"\n  KSRF relation:")
print(f"  f_π (framework) ≈ m_π/(√2 K*) = {f_pi_try:.1f} MeV (actual: 93 MeV, dev: {abs(f_pi_try-93)/93*100:.1f}%)")

# Try f_π = m_pion / √(2×H) = 140.9/√6
f_pi_try2 = m_pion / np.sqrt(2 * H)
print(f"  f_π (alt) ≈ m_π/√(2H) = {f_pi_try2:.1f} MeV (dev: {abs(f_pi_try2-93)/93*100:.1f}%)")

# Summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print("  The 400-1020 MeV window is where the leading-order eigenvalue")
print("  pattern breaks down. The reason is structural:")
print()
print("  1. KAONS: m_K involves the strange quark, which lives at the")
print("     interference node of the Koide pattern. The eigenvalue")
print("     spectroscopy gives band states; kaons require the chiral")
print("     dynamics that mix band and interference sectors.")
print()
print("  2. VECTOR MESONS (ρ, ω, φ): These are spin-1 qq̄ states.")
print("     The eigenvalue spectroscopy describes colour-singlet glue.")
print("     Vector mesons require the quark propagator (half-eigenvalue)")
print("     combined with the anti-quark propagator — a PRODUCT of")
print("     two half-eigenvalue modes, not a single eigenvalue power.")
print()
print("  3. The φ meson is ss̄ — two strange quarks. Its mass is")
print("     approximately 2×m_s(constituent) ≈ 1000 MeV, which is")
print("     a chiral-dynamics quantity, not a glue quantity.")
print()
print("  CONCLUSION: The 400-1020 MeV window requires the QUARK")
print("  propagator (half-eigenvalue × Koide), not just the glue")
print("  spectroscopy. This is consistent with the paper's identification")
print("  of the gap as a 'chiral dynamics' window.")
