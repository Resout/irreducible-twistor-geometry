#!/usr/bin/env python3
"""
Test the proposed correction: Δθ = -5/(39²) × α²

Where:
  39 = H × (H² + H + 1) = 3 × 13
  5 = dim(triv)
  α = fine structure constant

Also expressible as:
  Δθ = -dim(triv)³ / 195² × α²
  since 5³/195² = 125/38025 = 5/1521

This is O(α²) — exactly the Sumino 2-loop residual.
"""
import numpy as np

# Constants
H = 3
alpha = 7.2973525693e-3  # CODATA 2018 fine structure constant
m_proton = 938.272088     # MeV

# Experimental masses
m_e_exp = 0.51099895000
m_mu_exp = 105.6583755
m_tau_exp = 1776.86

# Framework
dim_triv = 5
dim_std = 10
dim_sign = 1
cyclotomic = H**2 + H + 1  # = 13

# The claimed correction
factor_39 = H * cyclotomic  # = 39
factor_1521 = factor_39**2  # = 1521

delta_theta_claimed = -(dim_triv / factor_1521) * alpha**2

print("=" * 70)
print("proposed correction: Δθ = -5/1521 × α²")
print("=" * 70)

print(f"\n  H = {H}")
print(f"  α = {alpha}")
print(f"  α² = {alpha**2:.10e}")
print(f"  dim(triv) = {dim_triv}")
print(f"  H × (H²+H+1) = {H} × {cyclotomic} = {factor_39}")
print(f"  39² = {factor_1521}")
print(f"  5/1521 = {dim_triv/factor_1521:.10e}")
print(f"  195 = 15 × 13 = {15 * 13}")
print(f"  5³/195² = {dim_triv**3}/{195**2} = {dim_triv**3/195**2:.10e}")
print(f"  5/1521 = {5/1521:.10e}  ← same")

print(f"\n  Δθ_claimed = -5/1521 × α² = {delta_theta_claimed:.10e}")

# Target from Brannen fit
theta_theory = 2.0 / 9.0
theta_brannen = 0.22222204717
delta_theta_target = theta_brannen - theta_theory

print(f"\n  Δθ_target (Brannen) = {delta_theta_target:.10e}")
print(f"  Δθ_claimed          = {delta_theta_claimed:.10e}")
print(f"  Difference          = {delta_theta_claimed - delta_theta_target:.4e}")
print(f"  Ratio               = {delta_theta_claimed / delta_theta_target:.8f}")

# Now compute lepton masses with corrected theta
print(f"\n{'='*70}")
print("LEPTON MASSES WITH CORRECTED θ")
print("=" * 70)

phases = np.array([0, 2*np.pi/3, 4*np.pi/3])

def koide_masses(theta, M):
    x = (1 + np.sqrt(2) * np.cos(theta + phases))**2
    return M * x

# θ corrected
theta_corrected = theta_theory + delta_theta_claimed
print(f"\n  θ = 2/9 + Δθ = {theta_corrected:.15f}")
print(f"  θ_Brannen     = {theta_brannen:.15f}")
print(f"  θ_residual    = {theta_corrected - theta_brannen:.4e}")

# M from 195 formula
x_mu_corr = (1 + np.sqrt(2) * np.cos(theta_corrected + 4*np.pi/3))**2
M_195 = m_proton / (H * (1 - (H-1) * x_mu_corr / 195))

masses_corr = koide_masses(theta_corrected, M_195)

print(f"\n  M(θ_corrected) = {M_195:.6f} MeV")
print(f"\n  {'Particle':>8}  {'Predicted':>14}  {'Experimental':>14}  {'Δ (eV)':>12}  {'σ':>8}")
print(f"  {'--------':>8}  {'----------':>14}  {'------------':>14}  {'------':>12}  {'---':>8}")

exp = [m_tau_exp, m_e_exp, m_mu_exp]
sig = [0.12, 0.00000000015, 0.0000023]
names = ['tau', 'electron', 'muon']
for i in range(3):
    delta_eV = (masses_corr[i] - exp[i]) * 1e6  # MeV to eV
    sigma = delta_eV / (sig[i] * 1e6) if sig[i] > 0 else 0
    print(f"  {names[i]:>8}  {masses_corr[i]:14.8f}  {exp[i]:14.8f}  {delta_eV:+12.2f}  {sigma:+8.1f}")

# Now with BOTH M and θ free, vs M_195 + corrected θ
print(f"\n{'='*70}")
print("COMPARISON: corrected vs best-fit")
print("=" * 70)

# Best fit (from our earlier computation)
from scipy.optimize import minimize

def chi2_both(params):
    theta, M = params
    pred = koide_masses(theta, M)
    sigmas = np.array([0.12, 0.00000000015, 0.0000023])
    return np.sum(((pred - np.array(exp)) / sigmas)**2)

res = minimize(chi2_both, [theta_theory, M_195], method='Nelder-Mead',
               options={'xatol': 1e-16, 'fatol': 1e-10, 'maxiter': 100000})
theta_best, M_best = res.x

print(f"  θ_corrected = {theta_corrected:.15f}")
print(f"  θ_best_fit  = {theta_best:.15f}")
print(f"  Residual    = {theta_corrected - theta_best:.4e}")
print()
print(f"  M_195(corr) = {M_195:.6f}")
print(f"  M_best_fit  = {M_best:.6f}")
print(f"  ΔM          = {M_195 - M_best:+.6f} MeV  ({(M_195-M_best)/M_best*1e6:+.1f} ppm)")

# The joint correction check
print(f"\n{'='*70}")
print("STRUCTURAL DECOMPOSITION")
print("=" * 70)
print(f"\n  The correction Δθ = -dim(triv)/(H×(H²+H+1))² × α²")
print(f"                    = -5/39² × α²")
print(f"                    = -dim(triv)³/195² × α²")
print(f"\n  195 appears in M formula (linearly): M = m_p/(H×(1-2x_μ/195))")
print(f"  195 appears in θ formula (quadratically): Δθ = -5³/195² × α²")
print(f"\n  This is O(α²) — consistent with Sumino cancellation of O(α)")
print(f"  The radiative correction is dressed by framework geometry (5/39²)")


# What about the JOINT correction to M?
print(f"\n{'='*70}")
print("JOINT M CORRECTION")
print("=" * 70)
print(f"\n  With θ_corrected, the 195 formula gives M = {M_195:.6f}")
print(f"  Best-fit M = {M_best:.6f}")
print(f"  Remaining M gap = {(M_195-M_best)/M_best*1e6:+.1f} ppm")
print(f"\n  Does the α² correction also fix M? Need M_corrected ≈ M_best")

# Try: what if M also gets an α² correction?
# M_corr = m_p / (H × (1 - 2x_μ/195)) × (1 + c × α²)
# What c would close the gap?
if abs(M_best - M_195) > 1e-10:
    c_needed = (M_best/M_195 - 1) / alpha**2
    print(f"  To close M gap: need (1 + c×α²) with c = {c_needed:.4f}")
    print(f"  That's c = {c_needed:.4f} ≈ ?")
    # Check some framework values
    candidates = {
        'dim(triv)/H': dim_triv/H,
        '1/K*': 30/7,
        'H²+1': H**2+1,
        'dim(std)/H': dim_std/H,
        '195/H²': 195/H**2,
        'cyclotomic': cyclotomic,
        '(H²+H+1)/H': cyclotomic/H,
        'dim(std+triv)/H': 15/H,
    }
    print(f"\n  Candidate framework values for c:")
    for name, val in sorted(candidates.items(), key=lambda x: abs(x[1]-c_needed)):
        print(f"    {name:>20} = {val:8.4f}  (ratio to needed: {val/c_needed:.4f})")
