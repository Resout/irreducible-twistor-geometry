#!/usr/bin/env python3
"""
Test joint (M, θ) closure with both α² corrections:
  θ = 2/9 - (5/1521) α²
  M = M_195 × (1 + (23/25) α²)

Does this close all three lepton masses?
"""
import numpy as np

H = 3
alpha = 7.2973525693e-3
m_proton = 938.272088

m_e_exp = 0.51099895000
m_mu_exp = 105.6583755
m_tau_exp = 1776.86

phases = np.array([0, 2*np.pi/3, 4*np.pi/3])

def koide_multipliers(theta):
    return (1 + np.sqrt(2) * np.cos(theta + phases))**2

# Corrected theta
theta_corr = 2.0/9.0 - (5.0/1521.0) * alpha**2

# M from 195 formula
x_mu = koide_multipliers(theta_corr)[2]
M_195 = m_proton / (H * (1 - (H-1) * x_mu / 195))

# M with 23/25 α² correction
M_corr = M_195 * (1 + (23.0/25.0) * alpha**2)

# Masses
masses_corr = M_corr * koide_multipliers(theta_corr)

print("=" * 70)
print("JOINT α² CLOSURE TEST")
print("=" * 70)

print(f"\n  θ = 2/9 - (5/1521)α² = {theta_corr:.15f}")
print(f"  M_195 = {M_195:.6f} MeV")
print(f"  M_corr = M_195 × (1 + 23/25 × α²) = {M_corr:.6f} MeV")
print(f"  ΔM = {M_corr - M_195:.6f} MeV = {(M_corr-M_195)/M_195*1e6:.1f} ppm")

# Best-fit M for comparison
from scipy.optimize import minimize
def chi2_both(params):
    theta, M = params
    pred = M * koide_multipliers(theta)
    sigmas = np.array([0.12, 0.00000000015, 0.0000023])
    return np.sum(((pred - np.array([m_tau_exp, m_e_exp, m_mu_exp])) / sigmas)**2)

res = minimize(chi2_both, [2.0/9.0, M_195], method='Nelder-Mead',
               options={'xatol': 1e-16, 'fatol': 1e-10, 'maxiter': 100000})
theta_best, M_best = res.x

print(f"\n  M_best_fit = {M_best:.6f} MeV")
print(f"  M_corr vs M_best: {(M_corr-M_best)/M_best*1e6:+.1f} ppm")

print(f"\n{'='*70}")
print("LEPTON MASSES")
print("=" * 70)

exp = [m_tau_exp, m_e_exp, m_mu_exp]
sig = [0.12, 0.00000000015, 0.0000023]
names = ['tau', 'electron', 'muon']

print(f"\n  {'Particle':>8}  {'Predicted':>14}  {'Experimental':>14}  {'Δ (eV)':>12}  {'σ':>8}")
print(f"  {'--------':>8}  {'----------':>14}  {'------------':>14}  {'------':>12}  {'---':>8}")
for i in range(3):
    delta_eV = (masses_corr[i] - exp[i]) * 1e6
    sigma = delta_eV / (sig[i] * 1e6) if sig[i] > 0 else 0
    print(f"  {names[i]:>8}  {masses_corr[i]:14.8f}  {exp[i]:14.8f}  {delta_eV:+12.2f}  {sigma:+8.1f}")

# Also test masses at (theta_best, M_best) for comparison
masses_best = M_best * koide_multipliers(theta_best)
print(f"\n  Best-fit comparison:")
for i in range(3):
    delta_eV = (masses_best[i] - exp[i]) * 1e6
    sigma = delta_eV / (sig[i] * 1e6) if sig[i] > 0 else 0
    print(f"  {names[i]:>8}  {masses_best[i]:14.8f}  {exp[i]:14.8f}  {delta_eV:+12.2f}  {sigma:+8.1f}")

# How much of the M gap does 23/25 close?
print(f"\n{'='*70}")
print("GAP ANALYSIS")
print("=" * 70)
M_gap_total = M_best - M_195
M_gap_after = M_best - M_corr
print(f"  M gap before 23/25 correction: {M_gap_total*1e6:.0f} eV = {M_gap_total/M_best*1e6:.1f} ppm")
print(f"  M gap after 23/25 correction:  {M_gap_after*1e6:.0f} eV = {M_gap_after/M_best*1e6:.1f} ppm")
print(f"  Fraction closed: {(1 - abs(M_gap_after/M_gap_total))*100:.1f}%")

# 23/25 decomposition
print(f"\n{'='*70}")
print("STRUCTURAL DECOMPOSITION OF 23/25")
print("=" * 70)
print(f"  23 = 30×(1-K*) = 30×23/30 = surviving information budget")
print(f"  25 = dim(triv)² = 5²")
print(f"  23/25 = {23/25:.4f}")
print(f"  23/25 × α² = {23/25 * alpha**2:.6e}")
print(f"  ΔM/M needed = {M_gap_total/M_195:.6e}")
print(f"  Ratio = {(23/25 * alpha**2) / (M_gap_total/M_195):.6f}")

# What tau mass does this predict?
print(f"\n{'='*70}")
print("TAU MASS PREDICTION (testable)")
print("=" * 70)
print(f"  Predicted m_tau = {masses_corr[0]:.4f} MeV")
print(f"  Current exp:     {m_tau_exp:.4f} ± 0.12 MeV")
print(f"  Deviation:       {masses_corr[0]-m_tau_exp:+.4f} MeV = {(masses_corr[0]-m_tau_exp)/0.12:+.1f}σ")
