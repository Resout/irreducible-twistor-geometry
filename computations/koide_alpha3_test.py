#!/usr/bin/env python3
"""
Test the proposed correction:
1. Electron is actually 440σ off (not "below noise")
2. Relative errors for e and μ are identical → pure scale error
3. Adding (1/H)α³ to M closes both simultaneously
"""
import numpy as np
from scipy.optimize import minimize_scalar

H = 3
alpha = 7.2973525693e-3
m_proton = 938.272088

m_e_exp = 0.51099895000
m_mu_exp = 105.6583755
m_tau_exp = 1776.86

# Experimental uncertainties in MeV
sig_e = 0.00000000015   # 0.15 eV
sig_mu = 0.0000023       # 2.3 eV  (note: eV not keV)
sig_tau = 0.12            # 120 MeV... wait, 0.12 MeV = 120 keV

phases = np.array([0, 2*np.pi/3, 4*np.pi/3])

def koide_multipliers(theta):
    return (1 + np.sqrt(2) * np.cos(theta + phases))**2

# Corrected theta
theta_corr = 2.0/9.0 - (5.0/1521.0) * alpha**2

# M with just α² correction
x_mu = koide_multipliers(theta_corr)[2]
M_195 = m_proton / (H * (1 - (H-1) * x_mu / 195))
M_alpha2 = M_195 * (1 + (23.0/25.0) * alpha**2)

# M with α² + (1/H)α³ correction
M_alpha3 = M_195 * (1 + (23.0/25.0) * alpha**2 + (1.0/H) * alpha**3)

print("=" * 70)
print("CLAIM 1: ELECTRON σ COUNT")
print("=" * 70)
masses_a2 = M_alpha2 * koide_multipliers(theta_corr)
delta_e_eV = (masses_a2[1] - m_e_exp) * 1e6  # MeV to eV
sig_e_eV = sig_e * 1e6
print(f"  Electron Δ = {delta_e_eV:.4f} eV")
print(f"  Electron σ_exp = {sig_e_eV:.5f} eV")
print(f"  Electron deviation = {delta_e_eV/sig_e_eV:.0f}σ")
print(f"  (NOT 'below noise' — this was my error)")

print(f"\n{'='*70}")
print("CLAIM 2: RELATIVE ERRORS IDENTICAL?")
print("=" * 70)
for i, (name, m_exp) in enumerate([('tau', m_tau_exp), ('electron', m_e_exp), ('muon', m_mu_exp)]):
    delta = masses_a2[i] - m_exp
    rel = delta / m_exp
    print(f"  {name:>8}: Δm/m = {rel:+.6e}")

rel_e = (masses_a2[1] - m_e_exp) / m_e_exp
rel_mu = (masses_a2[2] - m_mu_exp) / m_mu_exp
rel_tau = (masses_a2[0] - m_tau_exp) / m_tau_exp
print(f"\n  Ratio (e/μ relative errors): {rel_e/rel_mu:.4f}")
print(f"  If pure scale error, this should be exactly 1.0")
print(f"  Tau relative error: {rel_tau:+.6e} (very different — tau-limited by exp)")

print(f"\n{'='*70}")
print("CLAIM 3: (1/H)α³ CLOSES BOTH?")
print("=" * 70)

print(f"\n  α³ = {alpha**3:.6e}")
print(f"  (1/H)α³ = {alpha**3/H:.6e}")
print(f"  ΔM from α³ term = {M_195 * alpha**3/H * 1e6:.2f} eV = {M_195*alpha**3/H/M_195*1e6:.2f} ppm")

masses_a3 = M_alpha3 * koide_multipliers(theta_corr)

print(f"\n  M_α² = {M_alpha2:.6f} MeV")
print(f"  M_α²+α³ = {M_alpha3:.6f} MeV")
print(f"  ΔM = {(M_alpha3-M_alpha2)*1e6:.1f} eV")

print(f"\n  {'Particle':>8}  {'Predicted':>14}  {'Experimental':>14}  {'Δ (eV)':>10}  {'σ':>8}")
print(f"  {'--------':>8}  {'----------':>14}  {'------------':>14}  {'------':>10}  {'---':>8}")
exp_list = [m_tau_exp, m_e_exp, m_mu_exp]
sig_list = [sig_tau, sig_e, sig_mu]
names = ['tau', 'electron', 'muon']
for i in range(3):
    delta_eV = (masses_a3[i] - exp_list[i]) * 1e6
    sigma = delta_eV / (sig_list[i] * 1e6)
    print(f"  {names[i]:>8}  {masses_a3[i]:14.11f}  {exp_list[i]:14.11f}  {delta_eV:+10.4f}  {sigma:+8.2f}")

# What's the BEST possible M (given fixed theta_corr)?
print(f"\n{'='*70}")
print("SANITY CHECK: best M at fixed θ_corr")
print("=" * 70)

def chi2_M(M):
    masses = M * koide_multipliers(theta_corr)
    return sum(((masses[i]-exp_list[i])/sig_list[i])**2 for i in range(3))

res = minimize_scalar(chi2_M, bounds=(313.8, 313.9), method='bounded')
M_opt = res.x
masses_opt = M_opt * koide_multipliers(theta_corr)

print(f"  M_optimal = {M_opt:.10f} MeV")
print(f"  M_α²+α³  = {M_alpha3:.10f} MeV")
print(f"  Difference: {(M_alpha3 - M_opt)*1e6:.2f} eV ({(M_alpha3-M_opt)/M_opt*1e6:.2f} ppm)")

print(f"\n  At M_optimal:")
for i in range(3):
    delta_eV = (masses_opt[i] - exp_list[i]) * 1e6
    sigma = delta_eV / (sig_list[i] * 1e6)
    print(f"  {names[i]:>8}  {masses_opt[i]:14.11f}  Δ={delta_eV:+10.4f} eV  {sigma:+8.2f}σ")

# Is the relative error EXACTLY identical for e and mu?
print(f"\n{'='*70}")
print("DEEPER: IS IT REALLY A PURE SCALE ERROR?")
print("=" * 70)
# If it's a pure scale error, then for ALL three:
# m_predicted = m_true * (1 + epsilon)
# So m_pred/m_exp should be the same for all three
for i in range(3):
    ratio = masses_a2[i] / exp_list[i]
    print(f"  {names[i]:>8}: m_pred/m_exp = {ratio:.12f}  (1 + {(ratio-1):.6e})")

print(f"\n  If pure scale error, all three (1+ε) values would be identical.")
print(f"  e vs μ differ by {abs(rel_e - rel_mu):.3e} (relative)")
print(f"  This is {abs(rel_e - rel_mu)/abs(rel_mu)*100:.1f}% different — NOT identical")
