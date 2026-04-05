#!/usr/bin/env python3
"""
Koide θ feedback loop: does self-consistent M(θ) shift the equilibrium θ?

The system:
  1. θ → x_μ = (1 + √2 cos(θ + 4π/3))²
  2. x_μ → M = m_proton / (H × (1 - (H-1)×x_μ/195))
  3. M, θ → m_k = M × (1 + √2 cos(θ + 2πk/3))²

Question: when we solve for θ self-consistently (minimizing distance
to experimental masses WITH M = M(θ)), does the equilibrium shift
from 2/9 by anything near 1.75×10⁻⁷?
"""
import numpy as np
from scipy.optimize import minimize_scalar, brentq

# Constants
H = 3
m_proton = 938.272088  # MeV (CODATA 2018)

# Experimental lepton masses (MeV, CODATA/PDG)
m_e_exp = 0.51099895000   # ± 0.00000000015
m_mu_exp = 105.6583755     # ± 0.0000023
m_tau_exp = 1776.86         # ± 0.12

masses_exp = np.array([m_tau_exp, m_e_exp, m_mu_exp])  # k=0,1,2

# Phase offsets for k=0,1,2
phases = np.array([0, 2*np.pi/3, 4*np.pi/3])


def koide_multipliers(theta):
    """Return (1 + √2 cos(θ + 2πk/3))² for k=0,1,2."""
    return (1 + np.sqrt(2) * np.cos(theta + phases))**2


def M_from_theta(theta):
    """Self-consistent M via 195 formula."""
    x_mu = koide_multipliers(theta)[2]  # k=2 = muon
    return m_proton / (H * (1 - (H - 1) * x_mu / 195))


def masses_from_theta(theta):
    """Full mass prediction: M(θ) × multipliers(θ)."""
    M = M_from_theta(theta)
    return M * koide_multipliers(theta)


def masses_fixed_M(theta, M):
    """Mass prediction with fixed M (no feedback)."""
    return M * koide_multipliers(theta)


# ============================================================
# Method 1: Find θ that minimizes χ² to experimental masses
#            WITH self-consistent M(θ)
# ============================================================
def chi2_feedback(theta):
    """χ² with M = M(θ) feedback."""
    pred = masses_from_theta(theta)
    # Weight by inverse experimental uncertainty squared
    sigmas = np.array([0.12, 0.00000000015, 0.0000023])
    return np.sum(((pred - masses_exp) / sigmas)**2)


def chi2_fixed(theta, M):
    """χ² with fixed M."""
    pred = masses_fixed_M(theta, M)
    sigmas = np.array([0.12, 0.00000000015, 0.0000023])
    return np.sum(((pred - masses_exp) / sigmas)**2)


print("=" * 70)
print("KOIDE θ FEEDBACK ANALYSIS")
print("=" * 70)

theta_theory = 2.0 / 9.0
print(f"\nTheoretical θ = 2/9 = {theta_theory:.15f}")

# ============================================================
# Baseline: what M and masses does θ = 2/9 give?
# ============================================================
M_at_theory = M_from_theta(theta_theory)
masses_at_theory = masses_from_theta(theta_theory)
print(f"\nAt θ = 2/9:")
print(f"  M(θ) = {M_at_theory:.6f} MeV")
print(f"  m_tau = {masses_at_theory[0]:.4f} MeV  (exp: {m_tau_exp:.4f}, Δ = {masses_at_theory[0]-m_tau_exp:+.4f})")
print(f"  m_e   = {masses_at_theory[1]:.8f} MeV  (exp: {m_e_exp:.8f}, Δ = {masses_at_theory[1]-m_e_exp:+.8f})")
print(f"  m_mu  = {masses_at_theory[2]:.6f} MeV  (exp: {m_mu_exp:.6f}, Δ = {masses_at_theory[2]-m_mu_exp:+.6f})")

# ============================================================
# Fit 1: Optimize θ with M(θ) feedback
# ============================================================
print(f"\n{'='*70}")
print("FIT 1: θ optimized with M = M(θ) feedback")
print("=" * 70)

result_fb = minimize_scalar(chi2_feedback, bounds=(0.220, 0.225), method='bounded')
theta_fb = result_fb.x
M_fb = M_from_theta(theta_fb)
masses_fb = masses_from_theta(theta_fb)

print(f"  θ_fitted (feedback) = {theta_fb:.15f}")
print(f"  Δθ from 2/9         = {theta_fb - theta_theory:+.3e}")
print(f"  M(θ_fitted)         = {M_fb:.6f} MeV")
print(f"  m_tau = {masses_fb[0]:.4f}  (Δ = {masses_fb[0]-m_tau_exp:+.4f})")
print(f"  m_e   = {masses_fb[1]:.11f}  (Δ = {masses_fb[1]-m_e_exp:+.3e})")
print(f"  m_mu  = {masses_fb[2]:.7f}  (Δ = {masses_fb[2]-m_mu_exp:+.3e})")

# ============================================================
# Fit 2: Optimize θ with FIXED M (no feedback)
#         Use M derived at θ = 2/9 as fixed scale
# ============================================================
print(f"\n{'='*70}")
print("FIT 2: θ optimized with FIXED M (no feedback)")
print("=" * 70)

M_fixed = M_at_theory
result_nf = minimize_scalar(lambda t: chi2_fixed(t, M_fixed), bounds=(0.220, 0.225), method='bounded')
theta_nf = result_nf.x
masses_nf = masses_fixed_M(theta_nf, M_fixed)

print(f"  θ_fitted (no feedback) = {theta_nf:.15f}")
print(f"  Δθ from 2/9            = {theta_nf - theta_theory:+.3e}")
print(f"  M (fixed)              = {M_fixed:.6f} MeV")
print(f"  m_tau = {masses_nf[0]:.4f}  (Δ = {masses_nf[0]-m_tau_exp:+.4f})")
print(f"  m_e   = {masses_nf[1]:.11f}  (Δ = {masses_nf[1]-m_e_exp:+.3e})")
print(f"  m_mu  = {masses_nf[2]:.7f}  (Δ = {masses_nf[2]-m_mu_exp:+.3e})")

# ============================================================
# Fit 3: BOTH M and θ free (what M does experiment actually want?)
# ============================================================
print(f"\n{'='*70}")
print("FIT 3: θ optimized with FREE M (best possible)")
print("=" * 70)

from scipy.optimize import minimize

def chi2_both(params):
    theta, M = params
    pred = M * koide_multipliers(theta)
    sigmas = np.array([0.12, 0.00000000015, 0.0000023])
    return np.sum(((pred - masses_exp) / sigmas)**2)

res_both = minimize(chi2_both, [theta_theory, M_at_theory], method='Nelder-Mead',
                    options={'xatol': 1e-16, 'fatol': 1e-10, 'maxiter': 100000})
theta_free, M_free = res_both.x
masses_free = M_free * koide_multipliers(theta_free)

print(f"  θ_fitted (free)  = {theta_free:.15f}")
print(f"  M_fitted (free)  = {M_free:.6f} MeV")
print(f"  Δθ from 2/9      = {theta_free - theta_theory:+.3e}")
print(f"  ΔM from M(2/9)   = {M_free - M_at_theory:+.6f} MeV  ({(M_free-M_at_theory)/M_at_theory*1e6:+.2f} ppm)")
print(f"  m_tau = {masses_free[0]:.4f}  (Δ = {masses_free[0]-m_tau_exp:+.4f})")
print(f"  m_e   = {masses_free[1]:.11f}  (Δ = {masses_free[1]-m_e_exp:+.3e})")
print(f"  m_mu  = {masses_free[2]:.7f}  (Δ = {masses_free[2]-m_mu_exp:+.3e})")

# ============================================================
# Key comparison: does feedback shift θ?
# ============================================================
print(f"\n{'='*70}")
print("FEEDBACK EFFECT")
print("=" * 70)
print(f"  θ_fitted (with feedback):    {theta_fb:.15f}")
print(f"  θ_fitted (without feedback): {theta_nf:.15f}")
print(f"  θ_fitted (both free):        {theta_free:.15f}")
print(f"  θ = 2/9:                     {theta_theory:.15f}")
print()
shift = theta_fb - theta_nf
print(f"  Feedback-induced shift: {shift:+.3e} radians")
print(f"  Target anomaly:         {-1.75e-7:+.3e} radians")
if abs(shift) > 1e-12:
    print(f"  Ratio (shift/anomaly):  {shift / (-1.75e-7):.4f}")
    print(f"  Feedback explains {abs(shift/(-1.75e-7))*100:.2f}% of the anomaly")
else:
    print(f"  Feedback shift is negligible (< 1e-12)")

# ============================================================
# Jacobian: sensitivity of θ_fit to M
# ============================================================
print(f"\n{'='*70}")
print("SENSITIVITY ANALYSIS")
print("=" * 70)

# dθ/dM at the fitted point
dM = 0.001  # MeV
res_plus = minimize_scalar(lambda t: chi2_fixed(t, M_at_theory + dM), bounds=(0.220, 0.225), method='bounded')
res_minus = minimize_scalar(lambda t: chi2_fixed(t, M_at_theory - dM), bounds=(0.220, 0.225), method='bounded')
dtheta_dM = (res_plus.x - res_minus.x) / (2 * dM)

print(f"  dθ/dM = {dtheta_dM:.3e} rad/MeV")
print(f"  To get Δθ = -1.75e-7, need ΔM = {-1.75e-7 / dtheta_dM:.6f} MeV")
print(f"  That's {-1.75e-7 / dtheta_dM / M_at_theory * 1e6:.2f} ppm of M")

# x_mu sensitivity
print(f"\n  x_μ at θ=2/9: {koide_multipliers(theta_theory)[2]:.10f}")
print(f"  x_μ at θ_fit: {koide_multipliers(theta_free)[2]:.10f}")
print(f"  Δx_μ:         {koide_multipliers(theta_free)[2] - koide_multipliers(theta_theory)[2]:+.3e}")

# ============================================================
# The Jacobian of the full coupled system
# ============================================================
print(f"\n{'='*70}")
print("COUPLED SYSTEM JACOBIAN")
print("=" * 70)

# The coupled system: θ → x_μ(θ) → M(θ) → masses(θ)
# Compute d(masses)/dθ including the M(θ) chain

eps = 1e-10
masses_plus = masses_from_theta(theta_theory + eps)
masses_minus = masses_from_theta(theta_theory - eps)
dm_dtheta_coupled = (masses_plus - masses_minus) / (2 * eps)

masses_plus_nc = masses_fixed_M(theta_theory + eps, M_at_theory)
masses_minus_nc = masses_fixed_M(theta_theory - eps, M_at_theory)
dm_dtheta_uncoupled = (masses_plus_nc - masses_minus_nc) / (2 * eps)

print(f"  dm_k/dθ (with M feedback):")
labels = ['tau', 'e', 'mu']
for i in range(3):
    ratio = dm_dtheta_coupled[i] / dm_dtheta_uncoupled[i] if abs(dm_dtheta_uncoupled[i]) > 0 else 0
    print(f"    {labels[i]:>4}: coupled = {dm_dtheta_coupled[i]:+.4f}  uncoupled = {dm_dtheta_uncoupled[i]:+.4f}  ratio = {ratio:.6f}")

dM_dtheta = (M_from_theta(theta_theory + eps) - M_from_theta(theta_theory - eps)) / (2 * eps)
print(f"\n  dM/dθ = {dM_dtheta:.6f} MeV/rad")
print(f"  At θ=2/9, a shift of -1.75e-7 changes M by {dM_dtheta * (-1.75e-7):.3e} MeV")


if __name__ == '__main__':
    pass
