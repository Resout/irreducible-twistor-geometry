import numpy as np
from fractions import Fraction

H = 3
Kstar = 7/30
hE8 = 30  # h(E8) = H(H^2+1) = 30, derived from McKay correspondence

# ============================================================
# KOIDE FORMULA FROM THE FRAMEWORK
# ============================================================

print("="*70)
print("KOIDE FORMULA: LEPTON MASSES FROM H = 3")
print("="*70)

# Q = 2/3 (from S3 generation symmetry)
Q = 2/3

# theta_phys = 2/H^2 = 2/9
# derived as crystal angle minus E8 Coxeter correction
theta_crystal = (1 - Kstar) / H        # = 23/90
theta_E8_correction = 1 / hE8          # = 1/30 = 3/90
theta_phys = theta_crystal - theta_E8_correction  # = 20/90 = 2/9

print(f"Q = 2/3 = {Q:.10f}")
print(f"theta_crystal = (1-K*)/H = 23/90 = {theta_crystal:.10f}")
print(f"E8 correction = 1/h(E8) = 1/30 = {theta_E8_correction:.10f}")
print(f"theta_phys = 2/H^2 = 2/9 = {theta_phys:.10f}")

# Verify as exact fractions
Kstar_frac = Fraction(7, 30)
theta_frac = Fraction(1 - Kstar_frac.numerator/Kstar_frac.denominator) - Fraction(1, 30)
# Direct: (1-7/30)/3 - 1/30 = (23/30)/3 - 1/30 = 23/90 - 3/90 = 20/90 = 2/9
theta_exact = Fraction(23, 90) - Fraction(1, 30)
print(f"\nExact: theta = {theta_exact} = {float(theta_exact):.10f}")
print(f"2/H^2 = {Fraction(2, H**2)} = {float(Fraction(2, H**2)):.10f}")
print(f"Match (exact rational): {theta_exact == Fraction(2, H**2)}")

print(f"\nKey identity: 1-K* = 2/H + H/h(E8)")
print(f"  23/30 = 2/3 + 1/10 = {2/3 + 1/10:.10f}")
print(f"  Match: {Fraction(23,30) == Fraction(2,3) + Fraction(1,10)}")

# ============================================================
# MASS RATIOS FROM KOIDE PARAMETERIZATION
# ============================================================

print(f"\n" + "="*70)
print("MASS RATIOS")
print("="*70)

theta = theta_phys  # 2/9

x = []
for k in range(3):
    val = (1 + np.sqrt(2) * np.cos(theta + 2*np.pi*k/3))**2
    x.append(val)

x.sort()  # lightest first: e, mu, tau

print(f"\nKoide parameters at theta = 2/9:")
print(f"  x_e   = {x[0]:.10f}")
print(f"  x_mu  = {x[1]:.10f}")
print(f"  x_tau = {x[2]:.10f}")

# Measured masses (CODATA 2018)
m_e_meas = 0.51099895000  # MeV
m_mu_meas = 105.6583755   # MeV
m_tau_meas = 1776.86      # MeV

print(f"\nPredicted mass ratios vs measured:")
print(f"  m_tau/m_e = {x[2]/x[0]:.4f}  (measured: {m_tau_meas/m_e_meas:.4f}, dev: {abs(x[2]/x[0]-m_tau_meas/m_e_meas)/(m_tau_meas/m_e_meas)*100:.4f}%)")
print(f"  m_mu/m_e  = {x[1]/x[0]:.4f}  (measured: {m_mu_meas/m_e_meas:.4f}, dev: {abs(x[1]/x[0]-m_mu_meas/m_e_meas)/(m_mu_meas/m_e_meas)*100:.4f}%)")
print(f"  m_tau/m_mu = {x[2]/x[1]:.6f}  (measured: {m_tau_meas/m_mu_meas:.6f}, dev: {abs(x[2]/x[1]-m_tau_meas/m_mu_meas)/(m_tau_meas/m_mu_meas)*100:.4f}%)")

# Scale from tau mass
M_sq = m_tau_meas / x[2]
m_e_pred = M_sq * x[0]
m_mu_pred = M_sq * x[1]

print(f"\nAbsolute masses (m_tau = {m_tau_meas} MeV as input):")
print(f"  Scale M^2 = {M_sq:.6f} MeV")
print(f"  m_e  = {m_e_pred:.6f} MeV  (measured: {m_e_meas:.8f}, error: {abs(m_e_pred-m_e_meas)/m_e_meas*100:.4f}%)")
print(f"  m_mu = {m_mu_pred:.4f} MeV  (measured: {m_mu_meas:.7f}, error: {abs(m_mu_pred-m_mu_meas)/m_mu_meas*100:.4f}%)")

# Verify Q
Q_check = (m_e_pred + m_mu_pred + m_tau_meas) / (np.sqrt(m_e_pred) + np.sqrt(m_mu_pred) + np.sqrt(m_tau_meas))**2
print(f"\nKoide Q: {Q_check:.12f} (should be exactly 2/3 = {2/3:.12f})")
print(f"Error: {abs(Q_check - 2/3):.2e}")

# ============================================================
# QUARK TRIPLETS
# ============================================================

print(f"\n" + "="*70)
print("QUARK TRIPLET Q VALUES")
print("="*70)

triplets = {
    "Up-type  (u, c, t)": [2.16, 1270, 173000],
    "Down-type (d, s, b)": [4.70, 93.5, 4180],
    "Leptons (e, mu, tau)": [m_e_meas, m_mu_meas, m_tau_meas],
}
for name, masses in triplets.items():
    Q_t = sum(masses) / (sum(np.sqrt(m) for m in masses))**2
    print(f"  {name}: Q = {Q_t:.6f}  (2/3 = {2/3:.6f}, dev = {abs(Q_t-2/3)/(2/3)*100:.1f}%)")

print(f"\nOnly the lepton triplet satisfies Q = 2/3.")
