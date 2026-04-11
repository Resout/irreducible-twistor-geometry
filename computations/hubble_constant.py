"""
Hubble constant H_0 from H=3 informational geometry.

Framework inputs:
  - H = 3 (irreducibility dimension)
  - Omega_Lambda = (H-1)/H = 2/3
  - rho_Lambda = 2.59e-47 GeV^4  (derived from instanton hierarchy)
  - w = -1 exactly (cosmological constant)

Friedmann equation:  H_0^2 = 8*pi*G * rho_c / 3
  where rho_c = rho_Lambda / Omega_Lambda
"""

import numpy as np

# ============================================================
# Physical constants (CODATA 2018 / PDG)
# ============================================================
G_N       = 6.67430e-11      # Newton's constant [m^3 kg^-1 s^-2]
c_light   = 2.99792458e8     # speed of light [m/s]
hbar      = 1.054571817e-34  # reduced Planck constant [J s]
hbar_c_GeV_fm = 0.1973269804 # hbar*c in [GeV fm]
hbar_c    = hbar_c_GeV_fm * 1e-15  # hbar*c in [GeV m]
GeV_per_kg = 1.0 / 1.78266192e-27  # 1 kg in GeV/c^2  =>  1 GeV/c^2 = 1.78266192e-27 kg
kg_per_GeV = 1.78266192e-27  # 1 GeV/c^2 in kg
Mpc_in_m  = 3.0856775814913673e22  # 1 Mpc in meters
km_per_m  = 1e-3

# ============================================================
# Framework parameters
# ============================================================
H = 3
Omega_Lambda = (H - 1) / H  # = 2/3

rho_Lambda_nat_framework = 2.59e-47   # GeV^4 (framework prediction)
rho_Lambda_nat_planck    = 2.52e-47   # GeV^4 (Planck 2018 measured)

# ============================================================
# Unit conversion:  rho [GeV^4] (natural units) -> rho [kg/m^3] (SI)
#
# In natural units (hbar=c=1), energy density has dimension [GeV^4].
# To convert:
#   rho_SI = rho_nat * (1/hbar*c)^3 * (1 GeV/c^2)
#          = rho_nat / (hbar*c)^3 * kg_per_GeV
#
# where (hbar*c)^3 has units [GeV^3 m^3] and kg_per_GeV converts
# the remaining GeV to kg.
# ============================================================
def rho_nat_to_SI(rho_nat):
    """Convert energy density from GeV^4 (natural) to kg/m^3 (SI)."""
    return rho_nat / (hbar_c**3) * kg_per_GeV

# ============================================================
# Compute H_0
# ============================================================
def compute_H0(rho_Lambda_nat, label=""):
    rho_Lambda_SI = rho_nat_to_SI(rho_Lambda_nat)
    rho_c_SI = rho_Lambda_SI / Omega_Lambda   # rho_c = rho_Lambda / Omega_Lambda

    # Friedmann: H_0^2 = 8*pi*G*rho_c / 3
    H0_sq = 8 * np.pi * G_N * rho_c_SI / 3
    H0_per_s = np.sqrt(H0_sq)   # [s^-1]

    # Convert to km/s/Mpc
    H0_km_s_Mpc = H0_per_s * Mpc_in_m * km_per_m

    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    print(f"  rho_Lambda  = {rho_Lambda_nat:.2e} GeV^4  (natural units)")
    print(f"  rho_Lambda  = {rho_Lambda_SI:.4e} kg/m^3  (SI)")
    print(f"  Omega_Lambda = {Omega_Lambda:.6f}")
    print(f"  rho_c       = {rho_c_SI:.4e} kg/m^3  (SI)")
    print(f"  H_0         = {H0_per_s:.4e} s^-1")
    print(f"  H_0         = {H0_km_s_Mpc:.2f} km/s/Mpc")

    return H0_km_s_Mpc, rho_c_SI

print("╔══════════════════════════════════════════════════════════╗")
print("║   Hubble Constant from H=3 Informational Geometry      ║")
print("╚══════════════════════════════════════════════════════════╝")

print(f"\n  Framework dimension H = {H}")
print(f"  Dark energy fraction Omega_Lambda = (H-1)/H = {Omega_Lambda:.10f}")
print(f"  Equation of state w = -1 (exact, cosmological constant)")

# --- Framework prediction ---
H0_fw, rho_c_fw = compute_H0(rho_Lambda_nat_framework,
                               "FRAMEWORK: rho_Lambda = 2.59e-47 GeV^4")

# --- Planck measured ---
H0_pl, rho_c_pl = compute_H0(rho_Lambda_nat_planck,
                               "PLANCK MEASURED: rho_Lambda = 2.52e-47 GeV^4")

# ============================================================
# Sanity check: known critical density
# ============================================================
print(f"\n{'='*60}")
print(f"  SANITY CHECK: critical density")
print(f"{'='*60}")
rho_c_known = 9.47e-27  # kg/m^3  (for H_0 ~ 67-73)
print(f"  Known rho_c ~ 9.47e-27 kg/m^3  (for H_0 ~ 70 km/s/Mpc)")
print(f"  Framework rho_c = {rho_c_fw:.4e} kg/m^3")
print(f"  Planck    rho_c = {rho_c_pl:.4e} kg/m^3")
print(f"  Ratio framework/known = {rho_c_fw / rho_c_known:.4f}")

# ============================================================
# Comparison with measurements
# ============================================================
H0_CMB = 67.4;   sigma_CMB = 0.5    # Planck 2018 CMB
H0_SH0ES = 73.04; sigma_SH0ES = 1.04  # SH0ES 2022

print(f"\n{'='*60}")
print(f"  COMPARISON WITH OBSERVATIONS")
print(f"{'='*60}")

for label, H0_val in [("Framework (rho=2.59e-47)", H0_fw),
                       ("Planck rho (rho=2.52e-47)", H0_pl)]:
    dev_CMB = (H0_val - H0_CMB) / sigma_CMB
    dev_SH0ES = (H0_val - H0_SH0ES) / sigma_SH0ES

    print(f"\n  {label}:  H_0 = {H0_val:.2f} km/s/Mpc")
    print(f"    vs Planck CMB  (67.4 +/- 0.5):   deviation = {dev_CMB:+.2f} sigma")
    print(f"    vs SH0ES local (73.04 +/- 1.04):  deviation = {dev_SH0ES:+.2f} sigma")

    midpoint = (H0_CMB + H0_SH0ES) / 2
    print(f"    Midpoint of tension = {midpoint:.2f} km/s/Mpc")
    print(f"    Distance from midpoint = {abs(H0_val - midpoint):.2f} km/s/Mpc")

# ============================================================
# The Hubble tension
# ============================================================
print(f"\n{'='*60}")
print(f"  THE HUBBLE TENSION")
print(f"{'='*60}")
tension_sigma = (H0_SH0ES - H0_CMB) / np.sqrt(sigma_CMB**2 + sigma_SH0ES**2)
print(f"  CMB vs local tension: {tension_sigma:.1f} sigma")
print(f"  CMB:   {H0_CMB} +/- {sigma_CMB}")
print(f"  SH0ES: {H0_SH0ES} +/- {sigma_SH0ES}")
print(f"  Framework H_0 = {H0_fw:.2f} km/s/Mpc")
print(f"  Sits between the two? {H0_CMB < H0_fw < H0_SH0ES}")

# ============================================================
# Framework derivation chain summary
# ============================================================
print(f"\n{'='*60}")
print(f"  DERIVATION CHAIN")
print(f"{'='*60}")
print(f"  1. Irreducibility dimension H = 3")
print(f"  2. Dark energy fraction: Omega_Lambda = (H-1)/H = {Omega_Lambda}")
print(f"     (Observed: 0.6847 +/- 0.0073 [Planck 2018])")
obs_OL = 0.6847; sig_OL = 0.0073
print(f"     Deviation: {(Omega_Lambda - obs_OL)/sig_OL:+.1f} sigma")
print(f"  3. Vacuum energy rho_Lambda from instanton hierarchy")
print(f"  4. rho_c = rho_Lambda / Omega_Lambda = (3/2) rho_Lambda")
print(f"  5. Friedmann: H_0 = sqrt(8*pi*G*rho_c/3)")
print(f"  6. Result: H_0 = {H0_fw:.2f} km/s/Mpc")
print(f"\n  Free parameters used: 0 (H=3 is derived, rho_Lambda is derived)")
