"""
CKM Derivation Verification
============================
Verifies the algebraic identities underlying the CKM derivation
from the DS framework at H=3.

Key identity: (1 - alpha) / H = K*
where alpha = H/(H^2+1) is the channel discount (PMNS solar angle)
and K* = (H^2-H+1)/(H(H^2+1)) is the vacuum conflict.
"""

import numpy as np
from fractions import Fraction

H = 3

print("=" * 70)
print("CKM DERIVATION VERIFICATION")
print("=" * 70)

# ── Fundamental quantities ──
K_star = Fraction(H**2 - H + 1, H * (H**2 + 1))
alpha = Fraction(H, H**2 + 1)
sigma = -np.log(1 - float(K_star))

print(f"\nH = {H}")
print(f"K* = {K_star} = {float(K_star):.10f}")
print(f"alpha = {alpha} = {float(alpha):.10f}")
print(f"sigma = -ln(1 - K*) = {sigma:.10f}")

# ── KEY IDENTITY: (1 - alpha) / H = K* ──
complement = 1 - alpha
per_gen = complement / H

print(f"\n{'─' * 50}")
print("KEY IDENTITY: (1 - alpha) / H = K*")
print(f"  1 - alpha = {complement} = {float(complement):.10f}")
print(f"  (1 - alpha) / H = {per_gen} = {float(per_gen):.10f}")
print(f"  K* = {K_star} = {float(K_star):.10f}")
print(f"  MATCH: {per_gen == K_star}")

# ── Cyclotomic structure ──
Phi1 = H - 1  # x-1 at x=H... actually Phi_1(x)=x-1, but we use Phi_n as defined
# Let me use the paper's convention:
# Phi_4(H) = H^2 + 1 = 10
# Phi_6(H) = H^2 - H + 1 = 7
Phi4 = H**2 + 1
Phi6 = H**2 - H + 1

print(f"\n{'─' * 50}")
print("CYCLOTOMIC STRUCTURE")
print(f"  Phi_4(H) = H^2 + 1 = {Phi4}")
print(f"  Phi_6(H) = H^2 - H + 1 = {Phi6}")
print(f"  K* = Phi_6 / (H * Phi_4) = {Phi6}/{H * Phi4} = {Fraction(Phi6, H * Phi4)}")
print(f"  Cyclotomic Pythagoras: H + Phi_6 = Phi_4 => {H} + {Phi6} = {Phi4}: {H + Phi6 == Phi4}")
print(f"  cos^2(theta_12) = Phi_6/Phi_4 = {Fraction(Phi6, Phi4)} = {Phi6/Phi4:.4f}")
print(f"  sin^2(theta_12) = H/Phi_4 = {Fraction(H, Phi4)} = {H/Phi4:.4f}")

# ── Full CKM matrix ──
print(f"\n{'─' * 50}")
print("FULL CKM MATRIX")

A = Fraction(H, H + 1)
lam = K_star

# CKM elements
V_ud = 1 - lam**2 / 2
V_us = lam
V_ub = lam / (H + 1)**3
V_cd = lam
V_cs = 1 - lam**2 / 2
V_cb = A * lam**2
V_td = lam / H**3
V_ts = A * lam**2
V_tb = Fraction(1)

ckm = [
    ('V_ud', V_ud, 0.97435),
    ('V_us', V_us, 0.2253),
    ('V_ub', V_ub, 0.00356),
    ('V_cd', V_cd, 0.2210),
    ('V_cs', V_cs, 0.97345),
    ('V_cb', V_cb, 0.04082),
    ('V_td', V_td, 0.00859),
    ('V_ts', V_ts, 0.0415),
    ('V_tb', V_tb, 0.99915),
]

print(f"\n  lambda = K* = {lam} = {float(lam):.6f}")
print(f"  A = H/(H+1) = {A} = {float(A):.6f}")
print()

for name, pred, meas in ckm:
    dev = abs(float(pred) - meas) / meas * 100
    print(f"  |{name}| = {str(pred):>12s} = {float(pred):.6f}  (PDG: {meas:.5f}, dev: {dev:.2f}%)")

# Corner ratio
ratio_pred = Fraction((H+1)**3, H**3)
ratio_meas = 0.00859 / 0.00356
print(f"\n  |V_td|/|V_ub| = (H+1)^3/H^3 = {ratio_pred} = {float(ratio_pred):.4f}")
print(f"  Measured: {ratio_meas:.4f} (dev: {abs(float(ratio_pred) - ratio_meas)/ratio_meas*100:.1f}%)")

# ── CP phase and Jarlskog ──
print(f"\n{'─' * 50}")
print("CP PHASE AND JARLSKOG INVARIANT")

delta = np.pi - 2
print(f"  delta_CKM = pi - (H-1) = pi - 2 = {delta:.6f} rad")
print(f"  Measured: 1.144 +/- 0.027 rad (dev: {abs(delta - 1.144)/1.144*100:.2f}%)")

# Jarlskog
s12 = float(V_us)
s23 = float(V_cb) / float(V_us)**2 * float(V_us)**2  # A * lambda^2 / lambda ... let me just use s_ij
c12 = np.sqrt(1 - s12**2)
s13 = float(V_ub)
c13 = np.sqrt(1 - s13**2)
s23_ckm = float(V_cb)  # sin(theta_23) in CKM parametrization
c23_ckm = np.sqrt(1 - s23_ckm**2)

J = c12 * c13**2 * c23_ckm * s12 * s13 * s23_ckm * np.sin(delta)
print(f"  J = {J:.4e}")
print(f"  Measured: 3.08e-5 (dev: {abs(J - 3.08e-5)/3.08e-5*100:.1f}%)")

# ── PMNS angles ──
print(f"\n{'─' * 50}")
print("PMNS MIXING ANGLES")

sin2_12 = Fraction(H, H**2 + 1)
sin2_23 = Fraction(1, 2) + K_star / (2 * H)
sin2_13_exact = sigma / (H * (H + 1))

print(f"  sin^2(theta_12) = H/(H^2+1) = {sin2_12} = {float(sin2_12):.4f}  (PDG: 0.307, dev: {abs(float(sin2_12)-0.307)/0.307*100:.1f}%)")
print(f"  sin^2(theta_23) = 1/2 + K*/(2H) = {sin2_23} = {float(sin2_23):.4f}  (PDG: 0.546, dev: {abs(float(sin2_23)-0.546)/0.546*100:.1f}%)")
print(f"  sin^2(theta_13) = sigma/(H(H+1)) = {sin2_13_exact:.6f}  (PDG: 0.0220, dev: {abs(sin2_13_exact-0.0220)/0.0220*100:.1f}%)")

# PMNS CP phase
delta_pmns = -np.pi / (H - 1)
print(f"  delta_PMNS = -pi/(H-1) = -pi/2 = {delta_pmns:.6f} rad")

# ── Complementarity check ──
print(f"\n{'─' * 50}")
print("QUARK-LEPTON COMPLEMENTARITY")
print(f"  K* * H + sin^2(theta_12) = {float(K_star * H + sin2_12):.10f}")
print(f"  Should equal 1: {K_star * H + sin2_12 == 1}")
print(f"  cos^2(theta_12) / H = {float(Fraction(Phi6, Phi4) / H):.10f}")
print(f"  K* = {float(K_star):.10f}")
print(f"  Match: {Fraction(Phi6, Phi4) / H == K_star}")

# ── Nuclear quantities ──
print(f"\n{'─' * 50}")
print("NUCLEAR STRUCTURAL QUANTITIES")

m_p = 938.272  # MeV
m_pi = 139.570  # MeV
hbar_c = 197.327  # MeV fm
lambda_C = hbar_c / m_p
r_p_pred = (H + 1) * lambda_C
r_p_meas = 0.8414  # fm (muonic H)

print(f"  Compton wavelength = hbar/(m_p c) = {lambda_C:.5f} fm")
print(f"  r_p = (H+1) * lambda_C = {r_p_pred:.4f} fm")
print(f"  Measured: {r_p_meas} fm (dev: {abs(r_p_pred - r_p_meas)/r_p_meas*100:.2f}%)")

E_ope = m_pi**2 / (2 * m_p)
BA_pred = E_ope * float(1 - K_star)
print(f"\n  E_OPE = m_pi^2/(2m_N) = {E_ope:.2f} MeV")
print(f"  B/A = E_OPE * (1-K*) = {BA_pred:.2f} MeV")
print(f"  Measured: ~8.0 MeV (dev: {abs(BA_pred - 8.0)/8.0*100:.1f}%)")

mu_ratio = Fraction(H, H - 1)
print(f"\n  mu_p/|mu_n| = H/(H-1) = {mu_ratio} = {float(mu_ratio):.4f}")
print(f"  Measured: 1.4599 (dev: {abs(float(mu_ratio) - 1.4599)/1.4599*100:.1f}%)")

# ── GUT hierarchy ──
print(f"\n{'─' * 50}")
print("GUT HIERARCHY")

alpha_gut_inv = H * (H**2 - 1)
log_ratio = H * (H**2 + 1) + H  # h(E8) + H
M_Z = 91.1876  # GeV
M_GUT = M_Z * np.exp(log_ratio)
M_Pl = 1.221e19  # GeV

print(f"  1/alpha_GUT = H(H^2-1) = {alpha_gut_inv}")
print(f"  ln(M_GUT/M_Z) = h(E8) + H = {H*(H**2+1)} + {H} = {log_ratio}")
print(f"  M_GUT = M_Z * e^33 = {M_GUT:.3e} GeV")
print(f"  M_Pl / M_GUT = {M_Pl / M_GUT:.1f}")
print(f"  H^3 * (H^3 - 4) = {H**3 * (H**3 - 4)}")
print(f"  Match: {abs(M_Pl/M_GUT - H**3*(H**3-4)) / (H**3*(H**3-4)) * 100:.1f}%")

# Instanton action
S = Fraction(H**3, K_star)
print(f"\n  Instanton action S = H^3/K* = {S} = {float(S):.4f}")
print(f"  S + 1/K* = {S + 1/K_star} = {float(S + 1/K_star)}")
print(f"  5! = {120}")
print(f"  Match: {S + 1/K_star == 120}")

print(f"\n{'=' * 70}")
print("ALL VERIFICATIONS COMPLETE")
print(f"{'=' * 70}")
