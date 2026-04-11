#!/usr/bin/env python3
"""
Verify the two remaining conjectures:
1. conj:neutrinos — M_R = M_GUT/33, m_ν₃ = v²/(2M_R)
2. conj:ckm_full — all 9 CKM elements from K* and H
"""

from mpmath import mp, mpf, pi, sin, cos, sqrt, nstr, fabs, ln, exp

mp.dps = 30
H = 3
K = mpf(7)/30

print("=" * 70)
print("CONJECTURE 1: NEUTRINO MASS SCALE")
print("=" * 70)

# Framework parameters
v_multiplier = mpf((H+1)**2 * H**2)  # = 144
m0 = mpf("1.710")  # GeV
v = v_multiplier * m0 / 1000  # Convert: 144 × 1.710 MeV... no
# v is in units of m(0++): v/m(0++) = 144
# Physical v = 144 × m(0++) = 144 × 1710 MeV = 246240 MeV = 246.24 GeV
v_GeV = mpf(144) * mpf("1.710")  # GeV
print(f"  v = 144 × m(0++) = {nstr(v_GeV, 6)} GeV (actual: 246.22 GeV)")

# GUT scale
M_GUT = mpf("1.96e16")  # GeV (standard estimate)

# M_R = M_GUT / (h(E₈) + H) = M_GUT / 33
h_E8 = 30  # Coxeter number of E₈ = H(H²+1)
running_exp = h_E8 + H  # = 33
M_R = M_GUT / running_exp
print(f"\n  h(E₈) = {h_E8}, h(E₈)+H = {running_exp}")
print(f"  M_R = M_GUT/{running_exp} = {nstr(M_R, 4)} GeV")

# Seesaw: m_ν = v²/(2M_R)
m_nu = v_GeV**2 / (2 * M_R)
m_nu_eV = m_nu * mpf("1e9")  # GeV to eV
print(f"\n  m_ν₃ = v²/(2M_R) = {nstr(m_nu, 4)} GeV = {nstr(m_nu_eV, 4)} eV")
print(f"  Observed: √(Δm²_atm) ≈ 0.050 eV")
print(f"  Deviation: {nstr(fabs(m_nu_eV - mpf('0.050'))/mpf('0.050')*100, 3)}%")

# Can we express this entirely in framework terms?
# m_ν = v² / (2M_R) = (144 × m₀)² / (2 × M_GUT/33)
# = 144² × 33 × m₀² / (2 × M_GUT)
# The ratio m₀/M_GUT is the only non-framework quantity
print(f"\n  Framework form: m_ν = 144² × 33 / (2 × M_GUT) × m₀²")
print(f"  = {144**2 * 33 // 2} × m₀²/M_GUT")
print(f"  = {144**2 * 33 // 2} × ({nstr(m0, 4)} GeV)²/({nstr(M_GUT, 3)} GeV)")

# What determines M_GUT in the framework?
# The hierarchy tower: M_GUT/M_Z = exp(33 × something)
# Or: M_GUT = m(0++) × (some large multiplier)
# M_GUT/m(0++) = 1.96e16/1.710 = 1.146e16
print(f"\n  M_GUT/m(0++) = {nstr(M_GUT/m0, 4)}")
print(f"  = {nstr(M_GUT/m0, 4)} ← this multiplier is NOT derived from H and K*")

print("\n" + "=" * 70)
print("CONJECTURE 2: COMPLETE CKM MATRIX")
print("=" * 70)

# Framework CKM parameters
lam_W = K  # Wolfenstein λ = K*
A_W = mpf(H)/(H+1)  # Wolfenstein A = H/(H+1)
delta = pi - (H-1)  # CP phase = π - 2

# Mixing angles
s12 = lam_W  # sin(θ₁₂) = λ
s23 = A_W * lam_W**2  # sin(θ₂₃) = Aλ²
s13 = lam_W / (H+1)**3  # sin(θ₁₃) = λ/(H+1)³ = K*/(H+1)³

c12 = sqrt(1 - s12**2)
c23 = sqrt(1 - s23**2)
c13 = sqrt(1 - s13**2)

print(f"\n  Wolfenstein parameters:")
print(f"    λ = K* = {nstr(lam_W, 10)}")
print(f"    A = H/(H+1) = {nstr(A_W, 10)}")
print(f"    δ = π-2 = {nstr(delta, 10)} rad")

print(f"\n  Mixing angles:")
print(f"    sin(θ₁₂) = {nstr(s12, 10)}")
print(f"    sin(θ₂₃) = {nstr(s23, 10)}")
print(f"    sin(θ₁₃) = {nstr(s13, 10)}")

# Full CKM matrix (standard parameterization)
# V_ud = c12 c13
# V_us = s12 c13
# V_ub = s13 e^{-iδ}
# V_cd = -s12 c23 - c12 s23 s13 e^{iδ}
# V_cs = c12 c23 - s12 s23 s13 e^{iδ}
# V_cb = s23 c13
# V_td = s12 s23 - c12 c23 s13 e^{iδ}
# V_ts = -c12 s23 - s12 c23 s13 e^{iδ}
# V_tb = c23 c13

V_ud = c12 * c13
V_us = s12 * c13
V_ub = s13  # magnitude
V_cb = s23 * c13
V_tb = c23 * c13

# For V_cd, V_cs, V_td, V_ts need the phase
e_idelta_re = cos(delta)
e_idelta_im = sin(delta)

# V_cd = -s12*c23 - c12*s23*s13*(cos δ + i sin δ)
V_cd_re = -s12*c23 - c12*s23*s13*e_idelta_re
V_cd_im = -c12*s23*s13*e_idelta_im
V_cd = sqrt(V_cd_re**2 + V_cd_im**2)

V_cs_re = c12*c23 - s12*s23*s13*e_idelta_re
V_cs_im = -s12*s23*s13*e_idelta_im
V_cs = sqrt(V_cs_re**2 + V_cs_im**2)

V_td_re = s12*s23 - c12*c23*s13*e_idelta_re
V_td_im = -c12*c23*s13*e_idelta_im
V_td = sqrt(V_td_re**2 + V_td_im**2)

V_ts_re = -c12*s23 - s12*c23*s13*e_idelta_re
V_ts_im = -s12*c23*s13*e_idelta_im
V_ts = sqrt(V_ts_re**2 + V_ts_im**2)

# PDG values
pdg = {
    'V_ud': (V_ud, mpf('0.97435'), '0.97435'),
    'V_us': (V_us, mpf('0.2253'), '0.2253'),
    'V_ub': (V_ub, mpf('0.00356'), '0.00356'),
    'V_cd': (V_cd, mpf('0.2210'), '0.2210'),
    'V_cs': (V_cs, mpf('0.97345'), '0.97345'),
    'V_cb': (V_cb, mpf('0.04082'), '0.04082'),
    'V_td': (V_td, mpf('0.00859'), '0.00859'),
    'V_ts': (V_ts, mpf('0.0415'), '0.0415'),
    'V_tb': (V_tb, mpf('0.99915'), '0.99915'),
}

print(f"\n  {'Element':8s} {'Predicted':>12s} {'PDG':>12s} {'Error':>8s}")
print(f"  {'-'*8} {'-'*12} {'-'*12} {'-'*8}")
for name, (pred, obs, obs_str) in pdg.items():
    err = float(fabs(pred - obs)/obs * 100)
    print(f"  {name:8s} {float(pred):12.6f} {obs_str:>12s} {err:7.2f}%")

# Jarlskog invariant
J = c12 * c23 * c13**2 * s12 * s23 * s13 * sin(delta)
print(f"\n  Jarlskog invariant:")
print(f"    J = {nstr(J, 6)}")
print(f"    Observed: (3.08 ± 0.15) × 10⁻⁵")
print(f"    Ratio: {nstr(J / mpf('3.08e-5'), 4)}")

# Unitarity check: sum of squared magnitudes in each row
row1 = V_ud**2 + V_us**2 + V_ub**2
row2 = V_cd**2 + V_cs**2 + V_cb**2
row3 = V_td**2 + V_ts**2 + V_tb**2
print(f"\n  Unitarity check (rows should = 1):")
print(f"    Row 1: {nstr(row1, 10)}")
print(f"    Row 2: {nstr(row2, 10)}")
print(f"    Row 3: {nstr(row3, 10)}")

# What would it take to promote these?
print("\n" + "=" * 70)
print("WHAT'S NEEDED FOR PROMOTION TO CLASS A")
print("=" * 70)

print(f"""
  NEUTRINOS:
  - The Koide angle θ_ν = 7/60 is PROVED
  - The seesaw mechanism is IMPORTED (not derived)
  - M_R = M_GUT/33 needs: (a) M_GUT from framework, (b) why /33
  - The number 33 = h(E₈)+H IS structural
  - M_GUT = {nstr(M_GUT, 3)} GeV requires the running coupling equations
  - STATUS: the 33 is framework; the seesaw mechanism is external;
    M_GUT requires coupling unification (partially framework, partially SM)
  - VERDICT: cannot promote to Class A without deriving M_GUT

  CKM:
  - λ = K* = 7/30 is PROVED (equilibrium conflict)
  - A = H/(H+1) = 3/4 is PROVED (Schur uniformity)
  - δ = π-2 is VERIFIED (0.09σ), geometric interpretation clear
  - sin(θ₁₃) = K*/(H+1)³ is CONJECTURED (matches at 2.4%)
  - The formula says: 1→3 transition suppressed by mass-space volume
  - STATUS: 3 of 4 Wolfenstein parameters derived; θ₁₃ needs derivation
  - The Jarlskog J = {nstr(J, 4)} vs observed 3.08e-5 ({nstr(J/mpf('3.08e-5'), 3)})
""")

# Check: does the paper's formula |V_ub| = K*/(H+1)³ match our s13?
Vub_paper = K / (H+1)**3
print(f"  Paper's |V_ub| = K*/(H+1)³ = {nstr(Vub_paper, 8)}")
print(f"  Our s13 = K*/(H+1)³ = {nstr(s13, 8)}")
print(f"  These are identical: s13 IS the paper's V_ub formula")

# The key question: why (H+1)³ for V_ub and H³ for V_td?
Vtd_paper = K / H**3
ratio = Vtd_paper / Vub_paper
print(f"\n  Paper's |V_td| = K*/H³ = {nstr(Vtd_paper, 8)}")
print(f"  |V_td|/|V_ub| = (H+1)³/H³ = {nstr(ratio, 8)} = {float(ratio):.4f}")
print(f"  = {(H+1)**3}/{H**3} = 64/27")
print(f"  PDG: 0.00859/0.00356 = {0.00859/0.00356:.4f}")
print(f"  Framework: {float(ratio):.4f}")
print(f"  Error: {abs(float(ratio) - 0.00859/0.00356)/(0.00859/0.00356)*100:.1f}%")
