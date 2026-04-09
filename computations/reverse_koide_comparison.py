#!/usr/bin/env python3
"""
Reverse Koide: compare mass predictions from
  (A) algebraic path:   Q = 2/3,      θ = 2/9       (no Born floor)
  (B) eigenvalue path:  Q = λ₃/λ₀,   θ = K*·Δ₂/Δ₃  (through Born floor)

The Koide formula in Brannen parametrization:
  √m_k = M · (1 + r·cos(θ + 2πk/3))

where r = √(6Q - 2) and Q = (Σm)/(Σ√m)².
  - Q = 2/3 gives r = √2 (standard Koide)
  - Q ≠ 2/3 gives r ≠ √2 (floor-corrected)
"""

from mpmath import mp, mpf, sqrt, cos, pi, log, nstr, fabs

mp.dps = 50

# ============================================================
# 50-DIGIT EIGENVALUES (from a2_eigenvalue_50digit.py)
# ============================================================

lam0 = mpf('0.50216878356098934297723608630300571729693934175775')
lam1 = mpf('0.47447985206799007499882597657430324501036398090452')
lam2 = mpf('0.35265322403045475885006958862860908616028704218894')
lam3 = mpf('0.33442245925809154485648297121245740279265404305382')

K_star = mpf(7) / mpf(30)

# Mass gaps
Delta = [-log(l) for l in [lam0, lam1, lam2, lam3]]

# ============================================================
# PATH A: ALGEBRAIC (no Born floor)
# ============================================================

Q_alg = mpf(2) / mpf(3)
theta_alg = mpf(2) / mpf(9)
r_alg = sqrt(mpf(2))  # = sqrt(6·Q - 2) when Q = 2/3

# ============================================================
# PATH B: EIGENVALUE (through Born floor)
# ============================================================

Q_ev = lam3 / lam0
theta_ev = K_star * Delta[2] / Delta[3]
r_ev = sqrt(6*Q_ev - 2)

print("=" * 70)
print("REVERSE KOIDE: ALGEBRAIC vs EIGENVALUE PATH")
print("=" * 70)
print()
print("PATH A (algebraic, no floor):")
print(f"  Q   = 2/3          = {nstr(Q_alg, 20)}")
print(f"  θ   = 2/9          = {nstr(theta_alg, 20)}")
print(f"  r   = √2           = {nstr(r_alg, 20)}")
print()
print("PATH B (eigenvalue, through Born floor):")
print(f"  Q   = λ₃/λ₀        = {nstr(Q_ev, 20)}")
print(f"  θ   = K*·Δ₂/Δ₃     = {nstr(theta_ev, 20)}")
print(f"  r   = √(6Q-2)      = {nstr(r_ev, 20)}")
print()
print(f"  Q difference:  {float((Q_ev - Q_alg)/Q_alg * 100):.6f}%")
print(f"  θ difference:  {float((theta_ev - theta_alg)/theta_alg * 100):.6f}%")
print(f"  r difference:  {float((r_ev - r_alg)/r_alg * 100):.6f}%")
print()

# ============================================================
# KOIDE MASS FORMULA
# ============================================================

def koide_masses(r, theta, M_scale):
    """√m_k = M · (1 + r·cos(θ + 2πk/3)), k=0,1,2
    Compute all three, then sort by mass to identify τ, μ, e.
    """
    raw = []
    for k in range(3):
        alpha = theta + 2*pi*k/3
        sqrt_m = M_scale * (1 + r * cos(alpha))
        raw.append((k, sqrt_m**2))
    # Sort descending by mass: τ (heaviest), μ (middle), e (lightest)
    raw.sort(key=lambda x: -float(x[1]))
    return [x[1] for x in raw]  # [m_τ, m_μ, m_e]

# Experimental masses (MeV)
m_e_exp = mpf('0.51099895000')
m_mu_exp = mpf('105.6583755')
m_tau_exp = mpf('1776.86')

# For each path, determine M from the tau mass (or from m(0++))
# Actually, let's determine M by fitting to the sum of square roots
# (this uses ONE parameter: the overall scale)

# Method: fix M so that the tau mass is reproduced, then predict mu and e

def fit_and_predict(r, theta, label):
    """Fix M from tau mass, predict mu and e."""
    alpha_tau = theta  # k=0
    sqrt_m_tau_target = sqrt(m_tau_exp)
    M = sqrt_m_tau_target / (1 + r * cos(alpha_tau))

    masses = koide_masses(r, theta, M)
    m_tau, m_mu, m_e = masses

    print(f"  {label}:")
    print(f"    M (scale)  = {nstr(M, 15)} MeV^(1/2)")
    print()

    # Predictions
    print(f"    τ mass:  {float(m_tau):.4f} MeV  (input)")
    print(f"    μ mass:  {float(m_mu):.6f} MeV  (exp: {float(m_mu_exp):.6f})")
    print(f"    e mass:  {float(m_e):.8f} MeV  (exp: {float(m_e_exp):.8f})")
    print()

    mu_err = float((m_mu - m_mu_exp) / m_mu_exp * 100)
    e_err = float((m_e - m_e_exp) / m_e_exp * 100)

    print(f"    μ error: {mu_err:+.6f}%")
    print(f"    e error: {e_err:+.6f}%")
    print()

    # Koide ratio check
    sum_m = m_e + m_mu + m_tau
    sum_sqrt_m = sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau)
    Q_check = sum_m / sum_sqrt_m**2
    print(f"    Q check: {nstr(Q_check, 20)}")
    print()

    return m_e, m_mu, m_tau, mu_err, e_err

print("=" * 70)
print("MASS PREDICTIONS (τ mass as input)")
print("=" * 70)
print()

me_A, mmu_A, mtau_A, mu_err_A, e_err_A = fit_and_predict(r_alg, theta_alg, "PATH A (algebraic, Q=2/3, θ=2/9)")
me_B, mmu_B, mtau_B, mu_err_B, e_err_B = fit_and_predict(r_ev, theta_ev, "PATH B (eigenvalue, floor-corrected)")

print("=" * 70)
print("COMPARISON SUMMARY")
print("=" * 70)
print()
print(f"                    PATH A (no floor)    PATH B (floor)     Experiment")
print(f"  m_e (MeV)         {float(me_A):.8f}        {float(me_B):.8f}        {float(m_e_exp):.8f}")
print(f"  m_μ (MeV)         {float(mmu_A):.6f}          {float(mmu_B):.6f}          {float(m_mu_exp):.6f}")
print(f"  m_τ (MeV)         {float(mtau_A):.2f}            {float(mtau_B):.2f}            {float(m_tau_exp):.2f}")
print()
print(f"  μ error:          {mu_err_A:+.6f}%          {mu_err_B:+.6f}%")
print(f"  e error:          {e_err_A:+.6f}%          {e_err_B:+.6f}%")
print()

# Which is closer?
mu_A_abs = abs(mu_err_A)
mu_B_abs = abs(mu_err_B)
e_A_abs = abs(e_err_A)
e_B_abs = abs(e_err_B)

total_A = mu_A_abs + e_A_abs
total_B = mu_B_abs + e_B_abs

print(f"  Total |error|:    {total_A:.6f}%              {total_B:.6f}%")
print()
if total_B < total_A:
    print("  >>> PATH B (eigenvalue, through floor) gives BETTER predictions <<<")
elif total_A < total_B:
    print("  >>> PATH A (algebraic, no floor) gives BETTER predictions <<<")
else:
    print("  >>> Both paths give identical total error <<<")

print()
print("=" * 70)
print("THE KOIDE RATIO")
print("=" * 70)
print()

# Experimental Koide ratio from measured masses
sum_m_exp = m_e_exp + m_mu_exp + m_tau_exp
sum_sqrt_exp = sqrt(m_e_exp) + sqrt(m_mu_exp) + sqrt(m_tau_exp)
Q_exp = sum_m_exp / sum_sqrt_exp**2

print(f"  Q from experiment:     {nstr(Q_exp, 20)}")
print(f"  Q algebraic (2/3):     {nstr(Q_alg, 20)}")
print(f"  Q eigenvalue (λ₃/λ₀): {nstr(Q_ev, 20)}")
print()
print(f"  |Q_exp - 2/3|:         {float(fabs(Q_exp - Q_alg)):.10f}  ({float(fabs(Q_exp-Q_alg)/Q_alg*100):.6f}%)")
print(f"  |Q_exp - λ₃/λ₀|:      {float(fabs(Q_exp - Q_ev)):.10f}  ({float(fabs(Q_exp-Q_ev)/Q_ev*100):.6f}%)")
print()

if fabs(Q_exp - Q_ev) < fabs(Q_exp - Q_alg):
    print("  >>> Experimental Q is CLOSER to eigenvalue value <<<")
else:
    print("  >>> Experimental Q is CLOSER to algebraic 2/3 <<<")

print()

# ============================================================
# ALTERNATIVE: use BOTH eigenvalue parameters freely
# (don't force Brannen form, use eigenvalue ratios directly)
# ============================================================

print("=" * 70)
print("DIRECT EIGENVALUE MASS RATIOS (no Koide formula at all)")
print("=" * 70)
print()
print("If masses scale as eigenvalue powers, the dimensionless ratios are")
print("determined purely by the eigenvalue spectrum.")
print()

# The Koide formula says √m_k ∝ (1 + r·cos(θ+2πk/3))
# But more fundamentally, if the mass generation goes through the
# eigenvalue spectrum, the ratios might be:
#   m_τ/m_μ = f(λ₀, λ₁, λ₂, λ₃)
# where f is determined by the Z₃ generation structure

# In the paper's framework:
# m_k = M² · (1 + r·cos(θ_K + 2πk/3))²
# with M² = m(0++) · 11/60

# The mass scale from the paper
m_glueball = mpf('1712.5')  # MeV, predicted
M_sq = m_glueball * mpf(11) / mpf(60)

print(f"  m(0++) = {float(m_glueball)} MeV")
print(f"  M² = m(0++)·11/60 = {float(M_sq):.4f} MeV")
print(f"  M = √(M²) = {float(sqrt(M_sq)):.4f} MeV^(1/2)")
print()

# Path A masses (absolute, from M²)
M_A = sqrt(M_sq)
masses_A = koide_masses(r_alg, theta_alg, M_A)

# Path B masses (absolute, from M²)
M_B = sqrt(M_sq)
masses_B = koide_masses(r_ev, theta_ev, M_B)

print("  Absolute mass predictions from M² = m(0++)·11/60:")
print()
print(f"  {'':20s} PATH A        PATH B        Experiment")
print(f"  {'τ (MeV)':20s} {float(masses_A[0]):.2f}      {float(masses_B[0]):.2f}      {float(m_tau_exp):.2f}")
print(f"  {'μ (MeV)':20s} {float(masses_A[1]):.4f}      {float(masses_B[1]):.4f}      {float(m_mu_exp):.4f}")
print(f"  {'e (MeV)':20s} {float(masses_A[2]):.6f}      {float(masses_B[2]):.6f}      {float(m_e_exp):.6f}")
print()

for i, (name, m_exp) in enumerate([("τ", m_tau_exp), ("μ", m_mu_exp), ("e", m_e_exp)]):
    err_A = float((masses_A[i] - m_exp) / m_exp * 100)
    err_B = float((masses_B[i] - m_exp) / m_exp * 100)
    print(f"  {name} error:  A = {err_A:+.4f}%,  B = {err_B:+.4f}%")

print()
print("=" * 70)
print("DIMENSIONLESS RATIOS (the parameter-free test)")
print("=" * 70)
print()

# Pure ratios — no scale needed
ratio_tau_mu_A = masses_A[0] / masses_A[1]
ratio_tau_e_A = masses_A[0] / masses_A[2]
ratio_mu_e_A = masses_A[1] / masses_A[2]

ratio_tau_mu_B = masses_B[0] / masses_B[1]
ratio_tau_e_B = masses_B[0] / masses_B[2]
ratio_mu_e_B = masses_B[1] / masses_B[2]

ratio_tau_mu_exp = m_tau_exp / m_mu_exp
ratio_tau_e_exp = m_tau_exp / m_e_exp
ratio_mu_e_exp = m_mu_exp / m_e_exp

print(f"  {'Ratio':15s} {'PATH A':>15s} {'PATH B':>15s} {'Experiment':>15s}")
print(f"  {'-'*60}")
print(f"  {'m_τ/m_μ':15s} {float(ratio_tau_mu_A):15.6f} {float(ratio_tau_mu_B):15.6f} {float(ratio_tau_mu_exp):15.6f}")
print(f"  {'m_τ/m_e':15s} {float(ratio_tau_e_A):15.2f} {float(ratio_tau_e_B):15.2f} {float(ratio_tau_e_exp):15.2f}")
print(f"  {'m_μ/m_e':15s} {float(ratio_mu_e_A):15.4f} {float(ratio_mu_e_B):15.4f} {float(ratio_mu_e_exp):15.4f}")
print()

err_tm_A = float((ratio_tau_mu_A - ratio_tau_mu_exp)/ratio_tau_mu_exp * 100)
err_te_A = float((ratio_tau_e_A - ratio_tau_e_exp)/ratio_tau_e_exp * 100)
err_me_A = float((ratio_mu_e_A - ratio_mu_e_exp)/ratio_mu_e_exp * 100)

err_tm_B = float((ratio_tau_mu_B - ratio_tau_mu_exp)/ratio_tau_mu_exp * 100)
err_te_B = float((ratio_tau_e_B - ratio_tau_e_exp)/ratio_tau_e_exp * 100)
err_me_B = float((ratio_mu_e_B - ratio_mu_e_exp)/ratio_mu_e_exp * 100)

print(f"  Errors in ratios:")
print(f"  {'m_τ/m_μ':15s} {err_tm_A:+.6f}%     {err_tm_B:+.6f}%")
print(f"  {'m_τ/m_e':15s} {err_te_A:+.6f}%     {err_te_B:+.6f}%")
print(f"  {'m_μ/m_e':15s} {err_me_A:+.6f}%     {err_me_B:+.6f}%")
print()

total_ratio_A = abs(err_tm_A) + abs(err_te_A) + abs(err_me_A)
total_ratio_B = abs(err_tm_B) + abs(err_te_B) + abs(err_me_B)

print(f"  Total |ratio error|:  A = {total_ratio_A:.6f}%,  B = {total_ratio_B:.6f}%")
print()

if total_ratio_B < total_ratio_A:
    print("  >>> EIGENVALUE PATH (through Born floor) WINS on dimensionless ratios <<<")
elif total_ratio_A < total_ratio_B:
    print("  >>> ALGEBRAIC PATH (Q=2/3, θ=2/9) WINS on dimensionless ratios <<<")
else:
    print("  >>> TIE <<<")

print()
print("=" * 70)
