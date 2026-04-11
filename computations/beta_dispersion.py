#!/usr/bin/env python3
"""
Multi-site DS transfer operator: dispersion relation and momentum-dependent
running of coupling constants.

From the H=3 framework (Theorem thm:multisite), N sites on a 1D chain with
DS combination rule.  Linearising around the K*=7/30 fixed point in Fourier
space gives M_k = J_m + J_e cos(2πk/N).  Eigenvalues λ_α(k) encode the
dispersion relation;  spectral gaps Δ_α(k) = -ln|λ_α(k)| give decay rates.

Steps:
  1. Find the K*=7/30 fixed point at 40-digit precision.
  2. Compute J_m (state Jacobian) and J_e (evidence Jacobian) by central diffs.
  3. Dispersion relation for N=100 lattice.
  4. Extract running couplings and beta function.
  5. Compare to SM / MSSM one-loop coefficients.
  6. Prove H=3 uniqueness of β₀ = 11 = H²+H-1.
"""

from mpmath import (mp, mpf, matrix, sqrt, fabs, nstr, eig, chop,
                    findroot, eye, pi, cos, log, ln, re, im)

mp.dps = 40
H = 3
FLOOR = mpf(1) / mpf(H**3)

# ============================================================
# DS dynamics
# ============================================================

def ds_combine(m, e):
    s = m[:3]; theta = m[3]; ev = e[:3]; phi = e[3]
    s_pre = [s[i]*ev[i] + s[i]*phi + theta*ev[i] for i in range(3)]
    theta_pre = theta * phi
    total_pre = sum(s_pre) + theta_pre
    K = mpf(1) - total_pre
    if fabs(mpf(1) - K) < mpf(10)**(-35):
        return list(m), K
    denom = mpf(1) - K
    return [sp/denom for sp in s_pre] + [theta_pre/denom], K


def born_prob(m):
    L2sq = sum(x**2 for x in m)
    if L2sq < mpf(10)**(-60):
        return mpf(0)
    return m[3]**2 / L2sq


def enforce_floor(m):
    b = born_prob(m)
    if b >= FLOOR - mpf(10)**(-60):
        return list(m)
    S = sum(m[:3])
    if S < mpf(10)**(-60):
        return list(m)
    Sq = sum(x**2 for x in m[:3])
    # 26 = H³ - 1
    A_c = mpf(26) * S**2 - Sq
    B_c = mpf(2) * Sq
    C_c = -Sq
    disc = B_c**2 - 4*A_c*C_c
    if disc < 0:
        return list(m)
    t1 = (-B_c + sqrt(disc)) / (2*A_c)
    t2 = (-B_c - sqrt(disc)) / (2*A_c)
    cands = [t for t in [t1, t2] if mpf(0) < t < mpf(1)]
    if not cands:
        return list(m)
    t = min(cands, key=lambda x: fabs(x - m[3]))
    alpha_r = (mpf(1) - t) / S
    return [m[i]*alpha_r for i in range(3)] + [t]


def phi_map(m, e):
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds)


# ============================================================
# Step 1: Fixed point
# ============================================================

def eq_system(*p):
    s1, theta, w1, phi = p
    s2 = (mpf(1) - s1 - theta) / 2
    w2 = (mpf(1) - w1 - phi) / 2
    m = [s1, s2, s2, theta]
    e = [w1, w2, w2, phi]
    eq1 = theta**2 / (s1**2 + 2*s2**2 + theta**2) - FLOOR
    eq2 = phi**2 / (w1**2 + 2*w2**2 + phi**2) - FLOOR
    # K = sum of cross terms
    K_val = 2*s1*w2 + 2*s2*(w1 + w2)
    m_out = phi_map(m, e)
    return [eq1, eq2, K_val - mpf(7)/mpf(30), m_out[0] - s1]


print("="*70)
print("STEP 1: FINDING THE K*=7/30 FIXED POINT")
print("="*70)

sol = findroot(eq_system,
               [mpf('0.787'), mpf('0.155'), mpf('0.631'), mpf('0.128')])
s1s, ths, w1s, phis = sol
s2s = (mpf(1) - s1s - ths) / 2
w2s = (mpf(1) - w1s - phis) / 2
m_star = [s1s, s2s, s2s, ths]
e_star = [w1s, w2s, w2s, phis]

print(f"\nm* = ({nstr(s1s,15)}, {nstr(s2s,15)}, {nstr(s2s,15)}, {nstr(ths,15)})")
print(f"e* = ({nstr(w1s,15)}, {nstr(w2s,15)}, {nstr(w2s,15)}, {nstr(phis,15)})")

# Verify
m_check = phi_map(m_star, e_star)
resid = max(fabs(m_check[i] - m_star[i]) for i in range(4))
print(f"Fixed point residual: {nstr(resid, 6)}")

K_check = mpf(1) - sum(m_star[i]*e_star[i] for i in range(3)) - m_star[3]*e_star[3]
# Actually K = 1 - (s1*w1 + s2*w2 + s3*w3 + theta*phi)
# but Dempster K formula: K = sum over i≠j of products...
# Let's just compute it properly from cross terms
s_pre = [m_star[i]*e_star[i] + m_star[i]*e_star[3] + m_star[3]*e_star[i] for i in range(3)]
theta_pre = m_star[3]*e_star[3]
total_pre = sum(s_pre) + theta_pre
K_val = mpf(1) - total_pre
print(f"K at fixed point: {nstr(K_val, 15)}  (target: {nstr(mpf(7)/30, 15)})")
print(f"Born(m*): {nstr(born_prob(m_star), 15)}  (floor = {nstr(FLOOR, 15)})")


# ============================================================
# Step 2: Jacobians J_m and J_e by central finite differences
# ============================================================

print("\n" + "="*70)
print("STEP 2: JACOBIANS J_m AND J_e")
print("="*70)

eps = mpf(10)**(-20)

def compute_jacobian_m():
    """J_m[i,j] = d(Phi_i)/d(m_j) at (m*, e*)"""
    J = matrix(4, 4)
    for j in range(4):
        m_plus = list(m_star); m_plus[j] += eps
        m_minus = list(m_star); m_minus[j] -= eps
        out_p = phi_map(m_plus, e_star)
        out_m = phi_map(m_minus, e_star)
        for i in range(4):
            J[i, j] = (out_p[i] - out_m[i]) / (2*eps)
    return J

def compute_jacobian_e():
    """J_e[i,j] = d(Phi_i)/d(e_j) at (m*, e*)"""
    J = matrix(4, 4)
    for j in range(4):
        e_plus = list(e_star); e_plus[j] += eps
        e_minus = list(e_star); e_minus[j] -= eps
        out_p = phi_map(m_star, e_plus)
        out_m = phi_map(m_star, e_minus)
        for i in range(4):
            J[i, j] = (out_p[i] - out_m[i]) / (2*eps)
    return J

J_m = compute_jacobian_m()
J_e = compute_jacobian_e()

print("\nJ_m (state Jacobian):")
for i in range(4):
    row = [nstr(J_m[i,j], 10) for j in range(4)]
    print(f"  [{', '.join(row)}]")

print("\nJ_e (evidence Jacobian):")
for i in range(4):
    row = [nstr(J_e[i,j], 10) for j in range(4)]
    print(f"  [{', '.join(row)}]")

# On-site eigenvalues (k=0 check: M_0 = J_m + J_e)
M_0 = J_m + J_e
evals_0, _ = eig(M_0)
evals_0_real = sorted([re(chop(x)) for x in evals_0], key=lambda x: -fabs(x))
print(f"\nM_0 = J_m + J_e eigenvalues: {[nstr(x, 10) for x in evals_0_real]}")


# ============================================================
# Step 3: Dispersion relation
# ============================================================

print("\n" + "="*70)
print("STEP 3: DISPERSION RELATION (N=100 lattice)")
print("="*70)

N = 100
# We only need k = 0 .. N/2 by symmetry
n_modes = N // 2 + 1

# Store eigenvalues for each k
all_evals = []  # all_evals[k] = sorted list of 4 eigenvalues
all_gaps = []   # spectral gaps

for k in range(n_modes):
    cos_k = cos(2 * pi * k / N)
    M_k = J_m + J_e * cos_k
    evals_k, _ = eig(M_k)
    evals_sorted = sorted([chop(x) for x in evals_k], key=lambda x: -fabs(x))
    all_evals.append(evals_sorted)
    gaps = []
    for ev in evals_sorted:
        aev = fabs(ev)
        if aev > mpf(10)**(-30):
            gaps.append(-ln(aev))
        else:
            gaps.append(mpf('inf'))
    all_gaps.append(gaps)

print(f"\n{'k':>4s}  {'cos(2πk/N)':>12s}  {'|λ₀|':>12s}  {'|λ₁|':>12s}  {'|λ₂|':>12s}  {'|λ₃|':>12s}  {'Δ₀':>12s}  {'Δ₁':>12s}")
print("-" * 100)
for k in [0, 1, 2, 5, 10, 20, 25, 30, 40, 50]:
    cos_k = cos(2 * pi * k / N)
    evs = all_evals[k]
    gs = all_gaps[k]
    abs_evs = [fabs(e) for e in evs]
    print(f"{k:4d}  {nstr(float(cos_k), 6):>12s}  "
          f"{nstr(float(abs_evs[0]), 8):>12s}  "
          f"{nstr(float(abs_evs[1]), 8):>12s}  "
          f"{nstr(float(abs_evs[2]), 8):>12s}  "
          f"{nstr(float(abs_evs[3]), 8):>12s}  "
          f"{nstr(float(gs[0]), 8):>12s}  "
          f"{nstr(float(gs[1]), 8):>12s}")


# ============================================================
# Step 4: Running couplings and beta function
# ============================================================

print("\n" + "="*70)
print("STEP 4: RUNNING COUPLINGS AND BETA FUNCTION")
print("="*70)

# Option A: α(k) = |λ₀(k)|²  (dominant eigenvalue)
# Option B: α(k) = Tr(M_k²) / 4
# Option C: α(k) from the gauge mode eigenvalue

print("\n--- Option A: α(k) = |λ₀(k)|² ---")
alpha_A = []
for k in range(n_modes):
    a = fabs(all_evals[k][0])**2
    alpha_A.append(a)

print(f"  α(k=0)  = {nstr(float(alpha_A[0]), 10)}")
print(f"  α(k=25) = {nstr(float(alpha_A[25]), 10)}")
print(f"  α(k=50) = {nstr(float(alpha_A[50]), 10)}")

print("\n--- Option B: α(k) = Tr(M_k²) / (H+1) ---")
alpha_B = []
for k in range(n_modes):
    cos_k = cos(2 * pi * k / N)
    M_k = J_m + J_e * cos_k
    tr_sq = mpf(0)
    for i in range(4):
        for j in range(4):
            tr_sq += M_k[i,j] * M_k[j,i]
    # Actually Tr(M²) = sum_ij M_ij M_ji
    # But properly: Tr(M²) = sum_i (M²)_{ii} = sum_i sum_j M_ij M_ji
    alpha_B.append(tr_sq / (H + 1))
alpha_B_real = [re(x) for x in alpha_B]

print(f"  α(k=0)  = {nstr(float(alpha_B_real[0]), 10)}")
print(f"  α(k=25) = {nstr(float(alpha_B_real[25]), 10)}")
print(f"  α(k=50) = {nstr(float(alpha_B_real[50]), 10)}")

# Beta function: β = d(1/α)/d(ln k)
print("\n--- Beta function from Option A ---")
print(f"{'k':>4s}  {'α(k)':>14s}  {'1/α(k)':>14s}  {'β_eff':>14s}")
print("-" * 55)
for k in [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 49]:
    inv_a = 1 / alpha_A[k]
    if k > 1 and k < 50:
        inv_a_next = 1 / alpha_A[k+1]
        inv_a_prev = 1 / alpha_A[k-1]
        # Central difference for d(1/α)/d(ln k)
        beta_eff = (inv_a_next - inv_a_prev) / (ln(k+1) - ln(k-1))
        print(f"{k:4d}  {nstr(float(alpha_A[k]), 10):>14s}  "
              f"{nstr(float(inv_a), 10):>14s}  {nstr(float(beta_eff), 10):>14s}")
    else:
        print(f"{k:4d}  {nstr(float(alpha_A[k]), 10):>14s}  "
              f"{nstr(float(inv_a), 10):>14s}  {'---':>14s}")

print("\n--- Beta function from Option B ---")
print(f"{'k':>4s}  {'α(k)':>14s}  {'1/α(k)':>14s}  {'β_eff':>14s}")
print("-" * 55)
for k in [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 49]:
    inv_a = 1 / alpha_B_real[k]
    if k > 1 and k < 50:
        inv_a_next = 1 / alpha_B_real[k+1]
        inv_a_prev = 1 / alpha_B_real[k-1]
        beta_eff = (inv_a_next - inv_a_prev) / (ln(k+1) - ln(k-1))
        print(f"{k:4d}  {nstr(float(alpha_B_real[k]), 10):>14s}  "
              f"{nstr(float(inv_a), 10):>14s}  {nstr(float(beta_eff), 10):>14s}")
    else:
        print(f"{k:4d}  {nstr(float(alpha_B_real[k]), 10):>14s}  "
              f"{nstr(float(inv_a), 10):>14s}  {'---':>14s}")


# ============================================================
# Step 5: Compare to SM/MSSM
# ============================================================

print("\n" + "="*70)
print("STEP 5: COMPARISON TO SM / MSSM ONE-LOOP COEFFICIENTS")
print("="*70)

b3_SM = -7
b2_SM = mpf(-19)/6
b1_SM = mpf(41)/10

print(f"\nSM one-loop beta coefficients:")
print(f"  b₃ = {b3_SM}")
print(f"  b₂ = {nstr(float(b2_SM), 8)}")
print(f"  b₁ = {nstr(float(b1_SM), 8)}")

# Overall slope: fit 1/α vs ln(k) for mid-range k
# Use linear regression on k = 5..45
from mpmath import matrix as mpmatrix
sum_x = mpf(0); sum_y_A = mpf(0); sum_y_B = mpf(0)
sum_x2 = mpf(0); sum_xy_A = mpf(0); sum_xy_B = mpf(0)
n_pts = 0
for k in range(5, 46):
    x = ln(mpf(k))
    y_A = 1 / alpha_A[k]
    y_B = 1 / alpha_B_real[k]
    sum_x += x; sum_y_A += y_A; sum_y_B += y_B
    sum_x2 += x**2; sum_xy_A += x*y_A; sum_xy_B += x*y_B
    n_pts += 1

n = mpf(n_pts)
slope_A = (n*sum_xy_A - sum_x*sum_y_A) / (n*sum_x2 - sum_x**2)
slope_B = (n*sum_xy_B - sum_x*sum_y_B) / (n*sum_x2 - sum_x**2)

print(f"\nLinear fit slope of 1/α vs ln(k):")
print(f"  Option A (dominant eigenvalue²): slope = {nstr(float(slope_A), 10)}")
print(f"  Option B (Tr(M²)/(H+1)):        slope = {nstr(float(slope_B), 10)}")

# Normalize: in QCD, β₀ = -b₃ = 7 for SU(3) pure gauge gives d(1/α)/d(ln μ) = β₀/(2π)
beta0_QCD = 7  # pure gauge SU(3): 11*3/3 = 11, but b₃=-7 with quarks
print(f"\n  For comparison: β₀/(2π) = 11/(2π) = {nstr(float(11/(2*float(pi))), 8)} (pure gauge)")
print(f"  β₀/(2π) = 7/(2π)  = {nstr(float(7/(2*float(pi))), 8)} (with quarks)")

# Check ratio
if fabs(slope_A) > mpf(10)**(-20):
    ratio_A = slope_A / (11 / (2*pi))
    print(f"\n  slope_A / (11/(2π)) = {nstr(float(ratio_A), 10)}")
if fabs(slope_B) > mpf(10)**(-20):
    ratio_B = slope_B / (11 / (2*pi))
    print(f"  slope_B / (11/(2π)) = {nstr(float(ratio_B), 10)}")

# Per-eigenvalue slopes
print("\n--- Per-eigenvalue running ---")
for alpha_idx in range(4):
    alpha_mode = []
    for k in range(n_modes):
        alpha_mode.append(fabs(all_evals[k][alpha_idx])**2)

    sx = mpf(0); sy = mpf(0); sx2 = mpf(0); sxy = mpf(0); nn = 0
    for k in range(5, 46):
        if alpha_mode[k] > mpf(10)**(-30):
            x = ln(mpf(k))
            y = 1 / alpha_mode[k]
            sx += x; sy += y; sx2 += x**2; sxy += x*y; nn += 1
    if nn > 2:
        nf = mpf(nn)
        sl = (nf*sxy - sx*sy) / (nf*sx2 - sx**2)
        print(f"  Mode α={alpha_idx}: slope = {nstr(float(sl), 10)}, "
              f"ratio to 11/(2π) = {nstr(float(sl / (11/(2*pi))), 8)}")

# MSSM check
print(f"\nMSSM b₃ = -3.  DS b_eff contribution: see per-mode slopes above.")
print(f"SM-to-MSSM shift: Δb₃ = -3 - (-7) = +4 from superpartners.")


# ============================================================
# Step 6: H-specificity of β₀ = 11
# ============================================================

print("\n" + "="*70)
print("STEP 6: H=3 UNIQUENESS -- β₀ = 11 = H² + H - 1")
print("="*70)

print(f"\nPure SU(N) one-loop: β₀ = 11N/3.")
print(f"At N=H: β₀ = 11H/3.")
print(f"Massive gauge boson count at H: H²+H-1 (from framework).")
print(f"\nCondition: 11H/3 = H²+H-1")
print(f"  => 11H = 3H²+3H-3")
print(f"  => 3H²-8H-3 = 0")
print(f"  => (3H+1)(H-3) = 0")
print(f"  => H = 3  or  H = -1/3")
print(f"\nH=3 is the UNIQUE positive integer solution.")

print(f"\nVerification for H=1..10:")
print(f"{'H':>3s}  {'11H/3':>10s}  {'H²+H-1':>10s}  {'Match':>7s}")
print("-" * 38)
for h in range(1, 11):
    beta0 = mpf(11*h) / 3
    count = h**2 + h - 1
    match = "YES" if fabs(beta0 - count) < mpf(10)**(-10) else "no"
    print(f"{h:3d}  {nstr(float(beta0), 8):>10s}  {count:>10d}  {match:>7s}")

# Additional: verify the polynomial identity
print(f"\n3H²-8H-3 = 0 at H=3: {3*9 - 8*3 - 3} (should be 0)")
print(f"(3H+1)(H-3) at H=3: {(3*3+1)*(3-3)} (should be 0)")

# Physical interpretation
print(f"\nPhysical content:")
print(f"  At H=3: β₀ = 11*3/3 = 11")
print(f"  H²+H-1 = 9+3-1 = 11")
print(f"  This counts: dim(adjoint SU(3)) - dim(Cartan) + dim(S²) - 1")
print(f"            = 8 - 2 + 6 - 1 = 11  [gauge boson decomposition]")
print(f"  Or simply: H(H+1) - 1 = 3·4 - 1 = 11")

# ============================================================
# Step 7: IR and UV limits
# ============================================================

print("\n" + "="*70)
print("STEP 7: IR vs UV LIMITS")
print("="*70)

M_IR = J_m + J_e   # k=0
M_UV = J_m - J_e   # k=N/2

evals_IR, _ = eig(M_IR)
evals_UV, _ = eig(M_UV)

evals_IR_s = sorted([chop(x) for x in evals_IR], key=lambda x: -fabs(x))
evals_UV_s = sorted([chop(x) for x in evals_UV], key=lambda x: -fabs(x))

print(f"\nIR (k=0):  eigenvalues = {[nstr(float(re(x)), 10) for x in evals_IR_s]}")
print(f"UV (k=N/2): eigenvalues = {[nstr(float(re(x)), 10) for x in evals_UV_s]}")

print(f"\nIR spectral gaps:")
for i, ev in enumerate(evals_IR_s):
    aev = fabs(ev)
    if aev > mpf(10)**(-30):
        print(f"  Δ_{i}(IR) = {nstr(float(-ln(aev)), 10)}")

print(f"\nUV spectral gaps:")
for i, ev in enumerate(evals_UV_s):
    aev = fabs(ev)
    if aev > mpf(10)**(-30):
        print(f"  Δ_{i}(UV) = {nstr(float(-ln(aev)), 10)}")

# Gap between IR and UV for each mode
print(f"\nUV-IR gap shift per mode:")
for i in range(4):
    a_IR = fabs(evals_IR_s[i])
    a_UV = fabs(evals_UV_s[i])
    if a_IR > mpf(10)**(-30) and a_UV > mpf(10)**(-30):
        delta = -ln(a_UV) + ln(a_IR)
        print(f"  Mode {i}: Δ(UV) - Δ(IR) = {nstr(float(delta), 10)}")

print("\n" + "="*70)
print("COMPUTATION COMPLETE")
print("="*70)
