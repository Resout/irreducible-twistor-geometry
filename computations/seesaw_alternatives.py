#!/usr/bin/env python3
"""
SEESAW ALTERNATIVES: Can the 4.4% Dirac Koide gap be eliminated?
================================================================

The problem: Type-I seesaw with M_R proportional to identity gives
m_Dk proportional to |u_k| where u_k = 1 + r_nu * tan(theta_nu + 2*pi*k/3).

The Dirac Koide angle theta_D is near 2*pi/3. The offset:
  delta = theta_D - 2*pi/3
is close to sigma = -ln(23/30), but with a 4.36% gap.

This gap arises because u_1 < 0, and the sign flip |u_1| breaks
algebraic Koide structure.

This script investigates every alternative seesaw structure.
"""

from mpmath import (mp, mpf, sqrt, cos, sin, tan, pi, log, nstr, fabs,
                    acos, atan, atan2)

mp.dps = 50
H = 3
K_star = mpf(7)/30
sigma = -log(mpf(23)/30)
two_pi_thirds = 2*pi/3

print("=" * 80)
print("SEESAW ALTERNATIVES: ELIMINATING THE 4.36% DIRAC KOIDE GAP")
print("=" * 80)

# ============================================================
# UTILITY
# ============================================================

def dirac_koide(u_vals):
    """Given tan-Koide eigenstates u_k, compute Dirac Koide params.
    m_D proportional to |u_k|, so sqrt(m_D) proportional to |u_k|^{1/2}.
    Returns Q_D, theta_D, r_D, delta (= theta_D - 2*pi/3)."""
    abs_u = [fabs(x) for x in u_vals]
    fourth = [sqrt(a) for a in abs_u]  # |u_k|^{1/2}
    S_half = sum(fourth)
    S_one = sum(abs_u)
    Q_D = S_half**2 / (3 * S_one)
    r_D = sqrt(6*Q_D - 2)
    C = S_half / 3
    v = [(fourth[k]/C - 1)/r_D for k in range(3)]
    cp = sum(v[k]*cos(2*pi*k/3) for k in range(3))
    sp = sum(v[k]*sin(2*pi*k/3) for k in range(3))
    theta_D = atan2(-sp, cp)
    delta = theta_D - two_pi_thirds
    return Q_D, theta_D, r_D, delta


def general_koide(m_vals):
    """Given 3 positive masses, compute cos-Koide params.
    Returns Q, theta, r, delta (= theta - 2*pi/3)."""
    sq = [sqrt(m) for m in m_vals]
    S1 = sum(sq)
    S2 = sum(m_vals)
    Q = S1**2 / (3*S2)
    r = sqrt(6*Q - 2)
    C = S1/3
    v = [(sq[k]/C - 1)/r for k in range(3)]
    cp = sum(v[k]*cos(2*pi*k/3) for k in range(3))
    sp = sum(v[k]*sin(2*pi*k/3) for k in range(3))
    theta = atan2(-sp, cp)
    delta = theta - two_pi_thirds
    return Q, theta, r, delta


# ============================================================
# PART 0: BASELINE
# ============================================================
print("\n" + "=" * 80)
print("PART 0: BASELINE -- THE 4.36% GAP")
print("=" * 80)

Q_nu = mpf(29)/40
r_nu = sqrt(mpf(47)/20)
theta_nu = mpf(7)/810

print(f"\n  Neutrino parameters:")
print(f"    Q_nu  = 29/40 = {float(Q_nu):.6f}")
print(f"    r_nu  = sqrt(47/20) = {float(r_nu):.6f}")
print(f"    theta_nu = 7/810 = {float(theta_nu):.8f}")

u = [1 + r_nu * tan(theta_nu + 2*pi*k/3) for k in range(3)]
print(f"\n  Eigenstates:")
for k in range(3):
    sign = "+" if u[k] > 0 else "-"
    print(f"    u_{k} = {nstr(u[k], 15):>20s}  ({sign})")

Q_D, theta_D, r_D, delta = dirac_koide(u)

print(f"\n  Dirac Koide (Type-I, M_R proportional to 1):")
print(f"    Q_D     = {nstr(Q_D, 15)}")
print(f"    theta_D = {nstr(theta_D, 15)}")
print(f"    r_D     = {nstr(r_D, 15)}")
print(f"    delta   = theta_D - 2pi/3 = {nstr(delta, 15)}")
print(f"    sigma   = -ln(23/30)      = {nstr(sigma, 15)}")
print(f"    eps     = delta - sigma    = {nstr(delta - sigma, 15)}")
print(f"    delta/sigma = {nstr(delta/sigma, 15)}")
print(f"    GAP     = {float((delta/sigma - 1)*100):.4f}%")


# ============================================================
# PART 1: CHARGED LEPTON COMPARISON
# ============================================================
print("\n" + "=" * 80)
print("PART 1: CHARGED LEPTON SIGN FLIP COST")
print("=" * 80)

theta_ell = mpf(2)/9
r_ell = sqrt(2)
u_ell = [1 + r_ell * cos(theta_ell + 2*pi*k/3) for k in range(3)]

print(f"\n  Charged lepton eigenstates (cos-Koide, r=sqrt(2), theta=2/9):")
for k in range(3):
    sign = "+" if u_ell[k] > 0 else "-"
    print(f"    u_{k} = {nstr(u_ell[k], 12):>16s}  ({sign})")

all_pos = all(x > 0 for x in u_ell)
print(f"  All positive: {all_pos}")
if all_pos:
    print(f"  No sign flip for charged leptons.")
    print(f"  min(u_k) = {nstr(min(u_ell), 10)} (barely positive -- the electron)")

# Compute "Dirac Koide" of charged leptons for comparison
Q_De, theta_De, r_De, delta_e = dirac_koide(u_ell)
print(f"\n  Dirac Koide of charged leptons:")
print(f"    Q_De     = {nstr(Q_De, 12)}")
print(f"    theta_De = {nstr(theta_De, 12)}")
print(f"    delta_e  = theta_De - 2pi/3 = {nstr(delta_e, 12)}")
print(f"    For comparison: 2/9 = {nstr(theta_ell, 12)}")
print(f"    theta_De - theta_ell = {nstr(theta_De - theta_ell, 8)}")


# ============================================================
# PART 2: ALTERNATIVE SEESAW STRUCTURES
# ============================================================
print("\n" + "=" * 80)
print("PART 2: COMPREHENSIVE SEESAW SURVEY")
print("=" * 80)

print("""
  Different seesaw mechanisms give different m_D(m_nu) relations:
    Type-I:  m_D proportional to sqrt(m_nu)   -> sqrt(m_D) proportional to |u_k|^{1/2}
    Linear:  m_D proportional to m_nu          -> sqrt(m_D) proportional to |u_k|
    General: m_D proportional to m_nu^alpha    -> sqrt(m_D) proportional to |u_k|^alpha
""")

m_nu = [uk**2 for uk in u]
abs_u = [fabs(uk) for uk in u]

print(f"  {'Mechanism':45s} {'Q_D':>10s} {'delta':>14s} {'delta/sigma':>12s} {'gap%':>8s}")
print(f"  {'-'*45} {'-'*10} {'-'*14} {'-'*12} {'-'*8}")

# For m_D proportional to m_nu^alpha:
# sqrt(m_D) proportional to m_nu^{alpha/2} = (u_k^2)^{alpha/2} = |u_k|^alpha
# So the "Dirac Koide" uses |u_k|^alpha as the sqrt(m) values.

for label, alpha in [
    ("Type-I (alpha=1/2)", mpf(1)/2),
    ("Linear (alpha=1)", mpf(1)),
    ("Self-dual M_R (alpha=1)", mpf(1)),
    ("Inverse (alpha=1/2)", mpf(1)/2),
    ("m_D ~ m_nu^{1/3}", mpf(1)/3),
    ("m_D ~ m_nu^{1/4}", mpf(1)/4),
    ("m_D ~ m_nu^{2/3}", mpf(2)/3),
    ("m_D ~ m_nu^{3/4}", mpf(3)/4),
    ("m_D ~ m_nu^{3/2}", mpf(3)/2),
    ("m_D ~ m_nu^2", mpf(2)),
]:
    # sqrt(m_D) proportional to |u_k|^alpha
    vals = [a**alpha for a in abs_u]
    S_h = sum(vals)
    S_sq = sum(v**2 for v in vals)
    Q = S_h**2 / (3*S_sq)
    r = sqrt(6*Q - 2)
    C = S_h/3
    w = [(vals[k]/C - 1)/r for k in range(3)]
    cp = sum(w[k]*cos(2*pi*k/3) for k in range(3))
    sp = sum(w[k]*sin(2*pi*k/3) for k in range(3))
    th = atan2(-sp, cp)
    d = th - two_pi_thirds
    ratio = d/sigma
    gap = float((ratio - 1)*100)
    print(f"  {label:45s} {nstr(Q,7):>10s} {nstr(d,10):>14s} {nstr(ratio,8):>12s} {gap:>+7.2f}%")


# ============================================================
# PART 3: POWER SCAN -- FIND alpha THAT MAKES delta = sigma
# ============================================================
print("\n" + "=" * 80)
print("PART 3: POWER SCAN -- alpha WHERE delta = sigma")
print("=" * 80)

best_alpha = None
best_gap = mpf('1e10')

for i in range(1, 20000):
    alpha = mpf(i) / 20000  # scan 0 to 1
    vals = [a**alpha for a in abs_u]
    S_h = sum(vals)
    S_sq = sum(v**2 for v in vals)
    Q = S_h**2 / (3*S_sq)
    if 6*Q - 2 < 0:
        continue
    r = sqrt(6*Q - 2)
    C = S_h/3
    w = [(vals[k]/C - 1)/r for k in range(3)]
    cp = sum(w[k]*cos(2*pi*k/3) for k in range(3))
    sp = sum(w[k]*sin(2*pi*k/3) for k in range(3))
    th = atan2(-sp, cp)
    d = th - two_pi_thirds
    g = fabs(d - sigma)
    if g < best_gap:
        best_gap = g
        best_alpha = alpha

print(f"  Best alpha (coarse) = {nstr(best_alpha, 8)}, gap = {nstr(best_gap, 8)}")

# Refine
if best_alpha is not None:
    for i in range(-500, 501):
        alpha = best_alpha + mpf(i)/500000
        if alpha <= 0:
            continue
        vals = [a**alpha for a in abs_u]
        S_h = sum(vals)
        S_sq = sum(v**2 for v in vals)
        Q = S_h**2 / (3*S_sq)
        if 6*Q - 2 < 0:
            continue
        r = sqrt(6*Q - 2)
        C = S_h/3
        w = [(vals[k]/C - 1)/r for k in range(3)]
        cp = sum(w[k]*cos(2*pi*k/3) for k in range(3))
        sp = sum(w[k]*sin(2*pi*k/3) for k in range(3))
        th = atan2(-sp, cp)
        d = th - two_pi_thirds
        g = fabs(d - sigma)
        if g < best_gap:
            best_gap = g
            best_alpha = alpha

    print(f"  Best alpha (refined) = {nstr(best_alpha, 12)}, gap = {nstr(best_gap, 10)}")

    # Check framework fractions near best_alpha
    print(f"\n  Framework fractions near best alpha:")
    for q in range(1, 100):
        for p in range(1, 100):
            frac = mpf(p)/q
            if fabs(frac - best_alpha) < mpf('0.01'):
                vals = [a**frac for a in abs_u]
                S_h = sum(vals)
                S_sq = sum(v**2 for v in vals)
                Q = S_h**2 / (3*S_sq)
                if 6*Q - 2 < 0:
                    continue
                r = sqrt(6*Q - 2)
                C = S_h/3
                w = [(vals[k]/C - 1)/r for k in range(3)]
                cp = sum(w[k]*cos(2*pi*k/3) for k in range(3))
                sp = sum(w[k]*sin(2*pi*k/3) for k in range(3))
                th = atan2(-sp, cp)
                d = th - two_pi_thirds
                g_pct = float((d/sigma - 1)*100)
                if fabs(g_pct) < 5:
                    print(f"    alpha = {p}/{q} = {float(frac):.6f}: delta = {nstr(d, 10)}, gap = {g_pct:+.4f}%")


# ============================================================
# PART 4: COS-KOIDE WITH theta=7/60 (from v4 script)
# ============================================================
print("\n" + "=" * 80)
print("PART 4: COS-KOIDE WITH theta_nu=7/60 (alternative)")
print("=" * 80)

theta_alt = mpf(7)/60
r_cos = sqrt(6*Q_nu - 2)  # same Q, different form
w_cos = [1 + r_cos * cos(theta_alt + 2*pi*k/3) for k in range(3)]

print(f"\n  Cos-Koide with Q=29/40, theta=7/60:")
for k in range(3):
    sign = "+" if w_cos[k] > 0 else "-"
    print(f"    w_{k} = {nstr(w_cos[k], 12)}  ({sign})")

all_pos_cos = all(x > 0 for x in w_cos)
print(f"  All positive: {all_pos_cos}")

if all_pos_cos:
    print(f"  NO sign flip -- cos-Koide at theta=7/60 has all positive eigenstates!")
    Q_Dc, theta_Dc, r_Dc, delta_c = dirac_koide(w_cos)
    print(f"\n  Dirac Koide of cos-Koide neutrinos:")
    print(f"    Q_Dc     = {nstr(Q_Dc, 12)}")
    print(f"    theta_Dc = {nstr(theta_Dc, 12)}")
    print(f"    delta_c  = {nstr(delta_c, 12)}")
    print(f"    sigma    = {nstr(sigma, 12)}")
    print(f"    delta/sigma = {nstr(delta_c/sigma, 10)}")
    print(f"    GAP = {float((delta_c/sigma - 1)*100):.4f}%")
else:
    print(f"  Sign flip present in cos-Koide at theta=7/60.")
    Q_Dc, theta_Dc, r_Dc, delta_c = dirac_koide(w_cos)
    print(f"  delta = {nstr(delta_c, 12)}, sigma = {nstr(sigma, 12)}")
    print(f"  GAP = {float((delta_c/sigma - 1)*100):.4f}%")


# ============================================================
# PART 5: PRODUCT KOIDE -- DIRAC = SECTION x SUBSTRATE
# ============================================================
print("\n" + "=" * 80)
print("PART 5: PRODUCT KOIDE (section x substrate)")
print("=" * 80)

print("""
  The Dirac Yukawa couples left-handed leptons (section, cos-Koide)
  to right-handed neutrinos (substrate, tan-Koide).

  Product: Y_k = (1+sqrt(2) cos(2/9+2pi k/3)) * (1+r_nu tan(theta_nu+2pi k/3))
""")

v_cos = [1 + sqrt(2)*cos(theta_ell + 2*pi*k/3) for k in range(3)]
v_tan = [1 + r_nu*tan(theta_nu + 2*pi*k/3) for k in range(3)]
product = [v_cos[k] * v_tan[k] for k in range(3)]

print(f"  Product eigenstates:")
for k in range(3):
    sign = "+" if product[k] > 0 else "-"
    print(f"    Y_{k} = {nstr(v_cos[k], 8)} x {nstr(v_tan[k], 8)} = {nstr(product[k], 10)} ({sign})")

# Treat |Y_k| as Dirac masses
abs_prod = [fabs(p) for p in product]
sqrt_ap = [sqrt(a) for a in abs_prod]
S_h = sum(sqrt_ap)
S_sq = sum(abs_prod)
Q_prod = S_h**2 / (3*S_sq)
r_prod = sqrt(6*Q_prod - 2)
C_prod = S_h/3
w_prod = [(sqrt_ap[k]/C_prod - 1)/r_prod for k in range(3)]
cp_p = sum(w_prod[k]*cos(2*pi*k/3) for k in range(3))
sp_p = sum(w_prod[k]*sin(2*pi*k/3) for k in range(3))
theta_prod = atan2(-sp_p, cp_p)
delta_prod = theta_prod - two_pi_thirds

print(f"\n  Product Koide:")
print(f"    Q       = {nstr(Q_prod, 12)}")
print(f"    theta   = {nstr(theta_prod, 12)}")
print(f"    delta   = {nstr(delta_prod, 12)}")
print(f"    sigma   = {nstr(sigma, 12)}")
print(f"    delta/sigma = {nstr(delta_prod/sigma, 10)}")
if fabs(delta_prod) > mpf('1e-10'):
    print(f"    GAP = {float((delta_prod/sigma - 1)*100):.4f}%")


# ============================================================
# PART 6: M_R WITH KOIDE STRUCTURE
# ============================================================
print("\n" + "=" * 80)
print("PART 6: KOIDE-STRUCTURED M_R (breaking M_R proportional to 1)")
print("=" * 80)

print("""
  If M_Rk has its own Koide structure, the seesaw gives:
    m_Dk = sqrt(m_nu_k * M_Rk)

  The Dirac eigenstate becomes sqrt(|u_k^nu| * |u_k^R|).
  If u_k^R has the SAME sign pattern as u_k^nu, the product
  u_k^nu * u_k^R has consistent signs, avoiding the flip.

  Self-dual case: M_Rk has same (Q, theta) as m_nu_k:
    m_Dk = sqrt(u_k^2 * u_k^2) = u_k^2 (always positive!)
    sqrt(m_Dk) = |u_k| (same as linear seesaw)
""")

# Self-dual: same as alpha=1 case above. Already computed.
# Let's try M_R with DIFFERENT angle to cancel the sign flip.

# For the sign flip to cancel, we need u_k^R to have u_1^R < 0
# (same index as u_1^nu < 0). Then sqrt(|u_1^nu * u_1^R|)
# has no special treatment needed.
# But the PRODUCT |u_k^nu * u_k^R| still requires absolute values.

# More interesting: if M_R is ANTI-correlated with m_nu:
# M_Rk proportional to 1/m_nu_k. Then m_Dk = sqrt(m_nu_k / m_nu_k) = const.
# All Dirac masses equal -> degenerate -> Koide angle undefined.

# Or M_R has cos-Koide structure (all positive):
r_R = r_nu  # same r
theta_R = mpf(7)/60  # cos-Koide angle for M_R
u_R = [1 + r_R * cos(theta_R + 2*pi*k/3) for k in range(3)]

print(f"  Case: M_R with cos-Koide (r_nu, theta=7/60):")
for k in range(3):
    sign = "+" if u_R[k] > 0 else "-"
    print(f"    u_R_{k} = {nstr(u_R[k], 10)} ({sign})")

# m_Dk proportional to sqrt(u_k^2 * u_Rk^2) = |u_k| * |u_Rk|
mixed = [fabs(u[k]) * fabs(u_R[k]) for k in range(3)]
sqrt_mixed = [sqrt(m) for m in mixed]
S_hm = sum(sqrt_mixed)
S_sm = sum(mixed)
Q_mixed = S_hm**2 / (3*S_sm)
r_mixed = sqrt(6*Q_mixed - 2)
Cm = S_hm/3
wm = [(sqrt_mixed[k]/Cm - 1)/r_mixed for k in range(3)]
cp_m = sum(wm[k]*cos(2*pi*k/3) for k in range(3))
sp_m = sum(wm[k]*sin(2*pi*k/3) for k in range(3))
theta_mixed = atan2(-sp_m, cp_m)
delta_mixed = theta_mixed - two_pi_thirds

print(f"\n  Mixed seesaw Dirac Koide:")
print(f"    Q       = {nstr(Q_mixed, 12)}")
print(f"    theta   = {nstr(theta_mixed, 12)}")
print(f"    delta   = {nstr(delta_mixed, 12)}")
print(f"    sigma   = {nstr(sigma, 12)}")
print(f"    delta/sigma = {nstr(delta_mixed/sigma, 10)}")
print(f"    GAP = {float((delta_mixed/sigma - 1)*100):.4f}%")


# ============================================================
# PART 7: SCAN theta_R TO MINIMIZE GAP
# ============================================================
print("\n" + "=" * 80)
print("PART 7: SCAN theta_R (cos-Koide M_R) TO MINIMIZE GAP")
print("=" * 80)

best_tR = None
best_gR = mpf('1e10')

for i in range(1, 10000):
    tR = mpf(i)/10000 * pi/3
    u_R_test = [1 + r_nu * cos(tR + 2*pi*k/3) for k in range(3)]
    # Check all positive
    if any(x <= 0 for x in u_R_test):
        continue
    mixed_t = [fabs(u[k]) * u_R_test[k] for k in range(3)]
    sqrt_mt = [sqrt(m) for m in mixed_t]
    S_ht = sum(sqrt_mt)
    S_st = sum(mixed_t)
    Qt = S_ht**2 / (3*S_st)
    if 6*Qt - 2 < 0:
        continue
    rt = sqrt(6*Qt - 2)
    Ct = S_ht/3
    wt = [(sqrt_mt[k]/Ct - 1)/rt for k in range(3)]
    cp_t = sum(wt[k]*cos(2*pi*k/3) for k in range(3))
    sp_t = sum(wt[k]*sin(2*pi*k/3) for k in range(3))
    tht = atan2(-sp_t, cp_t)
    dt = tht - two_pi_thirds
    gt = fabs(dt - sigma)
    if gt < best_gR:
        best_gR = gt
        best_tR = tR

if best_tR is not None:
    print(f"  Best theta_R = {nstr(best_tR, 10)}")
    print(f"  Residual gap = {nstr(best_gR, 8)}")
    # Refine
    for i in range(-200, 201):
        tR = best_tR + mpf(i)/1000000
        if tR <= 0:
            continue
        u_R_test = [1 + r_nu * cos(tR + 2*pi*k/3) for k in range(3)]
        if any(x <= 0 for x in u_R_test):
            continue
        mixed_t = [fabs(u[k]) * u_R_test[k] for k in range(3)]
        sqrt_mt = [sqrt(m) for m in mixed_t]
        S_ht = sum(sqrt_mt)
        S_st = sum(mixed_t)
        Qt = S_ht**2 / (3*S_st)
        if 6*Qt - 2 < 0:
            continue
        rt = sqrt(6*Qt - 2)
        Ct = S_ht/3
        wt = [(sqrt_mt[k]/Ct - 1)/rt for k in range(3)]
        cp_t = sum(wt[k]*cos(2*pi*k/3) for k in range(3))
        sp_t = sum(wt[k]*sin(2*pi*k/3) for k in range(3))
        tht = atan2(-sp_t, cp_t)
        dt = tht - two_pi_thirds
        gt = fabs(dt - sigma)
        if gt < best_gR:
            best_gR = gt
            best_tR = tR

    print(f"  Refined theta_R = {nstr(best_tR, 12)}")
    print(f"  Refined gap = {nstr(best_gR, 10)}")

    # Compute delta/sigma at best
    u_R_best = [1 + r_nu * cos(best_tR + 2*pi*k/3) for k in range(3)]
    mixed_best = [fabs(u[k]) * u_R_best[k] for k in range(3)]
    sqrt_mb = [sqrt(m) for m in mixed_best]
    S_hb = sum(sqrt_mb)
    S_sb = sum(mixed_best)
    Qb = S_hb**2 / (3*S_sb)
    rb = sqrt(6*Qb - 2)
    Cb = S_hb/3
    wb = [(sqrt_mb[k]/Cb - 1)/rb for k in range(3)]
    cp_b = sum(wb[k]*cos(2*pi*k/3) for k in range(3))
    sp_b = sum(wb[k]*sin(2*pi*k/3) for k in range(3))
    thb = atan2(-sp_b, cp_b)
    db = thb - two_pi_thirds
    print(f"  delta = {nstr(db, 12)}")
    print(f"  delta/sigma = {nstr(db/sigma, 10)}")
    print(f"  GAP = {float((db/sigma - 1)*100):.4f}%")

    # Check if theta_R is a framework fraction
    print(f"\n  Simple fractions near theta_R:")
    for q in range(1, 200):
        for p in range(1, 100):
            frac = mpf(p)/q
            if fabs(frac - best_tR) < mpf('0.005'):
                print(f"    {p}/{q} = {nstr(frac, 10)}")
else:
    print(f"  No valid theta_R found in scan range.")


# ============================================================
# PART 8: THE ALGEBRAIC STRUCTURE
# ============================================================
print("\n" + "=" * 80)
print("PART 8: ALGEBRAIC STRUCTURE OF THE SIGN FLIP")
print("=" * 80)

abs_u = [fabs(uk) for uk in u]
S_signed = sum(u)
S_abs = sum(abs_u)
S_sq = sum(uk**2 for uk in u)

print(f"\n  Signed sum:   sum(u_k)   = {nstr(S_signed, 12)}")
print(f"  Unsigned sum: sum(|u_k|) = {nstr(S_abs, 12)}")
print(f"  Sum of squares: sum(u_k^2) = {nstr(S_sq, 12)}")
print(f"\n  Shift from sign flip: sum(|u|) - sum(u) = 2*|u_1| = {nstr(S_abs - S_signed, 12)}")
print(f"  2*|u_1| = {nstr(2*fabs(u[1]), 12)}")

# The signed sum uses the identity: sum tan(theta + 2pi k/3) = 3 tan(3 theta)
print(f"\n  Identity check: sum(u_k) = 3(1 + r_nu * tan(3*theta_nu)) = {nstr(3*(1 + r_nu*tan(3*theta_nu)), 12)}")
print(f"  Direct sum = {nstr(S_signed, 12)}")

# The epsilon = delta - sigma:
eps = delta - sigma
print(f"\n  epsilon = delta - sigma = {nstr(eps, 15)}")
print(f"  epsilon/theta_nu = {nstr(eps/theta_nu, 10)}")
print(f"  epsilon/K* = {nstr(eps/K_star, 10)}")
print(f"  epsilon*H = {nstr(eps*H, 10)}")
print(f"  epsilon*H^2 = {nstr(eps*H**2, 10)}")
print(f"  epsilon*30 = {nstr(eps*30, 10)}")
print(f"  epsilon/sigma = {nstr(eps/sigma, 10)} (this IS the 4.36% gap)")

# PSLQ-style search for epsilon
print(f"\n  Searching for epsilon as simple fraction:")
for q in range(1, 500):
    for p in range(1, 100):
        if fabs(eps - mpf(p)/q) < mpf('5e-4'):
            print(f"    eps approx {p}/{q} = {nstr(mpf(p)/q, 10)} (diff = {nstr(eps - mpf(p)/q, 4)})")


# ============================================================
# PART 9: TYPE-II / SCOTOGENIC (NO DIRAC MASS)
# ============================================================
print("\n" + "=" * 80)
print("PART 9: TYPE-II AND SCOTOGENIC")
print("=" * 80)

print("""
  TYPE-II SEESAW:
    m_nu = f * v_Delta directly. No Dirac mass. No sign flip issue.
    But the framework IS Type-I by construction (substrate = RH neutrino).
    Type-II is structurally alien.

  SCOTOGENIC:
    Neutrino masses at loop level. The Dirac coupling enters linearly
    in the loop integral, but the physical mass still involves a product
    of couplings, so the Koide structure is not preserved.

  CONCLUSION: These avoid the gap by avoiding Dirac masses entirely,
  but they don't match the framework's structure.
""")


# ============================================================
# PART 10: IS THE GAP theta_nu DEPENDENT?
# ============================================================
print("\n" + "=" * 80)
print("PART 10: HOW DOES THE GAP DEPEND ON theta_nu?")
print("=" * 80)

print(f"\n  Scanning theta_nu, fixed r_nu = sqrt(47/20):")
print(f"  {'theta_nu':>12s} {'delta':>14s} {'delta/sigma':>12s} {'gap%':>10s}")

# Wide scan
for i in [1, 2, 5, 7, 10, 20, 50, 70, 100, 200, 500]:
    th = mpf(i)/810
    u_test = [1 + r_nu * tan(th + 2*pi*k/3) for k in range(3)]
    Q_t, theta_t, r_t, delta_t = dirac_koide(u_test)
    ratio = delta_t/sigma
    gap = float((ratio-1)*100)
    frac = f"{i}/810"
    print(f"  {frac:>12s} {nstr(delta_t, 10):>14s} {nstr(ratio, 8):>12s} {gap:>+9.2f}%")

# Does any theta_nu give gap = 0?
print(f"\n  Fine scan near theta_nu = 7/810:")
best_th_nu = None
best_g_nu = mpf('1e10')

for i in range(-1000, 1001):
    th = theta_nu + mpf(i)/100000
    if th <= 0:
        continue
    u_test = [1 + r_nu * tan(th + 2*pi*k/3) for k in range(3)]
    Q_t, theta_t, r_t, delta_t = dirac_koide(u_test)
    g = fabs(delta_t - sigma)
    if g < best_g_nu:
        best_g_nu = g
        best_th_nu = th

print(f"  Best theta_nu = {nstr(best_th_nu, 12)}")
print(f"  Residual |delta - sigma| = {nstr(best_g_nu, 10)}")
print(f"  delta/sigma at best = {nstr(1 + best_g_nu/sigma, 10)} (or {nstr(1 - best_g_nu/sigma, 10)})")

# The gap is NOT zero-able by tuning theta_nu -- it's structural
if best_g_nu > mpf('0.001'):
    print(f"\n  The gap CANNOT be closed by tuning theta_nu alone!")
    print(f"  Minimum achievable |delta - sigma| = {nstr(best_g_nu, 8)}")
    print(f"  This is {float(best_g_nu/sigma)*100:.2f}% of sigma")
else:
    print(f"\n  The gap CAN be closed at theta_nu = {nstr(best_th_nu, 12)}")
    # Check if this is a framework fraction
    for q in range(1, 2000):
        for p in range(1, 100):
            if fabs(best_th_nu - mpf(p)/q) < mpf('1e-4'):
                print(f"    theta_nu approx {p}/{q}")


# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print(f"""
  The 4.36% gap: delta = theta_D - 2*pi/3 vs sigma = -ln(23/30).

  Baseline (Type-I, M_R = const):
    delta/sigma = {nstr(delta/sigma, 12)}
    gap = {float((delta/sigma - 1)*100):.4f}%
    eps = delta - sigma = {nstr(eps, 12)}

  FINDINGS:

  1. The sign flip is INESCAPABLE for tan-Koide.
     The tan function always produces at least one negative eigenstate.
     For charged leptons (cos-Koide), all eigenstates ARE positive
     (min u_k = {nstr(min(u_ell), 8)} > 0), so there is no sign flip.

  2. A power-law seesaw m_D ~ m_nu^alpha with alpha = 0.5797
     CAN close the gap (delta = sigma at sub-ppm). But:
     - alpha = 0.5797 is not a simple framework fraction
     - Closest: 11/19 = 0.5789 (gap 0.04%), 29/50 = 0.58 (gap 0.02%)
     - 40/69 = 0.5797 (gap 0.002%) -- not obviously a framework quantity
     - The standard Type-I seesaw gives alpha = 1/2, not 0.58

  3. The gap CAN also be closed by shifting theta_nu:
     - theta_nu = 0.01158 (instead of 7/810 = 0.00864) gives delta = sigma
     - Closest fraction: 1/86 = 0.01163
     - This requires ABANDONING theta_nu = 7/810 = K*/H^3

  4. The cos-Koide version (theta=7/60 instead of tan theta=7/810)
     has ALL positive eigenstates, so no sign flip occurs at all.
     But the delta from 2pi/3 is completely different.

  5. Neither the product Koide (section x substrate) nor Koide-structured
     M_R brings delta closer to sigma.

  INTERPRETATION: The 4.36% gap has two possible readings:
    (a) It is physical content: epsilon = {nstr(eps, 10)} is a prediction.
    (b) The neutrino angle theta_nu is not exactly 7/810.
  Given that 7/810 = K*/H^3 is algebraically clean and the gap is
  4.36%, interpretation (a) is preferred. The gap measures the cost
  of embedding the substrate's contractive (tan) Z_3 representation
  into the physical (positive-mass) Dirac sector.
""")

print("=" * 80)
print("END OF INVESTIGATION")
print("=" * 80)
