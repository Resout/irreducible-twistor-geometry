"""
Ward F+ coefficient at high precision.

Goal: determine whether |F+|^2 * det(M*)^2 = 1/H = 1/3 exactly,
and whether |F+|^2 * K* = 7/8 exactly (or some other clean fraction).

Uses mpmath for 100-digit precision throughout.
"""

from mpmath import mp, mpf, sqrt, matrix, nstr, fabs, eye

mp.dps = 100  # 100 decimal digits

ONE = mpf(1)
H = mpf(3)
FLOOR = ONE / H**3
TARGET_K = mpf(7) / mpf(30)

# Pauli matrices (as mpmath matrices)
sigma1 = matrix([[0, 1], [1, 0]])
sigma2 = matrix([[0, -1j], [1j, 0]])
sigma3 = matrix([[1, 0], [0, -1]])
I2 = eye(2)
basis = [sigma1, sigma2, sigma3, I2]

def to_matrix(m):
    s1, s2, s3, th = m
    return (th * I2 + s1 * sigma1 + s2 * sigma2 + s3 * sigma3) / sqrt(2)

def enforce_floor(m):
    s1, s2, s3, th = m
    ssq = s1**2 + s2**2 + s3**2
    born = th**2 / (ssq + th**2)
    if born >= FLOOR:
        return [s1, s2, s3, th]
    ss = s1 + s2 + s3
    sq = ssq
    r = sq / ss**2
    disc = (2*r)**2 + 4*(26 - r)*r
    t = (-2*r + sqrt(disc)) / (2*(26 - r))
    alpha = (ONE - t) / ss
    return [s1*alpha, s2*alpha, s3*alpha, t]

def ds_combine(m, e):
    s1, s2, s3, th = m
    e1, e2, e3, ph = e
    sn1 = s1*e1 + s1*ph + th*e1
    sn2 = s2*e2 + s2*ph + th*e2
    sn3 = s3*e3 + s3*ph + th*e3
    tn = th * ph
    K = s1*e2 + s1*e3 + s2*e1 + s2*e3 + s3*e1 + s3*e2
    d = ONE - K
    return [sn1/d, sn2/d, sn3/d, tn/d], K

def full_step(m, e):
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds), K

def make_evidence(p_dom):
    p_weak = (ONE - p_dom) / 2
    sc = ONE - FLOOR
    raw = [sqrt(p_dom * sc), sqrt(p_weak * sc), sqrt(p_weak * sc), sqrt(FLOOR)]
    tot = sum(raw)
    return [r / tot for r in raw]

print("="*70)
print("HIGH-PRECISION WARD F+ COEFFICIENT")
print(f"Working precision: {mp.dps} digits")
print("="*70)

# Secant method to find p_dom
print("\nFinding equilibrium via secant method...")

p0, p1 = mpf('0.932'), mpf('0.933')

def K_at_pdom(p):
    e = make_evidence(p)
    m = [mpf('0.4'), mpf('0.15'), mpf('0.15'), mpf('0.3')]
    for _ in range(10000):
        m, _ = full_step(m, e)
    _, K = ds_combine(m, e)
    return K

K0 = K_at_pdom(p0)
K1 = K_at_pdom(p1)
tol = mpf(10)**(-(mp.dps - 20))

for step in range(100):
    f0 = K0 - TARGET_K
    f1 = K1 - TARGET_K
    if fabs(f1 - f0) < mpf(10)**(-mp.dps):
        break
    p2 = p1 - f1 * (p1 - p0) / (f1 - f0)
    p0, K0 = p1, K1
    p1 = p2
    K1 = K_at_pdom(p1)
    if fabs(p1 - p0) < tol:
        break
    if step % 3 == 0:
        print(f"  Step {step}: |dp| = {nstr(fabs(p1-p0), 5)}")

p_dom = p1
e_star = make_evidence(p_dom)
m = [mpf('0.4'), mpf('0.15'), mpf('0.15'), mpf('0.3')]
for _ in range(10000):
    m, _ = full_step(m, e_star)
m_star = m
m_ds, K_star = ds_combine(m_star, e_star)

print(f"\n|K* - 7/30| = {nstr(fabs(K_star - TARGET_K), 10)}")

# Matrices
M_star = to_matrix(m_star)
M_ds = to_matrix(m_ds)
M_ds_inv = M_ds**(-1)
epsilon = M_star - M_ds

det_M_star = M_star[0,0]*M_star[1,1] - M_star[0,1]*M_star[1,0]
det_M_ds = M_ds[0,0]*M_ds[1,1] - M_ds[0,1]*M_ds[1,0]

print(f"det(M*) = {nstr(det_M_star, 30)}")
print(f"det(M_ds) = {nstr(det_M_ds, 30)}")
print(f"||epsilon|| = {nstr(sqrt(sum(fabs(epsilon[i,j])**2 for i in range(2) for j in range(2))), 20)}")

# Analytical Jacobians (same formulas, now at high precision)
s1, s2, s3, th = m_star
e1, e2, e3, ph = e_star
oneMinusK = ONE - K_star

S_vals = [s1*(e1+ph)+th*e1, s2*(e2+ph)+th*e2, s3*(e3+ph)+th*e3, th*ph]
N = [S_vals[i] / oneMinusK for i in range(4)]

# DS Jacobian
dSdm = [[mpf(0)]*4 for _ in range(4)]
for i in range(3):
    dSdm[i][i] = e_star[i] + ph
    dSdm[i][3] = e_star[i]
dSdm[3][3] = ph
dKdm = [e2+e3, e1+e3, e1+e2, mpf(0)]

J_DS = matrix(4, 4)
for i in range(4):
    for j in range(4):
        J_DS[i,j] = (dSdm[i][j] * oneMinusK + S_vals[i] * dKdm[j]) / oneMinusK**2

# Floor Jacobian
Ss = N[0]+N[1]+N[2]
Sq = N[0]**2+N[1]**2+N[2]**2
R = Sq/Ss**2
u = sqrt(mpf(26)*R)
t_val = (u - R)/(mpf(26) - R)
g = (ONE - t_val)/Ss

dtdR = (mpf(338) + mpf(13)*R - mpf(26)*u) / (u * (mpf(26) - R)**2)
dRdN = [2*(N[i]*Ss - Sq)/Ss**3 for i in range(3)]
dtdN = [dtdR * dRdN[i] for i in range(3)]
dgdN = [(-dtdN[i]*Ss - (ONE-t_val))/Ss**2 for i in range(3)]

J_floor = matrix(4, 4)
for i in range(3):
    for j in range(3):
        J_floor[i,j] = (g if i==j else mpf(0)) + N[i]*dgdN[j]
J_floor[3,0] = dtdN[0]; J_floor[3,1] = dtdN[1]; J_floor[3,2] = dtdN[2]; J_floor[3,3] = mpf(0)

J_correction = (J_floor - eye(4)) * J_DS

# Ward connection basis B[mu]
B = [matrix(2,2) for _ in range(4)]
for mu in range(4):
    for j in range(4):
        B_contrib = M_ds_inv * basis[j] * (J_correction[j, mu] / sqrt(2))
        B[mu] = B[mu] + B_contrib

# Incidence relation coefficients
sq2 = sqrt(2)
A1 = [mpf(0)]*4; C1 = [mpf(0)]*4
A2 = [mpf(0)]*4; C2 = [mpf(0)]*4

A1[0] = -1j/sq2; A1[3] = -1j/sq2
C1[1] = -1j/sq2; C1[2] = ONE/sq2
A2[1] = -1j/sq2; A2[2] = ONE/sq2
C2[0] = -1j/sq2; C2[3] = 1j/sq2

# Laurent coefficients
rho_m1 = matrix(2, 2)
for mu in range(4):
    coeff = C1[mu] + C2[mu]
    rho_m1 = rho_m1 + B[mu] * coeff

# F+ = rho_{-1} (the residue)
F_plus = rho_m1

# |F+|^2 = Tr(F+^dag F+)
F_plus_sq = mpf(0)
for i in range(2):
    for j in range(2):
        F_plus_sq += fabs(F_plus[i,j])**2

det_M_star_sq = fabs(det_M_star)**2

ratio_1 = F_plus_sq * det_M_star_sq
ratio_2 = F_plus_sq * K_star
ratio_3 = F_plus_sq * det_M_star_sq * H

print(f"\n{'='*70}")
print("HIGH-PRECISION RATIOS")
print(f"{'='*70}")
print(f"\n|F+|^2 = {nstr(F_plus_sq, 40)}")
print(f"det(M*)^2 = {nstr(det_M_star_sq, 40)}")
print(f"\n|F+|^2 * det(M*)^2 = {nstr(ratio_1, 40)}")
print(f"1/3                = {nstr(ONE/3, 40)}")
print(f"Difference         = {nstr(ratio_1 - ONE/3, 15)}")
print(f"Relative error     = {nstr(fabs(ratio_1 - ONE/3)/(ONE/3), 10)}")

print(f"\n|F+|^2 * det(M*)^2 * H = {nstr(ratio_3, 40)}")
print(f"1                       = 1")
print(f"Difference              = {nstr(ratio_3 - ONE, 15)}")

print(f"\n|F+|^2 * K* = {nstr(ratio_2, 40)}")
print(f"7/8         = {nstr(mpf(7)/8, 40)}")
print(f"Difference  = {nstr(ratio_2 - mpf(7)/8, 15)}")

# Check other simple fractions near ratio_1
print(f"\nNearby fractions for |F+|^2 * det(M*)^2 = {nstr(ratio_1, 20)}:")
for p in range(1, 50):
    for q in range(1, 100):
        frac = mpf(p)/mpf(q)
        if fabs(ratio_1 - frac) < mpf('0.001'):
            print(f"  {p}/{q} = {nstr(frac, 15)}, error = {nstr(fabs(ratio_1 - frac), 10)}")

print(f"\nNearby fractions for |F+|^2 * K* = {nstr(ratio_2, 20)}:")
for p in range(1, 50):
    for q in range(1, 100):
        frac = mpf(p)/mpf(q)
        if fabs(ratio_2 - frac) < mpf('0.005'):
            print(f"  {p}/{q} = {nstr(frac, 15)}, error = {nstr(fabs(ratio_2 - frac), 10)}")

# su(2) decomposition
tr_F = (F_plus[0,0] + F_plus[1,1]) / 2
su2_1 = (F_plus[0,1] + F_plus[1,0]) / 2
su2_2 = (F_plus[0,1] - F_plus[1,0]) / (2j)
su2_3 = (F_plus[0,0] - F_plus[1,1]) / 2

scalar_sq = fabs(tr_F)**2
gauge_sq = fabs(su2_1)**2 + fabs(su2_2)**2 + fabs(su2_3)**2
total = scalar_sq + gauge_sq

print(f"\nsu(2) decomposition:")
print(f"  Scalar fraction: {nstr(scalar_sq/total, 15)}")
print(f"  Gauge fraction:  {nstr(gauge_sq/total, 15)}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
