"""
RANK-1 FACTORISATION VERIFICATION
===================================

Verify that dbar(F) = column_vector * row_vector
where column = (-s1/S, -s2/S, -s3/S, 1)^T
and row = (dt/dz_bar_1, ..., dt/dz_bar_4)

Requirements: mpmath
"""

from mpmath import (mp, mpf, mpc, sqrt, matrix, eye, nstr, fabs, log)

PRECISION = 60
mp.dps = PRECISION
ONE = mpf(1); ZERO = mpf(0)
FLOOR = ONE / mpf(27)

sigma1 = matrix([[0, 1], [1, 0]])
sigma2 = matrix([[0, mpc(0,-1)], [mpc(0,1), 0]])
sigma3 = matrix([[1, 0], [0, -1]])
sq2 = sqrt(mpf(2))
I2 = eye(2)

def ds_combine(m, e):
    s1, s2, s3, th = m; e1, e2, e3, ph = e
    sn = [s1*e1+s1*ph+th*e1, s2*e2+s2*ph+th*e2, s3*e3+s3*ph+th*e3]
    tn = th*ph
    K = s1*e2+s1*e3+s2*e1+s2*e3+s3*e1+s3*e2
    d = ONE-K
    return [sn[0]/d, sn[1]/d, sn[2]/d, tn/d], K

def complex_floor_enforce(m):
    s1, s2, s3, th = m
    s_mag_sq = fabs(s1)**2 + fabs(s2)**2 + fabs(s3)**2
    th_mag = fabs(th)
    total_mag_sq = s_mag_sq + th_mag**2
    if total_mag_sq < mpf(10)**(-80): return list(m)
    born = th_mag**2 / total_mag_sq
    if born >= FLOOR - mpf(10)**(-30): return list(m)
    S = s1 + s2 + s3
    S_mag_sq = fabs(S)**2
    if S_mag_sq < mpf(10)**(-80): return list(m)
    phase = th / th_mag if th_mag > mpf(10)**(-80) else ONE
    R = s_mag_sq / S_mag_sq
    Re_ph = phase.real if hasattr(phase, 'real') else mpf(phase)
    ac = mpf(26) - R
    bc = mpf(2) * R * Re_ph
    cc = -R
    disc = bc**2 - 4*ac*cc
    r_new = (-bc + sqrt(disc)) / (2*ac)
    th_new = phase * r_new
    alpha = (ONE - th_new) / S
    return [s1*alpha, s2*alpha, s3*alpha, th_new]

def full_step(m, e):
    m_ds, K = ds_combine(m, e)
    return complex_floor_enforce(m_ds), K

def make_evidence(p_dom):
    pw = (ONE-p_dom)/2; sc = ONE-FLOOR
    raw = [sqrt(p_dom*sc), sqrt(pw*sc), sqrt(pw*sc), sqrt(FLOOR)]
    tot = sum(raw)
    return [r/tot for r in raw]

# Find equilibrium
e_star = make_evidence(mpf('0.932'))
m = [mpf('0.4'), mpf('0.15'), mpf('0.15'), mpf('0.3')]
for _ in range(100000):
    m2, _ = full_step(m, e_star)
    if max(fabs(m2[k]-m[k]) for k in range(4)) < mpf(10)**(-(PRECISION-10)):
        m = m2; break
    m = m2
m_star = m

print(f"Equilibrium: [{', '.join(nstr(x, 10) for x in m_star)}]")

# Compute dbar(Phi) at m* via Wirtinger finite differences
eps = mpf(10)**(-20)

J_anti = matrix(4, 4)
for alpha_idx in range(4):
    m_pr = list(m_star); m_mr = list(m_star)
    m_pi = list(m_star); m_mi = list(m_star)
    m_pr[alpha_idx] += eps
    m_mr[alpha_idx] -= eps
    m_pi[alpha_idx] += mpc(0, eps)
    m_mi[alpha_idx] -= mpc(0, eps)

    out_pr, _ = full_step(m_pr, e_star)
    out_mr, _ = full_step(m_mr, e_star)
    out_pi, _ = full_step(m_pi, e_star)
    out_mi, _ = full_step(m_mi, e_star)

    for j in range(4):
        d_re = (out_pr[j] - out_mr[j]) / (2*eps)
        d_im = (out_pi[j] - out_mi[j]) / (2*eps)
        J_anti[j, alpha_idx] = (d_re + mpc(0,1)*d_im) / 2

print(f"\ndbar(Phi) matrix:")
for i in range(4):
    row = [nstr(J_anti[i,j], 6) for j in range(4)]
    print(f"  [{', '.join(row)}]")

# Check factorisation: every row should be proportional to row 4
# column vector should be (-s1/S, -s2/S, -s3/S, 1)
S_val = m_star[0] + m_star[1] + m_star[2]
expected_col = [-m_star[0]/S_val, -m_star[1]/S_val, -m_star[2]/S_val, ONE]

print(f"\nExpected column vector (-s_i/S, 1):")
print(f"  [{', '.join(nstr(x, 10) for x in expected_col)}]")

# Extract actual column ratios from J_anti
# If J_anti = col * row, then J_anti[i,j] / J_anti[3,j] = col[i] / col[3] = col[i]
print(f"\nActual row ratios J_anti[i,:] / J_anti[3,:]:")
for i in range(4):
    # Use column 0 for the ratio
    if fabs(J_anti[3, 0]) > mpf(10)**(-40):
        ratio = J_anti[i, 0] / J_anti[3, 0]
        print(f"  row {i} / row 3 = {nstr(ratio, 10)}  (expected {nstr(expected_col[i], 10)})")

# Check all columns give same ratio
print(f"\nConsistency check across columns:")
for j in range(4):
    if fabs(J_anti[3, j]) > mpf(10)**(-40):
        ratios = [J_anti[i, j] / J_anti[3, j] for i in range(3)]
        expected = [-m_star[i]/S_val for i in range(3)]
        errors = [fabs(ratios[i] - expected[i]) for i in range(3)]
        max_err = max(errors)
        print(f"  col {j}: max ratio error = {nstr(max_err, 5)}")

# Final: compute rank via SVD proxy
# J^H J eigenvalues
JHJ = matrix(4, 4)
for i in range(4):
    for j in range(4):
        val = mpc(0, 0)
        for k in range(4):
            val += J_anti[k,i].conjugate() * J_anti[k,j]
        JHJ[i,j] = val

# Power iteration for top 2 singular values
def power_iter(M, n=300):
    v = matrix(4, 1)
    for i in range(4):
        v[i,0] = mpc(mpf(1)/(i+1), mpf(1)/(i+2))
    for _ in range(n):
        w = M * v
        norm = sqrt(sum(fabs(w[i,0])**2 for i in range(4)))
        v = w * (ONE/norm)
    Mv = M * v
    lam = sum(v[i,0].conjugate() * Mv[i,0] for i in range(4)).real
    return lam, v

sv1_sq, v1 = power_iter(JHJ)
# Deflate
JHJ2 = matrix(4, 4)
for i in range(4):
    for j in range(4):
        JHJ2[i,j] = JHJ[i,j] - sv1_sq * v1[i,0] * v1[j,0].conjugate()
sv2_sq, _ = power_iter(JHJ2)

sv1 = sqrt(fabs(sv1_sq))
sv2 = sqrt(fabs(sv2_sq))
print(f"\nSingular values:")
print(f"  sigma_1 = {nstr(sv1, 12)}")
print(f"  sigma_2 = {nstr(sv2, 12)}")
print(f"  sigma_2/sigma_1 = {nstr(sv2/sv1, 8)}")
if sv2/sv1 < mpf(10)**(-10):
    print(f"  RANK 1 CONFIRMED (sigma_2/sigma_1 = 10^{nstr(log(sv2/sv1)/log(mpf(10)), 4)})")
