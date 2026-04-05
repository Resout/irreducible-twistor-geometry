"""
Minimal polynomial of tr(J) — the trace of the projected Jacobian at K*=7/30.

If [Q(tr):Q] = 19, then the degree-19 comes from the TRACE, and λ₀ ∈ Q(tr).
If [Q(tr):Q] ≠ 19, the degree comes from the eigenvalue extraction step.

Uses 300 digits and maxcoeff 10^15 to handle large coefficients.
"""
import sys
try:
    from mpmath import (mp, mpf, sqrt, log, matrix, nstr, fabs,
                        pslq, chop, polyroots)
except ImportError:
    print("mpmath required"); sys.exit(1)

mp.dps = 300
ONE = mpf(1); FLOOR = ONE/mpf(27); K_STAR = mpf(7)/mpf(30)

def ds_mp(m, e):
    s1,s2,s3,th = m; e1,e2,e3,ph = e
    sn = [s1*e1+s1*ph+th*e1, s2*e2+s2*ph+th*e2, s3*e3+s3*ph+th*e3]
    tn = th*ph
    K = s1*e2+s1*e3+s2*e1+s2*e3+s3*e1+s3*e2
    d = ONE-K
    return [sn[0]/d, sn[1]/d, sn[2]/d, tn/d], K

def floor_mp(m):
    s1,s2,s3,th = m
    ssq = s1**2+s2**2+s3**2
    born = th**2/(ssq+th**2)
    if born < FLOOR:
        ss = s1+s2+s3
        r = ssq/ss**2
        ac = mpf(26)-r; bc = mpf(2)*r; cc = -r
        t = (-bc+sqrt(bc**2-mpf(4)*ac*cc))/(mpf(2)*ac)
        sc = (ONE-t)/ss
        return [s1*sc, s2*sc, s3*sc, t]
    return list(m)

def step_mp(m, e):
    md, K = ds_mp(m, e)
    return floor_mp(md), K

def make_ev_mp(p):
    pw = (ONE-p)/mpf(2); sc = ONE-FLOOR
    raw = [sqrt(p*sc), sqrt(pw*sc), sqrt(pw*sc), sqrt(FLOOR)]
    tot = sum(raw)
    return [r/tot for r in raw]

# Find equilibrium at 300 digits
print("Finding equilibrium at 300 digits...")
import numpy as np
from scipy.optimize import brentq

# Double precision seed
def ds_np(m, e):
    s, th = m[:3], m[3]
    se, ph = e[:3], e[3]
    s_new = s*se + s*ph + th*se
    th_new = th * ph
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i != j)
    return np.concatenate([s_new, [th_new]]) / (1 - K), K

def floor_np(m):
    s = m[:3]
    born = m[3]**2 / np.sum(m**2)
    if born < 1/27:
        ss = np.sum(s)
        r = np.sum(s**2)/ss**2
        a, b, c = 26-r, 2*r, -r
        t = (-b + np.sqrt(b**2-4*a*c))/(2*a)
        sc = (1-t)/ss
        return np.array([s[0]*sc, s[1]*sc, s[2]*sc, t])
    return m

def step_np(m, e):
    return floor_np(ds_np(m, e)[0])

def fp_np(e):
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(5000):
        m2 = step_np(m, e)
        if np.max(np.abs(m2-m)) < 1e-15: break
        m = m2
    return m2

def make_ev_np(p):
    pw=(1-p)/2; sc=1-1/27
    raw=np.array([np.sqrt(p*sc),np.sqrt(pw*sc),np.sqrt(pw*sc),np.sqrt(1/27)])
    return raw/np.sum(raw)

p_dom_dp = brentq(lambda p: ds_np(fp_np(make_ev_np(p)), make_ev_np(p))[1] - 7/30,
                   0.92, 0.94, xtol=1e-15)

# Refine at 300 digits
e_mp = make_ev_mp(mpf(p_dom_dp))
m_mp = [mpf(x) for x in fp_np(make_ev_np(p_dom_dp))]

for i in range(50000):
    m2, _ = step_mp(m_mp, e_mp)
    diff = max(fabs(m2[j]-m_mp[j]) for j in range(4))
    if diff < mpf(10)**(-280):
        print(f"  Converged at step {i}, residual = {nstr(diff, 5)}")
        break
    m_mp = m2
m_star = m2

_, K_check = ds_mp(m_star, e_mp)
print(f"  |K-7/30| = {nstr(fabs(K_check - K_STAR), 5)}")

# Jacobian at 300 digits
print("Computing Jacobian...")

def jac_4x4(m0, e0, eps_exp=-100):
    eps = mpf(10)**eps_exp
    J = matrix(4, 4)
    f0, _ = step_mp(list(m0), e0)
    for j in range(4):
        mp_ = list(m0); mp_[j] += eps
        fp, _ = step_mp(mp_, e0)
        for i in range(4):
            J[i,j] = (mpf(fp[i]) - mpf(f0[i])) / eps
    return J

J4 = jac_4x4(m_star, e_mp, -100)

V = matrix(4, 3)
for i in range(3): V[i,i] = ONE
for i in range(3): V[3,i] = -ONE
VtV = V.T * V
J3 = VtV**(-1) * V.T * J4 * V

tr_J = sum(J3[i,i] for i in range(3))
det_J = (J3[0,0]*(J3[1,1]*J3[2,2]-J3[1,2]*J3[2,1])
        -J3[0,1]*(J3[1,0]*J3[2,2]-J3[1,2]*J3[2,0])
        +J3[0,2]*(J3[1,0]*J3[2,1]-J3[1,1]*J3[2,0]))
cof = (J3[0,0]*J3[1,1]-J3[0,1]*J3[1,0]
      +J3[0,0]*J3[2,2]-J3[0,2]*J3[2,0]
      +J3[1,1]*J3[2,2]-J3[1,2]*J3[2,1])

lam0 = (tr_J + sqrt(tr_J**2 - 4*cof)) / 2
lam1 = (tr_J - sqrt(tr_J**2 - 4*cof)) / 2

print(f"  tr(J) = {nstr(tr_J, 60)}")
print(f"  cofsum = {nstr(cof, 60)}")
print(f"  det(J) = {nstr(det_J, 5)}")
print(f"  λ₀ = {nstr(lam0, 60)}")
print(f"  λ₁ = {nstr(lam1, 60)}")

# Verify the degree-19 polynomial
coeffs = [
    -4671418103,       # λ^0
     85202393725,      # λ^1
    -210759099887,     # λ^2
    -125920925463,     # λ^3
     65343281380,      # λ^4
    -215345191444,     # λ^5
     536180328934,     # λ^6
     10558534577,      # λ^7
    -458346126206,     # λ^8
     175793550949,     # λ^9
    -210848363196,     # λ^10
    -402921610991,     # λ^11
    -402784632820,     # λ^12
     833013126110,     # λ^13
    -300055480057,     # λ^14
     578938153168,     # λ^15
    -218313582412,     # λ^16
     22094975727,      # λ^17
     1088544475467,    # λ^18
    -497376766077,     # λ^19
]

P_lam0 = sum(mpf(c) * lam0**k for k, c in enumerate(coeffs))
P_lam1 = sum(mpf(c) * lam1**k for k, c in enumerate(coeffs))
print(f"\n  |P(λ₀)| = {nstr(fabs(P_lam0), 5)}")
print(f"  |P(λ₁)| = {nstr(fabs(P_lam1), 5)}")

# ============================================================
# PSLQ on tr(J)
# ============================================================
print("\n" + "=" * 70)
print("PSLQ: Minimal polynomial of tr(J)")
print("=" * 70)

for deg in range(1, 30):
    powers = [tr_J**k for k in range(deg+1)]
    rel = pslq(powers, maxcoeff=10**15)
    if rel is not None:
        check = sum(mpf(c)*tr_J**k for k,c in enumerate(rel))
        if fabs(check) < mpf(10)**(-50):
            print(f"\n  *** tr(J) has degree {deg} over Q ***")
            print(f"  |P_tr(tr)| = {nstr(fabs(check), 5)}")
            max_c = max(abs(c) for c in rel)
            print(f"  Max coefficient: {max_c}")
            if deg <= 25:
                print(f"  Coefficients: {rel}")
            break
    if deg % 5 == 0:
        print(f"  Checked degree {deg}...")
else:
    print("  No polynomial found up to degree 29")

# ============================================================
# PSLQ on cofsum
# ============================================================
print("\n" + "=" * 70)
print("PSLQ: Minimal polynomial of cofsum")
print("=" * 70)

for deg in range(1, 30):
    powers = [cof**k for k in range(deg+1)]
    rel = pslq(powers, maxcoeff=10**15)
    if rel is not None:
        check = sum(mpf(c)*cof**k for k,c in enumerate(rel))
        if fabs(check) < mpf(10)**(-50):
            print(f"\n  *** cofsum has degree {deg} over Q ***")
            print(f"  |P_cof(cof)| = {nstr(fabs(check), 5)}")
            max_c = max(abs(c) for c in rel)
            print(f"  Max coefficient: {max_c}")
            if deg <= 25:
                print(f"  Coefficients: {rel}")
            break
    if deg % 5 == 0:
        print(f"  Checked degree {deg}...")
else:
    print("  No polynomial found up to degree 29")

# ============================================================
# PSLQ on discriminant
# ============================================================
print("\n" + "=" * 70)
print("PSLQ: Minimal polynomial of disc = tr² - 4·cofsum")
print("=" * 70)

disc = tr_J**2 - 4*cof
print(f"  disc = {nstr(disc, 30)}")
print(f"  √disc = {nstr(sqrt(disc), 30)}")

for deg in range(1, 30):
    powers = [disc**k for k in range(deg+1)]
    rel = pslq(powers, maxcoeff=10**15)
    if rel is not None:
        check = sum(mpf(c)*disc**k for k,c in enumerate(rel))
        if fabs(check) < mpf(10)**(-50):
            print(f"\n  *** disc has degree {deg} over Q ***")
            print(f"  |P_disc(disc)| = {nstr(fabs(check), 5)}")
            max_c = max(abs(c) for c in rel)
            print(f"  Max coefficient: {max_c}")
            if deg <= 25:
                print(f"  Coefficients: {rel}")
            break
    if deg % 5 == 0:
        print(f"  Checked degree {deg}...")
else:
    print("  No polynomial found up to degree 29")

# ============================================================
# Check if λ₀ ∈ Q(tr) by testing if tr ∈ Q(λ₀)
# ============================================================
print("\n" + "=" * 70)
print("Is tr(J) ∈ Q(λ₀)?  (PSLQ on [1, λ₀, λ₀², ..., λ₀¹⁸, tr])")
print("=" * 70)

# If tr = Σ aₖ λ₀^k (rational combination), then tr ∈ Q(λ₀) and Q(tr) ⊆ Q(λ₀).
vec = [lam0**k for k in range(19)] + [tr_J]
rel = pslq(vec, maxcoeff=10**15)
if rel is not None:
    check = sum(mpf(c)*v for c,v in zip(rel, vec))
    if fabs(check) < mpf(10)**(-50):
        print(f"  YES: tr(J) ∈ Q(λ₀)")
        print(f"  tr = -(Σ aₖ λ₀^k) / a₁₉  where a₁₉ = {rel[19]}")
        print(f"  |check| = {nstr(fabs(check), 5)}")
    else:
        print(f"  PSLQ found relation but residual too large: {nstr(fabs(check), 5)}")
else:
    print(f"  NO: tr(J) ∉ Q(λ₀) (or coefficients too large)")

# Same for cofsum
print("\nIs cofsum ∈ Q(λ₀)?")
vec2 = [lam0**k for k in range(19)] + [cof]
rel2 = pslq(vec2, maxcoeff=10**15)
if rel2 is not None:
    check2 = sum(mpf(c)*v for c,v in zip(rel2, vec2))
    if fabs(check2) < mpf(10)**(-50):
        print(f"  YES: cofsum ∈ Q(λ₀)")
        print(f"  |check| = {nstr(fabs(check2), 5)}")
    else:
        print(f"  PSLQ found relation but residual too large: {nstr(fabs(check2), 5)}")
else:
    print(f"  NO: cofsum ∉ Q(λ₀) (or coefficients too large)")

print("\nDONE.")
