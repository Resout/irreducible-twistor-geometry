"""
Minimal Polynomial for λ₀ — the leading eigenvalue of the DS transfer operator.

Strategy:
1. Double precision: find exact equilibrium (p_dom, m*, e*) via scipy brentq
2. mpmath: refine the fixed point from the DP solution by iterating at 80 digits
3. mpmath: compute Jacobian at 80 digits
4. PSLQ to identify minimal polynomial
"""

import numpy as np
from scipy.optimize import brentq
import sys

try:
    from mpmath import (mp, mpf, sqrt, log, matrix, nstr, fabs, polyroots,
                        pslq, chop)
except ImportError:
    print("mpmath required: pip install mpmath")
    sys.exit(1)


# ============================================================
# Double-precision implementation (matching spectral_gap_computation.py)
# ============================================================
H = 3
FLOOR_np = 1.0 / 27.0

def ds_np(m, e):
    s, th = m[:3], m[3]
    se, ph = e[:3], e[3]
    s_new = s*se + s*ph + th*se
    th_new = th * ph
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i != j)
    out = np.concatenate([s_new, [th_new]]) / (1 - K)
    return out, K

def floor_np(m):
    s = m[:3]
    born = m[3]**2 / np.sum(m**2)
    if born < FLOOR_np:
        ss = np.sum(s)
        r = np.sum(s**2) / ss**2
        a, b, c = 26.0-r, 2.0*r, -r
        t = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)
        sc = (1.0-t)/ss
        return np.array([s[0]*sc, s[1]*sc, s[2]*sc, t])
    return m

def step_np(m, e):
    return floor_np(ds_np(m, e)[0])

def make_ev_np(p):
    pw = (1-p)/2; sc = 1-FLOOR_np
    raw = np.array([np.sqrt(p*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR_np)])
    return raw/np.sum(raw)

def fp_np(e):
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(5000):
        m2 = step_np(m, e)
        if np.max(np.abs(m2-m)) < 1e-15: break
        m = m2
    return m2

def K_eq_np(p):
    e = make_ev_np(p)
    m = fp_np(e)
    return ds_np(m, e)[1]


print("=" * 70)
print("MINIMAL POLYNOMIAL FOR λ₀")
print("=" * 70)

# Phase 1: exact equilibrium at double precision
print("\nPhase 1: Double-precision equilibrium...")
p_dom = brentq(lambda p: K_eq_np(p) - 7/30, 0.92, 0.94, xtol=1e-15)
e_dp = make_ev_np(p_dom)
m_dp = fp_np(e_dp)
K_dp = ds_np(m_dp, e_dp)[1]

print(f"  p_dom = {p_dom:.15f}")
print(f"  |K-7/30| = {abs(K_dp-7/30):.2e}")
print(f"  m* = {m_dp}")
print(f"  e* = {e_dp}")

# ============================================================
# Phase 2: mpmath refinement
# ============================================================
mp.dps = 80
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

print(f"\nPhase 2: mpmath refinement at {mp.dps} digits...")

# Convert DP solution to mpmath and refine
e_mp = make_ev_mp(mpf(p_dom))
m_mp = [mpf(x) for x in m_dp]

for i in range(20000):
    m2, _ = step_mp(m_mp, e_mp)
    diff = max(fabs(m2[j]-m_mp[j]) for j in range(4))
    if diff < mpf(10)**(-70):
        print(f"  Converged at step {i}, residual = {nstr(diff, 5)}")
        break
    m_mp = m2
m_star = m2

_, K_check = ds_mp(m_star, e_mp)
print(f"  |K-7/30| = {nstr(fabs(K_check - K_STAR), 5)}")
print(f"  Born(θ*) = {nstr(m_star[3]**2/sum(x**2 for x in m_star), 30)}")

# ============================================================
# Phase 3: Jacobian at high precision
# ============================================================
print(f"\nPhase 3: Jacobian computation...")

def jac_4x4_mp(m0, e0, eps_exp=-35):
    """Full 4×4 unconstrained Jacobian, matching the double-precision method."""
    eps = mpf(10)**eps_exp
    J = matrix(4, 4)
    f0, _ = step_mp(list(m0), e0)
    for j in range(4):
        mp_ = list(m0); mp_[j] += eps
        fp, _ = step_mp(mp_, e0)
        for i in range(4):
            J[i,j] = (mpf(fp[i]) - mpf(f0[i])) / eps
    return J

def project_3x3(J4):
    """Project 4×4 Jacobian to L₁=1 tangent space via V = [e₁-e₄, e₂-e₄, e₃-e₄]."""
    # V is 4×3
    V = matrix(4, 3)
    for i in range(3): V[i,i] = ONE
    for i in range(3): V[3,i] = -ONE

    VtV = V.T * V  # 3×3
    VtJV = V.T * J4 * V  # 3×3
    return VtV**(-1) * VtJV

J4_a = jac_4x4_mp(m_star, e_mp, -30)
J4_b = jac_4x4_mp(m_star, e_mp, -35)

J3_a = project_3x3(J4_a)
J3_b = project_3x3(J4_b)

d_ab = max(fabs(J3_a[i,j]-J3_b[i,j]) for i in range(3) for j in range(3))
a_ab = -int(float(mp.log10(d_ab))) if d_ab > 0 else 99
print(f"  Projected Jacobian agreement eps -30 vs -35: {a_ab} digits")

J = J3_b

# Characteristic polynomial
tr_J = J[0,0]+J[1,1]+J[2,2]
m01 = J[0,0]*J[1,1]-J[0,1]*J[1,0]
m02 = J[0,0]*J[2,2]-J[0,2]*J[2,0]
m12 = J[1,1]*J[2,2]-J[1,2]*J[2,1]
sm = m01+m02+m12
det_J = (J[0,0]*(J[1,1]*J[2,2]-J[1,2]*J[2,1])
        -J[0,1]*(J[1,0]*J[2,2]-J[1,2]*J[2,0])
        +J[0,2]*(J[1,0]*J[2,1]-J[1,1]*J[2,0]))

print(f"\n  tr(J) = {nstr(tr_J, 50)}")
print(f"  Σmin  = {nstr(sm, 50)}")
print(f"  det   = {nstr(det_J, 50)}")

# Eigenvalues
roots = polyroots([mpf(1), -tr_J, sm, -det_J])
roots = sorted(roots, key=lambda r: -fabs(r))

print(f"\n  Eigenvalues:")
for i,r in enumerate(roots):
    print(f"    λ_{i} = {nstr(chop(r), 50)}, |λ_{i}| = {nstr(fabs(r), 50)}")

lam0 = roots[0]
print(f"\n  Δ = -ln|λ₀| = {nstr(-log(fabs(lam0)), 50)}")

# ============================================================
# Phase 4: PSLQ
# ============================================================
print(f"\nPhase 4: PSLQ on λ₀...")

for deg in range(2, 40):
    powers = [lam0**k for k in range(deg+1)]
    rel = pslq(powers)
    if rel is not None:
        check = sum(mpf(c)*lam0**k for k,c in enumerate(rel))
        if fabs(check) < mpf(10)**(-20):
            print(f"\n  *** FOUND: degree {deg} ***")
            print(f"  Coefficients: {rel}")
            terms = []
            for k in range(deg,-1,-1):
                if rel[k] != 0:
                    terms.append(f"({rel[k]})λ^{k}" if k>0 else f"({rel[k]})")
            print(f"  {' + '.join(terms)} = 0")
            print(f"  |p(λ₀)| = {nstr(fabs(check), 5)}")
            # Check other eigenvalues
            for i,r in enumerate(roots[1:],1):
                v = sum(mpf(c)*r**k for k,c in enumerate(rel))
                print(f"  |p(λ_{i})| = {nstr(fabs(v), 5)}")
            break
else:
    print("  No integer polynomial found up to degree 39")

# Try with rational coefficients (denominators up to 30^4)
print(f"\nPhase 5: PSLQ on char poly coefficients...")
for name, val in [("tr(J)", tr_J), ("Σmin", sm), ("det(J)", det_J)]:
    # Try p/q for q up to 10^6
    best = (1, 0, fabs(val))
    for q in range(1, 2000001):
        p = int(float(val*q + 0.5))
        err = fabs(val - mpf(p)/mpf(q))
        if err < best[2]:
            best = (q, p, err)
        if err < mpf(10)**(-30):
            break
    print(f"  {name} = {best[1]}/{best[0]}  (err: {nstr(best[2], 5)})")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
