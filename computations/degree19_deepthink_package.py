"""
WHY DEGREE 19? — Self-contained computation package for high-precision analysis.

Context: The leading eigenvalue λ₀ = 0.28291... of the DS transfer operator at the
K*=7/30 equilibrium has minimal polynomial of degree 19 over Q, with Gal = S₁₉.
The evidence parameter v has degree 24 over Q. Why 19?

Conjecture: 19 = 2·dim(Sym²(C⁴)) - 1 = 2·(H²+1) - 1 at H=3.

This script needs to run at 500+ digit precision. It will:

Task A: Compute λ₀, λ₁, tr(J), cofsum to 400+ digits
  - Find p_dom by bisection at 500 digits (not from double-precision seed)
  - Use the ANALYTICAL Jacobian (chain rule J_floor · J_DS, no finite differences)
  - Verify the degree-19 polynomial P(λ₀) = 0 to 10^{-200}+

Task B: Find minimal polynomials of tr(J) and cofsum via PSLQ
  - If tr(J) has degree 10 = H²+1: confirms 19 = 2·10 - 1
  - If tr(J) has degree 19: confirms Q(tr) = Q(λ₀)
  - If tr(J) has degree d ≠ 10, 19: reveals different structure

Task C: Check Q(λ₀) = Q(tr)?
  - PSLQ on [1, λ₀, λ₀², ..., λ₀¹⁸, tr]: is tr a polynomial in λ₀?
  - PSLQ on [1, λ₀, λ₀², ..., λ₀¹⁸, cofsum]: is cofsum a polynomial in λ₀?
  - If both yes: Q(tr, cofsum) = Q(λ₀), degree 19 comes from the trace

Task D: Find minimal polynomial of λ₁
  - If deg(λ₁) = 19: same field, different root
  - If deg(λ₁) ≠ 19: different field entirely

Requirements: mpmath (pip install mpmath), scipy (for initial seed only)
Runtime estimate: 30-120 minutes depending on hardware (the bisection is the bottleneck)
"""

import sys
import time

try:
    from mpmath import (mp, mpf, sqrt, log, matrix, nstr, fabs,
                        pslq, chop, polyroots, eye)
except ImportError:
    print("ERROR: mpmath required. Install with: pip install mpmath")
    sys.exit(1)

# ============================================================
# CONFIGURATION
# ============================================================
PRECISION = 500          # working precision in decimal digits
PSLQ_MAX_COEFF = 10**15  # max integer coefficient for PSLQ
PSLQ_MAX_DEG = 30        # max degree to search
H = 3
TARGET_K = mpf(7) / mpf(30)

# The known degree-19 polynomial (from the paper, eq:minpoly)
P19_COEFFS = [  # ascending order: c₀, c₁, ..., c₁₉
    -4671418103,
     85202393725,
    -210759099887,
    -125920925463,
     65343281380,
    -215345191444,
     536180328934,
     10558534577,
    -458346126206,
     175793550949,
    -210848363196,
    -402921610991,
    -402784632820,
     833013126110,
    -300055480057,
     578938153168,
    -218313582412,
     22094975727,
     1088544475467,
    -497376766077,
]

# ============================================================
# DS FRAMEWORK AT ARBITRARY PRECISION
# ============================================================
mp.dps = PRECISION
ONE = mpf(1)
FLOOR = ONE / mpf(H**3)  # 1/27

def ds_combine(m, e):
    """DS combination. m, e are lists of 4 mpf values. Returns (m_out, K)."""
    s1, s2, s3, th = m
    e1, e2, e3, ph = e
    sn = [s1*e1 + s1*ph + th*e1,
          s2*e2 + s2*ph + th*e2,
          s3*e3 + s3*ph + th*e3]
    tn = th * ph
    K = s1*e2 + s1*e3 + s2*e1 + s2*e3 + s3*e1 + s3*e2
    d = ONE - K
    return [sn[0]/d, sn[1]/d, sn[2]/d, tn/d], K

def floor_enforce(m):
    """Born floor enforcement. Returns new m with Born(theta) >= 1/27."""
    s1, s2, s3, th = m
    ssq = s1**2 + s2**2 + s3**2
    born = th**2 / (ssq + th**2)
    if born < FLOOR:
        ss = s1 + s2 + s3
        r = ssq / ss**2
        a_c = mpf(26) - r
        b_c = mpf(2) * r
        c_c = -r
        disc = b_c**2 - mpf(4)*a_c*c_c
        t = (-b_c + sqrt(disc)) / (mpf(2) * a_c)
        sc = (ONE - t) / ss
        return [s1*sc, s2*sc, s3*sc, t]
    return list(m)

def full_step(m, e):
    """One DS step: combine then enforce floor."""
    m_ds, K = ds_combine(m, e)
    return floor_enforce(m_ds), K

def make_evidence(p_dom):
    """Construct evidence from dominant Born probability p_dom."""
    p_weak = (ONE - p_dom) / mpf(2)
    sc = ONE - FLOOR
    raw = [sqrt(p_dom * sc), sqrt(p_weak * sc),
           sqrt(p_weak * sc), sqrt(FLOOR)]
    tot = sum(raw)
    return [r / tot for r in raw]

def find_fixed_point(e, max_iter=100000):
    """Find DS fixed point by iteration."""
    m = [mpf('0.4'), mpf('0.15'), mpf('0.15'), mpf('0.3')]
    tol = mpf(10)**(-(PRECISION - 20))
    for i in range(max_iter):
        m2, _ = full_step(m, e)
        diff = max(fabs(m2[j] - m[j]) for j in range(4))
        if diff < tol:
            return m2, i
        m = m2
    return m2, max_iter

def K_at_pdom(p_dom):
    """Compute K at the fixed point for given p_dom."""
    e = make_evidence(p_dom)
    m, _ = find_fixed_point(e)
    _, K = ds_combine(m, e)
    return K


# ============================================================
# TASK A: HIGH-PRECISION EQUILIBRIUM AND EIGENVALUES
# ============================================================
def task_A():
    print("=" * 70)
    print("TASK A: HIGH-PRECISION EQUILIBRIUM AND EIGENVALUES")
    print(f"  Working precision: {PRECISION} digits")
    print("=" * 70)

    # Step 1: Find p_dom by bisection at full precision
    print("\n  Step 1: Bisecting for p_dom at full precision...")
    t0 = time.time()
    lo = mpf('0.932')
    hi = mpf('0.933')

    # Verify bracket
    K_lo = K_at_pdom(lo)
    K_hi = K_at_pdom(hi)
    print(f"    K(lo={nstr(lo,6)}) = {nstr(K_lo, 10)}")
    print(f"    K(hi={nstr(hi,6)}) = {nstr(K_hi, 10)}")
    print(f"    Target = {nstr(TARGET_K, 10)}")

    if not ((K_lo - TARGET_K) * (K_hi - TARGET_K) < 0):
        print("    ERROR: bracket doesn't straddle target. Widening...")
        lo = mpf('0.92')
        hi = mpf('0.94')

    bisect_tol = mpf(10)**(-(PRECISION - 30))
    n_bisect = 0
    while hi - lo > bisect_tol:
        mid = (lo + hi) / 2
        K_mid = K_at_pdom(mid)
        if K_mid > TARGET_K:
            lo = mid
        else:
            hi = mid
        n_bisect += 1
        if n_bisect % 100 == 0:
            print(f"    Bisection step {n_bisect}, interval width = {nstr(hi-lo, 5)}")

    p_dom = (lo + hi) / 2
    t1 = time.time()
    print(f"    Done: {n_bisect} bisection steps in {t1-t0:.1f}s")
    print(f"    p_dom = {nstr(p_dom, 60)}")

    # Step 2: Fixed point at full precision
    print("\n  Step 2: Computing fixed point...")
    e_star = make_evidence(p_dom)
    m_star, n_iter = find_fixed_point(e_star)
    _, K_star = ds_combine(m_star, e_star)
    print(f"    Converged in {n_iter} iterations")
    print(f"    |K - 7/30| = {nstr(fabs(K_star - TARGET_K), 5)}")
    print(f"    m* = [{nstr(m_star[0],40)}, {nstr(m_star[1],40)}, {nstr(m_star[2],40)}, {nstr(m_star[3],40)}]")

    # Step 3: ANALYTICAL Jacobian (no finite differences!)
    print("\n  Step 3: Analytical Jacobian (chain rule: J_floor @ J_DS)...")

    # --- J_DS ---
    # DS output (pre-floor): N_i = (s_i(e_i+phi) + theta*e_i) / (1-K)
    s1, s2, s3, th = m_star
    e1, e2, e3, ph = e_star

    oneMinusK = ONE - K_star
    N = [(s1*(e1+ph)+th*e1)/oneMinusK,
         (s2*(e2+ph)+th*e2)/oneMinusK,
         (s3*(e3+ph)+th*e3)/oneMinusK,
         (th*ph)/oneMinusK]

    # dN_i/dm_j before floor (4x4)
    # Numerator S_i and dS_i/dm_j
    S_vals = [s1*(e1+ph)+th*e1, s2*(e2+ph)+th*e2, s3*(e3+ph)+th*e3, th*ph]
    dSdm = matrix(4, 4)
    for i in range(3):
        dSdm[i, i] = e_star[i] + ph  # dS_i/ds_i
        dSdm[i, 3] = e_star[i]        # dS_i/dtheta
    dSdm[3, 3] = ph                    # dS_4/dtheta

    dKdm = [e2+e3, e1+e3, e1+e2, mpf(0)]

    J_DS = matrix(4, 4)
    for i in range(4):
        for j in range(4):
            J_DS[i,j] = (dSdm[i,j] * oneMinusK + S_vals[i] * dKdm[j]) / oneMinusK**2

    # --- J_floor ---
    Ss = N[0] + N[1] + N[2]
    Q_val = N[0]**2 + N[1]**2 + N[2]**2
    R = Q_val / Ss**2
    u = sqrt(mpf(26) * R)
    t = (u - R) / (mpf(26) - R)
    g = (ONE - t) / Ss

    # dtdR = (338 + 13R - 26*sqrt(26R)) / (sqrt(26R) * (26-R)^2)
    dtdR = (mpf(338) + mpf(13)*R - mpf(26)*u) / (u * (mpf(26) - R)**2)

    # dR/dN_j = 2*(N_j*Ss - Q) / Ss^3  (for j=0,1,2)
    dRdN = [mpf(2)*(N[j]*Ss - Q_val) / Ss**3 for j in range(3)]
    dtdN = [dtdR * dRdN[j] for j in range(3)]
    dgdN = [(-dtdN[j]*Ss - (ONE-t)) / Ss**2 for j in range(3)]

    J_floor = matrix(4, 4)
    for i in range(3):
        for j in range(3):
            J_floor[i,j] = (g if i == j else mpf(0)) + N[i] * dgdN[j]
        J_floor[i, 3] = mpf(0)
    J_floor[3, 0] = dtdN[0]
    J_floor[3, 1] = dtdN[1]
    J_floor[3, 2] = dtdN[2]
    J_floor[3, 3] = mpf(0)

    J_total = J_floor * J_DS

    # Project to L1=1 tangent space
    V = matrix(4, 3)
    for i in range(3): V[i,i] = ONE
    for i in range(3): V[3,i] = -ONE
    VtV = V.T * V
    J_proj = VtV**(-1) * V.T * J_total * V

    # Characteristic polynomial coefficients
    tr_J = sum(J_proj[i,i] for i in range(3))
    det_J = (J_proj[0,0]*(J_proj[1,1]*J_proj[2,2]-J_proj[1,2]*J_proj[2,1])
            -J_proj[0,1]*(J_proj[1,0]*J_proj[2,2]-J_proj[1,2]*J_proj[2,0])
            +J_proj[0,2]*(J_proj[1,0]*J_proj[2,1]-J_proj[1,1]*J_proj[2,0]))
    cof01 = J_proj[0,0]*J_proj[1,1] - J_proj[0,1]*J_proj[1,0]
    cof02 = J_proj[0,0]*J_proj[2,2] - J_proj[0,2]*J_proj[2,0]
    cof12 = J_proj[1,1]*J_proj[2,2] - J_proj[1,2]*J_proj[2,1]
    cofsum = cof01 + cof02 + cof12

    # Eigenvalues
    disc = tr_J**2 - mpf(4)*cofsum
    sqrt_disc = sqrt(disc)
    lam0 = (tr_J + sqrt_disc) / 2
    lam1 = (tr_J - sqrt_disc) / 2

    print(f"    tr(J)   = {nstr(tr_J, 80)}")
    print(f"    cofsum  = {nstr(cofsum, 80)}")
    print(f"    det(J)  = {nstr(det_J, 5)}")
    print(f"    disc    = {nstr(disc, 40)}")
    print(f"    λ₀      = {nstr(lam0, 80)}")
    print(f"    λ₁      = {nstr(lam1, 80)}")
    print(f"    Δ       = {nstr(-log(lam0), 40)}")

    # Verify degree-19 polynomial
    P_val = sum(mpf(c) * lam0**k for k, c in enumerate(P19_COEFFS))
    print(f"\n    |P₁₉(λ₀)| = {nstr(fabs(P_val), 5)}")
    P_val1 = sum(mpf(c) * lam1**k for k, c in enumerate(P19_COEFFS))
    print(f"    |P₁₉(λ₁)| = {nstr(fabs(P_val1), 5)}")

    return lam0, lam1, tr_J, cofsum, disc


# ============================================================
# TASK B: MINIMAL POLYNOMIALS VIA PSLQ
# ============================================================
def task_B(lam0, lam1, tr_J, cofsum, disc):
    print("\n" + "=" * 70)
    print("TASK B: MINIMAL POLYNOMIALS VIA PSLQ")
    print(f"  maxcoeff = {PSLQ_MAX_COEFF}, max degree = {PSLQ_MAX_DEG}")
    print("=" * 70)

    def find_minpoly(name, val, max_deg=PSLQ_MAX_DEG):
        print(f"\n  --- {name} ---")
        print(f"  Value = {nstr(val, 40)}")
        for deg in range(1, max_deg + 1):
            powers = [val**k for k in range(deg + 1)]
            rel = pslq(powers, maxcoeff=PSLQ_MAX_COEFF)
            if rel is not None:
                check = sum(mpf(c) * val**k for k, c in enumerate(rel))
                if fabs(check) < mpf(10)**(-PRECISION//2):
                    max_c = max(abs(c) for c in rel)
                    print(f"  *** FOUND: degree {deg}, max coeff = {max_c} ***")
                    print(f"  |P({name})| = {nstr(fabs(check), 5)}")
                    if deg <= 25:
                        print(f"  Coefficients: {rel}")
                    return deg, rel
            if deg % 5 == 0:
                print(f"    Checked up to degree {deg}...")
        print(f"  No polynomial found up to degree {max_deg}")
        return None, None

    results = {}

    # λ₀ (should find degree 19)
    d, c = find_minpoly("λ₀", lam0)
    results['lam0'] = (d, c)

    # λ₁
    d, c = find_minpoly("λ₁", lam1)
    results['lam1'] = (d, c)

    # tr(J)
    d, c = find_minpoly("tr(J)", tr_J)
    results['tr'] = (d, c)

    # cofsum
    d, c = find_minpoly("cofsum", cofsum)
    results['cof'] = (d, c)

    # discriminant
    d, c = find_minpoly("disc", disc)
    results['disc'] = (d, c)

    # √disc
    d, c = find_minpoly("√disc", sqrt(disc))
    results['sqrt_disc'] = (d, c)

    return results


# ============================================================
# TASK C: FIELD CONTAINMENT
# ============================================================
def task_C(lam0, lam1, tr_J, cofsum):
    print("\n" + "=" * 70)
    print("TASK C: FIELD CONTAINMENT — IS Q(tr) = Q(λ₀)?")
    print("=" * 70)

    def check_containment(name, val, basis_name, basis_val, basis_deg):
        """Check if val ∈ Q(basis_val) by PSLQ on [1, basis, basis², ..., basis^{d-1}, val]."""
        print(f"\n  Is {name} ∈ Q({basis_name})? (PSLQ with {basis_deg} powers)")
        vec = [basis_val**k for k in range(basis_deg)] + [val]
        rel = pslq(vec, maxcoeff=PSLQ_MAX_COEFF)
        if rel is not None:
            check = sum(mpf(c)*v for c, v in zip(rel, vec))
            if fabs(check) < mpf(10)**(-PRECISION//2):
                print(f"    YES: {name} = -(sum a_k {basis_name}^k) / a_last")
                print(f"    |residual| = {nstr(fabs(check), 5)}")
                print(f"    Denominator coefficient: {rel[-1]}")
                return True, rel
            else:
                print(f"    PSLQ found relation but residual = {nstr(fabs(check), 5)}")
                return False, None
        else:
            print(f"    NO (or coefficients > {PSLQ_MAX_COEFF})")
            return False, None

    # Assuming λ₀ has degree 19: check if tr, cofsum ∈ Q(λ₀)
    check_containment("tr(J)", tr_J, "λ₀", lam0, 19)
    check_containment("cofsum", cofsum, "λ₀", lam0, 19)
    check_containment("λ₁", lam1, "λ₀", lam0, 19)

    # Also check: is λ₀ ∈ Q(tr)?  (need to know degree of tr first)
    # This is attempted if Task B found the degree of tr


# ============================================================
# TASK D: SUMMARY AND DEGREE ANALYSIS
# ============================================================
def task_D(results):
    print("\n" + "=" * 70)
    print("TASK D: SUMMARY")
    print("=" * 70)

    print(f"\n  H = {H}")
    print(f"  H²+1 = {H**2+1} = dim(Sym²(C⁴))")
    print(f"  2(H²+1) - 1 = {2*(H**2+1)-1}")

    print(f"\n  Degrees found:")
    for name, (deg, _) in results.items():
        print(f"    {name:12s}: degree {deg if deg else '> ' + str(PSLQ_MAX_DEG)}")

    deg_lam0 = results.get('lam0', (None,))[0]
    deg_lam1 = results.get('lam1', (None,))[0]
    deg_tr = results.get('tr', (None,))[0]
    deg_cof = results.get('cof', (None,))[0]

    if deg_lam0 == 19:
        print(f"\n  ✓ λ₀ has degree 19 = 2(H²+1)-1. CONFIRMED.")
    if deg_lam1 is not None:
        if deg_lam1 == 19:
            print(f"  λ₁ also has degree 19 — same field, different root")
        else:
            print(f"  λ₁ has degree {deg_lam1} ≠ 19 — DIFFERENT field")
    if deg_tr is not None:
        if deg_tr == 10:
            print(f"  tr(J) has degree 10 = H²+1 = dim(Sym²(C⁴))")
            print(f"  → 19 = 2×10 - 1 confirmed: eigenvalue is quadratic over Sym²")
        elif deg_tr == 19:
            print(f"  tr(J) has degree 19 = same as λ₀ → Q(tr) = Q(λ₀)")
            print(f"  → degree comes from the TRACE, not the eigenvalue extraction")
        else:
            print(f"  tr(J) has degree {deg_tr} — unexpected!")
    if deg_cof is not None:
        print(f"  cofsum has degree {deg_cof}")

    print(f"\n  Key question answered:")
    if deg_tr == 10:
        print(f"  WHY 19? Because tr(J) has degree 10 = dim(Sym²(C⁴)), the eigenvalue")
        print(f"  equation is quadratic (×2 = 20), and m₂=m₃ symmetry kills one root (-1).")
    elif deg_tr == 19:
        print(f"  WHY 19? Because tr(J) itself has degree 19. The quadratic eigenvalue")
        print(f"  equation doesn't increase the degree (disc is a perfect square in Q(tr)).")
        print(f"  The degree comes from the algebraic interaction of the degree-24 evidence")
        print(f"  parameter with the floor's √(26R) operation.")
    else:
        print(f"  WHY 19? Still unclear — tr(J) degree not determined.")
        print(f"  Need higher precision or coefficient bound.")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print(f"Working precision: {PRECISION} digits")
    print(f"Expected runtime: 30-120 minutes\n")

    t_start = time.time()

    lam0, lam1, tr_J, cofsum, disc = task_A()
    results = task_B(lam0, lam1, tr_J, cofsum, disc)
    task_C(lam0, lam1, tr_J, cofsum)
    task_D(results)

    t_end = time.time()
    print(f"\nTotal runtime: {t_end - t_start:.1f}s ({(t_end-t_start)/60:.1f} minutes)")
    print("\nDONE.")
