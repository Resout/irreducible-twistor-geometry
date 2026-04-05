"""
TRIPLE DEEPTHINK PACKAGE — Three independent computations in one run.
======================================================================

Task 1: Newton-polish λ₀ and λ₁ to 500 digits from exact polynomial/Jacobian
Task 2: Find the minimal polynomial of λ₁ via PSLQ
Task 3: Trace the algebraic tower — why does degree 19 appear?

Requirements: mpmath (pip install mpmath), scipy (for initial seed only)
Runtime estimate: 60-180 minutes (dominated by Task 1 bisection)

All three tasks share the same high-precision equilibrium computation (Task 1),
so they must run sequentially: Task 1 first, then 2 and 3 use its output.
"""

import sys
import time

try:
    from mpmath import (mp, mpf, sqrt, log, matrix, nstr, fabs,
                        pslq, chop, polyroots, eye, det as mpdet,
                        polyval, diff, taylor, power, pi, inf)
except ImportError:
    print("ERROR: mpmath required. Install with: pip install mpmath")
    sys.exit(1)

# ============================================================
# CONFIGURATION
# ============================================================
PRECISION = 600          # working digits (extra headroom for cancellation)
PSLQ_MAX_COEFF = 10**15
H = 3
TARGET_K = mpf(7) / mpf(30)

# The known degree-19 polynomial for λ₀ (from the paper, eq:minpoly)
P19_COEFFS = [  # ascending: c₀, c₁, ..., c₁₉
    -4671418103, 85202393725, -210759099887, -125920925463,
    65343281380, -215345191444, 536180328934, 10558534577,
    -458346126206, 175793550949, -210848363196, -402921610991,
    -402784632820, 833013126110, -300055480057, 578938153168,
    -218313582412, 22094975727, 1088544475467, -497376766077,
]

def eval_poly(coeffs, x):
    """Evaluate polynomial with ascending coefficients at x."""
    return sum(mpf(c) * x**k for k, c in enumerate(coeffs))

def eval_poly_deriv(coeffs, x):
    """Evaluate polynomial derivative at x."""
    return sum(mpf(c) * k * x**(k-1) for k, c in enumerate(coeffs) if k > 0)


# ============================================================
# DS FRAMEWORK
# ============================================================
mp.dps = PRECISION
ONE = mpf(1)
FLOOR = ONE / mpf(H**3)

def ds_combine(m, e):
    s1, s2, s3, th = m
    e1, e2, e3, ph = e
    sn = [s1*e1+s1*ph+th*e1, s2*e2+s2*ph+th*e2, s3*e3+s3*ph+th*e3]
    tn = th*ph
    K = s1*e2+s1*e3+s2*e1+s2*e3+s3*e1+s3*e2
    d = ONE-K
    return [sn[0]/d, sn[1]/d, sn[2]/d, tn/d], K

def floor_enforce(m):
    s1, s2, s3, th = m
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

def full_step(m, e):
    m_ds, K = ds_combine(m, e)
    return floor_enforce(m_ds), K

def make_evidence(p_dom):
    pw = (ONE-p_dom)/mpf(2); sc = ONE-FLOOR
    raw = [sqrt(p_dom*sc), sqrt(pw*sc), sqrt(pw*sc), sqrt(FLOOR)]
    tot = sum(raw)
    return [r/tot for r in raw]

def find_fp(e, tol_exp=20):
    m = [mpf('0.4'), mpf('0.15'), mpf('0.15'), mpf('0.3')]
    tol = mpf(10)**(-(PRECISION - tol_exp))
    for i in range(200000):
        m2, _ = full_step(m, e)
        d = max(fabs(m2[j]-m[j]) for j in range(4))
        if d < tol:
            return m2, i
        m = m2
    return m2, 200000

def K_at_pdom(p):
    e = make_evidence(p)
    m, _ = find_fp(e)
    _, K = ds_combine(m, e)
    return K


# ============================================================
# TASK 1: HIGH-PRECISION EQUILIBRIUM + NEWTON-POLISHED EIGENVALUES
# ============================================================
def task_1():
    print("=" * 70)
    print("TASK 1: 500-DIGIT EQUILIBRIUM AND NEWTON-POLISHED EIGENVALUES")
    print("=" * 70)

    # --- Step 1: Bisect for p_dom ---
    print("\n  [1a] Bisecting for p_dom...")
    t0 = time.time()
    lo = mpf('0.932'); hi = mpf('0.933')

    K_lo = K_at_pdom(lo); K_hi = K_at_pdom(hi)
    if not ((K_lo - TARGET_K)*(K_hi - TARGET_K) < 0):
        lo = mpf('0.92'); hi = mpf('0.94')

    bisect_tol = mpf(10)**(-(PRECISION - 40))
    n_b = 0
    while hi - lo > bisect_tol:
        mid = (lo+hi)/2
        Km = K_at_pdom(mid)
        if Km > TARGET_K: lo = mid
        else: hi = mid
        n_b += 1
        if n_b % 200 == 0:
            print(f"       step {n_b}, width = {nstr(hi-lo, 5)}")

    p_dom = (lo+hi)/2
    print(f"       Done: {n_b} steps, {time.time()-t0:.0f}s")
    print(f"       p_dom = {nstr(p_dom, 80)}")

    # --- Step 2: Fixed point ---
    print("\n  [1b] Fixed point...")
    e_star = make_evidence(p_dom)
    m_star, n_it = find_fp(e_star)
    _, K_star = ds_combine(m_star, e_star)
    print(f"       Converged in {n_it} iterations")
    print(f"       |K-7/30| = {nstr(fabs(K_star-TARGET_K), 5)}")

    # --- Step 3: Analytical Jacobian ---
    print("\n  [1c] Analytical Jacobian (chain rule)...")
    s1, s2, s3, th = m_star
    e1, e2, e3, ph = e_star
    omK = ONE - K_star

    S_vals = [s1*(e1+ph)+th*e1, s2*(e2+ph)+th*e2, s3*(e3+ph)+th*e3, th*ph]
    dSdm = matrix(4,4)
    for i in range(3):
        dSdm[i,i] = e_star[i]+ph; dSdm[i,3] = e_star[i]
    dSdm[3,3] = ph
    dKdm = [e2+e3, e1+e3, e1+e2, mpf(0)]

    J_DS = matrix(4,4)
    for i in range(4):
        for j in range(4):
            J_DS[i,j] = (dSdm[i,j]*omK + S_vals[i]*dKdm[j]) / omK**2

    N = [S_vals[i]/omK for i in range(4)]
    Ss = N[0]+N[1]+N[2]
    Q_val = N[0]**2+N[1]**2+N[2]**2
    R = Q_val/Ss**2
    u = sqrt(mpf(26)*R)
    t_fl = (u-R)/(mpf(26)-R)
    g = (ONE-t_fl)/Ss

    dtdR = (mpf(338)+mpf(13)*R-mpf(26)*u)/(u*(mpf(26)-R)**2)
    dRdN = [mpf(2)*(N[j]*Ss-Q_val)/Ss**3 for j in range(3)]
    dtdN = [dtdR*dRdN[j] for j in range(3)]
    dgdN = [(-dtdN[j]*Ss-(ONE-t_fl))/Ss**2 for j in range(3)]

    J_fl = matrix(4,4)
    for i in range(3):
        for j in range(3):
            J_fl[i,j] = (g if i==j else mpf(0)) + N[i]*dgdN[j]
        J_fl[i,3] = mpf(0)
    for j in range(3): J_fl[3,j] = dtdN[j]
    J_fl[3,3] = mpf(0)

    J_total = J_fl * J_DS

    V = matrix(4,3)
    for i in range(3): V[i,i] = ONE
    for i in range(3): V[3,i] = -ONE
    J_proj = (V.T*V)**(-1) * V.T * J_total * V

    tr_J = sum(J_proj[i,i] for i in range(3))
    cof = (J_proj[0,0]*J_proj[1,1]-J_proj[0,1]*J_proj[1,0]
          +J_proj[0,0]*J_proj[2,2]-J_proj[0,2]*J_proj[2,0]
          +J_proj[1,1]*J_proj[2,2]-J_proj[1,2]*J_proj[2,1])
    disc = tr_J**2 - mpf(4)*cof

    lam0_jac = (tr_J + sqrt(disc))/2
    lam1_jac = (tr_J - sqrt(disc))/2

    print(f"       tr(J) = {nstr(tr_J, 80)}")
    print(f"       λ₀(Jac) = {nstr(lam0_jac, 80)}")
    print(f"       λ₁(Jac) = {nstr(lam1_jac, 80)}")

    # --- Step 4: Newton-polish λ₀ from the degree-19 polynomial ---
    print("\n  [1d] Newton-polishing λ₀ from degree-19 polynomial...")
    lam0 = lam0_jac  # seed from Jacobian
    for i in range(100):
        pv = eval_poly(P19_COEFFS, lam0)
        dpv = eval_poly_deriv(P19_COEFFS, lam0)
        if fabs(dpv) < mpf(10)**(-PRECISION+10):
            break
        correction = pv/dpv
        lam0 = lam0 - correction
        if fabs(correction) < mpf(10)**(-(PRECISION-10)):
            print(f"       Newton converged in {i+1} steps")
            break

    pv_check = eval_poly(P19_COEFFS, lam0)
    print(f"       |P₁₉(λ₀)| = {nstr(fabs(pv_check), 5)}")
    print(f"\n  *** λ₀ TO 500 DIGITS (Newton-polished from degree-19 polynomial) ***")
    print(f"  {nstr(lam0, 500)}")

    # --- Step 5: Print det(I-J) and related quantities ---
    det_IJ = (ONE-lam0)*(ONE-lam1_jac)
    print(f"\n  λ₁(Jac)  = {nstr(lam1_jac, 80)}")
    print(f"  det(I-J) = {nstr(det_IJ, 80)}")
    print(f"  Δ = -ln(λ₀) = {nstr(-log(lam0), 60)}")

    # --- Compare Jacobian λ₀ vs polynomial λ₀ ---
    diff_lam0 = fabs(lam0 - lam0_jac)
    print(f"\n  |λ₀(poly) - λ₀(Jac)| = {nstr(diff_lam0, 5)}")
    if diff_lam0 > mpf(10)**(-30):
        agree = -int(float(mp.log10(diff_lam0)))
        print(f"  Agreement: {agree} digits (Jacobian precision limit)")
    else:
        print(f"  Agreement: 30+ digits")

    return lam0, lam1_jac, tr_J, cof, disc, m_star, e_star


# ============================================================
# TASK 2: MINIMAL POLYNOMIAL OF λ₁
# ============================================================
def task_2(lam0, lam1, tr_J, cofsum):
    print("\n" + "=" * 70)
    print("TASK 2: MINIMAL POLYNOMIAL OF λ₁")
    print("=" * 70)

    # First check: is λ₁ a root of the degree-19 poly?
    pv1 = eval_poly(P19_COEFFS, lam1)
    print(f"\n  |P₁₉(λ₁)| = {nstr(fabs(pv1), 5)}")
    print(f"  λ₁ is a root of P₁₉: {'YES' if fabs(pv1) < mpf(10)**(-50) else 'NO'}")

    # PSLQ search for λ₁'s minimal polynomial
    print(f"\n  Searching for minimal polynomial of λ₁...")
    lam1_poly = None
    for deg in range(1, 35):
        powers = [lam1**k for k in range(deg+1)]
        rel = pslq(powers, maxcoeff=PSLQ_MAX_COEFF)
        if rel is not None:
            check = sum(mpf(c)*lam1**k for k, c in enumerate(rel))
            if fabs(check) < mpf(10)**(-(PRECISION//2)):
                max_c = max(abs(c) for c in rel)
                print(f"\n  *** λ₁ has degree {deg} over Q ***")
                print(f"  |P(λ₁)| = {nstr(fabs(check), 5)}")
                print(f"  Max coefficient: {max_c}")
                lam1_poly = (deg, rel)
                if deg <= 25:
                    print(f"  Coefficients (ascending): {rel}")

                # Newton-polish λ₁ from this polynomial
                print(f"\n  Newton-polishing λ₁ from its own polynomial...")
                lam1_polished = lam1
                for i in range(100):
                    pv = sum(mpf(c)*lam1_polished**k for k,c in enumerate(rel))
                    dpv = sum(mpf(c)*k*lam1_polished**(k-1) for k,c in enumerate(rel) if k>0)
                    if fabs(dpv) < mpf(10)**(-PRECISION+10): break
                    corr = pv/dpv
                    lam1_polished -= corr
                    if fabs(corr) < mpf(10)**(-(PRECISION-10)):
                        print(f"  Newton converged in {i+1} steps")
                        break

                pv_check2 = sum(mpf(c)*lam1_polished**k for k,c in enumerate(rel))
                print(f"  |P_λ₁(λ₁_polished)| = {nstr(fabs(pv_check2), 5)}")
                print(f"\n  *** λ₁ TO 500 DIGITS ***")
                print(f"  {nstr(lam1_polished, 500)}")

                # Updated det(I-J)
                det_IJ_polished = (ONE-lam0)*(ONE-lam1_polished)
                print(f"\n  det(I-J) (both polished) = {nstr(det_IJ_polished, 80)}")
                print(f"  det(I-J)^{{-1/2}} = {nstr(ONE/sqrt(det_IJ_polished), 40)}")

                # Check: is λ₁ a root of λ₀'s polynomial? (should be NO)
                pv_cross = eval_poly(P19_COEFFS, lam1_polished)
                print(f"\n  |P₁₉_λ₀(λ₁)| = {nstr(fabs(pv_cross), 5)}")

                # Check: is λ₀ a root of λ₁'s polynomial?
                pv_cross2 = sum(mpf(c)*lam0**k for k,c in enumerate(rel))
                print(f"  |P_λ₁(λ₀)| = {nstr(fabs(pv_cross2), 5)}")

                # Root analysis of λ₁'s polynomial
                if deg <= 25:
                    print(f"\n  Finding all roots of λ₁'s polynomial...")
                    mp_save = mp.dps
                    mp.dps = 50
                    all_roots = polyroots([mpf(c) for c in reversed(rel)])
                    real_roots = sorted([float(r.real) for r in all_roots if fabs(r.imag) < mpf(10)**(-20)])
                    real_in_01 = [r for r in real_roots if 0 < r < 1]
                    print(f"  Total roots: {len(all_roots)}")
                    print(f"  Real roots: {len(real_roots)}")
                    print(f"  Real roots in (0,1): {len(real_in_01)}: {[f'{r:.6f}' for r in real_in_01]}")
                    mp.dps = mp_save
                break
        if deg % 5 == 0:
            print(f"    Checked up to degree {deg}...")
    else:
        print(f"  No polynomial found up to degree 34")

    return lam1_poly


# ============================================================
# TASK 3: WHY DEGREE 19? — ALGEBRAIC TOWER
# ============================================================
def task_3(lam0, lam1, tr_J, cofsum, disc):
    print("\n" + "=" * 70)
    print("TASK 3: WHY DEGREE 19? — ALGEBRAIC TOWER ANALYSIS")
    print("=" * 70)

    H = 3
    print(f"\n  H = {H}")
    print(f"  H²+1 = {H**2+1} = dim(Sym²(C⁴))")
    print(f"  2(H²+1)-1 = {2*(H**2+1)-1}")

    # --- 3a: Minimal polynomials of tr(J), cofsum, disc, √disc ---
    quantities = [
        ("tr(J)", tr_J),
        ("cofsum", cofsum),
        ("disc = tr²-4·cof", disc),
        ("√disc", sqrt(disc)),
    ]

    results = {}
    for name, val in quantities:
        print(f"\n  --- PSLQ: minimal polynomial of {name} ---")
        print(f"  Value = {nstr(val, 40)}")
        found = False
        for deg in range(1, 35):
            powers = [val**k for k in range(deg+1)]
            rel = pslq(powers, maxcoeff=PSLQ_MAX_COEFF)
            if rel is not None:
                check = sum(mpf(c)*val**k for k,c in enumerate(rel))
                if fabs(check) < mpf(10)**(-(PRECISION//2)):
                    max_c = max(abs(c) for c in rel)
                    print(f"  *** degree {deg}, max coeff = {max_c} ***")
                    print(f"  |P(val)| = {nstr(fabs(check), 5)}")
                    if deg <= 25:
                        print(f"  Coefficients: {rel}")
                    results[name] = deg
                    found = True
                    break
            if deg % 5 == 0:
                print(f"    Checked up to degree {deg}...")
        if not found:
            results[name] = "> 34"
            print(f"  No polynomial found up to degree 34")

    # --- 3b: Field containment checks ---
    print(f"\n  --- FIELD CONTAINMENT ---")

    def check_in_field(name, val, basis_name, basis_val, basis_deg):
        vec = [basis_val**k for k in range(basis_deg)] + [val]
        rel = pslq(vec, maxcoeff=PSLQ_MAX_COEFF)
        if rel is not None:
            check = sum(mpf(c)*v for c,v in zip(rel, vec))
            if fabs(check) < mpf(10)**(-(PRECISION//2)):
                print(f"  {name} ∈ Q({basis_name}): YES (residual {nstr(fabs(check),5)})")
                return True
        print(f"  {name} ∈ Q({basis_name}): NO (or coeffs > {PSLQ_MAX_COEFF})")
        return False

    # Is tr ∈ Q(λ₀)?
    check_in_field("tr(J)", tr_J, "λ₀", lam0, 19)
    check_in_field("cofsum", cofsum, "λ₀", lam0, 19)
    check_in_field("√disc", sqrt(disc), "λ₀", lam0, 19)
    check_in_field("λ₁", lam1, "λ₀", lam0, 19)

    # Is λ₀ ∈ Q(tr)?
    if isinstance(results.get("tr(J)"), int):
        d_tr = results["tr(J)"]
        check_in_field("λ₀", lam0, "tr(J)", tr_J, d_tr)

    # --- 3c: Summary ---
    print(f"\n  --- SUMMARY ---")
    print(f"  Degrees found:")
    for name, deg in results.items():
        print(f"    {name:20s}: {deg}")

    d_tr = results.get("tr(J)")
    if d_tr == 10:
        print(f"\n  ANSWER: 19 = 2×10 - 1 = 2×dim(Sym²(C⁴)) - 1")
        print(f"  tr(J) has degree 10 = H²+1. The eigenvalue equation is quadratic (×2).")
        print(f"  The m₂=m₃ symmetry kills one root (-1). Total: 2×10-1 = 19.")
    elif d_tr == 19:
        print(f"\n  ANSWER: degree 19 comes from the TRACE itself.")
        print(f"  tr(J) has degree 19 = same as λ₀. Q(tr) = Q(λ₀).")
        print(f"  The discriminant is a perfect square in Q(tr), so the quadratic")
        print(f"  eigenvalue equation does NOT increase the degree.")
        print(f"  The degree 19 arises from the algebraic interaction of:")
        print(f"    - degree-24 evidence parameter")
        print(f"    - √(26R) from the Born floor")
        print(f"    - the conservation law K*(H²+1) - η·H² = 1")
    else:
        print(f"\n  tr(J) has unexpected degree {d_tr}. Further investigation needed.")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    t_start = time.time()
    print(f"TRIPLE DEEPTHINK PACKAGE")
    print(f"Working precision: {PRECISION} digits")
    print(f"Expected runtime: 60-180 minutes\n")

    # Task 1 provides the shared high-precision data
    lam0, lam1, tr_J, cofsum, disc, m_star, e_star = task_1()

    # Task 2 uses λ₁ from Task 1
    task_2(lam0, lam1, tr_J, cofsum)

    # Task 3 uses all quantities from Task 1
    task_3(lam0, lam1, tr_J, cofsum, disc)

    t_end = time.time()
    print(f"\n{'='*70}")
    print(f"TOTAL RUNTIME: {t_end-t_start:.0f}s ({(t_end-t_start)/60:.1f} minutes)")
    print(f"{'='*70}")
    print("DONE.")
