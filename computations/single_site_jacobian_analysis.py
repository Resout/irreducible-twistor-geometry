#!/usr/bin/env python3
"""
Analyse the single-site 4×4 Jacobian to understand the 1×1 (S₋) eigenvalue.

Key question: Is the S₋ eigenvalue exactly 1/H = 1/3?
"""

from mpmath import mp, mpf, matrix, sqrt, log, fabs, nstr, eig, chop
from mpmath import findroot, diff

mp.dps = 80
H = 3
FLOOR = mpf(1) / mpf(H**3)


# DS dynamics (same as before)
def ds_combine(m, e):
    s = m[:3]; theta = m[3]; ev = e[:3]; phi = e[3]
    s_pre = [s[i]*ev[i] + s[i]*phi + theta*ev[i] for i in range(3)]
    theta_pre = theta * phi
    total_pre = sum(s_pre) + theta_pre
    K = mpf(1) - total_pre
    if fabs(mpf(1) - K) < mpf(10)**(-70): return list(m), K
    denom = mpf(1) - K
    return [sp / denom for sp in s_pre] + [theta_pre / denom], K

def born_prob(m):
    L2sq = sum(x**2 for x in m)
    if L2sq < mpf(10)**(-60): return mpf(0)
    return m[3]**2 / L2sq

def enforce_floor(m, floor_val=None):
    if floor_val is None: floor_val = FLOOR
    b = born_prob(m)
    if b >= floor_val - mpf(10)**(-60): return list(m)
    S = sum(m[:3])
    if S < mpf(10)**(-60): return list(m)
    Sq = sum(x**2 for x in m[:3])
    A_c = mpf(26) * S**2 - Sq
    B_c = mpf(2) * Sq; C_c = -Sq
    disc = B_c**2 - 4*A_c*C_c
    if disc < 0: return list(m)
    t1 = (-B_c + sqrt(disc)) / (2*A_c)
    t2 = (-B_c - sqrt(disc)) / (2*A_c)
    cands = [t for t in [t1, t2] if mpf(0) < t < mpf(1)]
    if not cands: return list(m)
    t = min(cands, key=lambda x: fabs(x - m[3]))
    alpha = (mpf(1) - t) / S
    return [m[i]*alpha for i in range(3)] + [t]

def ds_step(m, e):
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds), K


def find_equilibrium():
    from scipy.optimize import fsolve
    def eq_rough(p):
        s1, theta, w1, phi = p
        s2 = (1.0 - s1 - theta) / 2.0; w2 = (1.0 - w1 - phi) / 2.0
        m = [s1, s2, s2, theta]; e = [w1, w2, w2, phi]
        L2m = s1**2 + 2*s2**2 + theta**2
        L2e = w1**2 + 2*w2**2 + phi**2
        K = s1*w2 + s1*w2 + s2*w1 + s2*w2 + s2*w1 + s2*w2
        m_mp = [mpf(x) for x in m]; e_mp = [mpf(x) for x in e]
        m_out, _ = ds_step(m_mp, e_mp)
        return [theta**2/L2m - 1/27, phi**2/L2e - 1/27, K - 7/30, float(m_out[0]) - s1]
    sol = fsolve(eq_rough, [0.787, 0.155, 0.631, 0.129])

    def eq_mp(*p):
        s1, theta, w1, phi = p
        s2 = (mpf(1)-s1-theta)/2; w2 = (mpf(1)-w1-phi)/2
        m = [s1,s2,s2,theta]; e = [w1,w2,w2,phi]
        eq1 = theta**2/(s1**2+2*s2**2+theta**2) - FLOOR
        eq2 = phi**2/(w1**2+2*w2**2+phi**2) - FLOOR
        K = s1*w2+s1*w2+s2*w1+s2*w2+s2*w1+s2*w2
        m_out, _ = ds_step(m, e)
        return [eq1, eq2, K - mpf(7)/30, m_out[0] - s1]

    x = findroot(eq_mp, [mpf(str(v)) for v in sol])
    s1,theta,w1,phi = x
    s2 = (mpf(1)-s1-theta)/2; w2 = (mpf(1)-w1-phi)/2
    return [s1,s2,s2,theta], [w1,w2,w2,phi]


def main():
    print("=" * 70)
    print("SINGLE-SITE JACOBIAN ANALYSIS")
    print("=" * 70)

    m_star, e_star = find_equilibrium()
    print(f"m* = [{', '.join(nstr(x,20) for x in m_star)}]")
    print(f"e* = [{', '.join(nstr(x,20) for x in e_star)}]")

    # ============================================================
    # 4×4 single-site Jacobian: ∂Φ/∂m at (m*, e*)
    # ============================================================
    print("\n--- Single-site 4×4 Jacobian ---")
    eps = mpf(10)**(-35)

    J4 = matrix(4, 4)
    for j in range(4):
        mp_p = list(m_star); mp_m = list(m_star)
        mp_p[j] += eps; mp_m[j] -= eps
        fp, _ = ds_step(mp_p, e_star)
        fm, _ = ds_step(mp_m, e_star)
        for i in range(4):
            J4[i,j] = (fp[i] - fm[i]) / (2*eps)

    print("\nJ4 =")
    for i in range(4):
        row = [nstr(J4[i,j], 12) for j in range(4)]
        print(f"  [{', '.join(row)}]")

    evals4, evecs4 = eig(J4)
    evals4 = [chop(e) for e in evals4]
    evals4_sorted = sorted(evals4, key=lambda x: -float(fabs(x)))

    print("\nEigenvalues of single-site Jacobian:")
    for k, ev in enumerate(evals4_sorted):
        print(f"  λ_{k} = {nstr(ev, 40)}  |λ| = {nstr(fabs(ev), 40)}")

    # ============================================================
    # s₂↔s₃ decomposition of J4
    # ============================================================
    print("\n--- s₂↔s₃ decomposition ---")

    # Change of basis: (s₁, s₂, s₃, θ) → (s₁, S₊, θ, S₋)
    T = matrix(4, 4)
    T[0,0] = 1
    T[1,1] = mpf(1)/sqrt(2); T[1,2] = mpf(1)/sqrt(2)
    T[2,3] = 1
    T[3,1] = mpf(1)/sqrt(2); T[3,2] = -mpf(1)/sqrt(2)

    Tinv = T**(-1)
    J4_rot = Tinv * J4 * T

    print("\nJ4 in (s₁, S₊, θ, S₋) basis:")
    for i in range(4):
        row = [nstr(J4_rot[i,j], 12) for j in range(4)]
        print(f"  [{', '.join(row)}]")

    # Check 3×3 ↔ 1×1 coupling
    coupling = max(fabs(J4_rot[i,3]) for i in range(3))
    coupling2 = max(fabs(J4_rot[3,j]) for j in range(3))
    print(f"\n||3→1 coupling|| = {nstr(coupling, 10)}")
    print(f"||1→3 coupling|| = {nstr(coupling2, 10)}")

    # 1×1 eigenvalue
    s_minus_eval = J4_rot[3,3]
    print(f"\nS₋ eigenvalue = {nstr(s_minus_eval, 40)}")
    print(f"1/H = {nstr(mpf(1)/H, 40)}")
    print(f"difference = {nstr(s_minus_eval - mpf(1)/H, 15)}")
    print(f"ratio = {nstr(s_minus_eval * H, 30)}")

    # 3×3 block
    J33 = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            J33[i,j] = J4_rot[i,j]

    evals33, _ = eig(J33)
    evals33 = [chop(e) for e in evals33]
    evals33_sorted = sorted(evals33, key=lambda x: -float(fabs(x)))
    print("\n3×3 block eigenvalues:")
    for k, ev in enumerate(evals33_sorted):
        print(f"  {nstr(ev, 40)}")

    # ============================================================
    # TEST: what ARE these eigenvalues in terms of framework numbers?
    # ============================================================
    print("\n--- Framework number identification ---")

    K_star = mpf(7) / 30
    targets = {
        "1/H": mpf(1)/H,
        "1/(H+1)": mpf(1)/(H+1),
        "K*": K_star,
        "1-K*": 1-K_star,
        "(H-1)/H²": mpf(H-1)/H**2,
        "K*/(H-1)": K_star/(H-1),
        "1/(H²-1)": mpf(1)/(H**2-1),
        "(H-1)/(H+1)": mpf(H-1)/(H+1),
        "1/H²": mpf(1)/H**2,
        "2/(H(H+1))": mpf(2)/(H*(H+1)),
        "K*/2": K_star/2,
        "1-2K*": 1-2*K_star,
        "H/(H²+1)": mpf(H)/(H**2+1),
        "1/(2H)": mpf(1)/(2*H),
        "(H-1)²/H²": mpf(H-1)**2/H**2,
        "K*²": K_star**2,
    }

    all_evals = list(evals4_sorted)
    for k, ev in enumerate(all_evals):
        if float(fabs(ev)) < 1e-10:
            continue
        print(f"\n  λ_{k} = {nstr(ev, 20)}")
        # Find closest framework number
        best_name = "???"
        best_err = mpf(1)
        for name, val in targets.items():
            err = fabs(fabs(ev) - val)
            if err < best_err:
                best_err = err
                best_name = name
        pct = float(best_err / fabs(ev) * 100)
        print(f"  closest: {best_name} = {nstr(targets[best_name], 15)} (err {pct:.4f}%)")

    # ============================================================
    # ANALYTIC S₋ EIGENVALUE
    # ============================================================
    print("\n\n--- Analytic S₋ computation ---")
    print("""
    At equilibrium, s₂* = s₃*, e₂* = e₃*.
    Under the perturbation δs₂ = +ε, δs₃ = -ε (pure S₋ mode):

    DS combination:
      δs₂' ∝ (e₂* + φ*)·δs₂ = (e₂* + φ*)·ε
      δs₃' ∝ (e₃* + φ*)·δs₃ = (e₂* + φ*)·(-ε)

    So S₋' = (e₂* + φ*)·S₋ / (1-K) in the DS step.

    Floor enforcement doesn't affect S₋ because it depends on
    S = s₁+s₂+s₃ and Sq = s₁²+s₂²+s₃² which are both
    symmetric in s₂, s₃. The S₋ perturbation doesn't change
    S or Sq to first order (δS = δs₂ + δs₃ = 0,
    δSq = 2s₂*δs₂ + 2s₃*δs₃ = 2s₂*(ε - ε) = 0).

    Therefore the S₋ eigenvalue is EXACTLY:
      λ_{S₋} = (e₂* + φ*) / (1 - K*)
    """)

    e2 = e_star[1]
    phi = e_star[3]
    lambda_Sminus_analytic = (e2 + phi) / (1 - K_star)

    print(f"  e₂* = {nstr(e2, 40)}")
    print(f"  φ*  = {nstr(phi, 40)}")
    print(f"  1-K* = {nstr(1 - K_star, 40)}")
    print(f"  (e₂* + φ*) / (1 - K*) = {nstr(lambda_Sminus_analytic, 40)}")
    print(f"  J4[3,3] (numerical)    = {nstr(s_minus_eval, 40)}")
    print(f"  difference = {nstr(fabs(lambda_Sminus_analytic - s_minus_eval), 15)}")

    # What IS e₂* + φ* in terms of framework numbers?
    sum_e2_phi = e2 + phi
    print(f"\n  e₂* + φ* = {nstr(sum_e2_phi, 40)}")
    print(f"  = {nstr(sum_e2_phi, 40)}")

    # Check: is (e₂+φ)/(1-K*) = 1/3?
    target_third = mpf(1) / 3
    ratio_to_third = lambda_Sminus_analytic / target_third
    print(f"\n  λ_{'{S₋}'} / (1/3) = {nstr(ratio_to_third, 30)}")
    print(f"  If exactly 1/3, need e₂* + φ* = (1-K*)/3 = 23/90")
    print(f"  23/90 = {nstr(mpf(23)/90, 40)}")
    print(f"  e₂* + φ* = {nstr(sum_e2_phi, 40)}")
    print(f"  difference = {nstr(sum_e2_phi - mpf(23)/90, 15)}")

    # ============================================================
    # ALSO: what does the single-site Jacobian's det and trace tell us?
    # ============================================================
    print("\n\n--- Trace and determinant of J4 ---")
    tr4 = sum(J4[i,i] for i in range(4))
    # det via eigenvalues
    det4 = evals4_sorted[0] * evals4_sorted[1]
    for i in range(2, len(evals4_sorted)):
        det4 *= evals4_sorted[i]

    print(f"  tr(J4) = {nstr(tr4, 30)}")
    print(f"  det(J4) = {nstr(det4, 15)}")
    print(f"  tr/H = {nstr(tr4/H, 15)}")
    print(f"  tr × H = {nstr(tr4*H, 15)}")

    # Sum of significant eigenvalues
    sig_sum = sum(fabs(e) for e in evals4_sorted if float(fabs(e)) > 1e-10)
    print(f"  Σ|λ_significant| = {nstr(sig_sum, 30)}")
    print(f"  vs 1 = {nstr(sig_sum - 1, 15)}")
    print(f"  vs K* = {nstr(sig_sum - K_star, 15)}")


if __name__ == "__main__":
    main()
