#!/usr/bin/env python3
"""
Trace the Born floor's effect on the S₋ eigenvalue.

The equilibrium sits ON the floor boundary. The DS combination
takes it below the floor, then floor enforcement brings it back.
The composite derivative in the S₋ direction = α × (e₂+φ)/(1-K*)
where α = (1-θ*)/S_DS_out is the floor projection factor.
"""

from mpmath import mp, mpf, matrix, sqrt, log, fabs, nstr, eig, chop, findroot

mp.dps = 80
H = 3
FLOOR = mpf(1) / mpf(H**3)


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
    A_c = mpf(26) * S**2 - Sq; B_c = mpf(2) * Sq; C_c = -Sq
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
    def eq_r(p):
        s1, theta, w1, phi = p
        s2 = (1-s1-theta)/2; w2 = (1-w1-phi)/2
        m = [s1,s2,s2,theta]; e = [w1,w2,w2,phi]
        L2m = s1**2+2*s2**2+theta**2; L2e = w1**2+2*w2**2+phi**2
        K = s1*w2+s1*w2+s2*w1+s2*w2+s2*w1+s2*w2
        m_mp = [mpf(x) for x in m]; e_mp = [mpf(x) for x in e]
        m_out, _ = ds_step(m_mp, e_mp)
        return [theta**2/L2m-1/27, phi**2/L2e-1/27, K-7/30, float(m_out[0])-s1]
    sol = fsolve(eq_r, [0.787, 0.155, 0.631, 0.129])
    def eq_mp(*p):
        s1,theta,w1,phi = p
        s2 = (mpf(1)-s1-theta)/2; w2 = (mpf(1)-w1-phi)/2
        m = [s1,s2,s2,theta]; e = [w1,w2,w2,phi]
        eq1 = theta**2/(s1**2+2*s2**2+theta**2) - FLOOR
        eq2 = phi**2/(w1**2+2*w2**2+phi**2) - FLOOR
        K = s1*w2+s1*w2+s2*w1+s2*w2+s2*w1+s2*w2
        m_out, _ = ds_step(m, e)
        return [eq1, eq2, K-mpf(7)/30, m_out[0]-s1]
    x = findroot(eq_mp, [mpf(str(v)) for v in sol])
    s1,theta,w1,phi = x
    s2 = (mpf(1)-s1-theta)/2; w2 = (mpf(1)-w1-phi)/2
    return [s1,s2,s2,theta], [w1,w2,w2,phi]


def main():
    print("=" * 70)
    print("FLOOR DERIVATIVE IN THE S₋ DIRECTION")
    print("=" * 70)

    m_star, e_star = find_equilibrium()
    s1, s2, s3, theta = m_star
    e1, e2, e3, phi = e_star
    S_star = s1 + s2 + s3  # = 1 - theta
    K_star = mpf(7)/30

    print(f"\nm* = ({nstr(s1,20)}, {nstr(s2,20)}, {nstr(s3,20)}, {nstr(theta,20)})")
    print(f"e* = ({nstr(e1,20)}, {nstr(e2,20)}, {nstr(e3,20)}, {nstr(phi,20)})")
    print(f"S* = 1-θ* = {nstr(S_star, 20)}")

    # Step 1: DS combination at equilibrium
    m_ds, K_val = ds_combine(m_star, e_star)
    print(f"\nDS combination output (before floor):")
    print(f"  m_DS = ({', '.join(nstr(x,20) for x in m_ds)})")
    print(f"  K = {nstr(K_val, 20)}")
    print(f"  Born(m_DS) = {nstr(born_prob(m_ds), 20)}")
    print(f"  Floor = {nstr(FLOOR, 20)}")
    print(f"  Floor active? {float(born_prob(m_ds)) < float(FLOOR)}")

    S_ds = sum(m_ds[:3])
    theta_ds = m_ds[3]
    print(f"\n  S_DS = {nstr(S_ds, 30)}")
    print(f"  θ_DS = {nstr(theta_ds, 30)}")

    # Step 2: Floor enforcement factor
    # After floor: new_si = (1-t)/S_DS * m_ds_si, new_theta = t
    # where t solves 26·t²·S_DS² = (1-t)²·Sq_DS
    m_final = enforce_floor(m_ds)
    print(f"\nFloor output:")
    print(f"  m_final = ({', '.join(nstr(x,20) for x in m_final)})")

    alpha_floor = m_final[0] / m_ds[0]  # rescaling factor
    t_floor = m_final[3]
    print(f"\n  Floor α = (1-t)/S_DS = {nstr(alpha_floor, 30)}")
    print(f"  Floor t = θ_final = {nstr(t_floor, 30)}")
    print(f"  (1-t)/S_DS = {nstr((1-t_floor)/S_ds, 30)}")

    # Step 3: S₋ eigenvalue decomposition
    print("\n" + "=" * 70)
    print("S₋ EIGENVALUE DECOMPOSITION")
    print("=" * 70)

    # The S₋ eigenvalue of the full map is:
    # λ_{S₋} = α_floor × λ_{S₋,DS}
    # where λ_{S₋,DS} = (e₂+φ)/(1-K*) is the DS part

    lambda_DS_Sminus = (e2 + phi) / (1 - K_star)
    print(f"\n  DS part: (e₂+φ)/(1-K*) = {nstr(lambda_DS_Sminus, 30)}")
    print(f"  Floor factor: α = {nstr(alpha_floor, 30)}")
    print(f"  Product: α × DS = {nstr(alpha_floor * lambda_DS_Sminus, 30)}")

    # Numerical S₋ eigenvalue from Jacobian
    eps = mpf(10)**(-35)
    # Perturb in S₋ direction: δs₂=+ε, δs₃=-ε
    mp_p = list(m_star); mp_p[1] += eps; mp_p[2] -= eps
    mp_m = list(m_star); mp_m[1] -= eps; mp_m[2] += eps
    fp, _ = ds_step(mp_p, e_star)
    fm, _ = ds_step(mp_m, e_star)

    # S₋ component of output = (out₂ - out₃)/2
    Sminus_out_p = (fp[1] - fp[2]) / 2
    Sminus_out_m = (fm[1] - fm[2]) / 2
    # Input S₋ = ε (since δs₂=+ε, δs₃=-ε → S₋ = ε)
    lambda_Sminus_numerical = (Sminus_out_p - Sminus_out_m) / (2*eps)

    print(f"\n  Direct numerical S₋ eigenvalue: {nstr(lambda_Sminus_numerical, 30)}")

    # Compare with J4[1,1] - J4[1,2]
    mp_p2 = list(m_star); mp_p2[1] += eps
    mp_m2 = list(m_star); mp_m2[1] -= eps
    fp2, _ = ds_step(mp_p2, e_star)
    fm2, _ = ds_step(mp_m2, e_star)
    J11 = (fp2[1] - fm2[1]) / (2*eps)
    J12 = (fp2[2] - fm2[2]) / (2*eps)  # wait, this is ∂out₂/∂s₂ and ∂out₃/∂s₂

    # Actually: J4[1,1] = ∂(out_s₂)/∂s₂ and J4[1,2] = ∂(out_s₂)/∂s₃
    mp_p3 = list(m_star); mp_p3[2] += eps
    mp_m3 = list(m_star); mp_m3[2] -= eps
    fp3, _ = ds_step(mp_p3, e_star)
    fm3, _ = ds_step(mp_m3, e_star)

    J_s2_s2 = (fp2[1] - fm2[1]) / (2*eps)  # ∂out₂/∂s₂
    J_s2_s3 = (fp3[1] - fm3[1]) / (2*eps)  # ∂out₂/∂s₃
    J_s3_s2 = (fp2[2] - fm2[2]) / (2*eps)  # ∂out₃/∂s₂
    J_s3_s3 = (fp3[2] - fm3[2]) / (2*eps)  # ∂out₃/∂s₃

    print(f"\n  J₂₂ = ∂out₂/∂s₂ = {nstr(J_s2_s2, 30)}")
    print(f"  J₂₃ = ∂out₂/∂s₃ = {nstr(J_s2_s3, 30)}")
    print(f"  J₃₂ = ∂out₃/∂s₂ = {nstr(J_s3_s2, 30)}")
    print(f"  J₃₃ = ∂out₃/∂s₃ = {nstr(J_s3_s3, 30)}")
    print(f"  J₂₂ - J₂₃ = {nstr(J_s2_s2 - J_s2_s3, 30)} (should = S₋ eigenvalue)")

    # ============================================================
    # WHAT DOES THE FLOOR FACTOR α EQUAL?
    # ============================================================
    print("\n" + "=" * 70)
    print("WHAT IS α IN FRAMEWORK TERMS?")
    print("=" * 70)

    print(f"\n  α = {nstr(alpha_floor, 40)}")
    print(f"  1/α = {nstr(1/alpha_floor, 40)}")

    # Test against framework numbers
    K = K_star
    tests = {
        "1": mpf(1),
        "(1-K)": 1-K,
        "(1-2K)/(1-K)": (1-2*K)/(1-K),
        "S*/S_ds": S_star/S_ds,
        "(H-1)/H": mpf(H-1)/H,
        "1-1/H": 1-mpf(1)/H,
        "(H²-1)/H²": mpf(H**2-1)/H**2,
        "θ*/θ_ds": theta/theta_ds,
        "1-K²": 1-K**2,
    }

    for name, val in tests.items():
        err = fabs(alpha_floor - val) / alpha_floor * 100
        if float(err) < 5:
            print(f"  α ≈ {name} = {nstr(val, 20)} (err {nstr(err, 5)}%)")

    # ============================================================
    # THE FULL S₋ EIGENVALUE
    # ============================================================
    print("\n" + "=" * 70)
    print("FULL S₋ EIGENVALUE = α × (e₂+φ)/(1-K*)")
    print("=" * 70)

    full_Sminus = alpha_floor * lambda_DS_Sminus
    print(f"\n  full S₋ = {nstr(full_Sminus, 40)}")
    print(f"  numerical = {nstr(lambda_Sminus_numerical, 40)}")
    print(f"  difference = {nstr(fabs(full_Sminus - lambda_Sminus_numerical), 15)}")

    # ============================================================
    # NOW: what about the s₁ direction (the other eigenvalue)?
    # ============================================================
    print("\n" + "=" * 70)
    print("THE 3×3 BLOCK: LARGEST EIGENVALUE")
    print("=" * 70)

    # The two single-site eigenvalues are ~0.283 and ~0.281
    # One is the S₋ mode, the other is from the 3×3 block
    # Let's identify them

    # S₋ eigenvalue:
    print(f"  S₋ eigenvalue = {nstr(lambda_Sminus_numerical, 30)}")

    # Full 4×4 eigenvalues
    J4 = matrix(4, 4)
    for j in range(4):
        sp = list(m_star); sm = list(m_star)
        sp[j] += eps; sm[j] -= eps
        fp, _ = ds_step(sp, e_star)
        fm, _ = ds_step(sm, e_star)
        for i in range(4):
            J4[i,j] = (fp[i] - fm[i]) / (2*eps)

    evals4, _ = eig(J4)
    evals4 = sorted([chop(e) for e in evals4], key=lambda x: -float(fabs(x)))

    print(f"  Full eigenvalues: {', '.join(nstr(fabs(e), 15) for e in evals4)}")
    print(f"  The S₋ eigenvalue matches λ₁ = {nstr(fabs(evals4[1]), 15)}")
    print(f"  The 3×3 dominant eigenvalue = λ₀ = {nstr(fabs(evals4[0]), 15)}")

    # ============================================================
    # HOW DO THE SINGLE-SITE EIGENVALUES BECOME THE 4 COUPLED EIGENVALUES?
    # ============================================================
    print("\n" + "=" * 70)
    print("SINGLE-SITE → COUPLED: HOW COUPLING SPLITS THE EIGENVALUES")
    print("=" * 70)

    # Single-site: λ_ss ≈ 0.283 (3×3) and λ_S₋ ≈ 0.281 (1×1)
    # Coupled: {0.502, 0.474, 0.353, 0.334}
    #
    # Under coupling g = K* = 7/30:
    # Antisymmetric (A-B): enhancement
    # Symmetric (A+B): suppression

    lambda_ss_0 = float(fabs(evals4[0]))  # 0.2829
    lambda_ss_1 = float(fabs(evals4[1]))  # 0.2813

    coupled = [0.50217, 0.47448, 0.35265, 0.33442]

    print(f"\n  Single-site eigenvalues: {lambda_ss_0:.6f}, {lambda_ss_1:.6f}")
    print(f"  Coupled eigenvalues: {coupled}")

    # Check: products
    print(f"\n  λ_ss_0 × λ_ss_1 = {lambda_ss_0 * lambda_ss_1:.8f}")
    print(f"  Coupled λ₀ × λ₃ = {coupled[0]*coupled[3]:.8f}")
    print(f"  Coupled λ₁ × λ₂ = {coupled[1]*coupled[2]:.8f}")

    # Ratio: coupled products / single-site product
    ss_prod = lambda_ss_0 * lambda_ss_1
    cp_03 = coupled[0]*coupled[3]
    cp_12 = coupled[1]*coupled[2]
    print(f"\n  λ₀λ₃ / (λ_ss₀ × λ_ss₁) = {cp_03/ss_prod:.6f}")
    print(f"  λ₁λ₂ / (λ_ss₀ × λ_ss₁) = {cp_12/ss_prod:.6f}")

    # Check: sums
    print(f"\n  λ₀ + λ₃ = {coupled[0]+coupled[3]:.6f}")
    print(f"  λ₁ + λ₂ = {coupled[1]+coupled[2]:.6f}")
    print(f"  2×λ_ss₀ = {2*lambda_ss_0:.6f}")
    print(f"  2×λ_ss₁ = {2*lambda_ss_1:.6f}")

    # So the coupled eigenvalues are NOT simply λ_ss ± g
    # Let's check what the coupling does more precisely
    print(f"\n  λ₀ - λ₃ = {coupled[0]-coupled[3]:.6f}")
    print(f"  λ₁ - λ₂ = {coupled[1]-coupled[2]:.6f}")
    print(f"  (λ₀-λ₃)/(λ₀+λ₃) = {(coupled[0]-coupled[3])/(coupled[0]+coupled[3]):.6f}")
    print(f"  (λ₁-λ₂)/(λ₁+λ₂) = {(coupled[1]-coupled[2])/(coupled[1]+coupled[2]):.6f}")
    print(f"  K* = {7/30:.6f}")

    # Check if the splits are related to K*
    split_0 = (coupled[0]-coupled[3])/(coupled[0]+coupled[3])
    split_1 = (coupled[1]-coupled[2])/(coupled[1]+coupled[2])
    print(f"\n  Relative split 0: {split_0:.6f}")
    print(f"  Relative split 1: {split_1:.6f}")
    print(f"  K* = {7/30:.6f}")
    print(f"  split₀/K* = {split_0/(7/30):.6f}")
    print(f"  split₁/K* = {split_1/(7/30):.6f}")


if __name__ == "__main__":
    main()
