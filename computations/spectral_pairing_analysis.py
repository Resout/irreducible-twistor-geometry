#!/usr/bin/env python3
"""
Spectral pairing analysis for the A₂ coupled Jacobian.

Questions to answer:
1. Is λ₀λ₃ = λ₁λ₂ exact, or approximate? (spectral pairing)
2. Why do rational power denominators q | H(H+1) = 12?
3. What is the block structure of the 8×8 Jacobian under symmetries?
"""

from mpmath import mp, mpf, matrix, sqrt, log, fabs, nstr, eig, pi
from mpmath import findroot, chop
import numpy as np

mp.dps = 80
H = 3
FLOOR = mpf(1) / mpf(H**3)


# ============================================================
# DS dynamics (from a2_eigenvalue_50digit.py)
# ============================================================

def ds_combine(m, e):
    s = m[:3]; theta = m[3]
    ev = e[:3]; phi = e[3]
    s_pre = [s[i]*ev[i] + s[i]*phi + theta*ev[i] for i in range(3)]
    theta_pre = theta * phi
    total_pre = sum(s_pre) + theta_pre
    K = mpf(1) - total_pre
    if fabs(mpf(1) - K) < mpf(10)**(-70):
        return list(m), K
    denom = mpf(1) - K
    out = [sp / denom for sp in s_pre] + [theta_pre / denom]
    return out, K

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
    A_coeff = mpf(26) * S**2 - Sq
    B_coeff = mpf(2) * Sq
    C_coeff = -Sq
    disc = B_coeff**2 - 4*A_coeff*C_coeff
    if disc < 0: return list(m)
    t1 = (-B_coeff + sqrt(disc)) / (2*A_coeff)
    t2 = (-B_coeff - sqrt(disc)) / (2*A_coeff)
    candidates = [t for t in [t1, t2] if mpf(0) < t < mpf(1)]
    if not candidates: return list(m)
    t = min(candidates, key=lambda x: fabs(x - m[3]))
    alpha = (mpf(1) - t) / S
    return [m[i]*alpha for i in range(3)] + [t]

def ds_step(m, e):
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds), K


# ============================================================
# Find equilibrium
# ============================================================

def find_single_site_equilibrium():
    from scipy.optimize import fsolve
    def equations_rough(params):
        s1, theta, w1, phi = params
        s2 = (1.0 - s1 - theta) / 2.0
        w2 = (1.0 - w1 - phi) / 2.0
        m = [s1, s2, s2, theta]; e = [w1, w2, w2, phi]
        L2m = s1**2 + 2*s2**2 + theta**2
        eq1 = theta**2 / L2m - 1.0/27.0
        L2e = w1**2 + 2*w2**2 + phi**2
        eq2 = phi**2 / L2e - 1.0/27.0
        K = s1*w2 + s1*w2 + s2*w1 + s2*w2 + s2*w1 + s2*w2
        eq3 = K - 7.0/30.0
        m_mp = [mpf(x) for x in m]; e_mp = [mpf(x) for x in e]
        m_out, _ = ds_step(m_mp, e_mp)
        eq4 = float(m_out[0]) - s1
        return [eq1, eq2, eq3, eq4]
    sol = fsolve(equations_rough, [0.787, 0.155, 0.631, 0.129])

    def equations_mp(*params):
        s1, theta, w1, phi = params
        s2 = (mpf(1) - s1 - theta) / 2
        w2 = (mpf(1) - w1 - phi) / 2
        m = [s1, s2, s2, theta]; e = [w1, w2, w2, phi]
        eq1 = theta**2 / (s1**2 + 2*s2**2 + theta**2) - FLOOR
        eq2 = phi**2 / (w1**2 + 2*w2**2 + phi**2) - FLOOR
        K = s1*w2 + s1*w2 + s2*w1 + s2*w2 + s2*w1 + s2*w2
        eq3 = K - mpf(7)/30
        m_out, _ = ds_step(m, e)
        eq4 = m_out[0] - s1
        return [eq1, eq2, eq3, eq4]

    x = findroot(equations_mp, [mpf(str(v)) for v in sol])
    s1, theta, w1, phi = x
    s2 = (mpf(1) - s1 - theta) / 2
    w2 = (mpf(1) - w1 - phi) / 2
    return [s1, s2, s2, theta], [w1, w2, w2, phi]


def coupled_evidence(m1, m2, e0, g):
    uniform = [mpf(1)/4]*4
    e1 = [e0[i] - g * (m2[i] - uniform[i]) for i in range(4)]
    e1 = [max(x, mpf(10)**(-50)) for x in e1]
    total = sum(e1)
    return [x / total for x in e1]

def coupled_A2_step(state, e0, g):
    m1 = state[:4]; m2 = state[4:]
    e1 = coupled_evidence(m1, m2, e0, g)
    e2 = coupled_evidence(m2, m1, e0, g)
    m1_new, _ = ds_step(m1, e1)
    m2_new, _ = ds_step(m2, e2)
    return m1_new + m2_new


# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():
    print("=" * 70)
    print("SPECTRAL PAIRING AND q|12 ANALYSIS")
    print("=" * 70)

    # Find equilibrium
    print("\n1. Finding equilibrium...")
    m_star, e_star = find_single_site_equilibrium()
    print(f"   m* = [{', '.join(nstr(x, 15) for x in m_star)}]")
    print(f"   e* = [{', '.join(nstr(x, 15) for x in e_star)}]")

    # Coupled fixed point
    g = mpf(7) / 30
    state = list(m_star) + list(m_star)
    for i in range(10000):
        state_new = coupled_A2_step(state, e_star, g)
        diff = max(fabs(state_new[j] - state[j]) for j in range(8))
        if diff < mpf(10)**(-70):
            print(f"   Coupled fixed point converged at step {i}")
            break
        state = state_new

    # 8×8 Jacobian
    print("\n2. Computing 8×8 Jacobian...")
    eps = mpf(10)**(-35)
    J = matrix(8, 8)
    for j in range(8):
        sp = list(state_new); sm = list(state_new)
        sp[j] += eps; sm[j] -= eps
        fp = coupled_A2_step(sp, e_star, g)
        fm = coupled_A2_step(sm, e_star, g)
        for i in range(8):
            J[i, j] = (fp[i] - fm[i]) / (2*eps)

    # ============================================================
    # SYMMETRY DECOMPOSITION
    # ============================================================
    print("\n3. Symmetry decomposition of the 8×8 Jacobian")

    # Node-swap symmetry: P swaps indices 0-3 with 4-7
    # Extract blocks: J = [A B; C D] where A,B,C,D are 4x4
    A_blk = matrix(4, 4); B_blk = matrix(4, 4)
    C_blk = matrix(4, 4); D_blk = matrix(4, 4)
    for i in range(4):
        for j in range(4):
            A_blk[i,j] = J[i,j]
            B_blk[i,j] = J[i,j+4]
            C_blk[i,j] = J[i+4,j]
            D_blk[i,j] = J[i+4,j+4]

    # Check A ≈ D and B ≈ C (node symmetry)
    AD_diff = max(fabs(A_blk[i,j] - D_blk[i,j]) for i in range(4) for j in range(4))
    BC_diff = max(fabs(B_blk[i,j] - C_blk[i,j]) for i in range(4) for j in range(4))
    print(f"   ||A - D|| = {nstr(AD_diff, 5)}  (should be ~0 for node symmetry)")
    print(f"   ||B - C|| = {nstr(BC_diff, 5)}  (should be ~0 for node symmetry)")

    # Symmetric sector: A + B
    # Antisymmetric sector: A - B
    Jsym = matrix(4, 4); Jasym = matrix(4, 4)
    for i in range(4):
        for j in range(4):
            Jsym[i,j] = A_blk[i,j] + B_blk[i,j]
            Jasym[i,j] = A_blk[i,j] - B_blk[i,j]

    evals_sym, _ = eig(Jsym)
    evals_asym, _ = eig(Jasym)

    print(f"\n   Symmetric sector eigenvalues (A+B):")
    for k, ev in enumerate(sorted(evals_sym, key=lambda x: -float(fabs(x)))):
        v = chop(ev)
        print(f"     μ_{k} = {nstr(v, 30)}  |μ| = {nstr(fabs(v), 30)}")

    print(f"\n   Antisymmetric sector eigenvalues (A-B):")
    for k, ev in enumerate(sorted(evals_asym, key=lambda x: -float(fabs(x)))):
        v = chop(ev)
        print(f"     ν_{k} = {nstr(v, 30)}  |ν| = {nstr(fabs(v), 30)}")

    # Full eigenvalues for comparison
    evals_full, evecs_full = eig(J)
    evals_sorted = sorted(evals_full, key=lambda x: -float(fabs(x)))
    evals_sorted = [chop(e) for e in evals_sorted]

    print(f"\n   Full 8×8 eigenvalues:")
    for k, ev in enumerate(evals_sorted):
        print(f"     λ_{k} = {nstr(ev, 40)}  |λ| = {nstr(fabs(ev), 40)}")

    # ============================================================
    # SPECTRAL PAIRING TEST
    # ============================================================
    print("\n" + "=" * 70)
    print("4. SPECTRAL PAIRING: λ₀λ₃ vs λ₁λ₂")
    print("=" * 70)

    # Get top 4 by magnitude (these are the physical eigenvalues)
    mags = [(fabs(ev), ev) for ev in evals_sorted]
    top4 = [ev for _, ev in mags[:4]]
    lam = [fabs(ev) for ev in top4]

    print(f"\n   Top 4 eigenvalues by magnitude:")
    for k in range(4):
        print(f"     |λ_{k}| = {nstr(lam[k], 50)}")

    prod_03 = lam[0] * lam[3]
    prod_12 = lam[1] * lam[2]
    pairing_ratio = prod_03 / prod_12
    pairing_gap = fabs(prod_03 - prod_12) / prod_12

    print(f"\n   λ₀·λ₃ = {nstr(prod_03, 50)}")
    print(f"   λ₁·λ₂ = {nstr(prod_12, 50)}")
    print(f"   ratio  = {nstr(pairing_ratio, 50)}")
    print(f"   gap    = {nstr(pairing_gap * 100, 20)}%")

    # Also check: are these from different sectors?
    print(f"\n   Check which sector each eigenvalue belongs to:")
    sym_mags = sorted([fabs(chop(e)) for e in evals_sym], reverse=True)
    asym_mags = sorted([fabs(chop(e)) for e in evals_asym], reverse=True)

    for k in range(4):
        # Find closest match in each sector
        d_sym = min(fabs(lam[k] - s) for s in sym_mags)
        d_asym = min(fabs(lam[k] - a) for a in asym_mags)
        sector = "SYM" if float(d_sym) < float(d_asym) else "ANTI"
        print(f"     |λ_{k}| = {nstr(lam[k], 15)} → {sector} (d_sym={nstr(d_sym,5)}, d_asym={nstr(d_asym,5)})")

    # ============================================================
    # CHARACTERISTIC POLYNOMIAL ANALYSIS
    # ============================================================
    print("\n" + "=" * 70)
    print("5. CHARACTERISTIC POLYNOMIAL COEFFICIENTS")
    print("=" * 70)

    # For the physical 4×4 block (antisymmetric sector typically)
    # Compute σ₁ = Σλᵢ, σ₂ = Σλᵢλⱼ, σ₃ = Σλᵢλⱼλₖ, σ₄ = Πλᵢ
    l = lam  # top 4 magnitudes

    sigma1 = l[0] + l[1] + l[2] + l[3]
    sigma2 = (l[0]*l[1] + l[0]*l[2] + l[0]*l[3] +
              l[1]*l[2] + l[1]*l[3] + l[2]*l[3])
    sigma3 = (l[0]*l[1]*l[2] + l[0]*l[1]*l[3] +
              l[0]*l[2]*l[3] + l[1]*l[2]*l[3])
    sigma4 = l[0] * l[1] * l[2] * l[3]

    print(f"\n   σ₁ = {nstr(sigma1, 40)}")
    print(f"   σ₂ = {nstr(sigma2, 40)}")
    print(f"   σ₃ = {nstr(sigma3, 40)}")
    print(f"   σ₄ = {nstr(sigma4, 40)}")

    # Test: is σ₄ = (λ₀λ₃)²? (consequence of exact pairing)
    print(f"\n   σ₄ = {nstr(sigma4, 40)}")
    print(f"   (λ₀λ₃)² = {nstr(prod_03**2, 40)}")
    print(f"   (λ₁λ₂)² = {nstr(prod_12**2, 40)}")

    # Test the discriminant of the pairing
    # If we write λᵢ = a_i * b_i where a is from one factor and b from another
    # The pairing λ₀λ₃ = λ₁λ₂ means the eigenvalues satisfy a "trace-det" relation

    # ============================================================
    # s₂ ↔ s₃ SYMMETRY IN EACH 4×4 BLOCK
    # ============================================================
    print("\n" + "=" * 70)
    print("6. s₂ ↔ s₃ SYMMETRY DECOMPOSITION")
    print("=" * 70)

    # The 4-component state is (s₁, s₂, s₃, θ)
    # s₂ ↔ s₃ symmetry means Jasym has block structure in the basis:
    #   (s₁, (s₂+s₃)/√2, θ, (s₂-s₃)/√2)
    # = 3×3 symmetric + 1×1 antisymmetric

    # Change of basis matrix
    T = matrix(4, 4)
    T[0,0] = 1  # s₁ → s₁
    T[1,1] = mpf(1)/sqrt(2); T[1,2] = mpf(1)/sqrt(2)  # (s₂+s₃)/√2
    T[2,3] = 1  # θ → θ
    T[3,1] = mpf(1)/sqrt(2); T[3,2] = -mpf(1)/sqrt(2)  # (s₂-s₃)/√2

    Tinv = T**(-1)

    for label, Jblock in [("Symmetric (A+B)", Jsym), ("Antisymmetric (A-B)", Jasym)]:
        J_rotated = Tinv * Jblock * T
        print(f"\n   {label} in symmetry-adapted basis:")
        print(f"   (basis: s₁, S₊, θ, S₋)")
        for i in range(4):
            row = [nstr(J_rotated[i,j], 8) for j in range(4)]
            print(f"     [{', '.join(row)}]")

        # Check off-diagonal coupling between 3×3 and 1×1 sectors
        coupling = max(fabs(J_rotated[3,0]), fabs(J_rotated[3,2]),
                      fabs(J_rotated[0,3]), fabs(J_rotated[2,3]),
                      fabs(J_rotated[1,3]), fabs(J_rotated[3,1]))
        print(f"   ||3×3 ↔ 1×1 coupling|| = {nstr(coupling, 10)}")

        # The 1×1 eigenvalue is J_rotated[3,3]
        print(f"   1×1 eigenvalue (S₋ sector) = {nstr(J_rotated[3,3], 30)}")

        # The 3×3 block
        J33 = matrix(3, 3)
        for i in range(3):
            for j in range(3):
                ii = [0,1,2][i]  # s₁, S₊, θ indices
                jj = [0,1,2][j]
                J33[i,j] = J_rotated[ii,jj]

        evals_33, _ = eig(J33)
        print(f"   3×3 eigenvalues:")
        for k, ev in enumerate(sorted(evals_33, key=lambda x: -float(fabs(x)))):
            v = chop(ev)
            print(f"     {nstr(v, 30)}")

    # ============================================================
    # MASS GAP PAIRING IN LOG SPACE
    # ============================================================
    print("\n" + "=" * 70)
    print("7. MASS GAP PAIRING (in log space)")
    print("=" * 70)

    Delta = [-log(l) for l in lam[:4]]
    sum_03 = Delta[0] + Delta[3]
    sum_12 = Delta[1] + Delta[2]

    print(f"\n   Δ₀ = {nstr(Delta[0], 30)}")
    print(f"   Δ₁ = {nstr(Delta[1], 30)}")
    print(f"   Δ₂ = {nstr(Delta[2], 30)}")
    print(f"   Δ₃ = {nstr(Delta[3], 30)}")
    print(f"\n   Δ₀ + Δ₃ = {nstr(sum_03, 40)}")
    print(f"   Δ₁ + Δ₂ = {nstr(sum_12, 40)}")
    print(f"   gap = {nstr(fabs(sum_03 - sum_12)/sum_12 * 100, 15)}%")

    # In log space: Δ₀+Δ₃ = -ln(λ₀λ₃), Δ₁+Δ₂ = -ln(λ₁λ₂)
    # So the pairing Δ₀+Δ₃ = Δ₁+Δ₂ iff λ₀λ₃ = λ₁λ₂
    print(f"\n   (Same as λ₀λ₃ = λ₁λ₂ test above)")

    # ============================================================
    # q | 12 INVESTIGATION
    # ============================================================
    print("\n" + "=" * 70)
    print("8. WHY q | H(H+1) = 12?")
    print("=" * 70)

    print("""
   The 4×4 Jacobian block has:
   - s₂↔s₃ symmetry → 3×3 (symmetric) + 1×1 (antisymmetric)
   - The 3×3 block acts on (s₁, S₊, θ) — the observable sector
   - The 1×1 block acts on S₋ — the colour difference mode

   Key structural numbers:
   - dim(SU(3)) = H²-1 = 8 generators
   - dim(SU(2)) = 3 generators
   - dim(U(1)) = 1 generator
   - Total: H²-1 + H + 1 = H(H+1) = 12

   Hypothesis: rational powers q | 12 arise because the eigenvalue
   spectrum is governed by the 3×3 block's action on a lattice
   graded by gauge generator counts. The denominators correspond to:
     q=1:  trivial
     q=2:  Z₂ (matter/antimatter, or chirality)
     q=3:  Z₃ ⊂ SU(3) (colour centre)
     q=4:  Z₄ = Z₂ × Z₂ (chirality × parity)
     q=6:  Z₆ = Z₃ × Z₂ (centre of SU(3)×SU(2))
     q=12: full Z₁₂ = generator lattice
""")

    # The 3×3 block's characteristic polynomial
    # If it factors over Z[1/12], that would explain q|12
    # Let's look at powers of the eigenvalues

    print("   Testing: do eigenvalue ratios simplify at q | 12?")
    for q in [1, 2, 3, 4, 6, 12]:
        for i in range(4):
            for j in range(i+1, 4):
                ratio = lam[i] / lam[j]
                rq = ratio**q
                # Check if rq is close to a simple rational p/q with small p
                rq_float = float(rq)
                # Try to identify as rational
                best_err = 1.0
                best_frac = ""
                for denom in range(1, 25):
                    numer = round(rq_float * denom)
                    if numer > 0:
                        err = abs(rq_float - numer/denom) / rq_float
                        if err < best_err:
                            best_err = err
                            best_frac = f"{numer}/{denom}"
                if best_err < 0.005:
                    print(f"     (λ_{i}/λ_{j})^{q} = {rq_float:.8f} ≈ {best_frac} (err {best_err*100:.4f}%)")

    # ============================================================
    # TRACE FORMULA TEST
    # ============================================================
    print("\n" + "=" * 70)
    print("9. TRACE AND DETERMINANT IDENTITIES")
    print("=" * 70)

    # For the 4×4 antisymmetric block:
    evals_a = sorted([chop(e) for e in evals_asym], key=lambda x: -float(fabs(x)))

    tr_a = sum(evals_a)
    det_a = evals_a[0] * evals_a[1] * evals_a[2] * evals_a[3]

    print(f"\n   Antisymmetric block (A-B):")
    print(f"   tr = {nstr(tr_a, 30)}")
    print(f"   det = {nstr(det_a, 30)}")

    # Check if determinant relates to K* or H
    print(f"\n   det vs framework constants:")
    print(f"   det / (K*)⁴ = {nstr(det_a / (mpf(7)/30)**4, 15)}")
    print(f"   det / (1/H³)² = {nstr(det_a / FLOOR**2, 15)}")
    print(f"   det × H⁴ = {nstr(det_a * 81, 15)}")
    print(f"   -ln(|det|) = {nstr(-log(fabs(det_a)), 30)}")
    print(f"   -ln(|det|) / ln(30/23) = {nstr(-log(fabs(det_a)) / log(mpf(30)/23), 15)}")

    # Similarly for symmetric block
    evals_s = sorted([chop(e) for e in evals_sym], key=lambda x: -float(fabs(x)))
    tr_s = sum(evals_s)
    det_s = evals_s[0] * evals_s[1] * evals_s[2] * evals_s[3]

    print(f"\n   Symmetric block (A+B):")
    print(f"   tr = {nstr(tr_s, 30)}")
    print(f"   det = {nstr(det_s, 30)}")

    # The spectral pairing in each sector
    print(f"\n   Spectral pairing in each sector:")
    for label, evs in [("Sym", evals_s), ("Anti", evals_a)]:
        mags = sorted([fabs(e) for e in evs], reverse=True)
        p03 = mags[0] * mags[3]
        p12 = mags[1] * mags[2]
        ratio = p03 / p12 if float(p12) > 1e-50 else mpf(0)
        print(f"     {label}: |e₀||e₃| / |e₁||e₂| = {nstr(ratio, 20)}")


if __name__ == "__main__":
    main()
