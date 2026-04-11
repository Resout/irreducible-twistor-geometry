#!/usr/bin/env python3
"""
Investigate Born(DS(m*)) — the "vacuum depth."

This is the Born probability of the DS combination output BEFORE
floor enforcement. It's a structural property of the fixed point,
not an intermediate step. It measures "how close to nothing" the
self-consistency requires.

Key questions:
1. What is Born_DS in terms of H and K*?
2. Does it relate to Λ (cosmological constant)?
3. What is θ_DS / θ* — the compression ratio?
4. What are all the derived quantities in this region?
"""

from mpmath import mp, mpf, matrix, sqrt, log, fabs, nstr, pi, exp, power, ln
from mpmath import findroot

mp.dps = 80
H = 3
FLOOR = mpf(1) / mpf(H**3)
K_star = mpf(7) / mpf(30)


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
    print("THE VACUUM DEPTH: Born(DS(m*))")
    print("=" * 70)

    m_star, e_star = find_equilibrium()
    s1, s2, s3, theta = m_star
    e1, e2, e3, phi = e_star

    print(f"\nEquilibrium:")
    print(f"  m* = ({nstr(s1,25)}, {nstr(s2,25)}, {nstr(s3,25)}, {nstr(theta,25)})")
    print(f"  e* = ({nstr(e1,25)}, {nstr(e2,25)}, {nstr(e3,25)}, {nstr(phi,25)})")
    print(f"  Born(m*) = {nstr(born_prob(m_star), 40)}")

    # DS combination output
    m_ds, K_val = ds_combine(m_star, e_star)
    s1_ds, s2_ds, s3_ds, theta_ds = m_ds

    Born_DS = born_prob(m_ds)
    print(f"\nDS output (before floor):")
    print(f"  m_DS = ({nstr(s1_ds,25)}, {nstr(s2_ds,25)}, {nstr(s3_ds,25)}, {nstr(theta_ds,25)})")
    print(f"  Born(m_DS) = {nstr(Born_DS, 50)}")
    print(f"  K = {nstr(K_val, 30)}")

    # ============================================================
    # ANALYTIC DECOMPOSITION OF Born_DS
    # ============================================================
    print("\n" + "=" * 70)
    print("ANALYTIC DECOMPOSITION")
    print("=" * 70)

    # θ_DS = θ*·φ* / (1-K*)
    theta_ds_analytic = theta * phi / (1 - K_star)
    print(f"\n  θ_DS = θ*·φ*/(1-K*)")
    print(f"    = {nstr(theta,15)} × {nstr(phi,15)} / {nstr(1-K_star,15)}")
    print(f"    = {nstr(theta_ds_analytic, 40)}")
    print(f"    check: {nstr(theta_ds, 40)}")

    # L2 norm squared of m_DS
    L2_ds = sum(x**2 for x in m_ds)
    print(f"\n  ||m_DS||² = {nstr(L2_ds, 40)}")

    # Born_DS = θ_DS² / ||m_DS||²
    # = (θ*φ*)² / [(1-K*)² × ||m_DS||²]
    print(f"\n  Born_DS = (θ*φ*)² / [(1-K*)² × ||m_DS||²]")
    print(f"         = {nstr(theta**2 * phi**2, 20)} / [{nstr((1-K_star)**2, 15)} × {nstr(L2_ds, 15)}]")
    print(f"         = {nstr(Born_DS, 40)}")

    # What is ||m_DS||²?
    # s_i_DS = (s_i*e_i + s_i*φ + θ*e_i) / (1-K)
    # For symmetric case s2=s3, e2=e3:
    S_ds = s1_ds + s2_ds + s3_ds
    Sq_ds = s1_ds**2 + s2_ds**2 + s3_ds**2

    print(f"\n  S_DS = {nstr(S_ds, 30)}")
    print(f"  Sq_DS = {nstr(Sq_ds, 30)}")
    print(f"  θ_DS² = {nstr(theta_ds**2, 30)}")
    print(f"  ||m_DS||² = Sq_DS + θ_DS² = {nstr(Sq_ds + theta_ds**2, 30)}")

    # ============================================================
    # BORN_DS IN TERMS OF FRAMEWORK NUMBERS
    # ============================================================
    print("\n" + "=" * 70)
    print("Born_DS vs FRAMEWORK NUMBERS")
    print("=" * 70)

    B = Born_DS
    print(f"\n  Born_DS = {nstr(B, 50)}")
    print(f"  1/Born_DS = {nstr(1/B, 30)}")

    # Test against many expressions
    tests = {
        # Simple
        "1/H⁶": mpf(1)/H**6,
        "1/H⁵": mpf(1)/H**5,
        "K*/H³": K_star/H**3,
        "K*²/H²": K_star**2/H**2,
        "(K*)³": K_star**3,
        "K*²/(H²-1)": K_star**2/(H**2-1),

        # Products involving θ* and φ*
        "θ*²φ*²": theta**2 * phi**2,

        # Born floor squared
        "(1/H³)²": FLOOR**2,
        "1/(H³(H³-1))": mpf(1)/(H**3*(H**3-1)),

        # Combinations
        "K*²/H³": K_star**2/H**3,
        "K*/(H⁴-1)": K_star/(H**4-1),
        "K*³/H": K_star**3/H,
        "1/(H²(H³-1))": mpf(1)/(H**2*(H**3-1)),
        "1/(H(H²-1)(H³-1))": mpf(1)/(H*(H**2-1)*(H**3-1)),

        # Floor × K*
        "K*/H³": K_star/H**3,
        "K*²/(H³-1)": K_star**2/(H**3-1),

        # Cosmological-scale candidates
        "exp(-H⁴)": exp(-mpf(H**4)),
        "exp(-H³)": exp(-mpf(H**3)),
        "K*^H": K_star**H,
        "K*^(H+1)": K_star**(H+1),
        "(1/H³)^H": FLOOR**H,
        "(K*/H)^H": (K_star/H)**H,

        # Deeper
        "1/(H²×h(E₈))": mpf(1)/(H**2*30),
        "K*/(H×h(E₈))": K_star/(H*30),
        "1/1236": mpf(1)/1236,
        "1/1237": mpf(1)/1237,
    }

    results = []
    for name, val in tests.items():
        if float(val) > 0:
            ratio = B / val
            err = float(fabs(ratio - 1) * 100)
            results.append((err, name, val, ratio))

    results.sort()
    print("\n  Closest matches:")
    for err, name, val, ratio in results[:15]:
        print(f"    {name:30s} = {nstr(val, 15):>20s}  ratio={nstr(ratio, 8):>12s}  err={err:.4f}%")

    # ============================================================
    # RECIPROCAL INVESTIGATION
    # ============================================================
    print("\n" + "=" * 70)
    print("1/Born_DS INVESTIGATION")
    print("=" * 70)

    invB = 1/B
    print(f"\n  1/Born_DS = {nstr(invB, 40)}")

    inv_tests = {
        "H⁶": mpf(H**6),
        "H⁵": mpf(H**5),
        "H³/K*": H**3/K_star,
        "H³(H³-1)": mpf(H**3*(H**3-1)),
        "H²(H³-1)": mpf(H**2*(H**3-1)),
        "H⁴/K*": mpf(H**4)/K_star,
        "H³/K*²": mpf(H**3)/K_star**2,
        "30²/K*": mpf(900)/K_star,
        "H²×120": mpf(H**2*120),
        "H²×h(E₈)²": mpf(H**2*900),
        "(H²-1)/K*²": mpf(H**2-1)/K_star**2,
        "H(H²-1)(H³-1)": mpf(H*(H**2-1)*(H**3-1)),
    }

    results2 = []
    for name, val in inv_tests.items():
        ratio = invB / val
        err = float(fabs(ratio - 1) * 100)
        results2.append((err, name, val, ratio))

    results2.sort()
    print("\n  Closest matches for 1/Born_DS:")
    for err, name, val, ratio in results2[:10]:
        print(f"    {name:30s} = {nstr(val, 15):>20s}  ratio={nstr(ratio, 8):>12s}  err={err:.4f}%")

    # ============================================================
    # θ COMPRESSION RATIO
    # ============================================================
    print("\n" + "=" * 70)
    print("θ COMPRESSION RATIO")
    print("=" * 70)

    comp = theta_ds / theta
    print(f"\n  θ_DS / θ* = {nstr(comp, 40)}")
    print(f"  = φ* / (1-K*) = {nstr(phi/(1-K_star), 40)}")

    comp_tests = {
        "1/H": mpf(1)/H,
        "K*": K_star,
        "φ*": phi,
        "1/(H+1)": mpf(1)/(H+1),
        "K*/(H-1)": K_star/(H-1),
        "1/H²": mpf(1)/H**2,
        "2/(H(H+1))": mpf(2)/(H*(H+1)),
    }

    results3 = []
    for name, val in comp_tests.items():
        ratio = comp / val
        err = float(fabs(ratio - 1) * 100)
        results3.append((err, name, val, ratio))

    results3.sort()
    for err, name, val, ratio in results3[:5]:
        print(f"    θ_DS/θ* ≈ {name:20s} = {nstr(val, 15):>15s}  err={err:.4f}%")

    # ============================================================
    # THE VACUUM DEPTH AS AN ENERGY
    # ============================================================
    print("\n" + "=" * 70)
    print("THE VACUUM DEPTH AS AN ENERGY")
    print("=" * 70)

    print(f"\n  Born_DS = {nstr(B, 20)}")
    print(f"  Born_floor = 1/H³ = {nstr(FLOOR, 20)}")
    print(f"  Ratio floor/depth = {nstr(FLOOR/B, 20)}")
    print(f"  Ratio depth/floor = {nstr(B/FLOOR, 20)}")

    # The "nothing depth" — how far below the floor
    depth = FLOOR - B
    print(f"\n  Floor - depth = {nstr(depth, 20)}")
    print(f"  depth/floor = {nstr(depth/FLOOR, 20)}")
    print(f"  1 - depth/floor = {nstr(B/FLOOR, 20)}")

    # ln of the depth
    print(f"\n  ln(Born_DS) = {nstr(ln(B), 30)}")
    print(f"  ln(1/27) = {nstr(ln(FLOOR), 30)}")
    print(f"  ln(Born_DS)/ln(1/27) = {nstr(ln(B)/ln(FLOOR), 20)}")
    print(f"  = ratio of how many e-foldings below 1")

    ratio_logs = ln(B) / ln(FLOOR)
    print(f"\n  ln(B)/ln(floor) = {nstr(ratio_logs, 40)}")
    # This means B = floor^ratio_logs, or equivalently (1/27)^ratio

    # What is ratio_logs as a framework number?
    log_tests = {
        "H": mpf(H),
        "H-1": mpf(H-1),
        "H+1": mpf(H+1),
        "H²/H": mpf(H),
        "(H+1)/H×H": mpf(H+1),
        "H²-1": mpf(H**2-1),
        "H(H-1)": mpf(H*(H-1)),
        "2H-1": mpf(2*H-1),
        "ln(1/K*)/ln(1/H³)": ln(1/K_star)/ln(mpf(H**3)),
        "π": pi,
    }

    results4 = []
    for name, val in log_tests.items():
        ratio = ratio_logs / val
        err = float(fabs(ratio - 1) * 100)
        results4.append((err, name, val, ratio))

    results4.sort()
    print("\n  ln(Born_DS)/ln(floor) closest matches:")
    for err, name, val, ratio in results4[:5]:
        print(f"    {name:30s}  err={err:.4f}%")

    # ============================================================
    # BORN_DS AS A POWER OF THE FLOOR
    # ============================================================
    print("\n" + "=" * 70)
    print("Born_DS AS A POWER OF THE FLOOR")
    print("=" * 70)

    # B = (1/27)^p → p = ln(B)/ln(1/27)
    p = ln(B) / ln(FLOOR)
    print(f"\n  Born_DS = (1/H³)^p where p = {nstr(p, 40)}")
    print(f"  If p = H+1/H = 10/3: (1/27)^(10/3) = {nstr(FLOOR**(mpf(10)/3), 20)}")
    print(f"  If p = π: (1/27)^π = {nstr(FLOOR**pi, 20)}")

    # ============================================================
    # COSMOLOGICAL CONSTANT CANDIDATES
    # ============================================================
    print("\n" + "=" * 70)
    print("COSMOLOGICAL CONSTANT INVESTIGATION")
    print("=" * 70)

    print(f"\n  Observed: Λ ≈ 2.888 × 10⁻¹²² (Planck units)")
    print(f"  ln(Λ) ≈ -280.5")

    # The framework's fundamental scale is m(0++) ≈ 1.71 GeV
    # Planck mass ≈ 1.22 × 10¹⁹ GeV
    # m(0++)/M_Planck ≈ 1.4 × 10⁻¹⁹
    m_ratio = mpf("1.4e-19")

    # Λ_obs in units of m(0++)⁴:
    # Λ_obs ≈ 2.888e-122 × M_Pl⁴
    # In m(0++)⁴ units: Λ = 2.888e-122 × (M_Pl/m(0++))⁴
    # = 2.888e-122 × (1/1.4e-19)⁴ = 2.888e-122 × 2.6e75 = 7.5e-47

    Lambda_glueball = mpf("7.5e-47")
    print(f"\n  Λ in glueball units (m(0++)⁴): ≈ {nstr(Lambda_glueball, 5)}")

    # Now: can Born_DS raised to some power give this?
    # Born_DS ≈ 8.09e-4
    # Born_DS^n = 7.5e-47 → n = ln(7.5e-47)/ln(8.09e-4) = -106.2/-7.12 ≈ 14.9

    n_lambda = ln(Lambda_glueball) / ln(B)
    print(f"  Born_DS^n = Λ → n = {nstr(n_lambda, 15)}")

    # What about Born_DS^(H²+H+1) = Born_DS^13?
    print(f"\n  Born_DS^13 = {nstr(B**13, 10)}")
    print(f"  Born_DS^14 = {nstr(B**14, 10)}")
    print(f"  Born_DS^15 = {nstr(B**15, 10)}")

    # What about floor-related expressions?
    print(f"\n  (1/H³)^(H⁴) = (1/27)^81 = {nstr(FLOOR**81, 10)}")
    print(f"  (1/H³)^(H³+H) = (1/27)^30 = {nstr(FLOOR**30, 10)}")
    print(f"  (1/H³)^(d_bos) = (1/27)^26 = {nstr(FLOOR**26, 10)}")

    # The key structural candidate: exp(-1/Born_DS)
    print(f"\n  exp(-1/Born_DS) = exp(-{nstr(1/B, 8)}) = {nstr(exp(-1/B), 10)}")
    print(f"  exp(-H³/Born_DS) = {nstr(exp(-H**3/B), 10)}")

    # ============================================================
    # WHAT IS Born_DS × H³? (Born_DS / floor)
    # ============================================================
    print("\n" + "=" * 70)
    print("Born_DS / FLOOR = Born_DS × H³")
    print("=" * 70)

    ratio_bf = B * H**3
    print(f"\n  Born_DS × H³ = {nstr(ratio_bf, 40)}")
    print(f"  = Born_DS / (1/H³)")

    bf_tests = {
        "K*²/H": K_star**2/H,
        "K*²/(H-1)": K_star**2/(H-1),
        "(φ*/(1-K*))²": (phi/(1-K_star))**2,
        "θ_DS²/θ*²×Born(m*)": (theta_ds/theta)**2 * FLOOR,
        "(φ/(1-K))²×Born_floor": (phi/(1-K_star))**2 * FLOOR,
        "K*²": K_star**2,
        "1/(H²+1)": mpf(1)/(H**2+1),
    }

    for name, val in bf_tests.items():
        err = float(fabs(ratio_bf/val - 1)*100)
        if err < 10:
            print(f"  ≈ {name} = {nstr(val, 15)} (err {err:.4f}%)")

    # ============================================================
    # THE FULL ANALYTIC FORM
    # ============================================================
    print("\n" + "=" * 70)
    print("FULL ANALYTIC BORN_DS")
    print("=" * 70)

    # Born_DS = θ_DS² / ||m_DS||²
    # θ_DS = θ*φ*/(1-K*)
    # ||m_DS||² = Σ(s_i_DS²) + θ_DS²

    # s1_DS = (s1*e1 + s1*φ + θ*e1)/(1-K)
    #       = s1*(e1+φ)/(1-K) + θ*e1/(1-K)
    s1_ds_check = (s1*(e1+phi) + theta*e1) / (1-K_star)
    s2_ds_check = (s2*(e2+phi) + theta*e2) / (1-K_star)
    print(f"  s1_DS check: {nstr(s1_ds_check, 20)} vs {nstr(s1_ds, 20)}")
    print(f"  s2_DS check: {nstr(s2_ds_check, 20)} vs {nstr(s2_ds, 20)}")

    # So Born_DS = [θ*φ*/(1-K*)]² / Σ[(s_i*(e_i+φ)+θ*e_i)/(1-K*)]² + [θ*φ*/(1-K*)]²
    # = (θ*φ*)² / {Σ[s_i*(e_i+φ)+θ*e_i]² + (θ*φ*)²}

    num = (theta*phi)**2
    denom_parts = [(s1*(e1+phi)+theta*e1)**2,
                   (s2*(e2+phi)+theta*e2)**2,
                   (s3*(e3+phi)+theta*e3)**2,
                   (theta*phi)**2]
    denom_val = sum(denom_parts)
    Born_check = num / denom_val
    print(f"\n  Analytic Born_DS = {nstr(Born_check, 40)}")
    print(f"  Numerical Born_DS = {nstr(Born_DS, 40)}")

    # Factor out θ*² from numerator and denominator:
    # num = θ*²φ*²
    # denom = Σ[s_i(e_i+φ) + θ*e_i]² + θ*²φ*²
    # = Σ[s_i(e_i+φ)]² + 2Σ[s_i(e_i+φ)·θ*e_i] + Σ[θ*e_i]² + θ*²φ*²
    # = Σs_i²(e_i+φ)² + 2θ*Σs_i·e_i(e_i+φ) + θ*²[Σe_i² + φ²]
    # = Σs_i²(e_i+φ)² + 2θ*Σs_i·e_i(e_i+φ) + θ*²·||e||²_L2

    # So Born_DS = φ*² / {Σ(s_i/θ*)²(e_i+φ)² + 2Σ(s_i/θ*)(e_i)(e_i+φ) + ||e||²_L2}

    # The ratio s_i/θ* at the floor: s1/θ = s1*H³ (since Born=1/27)
    # Actually Born = θ²/||m||² = 1/27 means ||m||²/θ² = 27 = H³
    # But s_i/θ is just the ratio of components

    print(f"\n  s1*/θ* = {nstr(s1/theta, 20)}")
    print(f"  s2*/θ* = {nstr(s2/theta, 20)}")
    print(f"  = {nstr(s1/theta, 20)} (dominant)")

    # The dominant term in the denominator is s1²(e1+φ)²
    dom_frac = (s1*(e1+phi)+theta*e1)**2 / denom_val
    print(f"\n  Fraction from s1 term: {nstr(dom_frac*100, 10)}%")

    # ============================================================
    # DOES THE DEPTH SCALE WITH THE FLOOR?
    # ============================================================
    print("\n" + "=" * 70)
    print("BORN_DS vs FLOOR ACROSS DIFFERENT H VALUES")
    print("=" * 70)

    print("\n  Testing: is Born_DS / floor^p constant for some p?")
    print("  (This requires solving the equilibrium for different H values)")

    # For H=3, we have Born_DS ≈ 8.09e-4, floor = 1/27 ≈ 3.70e-2
    # Born_DS/floor ≈ 0.02184
    # Born_DS/floor² ≈ 0.590
    # Born_DS/floor³ ≈ 15.9

    print(f"\n  H=3:")
    print(f"    floor = 1/{H**3} = {nstr(FLOOR, 10)}")
    print(f"    Born_DS = {nstr(B, 10)}")
    print(f"    Born_DS/floor = {nstr(B/FLOOR, 10)}")
    print(f"    Born_DS/floor² = {nstr(B/FLOOR**2, 10)}")
    print(f"    √(Born_DS/floor) = {nstr(sqrt(B/FLOOR), 10)}")
    print(f"    Born_DS^(1/2) = {nstr(sqrt(B), 10)}")

    # ============================================================
    # THE RATIO floor/Born_DS
    # ============================================================
    print("\n" + "=" * 70)
    print("floor/Born_DS = THE AMPLIFICATION FACTOR")
    print("=" * 70)

    amp = FLOOR / B
    print(f"\n  floor/Born_DS = {nstr(amp, 40)}")
    print(f"  = how much the floor amplifies the Born probability")

    amp_tests = {
        "H³-1 = 26": mpf(26),
        "H³ = 27": mpf(27),
        "H²(H+1)/2": mpf(H**2*(H+1)/2),
        "H(H²-1) = 24": mpf(24),
        "H⁴/2": mpf(H**4)/2,
        "5!/(H-1) = 60": mpf(60),
        "H²(H-1)²": mpf(H**2*(H-1)**2),
        "H(H+1)²-1 = 47": mpf(47),
        "H⁴-H²+1 = 73": mpf(73),
        "(H²-1)²/H": mpf((H**2-1)**2/H),
        "H²×5 = 45": mpf(45),
        "2H³+H²+H = 66": mpf(66),
        "H⁴-H = 78": mpf(78),
    }

    results5 = []
    for name, val in amp_tests.items():
        ratio = amp / val
        err = float(fabs(ratio - 1) * 100)
        results5.append((err, name, val, ratio))

    results5.sort()
    print("\n  Closest matches for floor/Born_DS:")
    for err, name, val, ratio in results5[:8]:
        print(f"    {name:30s} = {nstr(val, 10):>10s}  ratio={nstr(ratio, 8):>12s}  err={err:.2f}%")

    # ============================================================
    # sqrt(Born_DS) — THE AMPLITUDE
    # ============================================================
    print("\n" + "=" * 70)
    print("√Born_DS — THE AMPLITUDE OF NOTHING")
    print("=" * 70)

    sqrtB = sqrt(B)
    print(f"\n  √Born_DS = {nstr(sqrtB, 40)}")

    sqrt_tests = {
        "K*/H": K_star/H,
        "1/H²": mpf(1)/H**2,
        "K*/(H+1)": K_star/(H+1),
        "1/(H(H+1))": mpf(1)/(H*(H+1)),
        "φ*/(1-K*)×√(1/H³)": phi/(1-K_star)*sqrt(FLOOR),
        "θ_DS/θ*×√Born_floor": comp * sqrt(FLOOR),
    }

    for name, val in sqrt_tests.items():
        err = float(fabs(sqrtB/val - 1)*100)
        if err < 15:
            print(f"  ≈ {name} = {nstr(val, 15)} (err {err:.3f}%)")

    # ============================================================
    # THE DEPTH NUMBER ITSELF
    # ============================================================
    print("\n" + "=" * 70)
    print("SUMMARY: THE VACUUM DEPTH NUMBER")
    print("=" * 70)

    print(f"\n  Born_DS = {nstr(B, 50)}")
    print(f"  1/Born_DS = {nstr(1/B, 30)}")
    print(f"  floor/Born_DS = {nstr(amp, 20)}")
    print(f"  ln(Born_DS) = {nstr(ln(B), 30)}")
    print(f"  ln(Born_DS)/ln(floor) = {nstr(p, 20)}")
    print(f"\n  Born_DS = floor^{nstr(p, 15)}")
    print(f"  = (1/H³)^{nstr(p, 15)}")


if __name__ == "__main__":
    main()
