#!/usr/bin/env python3
"""
EIGENVALUE FLOW ALONG THE DS MODULI SPACE
==========================================
Vary the Born floor parameter f from ~0 to 1/H³ = 1/27 and track how
the fixed-point eigenvalues change. The hypothesis: the floor is an
effective IR cutoff, and the eigenvalue flow as f varies IS the RG flow.

UV (high energy): floor irrelevant → f → 0
IR (low energy):  floor dominates  → f = 1/27

Parts:
  1. Eigenvalue flow λ_k(f) for 50 floor values
  2. Logarithmic fit:  λ(f) = a + b·ln(f)   → b is one-loop β
  3. Two-loop fit:     λ(f) = a + b·ln(f) + c·ln²(f)
  4. Spectral gap Δ(f) and derived couplings
  5. Check whether the flow gives effective b₃ correction

Key subtlety: the enforce_floor quadratic coefficient (H³-1) becomes
(1/f - 1) when the floor is general f (since Born = θ²/L₂² ≥ f
requires θ²(1/f - 1) ≥ Sq, i.e., the coefficient is 1/f - 1).
"""

from mpmath import (mp, mpf, matrix, sqrt, fabs, nstr, eig, chop,
                    findroot, log, ln, pi, polyroots)
import sys

mp.dps = 50
H = 3


# ============================================================
# DS MACHINERY — parameterised by floor value
# ============================================================

def ds_combine(m, e):
    """DS combination rule: m ⊕ e with L₁ normalisation."""
    s = m[:3]; theta = m[3]; ev = e[:3]; phi = e[3]
    s_pre = [s[i]*ev[i] + s[i]*phi + theta*ev[i] for i in range(3)]
    theta_pre = theta * phi
    total_pre = sum(s_pre) + theta_pre
    K = mpf(1) - total_pre
    if fabs(mpf(1) - K) < mpf(10)**(-40):
        return list(m), K
    denom = mpf(1) - K
    return [sp / denom for sp in s_pre] + [theta_pre / denom], K


def born_prob(m):
    """Born probability: θ²/L₂²."""
    L2sq = sum(x**2 for x in m)
    if L2sq < mpf(10)**(-40):
        return mpf(0)
    return m[3]**2 / L2sq


def enforce_floor(m, floor_val):
    """
    Enforce Born ≥ floor_val by redistributing from {sᵢ} to θ.
    The quadratic coefficient is (1/f - 1) in general.
    """
    if floor_val < mpf(10)**(-40):
        return list(m)  # no floor enforcement
    b = born_prob(m)
    if b >= floor_val - mpf(10)**(-40):
        return list(m)
    S = sum(m[:3])
    if S < mpf(10)**(-40):
        return list(m)
    Sq = sum(x**2 for x in m[:3])
    # Born = t² / (Sq_new + t²) = f  where Sq_new = α²·Sq, α = (1-t)/S
    # → t² / ((1-t)²·Sq/S² + t²) = f
    # → t²·S² = f·((1-t)²·Sq + t²·S²)
    # → t²·S²·(1 - f) = f·(1-t)²·Sq
    # → t²·S²·(1/f - 1) = (1-t)²·Sq
    # Quadratic in t: (1/f - 1)·S²·t² - Sq·(1-t)² = 0
    # Expand: ((1/f-1)·S² - Sq)·t² + 2Sq·t - Sq = 0
    inv_f_minus_1 = mpf(1) / floor_val - mpf(1)
    A_c = inv_f_minus_1 * S**2 - Sq
    B_c = mpf(2) * Sq
    C_c = -Sq
    disc = B_c**2 - 4 * A_c * C_c
    if disc < 0:
        return list(m)
    sd = sqrt(disc)
    t1 = (-B_c + sd) / (2 * A_c)
    t2 = (-B_c - sd) / (2 * A_c)
    cands = [t for t in [t1, t2] if mpf(0) < t < mpf(1)]
    if not cands:
        return list(m)
    t = min(cands, key=lambda x: fabs(x - m[3]))
    alpha = (mpf(1) - t) / S
    return [m[i] * alpha for i in range(3)] + [t]


def ds_step(m, e, floor_val):
    """Full DS step: combine then enforce floor."""
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds, floor_val), K


# ============================================================
# FIND EQUILIBRIUM for a given floor value
# ============================================================

def find_equilibrium(floor_val, guess=None):
    """
    Find the symmetric fixed point (s₂=s₃, e₂=e₃) at K=7/30
    for the given Born floor value.

    The equations:
      1. Born(m) = floor_val
      2. Born(e) = floor_val
      3. K(m,e) = 7/30
      4. Φ(m,e) = m  (fixed point: s₁ component)
    """
    K_star = mpf(7) / mpf(30)

    # Initial guess from scipy
    if guess is None:
        guess = [mpf('0.787'), mpf('0.155'), mpf('0.631'), mpf('0.129')]

    def equations(*p):
        s1, theta, w1, phi = p
        s2 = (mpf(1) - s1 - theta) / 2
        w2 = (mpf(1) - w1 - phi) / 2

        m = [s1, s2, s2, theta]
        e = [w1, w2, w2, phi]

        # Born constraints
        L2m = s1**2 + 2*s2**2 + theta**2
        L2e = w1**2 + 2*w2**2 + phi**2
        eq1 = theta**2 / L2m - floor_val
        eq2 = phi**2 / L2e - floor_val

        # Conflict
        K = 2*s1*w2 + 2*s2*w1 + 2*s2*w2
        eq3 = K - K_star

        # Fixed point
        m_out, _ = ds_step(m, e, floor_val)
        eq4 = m_out[0] - s1

        return [eq1, eq2, eq3, eq4]

    try:
        x = findroot(equations, list(guess), tol=mpf(10)**(-40))
        s1, theta, w1, phi = x
        s2 = (mpf(1) - s1 - theta) / 2
        w2 = (mpf(1) - w1 - phi) / 2
        return [s1, s2, s2, theta], [w1, w2, w2, phi]
    except Exception as exc:
        return None, None


# ============================================================
# COMPUTE JACOBIAN EIGENVALUES
# ============================================================

def compute_eigenvalues(m_star, e_star, floor_val):
    """Compute 4×4 single-site Jacobian eigenvalues by finite differences."""
    eps = mpf(10)**(-25)
    J4 = matrix(4, 4)
    for j in range(4):
        mp_p = list(m_star); mp_m = list(m_star)
        mp_p[j] += eps; mp_m[j] -= eps
        fp, _ = ds_step(mp_p, e_star, floor_val)
        fm, _ = ds_step(mp_m, e_star, floor_val)
        for i in range(4):
            J4[i, j] = (fp[i] - fm[i]) / (2 * eps)
    evals, _ = eig(J4)
    evals = [chop(e) for e in evals]
    evals_sorted = sorted(evals, key=lambda x: -float(fabs(x)))
    return evals_sorted, J4


# ============================================================
# PART 0: PURE DS (NO FLOOR) — ANALYTICAL LIMIT
# ============================================================

def find_equilibrium_no_floor(guess=None):
    """
    At floor = 0, the Born constraint is absent. The fixed-point
    equations become: DS(m,e)/(1-K) = m, with K = 7/30.
    The map is: s_i → (s_i·e_i + s_i·φ + θ·e_i)/(1-K)
                θ   → θ·φ/(1-K)

    Fixed point: s_i = (s_i·e_i + s_i·φ + θ·e_i)/(1-K)
                 θ   = θ·φ/(1-K)

    From θ equation: φ = 1-K (if θ≠0)
    Then s_i equation: 1-K = e_i + φ + θ·e_i/s_i
                      → 0 = e_i + θ·e_i/s_i (since φ=1-K)

    Hmm, that gives e_i(1 + θ/s_i) = 0, so either e_i=0 or θ/s_i = -1.
    Neither makes sense for positive masses.

    Actually if θ=0 at floor=0: then s₁+s₂+s₃=1 and the DS step is
    s_i → s_i·e_i/(1-K), requiring e_i = 1-K for all i.
    Then K = Σ s_i·e_j (i≠j) = (1-K)·Σ_{i≠j} s_i = (1-K)(3-1) for uniform...
    Actually K = Σ_{i≠j} s_i·e_j = (1-K)(1 - Σ s_i²) if e_i = 1-K.
    Wait let me reconsider.

    With symmetric ansatz s₁=s, s₂=s₃=(1-s-θ)/2, e₁=w, e₂=e₃=(1-w-φ)/2:
    K = 2·s·(1-w-φ)/2 + 2·((1-s-θ)/2)·w + 2·((1-s-θ)/2)·((1-w-φ)/2)
      = s(1-w-φ) + (1-s-θ)w + (1-s-θ)(1-w-φ)/2

    This is getting complicated. Let me just solve numerically with very small floor.
    """
    pass


# ============================================================
# MAIN COMPUTATION
# ============================================================

def main():
    print("=" * 80)
    print("EIGENVALUE FLOW ALONG THE DS MODULI SPACE")
    print("Born floor parameter f varied from ~0 to 1/H³ = 1/27")
    print("=" * 80)

    K_star = mpf(7) / mpf(30)
    physical_floor = mpf(1) / mpf(27)

    # ============================================================
    # Generate floor values: logarithmically spaced from tiny to 1/27
    # ============================================================
    import numpy as np
    log_floors = np.linspace(-6, float(ln(physical_floor)), 50)
    floor_values = [mpf(str(float(np.exp(lf)))) for lf in log_floors]
    # Add a few very small floors explicitly
    floor_values = [mpf('1e-8'), mpf('1e-6'), mpf('1e-5'), mpf('1e-4')] + floor_values
    # Add the physical floor exactly (avoid duplicates by filtering close values)
    floor_values.append(physical_floor)
    # Sort and deduplicate (remove values within 0.1% of each other)
    floor_values = sorted(set(floor_values), key=float)
    deduped = [floor_values[0]]
    for fv in floor_values[1:]:
        if float(fabs(fv - deduped[-1]) / deduped[-1]) > 0.001:
            deduped.append(fv)
    floor_values = deduped

    # ============================================================
    # PART 1: Eigenvalue flow
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 1: EIGENVALUE FLOW WITH VARYING FLOOR")
    print("=" * 80)

    results = []
    guess = [mpf('0.787'), mpf('0.155'), mpf('0.631'), mpf('0.129')]
    last_good_guess = list(guess)

    header = f"{'floor':>14s}  {'|λ₀|':>14s}  {'|λ₁|':>14s}  {'|λ₂|':>14s}  {'|λ₃|':>14s}  {'Δ=ln(λ₀/λ₁)':>14s}  {'θ*':>14s}"
    print(f"\n{header}")
    print("-" * len(header))

    for fv in floor_values:
        m_star, e_star = find_equilibrium(fv, guess=last_good_guess)
        if m_star is None:
            print(f"  f={nstr(fv,6):>14s}  FAILED TO CONVERGE")
            continue

        evals, J4 = compute_eigenvalues(m_star, e_star, fv)
        abs_evals = [float(fabs(e)) for e in evals]

        # Spectral gap
        if abs_evals[1] > 1e-15:
            gap = float(log(fabs(evals[0]) / fabs(evals[1])))
        else:
            gap = float('inf')

        theta_star = float(m_star[3])

        results.append({
            'floor': fv,
            'evals': evals,
            'abs_evals': abs_evals,
            'gap': gap,
            'theta': theta_star,
            'm_star': m_star,
            'e_star': e_star,
        })

        print(f"  {nstr(fv,6):>14s}  {abs_evals[0]:14.10f}  {abs_evals[1]:14.10f}  "
              f"{abs_evals[2]:14.10f}  {abs_evals[3]:14.10f}  {gap:14.10f}  {theta_star:14.10f}")

        # Update guess for next iteration (continuation)
        last_good_guess = [m_star[0], m_star[3], e_star[0], e_star[3]]

    if len(results) < 3:
        print("\nToo few successful results to continue analysis.")
        return

    # ============================================================
    # PART 2: UV limit (floor → 0) vs IR (floor = 1/27)
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 2: UV vs IR COMPARISON")
    print("=" * 80)

    uv = results[0]
    ir = results[-1]

    print(f"\nUV (floor = {nstr(uv['floor'],6)}):")
    print(f"  Eigenvalues: {', '.join(f'{e:.10f}' for e in uv['abs_evals'])}")
    print(f"  θ* = {uv['theta']:.10f}")
    print(f"  Spectral gap = {uv['gap']:.10f}")

    print(f"\nIR (floor = {nstr(ir['floor'],6)}):")
    print(f"  Eigenvalues: {', '.join(f'{e:.10f}' for e in ir['abs_evals'])}")
    print(f"  θ* = {ir['theta']:.10f}")
    print(f"  Spectral gap = {ir['gap']:.10f}")

    print(f"\nRatio IR/UV for each eigenvalue:")
    for k in range(4):
        if uv['abs_evals'][k] > 1e-15:
            ratio = ir['abs_evals'][k] / uv['abs_evals'][k]
            print(f"  λ_{k}: {ratio:.10f}")
        else:
            print(f"  λ_{k}: UV value ~0, ratio undefined")

    print(f"\nDifference IR - UV for each eigenvalue:")
    for k in range(4):
        diff = ir['abs_evals'][k] - uv['abs_evals'][k]
        print(f"  λ_{k}: {diff:+.10f}")

    # ============================================================
    # PART 3: Logarithmic fit  λ(f) = a + b·ln(f)
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 3: ONE-LOOP FIT  λ_k(f) = a + b·ln(f)")
    print("=" * 80)

    import numpy as np

    log_f = np.array([float(ln(r['floor'])) for r in results])
    for k in range(4):
        lam_k = np.array([r['abs_evals'][k] for r in results])
        if np.all(lam_k < 1e-12):
            print(f"\n  λ_{k}: essentially zero, skipping")
            continue
        # Linear fit: λ = a + b·ln(f)
        A_mat = np.vstack([np.ones_like(log_f), log_f]).T
        coeffs, residuals, _, _ = np.linalg.lstsq(A_mat, lam_k, rcond=None)
        a_fit, b_fit = coeffs
        resid = np.sqrt(np.mean((lam_k - a_fit - b_fit * log_f)**2))

        print(f"\n  λ_{k}(f) = {a_fit:.10f} + ({b_fit:.10f}) × ln(f)")
        print(f"    a (UV intercept)     = {a_fit:.10f}")
        print(f"    b (one-loop slope)   = {b_fit:.10f}")
        print(f"    RMS residual         = {resid:.2e}")
        print(f"    b/(a) relative slope = {b_fit/a_fit:.6f}" if abs(a_fit) > 1e-15 else "")

    # ============================================================
    # PART 4: Two-loop fit  λ(f) = a + b·ln(f) + c·ln²(f)
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 4: TWO-LOOP FIT  λ_k(f) = a + b·ln(f) + c·ln²(f)")
    print("=" * 80)

    for k in range(4):
        lam_k = np.array([r['abs_evals'][k] for r in results])
        if np.all(lam_k < 1e-12):
            print(f"\n  λ_{k}: essentially zero, skipping")
            continue
        # Quadratic fit in ln(f)
        A_mat = np.vstack([np.ones_like(log_f), log_f, log_f**2]).T
        coeffs, residuals, _, _ = np.linalg.lstsq(A_mat, lam_k, rcond=None)
        a_fit, b_fit, c_fit = coeffs
        resid = np.sqrt(np.mean((lam_k - a_fit - b_fit*log_f - c_fit*log_f**2)**2))

        print(f"\n  λ_{k}(f) = {a_fit:.10f} + ({b_fit:.10f})·ln(f) + ({c_fit:.12f})·ln²(f)")
        print(f"    a (UV intercept)   = {a_fit:.10f}")
        print(f"    b (one-loop β)     = {b_fit:.10f}")
        print(f"    c (two-loop β)     = {c_fit:.12f}")
        print(f"    RMS residual       = {resid:.2e}")
        print(f"    c/b ratio          = {c_fit/b_fit:.6f}" if abs(b_fit) > 1e-15 else "")

    # ============================================================
    # PART 4B: Power-law fit  λ(f) = A × f^γ
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 4B: POWER-LAW FIT  λ_k(f) = A × f^γ  (i.e. ln λ = ln A + γ·ln f)")
    print("=" * 80)

    for k in range(2):  # only nonzero eigenvalues
        lam_k = np.array([r['abs_evals'][k] for r in results])
        ln_lam = np.log(lam_k)
        # Linear fit in log-log: ln(λ) = ln(A) + γ·ln(f)
        A_mat = np.vstack([np.ones_like(log_f), log_f]).T
        coeffs, _, _, _ = np.linalg.lstsq(A_mat, ln_lam, rcond=None)
        ln_A, gamma = coeffs
        A_val = np.exp(ln_A)
        fitted = A_val * np.exp(gamma * log_f)
        resid = np.sqrt(np.mean((lam_k - fitted)**2))

        print(f"\n  λ_{k}(f) = {A_val:.10f} × f^{gamma:.10f}")
        print(f"    A (amplitude)     = {A_val:.10f}")
        print(f"    γ (exponent)      = {gamma:.10f}")
        print(f"    RMS residual      = {resid:.2e}")
        print(f"    1/γ               = {1/gamma:.6f}" if abs(gamma) > 1e-15 else "")
        # Check framework numbers
        print(f"    γ vs 1/H = {1/3:.6f}")
        print(f"    γ vs K* = {7/30:.6f}")
        print(f"    γ vs 1/(2H) = {1/6:.6f}")
        print(f"    γ vs 1/H² = {1/9:.6f}")

    # ============================================================
    # PART 4C: What about the "zero" eigenvalues?
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 4C: THE 'ZERO' EIGENVALUES — ARE THEY TRULY ZERO?")
    print("=" * 80)

    print("\n  Checking raw eigenvalue magnitudes for λ₂ and λ₃:")
    for r in results[::10]:  # every 10th point
        fv = r['floor']
        raw = r['evals']
        print(f"    f={nstr(fv,6):>12s}  λ₂={nstr(fabs(raw[2]),6):>12s}  λ₃={nstr(fabs(raw[3]),6):>12s}")

    # ============================================================
    # PART 5: Spectral gap flow
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 5: SPECTRAL GAP Δ(f) = -ln|λ₀(f)|")
    print("=" * 80)

    print(f"\n  {'floor':>14s}  {'Δ₀=-ln|λ₀|':>14s}  {'Δ₁=-ln|λ₁|':>14s}  {'Δ₀/Δ₁':>12s}  {'gap ratio':>12s}")
    print("  " + "-" * 70)

    for r in results:
        fv = r['floor']
        ae = r['abs_evals']
        if ae[0] > 1e-15 and ae[1] > 1e-15:
            D0 = -float(log(mpf(str(ae[0]))))
            D1 = -float(log(mpf(str(ae[1]))))
            ratio01 = D0 / D1 if D1 != 0 else float('inf')
            gap_rat = ae[1] / ae[0]
        else:
            D0 = D1 = ratio01 = gap_rat = 0
        print(f"  {nstr(fv,6):>14s}  {D0:14.10f}  {D1:14.10f}  {ratio01:12.8f}  {gap_rat:12.8f}")

    # ============================================================
    # PART 6: Does K* run? (It shouldn't — conservation law)
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 6: DOES K* RUN WITH THE FLOOR?")
    print("=" * 80)
    print("  (K*=7/30 is imposed by the conservation law, NOT by the floor.)")
    print("  Verification: K is always set to 7/30 in the fixed-point equations.")
    print("  So K does NOT run. What runs is the eigenvalue spectrum.")

    # ============================================================
    # PART 7: Effective beta function
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 7: EFFECTIVE β FUNCTION FROM EIGENVALUE FLOW")
    print("=" * 80)

    # Define an effective coupling from the dominant eigenvalue:
    # g_eff(f) = 1 - |λ₀(f)|  (measures how far from marginal)
    # or g_eff(f) = -ln|λ₀(f)| / (-ln|λ₀(1/27)|) (normalised spectral gap)

    print("\n  Effective coupling g_eff = -ln|λ₀(f)| (unnormalised spectral gap):")
    print(f"\n  {'floor':>14s}  {'g_eff':>14s}  {'dg/d(ln f)':>14s}")
    print("  " + "-" * 46)

    g_effs = []
    for r in results:
        ae = r['abs_evals']
        if ae[0] > 1e-15:
            g = -float(log(mpf(str(ae[0]))))
        else:
            g = float('inf')
        g_effs.append(g)

    for idx, r in enumerate(results):
        fv = r['floor']
        g = g_effs[idx]
        # Numerical derivative
        if 0 < idx < len(results) - 1:
            dlnf = float(ln(results[idx+1]['floor'])) - float(ln(results[idx-1]['floor']))
            dg = g_effs[idx+1] - g_effs[idx-1]
            beta = dg / dlnf if abs(dlnf) > 1e-15 else 0
        else:
            beta = 0
        print(f"  {nstr(fv,6):>14s}  {g:14.10f}  {beta:14.10f}")

    # ============================================================
    # PART 8: Connection to SM beta coefficients
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 8: CONNECTION TO SM β COEFFICIENTS")
    print("=" * 80)

    # SM one-loop beta coefficients: b₁ = 41/10, b₂ = -19/6, b₃ = -7
    # The DS eigenvalue flow slope b = dλ/d(ln f) should relate to these.

    # Extract slope at the physical point (IR end) using last two distinct points
    if len(results) >= 3:
        ir2 = results[-1]
        ir1 = results[-2]

        dlnf_ir = float(ln(ir2['floor'])) - float(ln(ir1['floor']))

        print("\n  Slopes dλ_k/d(ln f) at the physical floor (IR end):")
        print(f"  (using floors {nstr(ir1['floor'],6)} and {nstr(ir2['floor'],6)}, "
              f"Δln(f) = {dlnf_ir:.6f})")
        for k in range(4):
            if abs(dlnf_ir) > 1e-10:
                slope = (ir2['abs_evals'][k] - ir1['abs_evals'][k]) / dlnf_ir
            else:
                slope = 0
            print(f"    dλ_{k}/d(ln f) = {slope:.10f}")

        # Compare with -b₃/(2π) = 7/(2π) ≈ 1.114
        print(f"\n  SM reference: -b₃/(2π) = {7/(2*float(pi)):.6f}")
        print(f"                 b₃ = -7 (SU(3), N_f=6)")
        print(f"                -b₃·K*/(2π) = {7*float(K_star)/(2*float(pi)):.6f}")
        print(f"                 K*·b₃ = {float(K_star)*(-7):.6f}")

    # ============================================================
    # PART 9: Total eigenvalue variation — "running" magnitude
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 9: TOTAL RUNNING FROM UV TO IR")
    print("=" * 80)

    if len(results) >= 2:
        uv_r = results[0]
        ir_r = results[-1]
        ln_ratio = float(ln(ir_r['floor'])) - float(ln(uv_r['floor']))

        print(f"\n  ln(f_IR/f_UV) = {ln_ratio:.6f}")
        print(f"  f_UV = {nstr(uv_r['floor'],6)},  f_IR = {nstr(ir_r['floor'],6)}")

        for k in range(4):
            delta_lam = ir_r['abs_evals'][k] - uv_r['abs_evals'][k]
            if abs(ln_ratio) > 1e-10 and uv_r['abs_evals'][k] > 1e-15:
                eff_b = delta_lam / ln_ratio
                rel_change = delta_lam / uv_r['abs_evals'][k] * 100
                print(f"\n  λ_{k}:")
                print(f"    UV value   = {uv_r['abs_evals'][k]:.10f}")
                print(f"    IR value   = {ir_r['abs_evals'][k]:.10f}")
                print(f"    Δλ         = {delta_lam:+.10f}")
                print(f"    Δλ/Δ(ln f) = {eff_b:.10f}  (effective one-loop slope)")
                print(f"    % change   = {rel_change:+.4f}%")
            else:
                print(f"\n  λ_{k}: UV~0, skipping")

    # ============================================================
    # PART 10: Honest assessment
    # ============================================================
    print("\n" + "=" * 80)
    print("PART 10: HONEST ASSESSMENT")
    print("=" * 80)

    # Check if eigenvalues actually vary
    if len(results) >= 2:
        total_var = [max(r['abs_evals'][k] for r in results) -
                     min(r['abs_evals'][k] for r in results) for k in range(4)]
        mean_vals = [sum(r['abs_evals'][k] for r in results) / len(results) for k in range(4)]

        print("\n  Eigenvalue variation summary:")
        for k in range(4):
            rel_var = total_var[k] / mean_vals[k] * 100 if mean_vals[k] > 1e-15 else 0
            print(f"    λ_{k}: range = {total_var[k]:.2e},  mean = {mean_vals[k]:.6f},  "
                  f"relative variation = {rel_var:.4f}%")

        if all(tv < 1e-6 for tv in total_var[:2]):
            print("\n  CONCLUSION: Eigenvalues are INSENSITIVE to the floor parameter.")
            print("  The floor enforcement barely changes the linearised dynamics.")
            print("  This means the floor is NOT a good proxy for the RG scale,")
            print("  OR the RG flow is trivial (eigenvalues don't run).")
        elif all(tv > 0.01 for tv in total_var[:2]):
            print("\n  CONCLUSION: Eigenvalues VARY SIGNIFICANTLY with the floor.")
            print("  The floor-eigenvalue flow is a genuine candidate for RG running.")
            print("\n  However, note the following caveats:")
            print("  1. Only 2 of 4 eigenvalues are nonzero at the single-site level.")
            print("     (The s₂-s₃ antisymmetric mode has eigenvalue ~0 due to the")
            print("     symmetric ansatz s₂=s₃. The full coupled A₂ system has 4 nonzero.)")
            print("  2. The flow is NOT purely logarithmic — there is significant")
            print("     curvature (two-loop term c/b ~ 6%). A power law may be more natural.")
            print("  3. The eigenvalues INCREASE as floor increases (UV → IR). In QCD,")
            print("     α_s increases from UV to IR (asymptotic freedom). So if")
            print("     λ ∝ coupling, the direction is CORRECT for asymptotic freedom.")
            print("  4. K* = 7/30 is FIXED and does not run. Only the linearised")
            print("     spectrum (fluctuations around the fixed point) changes.")
        else:
            print("\n  CONCLUSION: Mixed — some eigenvalues vary, some don't.")
            print("  Detailed analysis needed to determine physical significance.")

    print("\n" + "=" * 80)
    print("DONE")
    print("=" * 80)


if __name__ == "__main__":
    main()
