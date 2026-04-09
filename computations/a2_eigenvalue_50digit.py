#!/usr/bin/env python3
"""
A₂ coupled Jacobian eigenvalues to 50+ digit precision.
Tests whether the 0.11% gap is numerical artifact or real physics.

Ported from exotic_analysis/glueball_A2_spectrum.py (the correct dynamics)
into mpmath arbitrary-precision arithmetic.
"""

from mpmath import mp, mpf, matrix, sqrt, log, fabs, nstr, eig, mpf as M
from mpmath import findroot

# 80 internal digits → 50+ reliable output digits
mp.dps = 80

H = 3
FLOOR = mpf(1) / mpf(H**3)  # 1/27

# ============================================================
# DS COMBINATION AND BORN FLOOR (mpmath)
# ============================================================

def ds_combine(m, e):
    """DS combination: m=(s1,s2,s3,theta), e=(e1,e2,e3,phi)"""
    s = m[:3]
    theta = m[3]
    ev = e[:3]
    phi = e[3]

    # Pre-normalization output
    s_pre = [s[i]*ev[i] + s[i]*phi + theta*ev[i] for i in range(3)]
    theta_pre = theta * phi

    # Conflict
    total_pre = sum(s_pre) + theta_pre
    K = mpf(1) - total_pre

    if fabs(mpf(1) - K) < mpf(10)**(-70):
        return list(m), K  # degenerate

    denom = mpf(1) - K
    out = [sp / denom for sp in s_pre] + [theta_pre / denom]
    return out, K


def born_prob(m):
    """Born probability of ignorance: θ²/Σmᵢ²"""
    L2sq = sum(x**2 for x in m)
    if L2sq < mpf(10)**(-60):
        return mpf(0)
    return m[3]**2 / L2sq


def enforce_floor(m, floor_val=None):
    """Enforce Born floor analytically.
    Solve: 26·t²·S² = (1-t)²·Sq  for t (new theta).
    """
    if floor_val is None:
        floor_val = FLOOR

    b = born_prob(m)
    if b >= floor_val - mpf(10)**(-60):
        return list(m)

    S = sum(m[:3])
    if S < mpf(10)**(-60):
        return list(m)

    Sq = sum(x**2 for x in m[:3])

    # 26·t²·S² = (1-t)²·Sq
    # t²·(26·S² - Sq) + 2·t·Sq - Sq = 0
    A_coeff = mpf(26) * S**2 - Sq
    B_coeff = mpf(2) * Sq
    C_coeff = -Sq

    disc = B_coeff**2 - 4*A_coeff*C_coeff
    if disc < 0:
        return list(m)

    t1 = (-B_coeff + sqrt(disc)) / (2*A_coeff)
    t2 = (-B_coeff - sqrt(disc)) / (2*A_coeff)

    # Choose t in (0,1) closest to current theta
    candidates = [t for t in [t1, t2] if mpf(0) < t < mpf(1)]
    if not candidates:
        return list(m)

    t = min(candidates, key=lambda x: fabs(x - m[3]))

    # Rescale singletons
    alpha = (mpf(1) - t) / S
    out = [m[i]*alpha for i in range(3)] + [t]
    return out


def ds_step(m, e):
    """One full DS + floor step"""
    m_ds, K = ds_combine(m, e)
    m_floor = enforce_floor(m_ds)
    return m_floor, K


# ============================================================
# SINGLE-SITE EQUILIBRIUM: jointly solve (m*, e*)
# ============================================================

def find_single_site_equilibrium():
    """Find (m*, e*) by solving 4 equations in 4 unknowns:
    1. Born(m*) = 1/27
    2. Born(e*) = 1/27
    3. K(m*, e*) = 7/30
    4. m* = Φ(m*, e*)  [fixed point]

    Parameterized by (s1, theta, w1, phi) with s2=s3, w2=w3 by symmetry.
    """
    import numpy as np
    from scipy.optimize import fsolve as scipy_fsolve

    # First get a rough solution with scipy
    def equations_rough(params):
        s1, theta, w1, phi = params
        s2 = (1.0 - s1 - theta) / 2.0
        w2 = (1.0 - w1 - phi) / 2.0

        m = [s1, s2, s2, theta]
        e = [w1, w2, w2, phi]

        L2m = s1**2 + 2*s2**2 + theta**2
        eq1 = theta**2 / L2m - 1.0/27.0

        L2e = w1**2 + 2*w2**2 + phi**2
        eq2 = phi**2 / L2e - 1.0/27.0

        K = s1*w2 + s1*w2 + s2*w1 + s2*w2 + s2*w1 + s2*w2
        eq3 = K - 7.0/30.0

        # Use mpf for the ds_step (convert back)
        m_mp = [mpf(x) for x in m]
        e_mp = [mpf(x) for x in e]
        m_out, _ = ds_step(m_mp, e_mp)
        eq4 = float(m_out[0]) - s1

        return [eq1, eq2, eq3, eq4]

    x0 = [0.787, 0.155, 0.631, 0.129]
    sol_rough = scipy_fsolve(equations_rough, x0)
    print(f"  Rough solution (scipy): s1={sol_rough[0]:.8f}, θ={sol_rough[1]:.8f}, w1={sol_rough[2]:.8f}, φ={sol_rough[3]:.8f}")

    # Now refine with mpmath findroot
    def equations_mp(*params):
        s1, theta, w1, phi = params
        s2 = (mpf(1) - s1 - theta) / 2
        w2 = (mpf(1) - w1 - phi) / 2

        m = [s1, s2, s2, theta]
        e = [w1, w2, w2, phi]

        # Eq 1: Born(m) = 1/27
        L2m = s1**2 + 2*s2**2 + theta**2
        eq1 = theta**2 / L2m - FLOOR

        # Eq 2: Born(e) = 1/27
        L2e = w1**2 + 2*w2**2 + phi**2
        eq2 = phi**2 / L2e - FLOOR

        # Eq 3: K(m,e) = 7/30
        K = s1*w2 + s1*w2 + s2*w1 + s2*w2 + s2*w1 + s2*w2
        eq3 = K - mpf(7)/30

        # Eq 4: fixed point
        m_out, _ = ds_step(m, e)
        eq4 = m_out[0] - s1

        return [eq1, eq2, eq3, eq4]

    x_star = findroot(equations_mp, [mpf(str(x)) for x in sol_rough])

    s1, theta, w1, phi = x_star
    s2 = (mpf(1) - s1 - theta) / 2
    w2 = (mpf(1) - w1 - phi) / 2

    m_star = [s1, s2, s2, theta]
    e_star = [w1, w2, w2, phi]

    return m_star, e_star


# ============================================================
# COUPLED A₂ SYSTEM
# ============================================================

def coupled_evidence(m1, m2, e0, g):
    """Evidence for node 1 in A₂ system. A_{12} = -1."""
    uniform = [mpf(1)/4]*4
    e1 = [e0[i] - g * (m2[i] - uniform[i]) for i in range(4)]
    # Ensure positivity
    e1 = [max(x, mpf(10)**(-50)) for x in e1]
    # Renormalize
    total = sum(e1)
    e1 = [x / total for x in e1]
    return e1


def coupled_A2_step(state, e0, g):
    """One step of coupled A₂ system. state = m1 ++ m2 (8 elements)."""
    m1 = state[:4]
    m2 = state[4:]

    e1 = coupled_evidence(m1, m2, e0, g)
    e2 = coupled_evidence(m2, m1, e0, g)

    m1_new, K1 = ds_step(m1, e1)
    m2_new, K2 = ds_step(m2, e2)

    return m1_new + m2_new


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("A₂ COUPLED JACOBIAN: 50-DIGIT EIGENVALUE COMPUTATION")
    print("=" * 70)
    print(f"Working precision: {mp.dps} decimal digits")
    print()

    # Step 1: Find single-site equilibrium
    print("Step 1: Finding jointly-determined (m*, e*) equilibrium...")
    m_star, e_star = find_single_site_equilibrium()

    print(f"\n  m* = ({nstr(m_star[0],40)}, {nstr(m_star[1],40)}, {nstr(m_star[2],40)}, {nstr(m_star[3],40)})")
    print(f"  e* = ({nstr(e_star[0],40)}, {nstr(e_star[1],40)}, {nstr(e_star[2],40)}, {nstr(e_star[3],40)})")

    # Verify constraints
    born_m = born_prob(m_star)
    born_e = born_prob(e_star)
    K_check = (m_star[0]*e_star[1] + m_star[0]*e_star[2] +
               m_star[1]*e_star[0] + m_star[1]*e_star[2] +
               m_star[2]*e_star[0] + m_star[2]*e_star[1])
    m_check, _ = ds_step(m_star, e_star)
    fp_err = max(fabs(m_check[i] - m_star[i]) for i in range(4))

    print(f"\n  Born(m*) = {nstr(born_m, 40)}  (target: {nstr(FLOOR, 20)})")
    print(f"  Born(e*) = {nstr(born_e, 40)}")
    print(f"  K(m*,e*) = {nstr(K_check, 40)}  (target: {nstr(mpf(7)/30, 20)})")
    print(f"  ||Φ(m*)-m*|| = {nstr(fp_err, 10)}")
    print()

    # Step 2: Find coupled A₂ fixed point
    print("Step 2: Finding coupled A₂ fixed point...")
    g = mpf(7) / mpf(30)

    # Start at single-site equilibrium (both sites)
    state = list(m_star) + list(m_star)

    for iteration in range(10000):
        state_new = coupled_A2_step(state, e_star, g)
        diff = max(fabs(state_new[i] - state[i]) for i in range(8))
        if iteration < 5 or iteration % 1000 == 0:
            print(f"  Step {iteration}: diff = {nstr(diff, 5)}")
        if diff < mpf(10)**(-70):
            print(f"  Converged at step {iteration}")
            break
        state = state_new

    m1_eq = state_new[:4]
    m2_eq = state_new[4:]
    print(f"\n  m1* = ({nstr(m1_eq[0],30)}, {nstr(m1_eq[1],30)}, {nstr(m1_eq[2],30)}, {nstr(m1_eq[3],30)})")
    print(f"  m2* = ({nstr(m2_eq[0],30)}, {nstr(m2_eq[1],30)}, {nstr(m2_eq[2],30)}, {nstr(m2_eq[3],30)})")
    print()

    # Step 3: Compute 8×8 Jacobian
    print("Step 3: Computing 8×8 Jacobian at high precision...")

    # Use central differences for better accuracy
    eps = mpf(10)**(-35)  # perturbation (sweet spot for 80-digit precision)

    def coupled_map(s):
        return coupled_A2_step(s, e_star, g)

    J = matrix(8, 8)
    for j in range(8):
        s_plus = list(state_new)
        s_minus = list(state_new)
        s_plus[j] = s_plus[j] + eps
        s_minus[j] = s_minus[j] - eps
        f_plus = coupled_map(s_plus)
        f_minus = coupled_map(s_minus)
        for i in range(8):
            J[i, j] = (f_plus[i] - f_minus[i]) / (2*eps)

    print("  Jacobian computed (central differences).")
    print()

    # Step 4: Eigenvalues with mpmath
    print("Step 4: Computing eigenvalues with mpmath...")
    eigenvalues, eigenvectors = eig(J)

    # Sort by magnitude (descending)
    evals_data = []
    for k, ev in enumerate(eigenvalues):
        if hasattr(ev, 'imag') and fabs(ev.imag) > mpf(10)**(-40):
            mag = sqrt(ev.real**2 + ev.imag**2)
        else:
            mag = fabs(ev.real if hasattr(ev, 'real') else ev)
        evals_data.append((k, ev, mag))

    evals_data.sort(key=lambda x: -float(x[2]))

    print("\n  All 8 eigenvalues (sorted by |λ|):")
    print("  " + "-"*65)
    for rank, (k, ev, mag) in enumerate(evals_data):
        if hasattr(ev, 'imag') and fabs(ev.imag) > mpf(10)**(-40):
            print(f"    λ_{rank}: {nstr(ev.real, 30)} + {nstr(ev.imag, 30)}i")
            print(f"        |λ| = {nstr(mag, 50)}")
        else:
            val = ev.real if hasattr(ev, 'real') else ev
            print(f"    λ_{rank}: {nstr(val, 50)}")
            print(f"        |λ| = {nstr(mag, 50)}")
    print()

    # Step 5: Extract significant eigenvalues and compare
    sig = [(k, ev, mag) for k, ev, mag in evals_data if float(mag) > 1e-6]

    print("=" * 70)
    print("50-DIGIT EIGENVALUES vs 4-DIGIT ORIGINALS")
    print("=" * 70)
    print()

    expected = [
        (0.5022, "0++  ground state"),
        (0.4745, "1.082 angular mode"),
        (0.3527, "0-+  symmetric angular"),
        (0.3344, "0++* symmetric radial"),
    ]

    for i in range(min(4, len(sig))):
        _, ev_i, mag_i = sig[i]
        orig, label = expected[i]
        val = float(mag_i)
        diff = abs(val - orig)
        pct = 100 * diff / orig if orig > 0 else 0
        print(f"  λ_{i} ({label}):")
        print(f"    4-digit:  {orig}")
        print(f"    50-digit: {nstr(mag_i, 50)}")
        print(f"    Δ = {diff:.2e}  ({pct:.6f}%)")
        print()

    # Step 6: Mass gaps and bridge relations
    print("=" * 70)
    print("MASS GAPS AND BRIDGE RELATIONS")
    print("=" * 70)
    print()

    if len(sig) >= 4:
        lam = [sig[i][2] for i in range(4)]
        Delta = [-log(l) for l in lam]
        D0 = Delta[0]

        print("  Mass gaps Δ = -ln(λ):")
        for i in range(4):
            print(f"    Δ_{i} = {nstr(Delta[i], 50)}")
        print()

        print("  Mass ratios Δᵢ/Δ₀:")
        for i in range(4):
            ratio = Delta[i] / D0
            print(f"    Δ_{i}/Δ₀ = {nstr(ratio, 50)}")
        print()

        # Q bridge: λ₃/λ₀ vs 2/3
        Q = lam[3] / lam[0]
        Q_target = mpf(2) / mpf(3)
        Q_gap_pct = float(fabs(Q - Q_target) / Q_target * 100)

        print(f"  Q bridge (λ₃/λ₀):")
        print(f"    computed = {nstr(Q, 50)}")
        print(f"    target   = {nstr(Q_target, 50)}")
        print(f"    gap      = {Q_gap_pct:.10f}%")
        print()

        # Delta ratio: Δ₂/Δ₃ vs 20/21
        DR = Delta[2] / Delta[3]
        DR_target = mpf(20) / mpf(21)
        DR_gap_pct = float(fabs(DR - DR_target) / DR_target * 100)

        print(f"  Δ ratio (Δ₂/Δ₃):")
        print(f"    computed = {nstr(DR, 50)}")
        print(f"    target   = {nstr(DR_target, 50)}")
        print(f"    gap      = {DR_gap_pct:.10f}%")
        print()

        # λ₁/λ₀ ratio
        R10 = lam[1] / lam[0]
        print(f"  λ₁/λ₀ ratio:")
        print(f"    computed = {nstr(R10, 50)}")
        print()

        # λ₂/λ₃ ratio
        R23 = lam[2] / lam[3]
        print(f"  λ₂/λ₃ ratio:")
        print(f"    computed = {nstr(R23, 50)}")
        print()

        # Verdict
        print("=" * 70)
        if Q_gap_pct < 0.01 and DR_gap_pct < 0.01:
            print("VERDICT: GAPS CLOSED — discrepancy was numerical artifact.")
        elif Q_gap_pct < 0.05 or DR_gap_pct < 0.05:
            print("VERDICT: GAPS SIGNIFICANTLY REDUCED — mostly numerical.")
        elif abs(Q_gap_pct - 0.11) < 0.05:
            print("VERDICT: GAPS UNCHANGED at ~0.11% — this is REAL PHYSICS.")
            print("  The gap persists at 50 digits: it's not a precision artifact.")
        else:
            print("VERDICT: GAPS PERSIST — the discrepancy is real, not precision.")
        print(f"  Q gap:  {Q_gap_pct:.10f}%")
        print(f"  Δ gap:  {DR_gap_pct:.10f}%")
        print("=" * 70)


if __name__ == "__main__":
    main()
