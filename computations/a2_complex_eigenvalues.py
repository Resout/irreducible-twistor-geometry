#!/usr/bin/env python3
"""
A₂ coupled Jacobian eigenvalues on C⁴ (not R⁴).

The real Jacobian DΦ|_R mixes holomorphic and anti-holomorphic parts.
This computes the Wirtinger decomposition:
  ∂Φ/∂z  = (1/2)(∂Φ/∂x - i·∂Φ/∂y)   [holomorphic]
  ∂Φ/∂z̄  = (1/2)(∂Φ/∂x + i·∂Φ/∂y)   [anti-holomorphic]

and extracts eigenvalues of each separately.

Hypothesis: the holomorphic eigenvalues give Q = 2/3 exactly,
while the real Jacobian eigenvalues gave Q = 0.66596 (0.1% off).
"""

import numpy as np
from scipy.optimize import fsolve, brentq

H = 3
FLOOR = 1.0 / H**3

# ============================================================
# COMPLEX DS COMBINATION AND BORN FLOOR
# ============================================================

def ds_combine_complex(m, e):
    """DS combination for complex mass functions m, e in C^4."""
    s = m[:3]
    theta = m[3]
    se = e[:3]
    phi = e[3]

    s_new = s*se + s*phi + theta*se
    theta_new = theta * phi

    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i != j)

    total_pre = np.sum(s_new) + theta_new
    # K = 1 - total_pre (pre-normalization)

    denom = 1.0 - K
    if abs(denom) < 1e-15:
        return m.copy(), K

    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom

    # L1 normalize (complex: sum of components = 1)
    m_out = m_out / np.sum(m_out)

    return m_out, K


def born_prob_complex(m):
    """Born probability for complex mass function: |θ|²/Σ|m_i|²"""
    mod_sq = sum(abs(m[i])**2 for i in range(4))
    if mod_sq < 1e-30:
        return 0.0
    return abs(m[3])**2 / mod_sq


def enforce_floor_complex(m, floor_val=FLOOR):
    """Enforce Born floor on complex mass function.

    Preserves phase of θ, preserves ratios of s_i.
    Adjusts |θ| upward and rescales s_i to maintain Σm_i = 1.
    """
    b = born_prob_complex(m)
    if b >= floor_val - 1e-14:
        return m.copy()

    S = np.sum(m[:3])  # complex sum
    if abs(S) < 1e-15:
        return m.copy()

    theta = m[3]
    phase_theta = theta / abs(theta) if abs(theta) > 1e-30 else 1.0

    # s_i ratios preserved: s_i_new = s_i * α where α = (1 - θ_new)/S
    # Born constraint: |θ_new|² / (|α|²·Sq_mod + |θ_new|²) = 1/27
    # where Sq_mod = Σ|s_i|², θ_new = t·phase_θ, α = (1 - t·phase_θ)/S

    Sq_mod = sum(abs(m[i])**2 for i in range(3))

    def born_at_t(t):
        """Born probability when |θ_new| = t."""
        theta_new = t * phase_theta
        alpha = (1.0 - theta_new) / S
        sq_scaled = abs(alpha)**2 * Sq_mod
        return t**2 / (sq_scaled + t**2) - floor_val

    # Bisection: find t where Born = 1/27
    t_lo = abs(theta)  # current |θ| (Born too low)
    t_hi = 1.0

    # Make sure we bracket the root
    if born_at_t(t_hi) < 0:
        return m.copy()  # can't reach floor

    try:
        t_star = brentq(born_at_t, t_lo, t_hi, xtol=1e-15)
    except ValueError:
        # fallback: bisection manually
        for _ in range(100):
            t_mid = (t_lo + t_hi) / 2
            if born_at_t(t_mid) < 0:
                t_lo = t_mid
            else:
                t_hi = t_mid
        t_star = (t_lo + t_hi) / 2

    theta_new = t_star * phase_theta
    alpha = (1.0 - theta_new) / S

    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = m[:3] * alpha
    m_out[3] = theta_new

    return m_out


def ds_step_complex(m, e):
    """One full DS + floor step for complex mass functions."""
    m_ds, K = ds_combine_complex(m, e)
    m_floor = enforce_floor_complex(m_ds)
    return m_floor, K


# ============================================================
# SINGLE-SITE EQUILIBRIUM (real, from glueball_A2_spectrum.py)
# ============================================================

def find_equilibrium():
    """Find jointly-determined (m*, e*) with Born=1/27, K=7/30."""
    def enforce_floor_real(m):
        s, theta = m[:3], m[3]
        ssq = sum(si**2 for si in s)
        if theta**2 / (ssq + theta**2) >= FLOOR:
            return m.copy()
        ss = sum(s)
        if ss < 1e-15: return m.copy()
        r = ssq / ss**2
        disc = (2*r)**2 + 4*(26-r)*r
        t = (-2*r + np.sqrt(disc)) / (2*(26-r))
        alpha = (1-t) / ss
        return np.array([s[0]*alpha, s[1]*alpha, s[2]*alpha, t])

    def ds_step_real(m, e):
        s, theta = m[:3], m[3]
        se, phi = e[:3], e[3]
        s_new = s*se + s*phi + theta*se
        theta_new = theta*phi
        K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i != j)
        d = 1.0 - K
        m_out = np.zeros(4)
        m_out[:3] = s_new/d
        m_out[3] = theta_new/d
        m_out = m_out / np.sum(m_out)
        return enforce_floor_real(m_out), K

    def equations(params):
        s1, theta, w1, phi = params
        s2 = (1.0 - s1 - theta) / 2.0
        w2 = (1.0 - w1 - phi) / 2.0
        m = np.array([s1, s2, s2, theta])
        e = np.array([w1, w2, w2, phi])

        L2m = s1**2 + 2*s2**2 + theta**2
        eq1 = theta**2 / L2m - FLOOR

        L2e = w1**2 + 2*w2**2 + phi**2
        eq2 = phi**2 / L2e - FLOOR

        K = s1*w2 + s1*w2 + s2*w1 + s2*w2 + s2*w1 + s2*w2
        eq3 = K - 7.0/30

        m_out, _ = ds_step_real(m, e)
        eq4 = m_out[0] - s1

        return [eq1, eq2, eq3, eq4]

    sol = fsolve(equations, [0.787, 0.155, 0.631, 0.129])
    s1, theta, w1, phi = sol
    s2 = (1.0 - s1 - theta) / 2.0
    w2 = (1.0 - w1 - phi) / 2.0
    return np.array([s1, s2, s2, theta]), np.array([w1, w2, w2, phi])


# ============================================================
# A₂ COUPLED SYSTEM
# ============================================================

def coupled_evidence(m1, m2, e0, g):
    """Evidence for node 1 in A₂. Works for complex m."""
    uniform = np.array([0.25, 0.25, 0.25, 0.25], dtype=complex)
    e1 = e0.astype(complex) - g * (m2 - uniform)
    # Ensure real parts positive (clip near zero)
    for i in range(4):
        if e1[i].real < 1e-10:
            e1[i] = 1e-10 + e1[i].imag * 1j
    e1 = e1 / np.sum(e1)
    return e1


def coupled_A2_step(state, e0, g, use_complex=False):
    """One step of coupled A₂. state = m1 ++ m2 (8-vector)."""
    m1 = state[:4]
    m2 = state[4:]

    e1 = coupled_evidence(m1, m2, e0, g)
    e2 = coupled_evidence(m2, m1, e0, g)

    if use_complex:
        m1_new, _ = ds_step_complex(m1, e1)
        m2_new, _ = ds_step_complex(m2, e2)
    else:
        # Real step
        m1_new, _ = ds_step_complex(np.real(m1), np.real(e1))
        m2_new, _ = ds_step_complex(np.real(m2), np.real(e2))

    return np.concatenate([m1_new, m2_new])


# ============================================================
# WIRTINGER JACOBIANS
# ============================================================

def compute_wirtinger_jacobians(state, e0, g, eps=1e-8):
    """Compute holomorphic and anti-holomorphic Wirtinger Jacobians.

    ∂Φ/∂z_j  = (1/2)(∂Φ/∂x_j - i·∂Φ/∂y_j)
    ∂Φ/∂z̄_j = (1/2)(∂Φ/∂x_j + i·∂Φ/∂y_j)

    where x_j = Re(z_j), y_j = Im(z_j).
    """
    n = len(state)
    J_holo = np.zeros((n, n), dtype=complex)   # ∂Φ/∂z
    J_anti = np.zeros((n, n), dtype=complex)   # ∂Φ/∂z̄

    for j in range(n):
        # Real perturbation: ∂Φ/∂x_j
        s_px = state.astype(complex).copy()
        s_mx = state.astype(complex).copy()
        s_px[j] += eps
        s_mx[j] -= eps
        f_px = coupled_A2_step(s_px, e0, g, use_complex=True)
        f_mx = coupled_A2_step(s_mx, e0, g, use_complex=True)
        dfdx = (f_px - f_mx) / (2*eps)

        # Imaginary perturbation: ∂Φ/∂y_j
        s_py = state.astype(complex).copy()
        s_my = state.astype(complex).copy()
        s_py[j] += 1j*eps
        s_my[j] -= 1j*eps
        f_py = coupled_A2_step(s_py, e0, g, use_complex=True)
        f_my = coupled_A2_step(s_my, e0, g, use_complex=True)
        dfdy = (f_py - f_my) / (2*eps)

        # Wirtinger derivatives
        J_holo[:, j] = 0.5 * (dfdx - 1j * dfdy)
        J_anti[:, j] = 0.5 * (dfdx + 1j * dfdy)

    return J_holo, J_anti


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("A₂ EIGENVALUES ON C⁴: WIRTINGER DECOMPOSITION")
    print("=" * 70)
    print()

    # Step 1: Find equilibrium
    print("Step 1: Finding (m*, e*) equilibrium...")
    m_star, e_star = find_equilibrium()
    print(f"  m* = ({m_star[0]:.8f}, {m_star[1]:.8f}, {m_star[2]:.8f}, {m_star[3]:.8f})")
    print(f"  e* = ({e_star[0]:.8f}, {e_star[1]:.8f}, {e_star[2]:.8f}, {e_star[3]:.8f})")

    K_check = (m_star[0]*e_star[1] + m_star[0]*e_star[2] +
               m_star[1]*e_star[0] + m_star[1]*e_star[2] +
               m_star[2]*e_star[0] + m_star[2]*e_star[1])
    print(f"  K = {K_check:.10f} (target: {7/30:.10f})")
    print()

    # Step 2: Find coupled A₂ fixed point
    print("Step 2: Finding coupled A₂ fixed point...")
    g = 7.0/30
    state = np.concatenate([m_star, m_star])

    for i in range(5000):
        state_new = coupled_A2_step(state, e_star, g, use_complex=True)
        diff = np.max(np.abs(state_new - state))
        if i < 5:
            print(f"  Step {i}: diff = {diff:.2e}")
        if diff < 1e-14:
            print(f"  Converged at step {i}")
            break
        state = np.real(state_new)  # project to real for iteration

    state_eq = np.real(state_new)
    m1_eq = state_eq[:4]
    m2_eq = state_eq[4:]
    print(f"  m1* = ({m1_eq[0]:.8f}, {m1_eq[1]:.8f}, {m1_eq[2]:.8f}, {m1_eq[3]:.8f})")
    print()

    # Step 3: Real Jacobian (what we had before)
    print("Step 3: Real Jacobian eigenvalues (for comparison)...")
    eps_r = 1e-8
    f0 = coupled_A2_step(state_eq, e_star, g, use_complex=True)
    J_real = np.zeros((8, 8))
    for j in range(8):
        s_p = state_eq.copy()
        s_m = state_eq.copy()
        s_p[j] += eps_r
        s_m[j] -= eps_r
        fp = coupled_A2_step(s_p, e_star, g, use_complex=True)
        fm = coupled_A2_step(s_m, e_star, g, use_complex=True)
        J_real[:, j] = np.real((fp - fm) / (2*eps_r))

    evals_real = np.sort(np.abs(np.linalg.eigvals(J_real)))[::-1]
    print("  Real Jacobian |λ|:")
    for i in range(8):
        if evals_real[i] > 1e-10:
            print(f"    |λ_{i}| = {evals_real[i]:.10f}")
        else:
            print(f"    |λ_{i}| ≈ 0")

    Q_real = evals_real[3] / evals_real[0] if evals_real[0] > 0 else 0
    print(f"\n  Q_real = λ₃/λ₀ = {Q_real:.10f} (target 2/3 = {2/3:.10f})")
    print()

    # Step 4: Wirtinger Jacobians
    print("Step 4: Computing Wirtinger Jacobians (holomorphic + anti-holomorphic)...")
    J_holo, J_anti = compute_wirtinger_jacobians(state_eq, e_star, g, eps=1e-8)
    print("  Done.")
    print()

    # Eigenvalues of holomorphic Jacobian
    evals_holo_raw = np.linalg.eigvals(J_holo)
    evals_holo_mag = np.sort(np.abs(evals_holo_raw))[::-1]

    print("  Holomorphic Jacobian ∂Φ/∂z eigenvalues:")
    for i in range(8):
        ev = evals_holo_raw[np.argsort(-np.abs(evals_holo_raw))[i]]
        mag = abs(ev)
        if mag > 1e-10:
            print(f"    λ_{i} = {ev.real:+.10f} {ev.imag:+.10f}i,  |λ| = {mag:.10f}")
        else:
            print(f"    λ_{i} ≈ 0")

    # Eigenvalues of anti-holomorphic Jacobian
    evals_anti_raw = np.linalg.eigvals(J_anti)
    evals_anti_mag = np.sort(np.abs(evals_anti_raw))[::-1]

    print("\n  Anti-holomorphic Jacobian ∂Φ/∂z̄ eigenvalues:")
    for i in range(8):
        ev = evals_anti_raw[np.argsort(-np.abs(evals_anti_raw))[i]]
        mag = abs(ev)
        if mag > 1e-10:
            print(f"    λ_{i} = {ev.real:+.10f} {ev.imag:+.10f}i,  |λ| = {mag:.10f}")
        else:
            print(f"    λ_{i} ≈ 0")

    # Norms
    print(f"\n  ||∂Φ/∂z||_F  = {np.linalg.norm(J_holo, 'fro'):.6f}")
    print(f"  ||∂Φ/∂z̄||_F = {np.linalg.norm(J_anti, 'fro'):.6f}")

    # Ranks
    sv_holo = np.linalg.svd(J_holo, compute_uv=False)
    sv_anti = np.linalg.svd(J_anti, compute_uv=False)
    rank_holo = np.sum(sv_holo > 1e-8)
    rank_anti = np.sum(sv_anti > 1e-8)
    print(f"  rank(∂Φ/∂z)  = {rank_holo}")
    print(f"  rank(∂Φ/∂z̄) = {rank_anti}")
    print()

    # Step 5: THE COMPARISON
    print("=" * 70)
    print("THE COMPARISON: HOLOMORPHIC vs REAL vs ALGEBRAIC")
    print("=" * 70)
    print()

    # Extract top 4 nontrivial holomorphic eigenvalues
    holo_nontrivial = [m for m in evals_holo_mag if m > 1e-6]
    real_nontrivial = [m for m in evals_real if m > 1e-6]

    if len(holo_nontrivial) >= 4 and len(real_nontrivial) >= 4:
        Q_holo = holo_nontrivial[3] / holo_nontrivial[0]
        Q_real = real_nontrivial[3] / real_nontrivial[0]
        Q_target = 2.0/3

        print(f"  Q = λ₃/λ₀:")
        print(f"    Holomorphic (∂Φ/∂z):  {Q_holo:.10f}")
        print(f"    Real (DΦ|_R):          {Q_real:.10f}")
        print(f"    Algebraic (2/3):       {Q_target:.10f}")
        print()
        print(f"  Gaps from 2/3:")
        print(f"    Holomorphic: {abs(Q_holo - Q_target)/Q_target*100:.6f}%")
        print(f"    Real:        {abs(Q_real - Q_target)/Q_target*100:.6f}%")
        print()

        # Delta ratios
        Delta_holo = [-np.log(l) for l in holo_nontrivial[:4]]
        Delta_real = [-np.log(l) for l in real_nontrivial[:4]]

        DR_holo = Delta_holo[2] / Delta_holo[3] if Delta_holo[3] != 0 else 0
        DR_real = Delta_real[2] / Delta_real[3] if Delta_real[3] != 0 else 0
        DR_target = 20.0/21

        print(f"  Δ₂/Δ₃:")
        print(f"    Holomorphic: {DR_holo:.10f}")
        print(f"    Real:        {DR_real:.10f}")
        print(f"    Algebraic:   {DR_target:.10f}")
        print()
        print(f"  Gaps from 20/21:")
        print(f"    Holomorphic: {abs(DR_holo - DR_target)/DR_target*100:.6f}%")
        print(f"    Real:        {abs(DR_real - DR_target)/DR_target*100:.6f}%")
        print()

        # Eigenvalue comparison table
        print("  Eigenvalue comparison (|λ|):")
        print(f"  {'':5s} {'Holomorphic':>15s} {'Real':>15s} {'Difference':>15s}")
        print(f"  {'-'*50}")
        for i in range(4):
            h = holo_nontrivial[i]
            r = real_nontrivial[i]
            d = abs(h - r)
            print(f"    λ_{i}  {h:15.10f} {r:15.10f} {d:15.2e}")

    # Verdict
    print()
    print("=" * 70)
    if len(holo_nontrivial) >= 4:
        gap_holo = abs(Q_holo - Q_target)/Q_target*100
        gap_real = abs(Q_real - Q_target)/Q_target*100
        if gap_holo < gap_real / 2:
            print("VERDICT: Holomorphic eigenvalues CLOSER to Q=2/3.")
            print(f"  The 0.1% gap was contamination from the anti-holomorphic sector.")
            if gap_holo < 0.01:
                print(f"  Q_holo = {Q_holo:.10f} matches 2/3 to {gap_holo:.4f}%")
                print(f"  The Born floor correction lives in ∂̄Φ, not ∂Φ.")
        elif gap_holo > gap_real:
            print("VERDICT: Holomorphic eigenvalues FURTHER from Q=2/3.")
            print(f"  The real Jacobian was already the better estimate.")
        else:
            print("VERDICT: Similar gaps. The decomposition doesn't resolve the issue.")
        print(f"  Q_holo gap: {gap_holo:.6f}%")
        print(f"  Q_real gap: {gap_real:.6f}%")
    print("=" * 70)


if __name__ == "__main__":
    main()
