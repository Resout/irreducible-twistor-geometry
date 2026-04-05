#!/usr/bin/env python3
"""
SU(N) coupled mass gap: does the Lie-algebraic coupling preserve Δ > 0?

N-1 DS systems at H=3, one per simple root of su(N).
Coupled through the Cartan matrix A_ij.
Each system has K*=7/30, λ₀=0.2829 individually.

Question: does the COUPLED transfer operator on (CP³)^{N-1}
have spectral radius < 1?

Method: find the coupled fixed point, compute the coupled Jacobian,
extract eigenvalues. Test N = 2, 3, 4, 5, 6.
"""
import numpy as np
from scipy.optimize import brentq

H = 3
FLOOR = 1 / H**3


# ============================================================
# Single DS system (canonical floor from ds_gravity_k7_30.py)
# ============================================================
def ds_combine_real(m, e):
    """DS combination for real non-negative mass functions."""
    s, th = m[:3], m[3]
    se, the = e[:3], e[3]

    s_new = s * se + s * the + th * se
    th_new = th * the
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)

    denom = 1 - K
    if abs(denom) < 1e-15:
        return m.copy()

    out = np.zeros(4)
    out[:3] = s_new / denom
    out[3] = th_new / denom

    # Floor enforcement (correct quadratic formula)
    b = out[3]**2 / np.sum(out**2) if np.sum(out**2) > 0 else 1
    if b < FLOOR:
        abs_s = np.abs(out[:3])
        sum_s = np.sum(abs_s)
        sum_s_sq = np.sum(abs_s**2)
        if sum_s > 1e-15:
            r = sum_s_sq / sum_s**2
            # Solve (26-r)t² + 2rt - r = 0
            a_c = 26 - r
            b_c = 2 * r
            c_c = -r
            disc = b_c**2 - 4*a_c*c_c
            tn = (-b_c + np.sqrt(disc)) / (2*a_c)
            sc = (1.0 - tn) / sum_s
            out[:3] = abs_s * sc
            out[3] = tn

    return out


def find_single_equilibrium():
    """Find (m*, e*) for a single DS system at K*=7/30."""
    def K_at_eq(p_dom):
        p_w = (1 - p_dom) / 2
        sc = 1 - FLOOR
        raw = np.array([np.sqrt(p_dom*sc), np.sqrt(p_w*sc),
                        np.sqrt(p_w*sc), np.sqrt(FLOOR)])
        e = raw / np.sum(raw)
        m = np.array([0.4, 0.2, 0.2, 0.2])
        for _ in range(2000):
            m_new = ds_combine_real(m, e)
            if np.max(np.abs(m_new - m)) < 1e-15:
                break
            m = m_new
        s, th = m[:3], m[3]
        se, the = e[:3], e[3]
        K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
        return K

    # Find bracket by scanning
    lo, hi = 0.5, 0.99
    f_lo, f_hi = K_at_eq(lo) - 7/30, K_at_eq(hi) - 7/30
    for trial_lo in [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.92]:
        for trial_hi in [0.99, 0.98, 0.96, 0.94, 0.93]:
            fl = K_at_eq(trial_lo) - 7/30
            fh = K_at_eq(trial_hi) - 7/30
            if fl * fh < 0:
                lo, hi = trial_lo, trial_hi
                break
        else:
            continue
        break
    p_dom = brentq(lambda p: K_at_eq(p) - 7/30, lo, hi, xtol=1e-12)
    p_w = (1 - p_dom) / 2
    sc = 1 - FLOOR
    raw = np.array([np.sqrt(p_dom*sc), np.sqrt(p_w*sc),
                    np.sqrt(p_w*sc), np.sqrt(FLOOR)])
    e_star = raw / np.sum(raw)

    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(2000):
        m_new = ds_combine_real(m, e_star)
        if np.max(np.abs(m_new - m)) < 1e-15:
            break
        m = m_new

    return m, e_star


# ============================================================
# Cartan matrix for A_{N-1}
# ============================================================
def cartan_matrix(N):
    """Cartan matrix for su(N), type A_{N-1}."""
    r = N - 1  # rank
    A = np.zeros((r, r))
    for i in range(r):
        A[i, i] = 2
        if i > 0:
            A[i, i-1] = -1
        if i < r - 1:
            A[i, i+1] = -1
    return A


# ============================================================
# Coupled DS system
# ============================================================
def coupled_step(states, e_base, cartan, coupling=0.1):
    """
    One step of N-1 coupled DS systems.

    states: list of N-1 mass functions, each shape (4,)
    e_base: base evidence for each system
    cartan: (N-1)×(N-1) Cartan matrix
    coupling: strength of inter-system coupling

    Evidence for system i = e_base + coupling * sum_j A_ij * perturbation(state_j)
    """
    r = len(states)
    new_states = []

    for i in range(r):
        # Base evidence
        e_i = e_base.copy()

        # Coupling: modify evidence based on neighbors' states
        for j in range(r):
            if i != j and abs(cartan[i, j]) > 0.5:
                # Neighbor j influences system i's evidence
                # The perturbation is proportional to the state deviation from uniform
                uniform = np.ones(4) / 4
                perturbation = states[j] - uniform
                e_i = e_i + coupling * cartan[i, j] * perturbation

        # Ensure evidence is positive and normalized
        e_i = np.abs(e_i)
        e_i = np.maximum(e_i, 1e-10)
        e_i /= np.sum(e_i)

        new_states.append(ds_combine_real(states[i], e_i))

    return new_states


def find_coupled_equilibrium(N, coupling=0.1, n_iter=3000):
    """Find the coupled fixed point for SU(N)."""
    r = N - 1
    cartan = cartan_matrix(N)

    # Find single-system equilibrium for base evidence
    m_single, e_single = find_single_equilibrium()

    # Initialize all systems at the single-system equilibrium
    states = [m_single.copy() for _ in range(r)]

    for step in range(n_iter):
        new_states = coupled_step(states, e_single, cartan, coupling)
        max_diff = max(np.max(np.abs(new_states[i] - states[i])) for i in range(r))
        states = new_states
        if max_diff < 1e-14:
            break

    return states, e_single, cartan


def coupled_jacobian(states, e_base, cartan, coupling=0.1, eps=1e-7):
    """
    Compute the Jacobian of the coupled map at the fixed point.

    The state space is R^{4r} (r systems, 4 components each).
    The Jacobian is 4r × 4r.
    """
    r = len(states)
    dim = 4 * r

    # Flatten current state
    x0 = np.concatenate(states)

    def coupled_map(x_flat):
        ss = [x_flat[4*i:4*(i+1)] for i in range(r)]
        ss_new = coupled_step(ss, e_base, cartan, coupling)
        return np.concatenate(ss_new)

    f0 = coupled_map(x0)

    J = np.zeros((dim, dim))
    for j in range(dim):
        x_plus = x0.copy()
        x_minus = x0.copy()
        x_plus[j] += eps
        x_minus[j] -= eps
        J[:, j] = (coupled_map(x_plus) - coupled_map(x_minus)) / (2 * eps)

    return J


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("SU(N) COUPLED MASS GAP")
    print("=" * 70)

    # Single system reference
    m_single, e_single = find_single_equilibrium()
    print(f"\nSingle system: m* = {m_single}")
    print(f"Single system: e* = {e_single}")

    # Compute single-system Jacobian
    def single_map(m):
        return ds_combine_real(m, e_single)

    eps = 1e-7
    J_single = np.zeros((4, 4))
    for j in range(4):
        mp = m_single.copy(); mp[j] += eps
        mm = m_single.copy(); mm[j] -= eps
        J_single[:, j] = (single_map(mp) - single_map(mm)) / (2 * eps)

    evals_single = np.sort(np.abs(np.linalg.eigvals(J_single)))[::-1]
    print(f"Single system eigenvalues: {evals_single}")
    print(f"Single spectral radius: {evals_single[0]:.6f}")

    # Test coupled systems for N = 2, 3, 4, 5, 6
    print(f"\n{'='*70}")
    print(f"{'N':>3}  {'rank':>5}  {'dim':>5}  {'ρ(coupled)':>12}  {'ρ(single)':>12}  {'gap?':>6}  {'coupling iter':>14}")
    print(f"{'---':>3}  {'---':>5}  {'---':>5}  {'---':>12}  {'---':>12}  {'---':>6}  {'---':>14}")

    for coupling in [0.01, 0.05, 0.1]:
        print(f"\n  coupling = {coupling}")
        for N in range(2, 7):
            r = N - 1
            dim = 4 * r

            states, e_base, cartan = find_coupled_equilibrium(N, coupling)

            # Check convergence
            new_states = coupled_step(states, e_base, cartan, coupling)
            max_diff = max(np.max(np.abs(new_states[i] - states[i])) for i in range(r))

            J = coupled_jacobian(states, e_base, cartan, coupling)
            evals = np.sort(np.abs(np.linalg.eigvals(J)))[::-1]
            rho = evals[0]

            gap = "YES" if rho < 1 else "NO"
            print(f"  {N:3d}  {r:5d}  {dim:5d}  {rho:12.6f}  {evals_single[0]:12.6f}  {gap:>6}  {max_diff:.2e}")

    # Detailed eigenvalue spectrum for SU(3) and SU(5)
    print(f"\n{'='*70}")
    print("DETAILED EIGENVALUE SPECTRA (coupling = 0.05)")
    print("=" * 70)

    for N in [3, 5]:
        states, e_base, cartan = find_coupled_equilibrium(N, coupling=0.05)
        J = coupled_jacobian(states, e_base, cartan, coupling=0.05)
        evals = np.sort(np.abs(np.linalg.eigvals(J)))[::-1]
        print(f"\n  SU({N}): rank {N-1}, dim {4*(N-1)}")
        print(f"  Top 10 |eigenvalues|: {evals[:min(10,len(evals))]}")
        print(f"  Spectral radius: {evals[0]:.8f}")
        print(f"  Gap: Δ = {-np.log(evals[0]):.4f}" if evals[0] > 0 else "  Gap: ∞")


if __name__ == '__main__':
    np.random.seed(42)
    main()
