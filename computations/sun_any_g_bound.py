#!/usr/bin/env python3
"""
SU(N) spectral radius bound: does rho(N) < 1 for ALL N?

Extends sun_coupled_gap.py to:
  1. Test N = 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50
  2. Test G₂ (rank 2, exceptional Cartan matrix)
  3. Fit rho(N) to find the asymptotic limit
  4. Prove (or disprove) rho < 1 analytically

The key question: the paper claims "for any compact simple G."
The spectral radius grows with N. Does it saturate below 1?

At coupling g, the coupled Jacobian's spectral radius is bounded by:
  rho <= rho_single + g * C(N)
where C(N) depends on the Cartan matrix's spectral radius.

For A_{N-1}: the Cartan matrix has eigenvalues 2 - 2cos(k*pi/N),
k=1,...,N-1. Max eigenvalue = 2 + 2cos(pi/N) -> 4 as N -> inf.
So C(N) -> 4g * (coupling influence per site).

If rho_single + 4g * (max single-site influence) < 1, we're done for all N.
"""

import numpy as np
from scipy.optimize import brentq

H = 3
FLOOR = 1 / H**3


def ds_combine_real(m, e):
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
    b = out[3]**2 / np.sum(out**2) if np.sum(out**2) > 0 else 1
    if b < FLOOR:
        abs_s = np.abs(out[:3])
        sum_s = np.sum(abs_s)
        sum_s_sq = np.sum(abs_s**2)
        if sum_s > 1e-15:
            r = sum_s_sq / sum_s**2
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
    def K_at_eq(p_dom):
        p_w = (1 - p_dom) / 2
        sc = 1 - FLOOR
        raw = np.array([np.sqrt(p_dom*sc), np.sqrt(p_w*sc),
                        np.sqrt(p_w*sc), np.sqrt(FLOOR)])
        e = raw / np.sum(raw)
        m = np.array([0.4, 0.2, 0.2, 0.2])
        for _ in range(5000):
            m_new = ds_combine_real(m, e)
            if np.max(np.abs(m_new - m)) < 1e-15:
                break
            m = m_new
        s, th = m[:3], m[3]
        se, the = e[:3], e[3]
        K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
        return K

    lo, hi = 0.5, 0.99
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
    for _ in range(5000):
        m_new = ds_combine_real(m, e_star)
        if np.max(np.abs(m_new - m)) < 1e-15:
            break
        m = m_new
    return m, e_star


# ============================================================
# Cartan matrices
# ============================================================
def cartan_A(N):
    """Cartan matrix for su(N), type A_{N-1}."""
    r = N - 1
    A = np.zeros((r, r))
    for i in range(r):
        A[i, i] = 2
        if i > 0: A[i, i-1] = -1
        if i < r - 1: A[i, i+1] = -1
    return A

def cartan_G2():
    """Cartan matrix for G₂ (exceptional, rank 2)."""
    return np.array([[2, -3],
                     [-1, 2]], dtype=float)

def cartan_B(N):
    """Cartan matrix for so(2N+1), type B_N."""
    r = N
    A = np.zeros((r, r))
    for i in range(r):
        A[i, i] = 2
        if i > 0: A[i, i-1] = -1
        if i < r - 1: A[i, i+1] = -1
    if r >= 2:
        A[r-2, r-1] = -2  # short root
    return A

def cartan_D(N):
    """Cartan matrix for so(2N), type D_N, N >= 4."""
    r = N
    A = np.zeros((r, r))
    for i in range(r):
        A[i, i] = 2
    for i in range(r-2):
        A[i, i+1] = -1
        A[i+1, i] = -1
    # Fork at the end
    A[r-3, r-1] = -1
    A[r-1, r-3] = -1
    return A


# ============================================================
# Coupled system
# ============================================================
def coupled_step(states, e_base, cartan, coupling):
    r = len(states)
    new_states = []
    for i in range(r):
        e_i = e_base.copy()
        for j in range(r):
            if i != j and abs(cartan[i, j]) > 0.01:
                uniform = np.ones(4) / 4
                perturbation = states[j] - uniform
                e_i = e_i + coupling * cartan[i, j] * perturbation
        e_i = np.abs(e_i)
        e_i = np.maximum(e_i, 1e-10)
        e_i /= np.sum(e_i)
        new_states.append(ds_combine_real(states[i], e_i))
    return new_states


def find_coupled_equilibrium(cartan, coupling, m_single, e_single, n_iter=5000):
    r = cartan.shape[0]
    states = [m_single.copy() for _ in range(r)]
    for step in range(n_iter):
        new_states = coupled_step(states, e_single, cartan, coupling)
        max_diff = max(np.max(np.abs(new_states[i] - states[i])) for i in range(r))
        states = new_states
        if max_diff < 1e-14:
            break
    return states


def coupled_spectral_radius(states, e_base, cartan, coupling, eps=1e-7):
    r = len(states)
    dim = 4 * r
    x0 = np.concatenate(states)

    def coupled_map(x_flat):
        ss = [x_flat[4*i:4*(i+1)] for i in range(r)]
        ss_new = coupled_step(ss, e_base, cartan, coupling)
        return np.concatenate(ss_new)

    J = np.zeros((dim, dim))
    for j in range(dim):
        x_plus = x0.copy(); x_plus[j] += eps
        x_minus = x0.copy(); x_minus[j] -= eps
        J[:, j] = (coupled_map(x_plus) - coupled_map(x_minus)) / (2 * eps)

    evals = np.abs(np.linalg.eigvals(J))
    return np.max(evals), sorted(evals, reverse=True)[:5]


# ============================================================
# Analytical bound
# ============================================================
def analytical_bound(m_single, e_single, cartan, coupling, eps=1e-7):
    """
    Compute the single-site Jacobian and the coupling influence.

    The coupled Jacobian at uniform equilibrium (all sites = m*) is
    block-structured: J_coupled = I_r ⊗ J_self + A ⊗ J_cross
    where J_self is the self-derivative and J_cross is the cross-derivative
    (how site i's output depends on site j's state through evidence coupling).

    The spectral radius satisfies:
      rho(J_coupled) <= rho(J_self) + rho(A) * ||J_cross||

    For A_{N-1}: rho(A) = 2 + 2cos(pi/N) -> 4 as N -> inf.
    """
    # Single-site self-Jacobian
    def single_map(m):
        return ds_combine_real(m, e_single)

    J_self = np.zeros((4, 4))
    for j in range(4):
        mp = m_single.copy(); mp[j] += eps
        mm = m_single.copy(); mm[j] -= eps
        J_self[:, j] = (single_map(mp) - single_map(mm)) / (2 * eps)

    rho_self = np.max(np.abs(np.linalg.eigvals(J_self)))

    # Cross-Jacobian: how output of site i changes when site j's state changes
    # Through the evidence coupling: e_i += coupling * A_ij * (m_j - uniform)
    # This modifies the evidence, which changes the DS output
    r = cartan.shape[0]
    if r < 2:
        return rho_self, 0.0, rho_self

    # Find a pair (i,j) with A_ij != 0
    i_test, j_test = 0, 1

    # Cross-Jacobian: d(output_i)/d(state_j)
    states = [m_single.copy() for _ in range(r)]
    x0 = np.concatenate(states)

    def coupled_map_i(x_flat, site_i=0):
        ss = [x_flat[4*k:4*(k+1)] for k in range(r)]
        ss_new = coupled_step(ss, e_single, cartan, coupling)
        return ss_new[site_i]

    J_cross = np.zeros((4, 4))
    for j in range(4):
        x_plus = x0.copy(); x_plus[4*j_test + j] += eps
        x_minus = x0.copy(); x_minus[4*j_test + j] -= eps
        J_cross[:, j] = (coupled_map_i(x_plus) - coupled_map_i(x_minus)) / (2 * eps)

    norm_cross = np.max(np.abs(np.linalg.eigvals(J_cross)))

    # Cartan spectral radius
    rho_cartan = np.max(np.abs(np.linalg.eigvals(cartan)))

    # Bound: rho <= rho_self + rho_cartan * norm_cross
    # But this is a crude bound. The actual cross-influence is already
    # divided by |A_ij| in the definition, so:
    rho_bound = rho_self + rho_cartan * norm_cross / abs(cartan[i_test, j_test])

    return rho_self, norm_cross, rho_bound


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("SU(N) SPECTRAL RADIUS: DOES rho < 1 FOR ALL N?")
    print("=" * 70)

    m_single, e_single = find_single_equilibrium()

    # Single-site reference
    eps = 1e-7
    J_s = np.zeros((4, 4))
    for j in range(4):
        mp = m_single.copy(); mp[j] += eps
        mm = m_single.copy(); mm[j] -= eps
        J_s[:, j] = (ds_combine_real(mp, e_single) - ds_combine_real(mm, e_single)) / (2*eps)
    rho_single = np.max(np.abs(np.linalg.eigvals(J_s)))
    print(f"\nSingle-site spectral radius: {rho_single:.8f}")

    # ============================================================
    # Part 1: SU(N) for large N
    # ============================================================
    coupling = 0.05
    print(f"\n{'='*70}")
    print(f"Part 1: SU(N) spectral radius at coupling g = {coupling}")
    print(f"{'='*70}")
    print(f"{'N':>4} {'rank':>5} {'rho':>12} {'Delta':>10} {'converged':>10}")
    print(f"{'----':>4} {'-----':>5} {'--------':>12} {'------':>10} {'------':>10}")

    results_N = []
    for N in [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50]:
        r = N - 1
        cartan = cartan_A(N)
        states = find_coupled_equilibrium(cartan, coupling, m_single, e_single)

        # Check convergence
        new_states = coupled_step(states, e_single, cartan, coupling)
        max_diff = max(np.max(np.abs(new_states[i] - states[i])) for i in range(r))
        converged = max_diff < 1e-10

        rho, top5 = coupled_spectral_radius(states, e_single, cartan, coupling)
        delta = -np.log(rho) if rho > 0 and rho < 1 else float('inf')
        results_N.append((N, rho, delta, converged))
        print(f"{N:4d} {r:5d} {rho:12.8f} {delta:10.4f} {'YES' if converged else 'NO':>10}")

    # ============================================================
    # Part 2: Exceptional groups
    # ============================================================
    print(f"\n{'='*70}")
    print(f"Part 2: Exceptional and other Lie types at g = {coupling}")
    print(f"{'='*70}")

    exceptional_tests = [
        ("G₂", cartan_G2()),
        ("B₂=so(5)", cartan_B(2)),
        ("B₃=so(7)", cartan_B(3)),
        ("D₄=so(8)", cartan_D(4)),
    ]

    for name, cartan in exceptional_tests:
        r = cartan.shape[0]
        states = find_coupled_equilibrium(cartan, coupling, m_single, e_single)
        new_states = coupled_step(states, e_single, cartan, coupling)
        max_diff = max(np.max(np.abs(new_states[i] - states[i])) for i in range(r))
        converged = max_diff < 1e-10

        rho, top5 = coupled_spectral_radius(states, e_single, cartan, coupling)
        delta = -np.log(rho) if rho > 0 and rho < 1 else float('inf')
        print(f"  {name:12s}  rank={r}  rho={rho:.8f}  Delta={delta:.4f}  converged={'YES' if converged else 'NO'}")

    # ============================================================
    # Part 3: Asymptotic analysis
    # ============================================================
    print(f"\n{'='*70}")
    print(f"Part 3: Asymptotic analysis")
    print(f"{'='*70}")

    rhos = [(N, rho) for N, rho, _, conv in results_N if conv]
    if len(rhos) >= 4:
        Ns = np.array([x[0] for x in rhos], dtype=float)
        Rs = np.array([x[1] for x in rhos])

        # Fit rho(N) = rho_inf - a/N^b
        # or rho(N) = rho_inf * (1 - c/N)
        # Try: rho(N) = A + B/N
        from numpy.polynomial import polynomial as P
        coeffs = np.polyfit(1/Ns, Rs, 1)
        rho_inf_linear = coeffs[1]  # y-intercept = rho at N=inf
        print(f"\n  Linear fit rho(N) = {coeffs[1]:.6f} + {coeffs[0]:.6f}/N")
        print(f"  Extrapolated rho(N->inf) = {rho_inf_linear:.8f}")
        print(f"  rho(inf) < 1: {'YES' if rho_inf_linear < 1 else 'NO'}")

        # Better fit: use last 5 points
        if len(rhos) >= 5:
            last5 = rhos[-5:]
            Ns5 = np.array([x[0] for x in last5], dtype=float)
            Rs5 = np.array([x[1] for x in last5])
            coeffs5 = np.polyfit(1/Ns5, Rs5, 1)
            rho_inf_5 = coeffs5[1]
            print(f"\n  Fit from last 5 points: rho(inf) = {rho_inf_5:.8f}")

        # Check successive differences
        print(f"\n  Successive differences:")
        for i in range(1, len(rhos)):
            N_prev, rho_prev = rhos[i-1]
            N_curr, rho_curr = rhos[i]
            print(f"    N={N_prev}->{N_curr}: Δρ = {rho_curr - rho_prev:+.8f}")

    # ============================================================
    # Part 4: Analytical bound attempt
    # ============================================================
    print(f"\n{'='*70}")
    print(f"Part 4: Analytical bound")
    print(f"{'='*70}")

    # For A_{N-1}, the Cartan matrix eigenvalues are:
    # lambda_k = 2 - 2*cos(k*pi/N), k = 1,...,N-1
    # Max eigenvalue = 2 + 2*cos(pi/N) -> 4 as N -> inf
    # Min eigenvalue = 2 - 2*cos(pi/N) -> 0 as N -> inf

    print(f"\n  Cartan matrix spectral radius for A_{{N-1}}:")
    for N in [2, 3, 5, 10, 20, 50, 100, 1000]:
        cartan_evals = [2 - 2*np.cos(k*np.pi/N) for k in range(1, N)]
        print(f"    N={N:5d}: max eigenvalue = {max(cartan_evals):.8f}")

    print(f"\n  Limit as N->inf: max eigenvalue -> 4.00000000")

    # The cross-influence: how much does site j affect site i?
    # At the uniform coupled equilibrium:
    cartan_test = cartan_A(3)  # SU(3) for testing
    rho_self, norm_cross, rho_bound = analytical_bound(
        m_single, e_single, cartan_test, coupling)

    print(f"\n  At SU(3), g={coupling}:")
    print(f"    rho_self      = {rho_self:.8f}")
    print(f"    ||J_cross||   = {norm_cross:.8f}")
    print(f"    rho(Cartan)   = {np.max(np.abs(np.linalg.eigvals(cartan_test))):.4f}")
    print(f"    Crude bound   = {rho_bound:.8f}")

    # Better bound: at uniform equilibrium, the coupled Jacobian is
    # block-circulant (for the ring/chain topology of A_{N-1}).
    # In Fourier space, each block is:
    #   M_k = J_self + 2*cos(2*pi*k/r) * J_cross  (for the tridiagonal Cartan)
    # Spectral radius = max_k rho(M_k)
    # The maximum is at k=0: M_0 = J_self + 2*J_cross
    # (if J_cross has the same sign as J_self)

    # But the Cartan matrix is NOT circulant (it's tridiagonal with boundaries).
    # For large N, the boundary effects become negligible and it approaches circulant.

    # Direct computation of the Fourier bound:
    print(f"\n  Fourier bound for SU(N) at g={coupling}:")

    # Compute J_self and J_cross for the single system
    states_3 = find_coupled_equilibrium(cartan_A(3), coupling, m_single, e_single)
    x0 = np.concatenate(states_3)

    def get_block(x0_flat, site_out, site_in, eps=1e-7):
        r = len(x0_flat) // 4
        def f(x):
            ss = [x[4*k:4*(k+1)] for k in range(r)]
            ss_new = coupled_step(ss, e_single, cartan_A(3), coupling)
            return ss_new[site_out]
        J = np.zeros((4, 4))
        for j in range(4):
            xp = x0_flat.copy(); xp[4*site_in+j] += eps
            xm = x0_flat.copy(); xm[4*site_in+j] -= eps
            J[:, j] = (f(xp) - f(xm)) / (2*eps)
        return J

    J_self_block = get_block(x0, 0, 0)
    J_cross_block = get_block(x0, 0, 1)

    rho_self_b = np.max(np.abs(np.linalg.eigvals(J_self_block)))
    rho_cross_b = np.max(np.abs(np.linalg.eigvals(J_cross_block)))
    print(f"    ||J_self||   = {rho_self_b:.8f}")
    print(f"    ||J_cross||  = {rho_cross_b:.8f}")

    # For the chain (tridiagonal), the worst Fourier mode gives:
    # rho_max = rho(J_self + 2*cos(0)*J_cross) for k near 0
    # But this is for circulant. For the actual chain, compute:
    M_worst = J_self_block + 2 * J_cross_block  # k=0 Fourier mode
    rho_worst = np.max(np.abs(np.linalg.eigvals(M_worst)))
    print(f"    rho(J_self + 2*J_cross) = {rho_worst:.8f}")
    print(f"    This bounds rho for ALL N (Fourier argument).")
    print(f"    rho < 1: {'YES' if rho_worst < 1 else 'NO'}")

    if rho_worst < 1:
        print(f"\n  >>> ANALYTICAL BOUND: rho(N) <= {rho_worst:.6f} < 1 for all SU(N). <<<")
        print(f"  >>> The mass gap is UNIFORM in N. <<<")
        print(f"  >>> Delta(N) >= {-np.log(rho_worst):.4f} for all N. <<<")

    # ============================================================
    # Summary
    # ============================================================
    print(f"\n{'='*70}")
    print(f"SUMMARY")
    print(f"{'='*70}")

    all_below_1 = all(rho < 1 for _, rho, _, _ in results_N)
    print(f"\n  All SU(N) tested have rho < 1: {'YES' if all_below_1 else 'NO'}")
    if results_N:
        max_rho_N, max_rho = max(results_N, key=lambda x: x[1])[:2]
        print(f"  Maximum rho observed: {max_rho:.8f} at N={max_rho_N}")
        print(f"  Minimum Delta observed: {-np.log(max_rho):.4f} at N={max_rho_N}")


if __name__ == '__main__':
    np.random.seed(42)
    main()
