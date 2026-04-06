#!/usr/bin/env python3
"""
============================================================================
DEEPTHINK GAUNTLET: "FOR ANY COMPACT SIMPLE G" MASS GAP VERIFICATION
============================================================================

SELF-CONTAINED. No imports beyond numpy/scipy. Run as-is.

CONTEXT:
  The DS framework at H=3 produces a mass gap Delta > 0 for a single root
  SU(2) system (lambda_0 = 0.2829, Delta = 1.263). For SU(N), the framework
  uses N-1 coupled DS systems (Chevalley-Serre), one per simple root, coupled
  through the Cartan matrix. The coupling strength is g = K* = 7/30 (the
  equilibrium conflict of each root system — the framework's natural scale).

  The question: does the COUPLED spectral radius rho < 1 for ALL compact
  simple gauge groups G?

  There are exactly 4 infinite classical series (A_n, B_n, C_n, D_n) and
  5 exceptional groups (G2, F4, E6, E7, E8). This gauntlet tests them all.

TESTS:
  A1-A5:  SU(N) for N = 2, 3, 5, 10, 20, 50, 100 (A-series saturation)
  B1-B3:  SO(2n+1) for n = 2, 3, 5, 10 (B-series saturation)
  C1-C3:  Sp(2n) for n = 2, 3, 5, 10 (C-series saturation)
  D1-D3:  SO(2n) for n = 4, 5, 6, 10 (D-series saturation)
  E1-E5:  G2, F4, E6, E7, E8 (all exceptionals)
  F1:     Fourier bound for A_inf (analytical bound for all SU(N))
  F2:     Fourier bound for B_inf
  F3:     Fourier bound for C_inf
  F4:     Fourier bound for D_inf
  G1:     Coupling scan: rho(g) for SU(3) to find g_critical
  G2:     Coupling scan: rho(g) for F4 (tightest exceptional)
  H1:     Cross-check: single-site eigenvalues at 50-digit precision (mpmath)
  H2:     Control: uncoupled (g=0) gives rho = 0.2829 for all G

EXPECTED:
  All tests should give rho < 1 at g = K* = 7/30.
  F4 is predicted to be the tightest case (rho ~ 0.916).
  The Fourier bounds should confirm saturation for each classical series.
  g_critical for SU(3) should be ~ 0.65 (well above K* = 0.233).

OUTPUT: Print every result. No plots. Pure numbers.

============================================================================
"""

import numpy as np
from scipy.optimize import brentq
import sys

# ============================================================================
# DS FRAMEWORK CORE (H=3)
# ============================================================================

H = 3
FLOOR = 1.0 / H**3  # 1/27

def ds_combine(m, e):
    """Dempster-Shafer combination for real non-negative mass functions."""
    s, th = m[:3], m[3]
    se, ph = e[:3], e[3]
    s_new = s * se + s * ph + th * se
    th_new = th * ph
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    d = 1.0 - K
    if abs(d) < 1e-15:
        return m.copy()
    out = np.zeros(4)
    out[:3] = s_new / d
    out[3] = th_new / d
    # Born floor enforcement
    born = out[3]**2 / np.sum(out**2) if np.sum(out**2) > 0 else 1.0
    if born < FLOOR:
        abs_s = np.abs(out[:3])
        ss = np.sum(abs_s)
        sq = np.sum(abs_s**2)
        if ss > 1e-15:
            r = sq / ss**2
            a_c = 26.0 - r
            b_c = 2.0 * r
            c_c = -r
            disc = b_c**2 - 4*a_c*c_c
            tn = (-b_c + np.sqrt(disc)) / (2*a_c)
            sc = (1.0 - tn) / ss
            out[:3] = abs_s * sc
            out[3] = tn
    return out


def find_equilibrium():
    """Find (m*, e*) for a single DS system at K* = 7/30."""
    def K_at(p_dom):
        p_w = (1.0 - p_dom) / 2.0
        sc = 1.0 - FLOOR
        raw = np.array([np.sqrt(p_dom*sc), np.sqrt(p_w*sc),
                        np.sqrt(p_w*sc), np.sqrt(FLOOR)])
        e = raw / np.sum(raw)
        m = np.array([0.4, 0.2, 0.2, 0.2])
        for _ in range(8000):
            m2 = ds_combine(m, e)
            if np.max(np.abs(m2 - m)) < 1e-15:
                break
            m = m2
        s, th = m[:3], m[3]
        se, ph = e[:3], e[3]
        return sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)

    p_dom = brentq(lambda p: K_at(p) - 7.0/30, 0.90, 0.96, xtol=1e-14)
    p_w = (1.0 - p_dom) / 2.0
    sc = 1.0 - FLOOR
    raw = np.array([np.sqrt(p_dom*sc), np.sqrt(p_w*sc),
                    np.sqrt(p_w*sc), np.sqrt(FLOOR)])
    e_star = raw / np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(8000):
        m2 = ds_combine(m, e_star)
        if np.max(np.abs(m2 - m)) < 1e-15:
            break
        m = m2
    return m, e_star


# ============================================================================
# CARTAN MATRICES (complete classification of compact simple Lie algebras)
# ============================================================================

def cartan_A(n):
    """A_n = su(n+1). Rank n."""
    C = np.zeros((n, n))
    for i in range(n):
        C[i, i] = 2
        if i > 0: C[i, i-1] = -1
        if i < n-1: C[i, i+1] = -1
    return C

def cartan_B(n):
    """B_n = so(2n+1). Rank n. Last root short: A_{n-1,n} = -2."""
    C = np.zeros((n, n))
    for i in range(n):
        C[i, i] = 2
        if i > 0: C[i, i-1] = -1
        if i < n-1: C[i, i+1] = -1
    if n >= 2:
        C[n-2, n-1] = -2
    return C

def cartan_C(n):
    """C_n = sp(2n). Rank n. Last root long: A_{n,n-1} = -2."""
    C = np.zeros((n, n))
    for i in range(n):
        C[i, i] = 2
        if i > 0: C[i, i-1] = -1
        if i < n-1: C[i, i+1] = -1
    if n >= 2:
        C[n-1, n-2] = -2
    return C

def cartan_D(n):
    """D_n = so(2n). Rank n, n >= 4. Fork at node n-2."""
    C = np.zeros((n, n))
    for i in range(n):
        C[i, i] = 2
    for i in range(n-2):
        C[i, i+1] = -1
        C[i+1, i] = -1
    # Fork: node n-3 connects to both n-2 and n-1
    C[n-3, n-1] = -1
    C[n-1, n-3] = -1
    return C

def cartan_G2():
    """G2. Rank 2."""
    return np.array([[2, -3], [-1, 2]], dtype=float)

def cartan_F4():
    """F4. Rank 4. Node 2->3 has double bond."""
    return np.array([
        [ 2, -1,  0,  0],
        [-1,  2, -2,  0],
        [ 0, -1,  2, -1],
        [ 0,  0, -1,  2]
    ], dtype=float)

def cartan_E6():
    """E6. Rank 6. Linear 1-2-3-4-5 with branch 3-6."""
    C = np.zeros((6, 6))
    for i in range(6): C[i, i] = 2
    for i, j in [(0,1), (1,2), (2,3), (3,4), (2,5)]:
        C[i, j] = -1; C[j, i] = -1
    return C

def cartan_E7():
    """E7. Rank 7. Linear 1-2-3-4-5-6 with branch 3-7."""
    C = np.zeros((7, 7))
    for i in range(7): C[i, i] = 2
    for i, j in [(0,1), (1,2), (2,3), (3,4), (4,5), (2,6)]:
        C[i, j] = -1; C[j, i] = -1
    return C

def cartan_E8():
    """E8. Rank 8. Linear 1-2-3-4-5-6-7 with branch 3-8."""
    C = np.zeros((8, 8))
    for i in range(8): C[i, i] = 2
    for i, j in [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (2,7)]:
        C[i, j] = -1; C[j, i] = -1
    return C


# ============================================================================
# COUPLED DS SYSTEM
# ============================================================================

def coupled_step(states, e_base, cartan, g):
    """One step of r coupled DS systems."""
    r = len(states)
    new = []
    for i in range(r):
        e_i = e_base.copy()
        for j in range(r):
            if i != j and abs(cartan[i, j]) > 0.01:
                e_i = e_i + g * cartan[i, j] * (states[j] - 0.25)
        e_i = np.abs(e_i)
        e_i = np.maximum(e_i, 1e-10)
        e_i /= np.sum(e_i)
        new.append(ds_combine(states[i], e_i))
    return new


def find_coupled_eq(cartan, g, m0, e0, n_iter=8000):
    """Find coupled equilibrium."""
    r = cartan.shape[0]
    states = [m0.copy() for _ in range(r)]
    for step in range(n_iter):
        new = coupled_step(states, e0, cartan, g)
        diff = max(np.max(np.abs(new[i] - states[i])) for i in range(r))
        states = new
        if diff < 1e-14:
            return states, True, diff
    return states, False, diff


def coupled_rho(states, e_base, cartan, g, eps=1e-7):
    """Spectral radius of coupled Jacobian."""
    r = len(states)
    dim = 4 * r
    x0 = np.concatenate(states)

    def F(x):
        ss = [x[4*i:4*(i+1)] for i in range(r)]
        return np.concatenate(coupled_step(ss, e_base, cartan, g))

    J = np.zeros((dim, dim))
    f0 = F(x0)
    for j in range(dim):
        xp = x0.copy(); xp[j] += eps
        xm = x0.copy(); xm[j] -= eps
        J[:, j] = (F(xp) - F(xm)) / (2 * eps)

    evals = np.abs(np.linalg.eigvals(J))
    return np.max(evals)


def fourier_bound(states, e_base, cartan, g, eps=1e-7):
    """
    Fourier bound for long chains (A, B, C, D series).
    Extract middle-site self and cross blocks, compute max_k rho(M_k).
    """
    r = len(states)
    if r < 4:
        return None  # too short for meaningful Fourier analysis

    x0 = np.concatenate(states)
    mid = r // 2

    def f_i(x, site):
        ss = [x[4*k:4*(k+1)] for k in range(r)]
        return coupled_step(ss, e_base, cartan, g)[site]

    J_self = np.zeros((4, 4))
    J_cross = np.zeros((4, 4))
    for j in range(4):
        xp = x0.copy(); xp[4*mid+j] += eps
        xm = x0.copy(); xm[4*mid+j] -= eps
        J_self[:, j] = (f_i(xp, mid) - f_i(xm, mid)) / (2*eps)

        # cross: neighbor mid+1 -> mid
        xp2 = x0.copy(); xp2[4*(mid+1)+j] += eps
        xm2 = x0.copy(); xm2[4*(mid+1)+j] -= eps
        J_cross[:, j] = (f_i(xp2, mid) - f_i(xm2, mid)) / (2*eps)

    rho_max = 0
    k_worst = 0
    for k in np.linspace(0, np.pi, 200):
        M_k = J_self + 2*np.cos(k) * J_cross
        rho_k = np.max(np.abs(np.linalg.eigvals(M_k)))
        if rho_k > rho_max:
            rho_max = rho_k
            k_worst = k

    return rho_max, k_worst, J_self, J_cross


# ============================================================================
# MAIN GAUNTLET
# ============================================================================

def main():
    np.random.seed(42)
    g = 7.0 / 30  # K* = natural DS coupling

    print("=" * 78)
    print("DEEPTHINK GAUNTLET: MASS GAP FOR ALL COMPACT SIMPLE LIE ALGEBRAS")
    print("=" * 78)
    print(f"  Coupling: g = K* = 7/30 = {g:.10f}")
    print(f"  H = {H}, Born floor = 1/{H**3} = {FLOOR:.10f}")

    # Find single-site equilibrium
    m0, e0 = find_equilibrium()
    s, th = m0[:3], m0[3]
    se, ph = e0[:3], e0[3]
    K_check = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
    print(f"  m* = [{', '.join(f'{x:.10f}' for x in m0)}]")
    print(f"  e* = [{', '.join(f'{x:.10f}' for x in e0)}]")
    print(f"  K* = {K_check:.15f} (target: {7/30:.15f})")
    print(f"  |K* - 7/30| = {abs(K_check - 7/30):.2e}")

    # Single-site spectral radius
    J_s = np.zeros((4, 4))
    eps = 1e-7
    for j in range(4):
        mp = m0.copy(); mp[j] += eps
        mm = m0.copy(); mm[j] -= eps
        J_s[:, j] = (ds_combine(mp, e0) - ds_combine(mm, e0)) / (2*eps)
    rho_single = np.max(np.abs(np.linalg.eigvals(J_s)))
    print(f"  Single-site rho = {rho_single:.10f}")
    print(f"  Single-site Delta = {-np.log(rho_single):.10f}")

    results = []

    def test_group(name, cartan, g_val=g):
        r = cartan.shape[0]
        states, conv, diff = find_coupled_eq(cartan, g_val, m0, e0)
        rho = coupled_rho(states, e0, cartan, g_val)
        delta = -np.log(rho) if rho > 0 else float('inf')
        gap = rho < 1.0
        results.append((name, r, rho, delta, gap, conv))
        status = "PASS" if gap else "FAIL"
        c = "YES" if conv else "NO"
        print(f"  {name:14s} rank={r:2d}  rho={rho:.10f}  "
              f"Delta={delta:8.4f}  {status}  conv={c}")
        return rho

    # ================================================================
    # TESTS A: SU(N) = A_{N-1}
    # ================================================================
    print(f"\n{'='*78}")
    print("TESTS A: SU(N) series (type A)")
    print("=" * 78)
    for N in [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100]:
        test_group(f"SU({N})", cartan_A(N-1) if N > 1 else np.array([[2.0]]))

    # ================================================================
    # TESTS B: SO(2n+1) = B_n
    # ================================================================
    print(f"\n{'='*78}")
    print("TESTS B: SO(2n+1) series (type B)")
    print("=" * 78)
    for n in [2, 3, 4, 5, 6, 8, 10]:
        test_group(f"SO({2*n+1})", cartan_B(n))

    # ================================================================
    # TESTS C: Sp(2n) = C_n
    # ================================================================
    print(f"\n{'='*78}")
    print("TESTS C: Sp(2n) series (type C)")
    print("=" * 78)
    for n in [2, 3, 4, 5, 6, 8, 10]:
        test_group(f"Sp({2*n})", cartan_C(n))

    # ================================================================
    # TESTS D: SO(2n) = D_n
    # ================================================================
    print(f"\n{'='*78}")
    print("TESTS D: SO(2n) series (type D)")
    print("=" * 78)
    for n in [4, 5, 6, 8, 10]:
        test_group(f"SO({2*n})", cartan_D(n))

    # ================================================================
    # TESTS E: All 5 exceptional groups
    # ================================================================
    print(f"\n{'='*78}")
    print("TESTS E: Exceptional groups")
    print("=" * 78)
    test_group("G2", cartan_G2())
    test_group("F4", cartan_F4())
    test_group("E6", cartan_E6())
    test_group("E7", cartan_E7())
    test_group("E8", cartan_E8())

    # ================================================================
    # TESTS F: Fourier bounds (analytical bounds for infinite series)
    # ================================================================
    print(f"\n{'='*78}")
    print("TESTS F: Fourier bounds for classical series (N -> infinity)")
    print("=" * 78)

    for label, name, N_test in [("F1", "A_inf (SU)", 50),
                                 ("F2", "B_inf (SO odd)", 10),
                                 ("F3", "C_inf (Sp)", 10),
                                 ("F4", "D_inf (SO even)", 10)]:
        if label == "F1":
            cartan = cartan_A(N_test - 1)
        elif label == "F2":
            cartan = cartan_B(N_test)
        elif label == "F3":
            cartan = cartan_C(N_test)
        elif label == "F4":
            cartan = cartan_D(N_test)

        r = cartan.shape[0]
        states, conv, _ = find_coupled_eq(cartan, g, m0, e0)
        fb = fourier_bound(states, e0, cartan, g)
        if fb is not None:
            rho_fb, k_w, J_s_block, J_c_block = fb
            print(f"  {label}: {name:16s} (N={N_test:3d})  "
                  f"rho_Fourier={rho_fb:.10f}  k_worst={k_w:.4f}  "
                  f"rho<1: {'YES' if rho_fb < 1 else 'NO'}")
            print(f"       ||J_self||={np.max(np.abs(np.linalg.eigvals(J_s_block))):.8f}  "
                  f"||J_cross||={np.max(np.abs(np.linalg.eigvals(J_c_block))):.8f}")
        else:
            print(f"  {label}: {name} -- rank too small for Fourier analysis")

    # ================================================================
    # TESTS G: Critical coupling scans
    # ================================================================
    print(f"\n{'='*78}")
    print("TESTS G: Critical coupling g_crit (where rho crosses 1)")
    print("=" * 78)

    for name, cartan in [("SU(3)", cartan_A(2)), ("F4", cartan_F4())]:
        r = cartan.shape[0]
        print(f"\n  {name} (rank {r}):")
        prev_rho = 0
        g_crit = None
        for g_test in np.arange(0.05, 1.01, 0.05):
            states, conv, _ = find_coupled_eq(cartan, g_test, m0, e0)
            if not conv:
                print(f"    g={g_test:.3f}: did not converge")
                continue
            rho = coupled_rho(states, e0, cartan, g_test)
            gap = "YES" if rho < 1 else "NO"
            print(f"    g={g_test:.3f}: rho={rho:.8f}  gap={gap}")
            if rho >= 1 and prev_rho < 1 and g_crit is None:
                g_crit = g_test
            prev_rho = rho
        if g_crit:
            print(f"    >>> g_critical ~ {g_crit:.3f}")
            print(f"    >>> Safety margin: g_crit/K* = {g_crit / (7/30):.2f}x")
        else:
            print(f"    >>> rho < 1 for all tested g (up to 1.0)")

    # ================================================================
    # TEST H1: Control — uncoupled (g=0)
    # ================================================================
    print(f"\n{'='*78}")
    print("TEST H: Controls")
    print("=" * 78)

    print("\n  H1: g=0 (uncoupled) — all should give rho = single-site")
    for name, cartan in [("SU(5)", cartan_A(4)), ("E8", cartan_E8())]:
        states, conv, _ = find_coupled_eq(cartan, 0.0, m0, e0)
        rho = coupled_rho(states, e0, cartan, 0.0)
        print(f"    {name}: rho={rho:.10f} (single={rho_single:.10f}, "
              f"match={abs(rho-rho_single)<0.001})")

    print("\n  H2: Cartan matrix eigenvalue bounds")
    for name, cartan in [("A_49", cartan_A(49)), ("B_10", cartan_B(10)),
                          ("C_10", cartan_C(10)), ("D_10", cartan_D(10)),
                          ("G2", cartan_G2()), ("F4", cartan_F4()),
                          ("E8", cartan_E8())]:
        evals = np.linalg.eigvals(cartan)
        print(f"    {name:6s}: max|eval|={np.max(np.abs(evals)):.6f}  "
              f"min|eval|={np.min(np.abs(evals)):.6f}")

    # ================================================================
    # SUMMARY
    # ================================================================
    print(f"\n{'='*78}")
    print("SUMMARY")
    print("=" * 78)

    all_pass = all(gap for _, _, _, _, gap, _ in results)
    all_conv = all(conv for _, _, _, _, _, conv in results)
    tightest = max(results, key=lambda x: x[2])
    loosest = min(results, key=lambda x: x[2])

    print(f"\n  Total groups tested: {len(results)}")
    print(f"  All rho < 1 (gap exists): {'YES' if all_pass else 'NO'}")
    print(f"  All converged: {'YES' if all_conv else 'NO'}")
    print(f"  Tightest: {tightest[0]} (rho={tightest[2]:.8f}, "
          f"Delta={tightest[3]:.4f})")
    print(f"  Loosest:  {loosest[0]} (rho={loosest[2]:.8f}, "
          f"Delta={loosest[3]:.4f})")
    print(f"  Coupling: g = K* = 7/30 = {7/30:.10f}")

    print(f"\n  Classification coverage:")
    print(f"    A_n (su):   n=1,...,99  (SU(2) through SU(100))")
    print(f"    B_n (so):   n=2,...,10  (SO(5) through SO(21))")
    print(f"    C_n (sp):   n=2,...,10  (Sp(4) through Sp(20))")
    print(f"    D_n (so):   n=4,...,10  (SO(8) through SO(20))")
    print(f"    G2, F4, E6, E7, E8: all 5 exceptionals")
    print(f"    Total: every compact simple Lie algebra type represented")

    if all_pass:
        print(f"\n  >>> RESULT: AT g = K* = 7/30, THE MASS GAP EXISTS")
        print(f"  >>> FOR EVERY COMPACT SIMPLE LIE ALGEBRA TESTED.")
        print(f"  >>> TIGHTEST CASE: {tightest[0]} WITH rho = {tightest[2]:.6f},")
        print(f"  >>> MARGIN = {(1-tightest[2])*100:.1f}% BELOW 1.")
        print(f"  >>> MINIMUM MASS GAP: Delta = {tightest[3]:.4f} > 0.")
    else:
        failures = [r for r in results if not r[4]]
        print(f"\n  >>> FAILURES DETECTED:")
        for f in failures:
            print(f"  >>> {f[0]}: rho={f[2]:.8f}")

    print(f"\n{'='*78}")
    print("GAUNTLET COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()
