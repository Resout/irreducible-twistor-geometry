#!/usr/bin/env python3
"""
Non-triviality test: connected 4-point function at the K*=7/30 equilibrium.

Uses the CORRECT equilibrium (diverse evidence, not self-combination).
The equilibrium (m*, e*) is found by solving K(m*, e*) = 7/30 exactly.

A free (Gaussian) theory has G_n^conn = 0 for n >= 3.
Nonzero kurtosis excess proves the theory is interacting.
"""
import numpy as np
from scipy.optimize import brentq

H = 3
FLOOR = 1 / H**3  # 1/27


# ============================================================
# DS combination (from ds_gravity_k7_30.py)
# ============================================================
def ds_combine(m, e, apply_floor=True):
    s, th = m[:3], m[3]
    se, the = e[:3], e[3]

    s_new = s * se + s * the + th * se
    th_new = th * the
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)

    denom = 1 - K
    if abs(denom) < 1e-15:
        return m.copy(), K

    out = np.zeros(4, dtype=complex)
    out[:3] = s_new / denom
    out[3] = th_new / denom

    if apply_floor:
        b = np.abs(out[3])**2 / np.sum(np.abs(out)**2)
        if b < FLOOR:
            abs_s = np.abs(out[:3])
            ph_s = np.exp(1j * np.angle(out[:3]))
            ph_t = np.exp(1j * np.angle(out[3]))
            ss = np.sum(abs_s**2)
            tn = np.sqrt(FLOOR * ss / (1 - FLOOR))
            sc = (1.0 - tn) / np.sum(abs_s) if np.sum(abs_s) > 0 else 0
            out[:3] = ph_s * abs_s * sc
            out[3] = ph_t * tn

    return out, K


def born(m):
    return np.abs(m[3])**2 / np.sum(np.abs(m)**2)


def born_probs(m):
    sq = np.abs(m)**2
    return sq / np.sum(sq)


# ============================================================
# Find exact K*=7/30 equilibrium with diverse evidence
# ============================================================
def find_equilibrium():
    """Find (m*, e*) with K(m*, e*) = 7/30 exactly."""

    def K_at_eq(p_dom_val):
        p_w = (1 - p_dom_val) / 2
        sc = 1 - FLOOR
        raw = np.array([np.sqrt(p_dom_val * sc), np.sqrt(p_w * sc),
                        np.sqrt(p_w * sc), np.sqrt(FLOOR)])
        e = (raw / np.sum(raw)).astype(complex)
        m = np.array([0.4, 0.2, 0.2, 0.2], dtype=complex)
        for _ in range(2000):
            m_new, _ = ds_combine(m, e)
            if np.max(np.abs(m_new - m)) < 1e-15:
                break
            m = m_new
        _, K = ds_combine(m, e, apply_floor=False)
        return np.real(K)

    p_dom = brentq(lambda p: K_at_eq(p) - 7/30, 0.92, 0.94, xtol=1e-14)

    p_w = (1 - p_dom) / 2
    sc = 1 - FLOOR
    raw = np.array([np.sqrt(p_dom * sc), np.sqrt(p_w * sc),
                    np.sqrt(p_w * sc), np.sqrt(FLOOR)])
    e_star = (raw / np.sum(raw)).astype(complex)

    m = np.array([0.4, 0.2, 0.2, 0.2], dtype=complex)
    for _ in range(2000):
        m_new, _ = ds_combine(m, e_star)
        if np.max(np.abs(m_new - m)) < 1e-15:
            break
        m = m_new
    m_star = m_new

    return m_star, e_star


# ============================================================
# Main computation
# ============================================================
def main():
    print("=" * 70)
    print("NON-TRIVIALITY: 4-Point Function at K*=7/30 Equilibrium")
    print("=" * 70)

    m_star, e_star = find_equilibrium()
    _, K_check = ds_combine(m_star, e_star, apply_floor=False)
    p_star = born_probs(m_star)

    print(f"\nm* = ({', '.join(f'{np.real(x):.8f}' for x in m_star)})")
    print(f"e* = ({', '.join(f'{np.real(x):.8f}' for x in e_star)})")
    print(f"K(m*,e*) = {np.real(K_check):.15f}  (target: {7/30:.15f})")
    print(f"|K-7/30| = {abs(K_check - 7/30):.2e}")
    print(f"Born(θ*) = {np.real(born(m_star)):.10f}  (floor: {FLOOR:.10f})")
    print(f"Born probs = ({', '.join(f'{x:.6f}' for x in p_star)})")

    # Generate samples: perturb evidence around e*, combine with m*
    N = 500_000
    print(f"\nGenerating {N:,} samples: perturbed evidence around e*...")

    np.random.seed(42)

    p_out = np.zeros((N, H))  # singleton Born probs only

    for i in range(N):
        # Perturb evidence
        noise = np.random.randn(4) * 0.015
        e_pert = np.real(e_star) + noise
        e_pert = np.abs(e_pert)
        e_pert /= np.sum(e_pert)

        m_out, _ = ds_combine(m_star.copy(), e_pert.astype(complex))
        bp = born_probs(m_out)
        p_out[i] = np.real(bp[:H])

    # ================ TEST 1: Kurtosis ================
    print("\n" + "=" * 70)
    print("KURTOSIS EXCESS (Gaussian = 0)")
    print("=" * 70)

    se_kappa = np.sqrt(24.0 / N)
    for obs in range(H):
        dp = p_out[:, obs] - p_out[:, obs].mean()
        var = np.mean(dp**2)
        mu4 = np.mean(dp**4)
        kappa = mu4 / var**2 - 3
        sig = abs(kappa) / se_kappa
        print(f"  p_{obs}: κ = {kappa:+.6f} ± {se_kappa:.6f}  ({sig:.1f}σ)")

    # ================ TEST 2: Connected 4-point ================
    print("\n" + "=" * 70)
    print("CONNECTED 4-POINT (Wick subtracted)")
    print("=" * 70)

    fluct = p_out - p_out.mean(axis=0)

    G2 = np.zeros((H, H))
    for a in range(H):
        for b in range(H):
            G2[a, b] = np.mean(fluct[:, a] * fluct[:, b])

    print("\n  G2:")
    for a in range(H):
        print(f"    [{', '.join(f'{G2[a,b]:+.3e}' for b in range(H))}]")

    # (a,a,b,b) pairs
    print("\n  4-point (a,a,b,b):")
    for a in range(H):
        for b in range(a + 1, H):
            G4 = np.mean(fluct[:, a]**2 * fluct[:, b]**2)
            G4w = G2[a, a] * G2[b, b] + 2 * G2[a, b]**2
            G4c = G4 - G4w
            r = G4c / G4w if abs(G4w) > 1e-30 else 0
            print(f"    ({a}{a}{b}{b}): conn/Wick = {r:+.6f}")

    # (a,a,a,a) = kurtosis
    print("\n  4-point (a,a,a,a):")
    for a in range(H):
        G4 = np.mean(fluct[:, a]**4)
        G4w = 3 * G2[a, a]**2
        G4c = G4 - G4w
        r = G4c / G4w if abs(G4w) > 1e-30 else 0
        print(f"    ({a}{a}{a}{a}): conn/Wick = {r:+.6f}")

    # ================ TEST 3: Bootstrap ================
    print("\n" + "=" * 70)
    print("BOOTSTRAP (5000 resamples)")
    print("=" * 70)

    n_boot = 5000
    boot_kappa = np.zeros(n_boot)
    boot_ratio = np.zeros(n_boot)

    for bi in range(n_boot):
        idx = np.random.randint(0, N, size=N)
        bp = p_out[idx]
        fb = bp - bp.mean(axis=0)

        # Kurtosis of p_0
        dp = fb[:, 0]
        v = np.mean(dp**2)
        boot_kappa[bi] = np.mean(dp**4) / v**2 - 3 if v > 0 else 0

        # conn/wick (0,0,1,1)
        g2_00 = np.mean(fb[:, 0]**2)
        g2_11 = np.mean(fb[:, 1]**2)
        g2_01 = np.mean(fb[:, 0] * fb[:, 1])
        g4 = np.mean(fb[:, 0]**2 * fb[:, 1]**2)
        g4w = g2_00 * g2_11 + 2 * g2_01**2
        boot_ratio[bi] = (g4 - g4w) / g4w if abs(g4w) > 1e-30 else 0

    km, ks = boot_kappa.mean(), boot_kappa.std()
    rm, rs = boot_ratio.mean(), boot_ratio.std()

    print(f"  κ(p₀)         = {km:+.6f} ± {ks:.6f}  ({abs(km)/ks:.1f}σ)")
    print(f"  conn/Wick(0011)= {rm:+.6f} ± {rs:.6f}  ({abs(rm)/rs:.1f}σ)")

    # ================ ROBUSTNESS ================
    print("\n" + "=" * 70)
    print("ROBUSTNESS: varying σ")
    print("=" * 70)

    N_rob = 100_000
    for sigma in [0.005, 0.01, 0.015, 0.02, 0.03, 0.05]:
        p_rob = np.zeros((N_rob, H))
        for i in range(N_rob):
            noise = np.random.randn(4) * sigma
            e_pert = np.real(e_star) + noise
            e_pert = np.abs(e_pert)
            e_pert /= np.sum(e_pert)
            m_out, _ = ds_combine(m_star.copy(), e_pert.astype(complex))
            bp = born_probs(m_out)
            p_rob[i] = np.real(bp[:H])

        dp = p_rob[:, 0] - p_rob[:, 0].mean()
        v = np.mean(dp**2)
        kappa = np.mean(dp**4) / v**2 - 3 if v > 0 else 0
        se = np.sqrt(24.0 / N_rob)
        print(f"  σ={sigma:.3f}: κ(p₀)={kappa:+.6f} ({abs(kappa)/se:.1f}σ)")

    # ================ VERDICT ================
    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    print(f"  Equilibrium: K*=7/30 (verified to {abs(K_check-7/30):.1e})")
    print(f"  Kurtosis:    {km:+.6f} ± {ks:.6f} → {abs(km)/ks:.1f}σ from Gaussian")
    print(f"  conn/Wick:   {rm:+.6f} ± {rs:.6f} → {abs(rm)/rs:.1f}σ from Gaussian")

    if abs(km)/ks > 5:
        print(f"  DECISIVELY NON-GAUSSIAN: DS dynamics at K*=7/30 are interacting.")
    elif abs(km)/ks > 3:
        print(f"  NON-GAUSSIAN: evidence moderate.")
    else:
        print(f"  INCONCLUSIVE.")


if __name__ == '__main__':
    main()
