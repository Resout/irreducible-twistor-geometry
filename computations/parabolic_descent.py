"""
Parabolic Descent Computation
================================

Does DS combination with spatial neighbor-averaging produce diffusion?

Method:
1. Set up a 1D lattice with periodic boundary conditions
2. Initialize with a Born-floor mass function + small sinusoidal perturbation
3. Compare one DS step WITH vs WITHOUT spatial coupling
4. The DIFFERENCE is the spatial contribution
5. Check if this difference is proportional to -k^2 * sin(kx) (Laplacian signature)

If the effective coefficient nu > 0 for all perturbation directions and wavenumbers,
the spatial coupling is diffusive.

If nu is constant across wavenumbers, the operator is second-order (Laplacian).
"""

import numpy as np
import sys

H = 3


def ds_combine(m, e):
    """DS combination of mass function m with evidence e.

    m = (s1, s2, s3, theta), e = (e1, e2, e3, phi)

    Output components (pre-normalization):
        s_new,i = s_i * e_i + s_i * phi + theta * e_i
        theta_new = theta * phi
        K = sum_{i != j} s_i * e_j

    Normalized: m_out = (s_new, theta_new) / (1 - K)
    """
    s = m[:3]
    theta = m[3]
    ev = e[:3]
    phi = e[3]

    s_new = s * ev + s * phi + theta * ev
    theta_new = theta * phi

    K = 0.0
    for i in range(3):
        for j in range(3):
            if i != j:
                K += s[i] * ev[j]

    denom = 1.0 - K
    if abs(denom) < 1e-15:
        return m.copy(), K

    m_out = np.array([s_new[0], s_new[1], s_new[2], theta_new]) / denom
    return m_out, K


def born_enforce(m):
    """Enforce Born floor: Born(Theta) = theta^2 / sum(m_i^2) >= 1/27."""
    norm_sq = np.sum(m**2)
    if norm_sq < 1e-30:
        return m
    born = m[3]**2 / norm_sq
    if born < 1.0 / 27.0:
        s = m[:3]
        c = np.sqrt(np.sum(s**2) / 26.0)
        sign = 1.0 if m[3] >= 0 else -1.0
        m_new = np.array([s[0], s[1], s[2], sign * c])
        l1 = np.sum(np.abs(m_new))
        if l1 > 1e-30:
            m_new /= l1
        return m_new
    return m


def test_diffusion(m0, v, k, N=200, eps=1e-6):
    """Test spatial coupling for perturbation direction v at wavenumber k.

    Sets up lattice: m(x) = m0 + eps * sin(kx) * v  (L1-normalized)
    Runs one step with and without spatial coupling.
    Returns effective diffusion coefficient nu.
    """
    # Build lattice with sinusoidal perturbation
    lattice = np.zeros((N, 4))
    for x in range(N):
        lattice[x] = m0 + eps * np.sin(k * x) * v
        l1 = np.sum(np.abs(lattice[x]))
        if l1 > 1e-30:
            lattice[x] /= l1

    # One step WITH spatial coupling (evidence = average of neighbors)
    lat_s = np.zeros_like(lattice)
    for x in range(N):
        e = (lattice[(x + 1) % N] + lattice[(x - 1) % N]) / 2.0
        m_out, _ = ds_combine(lattice[x], e)
        m_out = born_enforce(m_out)
        lat_s[x] = m_out

    # One step WITHOUT spatial coupling (evidence = self)
    lat_ns = np.zeros_like(lattice)
    for x in range(N):
        m_out, _ = ds_combine(lattice[x], lattice[x])
        m_out = born_enforce(m_out)
        lat_ns[x] = m_out

    # Spatial contribution = difference
    spatial = lat_s - lat_ns

    # Project onto perturbation direction
    proj = np.array([np.dot(spatial[x], v) for x in range(N)])

    # Fit: proj = A * sin(kx)
    sin_kx = np.array([np.sin(k * x) for x in range(N)])
    denom = np.dot(sin_kx, sin_kx)
    if denom < 1e-30:
        return 0.0, 0.0
    A = np.dot(proj, sin_kx) / denom

    # For diffusion: spatial effect = nu * d^2/dx^2 [eps * sin(kx)] = -nu * k^2 * eps * sin(kx)
    # So A = -nu * k^2 * eps => nu = -A / (k^2 * eps)
    if abs(k**2 * eps) < 1e-30:
        return 0.0, 0.0
    nu = -A / (k**2 * eps)

    # Also compute residual (how well sin(kx) fits)
    fitted = A * sin_kx
    residual = np.sqrt(np.sum((proj - fitted)**2) / np.sum(proj**2 + 1e-30))

    return nu, residual


def main():
    print("=" * 70)
    print("PARABOLIC DESCENT COMPUTATION")
    print("Does DS spatial coupling produce diffusion?")
    print("=" * 70)
    print()

    # === BASE POINT 1: Symmetric, on Born floor ===
    theta0 = 1.0 / (3.0 * np.sqrt(26.0 / 3.0) + 1.0)
    s0 = theta0 * np.sqrt(26.0 / 3.0)
    m0_sym = np.array([s0, s0, s0, theta0])
    born0 = m0_sym[3]**2 / np.sum(m0_sym**2)

    print(f"BASE POINT 1: Symmetric, on Born floor")
    print(f"  m = [{s0:.6f}, {s0:.6f}, {s0:.6f}, {theta0:.6f}]")
    print(f"  L1 = {np.sum(m0_sym):.8f}, Born = {born0:.8f}, target = {1/27:.8f}")
    m_self, K0 = ds_combine(m0_sym, m0_sym)
    print(f"  Self-combination K = {K0:.6f}")
    print()

    # === BASE POINT 2: Asymmetric, on Born floor ===
    m0_asym = np.array([0.35, 0.28, 0.20, 0.17])
    m0_asym /= np.sum(m0_asym)
    m0_asym = born_enforce(m0_asym)
    born1 = m0_asym[3]**2 / np.sum(m0_asym**2)

    print(f"BASE POINT 2: Asymmetric, on Born floor")
    print(f"  m = [{m0_asym[0]:.6f}, {m0_asym[1]:.6f}, {m0_asym[2]:.6f}, {m0_asym[3]:.6f}]")
    print(f"  L1 = {np.sum(m0_asym):.8f}, Born = {born1:.8f}")
    print()

    # === BASE POINT 3: Inside Born floor (floor inactive) ===
    m0_inside = np.array([0.20, 0.20, 0.20, 0.40])
    m0_inside /= np.sum(m0_inside)
    born2 = m0_inside[3]**2 / np.sum(m0_inside**2)

    print(f"BASE POINT 3: Inside Born floor (floor inactive)")
    print(f"  m = [{m0_inside[0]:.6f}, {m0_inside[1]:.6f}, {m0_inside[2]:.6f}, {m0_inside[3]:.6f}]")
    print(f"  L1 = {np.sum(m0_inside):.8f}, Born = {born2:.8f} (> {1/27:.4f})")
    print()

    N = 200

    # Perturbation directions
    directions = {
        "s1":       np.array([1.0, 0.0, 0.0, 0.0]),
        "s2":       np.array([0.0, 1.0, 0.0, 0.0]),
        "s3":       np.array([0.0, 0.0, 1.0, 0.0]),
        "theta":    np.array([0.0, 0.0, 0.0, 1.0]),
        "s1-s2":    np.array([1.0, -1.0, 0.0, 0.0]) / np.sqrt(2),
        "s1-theta": np.array([1.0, 0.0, 0.0, -1.0]) / np.sqrt(2),
    }

    # Wavenumber modes
    modes = [1, 2, 3, 5, 10, 20]

    for bp_name, m0 in [("SYMMETRIC (Born boundary)", m0_sym),
                         ("ASYMMETRIC (Born boundary)", m0_asym),
                         ("INSIDE BORN FLOOR", m0_inside)]:
        print(f"--- {bp_name} ---")
        print(f"{'Direction':<12s}", end="")
        for n in modes:
            print(f" | n={n:>2d}", end="")
        print(" | const?")
        print("-" * (12 + 10 * len(modes) + 10))

        for dname, v in directions.items():
            nus = []
            for n in modes:
                k = 2 * np.pi * n / N
                nu, res = test_diffusion(m0, v, k, N=N, eps=1e-6)
                nus.append(nu)
            # Check constancy
            nu_arr = np.array(nus)
            if abs(np.mean(nu_arr)) > 1e-12:
                spread = np.std(nu_arr) / abs(np.mean(nu_arr))
            else:
                spread = 0.0

            line = f"{dname:<12s}"
            for nu in nus:
                line += f" | {nu:>+.4f}"
            sign = "YES" if spread < 0.1 else f"NO ({spread:.0%})"
            line += f" | {sign}"
            print(line)
        print()

    # === DETAILED: check k^2 scaling for one direction ===
    print("=" * 70)
    print("DETAILED k^2 SCALING TEST (symmetric base, s1 perturbation)")
    print("=" * 70)
    v = np.array([1.0, 0.0, 0.0, 0.0])
    print(f"{'n':>4s} {'k':>8s} {'k^2':>10s} {'nu':>12s} {'residual':>10s}")
    print("-" * 50)
    for n in [1, 2, 3, 5, 8, 10, 15, 20, 30, 40]:
        k = 2 * np.pi * n / N
        nu, res = test_diffusion(m0_sym, v, k, N=N, eps=1e-6)
        print(f"{n:>4d} {k:>8.4f} {k**2:>10.4f} {nu:>12.8f} {res:>10.4f}")

    print()
    print("=" * 70)
    print("VERDICT")
    print("=" * 70)
    # Collect all results for verdict
    all_positive = True
    all_constant = True
    for bp_name, m0 in [("sym", m0_sym), ("asym", m0_asym), ("inside", m0_inside)]:
        for dname, v in directions.items():
            nus = []
            for n in modes:
                k = 2 * np.pi * n / N
                nu, _ = test_diffusion(m0, v, k, N=N, eps=1e-6)
                nus.append(nu)
            nu_arr = np.array(nus)
            if any(nu_arr < -1e-10):
                all_positive = False
                print(f"  NEGATIVE nu found: {bp_name}/{dname}, nu = {min(nus):.6f}")
            if abs(np.mean(nu_arr)) > 1e-10:
                spread = np.std(nu_arr) / abs(np.mean(nu_arr))
                if spread > 0.2:
                    all_constant = False

    if all_positive:
        print("  All effective diffusion coefficients nu > 0: DIFFUSIVE")
    else:
        print("  Some nu < 0: ANTI-DIFFUSIVE in some directions")

    if all_constant:
        print("  nu approximately constant across wavenumbers: LAPLACIAN (2nd order)")
    else:
        print("  nu varies with wavenumber: HIGHER-ORDER operator")

    if all_positive and all_constant:
        print()
        print("  RESULT: DS spatial coupling IS a positive-definite Laplacian.")
        print("  The parabolic descent holds.")
    elif all_positive:
        print()
        print("  RESULT: Diffusive but not purely Laplacian.")
    else:
        print()
        print("  RESULT: NOT purely diffusive. Further analysis needed.")


if __name__ == "__main__":
    main()
