"""
Abelian Discrimination Test (Prediction 1)

The framework predicts:
  - SU(2): δ < 1 structurally, persists as N → ∞
  - U(1): δ → 1/H = 1/3 as N → ∞ (binning artefact only)

This is the single most important empirical test. If the curves
don't separate, the framework has a problem.

Method:
  1. Generate N samples from each theory (SU(2) lattice, U(1) lattice)
  2. Construct Wilson loops, bin into H=3
  3. Compute diagonal content δ of the conditional distribution C
  4. Repeat at N = 500, 1000, 2000, 5000, 10000, 50000
  5. Plot δ vs N for both theories

For SU(2): use Haar-random SU(2) matrices, Wilson action
For U(1): use uniform random phases, compact U(1) action
Wilson loops: plaquette traces in two orthogonal planes
"""

import numpy as np
from collections import Counter

H = 3


def random_su2():
    """Random SU(2) matrix via Haar measure (quaternion method)."""
    # Random unit quaternion
    x = np.random.randn(4)
    x /= np.linalg.norm(x)
    a, b, c, d = x
    return np.array([
        [a + 1j*b, c + 1j*d],
        [-c + 1j*d, a - 1j*b]
    ])


def random_u1():
    """Random U(1) phase."""
    return np.exp(2j * np.pi * np.random.random())


def su2_plaquette(links, x, y, mu, nu, L):
    """Compute SU(2) plaquette trace at position (x,y) in plane (mu, nu)."""
    # U_mu(x) * U_nu(x+mu) * U_mu(x+nu)^dag * U_nu(x)^dag
    def idx(px, py):
        return (px % L, py % L)

    if mu == 0 and nu == 1:
        U1 = links[0][idx(x, y)]
        U2 = links[1][idx(x+1, y)]
        U3 = links[0][idx(x, y+1)].conj().T
        U4 = links[1][idx(x, y)].conj().T
    else:
        U1 = links[1][idx(x, y)]
        U2 = links[0][idx(x, y+1)]
        U3 = links[1][idx(x+1, y)].conj().T
        U4 = links[0][idx(x, y)].conj().T

    return np.real(np.trace(U1 @ U2 @ U3 @ U4))


def u1_plaquette(links, x, y, mu, nu, L):
    """Compute U(1) plaquette phase at position (x,y)."""
    def idx(px, py):
        return (px % L, py % L)

    if mu == 0 and nu == 1:
        phase = (links[0][idx(x, y)] * links[1][idx(x+1, y)] *
                 np.conj(links[0][idx(x, y+1)]) * np.conj(links[1][idx(x, y)]))
    else:
        phase = (links[1][idx(x, y)] * links[0][idx(x, y+1)] *
                 np.conj(links[1][idx(x+1, y)]) * np.conj(links[0][idx(x, y)]))

    return np.real(phase)


def generate_su2_lattice(L, beta, n_therm=200, n_meas=1):
    """Generate thermalized SU(2) lattice config via heat bath."""
    # Initialize random
    links = [[random_su2() for _ in range(L)] for _ in range(2)]
    links = [
        {(x, y): random_su2() for x in range(L) for y in range(L)}
        for _ in range(2)
    ]

    # Simple Metropolis update
    for sweep in range(n_therm):
        for mu in range(2):
            for x in range(L):
                for y in range(L):
                    old = links[mu][(x, y)]
                    # Compute staple
                    old_action = 0
                    for nu in range(2):
                        if nu != mu:
                            old_action += su2_plaquette(links, x, y, mu, nu, L)
                            # Also the plaquette in the other direction
                            if mu == 0:
                                old_action += su2_plaquette(links, x-1, y, mu, nu, L)
                            else:
                                old_action += su2_plaquette(links, x, y-1, mu, nu, L)

                    # Propose new link
                    proposal = old @ random_su2()
                    # Make it close to old for better acceptance
                    eps = 0.3
                    delta = np.eye(2) + eps * (random_su2() - np.eye(2))
                    delta /= np.sqrt(np.real(np.linalg.det(delta)))
                    proposal = delta @ old

                    links[mu][(x, y)] = proposal
                    new_action = 0
                    for nu in range(2):
                        if nu != mu:
                            new_action += su2_plaquette(links, x, y, mu, nu, L)
                            if mu == 0:
                                new_action += su2_plaquette(links, x-1, y, mu, nu, L)
                            else:
                                new_action += su2_plaquette(links, x, y-1, mu, nu, L)

                    dS = -beta * (new_action - old_action)
                    if dS > 0 or np.random.random() < np.exp(dS):
                        pass  # accept
                    else:
                        links[mu][(x, y)] = old  # reject

    return links


def generate_u1_lattice(L, beta, n_therm=200):
    """Generate thermalized compact U(1) lattice config."""
    links = [
        {(x, y): random_u1() for x in range(L) for y in range(L)}
        for _ in range(2)
    ]

    for sweep in range(n_therm):
        for mu in range(2):
            for x in range(L):
                for y in range(L):
                    old = links[mu][(x, y)]
                    old_action = 0
                    for nu in range(2):
                        if nu != mu:
                            old_action += u1_plaquette(links, x, y, mu, nu, L)
                            if mu == 0:
                                old_action += u1_plaquette(links, x-1, y, mu, nu, L)
                            else:
                                old_action += u1_plaquette(links, x, y-1, mu, nu, L)

                    # Propose
                    proposal = old * np.exp(1j * np.random.uniform(-0.5, 0.5))
                    links[mu][(x, y)] = proposal

                    new_action = 0
                    for nu in range(2):
                        if nu != mu:
                            new_action += u1_plaquette(links, x, y, mu, nu, L)
                            if mu == 0:
                                new_action += u1_plaquette(links, x-1, y, mu, nu, L)
                            else:
                                new_action += u1_plaquette(links, x, y-1, mu, nu, L)

                    dS = -beta * (new_action - old_action)
                    if dS > 0 or np.random.random() < np.exp(dS):
                        pass
                    else:
                        links[mu][(x, y)] = old

    return links


def extract_observables(links, L, plaquette_fn):
    """Extract two Wilson loop observables from orthogonal planes."""
    obs1 = []  # plaquettes in (0,1) plane
    obs2 = []  # plaquettes at shifted positions

    for x in range(L):
        for y in range(L):
            obs1.append(plaquette_fn(links, x, y, 0, 1, L))
            # Second observable: plaquette at (x+1, y+1) — overlapping but different
            obs2.append(plaquette_fn(links, (x+1) % L, (y+1) % L, 0, 1, L))

    return np.array(obs1), np.array(obs2)


def compute_delta(obs1, obs2, N_subsample=None):
    """Compute diagonal content δ of the conditional distribution.

    δ = (1/H) * sum_i C_ii where C is the [H,H] conditional distribution.
    δ = 1: deterministic
    δ = 1/H: independent
    1/H < δ < 1: correlated
    """
    if N_subsample is not None and N_subsample < len(obs1):
        idx = np.random.choice(len(obs1), N_subsample, replace=False)
        obs1 = obs1[idx]
        obs2 = obs2[idx]

    # Bin into H=3 by equal-count terciles
    def bin_tercile(vals):
        sorted_idx = np.argsort(vals)
        bins = np.zeros(len(vals), dtype=int)
        n = len(vals)
        for rank, i in enumerate(sorted_idx):
            bins[i] = min(rank * H // n, H - 1)
        return bins

    b1 = bin_tercile(obs1)
    b2 = bin_tercile(obs2)

    # Build conditional distribution C
    C = np.zeros((H, H))
    for i in range(len(b1)):
        C[b1[i], b2[i]] += 1

    # Normalize rows
    row_sums = C.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    C = C / row_sums

    delta = np.trace(C) / H
    return delta, C


def run_test():
    """Run the abelian discrimination test."""
    L = 8  # lattice size
    beta_su2 = 2.3
    beta_u1 = 1.0
    n_configs = 100  # number of independent configurations

    sample_sizes = [500, 1000, 2000, 5000, 10000]

    print("=" * 70)
    print("ABELIAN DISCRIMINATION TEST")
    print("=" * 70)
    print(f"Lattice: {L}x{L}, SU(2) β={beta_su2}, U(1) β={beta_u1}")
    print(f"Configurations: {n_configs}")
    print(f"Sample sizes: {sample_sizes}")
    print()

    # Collect all observables
    print("Generating SU(2) configurations...")
    su2_obs1_all, su2_obs2_all = [], []
    for c in range(n_configs):
        links = generate_su2_lattice(L, beta_su2, n_therm=100)
        o1, o2 = extract_observables(links, L, su2_plaquette)
        su2_obs1_all.extend(o1)
        su2_obs2_all.extend(o2)
        if (c + 1) % 20 == 0:
            print(f"  SU(2): {c+1}/{n_configs} configs")

    su2_obs1_all = np.array(su2_obs1_all)
    su2_obs2_all = np.array(su2_obs2_all)
    print(f"  Total SU(2) samples: {len(su2_obs1_all)}")

    print("Generating U(1) configurations...")
    u1_obs1_all, u1_obs2_all = [], []
    for c in range(n_configs):
        links = generate_u1_lattice(L, beta_u1, n_therm=100)
        o1, o2 = extract_observables(links, L, u1_plaquette)
        u1_obs1_all.extend(o1)
        u1_obs2_all.extend(o2)
        if (c + 1) % 20 == 0:
            print(f"  U(1): {c+1}/{n_configs} configs")

    u1_obs1_all = np.array(u1_obs1_all)
    u1_obs2_all = np.array(u1_obs2_all)
    print(f"  Total U(1) samples: {len(u1_obs1_all)}")
    print()

    # Compute δ at each sample size
    print(f"{'N':>8s}  {'δ(SU2)':>10s}  {'δ(U1)':>10s}  {'SU2-1/3':>10s}  {'U1-1/3':>10s}  {'Separated?':>10s}")
    print("-" * 70)

    n_trials = 20  # repeat each measurement for error bars

    for N in sample_sizes:
        if N > len(su2_obs1_all):
            print(f"{N:>8d}  (insufficient samples)")
            continue

        su2_deltas = []
        u1_deltas = []

        for trial in range(n_trials):
            d_su2, _ = compute_delta(su2_obs1_all, su2_obs2_all, N_subsample=N)
            d_u1, _ = compute_delta(u1_obs1_all, u1_obs2_all, N_subsample=N)
            su2_deltas.append(d_su2)
            u1_deltas.append(d_u1)

        su2_mean = np.mean(su2_deltas)
        su2_std = np.std(su2_deltas)
        u1_mean = np.mean(u1_deltas)
        u1_std = np.std(u1_deltas)

        separated = "YES" if (su2_mean - su2_std) > (u1_mean + u1_std) else "no"

        print(f"{N:>8d}  {su2_mean:>8.4f}±{su2_std:.4f}  {u1_mean:>8.4f}±{u1_std:.4f}  "
              f"{su2_mean - 1/3:>+10.4f}  {u1_mean - 1/3:>+10.4f}  {separated:>10s}")

    print()
    print("Framework prediction:")
    print("  SU(2): δ should stay above 1/3, constant or increasing with N")
    print("  U(1):  δ should approach 1/3 = 0.3333 as N → ∞")
    print("  If both approach 1/3: framework has a problem")
    print("  If SU(2) stays high and U(1) drops: framework confirmed")


if __name__ == "__main__":
    np.random.seed(42)
    run_test()
