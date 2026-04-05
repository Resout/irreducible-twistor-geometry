"""
Correct Abelian Discrimination Test

The CORRECT test: observables that DON'T share links.

For SU(2): non-abelian gauge constraint propagates correlations through
plaquette coupling even between spatially separated observables.
δ > 1/H should PERSIST as N → ∞.

For U(1) with coprime charges: representations factorise regardless
of shared links. δ → 1/H as N → ∞.

Test design:
  1. Plaquettes at MAXIMUM SEPARATION on the lattice (no shared links)
     - P_near: plaquette at site (0,0,0,0) in (0,1) plane
     - P_far:  plaquette at site (L/2, L/2, L/2, L/2) in (2,3) plane
     These share ZERO links. Any correlation is from the gauge constraint.

  2. U(1) coprime test: Wilson loop with charge e=1 vs charge e=3
     at the SAME plaquette. These share links but the representations
     are coprime. Peter-Weyl orthogonality says they factorise.

  3. Measure δ (diagonal content of conditional distribution) vs N
     for all three cases. The framework predicts:
     - SU(2) separated: δ > 1/3, persists
     - U(1) separated: δ → 1/3
     - U(1) coprime: δ → 1/3
"""

import numpy as np
from time import time

H = 3
I2 = np.eye(2, dtype=complex)
sigma = [
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
]


def random_su2():
    a = np.random.randn(4)
    a /= np.linalg.norm(a)
    return a[0]*I2 + 1j*(a[1]*sigma[0] + a[2]*sigma[1] + a[3]*sigma[2])


def random_su2_near_identity(eps=0.3):
    """SU(2) element close to identity — better Metropolis acceptance."""
    a = np.array([1.0, 0, 0, 0]) + eps * np.random.randn(4)
    a /= np.linalg.norm(a)
    return a[0]*I2 + 1j*(a[1]*sigma[0] + a[2]*sigma[1] + a[3]*sigma[2])


# ═══════════════════════════════════════════════════════════════
#  U(1) Lattice
# ═══════════════════════════════════════════════════════════════

class U1Lattice:
    def __init__(self, L, beta):
        self.L = L
        self.beta = beta
        # links[x,y,z,t,mu] = phase angle
        self.links = np.random.uniform(0, 2*np.pi, (L, L, L, L, 4))

    def plaquette_phase(self, x, y, z, t, mu, nu):
        """Phase of plaquette at (x,y,z,t) in (mu,nu) plane."""
        L = self.L
        pos = [x, y, z, t]
        pos_mu = pos.copy(); pos_mu[mu] = (pos_mu[mu] + 1) % L
        pos_nu = pos.copy(); pos_nu[nu] = (pos_nu[nu] + 1) % L
        return (self.links[x, y, z, t, mu] +
                self.links[tuple(pos_mu)][nu] -
                self.links[tuple(pos_nu)][mu] -
                self.links[x, y, z, t, nu])

    def plaquette_trace(self, x, y, z, t, mu, nu, charge=1):
        """Re(exp(i * charge * plaquette_phase))."""
        return np.cos(charge * self.plaquette_phase(x, y, z, t, mu, nu))

    def sweep(self):
        L = self.L
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        for mu in range(4):
                            old = self.links[x, y, z, t, mu]
                            new = old + np.random.uniform(-0.5, 0.5)

                            dS = 0.0
                            for nu in range(4):
                                if nu == mu:
                                    continue
                                pos = [x, y, z, t]
                                pos_mu = pos.copy(); pos_mu[mu] = (pos_mu[mu]+1) % L
                                pos_nu = pos.copy(); pos_nu[nu] = (pos_nu[nu]+1) % L
                                pos_mnu = pos.copy(); pos_mnu[nu] = (pos_mnu[nu]-1) % L
                                pos_mu_mnu = pos.copy()
                                pos_mu_mnu[mu] = (pos_mu_mnu[mu]+1) % L
                                pos_mu_mnu[nu] = (pos_mu_mnu[nu]-1) % L

                                # Forward plaquette
                                p_old = old + self.links[tuple(pos_mu)][nu] - self.links[tuple(pos_nu)][mu] - self.links[x,y,z,t,nu]
                                p_new = new + self.links[tuple(pos_mu)][nu] - self.links[tuple(pos_nu)][mu] - self.links[x,y,z,t,nu]
                                dS -= self.beta * (np.cos(p_new) - np.cos(p_old))

                                # Backward plaquette
                                pb_old = self.links[tuple(pos_mnu)][nu] + old - self.links[tuple(pos_mu_mnu)][nu] - self.links[tuple(pos_mnu)][mu]
                                pb_new = self.links[tuple(pos_mnu)][nu] + new - self.links[tuple(pos_mu_mnu)][nu] - self.links[tuple(pos_mnu)][mu]
                                dS -= self.beta * (np.cos(pb_new) - np.cos(pb_old))

                            if dS < 0 or np.random.random() < np.exp(-dS):
                                self.links[x, y, z, t, mu] = new


# ═══════════════════════════════════════════════════════════════
#  SU(2) Lattice
# ═══════════════════════════════════════════════════════════════

class SU2Lattice:
    def __init__(self, L, beta):
        self.L = L
        self.beta = beta
        # links[x,y,z,t,mu] = 2x2 complex matrix
        self.links = np.zeros((L, L, L, L, 4, 2, 2), dtype=complex)
        for idx in np.ndindex(L, L, L, L):
            for mu in range(4):
                self.links[idx + (mu,)] = random_su2()

    def plaquette_trace(self, x, y, z, t, mu, nu):
        """Re(Tr(plaquette)) / 2, normalized to [-1, 1]."""
        L = self.L
        pos = (x, y, z, t)
        pos_mu = list(pos); pos_mu[mu] = (pos_mu[mu]+1) % L; pos_mu = tuple(pos_mu)
        pos_nu = list(pos); pos_nu[nu] = (pos_nu[nu]+1) % L; pos_nu = tuple(pos_nu)

        U1 = self.links[pos + (mu,)]
        U2 = self.links[pos_mu + (nu,)]
        U3 = self.links[pos_nu + (mu,)].conj().T
        U4 = self.links[pos + (nu,)].conj().T

        return np.real(np.trace(U1 @ U2 @ U3 @ U4)) / 2.0

    def staple(self, x, y, z, t, mu):
        """Sum of staples around link (x,y,z,t,mu)."""
        L = self.L
        pos = (x, y, z, t)
        S = np.zeros((2, 2), dtype=complex)

        for nu in range(4):
            if nu == mu:
                continue
            pos_mu = list(pos); pos_mu[mu] = (pos_mu[mu]+1) % L; pos_mu = tuple(pos_mu)
            pos_nu = list(pos); pos_nu[nu] = (pos_nu[nu]+1) % L; pos_nu = tuple(pos_nu)
            pos_mnu = list(pos); pos_mnu[nu] = (pos_mnu[nu]-1) % L; pos_mnu = tuple(pos_mnu)
            pos_mu_mnu = list(pos)
            pos_mu_mnu[mu] = (pos_mu_mnu[mu]+1) % L
            pos_mu_mnu[nu] = (pos_mu_mnu[nu]-1) % L
            pos_mu_mnu = tuple(pos_mu_mnu)

            # Forward staple
            S += (self.links[pos_mu + (nu,)] @
                  self.links[pos_nu + (mu,)].conj().T @
                  self.links[pos + (nu,)].conj().T)
            # Backward staple
            S += (self.links[pos_mu_mnu + (nu,)].conj().T @
                  self.links[pos_mnu + (mu,)].conj().T @
                  self.links[pos_mnu + (nu,)])
        return S

    def sweep(self, n_hits=8):
        """Metropolis sweep with multiple hits per link."""
        L = self.L
        accepted = 0
        total = 0
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        for mu in range(4):
                            S = self.staple(x, y, z, t, mu)
                            old_U = self.links[x, y, z, t, mu].copy()
                            old_action = np.real(np.trace(old_U @ S))

                            for _ in range(n_hits):
                                trial = random_su2_near_identity(0.3) @ old_U
                                trial_action = np.real(np.trace(trial @ S))
                                dS = self.beta / 2 * (trial_action - old_action)
                                total += 1
                                if dS > 0 or np.random.random() < np.exp(dS):
                                    self.links[x, y, z, t, mu] = trial
                                    old_U = trial
                                    old_action = trial_action
                                    accepted += 1

        return accepted / total if total > 0 else 0


# ═══════════════════════════════════════════════════════════════
#  Delta measurement
# ═══════════════════════════════════════════════════════════════

def compute_delta(obs1, obs2):
    """Diagonal content of the conditional distribution C.
    δ = 1: deterministic. δ = 1/H: independent."""
    n = len(obs1)
    # Equal-count tercile binning
    def bin_tercile(vals):
        order = np.argsort(vals)
        bins = np.zeros(n, dtype=int)
        for rank, i in enumerate(order):
            bins[i] = min(rank * H // n, H - 1)
        return bins

    b1 = bin_tercile(obs1)
    b2 = bin_tercile(obs2)

    C = np.zeros((H, H))
    for i in range(n):
        C[b1[i], b2[i]] += 1

    row_sums = C.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    C /= row_sums

    return np.trace(C) / H, C


# ═══════════════════════════════════════════════════════════════
#  Main test
# ═══════════════════════════════════════════════════════════════

def main():
    np.random.seed(42)
    L = 4  # 4^4 = 256 sites — manageable in pure Python
    n_therm = 50
    n_decorr = 5  # sweeps between measurements
    n_meas = 1000  # total measurements

    beta_u1 = 1.0
    beta_su2 = 2.2

    half = L // 2

    print("=" * 70)
    print("CORRECT ABELIAN DISCRIMINATION TEST")
    print("=" * 70)
    print(f"Lattice: {L}^4 = {L**4} sites")
    print(f"U(1) β = {beta_u1}, SU(2) β = {beta_su2}")
    print(f"Thermalisation: {n_therm} sweeps")
    print(f"Decorrelation: {n_decorr} sweeps between measurements")
    print(f"Measurements: {n_meas}")
    print()

    # ── U(1) ──
    print("--- U(1) ---")
    t0 = time()
    lat_u1 = U1Lattice(L, beta_u1)

    print(f"Thermalizing U(1)...")
    for s in range(n_therm):
        lat_u1.sweep()

    # Measure site-level plaquettes at maximum separation
    # O1: plaquette at (0,0,0,0) in (0,1) plane
    # O2: plaquette at (half,half,half,half) in (2,3) plane — NO shared links
    # O3: plaquette at (0,0,0,0) in (0,1) plane with charge 3 (coprime test)

    u1_O1 = []  # charge-1 plaquette at origin, (0,1) plane
    u1_O2 = []  # charge-1 plaquette at far site, (2,3) plane
    u1_O3 = []  # charge-3 plaquette at origin, (0,1) plane (coprime)

    # Also collect ALL site plaquettes for better statistics
    u1_near_all = []  # all plaquettes in (0,1) plane
    u1_far_all = []   # all plaquettes in (2,3) plane, shifted by (half,half,half,half)

    print(f"Measuring U(1)...")
    for m in range(n_meas):
        for _ in range(n_decorr):
            lat_u1.sweep()

        # Single-site observables
        u1_O1.append(lat_u1.plaquette_trace(0, 0, 0, 0, 0, 1, charge=1))
        u1_O2.append(lat_u1.plaquette_trace(half, half, half, half, 2, 3, charge=1))
        u1_O3.append(lat_u1.plaquette_trace(0, 0, 0, 0, 0, 1, charge=3))

        # All-site observables (more statistics per config)
        near_vals = []
        far_vals = []
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        near_vals.append(lat_u1.plaquette_trace(x, y, z, t, 0, 1))
                        fx = (x + half) % L
                        fy = (y + half) % L
                        fz = (z + half) % L
                        ft = (t + half) % L
                        far_vals.append(lat_u1.plaquette_trace(fx, fy, fz, ft, 2, 3))
        u1_near_all.extend(near_vals)
        u1_far_all.extend(far_vals)

        if (m + 1) % 200 == 0:
            dt = time() - t0
            print(f"  {m+1}/{n_meas} ({dt:.0f}s)")

    u1_O1 = np.array(u1_O1)
    u1_O2 = np.array(u1_O2)
    u1_O3 = np.array(u1_O3)
    u1_near_all = np.array(u1_near_all)
    u1_far_all = np.array(u1_far_all)

    dt = time() - t0
    print(f"  U(1) done in {dt:.0f}s")
    print()

    # ── SU(2) ──
    print("--- SU(2) ---")
    t0 = time()
    lat_su2 = SU2Lattice(L, beta_su2)

    print(f"Thermalizing SU(2)...")
    for s in range(n_therm):
        acc = lat_su2.sweep()
        if s % 10 == 0:
            print(f"  Sweep {s}, acceptance: {acc:.2f}")

    su2_near_all = []
    su2_far_all = []

    print(f"Measuring SU(2)...")
    for m in range(n_meas):
        for _ in range(n_decorr):
            lat_su2.sweep(n_hits=4)

        near_vals = []
        far_vals = []
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        near_vals.append(lat_su2.plaquette_trace(x, y, z, t, 0, 1))
                        fx = (x + half) % L
                        fy = (y + half) % L
                        fz = (z + half) % L
                        ft = (t + half) % L
                        far_vals.append(lat_su2.plaquette_trace(fx, fy, fz, ft, 2, 3))
        su2_near_all.extend(near_vals)
        su2_far_all.extend(far_vals)

        if (m + 1) % 100 == 0:
            dt = time() - t0
            print(f"  {m+1}/{n_meas} ({dt:.0f}s)")

    su2_near_all = np.array(su2_near_all)
    su2_far_all = np.array(su2_far_all)

    dt = time() - t0
    print(f"  SU(2) done in {dt:.0f}s")
    print()

    # ── Results ──
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    print()

    # Test 1: U(1) separated plaquettes (single-site, time series)
    d1, C1 = compute_delta(u1_O1, u1_O2)
    print(f"Test 1: U(1) single-site, separated, charge 1")
    print(f"  δ = {d1:.4f} (target 1/H = {1/H:.4f})")
    print(f"  C = \n{np.array2string(C1, precision=3)}")
    print()

    # Test 2: U(1) coprime (same site, charge 1 vs charge 3)
    d2, C2 = compute_delta(u1_O1, u1_O3)
    print(f"Test 2: U(1) coprime charges e=1 vs e=3 (same site)")
    print(f"  δ = {d2:.4f} (target 1/H = {1/H:.4f})")
    print(f"  C = \n{np.array2string(C2, precision=3)}")
    print()

    # Test 3: U(1) all-site separated (big N)
    # Sub-sample at various N to see scaling
    print(f"Test 3: U(1) separated, δ vs N (all-site)")
    for N in [500, 1000, 2000, 5000, 10000, 50000, len(u1_near_all)]:
        if N > len(u1_near_all):
            continue
        ds = []
        for trial in range(20):
            idx = np.random.choice(len(u1_near_all), N, replace=False)
            d, _ = compute_delta(u1_near_all[idx], u1_far_all[idx])
            ds.append(d)
        print(f"  N={N:>7d}: δ = {np.mean(ds):.4f} ± {np.std(ds):.4f}")
    print()

    # Test 4: SU(2) all-site separated (big N)
    print(f"Test 4: SU(2) separated, δ vs N (all-site)")
    for N in [500, 1000, 2000, 5000, 10000, 50000, len(su2_near_all)]:
        if N > len(su2_near_all):
            continue
        ds = []
        for trial in range(20):
            idx = np.random.choice(len(su2_near_all), N, replace=False)
            d, _ = compute_delta(su2_near_all[idx], su2_far_all[idx])
            ds.append(d)
        print(f"  N={N:>7d}: δ = {np.mean(ds):.4f} ± {np.std(ds):.4f}")
    print()

    # Summary
    print("=" * 70)
    print("FRAMEWORK PREDICTION vs OBSERVATION")
    print("=" * 70)
    print(f"  1/H = {1/H:.4f}")
    print()
    print(f"  SU(2) separated:  framework predicts δ > 1/H, persisting with N")
    print(f"  U(1) separated:   framework predicts δ → 1/H as N → ∞")
    print(f"  U(1) coprime:     framework predicts δ = 1/H exactly")
    print()
    print("  If SU(2) δ stays above U(1) δ as N grows: CONFIRMED")
    print("  If both converge to 1/H: FALSIFIED")


if __name__ == "__main__":
    main()
