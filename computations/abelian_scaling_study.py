"""
Abelian Discrimination Scaling Study

The definitive test: does U(1) coprime δ converge to 1/3 while
SU(2) cross-rep δ stays above 1/3 as the lattice grows?

Run the representation mixing test at L = 4, 6, 8.
Same β values, same observables, same analysis.
The ONLY thing that changes is the lattice size.

If the curves separate with increasing L: confirmed.
If both converge to 1/3: the discrimination fails.
If U(1) coprime plateaus above 1/3: Peter-Weyl argument has a hole.

This script is designed to run unattended for ~1 hour.
"""

import numpy as np
from time import time
import sys

H = 3
I2 = np.eye(2, dtype=complex)
sigma = [np.array([[0,1],[1,0]], dtype=complex),
         np.array([[0,-1j],[1j,0]], dtype=complex),
         np.array([[1,0],[0,-1]], dtype=complex)]

def random_su2():
    a = np.random.randn(4); a /= np.linalg.norm(a)
    return a[0]*I2 + 1j*(a[1]*sigma[0]+a[2]*sigma[1]+a[3]*sigma[2])

def near_id(eps=0.3):
    a = np.array([1.0,0,0,0]) + eps*np.random.randn(4)
    a /= np.linalg.norm(a)
    return a[0]*I2 + 1j*(a[1]*sigma[0]+a[2]*sigma[1]+a[3]*sigma[2])


class SU2Lattice:
    def __init__(self, L):
        self.L = L
        self.links = np.zeros((L,L,L,L,4,2,2), dtype=complex)
        for idx in np.ndindex(L,L,L,L):
            for mu in range(4):
                self.links[idx+(mu,)] = random_su2()

    def sweep(self, beta, n_hits=4):
        L = self.L
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        for mu in range(4):
                            pos = (x,y,z,t)
                            S = np.zeros((2,2), dtype=complex)
                            for nu in range(4):
                                if nu == mu: continue
                                pm = list(pos); pm[mu]=(pm[mu]+1)%L; pm=tuple(pm)
                                pn = list(pos); pn[nu]=(pn[nu]+1)%L; pn=tuple(pn)
                                pmn = list(pos); pmn[nu]=(pmn[nu]-1)%L; pmn=tuple(pmn)
                                pmun = list(pos); pmun[mu]=(pmun[mu]+1)%L; pmun[nu]=(pmun[nu]-1)%L; pmun=tuple(pmun)
                                S += self.links[pm+(nu,)] @ self.links[pn+(mu,)].conj().T @ self.links[pos+(nu,)].conj().T
                                S += self.links[pmun+(nu,)].conj().T @ self.links[pmn+(mu,)].conj().T @ self.links[pmn+(nu,)]
                            old = self.links[pos+(mu,)].copy()
                            oa = np.real(np.trace(old @ S))
                            for _ in range(n_hits):
                                trial = near_id() @ old
                                ta = np.real(np.trace(trial @ S))
                                if ta > oa or np.random.random() < np.exp(beta/2*(ta-oa)):
                                    self.links[pos+(mu,)] = trial; old = trial; oa = ta

    def measure(self):
        """Return (j12_01, j1_02, j12_02) for all sites."""
        L = self.L
        j12_01, j1_02, j12_02 = [], [], []
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        pos = (x,y,z,t)
                        pm0 = list(pos); pm0[0]=(pm0[0]+1)%L; pm0=tuple(pm0)
                        p1 = list(pos); p1[1]=(p1[1]+1)%L; p1=tuple(p1)
                        p2 = list(pos); p2[2]=(p2[2]+1)%L; p2=tuple(p2)

                        plaq01 = self.links[pos+(0,)] @ self.links[pm0+(1,)] @ self.links[p1+(0,)].conj().T @ self.links[pos+(1,)].conj().T
                        plaq02 = self.links[pos+(0,)] @ self.links[pm0+(2,)] @ self.links[p2+(0,)].conj().T @ self.links[pos+(2,)].conj().T

                        j12_01.append(np.real(np.trace(plaq01)) / 2)
                        tr02 = np.trace(plaq02)
                        j1_02.append((np.abs(tr02)**2 - 1) / 3)
                        j12_02.append(np.real(tr02) / 2)
        return np.array(j12_01), np.array(j1_02), np.array(j12_02)


class U1Lattice:
    def __init__(self, L):
        self.L = L
        self.links = np.random.uniform(0, 2*np.pi, (L,L,L,L,4))

    def sweep(self, beta):
        L = self.L
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        for mu in range(4):
                            old = self.links[x,y,z,t,mu]
                            new = old + np.random.uniform(-0.5, 0.5)
                            dS = 0.0
                            for nu in range(4):
                                if nu == mu: continue
                                pos = [x,y,z,t]
                                pm = pos.copy(); pm[mu]=(pm[mu]+1)%L
                                pn = pos.copy(); pn[nu]=(pn[nu]+1)%L
                                pmn = pos.copy(); pmn[nu]=(pmn[nu]-1)%L
                                pmun = pos.copy(); pmun[mu]=(pmun[mu]+1)%L; pmun[nu]=(pmun[nu]-1)%L
                                p_old = old + self.links[tuple(pm)][nu] - self.links[tuple(pn)][mu] - self.links[x,y,z,t,nu]
                                p_new = new + self.links[tuple(pm)][nu] - self.links[tuple(pn)][mu] - self.links[x,y,z,t,nu]
                                dS -= beta*(np.cos(p_new) - np.cos(p_old))
                                pb_old = self.links[tuple(pmn)][nu]+old-self.links[tuple(pmun)][nu]-self.links[tuple(pmn)][mu]
                                pb_new = self.links[tuple(pmn)][nu]+new-self.links[tuple(pmun)][nu]-self.links[tuple(pmn)][mu]
                                dS -= beta*(np.cos(pb_new) - np.cos(pb_old))
                            if dS < 0 or np.random.rand() < np.exp(-dS):
                                self.links[x,y,z,t,mu] = new

    def measure(self):
        """Return (e1_01, e3_02, e1_02) for all sites."""
        L = self.L
        e1_01, e3_02, e1_02 = [], [], []
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        pos = [x,y,z,t]
                        pm0 = pos.copy(); pm0[0]=(pm0[0]+1)%L
                        p1 = pos.copy(); p1[1]=(p1[1]+1)%L
                        p2 = pos.copy(); p2[2]=(p2[2]+1)%L

                        ph01 = self.links[x,y,z,t,0] + self.links[tuple(pm0)][1] - self.links[tuple(p1)][0] - self.links[x,y,z,t,1]
                        ph02 = self.links[x,y,z,t,0] + self.links[tuple(pm0)][2] - self.links[tuple(p2)][0] - self.links[x,y,z,t,2]

                        e1_01.append(np.cos(ph01))
                        e3_02.append(np.cos(3 * ph02))
                        e1_02.append(np.cos(ph02))
        return np.array(e1_01), np.array(e3_02), np.array(e1_02)


def compute_delta(o1, o2, N=None):
    if N is not None and N < len(o1):
        idx = np.random.choice(len(o1), N, replace=False)
        o1, o2 = o1[idx], o2[idx]
    n = len(o1)
    order1 = np.argsort(o1); order2 = np.argsort(o2)
    b1 = np.zeros(n, dtype=int); b2 = np.zeros(n, dtype=int)
    for rank, i in enumerate(order1): b1[i] = min(rank*H//n, H-1)
    for rank, i in enumerate(order2): b2[i] = min(rank*H//n, H-1)
    C = np.zeros((H,H))
    for i in range(n): C[b1[i], b2[i]] += 1
    rs = C.sum(1, keepdims=True); rs[rs==0]=1; C /= rs
    return np.trace(C)/H


def run_lattice_size(L, beta_su2, beta_u1, n_therm, n_meas, n_decorr):
    """Run the full test at one lattice size."""
    n_sites = L**4
    print(f"\n{'='*70}")
    print(f"  LATTICE SIZE L = {L}  ({n_sites} sites)")
    print(f"  SU(2) beta = {beta_su2}, U(1) beta = {beta_u1}")
    print(f"  Therm: {n_therm}, Meas: {n_meas}, Decorr: {n_decorr}")
    print(f"{'='*70}")

    # ── SU(2) ──
    print(f"\n  SU(2): thermalizing ({n_therm} sweeps)...")
    t0 = time()
    su2 = SU2Lattice(L)
    for s in range(n_therm):
        su2.sweep(beta_su2, n_hits=4)
        if (s+1) % 20 == 0:
            print(f"    Sweep {s+1}/{n_therm} ({time()-t0:.0f}s)")

    print(f"  SU(2): measuring ({n_meas} configs)...")
    su2_j12_01_all = []
    su2_j1_02_all = []
    su2_j12_02_all = []
    t_meas = time()
    for m in range(n_meas):
        for _ in range(n_decorr):
            su2.sweep(beta_su2, n_hits=3)
        j12_01, j1_02, j12_02 = su2.measure()
        su2_j12_01_all.extend(j12_01)
        su2_j1_02_all.extend(j1_02)
        su2_j12_02_all.extend(j12_02)
        if (m+1) % 50 == 0:
            elapsed = time() - t_meas
            rate = (m+1) / elapsed
            remaining = (n_meas - m - 1) / rate
            print(f"    Meas {m+1}/{n_meas} ({elapsed:.0f}s, ~{remaining:.0f}s remaining)")

    su2_j12_01 = np.array(su2_j12_01_all)
    su2_j1_02 = np.array(su2_j1_02_all)
    su2_j12_02 = np.array(su2_j12_02_all)
    dt_su2 = time() - t0
    print(f"  SU(2) done: {len(su2_j12_01)} samples in {dt_su2:.0f}s")

    # ── U(1) ──
    print(f"\n  U(1): thermalizing ({n_therm} sweeps)...")
    t0 = time()
    u1 = U1Lattice(L)
    for s in range(n_therm):
        u1.sweep(beta_u1)

    print(f"  U(1): measuring ({n_meas} configs)...")
    u1_e1_01_all = []
    u1_e3_02_all = []
    u1_e1_02_all = []
    t_meas = time()
    for m in range(n_meas):
        for _ in range(n_decorr):
            u1.sweep(beta_u1)
        e1_01, e3_02, e1_02 = u1.measure()
        u1_e1_01_all.extend(e1_01)
        u1_e3_02_all.extend(e3_02)
        u1_e1_02_all.extend(e1_02)
        if (m+1) % 50 == 0:
            elapsed = time() - t_meas
            rate = (m+1) / elapsed
            remaining = (n_meas - m - 1) / rate
            print(f"    Meas {m+1}/{n_meas} ({elapsed:.0f}s, ~{remaining:.0f}s remaining)")

    u1_e1_01 = np.array(u1_e1_01_all)
    u1_e3_02 = np.array(u1_e3_02_all)
    u1_e1_02 = np.array(u1_e1_02_all)
    dt_u1 = time() - t0
    print(f"  U(1) done: {len(u1_e1_01)} samples in {dt_u1:.0f}s")

    # ── Compute deltas ──
    n_trials = 30
    N_target = min(50000, len(su2_j12_01), len(u1_e1_01))

    results = {}
    for name, o1, o2 in [
        ("SU2_cross",  su2_j12_01, su2_j1_02),
        ("SU2_same",   su2_j12_01, su2_j12_02),
        ("U1_coprime", u1_e1_01,   u1_e3_02),
        ("U1_same",    u1_e1_01,   u1_e1_02),
    ]:
        ds = [compute_delta(o1, o2, N_target) for _ in range(n_trials)]
        m, s = np.mean(ds), np.std(ds)
        results[name] = (m, s)

    # Also full dataset
    for name, o1, o2 in [
        ("SU2_cross_full",  su2_j12_01, su2_j1_02),
        ("SU2_same_full",   su2_j12_01, su2_j12_02),
        ("U1_coprime_full", u1_e1_01,   u1_e3_02),
        ("U1_same_full",    u1_e1_01,   u1_e1_02),
    ]:
        d = compute_delta(o1, o2)
        results[name] = (d, 0.0)

    return results, dt_su2 + dt_u1


def main():
    np.random.seed(42)
    t_global = time()

    print("=" * 70)
    print("ABELIAN DISCRIMINATION SCALING STUDY")
    print("=" * 70)
    print()
    print("Does U(1) coprime delta converge to 1/3 while")
    print("SU(2) cross-rep delta stays above 1/3?")
    print()

    # Parameters per lattice size
    configs = [
        # (L, beta_su2, beta_u1, n_therm, n_meas, n_decorr)
        (4,  2.2, 1.0, 50,  400, 3),
        (6,  2.2, 1.0, 60,  150, 3),
        (8,  2.2, 1.0, 80,  60,  3),
    ]

    all_results = {}

    for L, beta_su2, beta_u1, n_therm, n_meas, n_decorr in configs:
        results, dt = run_lattice_size(L, beta_su2, beta_u1, n_therm, n_meas, n_decorr)
        all_results[L] = results

        # Print intermediate results
        print(f"\n  --- L={L} Results (N=50K subsampled, 30 trials) ---")
        for key in ["SU2_cross", "SU2_same", "U1_coprime", "U1_same"]:
            m, s = results[key]
            print(f"    {key:15s}: delta = {m:.4f} +/- {s:.4f}")
        print(f"    Full dataset:")
        for key in ["SU2_cross_full", "SU2_same_full", "U1_coprime_full", "U1_same_full"]:
            m, _ = results[key]
            label = key.replace("_full", "")
            print(f"    {label:15s}: delta = {m:.4f}")
        sys.stdout.flush()

    # ── Final summary ──
    print()
    print("=" * 70)
    print("SCALING SUMMARY")
    print("=" * 70)
    print()
    print(f"{'L':>3s}  {'SU2 cross':>12s}  {'SU2 same':>12s}  {'U1 coprime':>12s}  {'U1 same':>12s}  {'Gap':>8s}")
    print("-" * 70)

    for L in sorted(all_results.keys()):
        r = all_results[L]
        su2c = r["SU2_cross_full"][0]
        su2s = r["SU2_same_full"][0]
        u1c = r["U1_coprime_full"][0]
        u1s = r["U1_same_full"][0]
        gap = su2c - u1c
        print(f"{L:>3d}  {su2c:>12.4f}  {su2s:>12.4f}  {u1c:>12.4f}  {u1s:>12.4f}  {gap:>+8.4f}")

    print()
    print(f"1/H = {1/H:.4f}")
    print()

    # Trend analysis
    Ls = sorted(all_results.keys())
    su2_cross_vals = [all_results[L]["SU2_cross_full"][0] for L in Ls]
    u1_coprime_vals = [all_results[L]["U1_coprime_full"][0] for L in Ls]
    gaps = [s - u for s, u in zip(su2_cross_vals, u1_coprime_vals)]

    print("TREND ANALYSIS:")
    print(f"  SU(2) cross-rep trend: {su2_cross_vals}")
    if su2_cross_vals[-1] > 1/3 + 0.002:
        print(f"    -> Stays above 1/3 (last value {su2_cross_vals[-1]:.4f})")
    else:
        print(f"    -> Converging toward 1/3 (last value {su2_cross_vals[-1]:.4f})")

    print(f"  U(1) coprime trend: {u1_coprime_vals}")
    if u1_coprime_vals[-1] < u1_coprime_vals[0]:
        print(f"    -> Dropping with L (consistent with -> 1/3)")
    else:
        print(f"    -> Not dropping (may plateau above 1/3)")

    print(f"  Gap trend: {gaps}")
    if gaps[-1] > gaps[0]:
        print(f"    -> Gap INCREASING with L: discrimination strengthens")
    elif gaps[-1] > 0:
        print(f"    -> Gap positive but {'stable' if abs(gaps[-1]-gaps[0]) < 0.002 else 'shrinking'}")
    else:
        print(f"    -> Gap closed or negative: discrimination fails")

    print()
    dt_total = time() - t_global
    print(f"Total runtime: {dt_total:.0f}s ({dt_total/60:.1f} min)")
    print()

    if gaps[-1] > 0.003:
        print("VERDICT: Abelian discrimination is REAL and persists with lattice size.")
    elif gaps[-1] > 0:
        print("VERDICT: Abelian discrimination is present but small. Needs larger lattices.")
    else:
        print("VERDICT: No discrimination detected. The observable does not separate the theories.")


if __name__ == "__main__":
    main()
