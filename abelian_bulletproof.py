"""
Abelian Discrimination — Bulletproof Run

L = 4, 6, 8, 10, 12. Maximum statistics at each size.
The U(1) coprime δ must reach 1/3 and SU(2) cross-rep must stay above.
No ambiguity. No "directionally correct." Undeniable.
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
                        plaq01 = self.links[pos+(0,)] @ self.links[tuple(pm0)+(1,)] @ self.links[tuple(p1)+(0,)].conj().T @ self.links[pos+(1,)].conj().T
                        plaq02 = self.links[pos+(0,)] @ self.links[tuple(pm0)+(2,)] @ self.links[tuple(p2)+(0,)].conj().T @ self.links[pos+(2,)].conj().T
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


def run_one(L, beta_su2, beta_u1, n_therm, n_meas, n_decorr):
    n_sites = L**4
    t0 = time()

    print(f"\n{'='*70}")
    print(f"  L = {L}  ({n_sites} sites, {n_meas} configs)")
    print(f"{'='*70}")
    sys.stdout.flush()

    # SU(2)
    print(f"  SU(2) thermalizing...", end="", flush=True)
    su2 = SU2Lattice(L)
    for s in range(n_therm):
        su2.sweep(beta_su2, n_hits=4)
    print(f" done ({time()-t0:.0f}s)")

    t_m = time()
    print(f"  SU(2) measuring...", flush=True)
    su2_a, su2_b, su2_c = [], [], []
    for m in range(n_meas):
        for _ in range(n_decorr):
            su2.sweep(beta_su2, n_hits=3)
        a, b, c = su2.measure()
        su2_a.extend(a); su2_b.extend(b); su2_c.extend(c)
        if (m+1) % max(1, n_meas//5) == 0:
            el = time()-t_m; rm = el/(m+1)*(n_meas-m-1)
            print(f"    {m+1}/{n_meas} ({el:.0f}s, ~{rm:.0f}s left)")
            sys.stdout.flush()
    su2_a = np.array(su2_a); su2_b = np.array(su2_b); su2_c = np.array(su2_c)

    # U(1)
    t_u = time()
    print(f"  U(1) thermalizing...", end="", flush=True)
    u1 = U1Lattice(L)
    for s in range(n_therm):
        u1.sweep(beta_u1)
    print(f" done ({time()-t_u:.0f}s)")

    t_m = time()
    print(f"  U(1) measuring...", flush=True)
    u1_a, u1_b, u1_c = [], [], []
    for m in range(n_meas):
        for _ in range(n_decorr):
            u1.sweep(beta_u1)
        a, b, c = u1.measure()
        u1_a.extend(a); u1_b.extend(b); u1_c.extend(c)
        if (m+1) % max(1, n_meas//5) == 0:
            el = time()-t_m; rm = el/(m+1)*(n_meas-m-1)
            print(f"    {m+1}/{n_meas} ({el:.0f}s, ~{rm:.0f}s left)")
            sys.stdout.flush()
    u1_a = np.array(u1_a); u1_b = np.array(u1_b); u1_c = np.array(u1_c)

    # Compute deltas — full dataset
    d_su2_cross = compute_delta(su2_a, su2_b)
    d_su2_same = compute_delta(su2_a, su2_c)
    d_u1_coprime = compute_delta(u1_a, u1_b)
    d_u1_same = compute_delta(u1_a, u1_c)

    # Error bars via bootstrap
    n_boot = 50
    N_sub = min(50000, len(su2_a), len(u1_a))
    su2_cross_boots = [compute_delta(su2_a, su2_b, N_sub) for _ in range(n_boot)]
    u1_coprime_boots = [compute_delta(u1_a, u1_b, N_sub) for _ in range(n_boot)]

    su2_err = np.std(su2_cross_boots)
    u1_err = np.std(u1_coprime_boots)
    gap = d_su2_cross - d_u1_coprime
    gap_err = np.sqrt(su2_err**2 + u1_err**2)
    gap_sigma = gap / gap_err if gap_err > 0 else 0

    dt = time() - t0
    print(f"\n  L={L} RESULT ({dt:.0f}s, {len(su2_a)} samples):")
    print(f"    SU(2) cross:  {d_su2_cross:.4f} +/- {su2_err:.4f}")
    print(f"    SU(2) same:   {d_su2_same:.4f}")
    print(f"    U(1) coprime: {d_u1_coprime:.4f} +/- {u1_err:.4f}")
    print(f"    U(1) same:    {d_u1_same:.4f}")
    print(f"    GAP:          {gap:+.4f} +/- {gap_err:.4f} ({gap_sigma:.1f} sigma)")
    sys.stdout.flush()

    return {
        "su2_cross": d_su2_cross, "su2_cross_err": su2_err,
        "su2_same": d_su2_same,
        "u1_coprime": d_u1_coprime, "u1_coprime_err": u1_err,
        "u1_same": d_u1_same,
        "gap": gap, "gap_err": gap_err, "gap_sigma": gap_sigma,
        "samples": len(su2_a), "time": dt,
    }


def main():
    np.random.seed(42)
    t_global = time()

    print("=" * 70)
    print("ABELIAN DISCRIMINATION — BULLETPROOF RUN")
    print("=" * 70)
    print("L = 4, 6, 8, 10, 12")
    print("Maximum statistics. No ambiguity.")
    print()
    sys.stdout.flush()

    # Aggressive configs: more measurements at small L, fewer at large L
    # but more sites per measurement at large L, so total samples scale
    configs = [
        # (L, beta_su2, beta_u1, n_therm, n_meas, n_decorr)
        (4,   2.2, 1.0, 60,  500, 3),    # 128K samples
        (6,   2.2, 1.0, 80,  200, 3),    # 259K samples
        (8,   2.2, 1.0, 100, 80,  3),    # 328K samples
        (10,  2.2, 1.0, 100, 40,  3),    # 400K samples
        (12,  2.2, 1.0, 100, 20,  3),    # 414K samples
    ]

    results = {}
    for L, b2, b1, nt, nm, nd in configs:
        results[L] = run_one(L, b2, b1, nt, nm, nd)

    # Final table
    print()
    print("=" * 70)
    print("FINAL SCALING TABLE")
    print("=" * 70)
    print()
    print(f"{'L':>3s}  {'SU2 cross':>14s}  {'U1 coprime':>14s}  {'Gap':>14s}  {'Sig':>6s}  {'Samples':>8s}")
    print("-" * 70)
    for L in sorted(results.keys()):
        r = results[L]
        print(f"{L:>3d}  {r['su2_cross']:>8.4f}+/-{r['su2_cross_err']:.4f}"
              f"  {r['u1_coprime']:>8.4f}+/-{r['u1_coprime_err']:.4f}"
              f"  {r['gap']:>+8.4f}+/-{r['gap_err']:.4f}"
              f"  {r['gap_sigma']:>5.1f}s"
              f"  {r['samples']:>8d}")

    print()
    print(f"1/H = {1/H:.4f}")
    print()

    # Verdict
    Ls = sorted(results.keys())
    su2_vals = [results[L]["su2_cross"] for L in Ls]
    u1_vals = [results[L]["u1_coprime"] for L in Ls]
    gaps = [results[L]["gap"] for L in Ls]
    sigs = [results[L]["gap_sigma"] for L in Ls]

    all_positive = all(g > 0 for g in gaps)
    all_significant = all(s > 2 for s in sigs)
    u1_dropping = u1_vals[-1] < u1_vals[0]
    su2_stable = abs(su2_vals[-1] - su2_vals[0]) < 0.005
    u1_near_third = abs(u1_vals[-1] - 1/3) < 0.005

    print("ANALYSIS:")
    print(f"  All gaps positive: {all_positive}")
    print(f"  All gaps > 2 sigma: {all_significant}")
    print(f"  U(1) coprime dropping: {u1_dropping} ({u1_vals[0]:.4f} -> {u1_vals[-1]:.4f})")
    print(f"  SU(2) cross stable: {su2_stable} ({su2_vals[0]:.4f} -> {su2_vals[-1]:.4f})")
    print(f"  U(1) near 1/3: {u1_near_third} (last = {u1_vals[-1]:.4f})")
    print()

    if all_positive and u1_dropping and su2_stable:
        if all_significant:
            print("VERDICT: UNDENIABLE. The abelian discrimination is confirmed")
            print("at every lattice size with statistical significance.")
        else:
            print("VERDICT: CONFIRMED but some sizes lack statistical power.")
            print("The trend is clear; larger statistics would sharpen it.")
    elif all_positive:
        print("VERDICT: CONSISTENT. All gaps positive, trend present.")
    else:
        print("VERDICT: INCONCLUSIVE. Some gaps negative.")

    dt = time() - t_global
    print(f"\nTotal runtime: {dt:.0f}s ({dt/60:.1f} min, {dt/3600:.1f} hours)")


if __name__ == "__main__":
    main()
