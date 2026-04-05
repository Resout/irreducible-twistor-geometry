"""
Representation Mixing Test

THE discriminator: different representations at the same plaquette.

SU(2): j=1/2 trace = Re(Tr(U))/2 ranges in [-1, 1]
       j=1 trace = (|Tr(U)|^2 - 1)/3 ranges in [-1/3, 1]
       (this is the adjoint/spin-1 character)

       These are BOTH functions of the same matrix U, but they probe
       DIFFERENT representations. The Clebsch-Gordan decomposition
       1/2 ⊗ 1 = 1/2 ⊕ 3/2 means these representations MIX
       under the group action. The DS framework should detect this
       mixing as δ > 1/H.

U(1): charge 1 = Re(exp(iθ)) = cos(θ)
      charge 2 = Re(exp(2iθ)) = cos(2θ) = 2cos²(θ) - 1

      These are also functions of the same variable θ, but they're
      in DIFFERENT representations (charges). For U(1), different
      charges are orthogonal under Haar:
        ∫ exp(inθ) exp(-imθ) dθ/2π = δ_{nm}

      BUT: cos(θ) and cos(2θ) = 2cos²(θ)-1 are NOT independent!
      They're deterministically related. This isn't a representation
      orthogonality test — it's a functional dependence test.

      The CORRECT U(1) coprime test: two SEPARATE plaquettes, one
      measured with charge 1 and one with charge 3 (coprime).
      Same shared link, different charges.
      ∫ exp(iθ) * exp(3iθ) dθ/2π = ∫ exp(4iθ) dθ/2π = 0
      vs
      SU(2): ∫ Tr_{1/2}(U*A) * Tr_1(U*B) dU ≠ 0 generically
      (because 1/2 ⊗ 1 decomposes non-trivially)

Actually the clearest test:

Plaquettes in (0,1) and (0,2) sharing link U_0.
  SU(2): measure Tr_{1/2}(plaq_01) and Tr_1(plaq_02)
         Different representations at plaquettes sharing a link.
         The shared-link Haar integral mixes representations.

  U(1):  measure cos(1*plaq_01) and cos(3*plaq_02)
         Different charges at plaquettes sharing a link.
         Shared link phase: charges 1 and 3 give total charge 4.
         Haar: ∫ exp(iθ) * exp(3iθ) dθ = 0. Should decorrelate.

THIS is the test. Same shared link, different representations.
SU(2) mixes, U(1) doesn't (for coprime total charge).
"""

import numpy as np
from time import time

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

L = 4
np.random.seed(42)

# ── SU(2) ──
su2 = np.zeros((L,L,L,L,4,2,2), dtype=complex)
for idx in np.ndindex(L,L,L,L):
    for mu in range(4):
        su2[idx+(mu,)] = random_su2()

def su2_sweep(links, beta):
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
                            S += links[pm+(nu,)] @ links[pn+(mu,)].conj().T @ links[pos+(nu,)].conj().T
                            S += links[pmun+(nu,)].conj().T @ links[pmn+(mu,)].conj().T @ links[pmn+(nu,)]
                        old = links[pos+(mu,)].copy()
                        oa = np.real(np.trace(old @ S))
                        for _ in range(4):
                            trial = near_id() @ old
                            ta = np.real(np.trace(trial @ S))
                            if ta > oa or np.random.random() < np.exp(beta/2*(ta-oa)):
                                links[pos+(mu,)] = trial; old = trial; oa = ta

# ── U(1) ──
u1 = np.random.uniform(0, 2*np.pi, (L,L,L,L,4))

def u1_sweep(links, beta):
    for x in range(L):
        for y in range(L):
            for z in range(L):
                for t in range(L):
                    for mu in range(4):
                        old = links[x,y,z,t,mu]
                        new = old + np.random.uniform(-0.5, 0.5)
                        dS = 0.0
                        for nu in range(4):
                            if nu == mu: continue
                            pos = [x,y,z,t]
                            pm = pos.copy(); pm[mu]=(pm[mu]+1)%L
                            pn = pos.copy(); pn[nu]=(pn[nu]+1)%L
                            pmn = pos.copy(); pmn[nu]=(pmn[nu]-1)%L
                            pmun = pos.copy(); pmun[mu]=(pmun[mu]+1)%L; pmun[nu]=(pmun[nu]-1)%L
                            p_old = old + links[tuple(pm)][nu] - links[tuple(pn)][mu] - links[x,y,z,t,nu]
                            p_new = new + links[tuple(pm)][nu] - links[tuple(pn)][mu] - links[x,y,z,t,nu]
                            dS -= beta*(np.cos(p_new) - np.cos(p_old))
                            pb_old = links[tuple(pmn)][nu]+old-links[tuple(pmun)][nu]-links[tuple(pmn)][mu]
                            pb_new = links[tuple(pmn)][nu]+new-links[tuple(pmun)][nu]-links[tuple(pmn)][mu]
                            dS -= beta*(np.cos(pb_new) - np.cos(pb_old))
                        if dS < 0 or np.random.rand() < np.exp(-dS):
                            links[x,y,z,t,mu] = new


def compute_delta(o1, o2, N=None):
    if N and N < len(o1):
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


def main():
    n_therm = 50
    n_meas = 400
    n_decorr = 3

    print("=" * 70)
    print("REPRESENTATION MIXING TEST")
    print("=" * 70)
    print()
    print("SU(2): j=1/2 plaquette in (0,1) vs j=1 plaquette in (0,2)")
    print("       sharing direction-0 link. CG: 1/2 x 1 = 1/2 + 3/2")
    print()
    print("U(1):  charge-1 plaquette in (0,1) vs charge-3 in (0,2)")
    print("       sharing direction-0 link. Total charge 4 -> Haar kills it")
    print()

    # Thermalize
    print("Thermalizing...")
    t0 = time()
    for s in range(n_therm):
        su2_sweep(su2, 2.2)
        u1_sweep(u1, 1.0)
    print(f"  Done in {time()-t0:.0f}s")

    # Measure
    su2_j12 = []  # j=1/2 trace of (0,1) plaquette
    su2_j1 = []   # j=1 trace of (0,2) plaquette
    u1_e1 = []    # charge 1 of (0,1) plaquette
    u1_e3 = []    # charge 3 of (0,2) plaquette
    # Also same-rep controls
    su2_j12_02 = []  # j=1/2 trace of (0,2) plaquette (control: same rep)
    u1_e1_02 = []    # charge 1 of (0,2) plaquette (control: same charge)

    print(f"Measuring ({n_meas} configs)...")
    t0 = time()
    for m in range(n_meas):
        for _ in range(n_decorr):
            su2_sweep(su2, 2.2)
            u1_sweep(u1, 1.0)

        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        pos = (x,y,z,t)
                        pm0 = list(pos); pm0[0]=(pm0[0]+1)%L; pm0=tuple(pm0)
                        p1 = list(pos); p1[1]=(p1[1]+1)%L; p1=tuple(p1)
                        p2 = list(pos); p2[2]=(p2[2]+1)%L; p2=tuple(p2)

                        # SU(2) plaquettes
                        plaq01 = su2[pos+(0,)] @ su2[pm0+(1,)] @ su2[p1+(0,)].conj().T @ su2[pos+(1,)].conj().T
                        plaq02 = su2[pos+(0,)] @ su2[pm0+(2,)] @ su2[p2+(0,)].conj().T @ su2[pos+(2,)].conj().T

                        # j=1/2: Tr(U)/2
                        su2_j12.append(np.real(np.trace(plaq01)) / 2)
                        su2_j12_02.append(np.real(np.trace(plaq02)) / 2)

                        # j=1 (adjoint): (|Tr(U)|^2 - 1) / 3
                        tr02 = np.trace(plaq02)
                        su2_j1.append((np.abs(tr02)**2 - 1) / 3)

                        # U(1) plaquettes
                        ph01 = u1[x,y,z,t,0] + u1[pm0][1] - u1[p1][0] - u1[x,y,z,t,1]
                        ph02 = u1[x,y,z,t,0] + u1[pm0][2] - u1[p2][0] - u1[x,y,z,t,2]

                        # charge 1
                        u1_e1.append(np.cos(ph01))
                        u1_e1_02.append(np.cos(ph02))

                        # charge 3 on the (0,2) plaquette
                        u1_e3.append(np.cos(3 * ph02))

        if (m+1) % 100 == 0:
            print(f"  {m+1}/{n_meas} ({time()-t0:.0f}s)")

    su2_j12 = np.array(su2_j12)
    su2_j1 = np.array(su2_j1)
    su2_j12_02 = np.array(su2_j12_02)
    u1_e1 = np.array(u1_e1)
    u1_e3 = np.array(u1_e3)
    u1_e1_02 = np.array(u1_e1_02)

    print(f"  Done: {len(su2_j12)} samples in {time()-t0:.0f}s")
    print()

    # Results
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    print()

    Ns = [500, 1000, 2000, 5000, 10000, 50000, 100000]

    # Test A: SU(2) j=1/2 vs j=1 (different reps, shared link)
    print("TEST A: SU(2) j=1/2 (plaq 01) vs j=1 (plaq 02) — CROSS-REP")
    for N in Ns:
        if N > len(su2_j12): continue
        ds = [compute_delta(su2_j12, su2_j1, N) for _ in range(20)]
        print(f"  N={N:>7d}: delta = {np.mean(ds):.4f} +/- {np.std(ds):.4f}")
    print()

    # Test B: SU(2) j=1/2 vs j=1/2 (same rep, shared link — control)
    print("TEST B: SU(2) j=1/2 (plaq 01) vs j=1/2 (plaq 02) — SAME REP (control)")
    for N in Ns:
        if N > len(su2_j12): continue
        ds = [compute_delta(su2_j12, su2_j12_02, N) for _ in range(20)]
        print(f"  N={N:>7d}: delta = {np.mean(ds):.4f} +/- {np.std(ds):.4f}")
    print()

    # Test C: U(1) charge-1 vs charge-3 (coprime charges, shared link)
    print("TEST C: U(1) charge 1 (plaq 01) vs charge 3 (plaq 02) — COPRIME")
    for N in Ns:
        if N > len(u1_e1): continue
        ds = [compute_delta(u1_e1, u1_e3, N) for _ in range(20)]
        print(f"  N={N:>7d}: delta = {np.mean(ds):.4f} +/- {np.std(ds):.4f}")
    print()

    # Test D: U(1) charge-1 vs charge-1 (same charge, shared link — control)
    print("TEST D: U(1) charge 1 (plaq 01) vs charge 1 (plaq 02) — SAME CHARGE (control)")
    for N in Ns:
        if N > len(u1_e1): continue
        ds = [compute_delta(u1_e1, u1_e1_02, N) for _ in range(20)]
        print(f"  N={N:>7d}: delta = {np.mean(ds):.4f} +/- {np.std(ds):.4f}")
    print()

    # Summary
    d_su2_cross = compute_delta(su2_j12, su2_j1)
    d_su2_same = compute_delta(su2_j12, su2_j12_02)
    d_u1_coprime = compute_delta(u1_e1, u1_e3)
    d_u1_same = compute_delta(u1_e1, u1_e1_02)

    print("=" * 70)
    print("SUMMARY (full dataset)")
    print("=" * 70)
    print(f"  1/H = {1/H:.4f}")
    print()
    print(f"  SU(2) cross-rep (j=1/2 x j=1):  delta = {d_su2_cross:.4f}")
    print(f"  SU(2) same-rep  (j=1/2 x j=1/2): delta = {d_su2_same:.4f}")
    print(f"  U(1)  coprime   (e=1 x e=3):     delta = {d_u1_coprime:.4f}")
    print(f"  U(1)  same      (e=1 x e=1):     delta = {d_u1_same:.4f}")
    print()
    print("Framework predictions:")
    print("  SU(2) cross-rep: delta > 1/H (CG mixing)")
    print("  SU(2) same-rep:  delta > 1/H (shared link)")
    print("  U(1) coprime:    delta -> 1/H (orthogonal reps)")
    print("  U(1) same:       delta > 1/H (shared link)")
    print()

    if d_u1_coprime < d_u1_same - 0.005:
        print("*** U(1) coprime IS lower than U(1) same ***")
        print("    Coprime charges decorrelate — consistent with Peter-Weyl")
    if d_su2_cross > 1/H + 0.005:
        print("*** SU(2) cross-rep IS above 1/H ***")
        print("    Representation mixing detected — CG cross terms")
    if d_su2_cross > d_u1_coprime + 0.005:
        print("*** SU(2) cross-rep > U(1) coprime ***")
        print("    THIS IS THE ABELIAN DISCRIMINATION")


if __name__ == "__main__":
    main()
