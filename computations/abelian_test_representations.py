"""
Abelian Discrimination Test — Representation Version

The CORRECT discriminator: same links, different representations.

SU(2): Wilson loops in representation j=1/2 and j=1 at the SAME plaquette.
  The Clebsch-Gordan decomposition 1/2 ⊗ 1 = 1/2 ⊕ 3/2 has cross terms.
  The conditional distribution C should have δ > 1/H (structural correlation).

U(1): Wilson loops with charge e=1 and e=2 at the SAME plaquette.
  But gcd(1,2) = 2, not coprime. Try e=1 and e=3 (gcd=1, coprime).
  For coprime charges: the Peter-Weyl orthogonality says:
    ∫ χ_e1(U) χ_e2(U)* dU = δ_{e1,e2}
  So different charges are orthogonal under the Haar measure.
  BUT on a thermalized lattice, the Boltzmann weight breaks Haar — the
  action couples ALL representations through the plaquette action.

Actually, let me think more carefully.

For U(1): the plaquette is exp(i*phase). The charge-e Wilson loop is
Re(exp(i*e*phase)) = cos(e*phase). Two charges e1, e2 at the same
plaquette give cos(e1*θ) and cos(e2*θ) where θ is the plaquette angle.
These are ALWAYS correlated (deterministic functions of the same θ)
regardless of whether the gauge group is abelian or not.

For SU(2): the j=1/2 trace is Re(Tr(U))/2. The j=1 trace is
Re(Tr_adj(U))/3 = (|Tr(U)|^2 - 1)/3. Again, deterministic functions
of the same matrix U.

So representation-based observables at the SAME plaquette are always
correlated — they're functions of the same underlying variable. This
is NOT what the framework's abelian discrimination is about.

The ACTUAL discriminator the paper describes (Theorem thm:universal):
  Non-abelian: the tensor product R^(p) ⊗ R^(q) has non-trivial
  Clebsch-Gordan decomposition. Integration over SHARED LINKS produces
  cross terms that don't vanish.

  Abelian: all irreps are 1D. Tensor product is just multiplication
  of phases. No cross terms.

The test should be: two Wilson loops that SHARE some links but traverse
DIFFERENT paths. The shared-link integration creates the correlation,
and the Clebsch-Gordan structure determines whether it persists.

Concrete test:
  Loop A: plaquette in (0,1) plane at site x — uses links U_0(x), U_1(x+0), U_0(x+1)†, U_1(x)†
  Loop B: plaquette in (0,2) plane at site x — uses links U_0(x), U_2(x+0), U_0(x+2)†, U_2(x)†

  They SHARE link U_0(x). Integration over this shared link:
    SU(2): ∫ Tr(U*A) Tr(U*B) dU = Tr(A†B)/2 ≠ Tr(A†)Tr(B)/4 in general
    U(1):  ∫ exp(iθ)*a * exp(iθ)*b dθ = a*b * ∫ exp(2iθ) dθ = 0 (if coprime charges)

Wait — for U(1) with charge 1 in both loops:
  ∫ exp(iθ) * exp(iθ) dθ/2π = ∫ exp(2iθ) dθ/2π = 0

But both loops use U_0(x) in the SAME direction (not conjugate), so the
contribution from the shared link is exp(iθ) in both cases — same charge.
The product is exp(2iθ), which integrates to zero under Haar.

For SU(2) with j=1/2 in both:
  ∫ U_ij * U_kl dU = (1/2) δ_il δ_jk
  This gives a non-zero contribution to the plaquette-plaquette correlator.

THIS is the Peter-Weyl discrimination. Same charge, shared link, but the
SU(2) Haar integral gives non-zero cross terms while U(1) gives zero.

But on a THERMALIZED lattice, the Haar integration is weighted by exp(-βS).
The Boltzmann weight can restore correlations even for U(1). The question
is whether the connected correlator persists as the lattice thermalizes.

Let me just measure it directly: plaquettes in (0,1) and (0,2) at the
SAME SITE (sharing the direction-0 link) on a thermalized lattice, and
check δ scaling with N.

THIS is what gap7_abelian.py was already doing, and the result was:
  U(1) correlation: 0.628, connected: 0.5σ
  SU(2) correlation: 0.156, connected: 0.0σ

The U(1) was MORE correlated, which confused me. But that's because the
lattice-averaged plaquettes have autocorrelation from the Markov chain.
Let me use SITE-LEVEL plaquettes (not lattice averages) to get real
statistics and do the δ vs N properly.
"""

import numpy as np
from time import time

H = 3

# ═══════════════════════════════════════════════════════════════
#  Lattices (reuse from previous scripts)
# ═══════════════════════════════════════════════════════════════

I2 = np.eye(2, dtype=complex)
sigma = [np.array([[0,1],[1,0]],dtype=complex),
         np.array([[0,-1j],[1j,0]],dtype=complex),
         np.array([[1,0],[0,-1]],dtype=complex)]

def random_su2():
    a = np.random.randn(4); a /= np.linalg.norm(a)
    return a[0]*I2 + 1j*(a[1]*sigma[0]+a[2]*sigma[1]+a[3]*sigma[2])

def random_su2_near_id(eps=0.3):
    a = np.array([1.0,0,0,0]) + eps*np.random.randn(4)
    a /= np.linalg.norm(a)
    return a[0]*I2 + 1j*(a[1]*sigma[0]+a[2]*sigma[1]+a[3]*sigma[2])


def u1_sweep(links, beta, L):
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


def su2_sweep(links, beta, L):
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
                        old_U = links[pos+(mu,)].copy()
                        old_a = np.real(np.trace(old_U @ S))
                        for _ in range(4):
                            trial = random_su2_near_id(0.3) @ old_U
                            trial_a = np.real(np.trace(trial @ S))
                            dS = beta/2*(trial_a - old_a)
                            if dS > 0 or np.random.random() < np.exp(dS):
                                links[pos+(mu,)] = trial
                                old_U = trial; old_a = trial_a


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
    np.random.seed(42)
    L = 4
    n_therm = 50
    n_decorr = 3
    n_meas = 500

    print("=" * 70)
    print("ABELIAN DISCRIMINATION: SHARED-LINK PLAQUETTES")
    print("=" * 70)
    print(f"Lattice: {L}^4, n_meas={n_meas}, n_decorr={n_decorr}")
    print()
    print("Observable: plaquettes in (0,1) and (0,2) planes at SAME site")
    print("These share the direction-0 link. The shared-link integration")
    print("is where SU(2) Clebsch-Gordan cross terms appear.")
    print()

    # ── U(1) ──
    print("--- U(1) (β=1.0) ---")
    t0 = time()
    u1_links = np.random.uniform(0, 2*np.pi, (L,L,L,L,4))
    for s in range(n_therm):
        u1_sweep(u1_links, 1.0, L)

    u1_P01, u1_P02 = [], []
    for m in range(n_meas):
        for _ in range(n_decorr):
            u1_sweep(u1_links, 1.0, L)
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        pos = [x,y,z,t]
                        pm = pos.copy(); pm[0]=(pm[0]+1)%L
                        p1 = pos.copy(); p1[1]=(p1[1]+1)%L
                        p2 = pos.copy(); p2[2]=(p2[2]+1)%L
                        # (0,1) plaquette
                        ph01 = u1_links[x,y,z,t,0]+u1_links[tuple(pm)][1]-u1_links[tuple(p1)][0]-u1_links[x,y,z,t,1]
                        u1_P01.append(np.cos(ph01))
                        # (0,2) plaquette — shares link U_0(x)
                        ph02 = u1_links[x,y,z,t,0]+u1_links[tuple(pm)][2]-u1_links[tuple(p2)][0]-u1_links[x,y,z,t,2]
                        u1_P02.append(np.cos(ph02))
        if (m+1) % 100 == 0:
            print(f"  {m+1}/{n_meas} ({time()-t0:.0f}s)")

    u1_P01 = np.array(u1_P01); u1_P02 = np.array(u1_P02)
    print(f"  Done: {len(u1_P01)} samples in {time()-t0:.0f}s")
    print()

    # ── SU(2) ──
    print("--- SU(2) (β=2.2) ---")
    t0 = time()
    su2_links = np.zeros((L,L,L,L,4,2,2), dtype=complex)
    for idx in np.ndindex(L,L,L,L):
        for mu in range(4):
            su2_links[idx+(mu,)] = random_su2()
    for s in range(n_therm):
        su2_sweep(su2_links, 2.2, L)
        if s % 10 == 0:
            print(f"  Therm sweep {s} ({time()-t0:.0f}s)")

    su2_P01, su2_P02 = [], []
    for m in range(n_meas):
        for _ in range(n_decorr):
            su2_sweep(su2_links, 2.2, L)
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    for t in range(L):
                        pos = (x,y,z,t)
                        pm = list(pos); pm[0]=(pm[0]+1)%L; pm=tuple(pm)
                        p1 = list(pos); p1[1]=(p1[1]+1)%L; p1=tuple(p1)
                        p2 = list(pos); p2[2]=(p2[2]+1)%L; p2=tuple(p2)
                        # (0,1) plaquette
                        plaq01 = (su2_links[pos+(0,)] @ su2_links[pm+(1,)] @
                                  su2_links[p1+(0,)].conj().T @ su2_links[pos+(1,)].conj().T)
                        su2_P01.append(np.real(np.trace(plaq01))/2)
                        # (0,2) plaquette — shares link U_0(x)
                        plaq02 = (su2_links[pos+(0,)] @ su2_links[pm+(2,)] @
                                  su2_links[p2+(0,)].conj().T @ su2_links[pos+(2,)].conj().T)
                        su2_P02.append(np.real(np.trace(plaq02))/2)
        if (m+1) % 100 == 0:
            print(f"  {m+1}/{n_meas} ({time()-t0:.0f}s)")

    su2_P01 = np.array(su2_P01); su2_P02 = np.array(su2_P02)
    print(f"  Done: {len(su2_P01)} samples in {time()-t0:.0f}s")
    print()

    # ── Results ──
    print("=" * 70)
    print("δ vs N — SHARED-LINK PLAQUETTES")
    print("=" * 70)
    print()

    Ns = [500, 1000, 2000, 5000, 10000, 50000, 100000]

    print(f"{'N':>8s}  {'U(1) δ':>14s}  {'SU(2) δ':>14s}  {'Diff':>10s}  {'Sig':>6s}")
    print("-" * 65)

    for N in Ns:
        if N > min(len(u1_P01), len(su2_P01)):
            continue
        u1_ds, su2_ds = [], []
        for _ in range(30):
            u1_ds.append(compute_delta(u1_P01, u1_P02, N))
            su2_ds.append(compute_delta(su2_P01, su2_P02, N))

        u1_m, u1_s = np.mean(u1_ds), np.std(u1_ds)
        su2_m, su2_s = np.mean(su2_ds), np.std(su2_ds)
        diff = su2_m - u1_m
        sig = diff / np.sqrt(su2_s**2 + u1_s**2) if (su2_s + u1_s) > 0 else 0

        print(f"{N:>8d}  {u1_m:>8.4f}±{u1_s:.4f}  {su2_m:>8.4f}±{su2_s:.4f}  {diff:>+10.4f}  {sig:>5.1f}σ")

    print()
    print("1/H =", 1/H)
    print()

    # Full dataset
    N_full = min(len(u1_P01), len(su2_P01))
    d_u1_full = compute_delta(u1_P01[:N_full], u1_P02[:N_full])
    d_su2_full = compute_delta(su2_P01[:N_full], su2_P02[:N_full])
    print(f"Full dataset (N={N_full}):")
    print(f"  U(1)  δ = {d_u1_full:.6f}")
    print(f"  SU(2) δ = {d_su2_full:.6f}")
    print(f"  Diff    = {d_su2_full - d_u1_full:+.6f}")
    print()

    # Connected correlator
    corr_u1 = np.corrcoef(u1_P01[:N_full], u1_P02[:N_full])[0,1]
    corr_su2 = np.corrcoef(su2_P01[:N_full], su2_P02[:N_full])[0,1]
    print(f"Pearson correlation:")
    print(f"  U(1)  r = {corr_u1:.6f}")
    print(f"  SU(2) r = {corr_su2:.6f}")
    print()

    if d_su2_full > d_u1_full + 0.001:
        print("SU(2) has higher δ than U(1) at shared links.")
        print("This is the signature of non-abelian Clebsch-Gordan coupling.")
    elif d_u1_full > d_su2_full + 0.001:
        print("U(1) has higher δ than SU(2) — unexpected.")
        print("Both share links; the coupling mechanism differs.")
    else:
        print("δ values are indistinguishable within precision.")
        print("The shared-link correlation is similar for both theories.")


if __name__ == "__main__":
    main()
