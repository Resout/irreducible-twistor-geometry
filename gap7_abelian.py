"""
Gap 7: Abelian discrimination.
Test whether U(1) Wilson loops in different planes factorize (delta = 1)
when they share links, versus SU(2) where they don't.

The key: plaquettes in planes (0,1) and (0,2) share direction-0 links.
For U(1): integration over shared links factorizes (1D irreps).
For SU(2): Clebsch-Gordan cross terms prevent factorization.
"""
import numpy as np

# =====================================================
# PART 1: Compact U(1)
# =====================================================
print("=" * 60)
print("PART 1: Compact U(1)")
print("=" * 60)

L = 6
beta = 1.0
np.random.seed(42)
links_u1 = np.random.uniform(0, 2*np.pi, (L, L, L, L, 4))

def u1_sweep(links, beta, L):
    for idx in np.ndindex(L, L, L, L):
        for mu in range(4):
            old = links[idx][mu]
            new = old + np.random.uniform(-0.5, 0.5)
            dS = 0.0
            for nu in range(4):
                if nu == mu: continue
                pm = list(idx); pm[mu] = (pm[mu]+1)%L
                pn = list(idx); pn[nu] = (pn[nu]+1)%L
                pmn = list(idx); pmn[nu] = (pmn[nu]-1)%L
                pmun = list(idx); pmun[mu]=(pmun[mu]+1)%L; pmun[nu]=(pmun[nu]-1)%L
                p1o = old + links[tuple(pm)][nu] - links[tuple(pn)][mu] - links[idx][nu]
                p1n = new + links[tuple(pm)][nu] - links[tuple(pn)][mu] - links[idx][nu]
                dS -= beta*(np.cos(p1n)-np.cos(p1o))
                p2o = links[tuple(pmn)][nu]+old-links[tuple(pmun)][nu]-links[tuple(pmn)][mu]
                p2n = links[tuple(pmn)][nu]+new-links[tuple(pmun)][nu]-links[tuple(pmn)][mu]
                dS -= beta*(np.cos(p2n)-np.cos(p2o))
            if dS < 0 or np.random.rand() < np.exp(-dS):
                links[idx][mu] = new

print("Thermalizing...")
for s in range(100): u1_sweep(links_u1, beta, L)

# Measure plaquettes in (0,1) and (0,2) planes — share direction 0
N = 500
P01 = []
P02 = []

print("Measuring...")
for meas in range(N):
    for _ in range(3): u1_sweep(links_u1, beta, L)
    w01 = 0; w02 = 0; c = 0
    for idx in np.ndindex(L,L,L,L):
        p0 = list(idx); p0[0]=(p0[0]+1)%L
        p1 = list(idx); p1[1]=(p1[1]+1)%L
        p2 = list(idx); p2[2]=(p2[2]+1)%L
        # (0,1) plaquette
        f01 = links_u1[idx][0]+links_u1[tuple(p0)][1]-links_u1[tuple(p1)][0]-links_u1[idx][1]
        w01 += np.cos(f01)
        # (0,2) plaquette
        f02 = links_u1[idx][0]+links_u1[tuple(p0)][2]-links_u1[tuple(p2)][0]-links_u1[idx][2]
        w02 += np.cos(f02)
        c += 1
    P01.append(w01/c)
    P02.append(w02/c)

P01 = np.array(P01); P02 = np.array(P02)
corr_u1 = np.corrcoef(P01, P02)[0,1]
joint = np.mean(P01*P02)
prod = np.mean(P01)*np.mean(P02)
conn = abs(joint - prod)

# Bootstrap
boots = [abs(np.mean(P01[np.random.randint(0,N,N)]*P02[np.random.randint(0,N,N)])-np.mean(P01[np.random.randint(0,N,N)])*np.mean(P02[np.random.randint(0,N,N)])) for _ in range(1000)]
err = np.std(boots)

print(f"\nU(1) plaquettes (0,1) vs (0,2) [share dir-0 links]:")
print(f"  <P01> = {P01.mean():.6f}, <P02> = {P02.mean():.6f}")
print(f"  <P01*P02> = {joint:.6f}")
print(f"  <P01><P02> = {prod:.6f}")
print(f"  |connected| = {conn:.6e} +/- {err:.6e} ({conn/err:.1f} sigma)")
print(f"  correlation = {corr_u1:.6f}")

# =====================================================
# PART 2: SU(2) — for comparison
# =====================================================
print()
print("=" * 60)
print("PART 2: SU(2)")
print("=" * 60)

# SU(2) link variables: 2x2 unitary matrices with det=1
# Parametrize as U = a0*I + i*(a1*s1 + a2*s2 + a3*s3) with sum(a^2)=1

I2 = np.eye(2, dtype=complex)
sigma = [np.array([[0,1],[1,0]],dtype=complex),
         np.array([[0,-1j],[1j,0]],dtype=complex),
         np.array([[1,0],[0,-1]],dtype=complex)]

def random_su2():
    a = np.random.randn(4)
    a = a / np.linalg.norm(a)
    return a[0]*I2 + 1j*(a[1]*sigma[0]+a[2]*sigma[1]+a[3]*sigma[2])

def su2_staple(links, idx, mu, nu, L):
    pm = list(idx); pm[mu]=(pm[mu]+1)%L
    pn = list(idx); pn[nu]=(pn[nu]+1)%L
    fwd = links[tuple(pm)][nu] @ links[tuple(pn)][mu].conj().T @ links[idx][nu].conj().T
    pmn = list(idx); pmn[nu]=(pmn[nu]-1)%L
    pmun = list(idx); pmun[mu]=(pmun[mu]+1)%L; pmun[nu]=(pmun[nu]-1)%L
    bwd = links[tuple(pmun)][nu].conj().T @ links[tuple(pmn)][mu].conj().T @ links[tuple(pmn)][nu]
    return fwd + bwd

def su2_heat_bath(links, beta, L):
    for idx in np.ndindex(L,L,L,L):
        for mu in range(4):
            S = np.zeros((2,2), dtype=complex)
            for nu in range(4):
                if nu == mu: continue
                S += su2_staple(links, idx, mu, nu, L)
            # Creutz heat bath: generate new link proportional to exp(beta/2 * Re(Tr(U*S)))
            # Simplified: try a few random SU(2) elements and accept best
            best_U = links[idx][mu]
            best_action = np.real(np.trace(best_U @ S))
            for _ in range(10):
                trial = random_su2()
                trial_action = np.real(np.trace(trial @ S))
                if trial_action > best_action or np.random.rand() < np.exp(beta/2*(trial_action-best_action)):
                    best_U = trial
                    best_action = trial_action
            links[idx][mu] = best_U

# Initialize SU(2) lattice
links_su2 = np.zeros((L,L,L,L,4), dtype=object)
for idx in np.ndindex(L,L,L,L):
    for mu in range(4):
        links_su2[idx][mu] = random_su2()

# Use numpy array for SU(2) links
links_su2_arr = np.zeros((L,L,L,L,4,2,2), dtype=complex)
for idx in np.ndindex(L,L,L,L):
    for mu in range(4):
        links_su2_arr[idx+(mu,)] = random_su2()

def su2_sweep_arr(links, beta, L):
    for idx in np.ndindex(L,L,L,L):
        for mu in range(4):
            S = np.zeros((2,2), dtype=complex)
            for nu in range(4):
                if nu == mu: continue
                pm = list(idx); pm[mu]=(pm[mu]+1)%L
                pn = list(idx); pn[nu]=(pn[nu]+1)%L
                pmn = list(idx); pmn[nu]=(pmn[nu]-1)%L
                pmun = list(idx); pmun[mu]=(pmun[mu]+1)%L; pmun[nu]=(pmun[nu]-1)%L
                fwd = links[tuple(pm)+(nu,)] @ links[tuple(pn)+(mu,)].conj().T @ links[tuple(idx)+(nu,)].conj().T
                bwd = links[tuple(pmun)+(nu,)].conj().T @ links[tuple(pmn)+(mu,)].conj().T @ links[tuple(pmn)+(nu,)]
                S += fwd + bwd
            best_U = links[tuple(idx)+(mu,)]
            best_a = np.real(np.trace(best_U @ S))
            for _ in range(8):
                trial = random_su2()
                trial_a = np.real(np.trace(trial @ S))
                if trial_a > best_a or np.random.rand() < np.exp(beta/2*(trial_a - best_a)):
                    best_U = trial
                    best_a = trial_a
            links[tuple(idx)+(mu,)] = best_U

print("Thermalizing SU(2)...")
for s in range(30):
    su2_sweep_arr(links_su2_arr, 2.0, L)  # beta=2.0 for SU(2)
    if s % 10 == 0: print(f"  Sweep {s}")

N_su2 = 200
SP01 = []; SP02 = []

print("Measuring SU(2)...")
for meas in range(N_su2):
    su2_sweep_arr(links_su2_arr, 2.0, L)
    w01 = 0; w02 = 0; c = 0
    for idx in np.ndindex(L,L,L,L):
        pm = list(idx); pm[0]=(pm[0]+1)%L
        p1 = list(idx); p1[1]=(p1[1]+1)%L
        p2 = list(idx); p2[2]=(p2[2]+1)%L
        # (0,1) plaquette trace
        plaq01 = links_su2_arr[tuple(idx)+(0,)] @ links_su2_arr[tuple(pm)+(1,)] @ links_su2_arr[tuple(p1)+(0,)].conj().T @ links_su2_arr[tuple(idx)+(1,)].conj().T
        w01 += np.real(np.trace(plaq01)) / 2
        # (0,2) plaquette trace
        plaq02 = links_su2_arr[tuple(idx)+(0,)] @ links_su2_arr[tuple(pm)+(2,)] @ links_su2_arr[tuple(p2)+(0,)].conj().T @ links_su2_arr[tuple(idx)+(2,)].conj().T
        w02 += np.real(np.trace(plaq02)) / 2
        c += 1
    SP01.append(w01/c); SP02.append(w02/c)
    if (meas+1) % 50 == 0: print(f"  Measurement {meas+1}/{N_su2}")

SP01 = np.array(SP01); SP02 = np.array(SP02)
corr_su2 = np.corrcoef(SP01, SP02)[0,1]
joint_su2 = np.mean(SP01*SP02)
prod_su2 = np.mean(SP01)*np.mean(SP02)
conn_su2 = abs(joint_su2 - prod_su2)

boots_su2 = [abs(np.mean(SP01[np.random.randint(0,N_su2,N_su2)]*SP02[np.random.randint(0,N_su2,N_su2)])-np.mean(SP01[np.random.randint(0,N_su2,N_su2)])*np.mean(SP02[np.random.randint(0,N_su2,N_su2)])) for _ in range(1000)]
err_su2 = np.std(boots_su2)

print(f"\nSU(2) plaquettes (0,1) vs (0,2) [share dir-0 links]:")
print(f"  <P01> = {SP01.mean():.6f}, <P02> = {SP02.mean():.6f}")
print(f"  <P01*P02> = {joint_su2:.6f}")
print(f"  <P01><P02> = {prod_su2:.6f}")
print(f"  |connected| = {conn_su2:.6e} +/- {err_su2:.6e} ({conn_su2/err_su2:.1f} sigma)")
print(f"  correlation = {corr_su2:.6f}")

# =====================================================
# COMPARISON
# =====================================================
print()
print("=" * 60)
print("COMPARISON")
print("=" * 60)
print(f"  U(1)  correlation: {corr_u1:.6f}")
print(f"  SU(2) correlation: {corr_su2:.6f}")
print()
if abs(corr_su2) > 3*abs(corr_u1) and abs(corr_su2) > 0.1:
    print("  SU(2) has significantly stronger correlation than U(1)")
    print("  -> Non-abelian gauge constraint creates structural correlation")
    print("  -> Confirms: delta(SU(2)) < delta(U(1))")
elif abs(corr_u1) < 0.05 and abs(corr_su2) > 0.1:
    print("  U(1) factorizes, SU(2) does not")
    print("  -> Abelian discrimination confirmed")
else:
    print(f"  Both have similar correlation levels")
    print(f"  -> Need larger lattice or more statistics to discriminate")
