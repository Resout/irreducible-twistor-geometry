"""
Peter-Weyl Decomposition: the actual discriminator.

Two plaquettes sharing link U in direction 0:
  P01 = Tr(U * A) / 2   where A = staple in (0,1) plane
  P02 = Tr(U * B) / 2   where B = staple in (0,2) plane

Peter-Weyl under Haar measure on the shared link:
  SU(2): <Tr(UA) Tr(UB)>_U = (1/2) Tr(A†B)  -- non-zero generically
  U(1):  <cos(θ+α) cos(θ+β)>_θ = (1/2) cos(α-β)  -- also non-zero

Both have surviving cross terms. The difference is in the
REPRESENTATION CONTENT of the cross term:
  SU(2): Tr(A†B) contains all spin-j contributions
  U(1):  cos(α-β) is charge-0 only

Measure both on thermalized lattices and compare distributions.
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

# SU(2) lattice
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

# U(1) lattice
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

# Thermalize both
print("Thermalizing SU(2)...")
t0 = time()
for s in range(50): su2_sweep(su2, 2.2)
print(f"  Done in {time()-t0:.0f}s")

print("Thermalizing U(1)...")
t0 = time()
for s in range(50): u1_sweep(u1, 1.0)
print(f"  Done in {time()-t0:.0f}s")

# Measure
n_meas = 300
su2_trAB, u1_cosAB = [], []
su2_conn, u1_conn = [], []

print(f"Measuring ({n_meas} configs)...")
t0 = time()
for m in range(n_meas):
    for _ in range(3):
        su2_sweep(su2, 2.2)
        u1_sweep(u1, 1.0)

    su2_p01_cfg, su2_p02_cfg = [], []
    u1_p01_cfg, u1_p02_cfg = [], []

    for x in range(L):
        for y in range(L):
            for z in range(L):
                for t in range(L):
                    pos = (x,y,z,t)
                    pm0 = list(pos); pm0[0]=(pm0[0]+1)%L; pm0=tuple(pm0)
                    p1 = list(pos); p1[1]=(p1[1]+1)%L; p1=tuple(p1)
                    p2 = list(pos); p2[2]=(p2[2]+1)%L; p2=tuple(p2)

                    # SU(2)
                    A = su2[pm0+(1,)] @ su2[p1+(0,)].conj().T @ su2[pos+(1,)].conj().T
                    B = su2[pm0+(2,)] @ su2[p2+(0,)].conj().T @ su2[pos+(2,)].conj().T
                    su2_trAB.append(np.real(np.trace(A.conj().T @ B)) / 2)
                    p01 = np.real(np.trace(su2[pos+(0,)] @ A)) / 2
                    p02 = np.real(np.trace(su2[pos+(0,)] @ B)) / 2
                    su2_p01_cfg.append(p01)
                    su2_p02_cfg.append(p02)

                    # U(1)
                    alpha = u1[pm0][1] - u1[p1][0] - u1[x,y,z,t,1]
                    beta_ph = u1[pm0][2] - u1[p2][0] - u1[x,y,z,t,2]
                    u1_cosAB.append(np.cos(alpha - beta_ph) / 2)
                    u1_p01_cfg.append(np.cos(u1[x,y,z,t,0] + alpha))
                    u1_p02_cfg.append(np.cos(u1[x,y,z,t,0] + beta_ph))

    # Per-configuration connected correlator
    s01 = np.array(su2_p01_cfg); s02 = np.array(su2_p02_cfg)
    su2_conn.append(np.mean(s01*s02) - np.mean(s01)*np.mean(s02))
    u01 = np.array(u1_p01_cfg); u02 = np.array(u1_p02_cfg)
    u1_conn.append(np.mean(u01*u02) - np.mean(u01)*np.mean(u02))

    if (m+1) % 100 == 0:
        print(f"  {m+1}/{n_meas} ({time()-t0:.0f}s)")

su2_trAB = np.array(su2_trAB)
u1_cosAB = np.array(u1_cosAB)
su2_conn = np.array(su2_conn)
u1_conn = np.array(u1_conn)

print(f"  Done in {time()-t0:.0f}s")
print()

# Results
print("=" * 70)
print("PETER-WEYL CROSS TERM DISTRIBUTIONS")
print("=" * 70)
print()
print("SU(2): Tr(A'B)/2          U(1): cos(a-b)/2")
print(f"  Mean:   {su2_trAB.mean():>10.6f}         {u1_cosAB.mean():>10.6f}")
print(f"  Std:    {su2_trAB.std():>10.6f}         {u1_cosAB.std():>10.6f}")
sk_su2 = ((su2_trAB - su2_trAB.mean())**3).mean() / su2_trAB.std()**3
sk_u1 = ((u1_cosAB - u1_cosAB.mean())**3).mean() / u1_cosAB.std()**3
ku_su2 = ((su2_trAB - su2_trAB.mean())**4).mean() / su2_trAB.std()**4
ku_u1 = ((u1_cosAB - u1_cosAB.mean())**4).mean() / u1_cosAB.std()**4
print(f"  Skew:   {sk_su2:>10.4f}         {sk_u1:>10.4f}")
print(f"  Kurt:   {ku_su2:>10.4f}         {ku_u1:>10.4f}")
print()

print("PER-CONFIGURATION CONNECTED CORRELATOR <P01*P02> - <P01><P02>:")
print(f"  SU(2): mean = {su2_conn.mean():.6e}, std = {su2_conn.std():.6e}, sig = {su2_conn.mean()/su2_conn.std()*np.sqrt(len(su2_conn)):.1f}sigma")
print(f"  U(1):  mean = {u1_conn.mean():.6e}, std = {u1_conn.std():.6e}, sig = {u1_conn.mean()/u1_conn.std()*np.sqrt(len(u1_conn)):.1f}sigma")
print()

print("INTERPRETATION:")
print("  The cross term exists for BOTH theories (shared link).")
print("  The connected correlator measures the actual signal.")
if abs(su2_conn.mean()) > 2*su2_conn.std()/np.sqrt(len(su2_conn)):
    print(f"  SU(2): SIGNIFICANT connected correlator ({su2_conn.mean()/su2_conn.std()*np.sqrt(len(su2_conn)):.1f}sigma)")
else:
    print(f"  SU(2): connected correlator consistent with zero")
if abs(u1_conn.mean()) > 2*u1_conn.std()/np.sqrt(len(u1_conn)):
    print(f"  U(1):  SIGNIFICANT connected correlator ({u1_conn.mean()/u1_conn.std()*np.sqrt(len(u1_conn)):.1f}sigma)")
else:
    print(f"  U(1):  connected correlator consistent with zero")
