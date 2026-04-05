"""
NS Identification v4: Full su(2)-valued vorticity.

The descended vorticity is NOT a scalar — it's an su(2)-valued vector:
  B_i^a(x) = (1/2) eps_{ijk} Tr(sigma_a * F_{jk}(x)) / 2

That's 3 spatial × 3 su(2) = 9 components per site.

The NS nonlinear term involves TWO cross products:
  - Spatial: eps_{ijk} (from Hodge dual / curl)
  - su(2): eps_{abc} (from [sigma_a, sigma_b] = 2i eps_{abc} sigma_c)

A scalar projection (taking only sigma_3) kills the su(2) cross product.
The full 9-component system should show the coupling.

Strategy:
1. Extract B_i^a from lattice (9 components)
2. DS step with neighbours → delta_B_i^a
3. Subtract self-evidence reaction → spatial delta_B
4. Regress against:
   (a) Ordinary Laplacian of B (diffusion)
   (b) Commutator coupling: C_i^a = eps_{abc} sum_j a_j^b * (partial_j B_i^c)
       This is the gauge-covariant part that mixes su(2) indices
   (c) Full NS-like advection using Biot-Savart on each su(2) component
"""

import numpy as np
from itertools import product as iprod
from scipy.optimize import brentq
from numpy.linalg import lstsq, inv, norm
import sys

H = 3; FLOOR = 1.0/27.0

sigma = [
    np.array([[0,1],[1,0]], dtype=complex),
    np.array([[0,-1j],[1j,0]], dtype=complex),
    np.array([[1,0],[0,-1]], dtype=complex)
]

def enforce_floor(m):
    s = m[:3].copy(); th = m[3]
    born = th**2 / (np.sum(s**2) + th**2)
    if born >= FLOOR: return m.copy()
    S = np.sum(s); Sq = np.sum(s**2); R = Sq/S**2
    t = (np.sqrt(26*R) - R)/(26-R); alpha = (1-t)/S
    out = m.copy(); out[:3] = s*alpha; out[3] = t
    return out

def ds_combine(m, e):
    s,th = m[:3],m[3]; se,ph = e[:3],e[3]
    sn = s*se + s*ph + th*se; tn = th*ph
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
    d = 1-K; out = np.zeros(4); out[:3]=sn/d; out[3]=tn/d
    return enforce_floor(out), K

def to_matrix(m):
    s1,s2,s3,th = m
    return (th*np.eye(2,dtype=complex) + s1*sigma[0] + s2*sigma[1] + s3*sigma[2]) / np.sqrt(2)

def find_eq():
    def K_res(p):
        pw=(1-p)/2; sc=26./27.
        raw = np.array([np.sqrt(p*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR)])
        e = raw/np.sum(raw); m = np.array([0.4,0.2,0.2,0.2])
        for _ in range(3000): m,_ = ds_combine(m,e)
        K = sum(m[i]*e[j] for i in range(3) for j in range(3) if i!=j)
        return K - 7./30.
    p = brentq(K_res, 0.92, 0.94, xtol=1e-15)
    pw=(1-p)/2; sc=26./27.
    raw = np.array([np.sqrt(p*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR)])
    e = raw/np.sum(raw); m = np.array([0.4,0.2,0.2,0.2])
    for _ in range(3000): m,_ = ds_combine(m,e)
    return m, e

m_star, e_star = find_eq()
print(f"m* = {m_star}"); sys.stdout.flush()

N = 12; amp = 0.15
print(f"\nSetting up {N}^3 lattice, amp={amp}..."); sys.stdout.flush()

lattice = np.zeros((N,N,N,4))
for ix,iy,iz in iprod(range(N), repeat=3):
    x=2*np.pi*ix/N; y=2*np.pi*iy/N; z=2*np.pi*iz/N
    perturb = amp * np.array([np.sin(z)+np.cos(y), np.sin(x)+np.cos(z),
                              np.sin(y)+np.cos(x), 0.0])
    m = np.abs(m_star + perturb); m /= np.sum(m)
    lattice[ix,iy,iz] = enforce_floor(m)

# ============================================================
# Full su(2) vorticity extraction: B[ix,iy,iz,spatial,su2] = B_i^a
# ============================================================

def get_connection(lattice, ix, iy, iz, N):
    """Return a[d] = M^{-1}(M(x+d)-M(x-d))/2, as 2x2 matrices, plus M, M^{-1}."""
    M = to_matrix(lattice[ix,iy,iz]); Mi = inv(M)
    a = []
    for d in range(3):
        ip=[ix,iy,iz]; im=[ix,iy,iz]
        ip[d]=(ip[d]+1)%N; im[d]=(im[d]-1)%N
        a.append(Mi @ (to_matrix(lattice[tuple(ip)]) - to_matrix(lattice[tuple(im)])) / 2)
    return a, M, Mi

def extract_full_B(lattice, N):
    """Extract B_i^a = (1/2) eps_{ijk} Tr(sigma_a F_{jk}) / 2.
    Returns array of shape (N,N,N,3,3): [site_x,site_y,site_z,spatial_idx,su2_idx]."""
    B = np.zeros((N,N,N,3,3))

    for ix,iy,iz in iprod(range(N), repeat=3):
        a, M, Mi = get_connection(lattice, ix, iy, iz, N)

        # Curvature: F_{mu,nu} for (mu,nu) in {(0,1),(0,2),(1,2)}
        # F_{mn} = partial_m a_n - partial_n a_m + [a_m, a_n]
        F = {}
        for mu,nu in [(0,1),(0,2),(1,2)]:
            # partial_mu(a_nu): evaluate a_nu at x+mu and x-mu
            ip=[ix,iy,iz]; im=[ix,iy,iz]
            ip[mu]=(ip[mu]+1)%N; im[mu]=(im[mu]-1)%N
            a_nu_p = get_connection(lattice,*ip,N)[0][nu]
            a_nu_m = get_connection(lattice,*im,N)[0][nu]
            d_mu_a_nu = (a_nu_p - a_nu_m) / 2

            ip2=[ix,iy,iz]; im2=[ix,iy,iz]
            ip2[nu]=(ip2[nu]+1)%N; im2[nu]=(im2[nu]-1)%N
            a_mu_p = get_connection(lattice,*ip2,N)[0][mu]
            a_mu_m = get_connection(lattice,*im2,N)[0][mu]
            d_nu_a_mu = (a_mu_p - a_mu_m) / 2

            comm = a[mu] @ a[nu] - a[nu] @ a[mu]
            F[(mu,nu)] = d_mu_a_nu - d_nu_a_mu + comm

        # Hodge dual: B_0 = F_{12}, B_1 = F_{20} = -F_{02}, B_2 = F_{01}
        F_components = [F[(1,2)], -F[(0,2)], F[(0,1)]]

        for spatial in range(3):
            for su2_a in range(3):
                B[ix,iy,iz,spatial,su2_a] = np.real(np.trace(sigma[su2_a] @ F_components[spatial])) / 2

    return B

print("Extracting full su(2) vorticity B_i^a..."); sys.stdout.flush()
B0 = extract_full_B(lattice, N)
print(f"  B shape: {B0.shape}")
print(f"  max|B| = {np.max(np.abs(B0)):.8f}")
print(f"  per su(2) component: {[f'{np.max(np.abs(B0[:,:,:,:,a])):.6f}' for a in range(3)]}")
print(); sys.stdout.flush()

# ============================================================
# DS evolution: full step and self-evidence
# ============================================================

print("Computing self-evidence..."); sys.stdout.flush()
lattice_self = np.zeros_like(lattice)
for ix,iy,iz in iprod(range(N), repeat=3):
    lattice_self[ix,iy,iz], _ = ds_combine(lattice[ix,iy,iz], lattice[ix,iy,iz])

print("Computing full DS step..."); sys.stdout.flush()
lattice_full = np.zeros_like(lattice)
for ix,iy,iz in iprod(range(N), repeat=3):
    e_loc = np.zeros(4)
    for d in range(3):
        for s in [+1,-1]:
            idx=[ix,iy,iz]; idx[d]=(idx[d]+s)%N
            e_loc += lattice[tuple(idx)]
    e_loc /= 6.0; e_loc = enforce_floor(e_loc)
    lattice_full[ix,iy,iz], _ = ds_combine(lattice[ix,iy,iz], e_loc)

print("Extracting evolved B..."); sys.stdout.flush()
B_full = extract_full_B(lattice_full, N)
B_self = extract_full_B(lattice_self, N)

dB_total = B_full - B0
dB_self = B_self - B0
dB_spatial = dB_total - dB_self

print(f"  ||dB_total||   = {norm(dB_total):.8f}")
print(f"  ||dB_self||    = {norm(dB_self):.8f}")
print(f"  ||dB_spatial|| = {norm(dB_spatial):.8f}")
print(f"  Spatial: {norm(dB_spatial)/norm(dB_total)*100:.1f}%")
print(); sys.stdout.flush()

# ============================================================
# Compute candidate terms for the spatial evolution
# ============================================================

print("Computing Laplacian of B..."); sys.stdout.flush()
lapB = np.zeros_like(B0)
for d in range(3):
    for i in range(3):
        for a in range(3):
            for ix,iy,iz in iprod(range(N), repeat=3):
                ip=[ix,iy,iz]; ip[d]=(ip[d]+1)%N
                im=[ix,iy,iz]; im[d]=(im[d]-1)%N
                lapB[ix,iy,iz,i,a] += B0[tuple(ip)][i,a] + B0[tuple(im)][i,a] - 2*B0[ix,iy,iz,i,a]

print("Computing commutator coupling [a, dB]..."); sys.stdout.flush()
# C_i^a = eps_{abc} sum_j a_j^b * partial_j B_i^c
# This is the gauge-covariant derivative correction
# The su(2) structure constants: [sigma_a, sigma_b] = 2i eps_{abc} sigma_c
# So the commutator coupling in component form is:
# [a_j, partial_j B_i]^a = eps_{abc} a_j^b (partial_j B_i^c)
eps_tensor = np.zeros((3,3,3))
eps_tensor[0,1,2] = eps_tensor[1,2,0] = eps_tensor[2,0,1] = 1
eps_tensor[0,2,1] = eps_tensor[2,1,0] = eps_tensor[1,0,2] = -1

# First compute a_j^b at each site
print("  Extracting connection su(2) components..."); sys.stdout.flush()
a_comp = np.zeros((N,N,N,3,3))  # a_comp[...,j,b] = su(2) component b of connection in direction j
for ix,iy,iz in iprod(range(N), repeat=3):
    a_list, _, _ = get_connection(lattice, ix, iy, iz, N)
    for j in range(3):
        for b in range(3):
            a_comp[ix,iy,iz,j,b] = np.real(np.trace(sigma[b] @ a_list[j])) / 2

# Compute partial_j B_i^c
print("  Computing gradients of B..."); sys.stdout.flush()
dB_dx = np.zeros((N,N,N,3,3,3))  # dB_dx[...,j,i,c] = partial_j B_i^c
for j in range(3):
    for i in range(3):
        for c in range(3):
            for ix,iy,iz in iprod(range(N), repeat=3):
                ip=[ix,iy,iz]; ip[j]=(ip[j]+1)%N
                im=[ix,iy,iz]; im[j]=(im[j]-1)%N
                dB_dx[ix,iy,iz,j,i,c] = (B0[tuple(ip)][i,c] - B0[tuple(im)][i,c]) / 2

# Commutator: commB_i^a = sum_{j,b,c} eps_{abc} a_j^b (partial_j B_i^c)
print("  Contracting..."); sys.stdout.flush()
commB = np.zeros_like(B0)
for ix,iy,iz in iprod(range(N), repeat=3):
    for i in range(3):
        for a in range(3):
            val = 0.0
            for j in range(3):
                for b in range(3):
                    for c in range(3):
                        val += eps_tensor[a,b,c] * a_comp[ix,iy,iz,j,b] * dB_dx[ix,iy,iz,j,i,c]
            commB[ix,iy,iz,i,a] = val

print(f"  max|lapB|  = {np.max(np.abs(lapB)):.8f}")
print(f"  max|commB| = {np.max(np.abs(commB)):.8f}")
print(f"  ||lapB||   = {norm(lapB):.8f}")
print(f"  ||commB||  = {norm(commB):.8f}")
print(); sys.stdout.flush()

# ============================================================
# Regression: dB_spatial = nu * lapB + gamma * commB
# ============================================================

ds = dB_spatial.flatten()
la = lapB.flatten()
co = commB.flatten()

print("="*60)
print("REGRESSION: dB_spatial = nu * lapB + gamma * commB")
print("="*60)

# Full model
A_full = np.column_stack([la, co])
c_full, _, _, _ = lstsq(A_full, ds, rcond=None)
nu_fit, gamma_fit = c_full
pred_full = nu_fit*lapB + gamma_fit*commB
res_full = dB_spatial - pred_full
SS_res = np.sum(res_full**2)
SS_tot = np.sum((ds - np.mean(ds))**2)
R2_full = 1 - SS_res/SS_tot if SS_tot > 1e-30 else 0

# Diffusion only
c_diff, _, _, _ = lstsq(la.reshape(-1,1), ds, rcond=None)
pred_diff = c_diff[0]*lapB
R2_diff = 1 - np.sum((dB_spatial - pred_diff)**2)/SS_tot if SS_tot > 1e-30 else 0

# Commutator only
c_comm, _, _, _ = lstsq(co.reshape(-1,1), ds, rcond=None)
pred_comm = c_comm[0]*commB
R2_comm = 1 - np.sum((dB_spatial - pred_comm)**2)/SS_tot if SS_tot > 1e-30 else 0

print(f"  nu (diffusion)    = {nu_fit:.8f}")
print(f"  gamma (commutator)= {gamma_fit:.8f}")
print()
print(f"  R^2 (diff + comm) = {R2_full:.6f}")
print(f"  R^2 (diff only)   = {R2_diff:.6f}")
print(f"  R^2 (comm only)   = {R2_comm:.6f}")
print(f"  R^2 improvement from commutator: {R2_full - R2_diff:.6f}")
print()
print(f"  ||dB_spatial||    = {norm(dB_spatial):.8f}")
print(f"  ||nu*lapB||       = {norm(nu_fit*lapB):.8f}")
print(f"  ||gamma*commB||   = {norm(gamma_fit*commB):.8f}")
print(f"  ||residual||      = {norm(res_full):.8f}")
print()
print(f"  Commutator/diffusion ratio: {norm(gamma_fit*commB)/norm(nu_fit*lapB):.4f}")
print(f"  Commutator/residual ratio:  {norm(gamma_fit*commB)/norm(res_full):.4f}")
print()

# Check: does the su(2) coupling explain what the scalar projection missed?
print("="*60)
print("ANALYSIS: Is the commutator the missing NS term?")
print("="*60)
print()
print(f"The commutator [a, dB] couples su(2) components via eps_abc.")
print(f"This is the gauge-covariant derivative correction D_j B_i = partial_j B_i + [a_j, B_i].")
print(f"In the abelian limit (all a^b = 0), it vanishes -> pure diffusion.")
print(f"For non-abelian connections, it produces FIRST-DERIVATIVE bilinear terms")
print(f"of the form a * dB, which is exactly the advection-stretching structure.")
print()

if R2_full - R2_diff > 0.01:
    print(f"SIGNIFICANT: commutator explains {(R2_full-R2_diff)*100:.1f}% additional variance.")
    print(f"The su(2) coupling IS detectable and contributes to the spatial evolution.")
    if abs(gamma_fit - 1.0) < 0.5:
        print(f"gamma = {gamma_fit:.4f} is near 1.0 -> consistent with unit-coefficient NS coupling.")
    else:
        print(f"gamma = {gamma_fit:.4f} differs from 1.0 -> coefficient needs understanding.")
else:
    print(f"NOT YET SIGNIFICANT: commutator adds only {(R2_full-R2_diff)*100:.4f}% variance.")
    print(f"This could be because:")
    print(f"  - The perturbation is still too weak (commutator is quadratic in the field)")
    print(f"  - The lattice is too coarse for the derivative coupling to resolve")
    print(f"  - The gauge-covariant term is genuinely small at this equilibrium")

print()
print("DONE.")
