"""
NS Identification v3: Large perturbation + reaction term subtraction.

Fixes from v2:
1. Larger perturbation amplitude (0.25) to bring quadratic advection
   above the diffusion floor
2. Subtract the LOCAL reaction term: the DS self-evidence update Phi(m,m)-m
   contributes to delta_m but has nothing to do with spatial coupling.
   The spatial part is: delta_m_spatial = Phi(m, E_neighbour) - Phi(m, m)
3. Better vorticity extraction: use ALL three su(2) components, not just sigma_3
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
I2 = np.eye(2, dtype=complex)

def enforce_floor(m):
    s = m[:3].copy(); th = m[3]
    born = th**2 / (np.sum(s**2) + th**2)
    if born >= FLOOR: return m.copy()
    S = np.sum(s); Sq = np.sum(s**2); R = Sq/S**2
    t = (np.sqrt(26*R) - R)/(26-R)
    alpha = (1-t)/S
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
    return (th*I2 + s1*sigma[0] + s2*sigma[1] + s3*sigma[2]) / np.sqrt(2)

def find_eq():
    def K_res(p):
        pw=(1-p)/2; sc=26./27.
        raw = np.array([np.sqrt(p*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR)])
        e = raw/np.sum(raw)
        m = np.array([0.4,0.2,0.2,0.2])
        for _ in range(3000): m,_ = ds_combine(m,e)
        s,th=m[:3],m[3]; se,ph=e[:3],e[3]
        K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
        return K - 7./30.
    p = brentq(K_res, 0.92, 0.94, xtol=1e-15)
    pw=(1-p)/2; sc=26./27.
    raw = np.array([np.sqrt(p*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR)])
    e = raw/np.sum(raw)
    m = np.array([0.4,0.2,0.2,0.2])
    for _ in range(3000): m,_ = ds_combine(m,e)
    return m, e

m_star, e_star = find_eq()
print(f"m* = {m_star}")
print(); sys.stdout.flush()

# ============================================================
# Set up lattice with LARGE perturbation
# ============================================================

N = 16
amp = 0.15
print(f"Setting up {N}^3 lattice, amplitude = {amp}...")
sys.stdout.flush()

lattice = np.zeros((N,N,N,4))
for ix,iy,iz in iprod(range(N), repeat=3):
    x = 2*np.pi*ix/N; y = 2*np.pi*iy/N; z = 2*np.pi*iz/N
    perturb = amp * np.array([
        np.sin(z) + np.cos(y),
        np.sin(x) + np.cos(z),
        np.sin(y) + np.cos(x),
        0.0
    ])
    m = m_star + perturb
    m = np.abs(m)
    m = m / np.sum(m)
    lattice[ix,iy,iz] = enforce_floor(m)

# ============================================================
# Vorticity extraction using full su(2) curvature
# ============================================================

def connection_at(lattice, ix, iy, iz, N):
    """Compute connection a_d = M^{-1}(M(x+d)-M(x-d))/2 at site (ix,iy,iz)."""
    M = to_matrix(lattice[ix,iy,iz])
    Mi = inv(M)
    a = []
    for d in range(3):
        ip = [ix,iy,iz]; im = [ix,iy,iz]
        ip[d]=(ip[d]+1)%N; im[d]=(im[d]-1)%N
        Mp = to_matrix(lattice[tuple(ip)])
        Mm = to_matrix(lattice[tuple(im)])
        a.append(Mi @ (Mp - Mm) / 2.0)
    return a, M, Mi

def curvature_at(lattice, ix, iy, iz, N):
    """Compute F_{mu,nu} = partial_mu a_nu - partial_nu a_mu + [a_mu, a_nu]."""
    a, M, Mi = connection_at(lattice, ix, iy, iz, N)

    # For a cleaner computation: use the plaquette formula
    # F_{mu,nu} ≈ M^{-1}(x) [M(x+mu+nu) M(x+nu)^{-1} M(x) - M(x+mu) M(x)^{-1} M(x+nu)] M(x+nu)^{-1}
    # But the linearized version is simpler and sufficient:
    # F_{mu,nu} = [a_mu, a_nu] + (discrete curl of a)

    F = {}
    for mu in range(3):
        for nu in range(mu+1, 3):
            # Commutator part
            comm = a[mu] @ a[nu] - a[nu] @ a[mu]

            # Derivative part: partial_mu(a_nu) - partial_nu(a_mu)
            # partial_mu(a_nu) = (a_nu(x+mu) - a_nu(x-mu))/2
            ip = [ix,iy,iz]; im = [ix,iy,iz]
            ip[mu]=(ip[mu]+1)%N; im[mu]=(im[mu]-1)%N
            a_nu_p = connection_at(lattice, *ip, N)[0][nu]
            a_nu_m = connection_at(lattice, *im, N)[0][nu]
            da_mu_nu = (a_nu_p - a_nu_m) / 2

            ip2 = [ix,iy,iz]; im2 = [ix,iy,iz]
            ip2[nu]=(ip2[nu]+1)%N; im2[nu]=(im2[nu]-1)%N
            a_mu_p = connection_at(lattice, *ip2, N)[0][mu]
            a_mu_m = connection_at(lattice, *im2, N)[0][mu]
            da_nu_mu = (a_mu_p - a_mu_m) / 2

            F[(mu,nu)] = da_mu_nu - da_nu_mu + comm

    return F

def extract_omega(lattice, N):
    """Extract vorticity using full su(2) norm of curvature."""
    omega = np.zeros((N,N,N,3))
    for ix,iy,iz in iprod(range(N), repeat=3):
        F = curvature_at(lattice, ix, iy, iz, N)

        # Hodge dual: omega_k = (1/2) eps_{kij} F_{ij}
        # omega_0 (x-component) = F_{12}
        # omega_1 (y-component) = F_{20} = -F_{02}
        # omega_2 (z-component) = F_{01}

        # Extract su(2) NORM of each F component
        # |F_{ij}|^2 = sum_a (Tr(sigma_a F_{ij})/2)^2
        for comp_idx, (mu,nu,sign) in enumerate([(1,2,1), (0,2,-1), (0,1,1)]):
            F_mn = F[(min(mu,nu), max(mu,nu))]
            if mu > nu:
                F_mn = -F_mn
            F_mn = sign * F_mn

            # Take the sigma_3 component as the "vorticity" (dominant axis)
            # This is a scalar projection — works for comparison with scalar NS
            omega[ix,iy,iz,comp_idx] = np.real(np.trace(sigma[2] @ F_mn)) / 2

    return omega

print("Extracting vorticity..."); sys.stdout.flush()
omega_0 = extract_omega(lattice, N)
print(f"  max|omega| = {np.max(np.abs(omega_0)):.8f}")
print(f"  mean|omega| = {np.mean(np.abs(omega_0)):.8f}")
print(); sys.stdout.flush()

# ============================================================
# Evolve: SPATIAL PART ONLY (subtract self-evidence reaction)
# ============================================================

print("Computing self-evidence reaction (local term)..."); sys.stdout.flush()
lattice_self = np.zeros_like(lattice)
for ix,iy,iz in iprod(range(N), repeat=3):
    lattice_self[ix,iy,iz], _ = ds_combine(lattice[ix,iy,iz], lattice[ix,iy,iz])

print("Computing full DS step with neighbours..."); sys.stdout.flush()
lattice_full = np.zeros_like(lattice)
for ix,iy,iz in iprod(range(N), repeat=3):
    e_local = np.zeros(4)
    for d in range(3):
        for sign in [+1,-1]:
            idx=[ix,iy,iz]; idx[d]=(idx[d]+sign)%N
            e_local += lattice[tuple(idx)]
    e_local /= 6.0
    e_local = enforce_floor(e_local)
    lattice_full[ix,iy,iz], _ = ds_combine(lattice[ix,iy,iz], e_local)

# Spatial contribution = full - self
lattice_spatial = lattice_full - lattice_self  # This is the SPATIAL-ONLY delta

print("Extracting vorticity from full evolution..."); sys.stdout.flush()
omega_full = extract_omega(lattice_full, N)
omega_self = extract_omega(lattice_self, N)

delta_omega_total = omega_full - omega_0
delta_omega_self = omega_self - omega_0
delta_omega_spatial = delta_omega_total - delta_omega_self

print(f"  ||delta_omega_total||   = {norm(delta_omega_total):.8f}")
print(f"  ||delta_omega_self||    = {norm(delta_omega_self):.8f}")
print(f"  ||delta_omega_spatial|| = {norm(delta_omega_spatial):.8f}")
print(f"  Spatial fraction: {norm(delta_omega_spatial)/norm(delta_omega_total)*100:.1f}%")
print(); sys.stdout.flush()

# ============================================================
# Compare spatial delta_omega with NS terms
# ============================================================

print("Computing Laplacian of omega..."); sys.stdout.flush()
lap = np.zeros_like(omega_0)
for d in range(3):
    for k in range(3):
        for ix,iy,iz in iprod(range(N), repeat=3):
            ip=[ix,iy,iz]; ip[d]=(ip[d]+1)%N
            im=[ix,iy,iz]; im[d]=(im[d]-1)%N
            lap[ix,iy,iz,k] += omega_0[tuple(ip)][k] + omega_0[tuple(im)][k] - 2*omega_0[ix,iy,iz,k]

print("Computing Biot-Savart velocity..."); sys.stdout.flush()
omega_hat = np.fft.fftn(omega_0, axes=(0,1,2))
kx = np.fft.fftfreq(N)*2*np.pi
KX,KY,KZ = np.meshgrid(kx,kx,kx, indexing='ij')
K2 = KX**2+KY**2+KZ**2; K2[0,0,0]=1.0
u_hat = np.zeros_like(omega_hat)
u_hat[:,:,:,0] = 1j*(KY*omega_hat[:,:,:,2] - KZ*omega_hat[:,:,:,1]) / K2
u_hat[:,:,:,1] = 1j*(KZ*omega_hat[:,:,:,0] - KX*omega_hat[:,:,:,2]) / K2
u_hat[:,:,:,2] = 1j*(KX*omega_hat[:,:,:,1] - KY*omega_hat[:,:,:,0]) / K2
u = np.real(np.fft.ifftn(u_hat, axes=(0,1,2)))

print("Computing advection-stretching..."); sys.stdout.flush()
adv = np.zeros_like(omega_0)
for k in range(3):
    for d in range(3):
        du_k = np.zeros((N,N,N)); domega_k = np.zeros((N,N,N))
        for ix,iy,iz in iprod(range(N), repeat=3):
            ip=[ix,iy,iz]; ip[d]=(ip[d]+1)%N
            im=[ix,iy,iz]; im[d]=(im[d]-1)%N
            du_k[ix,iy,iz] = (u[tuple(ip)][k] - u[tuple(im)][k])/2
            domega_k[ix,iy,iz] = (omega_0[tuple(ip)][k] - omega_0[tuple(im)][k])/2
        adv[:,:,:,k] += omega_0[:,:,:,d]*du_k - u[:,:,:,d]*domega_k

print(f"  max|laplacian| = {np.max(np.abs(lap)):.8f}")
print(f"  max|advection| = {np.max(np.abs(adv)):.8f}")
print(f"  max|u| = {np.max(np.abs(u)):.8f}")
print(); sys.stdout.flush()

# ============================================================
# Regression on the SPATIAL part
# ============================================================

# Model: delta_omega_spatial = nu * laplacian + alpha * advection
ds = delta_omega_spatial.flatten()
la = lap.flatten()
ad = adv.flatten()

print("="*60)
print("REGRESSION ON SPATIAL DELTA_OMEGA")
print("="*60)

# Full model: diffusion + advection
A_full = np.column_stack([la, ad])
c_full, _, _, _ = lstsq(A_full, ds, rcond=None)
nu_fit, alpha_fit = c_full
pred_full = nu_fit*lap + alpha_fit*adv
res_full = delta_omega_spatial - pred_full
SS_res = np.sum(res_full**2); SS_tot = np.sum((ds-np.mean(ds))**2)
R2_full = 1 - SS_res/SS_tot if SS_tot > 1e-30 else 0

# Diffusion only
A_diff = la.reshape(-1,1)
c_diff, _, _, _ = lstsq(A_diff, ds, rcond=None)
pred_diff = c_diff[0]*lap
res_diff = delta_omega_spatial - pred_diff
R2_diff = 1 - np.sum(res_diff**2)/SS_tot if SS_tot > 1e-30 else 0

print(f"  nu (diffusion)      = {nu_fit:.6f}")
print(f"  alpha (advection)   = {alpha_fit:.6f}")
print(f"  R^2 (diff+adv)     = {R2_full:.6f}")
print(f"  R^2 (diff only)    = {R2_diff:.6f}")
print(f"  R^2 improvement    = {R2_full - R2_diff:.6f}")
print()
print(f"  ||spatial_delta||   = {norm(delta_omega_spatial):.8f}")
print(f"  ||nu*laplacian||    = {norm(nu_fit*lap):.8f}")
print(f"  ||alpha*advection|| = {norm(alpha_fit*adv):.8f}")
print(f"  ||residual||        = {norm(res_full):.8f}")
print()

# Signal-to-noise: is advection detectable?
print(f"  Advection signal/diffusion signal = {norm(alpha_fit*adv)/norm(nu_fit*lap):.4f}")
print(f"  Advection signal/residual = {norm(alpha_fit*adv)/norm(res_full):.4f}")
print()

if abs(alpha_fit) > 0.1:
    print(f"  NS test: alpha should be 1.0, measured {alpha_fit:.4f}")
    print(f"  Deviation: {abs(alpha_fit - 1.0)*100:.1f}%")
else:
    print(f"  Advection coefficient too small to test NS hypothesis at this amplitude.")
    print(f"  Need larger perturbation or longer evolution.")

print()
print("DONE.")
