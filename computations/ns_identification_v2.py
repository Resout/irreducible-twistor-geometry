"""
NS Identification v2: Fixed vorticity extraction.

The vorticity is the su(2) content of the curvature F = da + a^a,
extracted via Tr(sigma_k * F_{mu nu}) / 2, NOT via Tr(F) (which is zero
for su(2)-valued forms).

The curvature 2-form F_{mu nu} in 3D has 3 independent components.
In the Hodge dual: omega_k = (1/2) epsilon_{kij} F_{ij}.
Each F_{ij} is a 2x2 su(2) matrix. The physical vorticity components
are the su(2) components: omega_k^a = (1/2) epsilon_{kij} Tr(sigma_a F_{ij}).
"""

import numpy as np
from itertools import product as iprod
from scipy.optimize import brentq
from numpy.linalg import lstsq, inv

H = 3; FLOOR = 1.0/27.0

# Pauli matrices
sigma = [
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex)
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
    """Mass function -> 2x2 matrix."""
    s1,s2,s3,th = m
    return (th*I2 + s1*sigma[0] + s2*sigma[1] + s3*sigma[2]) / np.sqrt(2)

def su2_components_of(M):
    """Extract su(2) components: M = (a0*I + a_k*sigma_k)/sqrt(2)."""
    sq2 = np.sqrt(2)
    a0 = np.real(np.trace(M*sq2)) / 2
    a = np.array([np.real(np.trace(sigma[k] @ (M*sq2))) / 2 for k in range(3)])
    return a0, a

# Find equilibrium
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
print()

# ============================================================
# Set up lattice with non-trivial vorticity
# ============================================================

N = 12
print(f"Setting up {N}^3 lattice with vortex-tube perturbation...")

lattice = np.zeros((N, N, N, 4))
amp = 0.03  # perturbation amplitude

for ix, iy, iz in iprod(range(N), repeat=3):
    x = 2*np.pi*ix/N
    y = 2*np.pi*iy/N
    z = 2*np.pi*iz/N

    # Arnold-Beltrami-Childress (ABC) flow perturbation
    # This is a known solution of the Euler equations
    perturb = amp * np.array([
        np.sin(z) + np.cos(y),   # s1
        np.sin(x) + np.cos(z),   # s2
        np.sin(y) + np.cos(x),   # s3
        0.0
    ])

    m = m_star + perturb
    m = np.abs(m)  # keep positive
    m = m / np.sum(m)  # L1 = 1
    lattice[ix,iy,iz] = enforce_floor(m)

# ============================================================
# Extract vorticity from the connection curvature
# ============================================================

def extract_vorticity_su2(lattice, N):
    """Extract vorticity as the su(2) components of the curvature F = da + a^a.

    Connection: a_mu(x) = M(x)^{-1} * (M(x+mu) - M(x-mu)) / 2
    Curvature: F_{mu,nu} = a_nu(x+mu) - a_nu(x-mu))/2 - (same for mu<->nu) + [a_mu, a_nu]
    Vorticity: omega_k^a = (1/2) eps_{kij} Tr(sigma_a * F_{ij}) / 2
    """
    omega = np.zeros((N,N,N,3))  # 3 components of vorticity (using dominant su(2) direction)

    for ix,iy,iz in iprod(range(N), repeat=3):
        M = to_matrix(lattice[ix,iy,iz])
        Mi = inv(M)

        # Connection: a_d = M^{-1} * (M(x+d) - M(x-d)) / 2
        a = []
        for d in range(3):
            ip = [ix,iy,iz]; im = [ix,iy,iz]
            ip[d] = (ip[d]+1)%N; im[d] = (im[d]-1)%N
            Mp = to_matrix(lattice[tuple(ip)])
            Mm = to_matrix(lattice[tuple(im)])
            a.append(Mi @ (Mp - Mm) / 2.0)

        # Curvature F_{mu,nu} ≈ [a_mu, a_nu] + derivative terms
        # For leading order, use the commutator (the a^a term)
        # Plus the discrete derivative of a
        F = [[None]*3 for _ in range(3)]
        for mu in range(3):
            for nu in range(mu+1, 3):
                # Discrete: partial_mu a_nu - partial_nu a_mu
                ip_mu = [ix,iy,iz]; im_mu = [ix,iy,iz]
                ip_mu[mu] = (ip_mu[mu]+1)%N; im_mu[mu] = (im_mu[mu]-1)%N

                # a_nu at x+mu and x-mu
                M_pmu = to_matrix(lattice[tuple(ip_mu)])
                Mi_pmu = inv(M_pmu)
                ip2 = list(ip_mu); im2 = list(ip_mu)
                ip2[nu] = (ip2[nu]+1)%N; im2[nu] = (im2[nu]-1)%N
                a_nu_pmu = Mi_pmu @ (to_matrix(lattice[tuple(ip2)]) - to_matrix(lattice[tuple(im2)])) / 2

                M_mmu = to_matrix(lattice[tuple(im_mu)])
                Mi_mmu = inv(M_mmu)
                ip3 = list(im_mu); im3 = list(im_mu)
                ip3[nu] = (ip3[nu]+1)%N; im3[nu] = (im3[nu]-1)%N
                a_nu_mmu = Mi_mmu @ (to_matrix(lattice[tuple(ip3)]) - to_matrix(lattice[tuple(im3)])) / 2

                da_mu_nu = (a_nu_pmu - a_nu_mmu) / 2  # partial_mu a_nu

                # Similarly partial_nu a_mu
                ip_nu = [ix,iy,iz]; im_nu = [ix,iy,iz]
                ip_nu[nu] = (ip_nu[nu]+1)%N; im_nu[nu] = (im_nu[nu]-1)%N

                M_pnu = to_matrix(lattice[tuple(ip_nu)])
                Mi_pnu = inv(M_pnu)
                ip4 = list(ip_nu); im4 = list(ip_nu)
                ip4[mu] = (ip4[mu]+1)%N; im4[mu] = (im4[mu]-1)%N
                a_mu_pnu = Mi_pnu @ (to_matrix(lattice[tuple(ip4)]) - to_matrix(lattice[tuple(im4)])) / 2

                M_mnu = to_matrix(lattice[tuple(im_nu)])
                Mi_mnu = inv(M_mnu)
                ip5 = list(im_nu); im5 = list(im_nu)
                ip5[mu] = (ip5[mu]+1)%N; im5[mu] = (im5[mu]-1)%N
                a_mu_mnu = Mi_mnu @ (to_matrix(lattice[tuple(ip5)]) - to_matrix(lattice[tuple(im5)])) / 2

                da_nu_mu = (a_mu_pnu - a_mu_mnu) / 2  # partial_nu a_mu

                # Commutator
                comm = a[mu] @ a[nu] - a[nu] @ a[mu]

                F[mu][nu] = da_mu_nu - da_nu_mu + comm

        # Hodge dual: omega_k = (1/2) eps_{kij} F_{ij}
        # omega_0 = F_{12}, omega_1 = F_{20}, omega_2 = F_{01}
        F_01 = F[0][1]
        F_02 = F[0][2]
        F_12 = F[1][2]

        # Extract the DOMINANT su(2) component (sigma_3 direction at equilibrium)
        # Actually, take the full su(2) norm: |omega_k| = sqrt(sum_a (Tr(sigma_a F_{ij}))^2)
        # But for comparison with NS, we need a scalar vorticity.
        # Use: omega_k = Tr(sigma_3 * F_{ij}) / 2  (project onto dominant su(2) axis)
        # Actually for a proper vector, take all three:
        for a_idx in range(3):
            omega[ix,iy,iz,0] += np.real(np.trace(sigma[a_idx] @ F_12)) / 2 if a_idx == 2 else 0
            omega[ix,iy,iz,1] += np.real(np.trace(sigma[a_idx] @ F[0][2])) / 2 if a_idx == 2 else 0  # F_{20} = -F_{02}

        # Simpler: just take the sigma_3 component (dominant axis)
        omega[ix,iy,iz,0] = np.real(np.trace(sigma[2] @ F_12)) / 2     # omega_z from F_{xy}
        omega[ix,iy,iz,1] = -np.real(np.trace(sigma[2] @ F_02)) / 2    # omega_y from F_{zx}
        omega[ix,iy,iz,2] = np.real(np.trace(sigma[2] @ F_01)) / 2     # omega_x from F_{yz}

    return omega

print("Extracting vorticity (su(2) components of curvature)...")
omega_0 = extract_vorticity_su2(lattice, N)
print(f"  max|omega| = {np.max(np.abs(omega_0)):.8f}")
print(f"  mean|omega| = {np.mean(np.abs(omega_0)):.8f}")
print(f"  omega_0 nonzero: {np.sum(np.abs(omega_0) > 1e-15)}/{omega_0.size}")
print()

if np.max(np.abs(omega_0)) < 1e-12:
    print("ERROR: vorticity still zero. Aborting.")
    exit(1)

# Evolve one DS step
print("Evolving one DS step with neighbour averaging...")
lattice_new = np.zeros_like(lattice)
for ix,iy,iz in iprod(range(N), repeat=3):
    e_local = np.zeros(4)
    for d in range(3):
        for sign in [+1,-1]:
            idx = [ix,iy,iz]; idx[d] = (idx[d]+sign)%N
            e_local += lattice[tuple(idx)]
    e_local /= 6.0
    e_local = enforce_floor(e_local)
    m_new, _ = ds_combine(lattice[ix,iy,iz], e_local)
    lattice_new[ix,iy,iz] = m_new

print("Extracting evolved vorticity...")
omega_1 = extract_vorticity_su2(lattice_new, N)
delta_omega = omega_1 - omega_0
print(f"  max|delta_omega| = {np.max(np.abs(delta_omega)):.8f}")
print()

# Discrete Laplacian of omega
print("Computing Laplacian...")
lap_omega = np.zeros_like(omega_0)
for d in range(3):
    for k in range(3):
        for ix,iy,iz in iprod(range(N), repeat=3):
            ip = [ix,iy,iz]; ip[d]=(ip[d]+1)%N
            im = [ix,iy,iz]; im[d]=(im[d]-1)%N
            lap_omega[ix,iy,iz,k] += omega_0[tuple(ip)][k] + omega_0[tuple(im)][k] - 2*omega_0[ix,iy,iz,k]

# Biot-Savart via FFT
print("Computing velocity via Biot-Savart...")
omega_hat = np.fft.fftn(omega_0, axes=(0,1,2))
kx = np.fft.fftfreq(N)*2*np.pi; ky=kx.copy(); kz=kx.copy()
KX,KY,KZ = np.meshgrid(kx,ky,kz, indexing='ij')
K2 = KX**2+KY**2+KZ**2; K2[0,0,0]=1.0

u_hat = np.zeros_like(omega_hat)
u_hat[:,:,:,0] = 1j*(KY*omega_hat[:,:,:,2] - KZ*omega_hat[:,:,:,1]) / K2
u_hat[:,:,:,1] = 1j*(KZ*omega_hat[:,:,:,0] - KX*omega_hat[:,:,:,2]) / K2
u_hat[:,:,:,2] = 1j*(KX*omega_hat[:,:,:,1] - KY*omega_hat[:,:,:,0]) / K2
u = np.real(np.fft.ifftn(u_hat, axes=(0,1,2)))
print(f"  max|u| = {np.max(np.abs(u)):.8f}")

# Advection-stretching
print("Computing advection-stretching...")
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

print(f"  max|advection| = {np.max(np.abs(adv)):.8f}")
print()

# Regression
do = delta_omega.flatten()
la = lap_omega.flatten()
ad = adv.flatten()

A_mat = np.column_stack([la, ad, np.ones_like(do)])
coeffs, _, _, _ = lstsq(A_mat, do, rcond=None)
nu_fit, alpha_fit, const_fit = coeffs

pred = nu_fit*lap_omega + alpha_fit*adv + const_fit
res = delta_omega - pred

SS_res = np.sum(res**2)
SS_tot = np.sum((do - np.mean(do))**2)
R2 = 1 - SS_res/SS_tot if SS_tot > 0 else 0

# Diffusion-only
A_diff = np.column_stack([la, np.ones_like(do)])
c_diff, _, _, _ = lstsq(A_diff, do, rcond=None)
pred_diff = c_diff[0]*lap_omega + c_diff[1]
res_diff = delta_omega - pred_diff
R2_diff = 1 - np.sum(res_diff**2)/SS_tot if SS_tot > 0 else 0

print("="*60)
print("RESULTS: delta_omega = nu*laplacian + alpha*advection + const")
print("="*60)
print(f"  nu (diffusion)    = {nu_fit:.6f}")
print(f"  alpha (advection) = {alpha_fit:.6f}")
print(f"  const             = {const_fit:.2e}")
print(f"  R^2 (full model)  = {R2:.6f}")
print(f"  R^2 (diffusion only) = {R2_diff:.6f}")
print(f"  R^2 improvement   = {R2 - R2_diff:.6f}")
print()
print(f"  NS prediction: alpha = 1.0")
print(f"  Measured: alpha = {alpha_fit:.6f}")
if abs(alpha_fit) > 0.01:
    print(f"  Deviation: {abs(alpha_fit - 1.0)*100:.1f}%")
else:
    print(f"  WARNING: alpha ~ 0, advection term may not be detectable at this perturbation scale")
print()
print(f"  ||delta_omega|| = {np.sqrt(np.sum(delta_omega**2)):.8f}")
print(f"  ||laplacian||   = {np.sqrt(np.sum(lap_omega**2)):.8f}")
print(f"  ||advection||   = {np.sqrt(np.sum(adv**2)):.8f}")
print(f"  ||residual||    = {np.sqrt(np.sum(res**2)):.8f}")
print()
print("DONE.")
