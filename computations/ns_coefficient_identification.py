"""
NS Coefficient Identification.

We know the DS spatial evolution is:
  dm_spatial = J_eff * lap_m + Q(m, lap_m) + G(grad_m, grad_m)
with R^2 = 0.999.

Now extract the EXACT coefficient tensors and determine whether the
Penrose-transformed equation is Navier-Stokes.

The NS vorticity equation is:
  d_omega/dt = nu * nabla^2 omega + (omega . grad)u - (u . grad)omega

In component form (using omega_i, u_i on R^3):
  d_omega_i/dt = nu * nabla^2 omega_i + eps_{ijk} (omega_j * partial_k u_l - u_k * partial_l omega_j)

The key structural question: does the Q tensor (m_i * lap_m_j coupling)
have the antisymmetric structure eps_{ijk} needed for the cross product?

Strategy:
1. Extract the full Q_{ij} coefficient matrix from the regression
2. Decompose Q into symmetric + antisymmetric parts
3. Check if the antisymmetric part matches eps_{abc} (Levi-Civita)
4. Extract the G_{ijkl} tensor from gradient cross terms
5. Map everything through the Pauli embedding to vorticity components
6. Compare with NS
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
print(f"m* = {m_star}")
print(f"e* = {e_star}")
sys.stdout.flush()

# ============================================================
# Set up lattice
# ============================================================

N = 16; amp = 0.20
print(f"\nSetting up {N}^3 lattice, amp={amp}...")
sys.stdout.flush()

lattice = np.zeros((N,N,N,4))
for ix,iy,iz in iprod(range(N), repeat=3):
    x=2*np.pi*ix/N; y=2*np.pi*iy/N; z=2*np.pi*iz/N
    perturb = amp * np.array([np.sin(z)+np.cos(y), np.sin(x)+np.cos(z),
                              np.sin(y)+np.cos(x), 0.0])
    m = np.abs(m_star + perturb); m /= np.sum(m)
    lattice[ix,iy,iz] = enforce_floor(m)

# ============================================================
# Compute DS spatial evolution
# ============================================================

print("Computing DS spatial evolution..."); sys.stdout.flush()

dm_spatial = np.zeros_like(lattice)
for ix,iy,iz in iprod(range(N), repeat=3):
    m_self, _ = ds_combine(lattice[ix,iy,iz], lattice[ix,iy,iz])
    e_loc = np.zeros(4)
    for d in range(3):
        for s in [+1,-1]:
            idx=[ix,iy,iz]; idx[d]=(idx[d]+s)%N
            e_loc += lattice[tuple(idx)]
    e_loc /= 6.0; e_loc = enforce_floor(e_loc)
    m_full, _ = ds_combine(lattice[ix,iy,iz], e_loc)
    dm_spatial[ix,iy,iz] = m_full - m_self

# Spatial features
lap_m = np.zeros_like(lattice)
for d in range(3):
    for ix,iy,iz in iprod(range(N), repeat=3):
        ip=[ix,iy,iz]; ip[d]=(ip[d]+1)%N
        im=[ix,iy,iz]; im[d]=(im[d]-1)%N
        lap_m[ix,iy,iz] += lattice[tuple(ip)] + lattice[tuple(im)] - 2*lattice[ix,iy,iz]

grad_m = np.zeros((N,N,N,3,4))
for d in range(3):
    for ix,iy,iz in iprod(range(N), repeat=3):
        ip=[ix,iy,iz]; ip[d]=(ip[d]+1)%N
        im=[ix,iy,iz]; im[d]=(im[d]-1)%N
        grad_m[ix,iy,iz,d,:] = (lattice[tuple(ip)] - lattice[tuple(im)]) / 2

print(f"||dm_spatial|| = {norm(dm_spatial):.6f}")
print(f"||lap_m||      = {norm(lap_m):.6f}")
sys.stdout.flush()

# ============================================================
# Extract the Q tensor: dm_spatial_i = J_ij lap_m_j + Q_ijk m_j lap_m_k
# ============================================================

print("\n" + "="*60)
print("COEFFICIENT TENSOR EXTRACTION")
print("="*60)
sys.stdout.flush()

# For each output component i, regress against:
# - lap_m_j (4 linear features)
# - m_j * lap_m_k (16 quadratic features)
# - grad_m_{d,j} * grad_m_{d,k} summed over d (16 gradient cross features)

npts = N**3
DM = dm_spatial.reshape(npts, 4)
LAP = lap_m.reshape(npts, 4)
M = lattice.reshape(npts, 4)
G = grad_m.reshape(npts, 3, 4)

# Build feature matrices
linear_features = LAP  # shape (npts, 4)

quad_features = np.zeros((npts, 16))  # m_j * lap_m_k
for j in range(4):
    for k in range(4):
        quad_features[:, 4*j+k] = M[:,j] * LAP[:,k]

grad_features = np.zeros((npts, 16))  # sum_d grad_m_{d,j} * grad_m_{d,k}
for j in range(4):
    for k in range(4):
        grad_features[:, 4*j+k] = np.sum(G[:,:,j] * G[:,:,k], axis=1)

all_features = np.column_stack([linear_features, quad_features, grad_features])

# Extract coefficients for each output component
J = np.zeros((4,4))       # linear: J[i,j] = coeff of lap_m_j in dm_i
Q = np.zeros((4,4,4))     # quadratic: Q[i,j,k] = coeff of m_j*lap_m_k in dm_i
Gcoeff = np.zeros((4,4,4))  # gradient: Gcoeff[i,j,k] = coeff of sum_d(grad_j * grad_k) in dm_i

for i in range(4):
    c, _, _, _ = lstsq(all_features, DM[:,i], rcond=None)
    J[i,:] = c[:4]
    Q[i,:,:] = c[4:20].reshape(4,4)
    Gcoeff[i,:,:] = c[20:36].reshape(4,4)

# Verify R^2
pred_all = all_features @ np.zeros(36)  # placeholder
R2_comp = []
for i in range(4):
    c, _, _, _ = lstsq(all_features, DM[:,i], rcond=None)
    pred = all_features @ c
    ss_res = np.sum((DM[:,i] - pred)**2)
    ss_tot = np.sum((DM[:,i] - np.mean(DM[:,i]))**2)
    R2_comp.append(1 - ss_res/ss_tot if ss_tot > 1e-30 else 0)

print(f"\nR^2 per component: {[f'{r:.4f}' for r in R2_comp]}")
print(f"Mean R^2: {np.mean(R2_comp):.6f}")
print()

# ============================================================
# Analyze the J matrix (linear diffusion)
# ============================================================

print("="*60)
print("LINEAR TERM: J_eff (diffusion coefficient matrix)")
print("="*60)

print("\nJ_eff (dm_i = J_ij * lap_m_j):")
for i in range(4):
    print(f"  [{J[i,0]:+.8f} {J[i,1]:+.8f} {J[i,2]:+.8f} {J[i,3]:+.8f}]")

# Decompose J into components
# For scalar diffusion: J = nu * I (identity times viscosity)
# For anisotropic: J is diagonal but not proportional to I
print(f"\nDiagonal elements: {[f'{J[i,i]:.6f}' for i in range(4)]}")
print(f"Off-diagonal max:  {max(abs(J[i,j]) for i in range(4) for j in range(4) if i!=j):.6f}")
print()

# ============================================================
# Analyze the Q tensor (m * lap_m coupling)
# ============================================================

print("="*60)
print("QUADRATIC TERM: Q tensor (m_j * lap_m_k coupling)")
print("="*60)

# The Q tensor Q[i,j,k] gives dm_i from m_j * lap_m_k.
# For NS, this should encode the advection-stretching operator.
#
# In the Pauli embedding, m = (s1, s2, s3, theta) maps to M = (theta*I + s.sigma)/sqrt(2).
# The singleton components s1, s2, s3 map to the su(2) generators sigma_1, sigma_2, sigma_3.
# The theta component maps to the identity I.
#
# The cross product structure eps_{abc} requires:
# Q[a, b, c] ~ eps_{abc} for a,b,c in {0,1,2} (singleton indices)
# That is: Q[0,1,2] = -Q[0,2,1] = Q[1,2,0] = -Q[1,0,2] = Q[2,0,1] = -Q[2,1,0]
# and all other singleton-only Q components = 0.

print("\nQ tensor: singleton block Q[a,b,c] for a,b,c in {0,1,2}:")
print("(This should be proportional to eps_{abc} for NS)")
print()

# Extract singleton-singleton-singleton block
Q_sss = Q[:3,:3,:3]

# Decompose into symmetric and antisymmetric parts
Q_sym = np.zeros((3,3,3))
Q_anti = np.zeros((3,3,3))
for i in range(3):
    for j in range(3):
        for k in range(3):
            Q_sym[i,j,k] = (Q_sss[i,j,k] + Q_sss[i,k,j]) / 2
            Q_anti[i,j,k] = (Q_sss[i,j,k] - Q_sss[i,k,j]) / 2

eps = np.zeros((3,3,3))
eps[0,1,2] = eps[1,2,0] = eps[2,0,1] = 1
eps[0,2,1] = eps[2,1,0] = eps[1,0,2] = -1

# Project Q_anti onto eps tensor
# Q_anti[i,j,k] = lambda * eps[i,j,k] + remainder
# lambda = sum Q_anti * eps / sum eps * eps = sum Q_anti * eps / 6
lam = np.sum(Q_anti * eps) / 6
Q_eps = lam * eps
Q_remainder = Q_anti - Q_eps

print(f"Antisymmetric part of Q[a,b,c] (singleton block):")
print(f"  Projection onto eps_abc: lambda = {lam:.8f}")
print(f"  ||Q_anti||    = {norm(Q_anti):.8f}")
print(f"  ||lam * eps|| = {norm(Q_eps):.8f}")
print(f"  ||remainder|| = {norm(Q_remainder):.8f}")
print(f"  eps fraction:   {norm(Q_eps)**2 / norm(Q_anti)**2 * 100:.1f}%")
print()

print(f"Symmetric part of Q[a,b,c]:")
print(f"  ||Q_sym|| = {norm(Q_sym):.8f}")
print()

print(f"Full Q singleton block norms:")
print(f"  ||Q_sss||  = {norm(Q_sss):.8f}")
print(f"  ||Q_anti|| = {norm(Q_anti):.8f} ({norm(Q_anti)/norm(Q_sss)*100:.1f}%)")
print(f"  ||Q_sym||  = {norm(Q_sym):.8f} ({norm(Q_sym)/norm(Q_sss)*100:.1f}%)")
print()

# Show the actual Q_anti components
print("Q_anti[i,j,k] (should be ~ lambda * eps_{ijk}):")
for i in range(3):
    for j in range(3):
        for k in range(3):
            if abs(Q_anti[i,j,k]) > 1e-8:
                expected = lam * eps[i,j,k]
                print(f"  Q_anti[{i},{j},{k}] = {Q_anti[i,j,k]:+.8f}  (expected {expected:+.8f})")

print()

# ============================================================
# Theta coupling: how does theta participate?
# ============================================================

print("="*60)
print("THETA COUPLING: Q[i, 3, k] and Q[i, j, 3]")
print("="*60)
print()

print("Q[i, 3, k] (theta * lap_m_k contribution to dm_i):")
for i in range(4):
    print(f"  i={i}: [{Q[i,3,0]:+.6f} {Q[i,3,1]:+.6f} {Q[i,3,2]:+.6f} {Q[i,3,3]:+.6f}]")

print(f"\n||Q[:,:3,:3]|| (singleton block) = {norm(Q[:3,:3,:3]):.6f}")
print(f"||Q[:, 3,:]|| (theta-mixed)     = {norm(Q[:, 3,:]):.6f}")
print(f"||Q[:,:, 3]|| (theta-mixed)     = {norm(Q[:,:, 3]):.6f}")
print(f"||Q[3,:,:]|| (theta output)     = {norm(Q[3,:,:]):.6f}")
print()

# ============================================================
# Gradient tensor analysis
# ============================================================

print("="*60)
print("GRADIENT TERM: G tensor (grad_m . grad_m coupling)")
print("="*60)
print()

G_sss = Gcoeff[:3,:3,:3]
print(f"||G_sss|| (singleton block) = {norm(G_sss):.6f}")
print(f"||G_full||                  = {norm(Gcoeff):.6f}")
print(f"Gradient/Quadratic ratio:   = {norm(Gcoeff)/norm(Q):.4f}")
print()

# ============================================================
# Summary: what equation did we find?
# ============================================================

print("="*60)
print("SUMMARY: THE DESCENDED EQUATION")
print("="*60)
print()
print("dm/dt = J * nabla^2(m) + Q(m, nabla^2 m) + G(grad m, grad m)")
print()
print("where:")
print(f"  J = diffusion matrix (approximately diagonal)")
print(f"    diagonal: [{J[0,0]:.6f}, {J[1,1]:.6f}, {J[2,2]:.6f}, {J[3,3]:.6f}]")
print()
print(f"  Q = quadratic coupling:")
print(f"    singleton block (s_a * lap_s_b -> ds_c):")
print(f"      antisymmetric part: {norm(Q_anti)/norm(Q_sss)*100:.1f}% of total")
print(f"      eps_{'{abc}'} coefficient: lambda = {lam:.6f}")
if norm(Q_anti) > 0:
    print(f"      eps purity: {norm(Q_eps)/norm(Q_anti)*100:.1f}%")
print()
print(f"  G = gradient cross coupling:")
print(f"    total magnitude: {norm(Gcoeff):.6f}")
print(f"    relative to Q: {norm(Gcoeff)/norm(Q)*100:.1f}%")
print()

# The NS identification check:
# If the singleton block of Q is dominated by the eps_{abc} structure,
# then after the Penrose transform (which maps s_i to vorticity/velocity
# components), the equation has the advection-stretching form.
#
# The Biot-Savart relation u = (-Delta)^{-1}(curl omega) converts
# m_j * lap_m_k into omega_j * u_k, where the inverse Laplacian
# absorbs the "lap" factor.

if norm(Q_anti) > 0.3 * norm(Q_sss):
    print("FINDING: The quadratic coupling has significant antisymmetric content.")
    print(f"The eps_abc component accounts for {norm(Q_eps)/norm(Q_anti)*100:.1f}% of the antisymmetric part.")
    if norm(Q_eps) > 0.5 * norm(Q_anti):
        print("The antisymmetric part is dominated by the Levi-Civita structure.")
        print("This is the cross-product structure of Navier-Stokes advection.")
        print()
        print(f"After Penrose transform (s_i -> omega_i, with Biot-Savart u = curl^-1(omega)):")
        print(f"  dm_a/dt ~ ... + lambda * eps_abc * omega_b * u_c")
        print(f"          = ... + lambda * (omega x u)_a")
        print(f"          = ... + lambda * [(omega.grad)u - (u.grad)omega]_a")
        print(f"  with lambda = {lam:.6f}")
    else:
        print("But the eps_abc component does not dominate — the coupling is more complex.")
else:
    print("FINDING: The quadratic coupling is predominantly symmetric.")
    print("The Levi-Civita structure is subdominant.")
    print("The descended equation may not be standard NS.")

print()
print("DONE.")
