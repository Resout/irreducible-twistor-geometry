"""
NS Identification v5: Correct nonlinear term identification.

The key insight: the nonlinear term is NOT [a, dB] (static connection × gradient).
It is the cross-coupling between the CHANGE in connection (from DS spatial coupling)
and the existing curvature.

The DS step with neighbour evidence:
  E(x) = (1/6) sum_{y~x} m(y) = m(x) + (1/6) nabla^2 m(x) + ...

This produces a change in m: delta_m = Phi(m, E) - Phi(m, m)
The spatial part: delta_m_spatial ~ J_e * (1/6) nabla^2 m

In the matrix representation:
  delta_M = J_e_matrix * (1/6) nabla^2 M
  delta_a = M^{-1} d(delta_M) + [a, M^{-1} delta_M]
            ^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^
            second derivative         first × first
            (DIFFUSION)              (NONLINEAR)

The [a, M^{-1}delta_M] term: a has one spatial derivative, delta_M has one spatial
derivative (from nabla^2 M = sum of neighbours - 2M, which contains gradients).
The product is bilinear in first derivatives — the advection-stretching structure.

Actually, let me just directly measure: what fraction of dB_spatial is explained by
the FULL covariant Laplacian D^2 B = nabla^2 B + 2[a, dB] + [a,[a,B]]?

If R^2 is significantly higher than for the ordinary Laplacian, the covariant terms
are the missing physics.
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
# Instead of extracting vorticity from curvature, work directly
# with the mass functions and measure what DS does to them.
#
# The DIRECT approach: compare the 4-component delta_m_spatial
# at each site against spatial derivatives of m.
#
# delta_m_spatial(x) = Phi(m(x), E(x)) - Phi(m(x), m(x))
# where E(x) = (1/6) sum_{y~x} m(y) = m(x) + (1/6) lap_m(x) + ...
#
# To first order in E-m:
# delta_m_spatial ≈ J_e(m) * (1/6) lap_m(x)
#
# But to second order, there are QUADRATIC terms from the DS nonlinearity:
# delta_m_spatial ≈ J_e * (1/6) lap_m + (1/2) H_e * ((1/6) lap_m)^2 + ...
# where H_e is the Hessian of Phi w.r.t. evidence.
#
# The quadratic term H_e * (lap_m)^2 is the nonlinear coupling.
# In the connection language, this is the [a, phi] commutator.
#
# But actually: (1/6) lap_m = E - m = (1/6) sum (m(y) - m(x))
# The individual terms m(y) - m(x) are first-order in the gradient.
# Their squares couple to the Hessian.
#
# More precisely: E = m + (a^2/6)(d^2m/dx^2) = m + (a^2/6) nabla^2 m
# But this misses the CROSS terms between different directions.
# E = (1/6)(m(x+1,y,z) + m(x-1,y,z) + m(x,y+1,z) + ... )
# = m + (a^2/6) nabla^2 m + (a^4/360)(nabla^4 m) + ...
# But also: the DS combination is nonlinear, so:
# Phi(m, m+eps) = Phi(m,m) + J_e * eps + (1/2) * eps^T H_e eps + ...
# where eps = (1/6) nabla^2 m.
#
# The quadratic part: (1/2) (nabla^2 m)^T H_e (nabla^2 m) / 36
# This is O(amp^2 / N^4) — very small on a smooth lattice.
#
# Hmm. The DS nonlinearity couples DIFFERENT COMPONENTS of nabla^2 m.
# Let me just compute everything directly and see.
# ============================================================

print("\nDIRECT APPROACH: Measure DS spatial coupling structure")
print("="*60); sys.stdout.flush()

# Compute at each site:
# (a) delta_m_spatial = Phi(m, E_neighbour) - Phi(m, m)
# (b) nabla^2 m (discrete Laplacian)
# (c) directional second derivatives d^2m/dx_j^2

lap_m = np.zeros_like(lattice)
for d in range(3):
    for ix,iy,iz in iprod(range(N), repeat=3):
        ip=[ix,iy,iz]; ip[d]=(ip[d]+1)%N
        im=[ix,iy,iz]; im[d]=(im[d]-1)%N
        lap_m[ix,iy,iz] += lattice[tuple(ip)] + lattice[tuple(im)] - 2*lattice[ix,iy,iz]

# Self-evidence and full evolution
dm_spatial = np.zeros_like(lattice)
dm_self = np.zeros_like(lattice)
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
    dm_self[ix,iy,iz] = m_self - lattice[ix,iy,iz]

print(f"||dm_self||    = {norm(dm_self):.8f}")
print(f"||dm_spatial|| = {norm(dm_spatial):.8f}")
print(f"||lap_m||      = {norm(lap_m):.8f}")
print()

# Linear regression: dm_spatial = J * lap_m (4x4 effective Jacobian)
# At each site: dm_spatial_i = sum_j J_ij * lap_m_j
# Global fit: stack all sites

DM = dm_spatial.reshape(-1, 4)  # (N^3, 4)
LAP = lap_m.reshape(-1, 4)     # (N^3, 4)

# Fit each output component against all input components
J_eff = np.zeros((4,4))
R2_per_comp = []
for i in range(4):
    c, _, _, _ = lstsq(LAP, DM[:,i], rcond=None)
    J_eff[i,:] = c
    pred = LAP @ c
    ss_res = np.sum((DM[:,i] - pred)**2)
    ss_tot = np.sum((DM[:,i] - np.mean(DM[:,i]))**2)
    R2_per_comp.append(1 - ss_res/ss_tot if ss_tot > 1e-30 else 0)

print("Effective spatial Jacobian: dm_spatial = J_eff * lap_m")
print("J_eff (4x4):")
for i in range(4):
    print(f"  [{J_eff[i,0]:+.6f} {J_eff[i,1]:+.6f} {J_eff[i,2]:+.6f} {J_eff[i,3]:+.6f}]  R^2={R2_per_comp[i]:.4f}")
print(f"\nOverall R^2: {np.mean(R2_per_comp):.4f}")

# Now add QUADRATIC terms: cross products s_i * lap_s_j
# These represent the nonlinear coupling
print("\nAdding quadratic terms (cross-products of m components with Laplacian components)...")
sys.stdout.flush()

# Build quadratic features: m_i * lap_m_j for all (i,j) pairs = 16 features
M_vals = lattice.reshape(-1, 4)
L_vals = lap_m.reshape(-1, 4)
quad_features = np.zeros((N**3, 16))
for i in range(4):
    for j in range(4):
        quad_features[:, 4*i+j] = M_vals[:, i] * L_vals[:, j]

# Also add the gradient cross-terms: (dm/dx_mu)_i * (dm/dx_mu)_j
# These are the a * dB terms
grad_m = np.zeros((N,N,N,3,4))  # gradient of m in direction d
for d in range(3):
    for ix,iy,iz in iprod(range(N), repeat=3):
        ip=[ix,iy,iz]; ip[d]=(ip[d]+1)%N
        im=[ix,iy,iz]; im[d]=(im[d]-1)%N
        grad_m[ix,iy,iz,d,:] = (lattice[tuple(ip)] - lattice[tuple(im)]) / 2

# Cross-gradient features: sum_d (dm/dx_d)_i * (dm/dx_d)_j
grad_cross = np.zeros((N**3, 16))
G = grad_m.reshape(-1, 3, 4)
for i in range(4):
    for j in range(4):
        grad_cross[:, 4*i+j] = np.sum(G[:,:,i] * G[:,:,j], axis=1)

# Full regression: dm_spatial = J_eff * lap_m + Q1 * (m * lap_m) + Q2 * (grad * grad)
all_features = np.column_stack([LAP, quad_features, grad_cross])
print(f"Feature matrix: {all_features.shape[0]} samples x {all_features.shape[1]} features")

R2_full = []
for i in range(4):
    c, _, _, _ = lstsq(all_features, DM[:,i], rcond=None)
    pred = all_features @ c
    ss_res = np.sum((DM[:,i] - pred)**2)
    ss_tot = np.sum((DM[:,i] - np.mean(DM[:,i]))**2)
    R2_full.append(1 - ss_res/ss_tot if ss_tot > 1e-30 else 0)

print(f"\nR^2 with quadratic terms:")
for i in range(4):
    print(f"  Component {i}: linear R^2={R2_per_comp[i]:.4f}, full R^2={R2_full[i]:.4f}, improvement={R2_full[i]-R2_per_comp[i]:.4f}")

print(f"\nOverall: linear {np.mean(R2_per_comp):.4f}, full {np.mean(R2_full):.4f}")
print(f"Quadratic improvement: {np.mean(R2_full) - np.mean(R2_per_comp):.4f}")
print()

if np.mean(R2_full) - np.mean(R2_per_comp) > 0.01:
    print("SIGNIFICANT: Quadratic/cross-product terms improve the fit.")
    print("The nonlinear coupling is detectable in the mass function evolution.")
    print("This is the precursor of the NS advection-stretching term.")

    # Which quadratic terms dominate?
    # Refit with just linear + grad_cross (the a*dB terms)
    features_grad = np.column_stack([LAP, grad_cross])
    R2_grad = []
    for i in range(4):
        c, _, _, _ = lstsq(features_grad, DM[:,i], rcond=None)
        pred = features_grad @ c
        ss_res = np.sum((DM[:,i] - pred)**2)
        ss_tot = np.sum((DM[:,i] - np.mean(DM[:,i]))**2)
        R2_grad.append(1 - ss_res/ss_tot if ss_tot > 1e-30 else 0)

    features_quad = np.column_stack([LAP, quad_features])
    R2_quad = []
    for i in range(4):
        c, _, _, _ = lstsq(features_quad, DM[:,i], rcond=None)
        pred = features_quad @ c
        ss_res = np.sum((DM[:,i] - pred)**2)
        ss_tot = np.sum((DM[:,i] - np.mean(DM[:,i]))**2)
        R2_quad.append(1 - ss_res/ss_tot if ss_tot > 1e-30 else 0)

    print(f"\n  Linear only:           {np.mean(R2_per_comp):.4f}")
    print(f"  + gradient cross:      {np.mean(R2_grad):.4f} (improvement {np.mean(R2_grad)-np.mean(R2_per_comp):.4f})")
    print(f"  + m*lap_m cross:       {np.mean(R2_quad):.4f} (improvement {np.mean(R2_quad)-np.mean(R2_per_comp):.4f})")
    print(f"  + both:                {np.mean(R2_full):.4f} (improvement {np.mean(R2_full)-np.mean(R2_per_comp):.4f})")
else:
    print("Quadratic terms do not significantly improve the fit.")
    print("The spatial evolution is well-described by linear diffusion alone")
    print("at this perturbation amplitude.")

print("\nDONE.")
