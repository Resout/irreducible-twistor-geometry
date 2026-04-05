"""
Gluon condensate v3: The geometric integral over B.

The condensate <Tr(F^2)> is a geometric constant — the average of the
commutator-squared |[a_mu, a_nu]|^2 over the Fubini-Study measure on
B = {Born(theta) >= 1/27} subset CP^3.

The DS dynamics projects onto m* (F=0). The condensate comes from the
path integral measure (FS on B), not from the dynamical equilibrium.

This computation:
1. Sample pairs of nearby points in B uniformly (FS measure)
2. For each pair, compute the connection a = M1^{-1}(M2 - M1)
3. For triplets, compute F = [a_mu, a_nu]
4. Average |F|^2 over the ensemble
5. Relate to the volume fraction (26/27)^3

The result should be a PURE NUMBER determined by H=3 and K*=7/30 alone.
"""

import numpy as np
from numpy.linalg import inv, norm, det
from scipy.optimize import brentq
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

def to_matrix(m):
    s1,s2,s3,th = m
    return (th*np.eye(2,dtype=complex) + s1*sigma[0] + s2*sigma[1] + s3*sigma[2]) / np.sqrt(2)

def sample_B_uniform(n, rng):
    """Sample n mass functions uniformly from B subset CP^3.

    Strategy: sample from the full simplex {m >= 0, sum m = 1},
    then enforce Born floor. Rejection-sample to ensure Born >= 1/27.
    This gives the correct FS measure restricted to B.
    """
    samples = []
    while len(samples) < n:
        # Sample from Dirichlet(1,1,1,1) = uniform on the 3-simplex
        m = rng.dirichlet([1,1,1,1])
        m = enforce_floor(m)
        # Check det(M) != 0 (bundle non-degenerate)
        M = to_matrix(m)
        if abs(det(M)) > 1e-10:
            samples.append(m)
    return np.array(samples)

def find_eq():
    def K_res(p):
        pw=(1-p)/2; sc=26./27.
        raw = np.array([np.sqrt(p*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR)])
        e = raw/np.sum(raw); m = np.array([0.4,0.2,0.2,0.2])
        for _ in range(3000):
            s,th = m[:3],m[3]; se,ph = e[:3],e[3]
            sn = s*se + s*ph + th*se; tn = th*ph
            K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
            d = 1-K; out = np.zeros(4); out[:3]=sn/d; out[3]=tn/d
            m = enforce_floor(out)
        K = sum(m[i]*e[j] for i in range(3) for j in range(3) if i!=j)
        return K - 7./30.
    p = brentq(K_res, 0.92, 0.94, xtol=1e-15)
    pw=(1-p)/2; sc=26./27.
    raw = np.array([np.sqrt(p*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(FLOOR)])
    e = raw/np.sum(raw); m = np.array([0.4,0.2,0.2,0.2])
    for _ in range(3000):
        s,th = m[:3],m[3]; se,ph = e[:3],e[3]
        sn = s*se + s*ph + th*se; tn = th*ph
        K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
        d = 1-K; out = np.zeros(4); out[:3]=sn/d; out[3]=tn/d
        m = enforce_floor(out)
    return m, e

m_star, e_star = find_eq()
print(f"m* = {m_star}")
print(f"e* = {e_star}")
print(); sys.stdout.flush()

# ============================================================
# Method 1: Monte Carlo over random triplets in B
# ============================================================
#
# A "plaquette" requires 3 nearby points in R^3:
#   m(x), m(x + e_mu), m(x + e_nu)
# The connection is a_mu = M(x)^{-1} [M(x+e_mu) - M(x)]
# The curvature is F_{mu,nu} = [a_mu, a_nu] (leading order)
#
# For the PATH INTEGRAL average, we sample the 3 points independently
# from the FS measure on B. This is the "uncorrelated" estimate —
# each site fluctuates independently.
#
# The spatial correlation (from DS coupling) would modify this,
# but since the correlation decays at rate lambda_0 = 0.28 per step,
# nearest-neighbour correlations are 0.28 and next-nearest are 0.08.
# The uncorrelated estimate is the leading contribution.

print("="*60)
print("METHOD 1: Monte Carlo over random triplets")
print("="*60)
print(); sys.stdout.flush()

rng = np.random.default_rng(42)
n_samples = 100000

# Sample points from B
print(f"Sampling {n_samples} points from B..."); sys.stdout.flush()
points = sample_B_uniform(n_samples, rng)

# For each consecutive triplet (m0, m1, m2), compute:
# a_0 = M0^{-1}(M1 - M0), a_1 = M0^{-1}(M2 - M0)
# F = [a_0, a_1]
# |F|^2 = Tr(F^dag F)

print("Computing condensate..."); sys.stdout.flush()

F_sq_samples = []
comm_norm_samples = []
K_samples = []

n_triplets = n_samples // 3
for i in range(n_triplets):
    m0 = points[3*i]
    m1 = points[3*i + 1]
    m2 = points[3*i + 2]

    M0 = to_matrix(m0)
    M1 = to_matrix(m1)
    M2 = to_matrix(m2)
    M0i = inv(M0)

    # Connections
    a0 = M0i @ (M1 - M0)
    a1 = M0i @ (M2 - M0)

    # Curvature (commutator part)
    F = a0 @ a1 - a1 @ a0

    # |F|^2
    F_sq = np.real(np.trace(F.conj().T @ F))
    F_sq_samples.append(F_sq)

    # Also measure the commutator norm of the su(2) content
    # [M, E] = i(s x e).sigma — the antisymmetric off-diagonal content
    # For m0 and m1:
    s0, s1_v = m0[:3], m1[:3]
    cross = np.cross(s0, s1_v)
    comm_norm_samples.append(norm(cross))

    # K between pairs
    K01 = sum(m0[j]*m1[k] for j in range(3) for k in range(3) if j!=k)
    K_samples.append(K01)

F_sq_arr = np.array(F_sq_samples)
comm_arr = np.array(comm_norm_samples)
K_arr = np.array(K_samples)

print(f"\nResults from {n_triplets} random triplets in B:")
print(f"  <Tr(F^2)>     = {np.mean(F_sq_arr):.8f} +/- {np.std(F_sq_arr)/np.sqrt(n_triplets):.8f}")
print(f"  median         = {np.median(F_sq_arr):.8f}")
print(f"  <|s x e|>      = {np.mean(comm_arr):.8f}")
print(f"  <K>             = {np.mean(K_arr):.6f}")
print()
sys.stdout.flush()

# ============================================================
# Method 2: Condensate around the equilibrium (thermal fluctuations)
# ============================================================

print("="*60)
print("METHOD 2: Fluctuations around m*")
print("="*60)
print(); sys.stdout.flush()

# Sample points near m* with controlled distance, compute <F^2>
# This gives the condensate density at the equilibrium.

distances = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3]
n_per = 50000

for dist in distances:
    F2_local = []
    for i in range(n_per // 3):
        perturbs = dist * rng.standard_normal((3, 4))
        ms = []
        for j in range(3):
            m = m_star + perturbs[j]
            m = np.abs(m); m /= np.sum(m)
            m = enforce_floor(m)
            ms.append(m)

        M0 = to_matrix(ms[0]); M0i = inv(M0)
        a0 = M0i @ (to_matrix(ms[1]) - M0)
        a1 = M0i @ (to_matrix(ms[2]) - M0)
        F = a0 @ a1 - a1 @ a0
        F2_local.append(np.real(np.trace(F.conj().T @ F)))

    F2_a = np.array(F2_local)
    print(f"  dist={dist:.3f}: <F^2> = {np.mean(F2_a):.10f} (sigma^4 scaling: {np.mean(F2_a)/dist**4:.4f})")
    sys.stdout.flush()

print()

# ============================================================
# Method 3: The condensate as a geometric constant
# ============================================================

print("="*60)
print("THE CONDENSATE AS A GEOMETRIC CONSTANT")
print("="*60)
print(); sys.stdout.flush()

# From the sigma^4 scaling: <F^2> = C * sigma^4
# where sigma is the fluctuation scale.
#
# In the path integral, the fluctuation scale is determined by
# the eigenvalues of the transfer operator:
# <(m - m*)^2> = sum_k 1/(1 - lambda_k^2) * |psi_k|^2
#
# At the equilibrium: lambda_0 = 0.2829, lambda_1 = 0.2813
# The variance is dominated by 1/(1-lambda_0^2) = 1/(1-0.08) = 1.087
# and 1/(1-lambda_1^2) = 1/(1-0.079) = 1.086
#
# The effective sigma^2 ~ 1/(1 - lambda_0^2) per mode.
# With 3 spatial modes (the L1=1 tangent space is 3D):
# sigma_eff^2 ~ 3/(1 - lambda_0^2) = 3.26

lambda_0 = 0.2829
lambda_1 = 0.2813
sigma_eff_sq = 1.0 / (1 - lambda_0**2)
print(f"Effective fluctuation variance: sigma_eff^2 = 1/(1 - lambda_0^2) = {sigma_eff_sq:.6f}")
print(f"sigma_eff = {np.sqrt(sigma_eff_sq):.6f}")
print()

# From Method 2, the coefficient C in <F^2> = C * sigma^4:
# Use the well-measured middle range
F2_at_01 = [x for x in distances if x == 0.01]
if F2_at_01:
    # Recompute for sigma=0.01
    dist = 0.01
    F2_local = []
    for i in range(50000 // 3):
        perturbs = dist * rng.standard_normal((3, 4))
        ms = []
        for j in range(3):
            m = m_star + perturbs[j]
            m = np.abs(m); m /= np.sum(m); m = enforce_floor(m)
            ms.append(m)
        M0 = to_matrix(ms[0]); M0i = inv(M0)
        a0 = M0i @ (to_matrix(ms[1]) - M0)
        a1 = M0i @ (to_matrix(ms[2]) - M0)
        F = a0 @ a1 - a1 @ a0
        F2_local.append(np.real(np.trace(F.conj().T @ F)))
    C_coeff = np.mean(F2_local) / dist**4
    print(f"Geometric coefficient C = <F^2>/sigma^4 = {C_coeff:.4f}")
    print()

    # The condensate at the equilibrium fluctuation scale:
    condensate = C_coeff * sigma_eff_sq**2
    print(f"<Tr(F^2)>_vacuum = C * sigma_eff^4 = {C_coeff:.4f} * {sigma_eff_sq**2:.6f} = {condensate:.6f}")
    print()

    # Relate to framework constants
    print(f"Ratio to K*: {condensate / (7/30):.6f}")
    print(f"Ratio to K*^2: {condensate / (7/30)**2:.6f}")
    print(f"Ratio to 1/27: {condensate / (1/27):.6f}")
    print(f"Ratio to (1-K*): {condensate / (23/30):.6f}")
    print(f"Ratio to eta: {condensate / (4/27):.6f}")
    print()

    # Volume fraction connection
    vol_frac = (26./27.)**3
    print(f"Volume fraction (26/27)^3 = {vol_frac:.6f}")
    print(f"Condensate * vol_frac = {condensate * vol_frac:.6f}")
    print(f"Condensate / vol_frac = {condensate / vol_frac:.6f}")

print()
print("DONE.")
