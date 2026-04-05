"""
Gluon condensate <F^2> from the DS framework.

The equilibrium has F ≈ 0 (rank-1 dbar Phi, abelian connection).
But the path integral average <F^2> over fluctuations around m* should
be nonzero — this is the gluon condensate.

Strategy:
1. Sample mass functions from the Fubini-Study measure on B = {Born >= 1/27}
2. At each sample m, compute:
   - The transition function M(m)
   - The connection a = M^{-1} dM (using a small perturbation as "spatial gradient")
   - The curvature F = da + a ^ a
   - |F|^2 = Tr(F^dag F)
3. Average over samples: <F^2> = (1/N) sum |F_i|^2
4. Also measure the contribution from the a^a commutator term vs the da term

The fluctuations are sampled from the Fubini-Study measure restricted to B.
On B, the natural measure is the FS volume form.
"""

import numpy as np
from numpy.linalg import inv, norm, det
from itertools import product as iprod
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
M_star = to_matrix(m_star)
print(f"m* = {m_star}")
print(f"det(M*) = {det(M_star):.6f}")
print(); sys.stdout.flush()

# ============================================================
# Sampling strategy
# ============================================================
#
# The path integral measure is the Fubini-Study measure on B.
# For the condensate, we need to compute <F^2> where F is the
# curvature at each configuration.
#
# A "configuration" in the DS framework is a mass function m(x)
# at each spacetime point. The connection is a = M^{-1} dM,
# which requires spatial variation. At a SINGLE point, a = 0.
#
# The condensate comes from the COUPLED system: a lattice of
# mass functions, each interacting with neighbours via DS combination.
# F^2 at a site depends on the mass function AND its neighbours.
#
# Approach: set up a small lattice at thermal equilibrium (after
# many DS steps with the K*=7/30 evidence), then measure F^2.
# The "thermal" fluctuations around m* are the quantum fluctuations.

# ============================================================
# Phase 1: Measure <F^2> at the DS equilibrium on a lattice
# ============================================================

print("="*60)
print("PHASE 1: <F^2> on a thermalized lattice")
print("="*60)
print(); sys.stdout.flush()

N = 10  # lattice size

# Initialize near equilibrium with small random perturbations
np.random.seed(42)
lattice = np.zeros((N,N,N,4))
for ix,iy,iz in iprod(range(N), repeat=3):
    perturb = 0.01 * np.random.randn(4)
    perturb[3] = 0  # don't perturb theta directly
    m = m_star + perturb
    m = np.abs(m)
    m /= np.sum(m)
    lattice[ix,iy,iz] = enforce_floor(m)

# Thermalize: run DS with K*=7/30 evidence for many steps
print(f"Thermalizing {N}^3 lattice ({N**3} sites)...")
sys.stdout.flush()

n_therm = 50
for step in range(n_therm):
    lattice_new = np.zeros_like(lattice)
    for ix,iy,iz in iprod(range(N), repeat=3):
        # Evidence: average of neighbours
        e_loc = np.zeros(4)
        for d in range(3):
            for s in [+1,-1]:
                idx=[ix,iy,iz]; idx[d]=(idx[d]+s)%N
                e_loc += lattice[tuple(idx)]
        e_loc /= 6.0
        e_loc = enforce_floor(e_loc)
        lattice_new[ix,iy,iz], _ = ds_combine(lattice[ix,iy,iz], e_loc)
    lattice = lattice_new

    if (step+1) % 10 == 0:
        # Check convergence: max deviation from m*
        dev = np.max([norm(lattice[ix,iy,iz] - m_star) for ix,iy,iz in iprod(range(N), repeat=3)])
        print(f"  Step {step+1}: max|m - m*| = {dev:.6f}")
        sys.stdout.flush()

print("Thermalization complete.")
print(); sys.stdout.flush()

# ============================================================
# Phase 2: Compute F^2 at each site
# ============================================================

print("Computing F^2 at each lattice site...")
sys.stdout.flush()

def compute_F_squared(lattice, N):
    """Compute |F|^2 = sum_{mu<nu} Tr(F_{mu,nu}^dag F_{mu,nu}) at each site.

    Connection: a_mu(x) = M(x)^{-1} [M(x+mu) - M(x-mu)] / 2
    Curvature: F_{mu,nu} = [a_mu, a_nu] + (discrete curl of a)
    |F|^2 = sum over plaquettes of |F_{mu,nu}|^2

    Returns: F_sq (scalar at each site), F_comm_sq (commutator part),
             F_da_sq (derivative part), plaquette details
    """
    F_sq = np.zeros((N,N,N))
    F_comm_sq = np.zeros((N,N,N))
    F_da_sq = np.zeros((N,N,N))

    for ix,iy,iz in iprod(range(N), repeat=3):
        M = to_matrix(lattice[ix,iy,iz])
        Mi = inv(M)

        # Connection at this site
        a = []
        for d in range(3):
            ip=[ix,iy,iz]; im=[ix,iy,iz]
            ip[d]=(ip[d]+1)%N; im[d]=(im[d]-1)%N
            Mp = to_matrix(lattice[tuple(ip)])
            Mm = to_matrix(lattice[tuple(im)])
            a.append(Mi @ (Mp - Mm) / 2.0)

        for mu in range(3):
            for nu in range(mu+1, 3):
                # Commutator part: [a_mu, a_nu]
                comm = a[mu] @ a[nu] - a[nu] @ a[mu]

                # Derivative part: partial_mu a_nu - partial_nu a_mu
                # partial_mu a_nu = (a_nu(x+mu) - a_nu(x-mu)) / 2
                ip=[ix,iy,iz]; im=[ix,iy,iz]
                ip[mu]=(ip[mu]+1)%N; im[mu]=(im[mu]-1)%N

                M_p = to_matrix(lattice[tuple(ip)])
                Mi_p = inv(M_p)
                ip2=list(ip); im2=list(ip)
                ip2[nu]=(ip2[nu]+1)%N; im2[nu]=(im2[nu]-1)%N
                a_nu_p = Mi_p @ (to_matrix(lattice[tuple(ip2)]) - to_matrix(lattice[tuple(im2)])) / 2

                M_m = to_matrix(lattice[tuple(im)])
                Mi_m = inv(M_m)
                ip3=list(im); im3=list(im)
                ip3[nu]=(ip3[nu]+1)%N; im3[nu]=(im3[nu]-1)%N
                a_nu_m = Mi_m @ (to_matrix(lattice[tuple(ip3)]) - to_matrix(lattice[tuple(im3)])) / 2

                da_mu_nu = (a_nu_p - a_nu_m) / 2

                ip4=[ix,iy,iz]; im4=[ix,iy,iz]
                ip4[nu]=(ip4[nu]+1)%N; im4[nu]=(im4[nu]-1)%N

                M_p2 = to_matrix(lattice[tuple(ip4)])
                Mi_p2 = inv(M_p2)
                ip5=list(ip4); im5=list(ip4)
                ip5[mu]=(ip5[mu]+1)%N; im5[mu]=(im5[mu]-1)%N
                a_mu_p2 = Mi_p2 @ (to_matrix(lattice[tuple(ip5)]) - to_matrix(lattice[tuple(im5)])) / 2

                M_m2 = to_matrix(lattice[tuple(im4)])
                Mi_m2 = inv(M_m2)
                ip6=list(im4); im6=list(im4)
                ip6[mu]=(ip6[mu]+1)%N; im6[mu]=(im6[mu]-1)%N
                a_mu_m2 = Mi_m2 @ (to_matrix(lattice[tuple(ip6)]) - to_matrix(lattice[tuple(im6)])) / 2

                da_nu_mu = (a_mu_p2 - a_mu_m2) / 2

                da = da_mu_nu - da_nu_mu

                F_mn = da + comm

                # |F|^2 = Tr(F^dag F) (real, non-negative)
                F_sq[ix,iy,iz] += np.real(np.trace(F_mn.conj().T @ F_mn))
                F_comm_sq[ix,iy,iz] += np.real(np.trace(comm.conj().T @ comm))
                F_da_sq[ix,iy,iz] += np.real(np.trace(da.conj().T @ da))

    return F_sq, F_comm_sq, F_da_sq

F_sq, F_comm_sq, F_da_sq = compute_F_squared(lattice, N)

print(f"\n<F^2> (full curvature):")
print(f"  Mean:   {np.mean(F_sq):.10f}")
print(f"  Std:    {np.std(F_sq):.10f}")
print(f"  Min:    {np.min(F_sq):.10f}")
print(f"  Max:    {np.max(F_sq):.10f}")
print()
print(f"<|[a,a]|^2> (commutator / non-abelian part):")
print(f"  Mean:   {np.mean(F_comm_sq):.10f}")
print(f"  Fraction of <F^2>: {np.mean(F_comm_sq)/np.mean(F_sq)*100:.2f}%")
print()
print(f"<|da|^2> (derivative / abelian part):")
print(f"  Mean:   {np.mean(F_da_sq):.10f}")
print(f"  Fraction of <F^2>: {np.mean(F_da_sq)/np.mean(F_sq)*100:.2f}%")
print()
sys.stdout.flush()

# ============================================================
# Phase 3: Measure at the EXACT equilibrium (should be ~0)
# ============================================================

print("="*60)
print("PHASE 2: F^2 at exact equilibrium (should be ~0)")
print("="*60)
print()

lattice_eq = np.zeros((N,N,N,4))
for ix,iy,iz in iprod(range(N), repeat=3):
    lattice_eq[ix,iy,iz] = m_star.copy()

F_sq_eq, F_comm_eq, F_da_eq = compute_F_squared(lattice_eq, N)
print(f"<F^2> at uniform equilibrium: {np.mean(F_sq_eq):.2e}")
print(f"  (Should be ~0 since all sites identical -> a = 0)")
print()
sys.stdout.flush()

# ============================================================
# Phase 4: Measure with controlled fluctuation amplitudes
# ============================================================

print("="*60)
print("PHASE 3: <F^2> vs fluctuation amplitude")
print("="*60)
print()

amplitudes = [0.001, 0.003, 0.01, 0.03, 0.05, 0.1]
F2_vs_amp = []

for amp in amplitudes:
    np.random.seed(123)
    lat = np.zeros((N,N,N,4))
    for ix,iy,iz in iprod(range(N), repeat=3):
        perturb = amp * np.random.randn(4)
        perturb[3] = 0
        m = m_star + perturb
        m = np.abs(m); m /= np.sum(m)
        lat[ix,iy,iz] = enforce_floor(m)

    # Thermalize briefly
    for step in range(20):
        lat_new = np.zeros_like(lat)
        for ix,iy,iz in iprod(range(N), repeat=3):
            e_loc = np.zeros(4)
            for d in range(3):
                for s in [+1,-1]:
                    idx=[ix,iy,iz]; idx[d]=(idx[d]+s)%N
                    e_loc += lat[tuple(idx)]
            e_loc /= 6.0; e_loc = enforce_floor(e_loc)
            lat_new[ix,iy,iz], _ = ds_combine(lat[ix,iy,iz], e_loc)
        lat = lat_new

    Fs, Fc, Fd = compute_F_squared(lat, N)
    F2_mean = np.mean(Fs)
    Fc_mean = np.mean(Fc)
    F2_vs_amp.append((amp, F2_mean, Fc_mean))
    print(f"  amp={amp:.4f}: <F^2>={F2_mean:.8f}, <[a,a]^2>={Fc_mean:.8f} ({Fc_mean/F2_mean*100:.1f}% non-abelian)")
    sys.stdout.flush()

print()

# Check scaling: does <F^2> scale as amp^2 (linear response) or amp^4 (quadratic)?
amps = np.array([x[0] for x in F2_vs_amp])
F2s = np.array([x[1] for x in F2_vs_amp])

# Fit log-log: log(F2) = alpha * log(amp) + const
mask = F2s > 1e-20
if np.sum(mask) >= 2:
    log_amps = np.log(amps[mask])
    log_F2 = np.log(F2s[mask])
    alpha, const = np.polyfit(log_amps, log_F2, 1)
    print(f"Scaling: <F^2> ~ amp^{alpha:.2f}")
    print(f"  (amp^2 = linear response / Gaussian fluctuations)")
    print(f"  (amp^4 = two-gradient coupling)")
    print()

# ============================================================
# Phase 5: The condensate in framework units
# ============================================================

print("="*60)
print("THE GLUON CONDENSATE")
print("="*60)
print()

# The thermalized <F^2> is the condensate.
# In the DS framework, everything is dimensionless (natural units where
# the mass function has L1 = 1). The condensate is a pure number.

print(f"<Tr(F^2)> at thermal equilibrium = {np.mean(F_sq):.8f}")
print(f"  from {N}^3 = {N**3} sites, {n_therm} thermalization steps")
print(f"  initial perturbation sigma = 0.01")
print()
print(f"Decomposition:")
print(f"  da part (abelian):      {np.mean(F_da_sq):.8f} ({np.mean(F_da_sq)/np.mean(F_sq)*100:.1f}%)")
print(f"  [a,a] part (non-abelian): {np.mean(F_comm_sq):.8f} ({np.mean(F_comm_sq)/np.mean(F_sq)*100:.1f}%)")
print()

# Relate to K*
# The paper shows K measures the symmetric off-diagonal content
# and F involves the antisymmetric content.
# At equilibrium K* = 7/30, but F ≈ 0.
# The fluctuation <F^2> should be related to the variance of K around K*.

# Compute K variance
K_values = []
for ix,iy,iz in iprod(range(N), repeat=3):
    m = lattice[ix,iy,iz]
    e_loc = np.zeros(4)
    for d in range(3):
        for s in [+1,-1]:
            idx=[ix,iy,iz]; idx[d]=(idx[d]+s)%N
            e_loc += lattice[tuple(idx)]
    e_loc /= 6.0; e_loc = enforce_floor(e_loc)
    K = sum(m[i]*e_loc[j] for i in range(3) for j in range(3) if i!=j)
    K_values.append(K)

K_arr = np.array(K_values)
print(f"K statistics at thermal equilibrium:")
print(f"  <K>       = {np.mean(K_arr):.6f} (K* = {7/30:.6f})")
print(f"  std(K)    = {np.std(K_arr):.6f}")
print(f"  <(K-K*)^2> = {np.mean((K_arr - 7/30)**2):.8f}")
print()

# The volume fraction
vol_frac = (26/27)**3
print(f"Volume fraction of B in CP^3: (26/27)^3 = {vol_frac:.6f}")
print(f"<F^2> * vol_frac = {np.mean(F_sq) * vol_frac:.8f}")
print(f"<F^2> / K* = {np.mean(F_sq) / (7/30):.8f}")
print(f"<F^2> / (K*)^2 = {np.mean(F_sq) / (7/30)**2:.8f}")
print()

print("DONE.")
