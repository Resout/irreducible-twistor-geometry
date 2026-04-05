"""
Koopman operator compactness argument for uniform-in-N spectral gap.

Question: On the N-site coupled system B^N, does the Koopman operator
have a spectral gap that is bounded below independently of N?

Approach:
  1. Single site: Koopman operator T on L²(B, μ_FS). B compact, Φ smooth.
     Is T compact? If so, discrete spectrum, gap = -ln|λ₀|.

  2. N sites (uncoupled): Koopman on L²(B^N) = L²(B)^⊗N.
     T_N = T⊗T⊗...⊗T. Spectral gap of tensor product = min of individual gaps.
     (Because eigenvalues of T^⊗N are products of individual eigenvalues,
      and the largest non-trivial product eigenvalue is λ₀·1·1·...·1 = λ₀.)

  3. N sites (coupled): T_N^coupled ≠ T^⊗N. The coupling modifies the operator.
     The question: does coupling change the spectral gap?

  4. Key insight: coupling is NEAREST-NEIGHBOR and the coupling Jacobian has
     range contained in T_p(S) (proven in the paper). The coupling is a
     PERTURBATION of the uncoupled operator. If the perturbation is bounded
     and the gap of the uncoupled system is Δ, then by stability of spectral
     gaps under perturbation, the coupled gap is ≥ Δ - ||perturbation||.

Let's compute this.
"""

import numpy as np
from scipy.optimize import brentq

H, FLOOR = 3, 1.0/27

def enforce_floor(m):
    s = m[:3].copy()
    lo, hi = m[3], 1.0
    for _ in range(50):
        mid = (lo+hi)/2
        ss = (1-mid)/np.sum(s)*s if np.sum(s)>0 else s
        if mid**2/(np.sum(ss**2)+mid**2) < FLOOR: lo=mid
        else: hi=mid
    tn = (lo+hi)/2
    ss = (1-tn)/np.sum(s)*s if np.sum(s)>0 else s
    return np.array([ss[0],ss[1],ss[2],tn])

def ds_combine(m, e):
    s,th = m[:3],m[3]; se,ph = e[:3],e[3]
    sn = s*se + s*ph + th*se; tn = th*ph
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
    mo = np.append(sn,tn)/(1-K)
    if mo[3]**2/np.sum(mo**2) < FLOOR: mo = enforce_floor(mo)
    return mo

# K*=7/30 equilibrium
def K_eq(pd):
    pw=(1-pd)/2; sc=1-FLOOR
    raw=np.array([np.sqrt(pd*sc),np.sqrt(pw*sc),np.sqrt(pw*sc),np.sqrt(FLOOR)])
    e=raw/np.sum(raw); m=np.array([.4,.2,.2,.2])
    for _ in range(500):
        mn=ds_combine(m,e)
        if np.max(np.abs(mn-m))<1e-15: break
        m=mn
    s,th=m[:3],m[3]; se,ph=e[:3],e[3]
    return sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)

pd=brentq(lambda p: K_eq(p)-7/30, .9, .95)
sc=1-FLOOR; pw=(1-pd)/2
raw=np.array([np.sqrt(pd*sc),np.sqrt(pw*sc),np.sqrt(pw*sc),np.sqrt(FLOOR)])
e_star=raw/np.sum(raw)
m=np.array([.4,.2,.2,.2])
for _ in range(500):
    mn=ds_combine(m,e_star)
    if np.max(np.abs(mn-m))<1e-15: break
    m=mn
m_star=mn

print(f"m* = [{m_star[0]:.6f}, {m_star[1]:.6f}, {m_star[2]:.6f}, {m_star[3]:.6f}]")

# ============================================================
# PART 1: Single-site Jacobian
# ============================================================
print("\n" + "=" * 70)
print("PART 1: SINGLE-SITE JACOBIAN")
print("=" * 70)

eps = 1e-8
def jacobian(m, e, eps=1e-8):
    J = np.zeros((4,4))
    f0 = ds_combine(m, e)
    for j in range(4):
        mp = m.copy(); mp[j] += eps
        J[:,j] = (ds_combine(mp, e) - f0) / eps
    return J

J_local = jacobian(m_star, e_star)

# Project to L1=1 tangent space
V = np.array([[1,0,0,-1],[0,1,0,-1],[0,0,1,-1]], dtype=float).T
VTV_inv = np.linalg.inv(V.T @ V)
J_local_proj = VTV_inv @ V.T @ J_local @ V

evals_local = np.sort(np.abs(np.linalg.eigvals(J_local_proj)))[::-1]
print(f"  Single-site eigenvalues: {evals_local}")
print(f"  λ₀ = {evals_local[0]:.10f}")
print(f"  Δ_local = {-np.log(evals_local[0]):.6f}")

lambda_0 = evals_local[0]

# ============================================================
# PART 2: Coupling Jacobian — derivative of DS step w.r.t. evidence
# ============================================================
print("\n" + "=" * 70)
print("PART 2: COUPLING JACOBIAN (dΦ/de at equilibrium)")
print("=" * 70)

def jacobian_wrt_evidence(m, e, eps=1e-8):
    """Derivative of DS(m,e) with respect to e."""
    J = np.zeros((4,4))
    f0 = ds_combine(m, e)
    for j in range(4):
        ep = e.copy(); ep[j] += eps
        ep = ep / np.sum(ep)  # renormalize evidence
        J[:,j] = (ds_combine(m, ep) - f0) / eps
    return J

# At equilibrium, what is the evidence? For the coupled system,
# evidence at site i = average of neighbors' states
# At coupled equilibrium, all sites = m*, so evidence = m*

# But we computed the single-site gap with evidence e* ≠ m*.
# The COUPLED system uses evidence = neighbor = m* at equilibrium.
# Let's compute both Jacobians.

# Jacobian of DS(m*, m*) w.r.t. m (self-evidence)
J_self = jacobian(m_star, m_star)
J_self_proj = VTV_inv @ V.T @ J_self @ V
evals_self = np.sort(np.abs(np.linalg.eigvals(J_self_proj)))[::-1]
print(f"  Self-evidence Jacobian eigenvalues: {evals_self}")
print(f"  λ₀(self) = {evals_self[0]:.10f}")
if evals_self[0] < 1:
    print(f"  Δ_self = {-np.log(evals_self[0]):.6f}")
else:
    print(f"  Δ_self = NO GAP (λ₀ ≥ 1)")

# Jacobian of DS(m*, e*) w.r.t. m (external evidence)
print(f"\n  External evidence Jacobian eigenvalues: {evals_local}")
print(f"  λ₀(external) = {evals_local[0]:.10f}")
print(f"  Δ_external = {-np.log(evals_local[0]):.6f}")

# K at self-combination
s = m_star[:3]; th = m_star[3]
K_self = sum(s[i]*s[j] for i in range(3) for j in range(3) if i!=j)
print(f"\n  K(m*, m*) = {K_self:.6f} (self-evidence)")
print(f"  K(m*, e*) = {7/30:.6f} (external evidence)")

# ============================================================
# PART 3: The actual coupled system at equilibrium
# ============================================================
print("\n" + "=" * 70)
print("PART 3: COUPLED SYSTEM — WHAT ACTUALLY HAPPENS")
print("=" * 70)

def coupled_equilibrium(N, n_sweeps=500):
    rng = np.random.default_rng(42 + N)
    lattice = np.zeros((N, 4))
    for i in range(N):
        m = rng.dirichlet([1,1,1,1])
        while m[3]**2/np.sum(m**2) < FLOOR:
            m = rng.dirichlet([1,1,1,1])
        lattice[i] = m
    for _ in range(n_sweeps):
        new = np.zeros_like(lattice)
        for i in range(N):
            left = lattice[(i-1)%N]; right = lattice[(i+1)%N]
            e = (left+right)/2; e = e/np.sum(e)
            new[i] = ds_combine(lattice[i], e)
        lattice = new
    return lattice

for N in [2, 4, 8]:
    eq = coupled_equilibrium(N, 500)

    # What is the equilibrium state?
    print(f"\n  N={N}: equilibrium states:")
    for i in range(min(N, 4)):
        print(f"    site {i}: [{eq[i,0]:.6f}, {eq[i,1]:.6f}, {eq[i,2]:.6f}, {eq[i,3]:.6f}]")

    # Are they all the same?
    spread = np.max([np.linalg.norm(eq[i] - eq[0]) for i in range(N)])
    print(f"    Max spread: {spread:.6e}")

    # What is K at equilibrium?
    for i in range(min(N, 2)):
        left = eq[(i-1)%N]; right = eq[(i+1)%N]
        e = (left+right)/2; e = e/np.sum(e)
        s = eq[i,:3]; se = e[:3]
        K = sum(s[ii]*se[jj] for ii in range(3) for jj in range(3) if ii!=jj)
        print(f"    K at site {i}: {K:.6f}")

    # Compute the coupled Jacobian at equilibrium
    # Perturb site 0, measure response at site 0
    J_coupled = np.zeros((4,4))
    f0 = eq.copy()
    new_f0 = np.zeros_like(f0)
    for i in range(N):
        left = f0[(i-1)%N]; right = f0[(i+1)%N]
        e = (left+right)/2; e = e/np.sum(e)
        new_f0[i] = ds_combine(f0[i], e)

    for j in range(4):
        f_pert = eq.copy()
        f_pert[0, j] += eps
        new_pert = np.zeros_like(f_pert)
        for i in range(N):
            left = f_pert[(i-1)%N]; right = f_pert[(i+1)%N]
            e = (left+right)/2; e = e/np.sum(e)
            new_pert[i] = ds_combine(f_pert[i], e)
        J_coupled[:, j] = (new_pert[0] - new_f0[0]) / eps

    J_coupled_proj = VTV_inv @ V.T @ J_coupled @ V
    evals_coupled = np.sort(np.abs(np.linalg.eigvals(J_coupled_proj)))[::-1]
    print(f"    Coupled Jacobian eigenvalues at site 0: {evals_coupled}")
    if evals_coupled[0] < 1 and evals_coupled[0] > 0:
        print(f"    Δ_coupled = {-np.log(evals_coupled[0]):.6f}")

# ============================================================
# PART 4: Key question — does the coupled equilibrium have a gap?
# ============================================================
print("\n" + "=" * 70)
print("PART 4: COUPLED GAP vs SINGLE-SITE GAP")
print("=" * 70)

print(f"\n  {'N':>4}  {'λ₀(coupled)':>14}  {'Δ(coupled)':>12}  {'Δ(single)':>12}  {'ratio':>8}")

for N in [2, 3, 4, 6, 8, 12, 16]:
    eq = coupled_equilibrium(N, 500)

    # Coupled Jacobian at site 0
    f0_full = eq.copy()
    new_f0 = np.zeros_like(f0_full)
    for i in range(N):
        left = f0_full[(i-1)%N]; right = f0_full[(i+1)%N]
        e = (left+right)/2; e = e/np.sum(e)
        new_f0[i] = ds_combine(f0_full[i], e)

    J_c = np.zeros((4,4))
    for j in range(4):
        fp = eq.copy()
        fp[0,j] += eps
        new_p = np.zeros_like(fp)
        for i in range(N):
            left = fp[(i-1)%N]; right = fp[(i+1)%N]
            e = (left+right)/2; e = e/np.sum(e)
            new_p[i] = ds_combine(fp[i], e)
        J_c[:,j] = (new_p[0] - new_f0[0]) / eps

    J_c_proj = VTV_inv @ V.T @ J_c @ V
    ev = np.sort(np.abs(np.linalg.eigvals(J_c_proj)))[::-1]
    lam = ev[0]
    if lam > 0 and lam < 1:
        gap = -np.log(lam)
        ratio = gap / (-np.log(lambda_0))
        print(f"  {N:>4}  {lam:>14.10f}  {gap:>12.6f}  {-np.log(lambda_0):>12.6f}  {ratio:>8.4f}")
    else:
        print(f"  {N:>4}  {lam:>14.10f}  {'N/A':>12}  {-np.log(lambda_0):>12.6f}  {'N/A':>8}")
