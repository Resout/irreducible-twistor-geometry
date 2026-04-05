"""
Verify UV behaviour: compute the spacetime correlator by explicit
fiber integration of the DS construction.

No assumed Penrose kernels. We directly:
1. Set up CP³ → S⁴ twistor fibration
2. Define the DS observable O(x) = Born_i(m*|_{L_x}) on each fiber
3. Perturb the equilibrium, propagate via DS dynamics
4. Compute ⟨δO(x) δO(y)⟩ by integrating fluctuations over B ⊂ CP³
5. Measure how this scales with |x-y|

The question: does it go to a constant (as equation 38 suggests),
or does the fiber integration introduce a power law?
"""
import numpy as np
from scipy.linalg import expm

H = 3; FLOOR = 1.0/27.0

# ================================================================
# The DS equilibrium at K*=7/30
# ================================================================

def ds_combine(m, e):
    s, th = m[:3], m[3]; se, ph = e[:3], e[3]
    sn = s*se + s*ph + th*se; tn = th*ph
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i!=j)
    d = 1.0 - K
    out = np.zeros(4); out[:3] = sn/d; out[3] = tn/d
    born = out[3]**2/np.sum(out**2)
    if born < FLOOR:
        ss = np.sum(np.abs(out[:3])); ssq = np.sum(out[:3]**2)
        r = ssq/ss**2; t = (np.sqrt(26*r)-r)/(26-r)
        out[:3] = np.abs(out[:3])*(1-t)/ss; out[3] = t
    return out, K

def make_ev(p):
    pw = (1-p)/2; sc = 26.0/27.0
    raw = np.array([np.sqrt(p*sc), np.sqrt(pw*sc), np.sqrt(pw*sc), np.sqrt(1.0/27.0)])
    return raw / np.sum(raw)

def find_equilibrium():
    from scipy.optimize import brentq
    def K_of_p(p):
        e = make_ev(p)
        m = np.array([0.4, 0.2, 0.2, 0.2])
        for _ in range(2000):
            m2, K = ds_combine(m, e)
            if np.max(np.abs(m2-m)) < 1e-15: break
            m = m2
        _, K = ds_combine(m2, e)
        return K - 7.0/30.0, m2, e

    p = brentq(lambda p: K_of_p(p)[0], 0.92, 0.94, xtol=1e-15)
    _, m_star, e_star = K_of_p(p)
    return m_star, e_star

m_star, e_star = find_equilibrium()
_, K_check = ds_combine(m_star, e_star)
print(f"Equilibrium: K* = {K_check:.10f} (target {7/30:.10f})")
print(f"  m* = {m_star}")
print(f"  e* = {e_star}")

# Born probabilities at equilibrium
born_star = m_star**2 / np.sum(m_star**2)
print(f"  Born(m*) = {born_star}")

# ================================================================
# CP³ → S⁴ twistor fibration
#
# Represent CP³ via homogeneous coordinates Z = (Z⁰,Z¹,Z²,Z³) ∈ C⁴.
# The twistor fibration π: CP³ → S⁴ is defined by the Penrose
# incidence relation. For our purposes:
#
# A point x ∈ S⁴ ↪ R⁵ (unit 5-vector) determines a CP¹ in CP³.
# We use the quaternionic description:
#   x ↔ 2×2 matrix X = x⁰I + i(x¹σ₁+x²σ₂+x³σ₃)
#   The fiber L_x = {Z = (ω, Xω) : ω ∈ CP¹}
#
# For two points x, y ∈ R⁴ (via stereographic projection from S⁴):
#   X_x = x_μ σ_μ, X_y = y_μ σ_μ
#   d_FS(L_x, L_y) ≈ |x-y|/2 for |x-y| small
# ================================================================

sigma = [np.eye(2, dtype=complex),
         np.array([[0,1],[1,0]], dtype=complex),
         np.array([[0,-1j],[1j,0]], dtype=complex),
         np.array([[1,0],[0,-1]], dtype=complex)]

def spacetime_to_twistor_matrix(x):
    """Convert R⁴ point x to 2×2 matrix X = x_μ σ_μ."""
    return sum(x[mu] * sigma[mu] for mu in range(4))

def fiber_point(x, omega):
    """
    A point on the CP¹ fiber over x ∈ R⁴.
    omega ∈ CP¹ parameterised as (cos(θ), e^{iφ} sin(θ)).
    Returns Z ∈ CP³ as (ω, X·ω) ∈ C⁴.
    """
    X = spacetime_to_twistor_matrix(x)
    pi = X @ omega
    return np.concatenate([omega, pi])

def fs_distance_cp3(Z, W):
    """Fubini-Study distance between Z, W ∈ CP³."""
    inner = np.abs(np.vdot(Z, W))**2
    norm_Z = np.vdot(Z, Z).real
    norm_W = np.vdot(W, W).real
    cos2 = inner / (norm_Z * norm_W)
    cos2 = min(cos2, 1.0)  # numerical safety
    return np.arccos(np.sqrt(cos2))

# ================================================================
# The DS observable on a fiber
#
# At each spacetime point x, the observable is defined by restricting
# the equilibrium mass function to the fiber L_x.
#
# The mass function m* lives at a specific point in CP³. Its "value"
# at a fiber point Z ∈ L_x is determined by the DS propagator:
#   O(Z) = Born_1(T^n(m*, e*))
# where n = d_FS(m*, Z) / δ is the number of DS steps.
#
# The spacetime observable is the fiber average:
#   O(x) = ∮_{L_x} O(Z) dZ / Vol(L_x)
# ================================================================

# DS propagator: n steps from m_star with evidence e_star
Delta = -np.log(0.28291)  # spectral gap
delta_step = 0.01  # FS displacement per DS step (normalisation)

def ds_propagated_born(d_fs):
    """
    Born probability of hypothesis 1 after propagating d_fs in FS distance.
    Uses the spectral decomposition: Born_1(n steps) ≈ Born_1* + c₁ λ₀^n
    where n = d_fs / δ.
    """
    n_steps = d_fs / delta_step
    # At equilibrium: Born_1 = born_star[0]
    # Fluctuation decays as λ₀^n
    return born_star[0]  # At equilibrium, no fluctuation — need to perturb

def compute_fiber_observable(x, m_point, n_fiber_samples=100):
    """
    Compute ⟨O⟩ on the fiber L_x by sampling points on CP¹.

    For each point Z on the fiber, compute d_FS(m_point, Z)
    and propagate the DS dynamics that many steps.
    The observable is the Born probability after propagation.
    """
    # Sample CP¹: parameterise omega = (cos θ, e^{iφ} sin θ)
    thetas = np.linspace(0, np.pi, n_fiber_samples, endpoint=False)
    phis = np.linspace(0, 2*np.pi, n_fiber_samples, endpoint=False)

    values = []
    for th in thetas:
        omega = np.array([np.cos(th), np.exp(1j*phis[0])*np.sin(th)])
        Z = fiber_point(x, omega)
        d = fs_distance_cp3(m_point, Z)
        # DS propagator: after n = d/δ steps, fluctuation decays as λ₀^n
        n = d / delta_step
        # The connected correlator contribution from this point
        values.append(0.28291**n)  # leading eigenvalue decay

    return np.mean(values)

# ================================================================
# Direct computation: spacetime correlator vs separation
#
# Instead of the full fiber integral, use the KNOWN spectral
# decomposition to compute the exact correlator:
#
# S₂(x,y) = Σ_k |c_k|² λ_k^{d_FS(L_x, L_y)/δ}
#
# But the question is: does the PENROSE TRANSFORM of this quantity
# (the fiber integral) introduce additional structure?
#
# The fiber-averaged observable at x is:
#   ⟨O⟩_x = (1/Vol) ∮_{L_x} Born_1(m(Z)) dZ
#
# For a perturbation δm at the equilibrium point m* ∈ CP³:
#   δO(x) = (1/Vol) ∮_{L_x} [∂Born_1/∂m · G(m*, Z)] dZ
#
# where G(m*, Z) = Σ_k c_k v_k λ_k^{d(m*,Z)/δ} is the Green function.
#
# The two-point function:
#   ⟨δO(x) δO(y)⟩ = (1/Vol²) ∮_{L_x} ∮_{L_y} Σ_k |c_k|² λ_k^{(d(Z)+d(W))/δ} dZ dW
#
# For Z ∈ L_x and W ∈ L_y, the path from m* to Z and from m* to W
# are INDEPENDENT. So d(Z) + d(W) ≥ d(Z,W) by triangle inequality.
#
# KEY INSIGHT: The FS distance from m* to points on the fiber L_x
# depends on WHERE on the fiber. This variation is the mechanism
# that could introduce non-trivial scaling.
# ================================================================

print("\n" + "="*65)
print("FIBER DISTANCE VARIATION")
print("="*65)
print("\nFor each spacetime point x, the fiber L_x is a CP¹ in CP³.")
print("The distance from the equilibrium m* to points on L_x VARIES")
print("as you move around the fiber. This variation is key.\n")

# Place m* at the "north pole" of CP³
# m* as a point in CP³: proportional to (m*₁, m*₂, m*₃, m*₄) ∈ C⁴
m_point = m_star.astype(complex)  # treat as homogeneous coords in CP³

# Compute fiber distances for various spacetime separations
for r in [0.01, 0.05, 0.1, 0.5, 1.0, 2.0]:
    x = np.array([1.0, r, 0, 0])  # point at distance ~r from origin

    # Sample the fiber L_x
    n_samples = 200
    distances = []
    for i in range(n_samples):
        th = np.pi * i / n_samples
        phi = 0  # one slice suffices for statistics
        omega = np.array([np.cos(th), np.exp(1j*phi)*np.sin(th)])
        Z = fiber_point(x, omega)
        d = fs_distance_cp3(m_point, Z)
        distances.append(d)

    distances = np.array(distances)
    print(f"  r = {r:.4f}: d_FS range [{distances.min():.4f}, {distances.max():.4f}], "
          f"mean = {distances.mean():.4f}, std = {distances.std():.4f}")

# ================================================================
# THE CRITICAL TEST: correlator vs spacetime separation
# ================================================================
print("\n" + "="*65)
print("SPACETIME CORRELATOR: FIBER-INTEGRATED DS PROPAGATOR")
print("="*65)

def spacetime_correlator(r, n_fiber=200):
    """
    Compute ⟨δO(0) δO(r)⟩ by fiber integration.

    O(x) = (1/Vol) ∮_{L_x} Σ_k c_k λ_k^{d(m*,Z)/δ} dZ

    For the connected correlator with a single dominant mode:
    ⟨δO(0)δO(r)⟩ = |c₀|² × [∮_{L_0} λ₀^{d(m*,Z)/δ} dZ]
                           × [∮_{L_r} λ₀^{d(m*,W)/δ} dW]

    Each fiber integral independently samples the propagator.
    """
    lambda0 = 0.28291

    # Fiber over origin: x = (1, 0, 0, 0)
    x0 = np.array([1.0, 0, 0, 0])
    # Fiber over point at distance r: x = (1, r, 0, 0)
    xr = np.array([1.0, r, 0, 0])

    # Integrate over both fibers
    integral_0 = 0
    integral_r = 0

    for i in range(n_fiber):
        th = np.pi * i / n_fiber
        omega = np.array([np.cos(th), np.sin(th) + 0j])

        # Fiber over origin
        Z0 = fiber_point(x0, omega)
        d0 = fs_distance_cp3(m_point, Z0)
        integral_0 += lambda0**(d0/delta_step)

        # Fiber over r
        Zr = fiber_point(xr, omega)
        dr = fs_distance_cp3(m_point, Zr)
        integral_r += lambda0**(dr/delta_step)

    integral_0 /= n_fiber
    integral_r /= n_fiber

    # Connected correlator = product of fiber integrals
    # (this is the factorised form; the connected part subtracts ⟨O⟩²)
    return integral_0 * integral_r

print(f"\n{'r':>10s}  {'C(r)':>14s}  {'log₁₀(C)':>10s}  {'slope':>8s}")
print(f"  {'-'*50}")

separations = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
log_r = []
log_C = []

for r in separations:
    C = spacetime_correlator(r, n_fiber=500)
    lr = np.log10(r)
    lc = np.log10(abs(C)) if C > 0 else float('nan')
    log_r.append(lr)
    log_C.append(lc)

    if len(log_r) >= 2:
        slope = (log_C[-1] - log_C[-2]) / (log_r[-1] - log_r[-2])
    else:
        slope = float('nan')

    print(f"  {r:8.4f}  {C:14.6e}  {lc:10.4f}  {slope:8.3f}")

# Fit power law at short distances
if len(log_r) >= 3:
    coeffs_short = np.polyfit(log_r[:4], log_C[:4], 1)
    coeffs_long = np.polyfit(log_r[-3:], log_C[-3:], 1)
    print(f"\nShort-distance power law: C(r) ~ r^{coeffs_short[0]:.3f}")
    print(f"Long-distance power law:  C(r) ~ r^{coeffs_long[0]:.3f}")
    print(f"\nIf short-distance slope ≈ 0: fiber correlator is FLAT (no UV singularity)")
    print(f"If short-distance slope < 0: fiber integration introduces power law")
    print(f"If short-distance slope ≈ -8: matches tr(F²) conformal dimension")

print(f"\n" + "="*65)
print("INTERPRETATION")
print("="*65)
print("""
What we're measuring: the DS spacetime correlator obtained by
integrating the Born probability over the CP¹ fiber at each point.

Case A: If C(r) → constant as r → 0
  → The DS fiber correlator IS the spacetime correlator
  → No UV singularity from the fiber integration
  → The remark's claim about "Penrose kernel" is WRONG
  → The UV physics must come from somewhere else

Case B: If C(r) ~ r^{-p} as r → 0 for some p > 0
  → The fiber integration introduces a power-law singularity
  → p should match the conformal dimension of the observable
  → The remark's mechanism is CORRECT
""")
