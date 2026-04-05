"""
BIRKHOFF FACTORISATION AT THE K*=7/30 EQUILIBRIUM
===================================================

The single most valuable computation remaining in the framework.

Context: The Maurer-Cartan form M⁻¹dM is identically flat. The PHYSICAL
gauge connection is the Ward-reconstructed connection h₊⁻¹∂_μ h₊, obtained
by Birkhoff factorisation of the transition function on twistor lines.

At the equilibrium m*, the transition function on each twistor line
L_x ≅ CP¹ (parameterised by ζ) is:

    M_tilde(ζ, ζ̄) = M₀(ζ) · (I + M₀⁻¹ · ε(ζ, ζ̄))

where M₀ is the holomorphic part (constant on fibres at equilibrium)
and ε is the anti-holomorphic correction from the Born floor.

Because ||ε|| ≤ C(H) = 0.344 < 1, the Birkhoff factorisation reduces
to a convergent Neumann series. At leading order:

    h₊(ζ) = M₀ · (I + P₊[M₀⁻¹ ε])

where P₊ is the Cauchy projector (keeps the holomorphic part):

    P₊[f](ζ) = (1/2πi) ∮_{|z|=1} f(z)/(z - ζ) dz

The physical connection is A_μ = h₊⁻¹ ∂_μ h₊.
The physical curvature is F = dA + A∧A.

This script:
  1. Constructs M₀ and ε at the equilibrium
  2. Computes the Cauchy projector integral numerically on CP¹
  3. Extracts h₊ at leading order
  4. Computes the Ward-reconstructed connection A
  5. Computes the physical curvature F_Ward
  6. Reports: is F_Ward = 0 or ≠ 0 at the equilibrium?

If F_Ward ≠ 0: the vacuum has classical curvature (condensate).
If F_Ward = 0: the vacuum is classically flat (perturbative vacuum).

Requirements: mpmath
Runtime: ~1-5 minutes
"""

import sys
import time
import numpy as np

try:
    from mpmath import (mp, mpf, mpc, sqrt, log, matrix, nstr, fabs,
                        pi, exp, cos, sin, quad, eye, inf, j as J_imag)
except ImportError:
    print("ERROR: mpmath required"); sys.exit(1)

# ============================================================
# PRECISION FIRST
# ============================================================
PRECISION = 100  # 100 digits is plenty for this integral computation
mp.dps = PRECISION
ONE = mpf(1)
ZERO = mpf(0)
FLOOR = ONE / mpf(27)
TARGET_K = mpf(7) / mpf(30)
I2 = eye(2)

# Pauli matrices (mpmath matrix format)
sigma1 = matrix([[0, 1], [1, 0]])
sigma2 = matrix([[0, mpc(0, -1)], [mpc(0, 1), 0]])
sigma3 = matrix([[1, 0], [0, -1]])
sigmas = [sigma1, sigma2, sigma3]


def mass_to_M(m):
    """Mass function m = [s1,s2,s3,theta] -> 2x2 matrix M = (theta*I + s·sigma)/sqrt(2)."""
    s1, s2, s3, th = m
    sq2 = sqrt(mpf(2))
    return (th * I2 + s1 * sigma1 + s2 * sigma2 + s3 * sigma3) / sq2


def mat_inv_2x2(M):
    """Inverse of a 2x2 mpmath matrix."""
    a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
    det = a*d - b*c
    return matrix([[d/det, -b/det], [-c/det, a/det]])


def mat_mul(A, B):
    """2x2 matrix multiply."""
    return A * B


def mat_norm_sq(M):
    """Frobenius norm squared of 2x2 matrix."""
    return sum(fabs(M[i,k])**2 for i in range(2) for k in range(2))


def mat_norm(M):
    return sqrt(mat_norm_sq(M))


# ============================================================
# DS FRAMEWORK (same as before, precision-safe)
# ============================================================
def ds_combine(m, e):
    s1, s2, s3, th = m; e1, e2, e3, ph = e
    sn = [s1*e1+s1*ph+th*e1, s2*e2+s2*ph+th*e2, s3*e3+s3*ph+th*e3]
    tn = th*ph
    K = s1*e2+s1*e3+s2*e1+s2*e3+s3*e1+s3*e2
    d = ONE-K
    return [sn[0]/d, sn[1]/d, sn[2]/d, tn/d], K

def floor_enforce(m):
    s1, s2, s3, th = m
    ssq = s1**2+s2**2+s3**2
    born = th**2/(ssq+th**2)
    if born < FLOOR:
        ss = s1+s2+s3
        r = ssq/ss**2
        ac = mpf(26)-r; bc = mpf(2)*r; cc = -r
        t = (-bc+sqrt(bc**2-mpf(4)*ac*cc))/(mpf(2)*ac)
        sc = (ONE-t)/ss
        return [s1*sc, s2*sc, s3*sc, t]
    return list(m)

def full_step(m, e):
    m_ds, K = ds_combine(m, e)
    return floor_enforce(m_ds), K

def make_evidence(p_dom):
    pw = (ONE-p_dom)/mpf(2); sc = ONE-FLOOR
    raw = [sqrt(p_dom*sc), sqrt(pw*sc), sqrt(pw*sc), sqrt(FLOOR)]
    tot = sum(raw)
    return [r/tot for r in raw]

_last_m = [mpf('0.4'), mpf('0.15'), mpf('0.15'), mpf('0.3')]
def find_fp(e):
    global _last_m
    m = list(_last_m)
    tol = mpf(10)**(-(PRECISION-20))
    for i in range(100000):
        m2, _ = full_step(m, e)
        d = max(fabs(m2[k]-m[k]) for k in range(4))
        if d < tol:
            _last_m = m2; return m2, i
        m = m2
    _last_m = m2; return m2, 100000

def K_at_pdom(p):
    e = make_evidence(p)
    m, _ = find_fp(e)
    _, K = ds_combine(m, e)
    return K


# ============================================================
# FIND EQUILIBRIUM
# ============================================================
print("=" * 70)
print("BIRKHOFF FACTORISATION AT K*=7/30 EQUILIBRIUM")
print(f"Working precision: {PRECISION} digits")
print("=" * 70)

print("\n  Finding equilibrium (Secant method)...")
t0 = time.time()
p0, p1 = mpf('0.932'), mpf('0.933')
K0, K1 = K_at_pdom(p0), K_at_pdom(p1)
tol = mpf(10)**(-(PRECISION-30))
for n_sec in range(50):
    f0, f1 = K0-TARGET_K, K1-TARGET_K
    if fabs(f1-f0) < mpf(10)**(-PRECISION+5): break
    p2 = p1 - f1*(p1-p0)/(f1-f0)
    p0, K0 = p1, K1
    p1, K1 = p2, K_at_pdom(p2)
    if fabs(p1-p0) < tol: break

p_dom = p1
e_star = make_evidence(p_dom)
m_star, _ = find_fp(e_star)
_, K_star = ds_combine(m_star, e_star)
print(f"  Done in {time.time()-t0:.1f}s, |K-7/30| = {nstr(fabs(K_star-TARGET_K),5)}")

M_star = mass_to_M(m_star)
M_star_inv = mat_inv_2x2(M_star)
print(f"  m* = [{nstr(m_star[0],15)}, {nstr(m_star[1],15)}, {nstr(m_star[2],15)}, {nstr(m_star[3],15)}]")
print(f"  det(M*) = {nstr(M_star[0,0]*M_star[1,1]-M_star[0,1]*M_star[1,0], 15)}")


# ============================================================
# CONSTRUCT ε(ζ̄) ON THE TWISTOR LINE
# ============================================================
print("\n" + "=" * 70)
print("STEP 1: Construct the non-holomorphic perturbation ε(ζ̄)")
print("=" * 70)

# At the equilibrium, the transition function on a twistor line L_x is:
#   M_tilde(ζ, ζ̄) = M(Φ(m(ζ, ζ̄)))
# where m depends on the fibre coordinate ζ through the twistor incidence.
#
# At equilibrium with fibre uniformity: m(x) = m* for all x.
# The HOLOMORPHIC part: M₀ = M(m*) = constant on L_x.
# The ANTI-HOLOMORPHIC part: comes from the floor's ζ̄-dependence.
#
# On the twistor line, parameterise by ζ ∈ CP¹ = S².
# The mass function at fibre position ζ has an anti-holomorphic correction
# from the DS dynamics: m → Φ(m) where Φ includes the floor.
#
# The perturbation ε = M₀⁻¹ · (M(Φ(m)) - M₀) encodes the floor's
# anti-holomorphic contribution.
#
# At the equilibrium: Φ(m*) = m*, so naively ε = 0. BUT the key is that
# ε depends on the FIBRE DIRECTION, not on the base point. Different
# fibre directions (different ζ on L_x) see different projections of
# the 4D mass function onto the 2D twistor line.
#
# The standard twistor incidence: a point Z = (ω^A, π_A') ∈ CP³
# with ω^A = x^{AA'} π_A'. On the twistor line L_x (fixed x),
# π_A' = (1, ζ) parameterises the fibre. The mass function decomposes as:
#
#   m(ζ) = m_+ · (1, ζ)^T  (schematically: the restriction to the line)
#
# For the su(2) embedding with m = (s1,s2,s3,θ):
#   M(ζ) = (θI + s·σ)/√2
# evaluated at the mass function restricted to the twistor line through
# the incidence relation.
#
# At fibre uniformity (equilibrium), the mass function doesn't depend
# on ζ, so M₀(ζ) = M* = const. The anti-holomorphic correction comes
# from the DYNAMICS: applying Φ at each ζ produces ζ̄-dependent output
# because the floor enforcement acts differently at different fibre
# positions (the Born probability depends on which components are
# "visible" from that fibre direction).

# To compute ε(ζ̄), we need the anti-holomorphic Jacobian dbar(Φ)
# projected onto the fibre direction. On the twistor line, the fibre
# coordinate ζ corresponds to a specific direction in CP³.

# The anti-holomorphic part of the floor correction at m* is:
# ε_ij = M⁻¹ · (∂M/∂m_k) · (∂̄Φ)_k  evaluated along the fibre

# We compute ∂̄Φ at m* using Wirtinger finite differences
print("  Computing anti-holomorphic Jacobian ∂̄Φ at m*...")

eps_w = mpf(10)**(-30)  # Wirtinger step
J_anti = matrix(4, 4)   # (∂̄Φ)_{j,alpha} = ∂Φ_j/∂m̄_alpha

for alpha in range(4):
    e_alpha_r = [ZERO]*4; e_alpha_r[alpha] = ONE
    e_alpha_i = [ZERO]*4; e_alpha_i[alpha] = ONE

    # Real perturbation
    m_pr = [m_star[k] + eps_w * e_alpha_r[k] for k in range(4)]
    m_mr = [m_star[k] - eps_w * e_alpha_r[k] for k in range(4)]
    phi_pr, _ = full_step(m_pr, e_star)
    phi_mr, _ = full_step(m_mr, e_star)
    df_dx = [(phi_pr[k]-phi_mr[k])/(2*eps_w) for k in range(4)]

    # Imaginary perturbation (complexify the mass function)
    # For real m* perturbed by i*eps in direction alpha:
    # The DS map is polynomial (holomorphic in m), and the floor involves |theta|.
    # For small imaginary perturbation, use the Wirtinger relation:
    # At real m: ∂̄Φ/∂m̄_α = (1/2)(∂Φ/∂x_α + i·∂Φ/∂y_α)
    # For real-valued Φ at real inputs: ∂Φ/∂y_α comes from the floor's
    # dependence on |theta| = sqrt(theta·theta_bar).
    #
    # For the Wirtinger derivative at real inputs:
    # ∂̄Φ = (1/2) dΦ/dx (since ∂Φ/∂y = i·(∂Φ_hol/∂x - ∂Φ/∂x) and
    # for the floor: the imaginary response is determined by the floor's
    # |·| dependence).
    #
    # Actually, for a REAL map f: R^n -> R^n, the Wirtinger anti-holomorphic
    # derivative at real inputs is: ∂̄f/∂z̄_α = (1/2) ∂f/∂x_α
    # This is because ∂f/∂y_α = 0 for a map that takes reals to reals
    # when evaluated at real inputs... NO, that's wrong. The floor involves
    # |theta| which IS differentiable with respect to Im(theta).
    #
    # For the correct computation, we need to extend Φ to complex inputs.
    # Since we're at 100 digits, let's use the complex floor enforcement.

    # Complex perturbation in imaginary direction
    # m* + i·eps·e_alpha: complexify
    m_pi = [mpc(m_star[k], eps_w * e_alpha_i[k]) for k in range(4)]
    m_mi = [mpc(m_star[k], -eps_w * e_alpha_i[k]) for k in range(4)]

    # Need complex-capable floor enforcement
    phi_pi = complex_full_step(m_pi, e_star)
    phi_mi = complex_full_step(m_mi, e_star)
    df_dy = [(phi_pi[k]-phi_mi[k])/(2*eps_w) for k in range(4)]

    for k in range(4):
        J_anti[k, alpha] = (df_dx[k] + mpc(0, 1)*df_dy[k]) / 2

print(f"  ||∂̄Φ|| = {nstr(sqrt(sum(fabs(J_anti[i,k])**2 for i in range(4) for k in range(4))), 10)}")

# Check rank (should be ~1 from the rank-1 finding)
# Compute SVD via eigenvalues of J_anti^H · J_anti
JHJ = J_anti.T.conjugate() * J_anti
eigs_JHJ = mp.eigvals(JHJ) if hasattr(mp, 'eigvals') else None
if eigs_JHJ:
    svals = sorted([sqrt(fabs(e)) for e in eigs_JHJ], reverse=True)
    print(f"  Singular values of ∂̄Φ: {[nstr(s,6) for s in svals]}")
    print(f"  Effective rank: {sum(1 for s in svals if s > mpf(10)**(-10))}")


# ============================================================
# Complex floor enforcement
# ============================================================
def complex_floor_enforce(m):
    """Floor enforcement for complex mass functions."""
    s1, s2, s3, th = m
    s_mag_sq = fabs(s1)**2 + fabs(s2)**2 + fabs(s3)**2
    th_mag = fabs(th)
    total_mag_sq = s_mag_sq + th_mag**2
    if total_mag_sq < mpf(10)**(-50):
        return list(m)
    born = th_mag**2 / total_mag_sq
    if born >= FLOOR - mpf(10)**(-20):
        return list(m)
    # Floor active
    S = s1 + s2 + s3
    S_mag_sq = fabs(S)**2
    if S_mag_sq < mpf(10)**(-50):
        return list(m)
    phase = th / th_mag if th_mag > mpf(10)**(-50) else ONE
    R = s_mag_sq / S_mag_sq
    Re_phase = mpf(phase.real) if hasattr(phase, 'real') else phase
    ac = mpf(26) - R
    bc = mpf(2) * R * Re_phase
    cc = -R
    disc = bc**2 - 4*ac*cc
    r_new = (-bc + sqrt(disc)) / (2*ac)
    if r_new <= 0:
        r_new = (-bc - sqrt(disc)) / (2*ac)
    th_new = phase * r_new
    alpha = (ONE - th_new) / S
    return [s1*alpha, s2*alpha, s3*alpha, th_new]

def complex_ds_combine(m, e):
    s1, s2, s3, th = m; e1, e2, e3, ph = e
    sn = [s1*e1+s1*ph+th*e1, s2*e2+s2*ph+th*e2, s3*e3+s3*ph+th*e3]
    tn = th*ph
    K = s1*e2+s1*e3+s2*e1+s2*e3+s3*e1+s3*e2
    d = ONE-K
    return [sn[0]/d, sn[1]/d, sn[2]/d, tn/d], K

def complex_full_step(m, e):
    m_ds, K = complex_ds_combine(m, e)
    return complex_floor_enforce(m_ds)


# ============================================================
# STEP 2: Project ε onto the twistor line
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Construct ε on the twistor line CP¹")
print("=" * 70)

# On the twistor line L_x parameterised by ζ = e^{iφ} ∈ S¹ ⊂ CP¹,
# the fibre direction corresponds to a specific tangent vector in CP³.
#
# The standard twistor fibration uses homogeneous coords [Z₀:Z₁:Z₂:Z₃].
# Identifying (s₁,s₂,s₃,θ) ↔ (Z₀,Z₁,Z₂,Z₃), the twistor line through
# the origin of S⁴ is parameterised by [s₁:s₂] free, [s₃:θ] = [ζ:1].
#
# At the equilibrium, the perturbation in the fibre direction (varying ζ)
# generates an anti-holomorphic response through ∂̄Φ.
#
# The fibre tangent vector at ζ = e^{iφ} is:
#   v_fibre(φ) = (0, 0, cos(φ), sin(φ)) (real part)
#              + (0, 0, -sin(φ), cos(φ)) (imaginary part)
# This is the direction in which ζ varies on the twistor line.

# The perturbation ε on the twistor line is:
#   ε(ζ̄) = M*⁻¹ · dM · ∂̄Φ · v_fibre(ζ̄)
# where v_fibre depends on ζ̄ (the anti-holomorphic fibre coordinate).

# For the Cauchy projector, we need ε as a function of ζ on the unit circle.
# Parameterise ζ = e^{iφ}, so ζ̄ = e^{-iφ}.

# dM/dm_j are constant: σ_j/√2 for j=0,1,2 and I/√2 for j=3
sq2 = sqrt(mpf(2))
dM = [sigma1/sq2, sigma2/sq2, sigma3/sq2, I2/sq2]

def epsilon_on_circle(phi_val):
    """Compute ε(ζ̄) at ζ = e^{iφ} on the twistor line.

    The fibre direction at angle φ is v = (0, 0, e^{-iφ}, i·e^{-iφ})
    (the ∂/∂ζ̄ direction in the (s₃, θ) plane).

    ε = M*⁻¹ · Σ_j dM_j · Σ_α (∂̄Φ)_{j,α} · v_α
    """
    zeta_bar = exp(mpc(0, -phi_val))

    # The fibre direction ∂/∂ζ̄ in the (s₃,θ) coordinates:
    # ζ parameterises [s₃:θ] = [ζ:1], so varying ζ̄ gives
    # v_2 = ζ̄ (for s₃ component, index 2) and v_3 = 0 for θ
    # Actually, for the anti-holomorphic variation:
    # v = (0, 0, 1, 0) in the ζ̄ direction (only s₃ responds)
    # But this depends on the specific incidence relation.
    #
    # More carefully: on the twistor line, ζ̄ enters through the
    # anti-holomorphic part of the incidence relation. The perturbation
    # of the mass function is:
    #   δm_α = (∂m_α/∂ζ̄) · δζ̄
    # At the equilibrium with fibre uniformity, ∂m/∂ζ̄ = 0 (m is constant).
    # The anti-holomorphic content comes from the DYNAMICS (∂̄Φ), not from
    # the incidence relation directly.
    #
    # The correct construction: ε(ζ̄) captures how the floor enforcement
    # breaks the holomorphic structure. On the twistor line, the transition
    # function M(ζ) = M* (constant). The floor makes M depend on ζ̄ through:
    #   M_floor(ζ, ζ̄) = M(Φ(m*(ζ, ζ̄)))
    # At the equilibrium: m* doesn't depend on ζ, but the FLOOR RESPONSE
    # does depend on the fibre direction because the floor enforcement
    # treats different mass components differently.
    #
    # The perturbation is:
    #   ε(ζ̄) = M*⁻¹ · (M_floor - M*) ≈ M*⁻¹ · dM · ∂̄Φ · v(ζ̄)
    # where v(ζ̄) is the anti-holomorphic fibre tangent.
    #
    # On the standard twistor line with π_A' = (1, ζ):
    # The anti-holomorphic tangent is ∂/∂ζ̄ which maps to the direction
    # (0, 0, 1, 0) rotated by ζ̄ in the (s₃, θ) plane.
    # For a general point on S⁴, the fibre tangent rotates the mass
    # components in a ζ-dependent way.

    # For the SIMPLEST case (origin of S⁴):
    # v(ζ̄) = (0, 0, ζ̄, 0) in the (s₁, s₂, s₃, θ) basis
    v = [ZERO, ZERO, zeta_bar, ZERO]

    # ε = M*⁻¹ · Σ_j (dM_j · Σ_α J_anti[j,α] · v[α])
    inner = [sum(J_anti[kk, alpha] * v[alpha] for alpha in range(4)) for kk in range(4)]
    dM_contribution = sum(dM[kk] * inner[kk] for kk in range(4))
    return mat_mul(M_star_inv, dM_contribution)


# Sample ε on the unit circle
N_sample = 64
print(f"\n  Sampling ε(ζ̄) at {N_sample} points on the unit circle...")

eps_samples = []
eps_norms = []
for i in range(N_sample):
    phi = 2 * pi * mpf(i) / mpf(N_sample)
    eps_val = epsilon_on_circle(phi)
    eps_samples.append(eps_val)
    eps_norms.append(mat_norm(eps_val))

max_eps = max(eps_norms)
mean_eps = sum(eps_norms) / len(eps_norms)
print(f"  max ||ε|| over circle = {nstr(max_eps, 10)}")
print(f"  mean ||ε|| = {nstr(mean_eps, 10)}")
print(f"  Neumann series converges: {max_eps < ONE}")


# ============================================================
# STEP 3: Cauchy projector — extract h₊
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Cauchy projector P₊[M₀⁻¹ ε]")
print("=" * 70)

# P₊[f](ζ) = (1/2πi) ∮ f(z)/(z-ζ) dz
#
# For ζ INSIDE the unit circle, the Cauchy integral picks up the
# holomorphic (positive Fourier) part of f.
#
# We compute P₊[ε](ζ) at a specific interior point ζ₀.
# Then h₊(ζ₀) = M₀ · (I + P₊[ε](ζ₀)) at leading order.

# Compute the Fourier decomposition of ε(φ)
# ε(φ) = Σ_n ε_n e^{inφ}
# P₊ keeps n ≥ 0 terms

print(f"  Computing Fourier modes of ε...")
n_modes = N_sample // 2
fourier_modes = []  # List of 2x2 matrices, one per Fourier mode

for n in range(-n_modes, n_modes + 1):
    # ε_n = (1/2π) ∫₀²π ε(φ) e^{-inφ} dφ
    # Numerical integration via trapezoidal rule (exponentially convergent for smooth periodic)
    coeff = matrix(2, 2)
    for i in range(N_sample):
        phi = 2 * pi * mpf(i) / mpf(N_sample)
        phase = exp(mpc(0, -n * float(phi)))
        for r in range(2):
            for c in range(2):
                coeff[r, c] += eps_samples[i][r, c] * phase
    coeff = coeff / mpf(N_sample)
    fourier_modes.append((n, coeff))

# Print mode norms
print(f"\n  Fourier mode norms:")
for n, coeff in fourier_modes:
    nm = mat_norm(coeff)
    if nm > mpf(10)**(-20):
        print(f"    n={n:+3d}: ||ε_n|| = {nstr(nm, 8)}")

# P₊ keeps modes with n ≥ 0
# P₊[ε] = Σ_{n≥0} ε_n e^{inφ}
# At ζ = 0 (centre of disc): P₊[ε](0) = ε₀ (just the zero mode)

print(f"\n  P₊[ε](ζ=0) = ε₀ (zero Fourier mode):")
eps_0 = [c for n, c in fourier_modes if n == 0][0]
for r in range(2):
    print(f"    [{nstr(eps_0[r,0], 12)}, {nstr(eps_0[r,1], 12)}]")
print(f"  ||ε₀|| = {nstr(mat_norm(eps_0), 10)}")


# ============================================================
# STEP 4: Compute h₊ and the Ward-reconstructed connection
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Ward-reconstructed connection")
print("=" * 70)

# At leading order: h₊(ζ) = M₀ · (I + P₊[M₀⁻¹ ε](ζ))
# Since M₀⁻¹ ε is what we Fourier-analysed, P₊[M₀⁻¹ ε] is the positive part.

# At ζ = 0: h₊ = M* · (I + ε₀)
h_plus_0 = mat_mul(M_star, I2 + eps_0)
print(f"  h₊(0) = M* · (I + ε₀):")
for r in range(2):
    print(f"    [{nstr(h_plus_0[r,0], 12)}, {nstr(h_plus_0[r,1], 12)}]")

# The physical connection involves ∂_μ h₊, which requires the spacetime
# dependence. At the equilibrium with fibre uniformity, the spacetime
# derivative ∂_μ h₊ comes from how h₊ varies with the base point x.
#
# Since m* is the same at all x, the variation of h₊ with x comes from
# the variation of the twistor line L_x with x. Different x gives
# different twistor lines, hence different Cauchy integrals.
#
# The spacetime derivative of h₊ at x₀ is:
#   ∂_μ h₊ = M* · ∂_μ P₊[ε](ζ₀)
# where ∂_μ acts by changing which twistor line we're on.
#
# For the Penrose transform, the spacetime derivative maps to
# the fibre derivative via the incidence relation:
#   ∂/∂x^{AA'} ↔ π_{A'} ∂/∂ω^A
#
# This means ∂_μ h₊ involves the FIBRE DERIVATIVE of the positive
# Fourier modes, weighted by the spinor π_{A'} = (1, ζ).

# The physical curvature F is:
#   F_{μν} = ∂_μ A_ν - ∂_ν A_μ + [A_μ, A_ν]
# where A_μ = h₊⁻¹ ∂_μ h₊.
#
# At leading order in ε:
#   A_μ ≈ M*⁻¹ ∂_μ(M* ε₊) = ∂_μ ε₊ + [M*⁻¹ ∂_μ M*, ε₊]
#       = ∂_μ ε₊  (since M* is constant in x)
#
# So A ≈ ∂_μ ε₊ and F ≈ ∂_μ ∂_ν ε₊ - ∂_ν ∂_μ ε₊ + [∂_μ ε₊, ∂_ν ε₊]
#       = [∂_μ ε₊, ∂_ν ε₊] (since partial derivatives commute)
#
# This is NONZERO if the spacetime derivatives of ε₊ don't commute!

print(f"\n  Leading-order analysis:")
print(f"  At O(ε): A_μ ≈ ∂_μ(P₊[ε]) and F ≈ [∂_μ ε₊, ∂_ν ε₊]")
print(f"  The curvature is the COMMUTATOR of spacetime derivatives of the")
print(f"  holomorphic projection of the floor correction.")
print(f"  This is nonzero whenever ε₊ has non-commuting spacetime derivatives.")

# Check: do the Fourier modes give non-commuting contributions?
# The n=0 mode is constant → no spacetime derivative → no contribution to F.
# The n≥1 modes vary along the fibre → via Penrose transform → spacetime variation.
# F ≠ 0 requires at least TWO non-commuting n≥1 modes.

positive_modes = [(n, c) for n, c in fourier_modes if n > 0 and mat_norm(c) > mpf(10)**(-15)]
print(f"\n  Positive Fourier modes (n>0):")
for n, c in positive_modes:
    print(f"    n={n}: ||ε_n|| = {nstr(mat_norm(c), 8)}")

if len(positive_modes) >= 2:
    n1, c1 = positive_modes[0]
    n2, c2 = positive_modes[1]
    comm = mat_mul(c1, c2) - mat_mul(c2, c1)  # [ε_{n1}, ε_{n2}]
    print(f"\n  [ε_{{{n1}}}, ε_{{{n2}}}]:")
    for r in range(2):
        print(f"    [{nstr(comm[r,0], 10)}, {nstr(comm[r,1], 10)}]")
    print(f"  ||[ε_{{{n1}}}, ε_{{{n2}}}]|| = {nstr(mat_norm(comm), 10)}")
    if mat_norm(comm) > mpf(10)**(-10):
        print(f"\n  *** F_Ward ≠ 0: The Ward-reconstructed curvature is NONZERO ***")
        print(f"  The vacuum has classical gauge curvature.")
    else:
        print(f"\n  [ε_{{{n1}}}, ε_{{{n2}}}] ≈ 0: modes commute.")
        print(f"  Need higher-order analysis or more modes.")
elif len(positive_modes) == 1:
    print(f"\n  Only one positive mode → F = 0 at leading order (abelian).")
    print(f"  But the rank-1 ∂̄Φ predicts exactly this.")
    print(f"  F_Ward at this order is ZERO. Need O(ε²) or multi-site coupling.")
else:
    print(f"\n  No positive modes → ε is purely anti-holomorphic → P₊[ε] = ε₀ = const.")
    print(f"  F_Ward = 0 at leading order.")

print(f"\n  Total ||ε₊|| (positive part) = {nstr(sqrt(sum(mat_norm_sq(c) for n,c in positive_modes)), 10)}")
print(f"  Total ||ε₋|| (negative part) = {nstr(sqrt(sum(mat_norm_sq(c) for n,c in fourier_modes if n < 0 and mat_norm(c) > mpf(10)**(-15))), 10)}")
print(f"  ||ε₀|| (zero mode) = {nstr(mat_norm(eps_0), 10)}")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
  At the K*=7/30 equilibrium:
  - M₀ = M* (constant transition function)
  - ε encodes the Born floor's anti-holomorphic correction
  - ||ε|| ≤ {nstr(max_eps, 6)} on the twistor line
  - Cauchy projector decomposes ε into holomorphic (P₊) and anti-holomorphic (P₋) parts
  - The Ward-reconstructed connection A = h₊⁻¹ dh₊ has curvature
    determined by the commutator structure of the positive Fourier modes

  The physical curvature F_Ward is:
  - NONZERO if the positive Fourier modes of ε don't commute (non-abelian vacuum)
  - ZERO if they commute or if only one mode exists (abelian/flat vacuum)

  This resolves Gap 7: whether F_Ward = 0 or ≠ 0 at equilibrium.
""")

t_total = time.time() - t0
print(f"Total runtime: {t_total:.1f}s")
print("DONE.")
