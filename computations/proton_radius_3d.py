#!/usr/bin/env python3
"""
PROTON CHARGE RADIUS FROM 3D LATTICE OF DS SITES
=================================================

On a 3D cubic lattice, each site i carries mass function m_i = (s1, s2, s3, theta).
The evidence at each site is the average of its 6 nearest neighbours:
  e_i = (1/6) sum_{j nn i} m_j

Linearized dynamics around equilibrium:
  delta_m_i(t+1) = J_m . delta_m_i(t) + J_e . (1/6) sum_{j nn} delta_m_j(t)

In Fourier space:
  delta_m(k, t+1) = M(k) . delta_m(k, t)
  M(k) = J_m + (J_e/3)(cos kx + cos ky + cos kz)

SPATIAL CORRELATOR:
The EQUAL-TIME spatial correlator in the steady state is:
  C(r) = <delta_m(0,t) delta_m(r,t)^T>
This satisfies C = M C M^T + noise (discrete Lyapunov equation).

The charge-charge correlator is G(r) = c^T C(r) c where c = (1,-1,-1,0).

In Fourier space: G_tilde(k) = c^T S(k) c where S(k) solves the k-space Lyapunov eq.

For the SPATIAL POLE MASS (the physical proton mass from the spatial decay):
  G(r) ~ exp(-m_p * |r|) for large |r|
  m_p is determined by the analytic continuation of M(k) to imaginary k.

The proton mass in lattice units is the solution to:
  det(I - M(iq,0,0)^2) = 0  for the charge channel
i.e., the point where the largest eigenvalue of M reaches 1.
  lambda_charge(gamma) = 1 where gamma = (cosh(q) + 2)/3
  m_p * a = q

The form factor is obtained from the q-dependence of the charge propagator.
"""

import numpy as np
from mpmath import mp, mpf, matrix, sqrt, fabs, nstr, eig, chop, findroot, log, ln
from scipy.optimize import brentq

mp.dps = 50
H = 3
FLOOR = mpf(1) / mpf(H**3)
K_STAR = mpf(7) / 30

# Physical constants
hbar_c = 197.3269804   # MeV fm
m_proton = 938.272046   # MeV
lambda_C = hbar_c / m_proton
r_p_exp = 0.84087       # fm (muonic hydrogen)

print("=" * 80)
print("PROTON CHARGE RADIUS FROM 3D DS LATTICE")
print("=" * 80)
print(f"\nTarget: r_p(muonic) = {r_p_exp} fm")
print(f"        r_p = (H+1)*lambda_C = {(H+1)*lambda_C:.6f} fm")


# ============================================================
# DS DYNAMICS
# ============================================================

def ds_combine(m, e):
    s = m[:3]; theta = m[3]
    ev = e[:3]; phi = e[3]
    s_pre = [s[i]*ev[i] + s[i]*phi + theta*ev[i] for i in range(3)]
    theta_pre = theta * phi
    total_pre = sum(s_pre) + theta_pre
    K = mpf(1) - total_pre
    if fabs(mpf(1) - K) < mpf(10)**(-70):
        return list(m), K
    denom = mpf(1) - K
    return [sp / denom for sp in s_pre] + [theta_pre / denom], K

def born_prob(m):
    L2sq = sum(x**2 for x in m)
    if L2sq < mpf(10)**(-60): return mpf(0)
    return m[3]**2 / L2sq

def enforce_floor(m, floor_val=None):
    if floor_val is None: floor_val = FLOOR
    b = born_prob(m)
    if b >= floor_val - mpf(10)**(-60): return list(m)
    S = sum(m[:3])
    if S < mpf(10)**(-60): return list(m)
    Sq = sum(x**2 for x in m[:3])
    A_c = mpf(26) * S**2 - Sq
    B_c = mpf(2) * Sq
    C_c = -Sq
    disc = B_c**2 - 4*A_c*C_c
    if disc < 0: return list(m)
    t1 = (-B_c + sqrt(disc)) / (2*A_c)
    t2 = (-B_c - sqrt(disc)) / (2*A_c)
    cands = [t for t in [t1, t2] if mpf(0) < t < mpf(1)]
    if not cands: return list(m)
    t = min(cands, key=lambda x: fabs(x - m[3]))
    alpha = (mpf(1) - t) / S
    return [m[i]*alpha for i in range(3)] + [t]

def ds_step(m, e):
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds), K


# ============================================================
# STEP 1: FIND EQUILIBRIUM
# ============================================================
print("\n" + "=" * 80)
print("STEP 1: EQUILIBRIUM AND JACOBIANS")
print("=" * 80)

def find_equilibrium():
    from scipy.optimize import fsolve
    def eq_rough(p):
        s1, theta, w1, phi = p
        s2 = (1.0 - s1 - theta) / 2.0
        w2 = (1.0 - w1 - phi) / 2.0
        m = [s1, s2, s2, theta]; e = [w1, w2, w2, phi]
        L2m = s1**2 + 2*s2**2 + theta**2
        m_mp = [mpf(x) for x in m]; e_mp = [mpf(x) for x in e]
        m_out, _ = ds_step(m_mp, e_mp)
        return [theta**2/L2m - 1/27, phi**2/(w1**2+2*w2**2+phi**2) - 1/27,
                s1*w2+s1*w2+s2*w1+s2*w2+s2*w1+s2*w2 - 7/30,
                float(m_out[0]) - s1]
    sol = fsolve(eq_rough, [0.787, 0.155, 0.631, 0.129])

    def eq_mp(*p):
        s1, theta, w1, phi = p
        s2 = (mpf(1)-s1-theta)/2; w2 = (mpf(1)-w1-phi)/2
        m = [s1,s2,s2,theta]; e = [w1,w2,w2,phi]
        eq1 = theta**2/(s1**2+2*s2**2+theta**2) - FLOOR
        eq2 = phi**2/(w1**2+2*w2**2+phi**2) - FLOOR
        K = s1*w2+s1*w2+s2*w1+s2*w2+s2*w1+s2*w2
        m_out, _ = ds_step(m, e)
        return [eq1, eq2, K - mpf(7)/30, m_out[0] - s1]

    x = findroot(eq_mp, [mpf(str(v)) for v in sol])
    s1, theta, w1, phi = x
    s2 = (mpf(1)-s1-theta)/2; w2 = (mpf(1)-w1-phi)/2
    return [s1,s2,s2,theta], [w1,w2,w2,phi]

m_star, e_star = find_equilibrium()
print(f"  m* = [{', '.join(nstr(x, 15) for x in m_star)}]")
print(f"  e* = [{', '.join(nstr(x, 15) for x in e_star)}]")

# Jacobians
eps = mpf(10)**(-35)

J_m = matrix(4, 4)
for j in range(4):
    mp_p = list(m_star); mp_m = list(m_star)
    mp_p[j] += eps; mp_m[j] -= eps
    fp, _ = ds_step(mp_p, e_star)
    fm, _ = ds_step(mp_m, e_star)
    for i in range(4):
        J_m[i,j] = (fp[i] - fm[i]) / (2*eps)

J_e = matrix(4, 4)
for j in range(4):
    ep_p = list(e_star); ep_m = list(e_star)
    ep_p[j] += eps; ep_m[j] -= eps
    fp, _ = ds_step(m_star, ep_p)
    fm, _ = ds_step(m_star, ep_m)
    for i in range(4):
        J_e[i,j] = (fp[i] - fm[i]) / (2*eps)

# Convert to numpy
J_m_np = np.array([[float(J_m[i,j]) for j in range(4)] for i in range(4)])
J_e_np = np.array([[float(J_e[i,j]) for j in range(4)] for i in range(4)])

evals_m_np = np.linalg.eigvals(J_m_np)
idx = np.argsort(-np.abs(evals_m_np))
print(f"\n  J_m eigenvalues:")
for k in range(4):
    print(f"    lambda_{k}^(m) = {evals_m_np[idx[k]]:.12f}")

evals_e_np = np.linalg.eigvals(J_e_np)
idx_e = np.argsort(-np.abs(evals_e_np))
print(f"  J_e eigenvalues:")
for k in range(4):
    print(f"    lambda_{k}^(e) = {evals_e_np[idx_e[k]]:.12f}")


# ============================================================
# STEP 2: 3D TRANSFER MATRIX
# ============================================================
print("\n" + "=" * 80)
print("STEP 2: 3D TRANSFER MATRIX M(k)")
print("=" * 80)

def Mk_gamma(gamma, J_m_np, J_e_np):
    """Transfer matrix for structure factor gamma = (cos kx + cos ky + cos kz)/3."""
    return J_m_np + J_e_np * gamma

c_vec = np.array([1.0, -1.0, -1.0, 0.0])

# Eigenvalues at k=0 (gamma=1) and k=(pi,pi,pi) (gamma=-1)
for label, gam in [("k=0 (gamma=1)", 1.0), ("k=(pi,pi,pi) (gamma=-1)", -1.0)]:
    M = Mk_gamma(gam, J_m_np, J_e_np)
    evals, evecs = np.linalg.eig(M)
    idx = np.argsort(-np.abs(evals))
    print(f"\n  {label}:")
    for k in range(4):
        v = evecs[:, idx[k]]; v = v / np.linalg.norm(v)
        cv = np.dot(c_vec, v)
        print(f"    lambda_{k} = {evals[idx[k]]:.10f}  |c.v|^2 = {abs(cv)**2:.6f}")


# ============================================================
# STEP 3: SPATIAL MASS FROM POLE OF PROPAGATOR
# ============================================================
print("\n" + "=" * 80)
print("STEP 3: SPATIAL MASS FROM POLE")
print("=" * 80)

print("""
  The spatial propagator has a pole when the largest eigenvalue of M
  (in the charge channel) reaches 1.

  For k = (iq, 0, 0): gamma = (cosh(q) + 2)/3
  Find q_0 such that lambda_charge(gamma(q_0)) = 1.
  Then m_p * a = q_0.
""")

def charge_eigenvalue(gamma):
    """Return the eigenvalue with largest charge projection."""
    M = Mk_gamma(gamma, J_m_np, J_e_np)
    evals, evecs = np.linalg.eig(M)
    best_idx = -1
    best_cv2 = -1
    for alpha in range(4):
        v = evecs[:, alpha]; v = v / np.linalg.norm(v)
        cv2 = abs(np.dot(c_vec, v))**2
        if cv2 > best_cv2:
            best_cv2 = cv2
            best_idx = alpha
    return float(np.real(evals[best_idx]))

def dominant_eigenvalue(gamma):
    """Return the largest eigenvalue (by magnitude)."""
    M = Mk_gamma(gamma, J_m_np, J_e_np)
    evals = np.linalg.eigvals(M)
    return float(np.max(np.abs(evals)))

# Scan to find where charge eigenvalue crosses 1
print(f"  Scanning charge eigenvalue vs gamma:")
gammas = np.arange(1.0, 5.0, 0.25)
for g in gammas:
    lam_c = charge_eigenvalue(g)
    q = np.arccosh(3*g - 2) if 3*g - 2 >= 1 else 0
    print(f"    gamma={g:.2f}: lambda_charge = {lam_c:.8f}, q = {q:.4f}")

# Find exact gamma where lambda_charge = 1
def charge_eig_minus_1(gamma):
    return charge_eigenvalue(gamma) - 1.0

gamma_pole = brentq(charge_eig_minus_1, 1.0, 5.0)
q_pole = np.arccosh(3*gamma_pole - 2)
m_pa = q_pole  # proton mass in lattice units

print(f"\n  Pole found at gamma = {gamma_pole:.10f}")
print(f"  cosh(q) = {3*gamma_pole - 2:.10f}")
print(f"  q_0 = m_p * a = {q_pole:.10f}")

# Lattice spacing
a_fm = q_pole * hbar_c / m_proton
print(f"\n  Lattice spacing: a = q_0 * hbar_c / m_p = {a_fm:.6f} fm")
print(f"  (= {q_pole:.6f} * {hbar_c:.4f} / {m_proton:.3f})")


# ============================================================
# STEP 4: FORM FACTOR FROM CHARGE PROPAGATOR
# ============================================================
print("\n" + "=" * 80)
print("STEP 4: CHARGE FORM FACTOR")
print("=" * 80)

print("""
  The charge propagator in Fourier space is:
    G_tilde(k) = sum_alpha |c.v_alpha(k)|^2 * P_alpha(k)

  where P_alpha is the propagator for mode alpha.

  For a SPATIAL propagator (static correlator), the propagator of mode alpha
  at momentum k is:
    P_alpha(k) = 1 / (1 - lambda_alpha(k))  [one-sided geometric sum]
  or equivalently for the full (time-summed) correlator:
    P_alpha(k) = 1 / (1 - lambda_alpha(k)^2) = 1/[(1-lam)(1+lam)]

  Actually, the equal-time correlator from the Lyapunov equation is more
  subtle. Let me use the RESOLVENT approach.

  The spatial propagator for a particle of mass M on a cubic lattice:
    G(k) = 1 / (M^2 + 4 sum_mu sin^2(k_mu/2))   [continuum: 1/(M^2+k^2)]

  In our system, the "mass" for each mode is determined by the eigenvalue
  at k=0:
    M_alpha^2 = 4(1 - lambda_alpha(0)) / |d lambda/d(gamma)|

  But this mixes lattice and transfer-matrix concepts. Let me instead
  directly compute the spatial propagator.

  For a transfer matrix T acting in time, the equal-time correlator
  in the steady state satisfies the Lyapunov equation:
    C = T C T^T + Sigma
  where Sigma is the noise covariance.

  In Fourier space (for k-dependent T(k)):
    C(k) = T(k) C(k) T(k)^T + Sigma(k)

  This gives: C(k) = (I - T(k) kron T(k))^{-1} Sigma

  For our purposes, we can take a simpler route: the SPECTRAL REPRESENTATION.

  Each eigenmode alpha of T(k) contributes to the correlator:
    C_alpha(k) = sigma_alpha^2 / (1 - lambda_alpha(k)^2)

  The charge correlator is then:
    G(k) = sum_alpha |c.v_alpha(k)|^2 * sigma_alpha^2 / (1 - lambda_alpha(k)^2)

  The noise sigma_alpha^2 depends on details of the stochastic forcing.
  For a UNIVERSAL form factor (independent of noise), we look at:
    G_E(Q^2) = G(Q) / G(0)

  If sigma_alpha is the same for all modes, or if one mode dominates,
  the noise drops out of the ratio.

  DOMINANT MODE APPROXIMATION:
  At k=0, mode 0 has lambda_0 = 0.542 and |c.v_0|^2 = 2.75 (dominant).
  Mode 1 has |c.v_1|^2 = 0 (zero charge projection).
  So the charge form factor is determined by mode 0 alone:

    G_E(Q^2) = [1 - lambda_0(0)^2] / [1 - lambda_0(Q)^2]

  where lambda_0(Q) is the eigenvalue at spatial momentum Q.
""")

# Compute the charge form factor using the dominant mode approximation
def form_factor_dominant(q, direction='100'):
    """
    Compute G_E(q) = [1 - lam(0)^2] / [1 - lam(q)^2]
    for momentum q along the specified direction.
    """
    if direction == '100':
        gamma_q = (np.cos(q) + 2) / 3.0
    elif direction == '111':
        gamma_q = np.cos(q / np.sqrt(3))
    else:
        gamma_q = (np.cos(q) + 2) / 3.0

    lam_q = charge_eigenvalue(gamma_q)
    lam_0 = charge_eigenvalue(1.0)
    return (1 - lam_0**2) / (1 - lam_q**2)

# Compute full form factor using all modes
def form_factor_full(q):
    """
    G_E(q) using all 4 modes.
    q is along (1,0,0), so gamma = (cos q + 2)/3.
    """
    gamma_0 = 1.0
    gamma_q = (np.cos(q) + 2) / 3.0

    G0 = 0.0
    Gq = 0.0

    for gamma, G_ref in [(gamma_0, None), (gamma_q, None)]:
        M = Mk_gamma(gamma, J_m_np, J_e_np)
        evals, evecs = np.linalg.eig(M)
        G_val = 0.0
        for alpha in range(4):
            v = evecs[:, alpha]; v = v / np.linalg.norm(v)
            cv2 = abs(np.dot(c_vec, v))**2
            lam = evals[alpha]
            denom = 1.0 - float(np.real(lam))**2
            if abs(denom) > 1e-15:
                G_val += cv2 / denom
        if gamma == gamma_0:
            G0 = G_val
        else:
            Gq = G_val

    # Recompute properly
    M0 = Mk_gamma(gamma_0, J_m_np, J_e_np)
    evals0, evecs0 = np.linalg.eig(M0)
    G0 = 0.0
    for alpha in range(4):
        v = evecs0[:, alpha]; v = v / np.linalg.norm(v)
        cv2 = abs(np.dot(c_vec, v))**2
        lam = float(np.real(evals0[alpha]))
        denom = 1.0 - lam**2
        if abs(denom) > 1e-15:
            G0 += cv2 / denom

    Mq = Mk_gamma(gamma_q, J_m_np, J_e_np)
    evalsq, evecsq = np.linalg.eig(Mq)
    Gq = 0.0
    for alpha in range(4):
        v = evecsq[:, alpha]; v = v / np.linalg.norm(v)
        cv2 = abs(np.dot(c_vec, v))**2
        lam = float(np.real(evalsq[alpha]))
        denom = 1.0 - lam**2
        if abs(denom) > 1e-15:
            Gq += cv2 / denom

    return Gq / G0

# Print form factor values
print(f"\n  Form factor G_E(q) along (q,0,0):")
print(f"  {'q':>8} {'q_phys(fm^-1)':>14} {'GE_dom':>12} {'GE_full':>12}")
q_scan = np.linspace(0, 2.5, 50)
for q in q_scan:
    q_phys = q / a_fm  # physical momentum in fm^-1
    GE_d = form_factor_dominant(q)
    GE_f = form_factor_full(q)
    if q < 0.51 or abs(q - 1.0) < 0.06 or abs(q - 2.0) < 0.06:
        print(f"  {q:>8.4f} {q_phys:>14.4f} {GE_d:>12.8f} {GE_f:>12.8f}")


# ============================================================
# STEP 5: CHARGE RADIUS FROM FORM FACTOR SLOPE
# ============================================================
print("\n" + "=" * 80)
print("STEP 5: CHARGE RADIUS EXTRACTION")
print("=" * 80)

# Method A: Numerical derivative of G_E(q) at q=0
# r_p^2 = -6 * dG_E/d(Q_phys^2) |_0
# Q_phys = q/a, so Q_phys^2 = q^2/a^2
# dG_E/d(Q^2) = a^2 * dG_E/d(q^2)
# r_p^2 = -6 * a^2 * dG_E/d(q^2)|_0

# For the dominant mode:
# G_E(q) = (1-lam0^2)/(1-lam(q)^2)
# where lam(q) depends on gamma(q) = (cos q + 2)/3 ≈ 1 - q^2/6

# d G_E / d(q^2) = (1-lam0^2) * d/d(q^2) [1/(1-lam^2)]
#                = (1-lam0^2) * 2*lam * (dlam/d(q^2)) / (1-lam^2)^2
# At q=0: lam = lam0, and
#   dlam/d(q^2) = (dlam/dgamma) * (dgamma/d(q^2))
#               = (dlam/dgamma) * (-1/6)

# Compute dlam/dgamma at gamma=1
dg = 1e-6
lam_p = charge_eigenvalue(1.0 + dg)
lam_m = charge_eigenvalue(1.0 - dg)
dlam_dg = (lam_p - lam_m) / (2*dg)

lam0 = charge_eigenvalue(1.0)
dlam_dq2 = dlam_dg * (-1.0/6.0)

dGE_dq2 = (1-lam0**2) * 2*lam0 * dlam_dq2 / (1-lam0**2)**2
dGE_dq2 = 2*lam0 * dlam_dq2 / (1-lam0**2)

print(f"\n  Dominant mode analysis:")
print(f"    lambda_0(k=0) = {lam0:.10f}")
print(f"    1 - lambda_0^2 = {1-lam0**2:.10f}")
print(f"    d(lambda)/d(gamma) = {dlam_dg:.10f}")
print(f"    d(lambda)/d(q^2) = dlam_dg * (-1/6) = {dlam_dq2:.10f}")
print(f"    dG_E/d(q^2)|_0 = {dGE_dq2:.10f}")

r_sq_latt = -6 * dGE_dq2  # in lattice units (a^2)
r_sq_phys = r_sq_latt * a_fm**2  # in fm^2

print(f"\n    r^2 (lattice) = -6 * dG_E/d(q^2) = {r_sq_latt:.10f} (in units of a^2)")
print(f"    r^2 (physical) = r^2_latt * a^2 = {r_sq_phys:.10f} fm^2")
print(f"    r_p = {np.sqrt(abs(r_sq_phys)):.6f} fm")

# Method B: Numerical derivative directly
print(f"\n  Numerical check (central difference on G_E):")
for step in [0.01, 0.001, 0.0001]:
    GE_p = form_factor_dominant(step)
    GE_m = form_factor_dominant(-step)  # same as +step by symmetry
    d2GE = (GE_p - 2.0 + GE_m) / step**2
    # G_E(q) ≈ 1 + (1/2) G_E'' q^2, so dG_E/d(q^2) = G_E''/2
    # But for isotropic: dG_E/d(|k|^2) = (1/2) d^2G_E/dq^2 when scanning along one axis? No.
    # Along (q,0,0): |k|^2 = q^2, so d/d(|k|^2) = d/d(q^2)
    # And G_E(q) ≈ 1 + dG_E/d(q^2)|_0 * q^2 + ...
    # The Taylor: G_E(q) = 1 + (1/2) G_E''(0) q^2 + ...
    # So dG_E/d(q^2)|_0 = (1/2) G_E''(0)
    slope = d2GE / 2.0
    r2 = -6 * slope * a_fm**2
    print(f"    step={step}: G_E''={d2GE:.10f}, dGE/d(q^2)={slope:.10f}, r^2={r2:.6f} fm^2, r={np.sqrt(abs(r2)):.6f} fm")

# Method C: Full mode form factor numerical derivative
print(f"\n  Full-mode numerical check:")
for step in [0.01, 0.001, 0.0001]:
    GE_p = form_factor_full(step)
    GE_m = form_factor_full(-step)
    d2GE = (GE_p - 2.0 + GE_m) / step**2
    slope = d2GE / 2.0
    r2 = -6 * slope * a_fm**2
    print(f"    step={step}: dGE/d(q^2)={slope:.10f}, r={np.sqrt(abs(r2)):.6f} fm")

# Method D: Fit G_E to dipole form
print(f"\n  Dipole fit: G_E(q) = 1/(1 + q^2/Lambda^2)^2")
# For small q, dipole gives: G_E ≈ 1 - 2q^2/Lambda^2
# So dG_E/d(q^2) = -2/Lambda^2
# r^2 = -6 * (-2/Lambda^2) = 12/Lambda^2
# Lambda^2 (lattice) is found by fitting

q_fit_range = np.linspace(0.01, 1.0, 100)
GE_fit_data = np.array([form_factor_dominant(q) for q in q_fit_range])
# Fit: G_E = (1 + q^2/L2)^{-2}
# => (G_E)^{-1/2} = 1 + q^2/L2
# => (G_E^{-1/2} - 1) = q^2/L2
y_fit = GE_fit_data**(-0.5) - 1
x_fit = q_fit_range**2
# Linear fit: y = (1/L2) * x
L2_inv = np.polyfit(x_fit[:20], y_fit[:20], 1)[0]  # use small-q region
Lambda2_latt = 1.0 / L2_inv
r2_dipole_latt = 12.0 / Lambda2_latt  # Hmm, wait: for G_E = 1/(1+q^2/L2)^2,
                                        # dG_E/d(q^2)|_0 = -2/L2
                                        # r^2/a^2 = -6*(-2/L2) = 12/L2
r2_dipole_phys = r2_dipole_latt * a_fm**2

print(f"    Lambda^2 (lattice) = {Lambda2_latt:.6f}")
print(f"    r^2/a^2 (dipole) = 12/Lambda^2 = {r2_dipole_latt:.6f}")
print(f"    r_p (dipole) = {np.sqrt(abs(r2_dipole_phys)):.6f} fm")


# ============================================================
# STEP 6: ALTERNATIVE LATTICE SPACING METHODS
# ============================================================
print("\n" + "=" * 80)
print("STEP 6: ALTERNATIVE LATTICE SPACING METHODS")
print("=" * 80)

# Method I: From the pole (Step 3)
print(f"\n  Method I: Pole of spatial propagator")
print(f"    m_p * a = q_pole = {q_pole:.8f}")
print(f"    a = {a_fm:.6f} fm")

# Method II: From the paper's Delta_0
Delta_0_paper = 0.6888  # = -ln(0.5022)
a_paper = Delta_0_paper * hbar_c / m_proton  # This gives the mass of 0++ state
# Actually Delta_0 corresponds to the 0++ glueball, not the proton
# Proton mass: m_p = Delta_1/2 * scale = -ln(0.4745^{1/2}) * 1710/Delta_0
# In lattice units: m_p * a_paper_proton = -ln(0.4745^{1/2})
m_pa_paper = -np.log(0.4745**0.5)
a_paper_p = m_pa_paper * hbar_c / m_proton
print(f"\n  Method II: Paper proton eigenvalue")
print(f"    m_p * a = -ln(lambda_1^{{1/2}}) = {m_pa_paper:.8f}")
print(f"    a = {a_paper_p:.6f} fm")

# Method III: From lattice lambda_0 at k=0
# The 3D lattice eigenvalue at k=0 is lambda_0 = 0.5416
# The 1D (paper) eigenvalue is lambda_0 = 0.5022
# These differ because the 3D combines J_m + J_e (full coupling)
# The 1D A_2 had a different coupling structure
lam0_3d = charge_eigenvalue(1.0)
m_pa_3d = -np.log(lam0_3d)
a_3d = m_pa_3d * hbar_c / m_proton
print(f"\n  Method III: 3D k=0 eigenvalue")
print(f"    lambda_0(k=0) = {lam0_3d:.8f}")
print(f"    m_p * a = -ln(lambda_0) = {m_pa_3d:.8f}")
print(f"    a = {a_3d:.6f} fm")

# For the proton, the mass in lattice units is given by the pole (Method I).
# The k=0 eigenvalue gives the correlation time, not the spatial mass.
# The pole q_0 = arccosh(3*gamma_pole - 2) gives the spatial mass.


# ============================================================
# STEP 7: FINAL RESULTS AND COMPARISON
# ============================================================
print("\n" + "=" * 80)
print("STEP 7: FINAL RESULTS")
print("=" * 80)

# Using pole mass (Method I) for the lattice spacing:
# r^2 in lattice units (from dominant mode analytic)
r_latt = np.sqrt(abs(r_sq_latt))
r_phys = r_latt * a_fm

print(f"\n  DOMINANT MODE RESULT:")
print(f"    r/a = {r_latt:.6f}")
print(f"    a   = {a_fm:.6f} fm  (from spatial pole)")
print(f"    r_p = {r_phys:.6f} fm")
print(f"    r_p(exp) = {r_p_exp:.5f} fm")
print(f"    Ratio: {r_phys/r_p_exp:.4f}")

# Also try with paper lattice spacing
r_phys_paper = r_latt * a_paper_p
print(f"\n  WITH PAPER LATTICE SPACING:")
print(f"    a   = {a_paper_p:.6f} fm")
print(f"    r_p = {r_phys_paper:.6f} fm")
print(f"    Ratio: {r_phys_paper/r_p_exp:.4f}")

# What (H+1) * lambda_C gives in our framework
print(f"\n  ANALYTIC PREDICTION:")
print(f"    r_p = (H+1) * hbar_c / m_p = {(H+1)*lambda_C:.6f} fm")
print(f"    Ratio to exp: {(H+1)*lambda_C/r_p_exp:.6f}")

# Can we understand the discrepancy?
print(f"\n  DIAGNOSTIC:")
print(f"    The charge radius in lattice units: r/a = {r_latt:.6f}")
print(f"    The proton mass in lattice units: m_p*a = {q_pole:.6f}")
print(f"    Product r * m_p = r/a * m_p*a = {r_latt * q_pole:.6f}")
print(f"    Compare: (H+1) = {H+1}")
print(f"    Ratio: {r_latt * q_pole / (H+1):.6f}")
print(f"    The dimensionless product r_p * m_p / hbar_c should be H+1 = 4")
print(f"    Our lattice gives: {r_latt * q_pole:.4f}")

# The product r*m in lattice units is dimensionless.
# Experimentally: r_p * m_p / hbar_c = 0.84087 * 938.272 / 197.327 = 3.998
product_exp = r_p_exp * m_proton / hbar_c
print(f"    Experimental: r_p * m_p / hbar_c = {product_exp:.4f}")

# From the pole analysis:
# m_p * a = q_pole where cosh(q_pole) = 3*gamma_pole - 2
# The charge radius depends on how fast the form factor falls off.
# dG_E/d(q^2) is controlled by d(lambda)/d(gamma) and 1-lambda^2.

# More detailed: what r/a SHOULD be to match experiment
r_over_a_needed = r_p_exp / a_fm
print(f"\n    To match r_p = {r_p_exp} fm with a = {a_fm:.4f} fm: need r/a = {r_over_a_needed:.4f}")
print(f"    Current r/a = {r_latt:.4f}")

# Alternative: what a should be if r/a is correct
a_needed = r_p_exp / r_latt
print(f"    To match r_p = {r_p_exp} fm with r/a = {r_latt:.4f}: need a = {a_needed:.4f} fm")
print(f"    That would require m_p*a = {a_needed * m_proton / hbar_c:.4f}")


# ============================================================
# STEP 8: FORM FACTOR TABLE
# ============================================================
print("\n" + "=" * 80)
print("STEP 8: FORM FACTOR TABLE")
print("=" * 80)

print(f"\n  {'Q^2 (fm^-2)':>12} {'Q^2 (GeV^2)':>12} {'GE_DS':>10} {'GE_dipole':>10}")
print(f"  {'-'*12} {'-'*12} {'-'*10} {'-'*10}")

# Standard dipole: G_E = 1/(1 + Q^2/0.71)^2 with Q in GeV
for Q2_gev in [0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
    Q2_fm = Q2_gev / (hbar_c * 1e-3)**2 * 1e-6  # GeV^2 to fm^-2 ... no
    # Q in GeV, Q^2 in GeV^2
    # Q in fm^-1: Q_fm = Q_GeV * 1000 / hbar_c  (since hbar_c = 197.3 MeV fm)
    # So Q^2_fm = Q^2_GeV * (1000/197.3)^2
    Q2_fm = Q2_gev * (1000.0/hbar_c)**2
    # Lattice momentum: q = Q_phys * a
    q_latt = np.sqrt(Q2_fm) * a_fm if Q2_fm > 0 else 0
    GE_ds = form_factor_dominant(q_latt) if q_latt < np.pi else float('nan')
    GE_dip = 1.0 / (1 + Q2_gev / 0.71)**2
    print(f"  {Q2_fm:>12.4f} {Q2_gev:>12.4f} {GE_ds:>10.6f} {GE_dip:>10.6f}")


# ============================================================
# STEP 9: ANALYTIC UNDERSTANDING
# ============================================================
print("\n" + "=" * 80)
print("STEP 9: ANALYTIC UNDERSTANDING")
print("=" * 80)

print(f"""
  KEY QUANTITIES:
    lambda_0 = {lam0:.10f}  (charge eigenvalue at k=0)
    d(lambda)/d(gamma) = {dlam_dg:.10f}
    1 - lambda_0^2 = {1-lam0**2:.10f}

  The form factor slope in the dominant mode approximation:
    dG_E/d(q^2)|_0 = (2*lam0*dlam_dq2) / (1-lam0^2)
                   = (2*{lam0:.6f}*{dlam_dq2:.6f}) / {1-lam0**2:.6f}
                   = {dGE_dq2:.10f}

  Charge radius in lattice units:
    r^2/a^2 = -6 * dG_E/d(q^2) = {r_sq_latt:.10f}

  Spatial pole mass:
    m_p * a = {q_pole:.10f}

  Dimensionless product:
    r_p * m_p = sqrt(r^2/a^2) * (m_p*a) = {r_latt:.6f} * {q_pole:.6f} = {r_latt*q_pole:.6f}
    (Should be ~4 for experiment)

  THE DOMINANT MODE FORMULA gives:
    r^2 * m^2 = r^2/a^2 * (m*a)^2

  where m*a = arccosh(3*gamma_pole - 2) and gamma_pole satisfies lambda(gamma_pole) = 1.

  For a LINEAR eigenvalue: lambda(gamma) = lambda_0 + dlam_dg * (gamma - 1)
    gamma_pole = 1 + (1 - lambda_0) / dlam_dg
    cosh(q_pole) = 3*gamma_pole - 2 = 1 + 3*(1-lambda_0)/dlam_dg

  And:
    r^2/a^2 = -6 * (2*lam0 * dlam_dg * (-1/6)) / (1-lam0^2)
            = 2*lam0*dlam_dg / (1-lam0^2)
            = 2*lam0*dlam_dg / ((1-lam0)(1+lam0))

  Let's check with numbers:
    gamma_pole = 1 + (1 - {lam0:.6f}) / {dlam_dg:.6f} = {1 + (1-lam0)/dlam_dg:.6f}
    vs actual: {gamma_pole:.6f}
""")

gamma_pole_lin = 1 + (1 - lam0) / dlam_dg
print(f"  Linear approximation for gamma_pole: {gamma_pole_lin:.6f} vs actual: {gamma_pole:.6f}")
print(f"  (Error: {abs(gamma_pole_lin - gamma_pole)/gamma_pole * 100:.2f}%)")

# Full analytic expression for r*m
r2_a2 = 2*lam0*dlam_dg / (1-lam0**2)
cosh_q = 1 + 3*(1-lam0)/dlam_dg
q_lin = np.arccosh(cosh_q)
rm_product = np.sqrt(r2_a2) * q_lin
print(f"\n  r^2/a^2 (analytic) = {r2_a2:.10f}")
print(f"  vs numerical: {r_sq_latt:.10f}")
print(f"  m*a (linear approx) = {q_lin:.10f}")
print(f"  vs actual: {q_pole:.10f}")
print(f"  r*m product (analytic) = {rm_product:.6f}")
print(f"  vs needed: {H+1}")

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)
