"""
Black hole maximum density from H=3 informational geometry.

The Born floor halts gravitational collapse at a maximum density
determined by the fibre curvature of the map Phi = floor . DS
at the K*=7/30 fixed point.

rho_max / rho_Planck = H^2 * ||dbar Phi|| / (8 pi)

where ||dbar Phi|| = operator norm of the second derivative (Hessian)
of Phi at the fixed point m*.
"""

from mpmath import (mp, mpf, matrix, sqrt, fabs, nstr, eig, chop,
                    findroot, eye, svd, pi, log10)

mp.dps = 40
H = 3
FLOOR = mpf(1) / mpf(H**3)   # 1/27

print("=" * 70)
print("BLACK HOLE MAXIMUM DENSITY FROM H=3 INFORMATIONAL GEOMETRY")
print("=" * 70)
print(f"\nH = {H}")
print(f"Born floor = 1/H^3 = 1/{H**3} = {nstr(FLOOR, 15)}")
print(f"K* = 7/30 = {nstr(mpf(7)/mpf(30), 15)}")

# ─────────────────────────────────────────────────────────────────────
# DS combination and Born floor enforcement
# ─────────────────────────────────────────────────────────────────────

def ds_combine(m, e):
    """Dempster-Shafer combination of mass functions m and e."""
    s = m[:3]; theta = m[3]; ev = e[:3]; phi = e[3]
    s_pre = [s[i]*ev[i] + s[i]*phi + theta*ev[i] for i in range(3)]
    theta_pre = theta * phi
    total_pre = sum(s_pre) + theta_pre
    K = mpf(1) - total_pre
    denom = mpf(1) - K
    return [sp/denom for sp in s_pre] + [theta_pre/denom], K

def born_prob(m):
    """Born probability of the ignorance component."""
    L2sq = sum(x**2 for x in m)
    if L2sq < mpf(10)**(-60):
        return mpf(0)
    return m[3]**2 / L2sq

def enforce_floor(m):
    """Enforce Born floor: redistribute if theta Born prob < 1/27."""
    b = born_prob(m)
    if b >= FLOOR - mpf(10)**(-60):
        return list(m)
    S = sum(m[:3])
    Sq = sum(x**2 for x in m[:3])
    # Solve for t such that t^2 / (Sq_rescaled + t^2) = 1/27
    # where rescaled s_i = m_i * (1-t)/S
    # => Sq_rescaled = Sq * ((1-t)/S)^2
    # => t^2 / (Sq*((1-t)/S)^2 + t^2) = 1/27
    # => 27*t^2 = Sq*((1-t)/S)^2 + t^2
    # => 26*t^2 = Sq*(1-t)^2/S^2
    # => 26*t^2*S^2 = Sq*(1-2t+t^2)
    # => (26*S^2 - Sq)*t^2 + 2*Sq*t - Sq = 0
    A_c = mpf(26) * S**2 - Sq
    B_c = mpf(2) * Sq
    C_c = -Sq
    disc = B_c**2 - 4 * A_c * C_c
    t1 = (-B_c + sqrt(disc)) / (2 * A_c)
    t2 = (-B_c - sqrt(disc)) / (2 * A_c)
    cands = [t for t in [t1, t2] if mpf(0) < t < mpf(1)]
    t = min(cands, key=lambda x: fabs(x - m[3]))
    alpha_r = (mpf(1) - t) / S
    return [m[i] * alpha_r for i in range(3)] + [t]

def phi_map(m, e):
    """Full map: DS combine then enforce Born floor."""
    m_ds, K = ds_combine(m, e)
    return enforce_floor(m_ds)

# ─────────────────────────────────────────────────────────────────────
# Step 1: Find the fixed point m* at K* = 7/30
# ─────────────────────────────────────────────────────────────────────
print("\n" + "─" * 70)
print("STEP 1: Finding fixed point m* at K* = 7/30")
print("─" * 70)

def eq_system(*p):
    s1, theta, w1, phi = p
    s2 = (mpf(1) - s1 - theta) / 2
    w2 = (mpf(1) - w1 - phi) / 2
    m = [s1, s2, s2, theta]
    e = [w1, w2, w2, phi]
    # Born floor saturated for both state and evidence
    eq1 = theta**2 / (s1**2 + 2*s2**2 + theta**2) - FLOOR
    eq2 = phi**2 / (w1**2 + 2*w2**2 + phi**2) - FLOOR
    # Conflict = 7/30
    m_ds, K_val = ds_combine(m, e)
    eq3 = K_val - mpf(7) / mpf(30)
    # Fixed point condition: Phi(m,e) = m
    m_out = phi_map(m, e)
    eq4 = m_out[0] - s1
    return [eq1, eq2, eq3, eq4]

sol = findroot(eq_system,
               [mpf('0.787'), mpf('0.155'), mpf('0.631'), mpf('0.128')])
s1s, ths, w1s, phis = sol
s2s = (mpf(1) - s1s - ths) / 2
w2s = (mpf(1) - w1s - phis) / 2
m_star = [s1s, s2s, s2s, ths]
e_star = [w1s, w2s, w2s, phis]

print(f"\nm* = ({nstr(s1s,15)}, {nstr(s2s,15)}, {nstr(s2s,15)}, {nstr(ths,15)})")
print(f"e* = ({nstr(w1s,15)}, {nstr(w2s,15)}, {nstr(w2s,15)}, {nstr(phis,15)})")
print(f"Sum m* = {nstr(sum(m_star),15)}")
print(f"Sum e* = {nstr(sum(e_star),15)}")

# Verify
m_out = phi_map(m_star, e_star)
_, K_check = ds_combine(m_star, e_star)
b_check = born_prob(m_star)
print(f"\nVerification:")
print(f"  K = {nstr(K_check, 15)}  (target: {nstr(mpf(7)/30, 15)})")
print(f"  Born(m*) = {nstr(b_check, 15)}  (target: {nstr(FLOOR, 15)})")
print(f"  Phi(m*) = ({nstr(m_out[0],12)}, {nstr(m_out[1],12)}, {nstr(m_out[2],12)}, {nstr(m_out[3],12)})")
print(f"  |Phi(m*) - m*| = {nstr(sqrt(sum((m_out[i]-m_star[i])**2 for i in range(4))), 6)}")

# ─────────────────────────────────────────────────────────────────────
# Step 2: Compute Jacobian J (4x4) via central differences
# ─────────────────────────────────────────────────────────────────────
print("\n" + "─" * 70)
print("STEP 2: Jacobian of Phi at m*")
print("─" * 70)

def phi_vec(m_list):
    """Phi as a function of m only, with e fixed at e*."""
    return phi_map(list(m_list), e_star)

eps_J = mpf(10)**(-20)
J = matrix(4, 4)
for j in range(4):
    m_plus = list(m_star)
    m_minus = list(m_star)
    m_plus[j] += eps_J
    m_minus[j] -= eps_J
    f_plus = phi_vec(m_plus)
    f_minus = phi_vec(m_minus)
    for i in range(4):
        J[i, j] = (f_plus[i] - f_minus[i]) / (2 * eps_J)

print("\nJacobian J = dPhi/dm at m*:")
for i in range(4):
    row = [nstr(J[i, j], 10) for j in range(4)]
    print(f"  [{', '.join(row)}]")

# Eigenvalues of J
evals_J = eig(J, left=False, right=False)
print("\nEigenvalues of J:")
for ev in evals_J:
    ev_c = chop(ev, tol=mpf(10)**(-15))
    print(f"  {nstr(ev_c, 12)}  (|lambda| = {nstr(fabs(ev_c), 12)})")

spectral_radius = max(fabs(chop(ev, tol=mpf(10)**(-15))) for ev in evals_J)
print(f"\nSpectral radius of J = {nstr(spectral_radius, 15)}")

# ─────────────────────────────────────────────────────────────────────
# Step 3: Compute Hessian tensors H^(k) via finite differences
# ─────────────────────────────────────────────────────────────────────
print("\n" + "─" * 70)
print("STEP 3: Hessian tensors (second derivatives) of Phi at m*")
print("─" * 70)

eps_H = mpf(10)**(-12)

# H^(k)_{ij} = d^2 Phi_k / dm_i dm_j
# Use 4-point central difference:
# f_{++} - f_{+-} - f_{-+} + f_{--}  /  (4 * eps^2)
Hess = []  # List of 4 matrices, each 4x4
for k in range(4):
    Hk = matrix(4, 4)
    for i in range(4):
        for j in range(4):
            m_pp = list(m_star); m_pp[i] += eps_H; m_pp[j] += eps_H
            m_pm = list(m_star); m_pm[i] += eps_H; m_pm[j] -= eps_H
            m_mp = list(m_star); m_mp[i] -= eps_H; m_mp[j] += eps_H
            m_mm = list(m_star); m_mm[i] -= eps_H; m_mm[j] -= eps_H
            f_pp = phi_vec(m_pp)[k]
            f_pm = phi_vec(m_pm)[k]
            f_mp = phi_vec(m_mp)[k]
            f_mm = phi_vec(m_mm)[k]
            Hk[i, j] = (f_pp - f_pm - f_mp + f_mm) / (4 * eps_H**2)
    Hess.append(Hk)

for k in range(4):
    print(f"\nH^({k}) = d^2 Phi_{k} / dm_i dm_j:")
    for i in range(4):
        row = [nstr(Hess[k][i, j], 8) for j in range(4)]
        print(f"  [{', '.join(row)}]")

# ─────────────────────────────────────────────────────────────────────
# Step 4: Operator norm ||dbar Phi|| = sqrt(max eig of sum H^(k)^T H^(k))
# ─────────────────────────────────────────────────────────────────────
print("\n" + "─" * 70)
print("STEP 4: Fibre curvature ||dbar Phi||")
print("─" * 70)

# Build M = sum_k H^(k)^T H^(k)
M = matrix(4, 4)
for i in range(4):
    for j in range(4):
        M[i, j] = mpf(0)

for k in range(4):
    HkT = Hess[k].T
    prod = HkT * Hess[k]
    for i in range(4):
        for j in range(4):
            M[i, j] += prod[i, j]

print("\nM = sum_k H^(k)^T H^(k):")
for i in range(4):
    row = [nstr(M[i, j], 8) for j in range(4)]
    print(f"  [{', '.join(row)}]")

evals_M = eig(M, left=False, right=False)
print("\nEigenvalues of M:")
for ev in evals_M:
    ev_r = chop(ev, tol=mpf(10)**(-10))
    print(f"  {nstr(ev_r, 15)}")

max_eig_M = max(chop(ev, tol=mpf(10)**(-10)) for ev in evals_M)
dbar_phi = sqrt(max_eig_M)

print(f"\nmax eigenvalue of M = {nstr(max_eig_M, 15)}")
print(f"||dbar Phi|| = sqrt(max eig) = {nstr(dbar_phi, 15)}")

# ─────────────────────────────────────────────────────────────────────
# Step 5: Maximum density
# ─────────────────────────────────────────────────────────────────────
print("\n" + "─" * 70)
print("STEP 5: Maximum black hole density")
print("─" * 70)

rho_ratio = H**2 * dbar_phi / (8 * pi)

print(f"\nrho_max / rho_Planck = H^2 * ||dbar Phi|| / (8 pi)")
print(f"                     = {H**2} * {nstr(dbar_phi, 10)} / (8 pi)")
print(f"                     = {nstr(mpf(H**2) * dbar_phi, 12)} / {nstr(8*pi, 12)}")
print(f"                     = {nstr(rho_ratio, 15)}")

print(f"\n{'=' * 70}")
print(f"RESULT:  rho_max = {nstr(rho_ratio, 10)} rho_Planck")
print(f"{'=' * 70}")

# Comparison with LQC
rho_LQC = mpf('0.41')
print(f"\nComparison with Loop Quantum Cosmology:")
print(f"  LQC:       rho_max = {nstr(rho_LQC, 4)} rho_Planck")
print(f"  This work: rho_max = {nstr(rho_ratio, 10)} rho_Planck")
print(f"  Ratio (this/LQC) = {nstr(rho_ratio / rho_LQC, 10)}")
print(f"  Difference = {nstr(fabs(rho_ratio - rho_LQC) / rho_LQC * 100, 6)}%")

# Physical density
rho_Planck_SI = mpf('5.155e96')  # kg/m^3
rho_max_SI = rho_ratio * rho_Planck_SI
print(f"\nIn SI units (rho_Planck = 5.155e96 kg/m^3):")
print(f"  rho_max = {nstr(rho_max_SI, 6)} kg/m^3")
print(f"  log10(rho_max) = {nstr(log10(rho_max_SI), 6)}")
