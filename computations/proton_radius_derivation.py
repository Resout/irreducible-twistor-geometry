#!/usr/bin/env python3
"""
PROTON CHARGE RADIUS DERIVATION FROM H=3 AXIOMS

Observation: r_p = (H+1) × ℏ/(m_p c) = 4 × 0.2103 fm = 0.8412 fm
             CODATA 2018: 0.8414(19) fm,  muonic hydrogen: 0.84087(39) fm

Goal: derive r_p² = -6 dG_E/dq²|₀ = (H+1)² (ℏ/m_p c)² from the DS dynamics.

The transfer-operator approach:
  The A₂ coupled system gives an 8×8 Jacobian with structure
  [[J_self, J_cross], [J_cross, J_self]] (by exchange symmetry).
  Fourier transform: T(k) = J_self + J_cross e^{ik} for k ∈ [0,2π].

  At k=0: T(0) = J_self + J_cross (bonding modes)
  At k=π: T(π) = J_self - J_cross (antibonding modes)

  The A₂ eigenvalues split: two are bonding (at k=0), two antibonding (at k=π).
  Each "band" connects a k=0 eigenvalue to a k=π eigenvalue as k varies.

  The proton form factor G_E(q²) relates to the momentum-dependent eigenvalue:
  r_p² = -6 a² d²ln|λ_proton(k)|/dk² |_{k=k₀}
  where k₀ is the momentum at which the proton band sits.
"""

import numpy as np
from mpmath import mp, mpf, matrix, sqrt, log, ln, fabs, nstr, eig, chop, exp, pi
from mpmath import findroot, mpc

mp.dps = 50
H = 3
FLOOR = mpf(1) / mpf(H**3)
K_STAR = mpf(7) / 30

# Physical constants
hbar_c = 197.3269804   # MeV·fm
m_proton = 938.272046   # MeV
lambda_C = hbar_c / m_proton  # proton Compton wavelength in fm
r_p_codata = 0.8414     # fm (CODATA 2018)
r_p_muonic = 0.84087    # fm (muonic hydrogen)
r_p_target = (H + 1) * lambda_C  # = 4 × λ_C

print("=" * 80)
print("PROTON CHARGE RADIUS DERIVATION FROM H=3 FRAMEWORK")
print("=" * 80)
print(f"\nPhysical constants:")
print(f"  ℏc = {hbar_c} MeV·fm")
print(f"  m_p = {m_proton} MeV")
print(f"  λ_C = ℏ/(m_p c) = {lambda_C:.6f} fm")
print(f"  r_p(CODATA 2018) = {r_p_codata} fm")
print(f"  r_p(muonic)      = {r_p_muonic} fm")
print(f"  r_p = (H+1)×λ_C  = {r_p_target:.6f} fm")
print(f"  Error vs CODATA:  {abs(1-r_p_target/r_p_codata)*100:.3f}%")
print(f"  Error vs muonic:  {abs(1-r_p_target/r_p_muonic)*100:.3f}%")


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
    out = [sp / denom for sp in s_pre] + [theta_pre / denom]
    return out, K

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
# 1. FIND EQUILIBRIUM
# ============================================================

def find_equilibrium():
    from scipy.optimize import fsolve
    def eq_rough(p):
        s1, theta, w1, phi = p
        s2 = (1.0 - s1 - theta) / 2.0
        w2 = (1.0 - w1 - phi) / 2.0
        m = [s1, s2, s2, theta]; e = [w1, w2, w2, phi]
        L2m = s1**2 + 2*s2**2 + theta**2
        L2e = w1**2 + 2*w2**2 + phi**2
        m_mp = [mpf(x) for x in m]; e_mp = [mpf(x) for x in e]
        m_out, _ = ds_step(m_mp, e_mp)
        return [theta**2/L2m - 1/27, phi**2/L2e - 1/27,
                s1*w2 + s1*w2 + s2*w1 + s2*w2 + s2*w1 + s2*w2 - 7/30,
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


print("\n" + "=" * 80)
print("STEP 1: EQUILIBRIUM AND JACOBIANS")
print("=" * 80)

m_star, e_star = find_equilibrium()
print(f"\n  m* = [{', '.join(nstr(x, 15) for x in m_star)}]")
print(f"  e* = [{', '.join(nstr(x, 15) for x in e_star)}]")

# Single-site 4×4 Jacobian
eps = mpf(10)**(-35)
J_m = matrix(4, 4)
for j in range(4):
    mp_p = list(m_star); mp_m = list(m_star)
    mp_p[j] += eps; mp_m[j] -= eps
    fp, _ = ds_step(mp_p, e_star)
    fm, _ = ds_step(mp_m, e_star)
    for i in range(4):
        J_m[i,j] = (fp[i] - fm[i]) / (2*eps)

evals_m, _ = eig(J_m)
evals_m = sorted([chop(e) for e in evals_m], key=lambda x: -float(fabs(x)))
print(f"\n  Single-site J_m eigenvalues:")
for k, ev in enumerate(evals_m):
    print(f"    λ_{k}^(ss) = {nstr(ev, 18)}  |λ| = {nstr(fabs(ev), 12)}")


# ============================================================
# 2. COUPLED A₂ SYSTEM
# ============================================================

print("\n" + "=" * 80)
print("STEP 2: COUPLED A₂ SYSTEM")
print("=" * 80)

def coupled_evidence(m1, m2, e0, g):
    uniform = [mpf(1)/4]*4
    e1 = [e0[i] - g * (m2[i] - uniform[i]) for i in range(4)]
    e1 = [max(x, mpf(10)**(-50)) for x in e1]
    total = sum(e1)
    return [x / total for x in e1]

def coupled_A2_step(state, e0, g):
    m1 = state[:4]; m2 = state[4:]
    e1 = coupled_evidence(m1, m2, e0, g)
    e2 = coupled_evidence(m2, m1, e0, g)
    m1_new, _ = ds_step(m1, e1)
    m2_new, _ = ds_step(m2, e2)
    return m1_new + m2_new

g = K_STAR
state = list(m_star) + list(m_star)
for i in range(30000):
    state_new = coupled_A2_step(state, e_star, g)
    diff_val = max(fabs(state_new[j] - state[j]) for j in range(8))
    if diff_val < mpf(10)**(-45):
        print(f"  Converged at step {i}")
        break
    state = state_new

state_star = state_new
m1_star = state_star[:4]
m2_star = state_star[4:]
print(f"  m1* = [{', '.join(nstr(x, 12) for x in m1_star)}]")
print(f"  m2* = [{', '.join(nstr(x, 12) for x in m2_star)}]")

# 8×8 Jacobian
J8 = matrix(8, 8)
for j in range(8):
    sp = list(state_star); sm = list(state_star)
    sp[j] += eps; sm[j] -= eps
    fp = coupled_A2_step(sp, e_star, g)
    fm = coupled_A2_step(sm, e_star, g)
    for i in range(8):
        J8[i,j] = (fp[i] - fm[i]) / (2*eps)

evals8, _ = eig(J8)
evals8 = [chop(e) for e in evals8]
evals8_sorted = sorted(evals8, key=lambda x: -float(fabs(x)))

print(f"\n  A₂ eigenvalues:")
for k, ev in enumerate(evals8_sorted):
    print(f"    Λ_{k} = {nstr(ev, 18):>25s}  |Λ| = {nstr(fabs(ev), 12)}")

# Extract J_self and J_cross from the 8×8 block structure
J_self = matrix(4, 4)
J_cross = matrix(4, 4)
for i in range(4):
    for j in range(4):
        J_self[i,j] = J8[i,j]
        J_cross[i,j] = J8[i,j+4]

# Verify block structure (should be [[A,B],[B,A]] by symmetry)
sym_err = max(fabs(J8[i+4,j+4] - J8[i,j]) for i in range(4) for j in range(4))
cross_err = max(fabs(J8[i+4,j] - J8[i,j+4]) for i in range(4) for j in range(4))
print(f"\n  Block symmetry check:")
print(f"    ||J_self(1,1) - J_self(2,2)|| = {float(sym_err):.2e}")
print(f"    ||J_cross(2,1) - J_cross(1,2)|| = {float(cross_err):.2e}")

# ============================================================
# 3. FOURIER DECOMPOSITION: T(k) = J_self + J_cross e^{ik}
# ============================================================

print("\n" + "=" * 80)
print("STEP 3: DISPERSION RELATION T(k)")
print("=" * 80)

# Eigenvalues at k=0 and k=π
T0 = J_self + J_cross
Tpi = J_self - J_cross
evals_T0, _ = eig(T0)
evals_Tpi, _ = eig(Tpi)
evals_T0 = sorted([chop(e) for e in evals_T0], key=lambda x: -float(fabs(x)))
evals_Tpi = sorted([chop(e) for e in evals_Tpi], key=lambda x: -float(fabs(x)))

print(f"\n  T(k=0) eigenvalues (bonding):")
for e in evals_T0:
    print(f"    {nstr(e, 15)}")
print(f"  T(k=π) eigenvalues (antibonding):")
for e in evals_Tpi:
    print(f"    {nstr(e, 15)}")

# The A₂ eigenvalues are the FULL SET from both k=0 and k=π
all_boundary = sorted([float(fabs(e)) for e in evals_T0 + evals_Tpi], reverse=True)
print(f"\n  All boundary eigenvalues (sorted): {[f'{x:.6f}' for x in all_boundary[:6]]}")
print(f"  Paper A₂ eigenvalues:              [0.502200, 0.474500, 0.352700, 0.334400]")

# Compute full dispersion relation with COMPLEX transfer matrix
n_k = 1001
k_grid = np.linspace(0, np.pi, n_k)
dk = k_grid[1] - k_grid[0]

# For each k, compute eigenvalues of T(k) = J_self + J_cross * e^{ik}
# and track |λ(k)| for each band
bands = [[] for _ in range(4)]

for k in k_grid:
    Tk = matrix(4, 4)
    ck = mpf(str(np.cos(k)))
    sk = mpf(str(np.sin(k)))
    for i in range(4):
        for j in range(4):
            Tk[i,j] = J_self[i,j] + J_cross[i,j] * mpc(ck, sk)
    evals_k, _ = eig(Tk)
    mods = sorted([float(fabs(e)) for e in evals_k], reverse=True)
    for b in range(4):
        bands[b].append(mods[b])

print(f"\n  Dispersion computed on {n_k} points in [0,π].")
print(f"  Band endpoints:")
for b in range(4):
    print(f"    Band {b}: λ(0)={bands[b][0]:.8f} → λ(π)={bands[b][-1]:.8f}")
    print(f"             bandwidth = {abs(bands[b][-1]-bands[b][0]):.8f}")

# Identify which band is the proton band
# The proton is λ₁^{1/2} with |λ₁| ≈ 0.4745
# This eigenvalue should appear at one of the band edges
print(f"\n  Identifying proton band:")
for b in range(4):
    for label, val in [("k=0", bands[b][0]), ("k=π", bands[b][-1])]:
        if abs(val - 0.4745) < 0.02:
            print(f"    Band {b} at {label}: {val:.6f} ≈ λ₁=0.4745")

# The proton eigenvalue 0.4745 appears at k=π for the antibonding mode
# The corresponding band runs from some value at k=0 to 0.4745 at k=π
# Or from 0.4745 at k=0 to some value at k=π

# Let's find which band has 0.4745 at either endpoint
proton_band = None
proton_at_kpi = False
for b in range(4):
    if abs(bands[b][-1] - 0.4745) < 0.02:
        proton_band = b
        proton_at_kpi = True
        break
    if abs(bands[b][0] - 0.4745) < 0.02:
        proton_band = b
        proton_at_kpi = False
        break

if proton_band is None:
    # Fallback: find band closest to 0.4745 at any point
    min_dist = 1e10
    for b in range(4):
        for val in bands[b]:
            d = abs(val - 0.4745)
            if d < min_dist:
                min_dist = d
                proton_band = b
    print(f"    Fallback: using band {proton_band}")
else:
    loc = "k=π" if proton_at_kpi else "k=0"
    print(f"    Proton band = {proton_band} (0.4745 at {loc})")


# ============================================================
# 4. CURVATURE AND CHARGE RADIUS
# ============================================================

print("\n" + "=" * 80)
print("STEP 4: CURVATURE → CHARGE RADIUS")
print("=" * 80)

# The charge radius comes from the curvature of the form factor.
# For a state at the Brillouin zone edge k=π, the relevant curvature
# is d²ln|λ|/dk² evaluated at k=π.
# For a state at k=0, evaluate at k=0.

# Compute curvatures at BOTH k=0 and k=π for all bands
print(f"\n  Curvatures d²ln|λ|/dk²:")

curvatures_k0 = []
curvatures_kpi = []

for b in range(4):
    ln_b = [np.log(max(x, 1e-30)) for x in bands[b]]

    # At k=0 (index 0), using symmetry λ(-k)=λ(k): f[-j]=f[j]
    # d²f/dk² = (-2f[2] + 32f[1] - 30f[0]) / (12dk²)
    d2_k0 = (-2*ln_b[2] + 32*ln_b[1] - 30*ln_b[0]) / (12 * dk**2)
    curvatures_k0.append(d2_k0)

    # At k=π (last index), using symmetry about π: f[N+j]=f[N-j]
    N = len(ln_b) - 1
    d2_kpi = (-2*ln_b[N-2] + 32*ln_b[N-1] - 30*ln_b[N]) / (12 * dk**2)
    curvatures_kpi.append(d2_kpi)

    print(f"    Band {b}: at k=0: {d2_k0:+.8f},  at k=π: {d2_kpi:+.8f}")

# Framework mass and lattice parameters
m_0pp = 1710.0  # MeV
lam_paper = [0.5022, 0.4745, 0.3527, 0.3344]
Delta_0_paper = -np.log(lam_paper[0])  # 0.6888
scale = m_0pp / Delta_0_paper  # 2482.6 MeV

# Lattice spacing: a = ℏc / scale
a = hbar_c / scale
print(f"\n  Scale = m₀/Δ₀ = {scale:.2f} MeV")
print(f"  Lattice spacing a = ℏc/scale = {a:.6f} fm")
print(f"  λ_C = ℏ/(m_p c) = {lambda_C:.6f} fm")
print(f"  a/λ_C = {a/lambda_C:.6f}")

# Compute radius for each band at relevant endpoint
print(f"\n  Physical radii (r = √(6 a² |curv|)):")
for b in range(4):
    # At k=0
    r_sq_0 = 6 * a**2 * abs(curvatures_k0[b])
    r_0 = np.sqrt(max(r_sq_0, 0))
    # At k=π
    r_sq_pi = 6 * a**2 * abs(curvatures_kpi[b])
    r_pi = np.sqrt(max(r_sq_pi, 0))

    print(f"    Band {b}: r(k=0) = {r_0:.6f} fm [{r_0/lambda_C:.4f} λ_C],"
          f"  r(k=π) = {r_pi:.6f} fm [{r_pi/lambda_C:.4f} λ_C]")

# The proton band curvature at the relevant k-point
if proton_band is not None:
    if proton_at_kpi:
        curv_proton = curvatures_kpi[proton_band]
        k_label = "k=π"
    else:
        curv_proton = curvatures_k0[proton_band]
        k_label = "k=0"

    r_sq_proton = 6 * a**2 * abs(curv_proton)
    r_proton = np.sqrt(max(r_sq_proton, 0))

    print(f"\n  === PROTON (band {proton_band}, {k_label}) ===")
    print(f"  d²ln|λ|/dk² = {curv_proton:.8f}")
    print(f"  r_p = {r_proton:.6f} fm")
    print(f"  r_p / λ_C = {r_proton/lambda_C:.6f}")
    print(f"  r_p / (H+1)λ_C = {r_proton/r_p_target:.6f}")
    print(f"  r_p / r_CODATA = {r_proton/r_p_codata:.6f}")


# ============================================================
# 5. ALTERNATIVE LATTICE SPACING CONVENTIONS
# ============================================================

print("\n" + "=" * 80)
print("STEP 5: EXPLORING LATTICE SPACING CONVENTIONS")
print("=" * 80)

# The key question: what is the correct physical lattice spacing?
# Convention 1: a = ℏc / (m₀/Δ₀) = ℏc·Δ₀/m₀ [used above]
# Convention 2: a = ℏc / m₀  [the glueball Compton wavelength]
# Convention 3: a = ℏc / (m_p × gap_proton) [proton correlation length]
# Convention 4: a set by requiring proton mass comes out right

a_conventions = {
    "a₁ = ℏcΔ₀/m₀": hbar_c * Delta_0_paper / m_0pp,
    "a₂ = ℏc/m₀": hbar_c / m_0pp,
    "a₃ = ℏc/m_p": hbar_c / m_proton,  # = λ_C
    "a₄ = ℏc/(m₀√Δ₀)": hbar_c / (m_0pp * np.sqrt(Delta_0_paper)),
}

print(f"\n  Lattice spacing conventions:")
for name, a_val in a_conventions.items():
    print(f"    {name} = {a_val:.6f} fm")

# For each convention, what radius does the proton band give?
if proton_band is not None:
    print(f"\n  Proton radius under each convention (band {proton_band}, {k_label}):")
    for name, a_val in a_conventions.items():
        r_val = np.sqrt(6 * a_val**2 * abs(curv_proton))
        print(f"    {name}: r_p = {r_val:.6f} fm, r_p/λ_C = {r_val/lambda_C:.4f}")

# What lattice spacing is REQUIRED to get r_p = (H+1)λ_C?
if proton_band is not None and abs(curv_proton) > 1e-10:
    a_needed = r_p_target / np.sqrt(6 * abs(curv_proton))
    scale_needed = hbar_c / a_needed
    print(f"\n  Required lattice spacing for r_p = (H+1)λ_C:")
    print(f"    a_needed = {a_needed:.6f} fm")
    print(f"    scale_needed = ℏc/a = {scale_needed:.2f} MeV")
    print(f"    Ratio a_needed/a₁ = {a_needed/a:.6f}")

    # Check if scale_needed matches any framework quantity
    tests = {
        "m₀/Δ₀": scale,
        "m₀": m_0pp,
        "m_p": m_proton,
        "m₀/(H+1)": m_0pp/(H+1),
        "m₀/H": m_0pp/H,
        "m₀×Δ₀": m_0pp * Delta_0_paper,
        "m₀√(Δ₀)": m_0pp * np.sqrt(Delta_0_paper),
        "m₀/√(Δ₀)": m_0pp / np.sqrt(Delta_0_paper),
        "m_p/Δ₀": m_proton / Delta_0_paper,
        "m_p×(H+1)": m_proton * (H+1),
        "m₀×K*": m_0pp * 7/30,
    }

    print(f"\n  Scale needed ({scale_needed:.2f} MeV) vs framework scales:")
    results = [(abs(scale_needed - v)/scale_needed*100, name, v) for name, v in tests.items()]
    results.sort()
    for err, name, val in results[:6]:
        print(f"    {name:20s} = {val:.2f} MeV  ({err:.2f}% off)")


# ============================================================
# 6. THE MULTI-SITE CHAIN APPROACH
# ============================================================

print("\n" + "=" * 80)
print("STEP 6: FULL CHAIN CORRELATION FUNCTION")
print("=" * 80)

print("""
  Instead of the Fourier-space approach, compute the position-space
  correlator directly. For a long chain of N sites, the 2-point
  correlator ⟨ρ(0)ρ(n)⟩ decays as λⁿ where λ is the largest
  eigenvalue of the transfer matrix in the relevant channel.

  The charge radius is:
    r² = Σ_n n² ⟨ρ(0)ρ(n)⟩ / Σ_n ⟨ρ(0)ρ(n)⟩

  For ⟨ρ(0)ρ(n)⟩ = C·λⁿ:
    r² = a² × Σ n² λⁿ / Σ λⁿ = a² × λ(1+λ)/(1-λ)²

  This is the EXACT result for a geometric correlator on a lattice.
""")

# For the proton, the relevant eigenvalue is the one from the A₂ system
# The paper says the proton corresponds to λ₁^{1/2}, so |λ₁| ≈ 0.4745
lam_proton = 0.4745
lam_glueball = 0.5022

# Geometric sum: r² = a² × λ(1+λ)/(1-λ)²
for name, lam_val in [("glueball λ₀=0.5022", lam_glueball),
                       ("proton λ₁=0.4745", lam_proton)]:
    r_sq_geom = lam_val * (1 + lam_val) / (1 - lam_val)**2
    print(f"\n  {name}:")
    print(f"    r²/a² = λ(1+λ)/(1-λ)² = {r_sq_geom:.6f}")
    print(f"    r/a = {np.sqrt(r_sq_geom):.6f}")

    for a_name, a_val in a_conventions.items():
        r_val = np.sqrt(r_sq_geom) * a_val
        print(f"      with {a_name}: r = {r_val:.6f} fm, r/λ_C = {r_val/lambda_C:.4f}")

# What about λ₁^{1/2}?  The proton is a half-power state.
# If the correlator goes as (λ₁^{1/2})^n = λ₁^{n/2}:
lam_proton_half = np.sqrt(lam_proton)
r_sq_half = lam_proton_half * (1 + lam_proton_half) / (1 - lam_proton_half)**2
print(f"\n  Proton as √λ₁ = {lam_proton_half:.6f}:")
print(f"    r²/a² = {r_sq_half:.6f}")
print(f"    r/a = {np.sqrt(r_sq_half):.6f}")
for a_name, a_val in a_conventions.items():
    r_val = np.sqrt(r_sq_half) * a_val
    print(f"      with {a_name}: r = {r_val:.6f} fm, r/λ_C = {r_val/lambda_C:.4f}")


# ============================================================
# 7. THE DIRECT STRUCTURAL ARGUMENT
# ============================================================

print("\n" + "=" * 80)
print("STEP 7: STRUCTURAL DERIVATION ATTEMPT")
print("=" * 80)

print("""
  The proton mass in the framework is:
    m_p = -ln(λ₁^{1/2}) × m₀/Δ₀ = ½|ln λ₁| × scale

  The Compton wavelength is:
    λ_C = ℏ/(m_p c) = ℏc/m_p = 2ℏcΔ₀/(m₀|ln λ₁|)

  For the charge radius, we need the spatial extent of the charge
  distribution. In the framework, this is set by:

  1. The transfer-matrix correlation length ξ = a/|ln λ₁|
     where a is the lattice spacing.

  2. The charge distribution involves H+1 = 4 real degrees of freedom
     (the real dimension of the spinor S = C²).

  3. The form factor probes the overlap of the H+1 components.
""")

# Analytical formulas
gap_p = 0.5 * abs(np.log(lam_proton))  # effective proton gap (half-power)
m_p_pred = gap_p * scale
print(f"  Proton mass: m_p = ½|ln λ₁| × scale = {m_p_pred:.2f} MeV (PDG: {m_proton:.2f})")

# Correlation length with the "natural" lattice spacing a = ℏc/scale
xi_p = a / (2 * abs(np.log(lam_proton)))  # factor of 2 for λ^{1/2}
print(f"  ξ_p = a/(2|ln λ₁|) = {xi_p:.6f} fm")
print(f"  ξ_p / λ_C = {xi_p / lambda_C:.6f}")

# Check: is r_p / ξ_p = some integer or framework number?
print(f"\n  r_target / ξ_p = {r_p_target / xi_p:.6f}")
print(f"  This should equal a factor from the charge structure.")

# Test: the charge structure factor
# For a dipole form factor with correlation length ξ:
# r² = 12ξ² → r = 2√3 ξ
r_dipole = 2 * np.sqrt(3) * xi_p
print(f"\n  Dipole: r = 2√3 × ξ = {r_dipole:.6f} fm, ratio to target = {r_dipole/r_p_target:.4f}")

# For a single pole: r² = 6ξ² → r = √6 ξ
r_monopole = np.sqrt(6) * xi_p
print(f"  Monopole: r = √6 × ξ = {r_monopole:.6f} fm, ratio = {r_monopole/r_p_target:.4f}")

# What multiplier of ξ gives r_target?
mult = r_p_target / xi_p
print(f"\n  Required: r_target = {mult:.6f} × ξ_p")
# Check against dim_R(S)² = (H+1)² = 16? Or...
print(f"  (H+1)² = {(H+1)**2}")
print(f"  2(H+1) = {2*(H+1)}")
print(f"  √((H+1)³) = {np.sqrt((H+1)**3):.6f}")
print(f"  (H+1)√(H+1) = {(H+1)*np.sqrt(H+1):.6f}")

# Actually compute what r_p / λ_C IS in terms of the eigenvalue
# r_p = (H+1) λ_C = (H+1) ℏc / m_p = (H+1) ℏc / (½|ln λ₁| × m₀/Δ₀)
# = (H+1) × 2Δ₀ / |ln λ₁| × ℏc/m₀
# = (H+1) × 2Δ₀ / |ln λ₁| × a/Δ₀     [since a = ℏc × Δ₀/m₀... NO, a = ℏc/scale = ℏc Δ₀/m₀]
# = (H+1) × 2 × a / |ln λ₁|
# = (H+1) × 2 × ξ_p   (since ξ_p = a / |ln λ₁| ... wait)

# Let me be precise:
# a = ℏc / scale = ℏc Δ₀ / m₀
# m_p = ½|ln λ₁| × scale = ½|ln λ₁| × m₀/Δ₀
# λ_C = ℏc / m_p = ℏc × 2Δ₀ / (m₀ |ln λ₁|) = 2a / |ln λ₁|
# r_p = (H+1) λ_C = 2(H+1) a / |ln λ₁|

r_check = 2 * (H+1) * a / abs(np.log(lam_proton))
print(f"\n  CROSS-CHECK: r_p = 2(H+1) a / |ln λ₁|")
print(f"  = 2 × {H+1} × {a:.6f} / {abs(np.log(lam_proton)):.6f}")
print(f"  = {r_check:.6f} fm")
print(f"  Target = {r_p_target:.6f} fm")
print(f"  Match: {'YES' if abs(r_check - r_p_target) < 1e-4 else 'NO'}")

# So the question reduces to:
# WHY does the form factor give r² = (H+1)² λ_C² = 4(H+1)² a² / |ln λ₁|²?
# In the geometric-correlator formula r²/a² = λ_eff(1+λ_eff)/(1-λ_eff)²,
# what value of λ_eff gives r/a = 2(H+1)/|ln λ₁|?

target_r_over_a = 2*(H+1) / abs(np.log(lam_proton))
print(f"\n  Target r/a = 2(H+1)/|ln λ₁| = {target_r_over_a:.6f}")

# Solve λ(1+λ)/(1-λ)² = (r/a)²
# Let R = r/a. Then λ(1+λ)/(1-λ)² = R²
# λ + λ² = R²(1-2λ+λ²) = R² - 2R²λ + R²λ²
# (R²-1)λ² - (2R²+1)λ + R² = 0
R = target_r_over_a
A_coeff = R**2 - 1
B_coeff = -(2*R**2 + 1)
C_coeff = R**2
disc = B_coeff**2 - 4*A_coeff*C_coeff
if disc >= 0:
    lam_needed_1 = (-B_coeff - np.sqrt(disc)) / (2*A_coeff)
    lam_needed_2 = (-B_coeff + np.sqrt(disc)) / (2*A_coeff)
    print(f"  Geometric correlator needs λ_eff = {lam_needed_1:.8f} or {lam_needed_2:.8f}")
    if 0 < lam_needed_1 < 1:
        print(f"  Physical solution: λ_eff = {lam_needed_1:.8f}")
        print(f"  Compare: λ₁ = {lam_proton}, √λ₁ = {np.sqrt(lam_proton):.8f}")
        print(f"  Compare: λ₁^{1/(H+1)} = {lam_proton**(1/(H+1)):.8f}")


# ============================================================
# 8. DIMENSIONAL ANALYSIS AND THE H+1 FACTOR
# ============================================================

print("\n" + "=" * 80)
print("STEP 8: WHY (H+1) — THE DIMENSIONAL ARGUMENT")
print("=" * 80)

print(f"""
  The proton is a baryon: a bound state of H = 3 quarks.
  In the framework, quarks are sections s_i ∈ S = C² (fundamental spinor).

  The charge distribution of the proton is the probability density
  of the "up-type" charge over the real spatial extent of the bound state.

  KEY OBSERVATION:
  The mass function lives in C^{{H+1}} = C⁴ (H sections + θ substrate).
  The REAL dimension of C^{{H+1}} is 2(H+1) = 8.
  The L₁=1 constraint removes 1 real dimension → 7.
  The Born floor constrains 1 more → 6 effective real dimensions.

  But the CHARGE direction (s₁ - s₂ - s₃) lives in a 1D subspace
  of the 3D section space. The form factor probes the projection
  of the 6D mass-function wavepacket onto the 3D spatial subspace.

  The RMS radius of a d-dimensional isotropic distribution with
  correlation length ξ in each direction is:
    r_RMS = √d × ξ

  Here ξ = λ_C (the Compton wavelength) and d = (H+1)² = 16? No.

  SIMPLER: The factor (H+1) = 4 counts the REAL spinor components
  of S = C². Each component contributes independently to the
  charge radius. The proton samples ALL H+1 real directions of S.

  This gives r_p = (H+1) × (minimum localization scale) = (H+1) × λ_C.
""")

# Numerical verification
print(f"  r_p = (H+1) × ℏ/(m_p c)")
print(f"      = {H+1} × {lambda_C:.6f} fm")
print(f"      = {r_p_target:.6f} fm")
print(f"  CODATA 2018: {r_p_codata} fm")
print(f"  Muonic H:    {r_p_muonic} fm")
print(f"  Error:       {abs(1-r_p_target/r_p_muonic)*100:.3f}% (vs muonic)")


# ============================================================
# 9. THE FORM FACTOR FROM THE SPECTRAL REPRESENTATION
# ============================================================

print("\n" + "=" * 80)
print("STEP 9: SPECTRAL REPRESENTATION OF G_E(q²)")
print("=" * 80)

print("""
  In the spectral representation, the proton electromagnetic form factor is:

    G_E(q²) = ∫ ρ(s)/(s + q²) ds

  where ρ(s) is the spectral function. For a dipole:

    G_E(q²) = (Λ²/(Λ² + q²))²  with  r² = 12/Λ²

  The dipole mass Λ sets the charge radius.

  In the DS framework, the spectral function is determined by the
  eigenvalue spectrum. The proton state λ₁^{1/2} has:
    mass gap = ½|ln λ₁| in lattice units
    Λ_dipole = ? (in principle from the spectral function of the transfer matrix)

  For the identification to work as a derivation, we need:
    Λ = m_p/(H+1) = m_p/4  (giving r = (H+1)λ_C from dipole formula)
""")

# Dipole mass
Lambda_dipole = np.sqrt(12) / r_p_target  # in fm⁻¹
Lambda_dipole_MeV = Lambda_dipole * hbar_c
print(f"  If r_p = {r_p_target:.4f} fm (dipole), then Λ = √12/r_p:")
print(f"    Λ = {Lambda_dipole:.4f} fm⁻¹ = {Lambda_dipole_MeV:.2f} MeV")
print(f"    m_p / Λ = {m_proton / Lambda_dipole_MeV:.6f}")
print(f"    This should be (H+1)/√12 = {(H+1)/np.sqrt(12):.6f} for dipole")

# What form factor gives r² = (H+1)²λ_C² exactly?
# Monopole: G = 1/(1 + q²r₀²) → r² = 6r₀² → need r₀ = (H+1)λ_C/√6
r0_monopole = r_p_target / np.sqrt(6)
Lambda_monopole = 1/r0_monopole * hbar_c  # MeV
print(f"\n  If monopole: G = 1/(1+q²r₀²), r₀ = r_p/√6 = {r0_monopole:.4f} fm")
print(f"    Λ_monopole = ℏc/r₀ = {Lambda_monopole:.2f} MeV")
print(f"    Λ/m_p = {Lambda_monopole/m_proton:.6f}")
print(f"    m_p/Λ = {m_proton/Lambda_monopole:.6f}")

# Gaussian: G = exp(-q²σ²/6) → r² = σ²
sigma_gauss = r_p_target  # fm
Lambda_gauss = hbar_c / sigma_gauss
print(f"\n  If Gaussian: G = exp(-q²r²/6), σ = r_p = {sigma_gauss:.4f} fm")
print(f"    Λ_gauss = ℏc/σ = {Lambda_gauss:.2f} MeV")
print(f"    Λ/m_p = {Lambda_gauss/m_proton:.6f}")
print(f"    m_p/Λ = {m_proton/Lambda_gauss:.6f}")
print(f"    (H+1) = {H+1}")
print(f"    → m_p/Λ_gauss = m_p × r_p / ℏc = m_p × (H+1)λ_C / ℏc = (H+1) ✓")


# ============================================================
# 10. SUMMARY AND HONEST ASSESSMENT
# ============================================================

print("\n" + "=" * 80)
print("SUMMARY AND HONEST ASSESSMENT")
print("=" * 80)

print(f"""
  THE IDENTIFICATION:
    r_p = (H+1) × ℏ/(m_p c) = 4 × {lambda_C:.6f} fm = {r_p_target:.6f} fm
    CODATA 2018: {r_p_codata} fm  →  error {abs(1-r_p_target/r_p_codata)*100:.3f}%
    Muonic H:    {r_p_muonic} fm  →  error {abs(1-r_p_target/r_p_muonic)*100:.3f}%

  WHAT WE ESTABLISHED:
    1. (H+1) = dim_R(S) = dim_R(C²) = 4 is a FRAMEWORK quantity.
    2. The transfer-operator T(k) = J_self + J_cross·e^{{ik}} produces
       dispersing bands, but the lattice-spacing conversion introduces
       model dependence that prevents a clean derivation from the
       transfer operator alone.
    3. The geometric correlator formula r²/a² = λ(1+λ)/(1-λ)² requires
       λ_eff ≈ {lam_needed_1:.4f} to match — this is neither λ₁ nor √λ₁,
       suggesting the simple 1D chain model is too crude.
    4. The identification r_p = (H+1)λ_C is equivalent to saying the
       proton form-factor scale equals m_p/(H+1), i.e., the proton's
       charge is "spread" over dim_R(S) = (H+1) Compton wavelengths.

  THE STRUCTURAL ARGUMENT (strongest):
    The proton is a bound state in the mass space C^{{H+1}}.
    The spinor S = C² has real dimension dim_R(S) = 2·dim_C(S) = 2·2 = 4 = H+1.
    The electromagnetic coupling probes the REAL charge distribution,
    which extends over all dim_R(S) independent real directions of S.
    Each direction contributes one Compton wavelength to the RMS radius:
      r_p = dim_R(S) × λ_C = (H+1) × ℏ/(m_p c)

  STATUS: Numerically excellent identification (0.04%).
    The dim_R(S) argument is structural but not yet a rigorous derivation
    from the form factor. A full derivation requires computing G_E(q²)
    from the multi-site transfer operator on a 3D lattice (not just 1D).
""")
