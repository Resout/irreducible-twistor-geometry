#!/usr/bin/env python3
"""
RIGOROUS VERIFICATION OF THE SQRT(2) MASS RATIO CLAIM

The claim: m(2++)/m(0++) = sqrt(2) exactly, from H=3 geometry.

This script tests every step:
1. Is s2 = s3 exactly at equilibrium (S2 symmetry of m*)?
2. Is the tensor ratio |T2|^2/|T0|^2 = 2 identically for (v1,v2,v2) vectors?
3. Does the mass ratio equal the sqrt of the coupling ratio? (coupling vs log-decay)
4. What are the actual Jacobian eigenvalues and what mass ratio do they give?
5. Is the U(1) photon argument correct (K*=0 for abelian)?
6. Is the SU(2) orbit argument for the graviton correct?

We try to BREAK each step.
"""

import numpy as np
from scipy.optimize import brentq

H = 3
FLOOR = 1.0 / H**3
K_STAR = 7.0 / 30

# ============================================================
# Core dynamics
# ============================================================
def ds_combine(m, e):
    s, th = m[:3], m[3]
    se, ph = e[:3], e[3]
    s_new = s * se + s * ph + th * se
    th_new = th * ph
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)
    d = 1.0 - K
    if abs(d) < 1e-15:
        return m.copy()
    out = np.zeros(4)
    out[:3] = s_new / d
    out[3] = th_new / d
    born = out[3]**2 / np.sum(out**2) if np.sum(out**2) > 0 else 1.0
    if born < FLOOR:
        abs_s = np.abs(out[:3])
        ss = np.sum(abs_s)
        sq = np.sum(abs_s**2)
        if ss > 1e-15:
            r = sq / ss**2
            a_c = 26.0 - r
            b_c = 2.0 * r
            c_c = -r
            disc = b_c**2 - 4*a_c*c_c
            tn = (-b_c + np.sqrt(disc)) / (2*a_c)
            sc = (1.0 - tn) / ss
            out[:3] = abs_s * sc
            out[3] = tn
    return out

def find_eq():
    def K_at(p_dom):
        p_w = (1.0 - p_dom) / 2.0
        sc = 1.0 - FLOOR
        raw = np.array([np.sqrt(p_dom*sc), np.sqrt(p_w*sc),
                        np.sqrt(p_w*sc), np.sqrt(FLOOR)])
        e = raw / np.sum(raw)
        m = np.array([0.4, 0.2, 0.2, 0.2])
        for _ in range(10000):
            m2 = ds_combine(m, e)
            if np.max(np.abs(m2 - m)) < 1e-15:
                break
            m = m2
        s = m[:3]; se = e[:3]
        return sum(s[i]*se[j] for i in range(3) for j in range(3) if i != j)

    p_dom = brentq(lambda p: K_at(p) - K_STAR, 0.90, 0.96, xtol=1e-14)
    p_w = (1.0 - p_dom) / 2.0
    sc = 1.0 - FLOOR
    raw = np.array([np.sqrt(p_dom*sc), np.sqrt(p_w*sc),
                    np.sqrt(p_w*sc), np.sqrt(FLOOR)])
    e_star = raw / np.sum(raw)
    m = np.array([0.4, 0.2, 0.2, 0.2])
    for _ in range(10000):
        m2 = ds_combine(m, e_star)
        if np.max(np.abs(m2 - m)) < 1e-15:
            break
        m = m2
    return m, e_star

def jacobian(m_star, e_star, eps=1e-8):
    J = np.zeros((4, 4))
    for j in range(4):
        mp = m_star.copy(); mp[j] += eps
        mm = m_star.copy(); mm[j] -= eps
        J[:, j] = (ds_combine(mp, e_star) - ds_combine(mm, e_star)) / (2*eps)
    return J

# ============================================================
# STEP 1: Is s2 = s3 exactly?
# ============================================================
print("=" * 70)
print("STEP 1: S2 SYMMETRY OF EQUILIBRIUM")
print("=" * 70)

m_star, e_star = find_eq()
print(f"m* = {m_star}")
print(f"e* = {e_star}")
print(f"s1 = {m_star[0]:.15f}")
print(f"s2 = {m_star[1]:.15f}")
print(f"s3 = {m_star[2]:.15f}")
print(f"s2 - s3 = {m_star[1] - m_star[2]:.3e}")
print(f"e2 - e3 = {e_star[1] - e_star[2]:.3e}")

# Is this by construction (the evidence was symmetric in 2,3)?
print(f"\nThe evidence was constructed with p_w = (1-p_dom)/2 for BOTH components 2 and 3.")
print(f"So e2=e3 by construction. DS combination preserves s2=s3 if e2=e3 and initial s2=s3.")
print(f"Starting from s2=s3=0.2, with e2=e3: YES, s2=s3 is maintained by symmetry.")
print(f"This is exact: not numerical coincidence, but a consequence of the symmetric initialisation.")
print(f"\nBUT: Is this a PHYSICAL requirement or just a consequence of how we search for the equilibrium?")
print(f"The equilibrium WITH s2≠s3 would also have K*=7/30 — it's on the moduli curve.")
print(f"The symmetric equilibrium (s2=s3) is the one the symmetric initial condition converges to.")
print(f"There is a 1D moduli space at K*=7/30. The symmetric point is ONE point on it.")

# ============================================================
# STEP 2: Is |T2|^2/|T0|^2 = 2 for (v1,v2,v2) vectors?
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: TENSOR RATIO IDENTITY")
print("=" * 70)

print("""
The claim: for δm = (v1, v2, v2, 0) [no θ component],
the ratio |T2|^2/|T0|^2 = 2 identically.

T0 = trace of the Pauli matrix: (1/sqrt(2)) * sum_i delta_s_i * I
   = 0 (since Tr(sigma_i) = 0, so trace part is the theta component)

Wait -- let me think about this carefully.

Under the SO(4) decomposition of dbar(Phi):
  (1,1) sector = trace part = pure delta_theta
  (3,3) sector = symmetric traceless = anisotropic delta_s
  (3,1)+(1,3) = antisymmetric = gauge

The Pauli embedding: M = (theta*I + s1*sigma1 + s2*sigma2 + s3*sigma3)/sqrt(2)

A perturbation delta_m = (delta_s1, delta_s2, delta_s3, delta_theta):
  delta_M = (delta_theta * I + delta_s . sigma) / sqrt(2)

T0 = scalar sector = (1/4) tr(delta_M) * I * sqrt(2) = delta_theta * I / sqrt(2)
T2 = symmetric traceless = (delta_M + delta_M^T)/2 - trace part

For delta_m = (v1, v2, v2, 0):
  delta_M = (v1*sigma1 + v2*sigma2 + v2*sigma3) / sqrt(2)

sigma1 = [[0,1],[1,0]], sigma2 = [[0,-i],[i,0]], sigma3 = [[1,0],[0,-1]]

delta_M = [[v2/sqrt(2) * (-i+1) ... ]]

Hmm, the (v1,v2,v2) vector refers to which section components?
The claim says it's the FLOOR KICK that has this form.
Let me compute the floor kick directly.
""")

# Compute what the floor does to a small perturbation of m*
eps_test = 1e-6
J_full = jacobian(m_star, e_star)

print("Full Jacobian eigenvalues:")
evals, evecs = np.linalg.eig(J_full)
for i, (ev, vec) in enumerate(sorted(zip(evals, evecs.T), key=lambda x: -abs(x[0]))):
    print(f"  lambda_{i} = {ev:.10f}  |lambda| = {abs(ev):.10f}")
    print(f"  eigvec = [{', '.join(f'{v:.6f}' for v in vec.real)}]")

# The two physical eigenvalues
lam0 = sorted(evals, key=lambda x: -abs(x))[0]
lam1 = sorted(evals, key=lambda x: -abs(x))[1]
print(f"\nlambda_0 = {lam0:.10f}")
print(f"lambda_1 = {lam1:.10f}")
print(f"-ln(lambda_0) = {-np.log(abs(lam0)):.10f}")
print(f"-ln(lambda_1) = {-np.log(abs(lam1)):.10f}")
print(f"ratio of Deltas = {-np.log(abs(lam1)) / (-np.log(abs(lam0))):.10f}")
print(f"sqrt(2) = {np.sqrt(2):.10f}")
print(f"ratio / sqrt(2) = {(-np.log(abs(lam1)) / (-np.log(abs(lam0)))) / np.sqrt(2):.10f}")

# ============================================================
# STEP 3: The crucial step -- coupling ratio vs mass ratio
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: COUPLING RATIO vs MASS RATIO")
print("=" * 70)

print("""
The claim says |T2|^2/|T0|^2 = 2, hence m(2++)/m(0++) = sqrt(2).

This step requires: mass ratio = sqrt(coupling ratio).

Is this true? Mass = -ln(lambda). If the couplings (eigenvalues) satisfy
lambda_1^2 = lambda_0... then ln(lambda_1) = ln(lambda_0)/2, so
Delta_1 = Delta_0/2. That would give ratio 1/2, not sqrt(2).

If instead lambda_1 = lambda_0^(1/sqrt(2))... that's not natural.

The ACTUAL question is: what determines the eigenvalue ratio?

Let me compute: lambda_1 / lambda_0 and (lambda_1/lambda_0)^2
""")

ratio_lam = abs(lam1) / abs(lam0)
print(f"lambda_1 / lambda_0 = {ratio_lam:.10f}")
print(f"(lambda_1/lambda_0)^2 = {ratio_lam**2:.10f}")
print(f"sqrt(lambda_1/lambda_0) = {np.sqrt(ratio_lam):.10f}")
print(f"Delta_1/Delta_0 = {(-np.log(abs(lam1)))/(-np.log(abs(lam0))):.10f}")
print(f"sqrt(2) = {np.sqrt(2):.10f}")
print(f"Difference from sqrt(2): {abs((-np.log(abs(lam1)))/(-np.log(abs(lam0))) - np.sqrt(2)):.4e}")

# ============================================================
# STEP 4: The SO(4) structure of eigenvectors
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: EIGENVECTOR STRUCTURE (what modes do lambda_0 and lambda_1 correspond to?)")
print("=" * 70)

evals_sorted = sorted(zip(np.abs(evals), evals, evecs.T), key=lambda x: -x[0])
for i, (mag, ev, vec) in enumerate(evals_sorted[:4]):
    s_content = np.sum(np.abs(vec[:3])**2) / (np.sum(np.abs(vec)**2) + 1e-30)
    th_content = np.abs(vec[3])**2 / (np.sum(np.abs(vec)**2) + 1e-30)
    s = vec[:3].real
    if np.sum(s**2) > 1e-20:
        s_unit = s / np.sqrt(np.sum(s**2))
        iso = np.dot(s_unit, np.array([1,1,1])/np.sqrt(3))**2
    else:
        iso = 0
    print(f"mode {i}: |lambda|={mag:.8f}  s%={s_content*100:.1f}  theta%={th_content*100:.1f}  "
          f"isotropy={iso:.4f}  vec={vec.real}")

# ============================================================
# STEP 5: What IS the scalar glueball eigenvalue?
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: SCALAR GLUEBALL -- ISOTROPIC (0++) MODE")
print("=" * 70)

# Project J onto the isotropic direction (1,1,1,0)/sqrt(3)
iso_dir = np.array([1.0, 1.0, 1.0, 0.0]) / np.sqrt(3)
aniso_dir = np.array([1.0, -0.5, -0.5, 0.0]) / np.sqrt(1.5)  # traceless s1-dominant
aniso2_dir = np.array([0.0, 1.0, -1.0, 0.0]) / np.sqrt(2)   # s2-s3

J_iso = J_full @ iso_dir
J_aniso = J_full @ aniso_dir
J_aniso2 = J_full @ aniso2_dir

print(f"Isotropic direction (1,1,1,0)/sqrt(3):")
print(f"  J @ iso_dir = {J_iso}")
print(f"  Rayleigh quotient = {np.dot(J_iso, iso_dir):.8f}")

print(f"\nAnisotropic direction (1,-1/2,-1/2,0) [s1-dominant]:")
print(f"  J @ aniso = {J_aniso}")
print(f"  Rayleigh quotient = {np.dot(J_aniso, aniso_dir):.8f}")

print(f"\nAnisotropic direction (0,1,-1,0)/sqrt(2) [s2-s3]:")
print(f"  J @ aniso2 = {J_aniso2}")
print(f"  Rayleigh quotient = {np.dot(J_aniso2, aniso2_dir):.8f}")

# ============================================================
# STEP 6: SU(2) orbit -- what do rotations of s do?
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: SU(2) ORBIT -- does M -> UMU† preserve Born and K?")
print("=" * 70)

# Apply a rotation in the s2-s3 plane (s1 fixed, s2->s2*cos+s3*sin, s3->-s2*sin+s3*cos)
def rotate_m(m, angle):
    """Rotate s2,s3 components by angle (leave s1,theta fixed)."""
    m_rot = m.copy()
    c, s_a = np.cos(angle), np.sin(angle)
    m_rot[1] = c * m[1] - s_a * m[2]
    m_rot[2] = s_a * m[1] + c * m[2]
    return m_rot

print("Effect of rotating s on S^2 (preserving |s|, theta):")
print(f"Original m* = {m_star}")

for angle in [0.0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, np.pi]:
    m_rot = rotate_m(m_star, angle)
    born_rot = m_rot[3]**2 / np.sum(m_rot**2)
    K_rot = sum(m_rot[i]*e_star[j] for i in range(3) for j in range(3) if i != j)
    print(f"  angle={angle/np.pi:.3f}pi: Born={born_rot:.8f}  K={K_rot:.8f}  "
          f"K-K*={K_rot-K_STAR:.2e}  Born-1/27={born_rot-1/27:.2e}")

print(f"\nKey question: does K change when we rotate m but keep e* fixed?")
print(f"If yes: the 'graviton mode' DOES change K -> it's not flat -> it gets a mass.")
print(f"If no: the rotation is truly a flat direction -> massless.")

# ============================================================
# STEP 7: What happens when we APPLY THE DYNAMICS to a rotated m?
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: DYNAMICS ON THE SU(2) ORBIT")
print("=" * 70)

print("Apply Phi = DS+floor to rotated m* with FIXED e*:")
print("If Phi(m_rot, e*) = m_rot (rotated equilibrium), the orbit is flat.")
print("If Phi(m_rot, e*) != m_rot, the orbit gets pulled back to m*.")

for angle in [np.pi/8, np.pi/4, np.pi/2]:
    m_rot = rotate_m(m_star, angle)
    m_out = ds_combine(m_rot, e_star)
    diff = m_out - m_rot
    diff_from_mstar = m_out - m_star
    print(f"\n  angle={angle/np.pi:.3f}pi:")
    print(f"    m_rot = {m_rot}")
    print(f"    Phi(m_rot, e*) = {m_out}")
    print(f"    Phi - m_rot = {diff}")
    print(f"    ||Phi - m_rot|| = {np.linalg.norm(diff):.6e}")
    print(f"    ||Phi - m*|| = {np.linalg.norm(diff_from_mstar):.6e}")

print(f"""
CRITICAL: The SU(2) orbit argument requires that when m is rotated
but e is NOT rotated (e stays at e*), Phi maps back to m*.
This would mean the rotated state is NOT an equilibrium -- it gets
pulled back. That means the rotational mode IS massive.

For the graviton to be massless, we need either:
(a) e rotates WITH m (self-evidence: e=m), or
(b) The Jacobian eigenvalue for the rotational mode is exactly 1.

The paper's argument is about BASE VARIATIONS: smooth changes of
m*(x) across S^4 that preserve K*=7/30 everywhere. This requires
BOTH m and e to change coherently at every point. A single-site
rotation with fixed e is NOT the right test.
""")

# ============================================================
# STEP 8: The actual mass ratio from the Jacobian
# ============================================================
print("\n" + "=" * 70)
print("STEP 8: ACTUAL DELTA RATIO FROM JACOBIAN")
print("=" * 70)

d0 = -np.log(abs(lam0))
d1 = -np.log(abs(lam1))
print(f"Delta_0 = {d0:.10f}")
print(f"Delta_1 = {d1:.10f}")
print(f"Delta_1 / Delta_0 = {d1/d0:.10f}")
print(f"sqrt(2) = {np.sqrt(2):.10f}")
print(f"Deviation: {abs(d1/d0 - np.sqrt(2)):.4e}")
print(f"Deviation in %: {abs(d1/d0 - np.sqrt(2))/np.sqrt(2)*100:.4f}%")

print(f"""
The ratio Delta_1/Delta_0 = {d1/d0:.8f}
This is {'close to' if abs(d1/d0 - np.sqrt(2)) < 0.01 else 'NOT close to'} sqrt(2) = {np.sqrt(2):.8f}

The deviation is {abs(d1/d0 - np.sqrt(2)):.4e}.

This is the single-site SU(2) (rank 1) result.
The claim is for SU(3) (rank 2, A2 Cartan matrix).
Let me check SU(3).
""")

# ============================================================
# STEP 9: SU(3) computation
# ============================================================
print("\n" + "=" * 70)
print("STEP 9: SU(3) COUPLED JACOBIAN EIGENVALUES")
print("=" * 70)

cartan_su3 = np.array([[2, -1], [-1, 2]], dtype=float)

def coupled_step_su3(states, cartan, g, e_base):
    r = cartan.shape[0]
    new = []
    for i in range(r):
        e_i = e_base.copy()
        for j in range(r):
            if i != j and abs(cartan[i,j]) > 0.01:
                e_i = e_i + g * cartan[i,j] * (states[j] - 0.25)
        e_i = np.abs(e_i); e_i = np.maximum(e_i, 1e-10); e_i /= np.sum(e_i)
        new.append(ds_combine(states[i], e_i))
    return new

def find_coupled_eq_su3(n_iter=20000):
    states = [m_star.copy(), m_star.copy()]
    for _ in range(n_iter):
        new = coupled_step_su3(states, cartan_su3, K_STAR, e_star)
        diff = max(np.max(np.abs(new[i] - states[i])) for i in range(2))
        states = new
        if diff < 1e-15:
            break
    return states

states_su3 = find_coupled_eq_su3()
print(f"SU(3) equilibrium:")
for i, s in enumerate(states_su3):
    print(f"  site {i}: {s}")

# Full Jacobian of SU(3)
x0 = np.concatenate(states_su3)
def F_su3(x):
    ss = [x[:4], x[4:]]
    return np.concatenate(coupled_step_su3(ss, cartan_su3, K_STAR, e_star))

eps = 1e-7
J_su3 = np.zeros((8, 8))
for j in range(8):
    xp = x0.copy(); xp[j] += eps
    xm = x0.copy(); xm[j] -= eps
    J_su3[:, j] = (F_su3(xp) - F_su3(xm)) / (2*eps)

evals_su3 = np.linalg.eigvals(J_su3)
physical_su3 = sorted([e for e in evals_su3 if 1e-8 < abs(e) < 1-1e-8],
                      key=lambda x: -abs(x))

print(f"\nSU(3) physical eigenvalues:")
for i, ev in enumerate(physical_su3):
    delta = -np.log(abs(ev))
    print(f"  lambda_{i} = {ev:.8f}  |lambda| = {abs(ev):.8f}  Delta = {delta:.6f}")

if len(physical_su3) >= 2:
    d0_su3 = -np.log(abs(physical_su3[0]))
    d1_su3 = -np.log(abs(physical_su3[1]))
    print(f"\nDelta_1/Delta_0 = {d1_su3/d0_su3:.10f}")
    print(f"sqrt(2)         = {np.sqrt(2):.10f}")
    print(f"Deviation       = {abs(d1_su3/d0_su3 - np.sqrt(2)):.4e}")
    print(f"Deviation in %  = {abs(d1_su3/d0_su3 - np.sqrt(2))/np.sqrt(2)*100:.4f}%")

    # Extended table
    print(f"\nFull mass ratio table (SU(3)):")
    ref = d0_su3
    for i, ev in enumerate(physical_su3):
        d = -np.log(abs(ev))
        print(f"  mode {i}: m/m0 = {d/ref:.6f}")

# ============================================================
# STEP 10: Is the T0/T2 argument actually about Jacobian eigenvectors?
# ============================================================
print("\n" + "=" * 70)
print("STEP 10: WHAT IS T0 AND T2 EXACTLY?")
print("=" * 70)

print("""
The claim says: for the floor kick vector delta_m = (v1, v2, v2, 0),
|T2|^2/|T0|^2 = 2.

Let me understand what T0 and T2 are in this context.

Under SO(4) decomposition of a 4-vector delta_m = (ds1, ds2, ds3, dtheta):
  Embedding into M_2(C): delta_M = (dtheta*I + ds1*sigma1 + ds2*sigma2 + ds3*sigma3)/sqrt(2)

  T0 = (1,1) part = (1/2) tr(delta_M) * I = dtheta * I / sqrt(2)
  T2 = (3,3) part = symmetric traceless = delta_M - T0 - antisymmetric

For delta_m = (v1, v2, v2, 0) [no theta component]:
  T0 = 0 (no theta)
  T2 = (v1*sigma1 + v2*sigma2 + v2*sigma3)/sqrt(2) -- but wait, this is ALL of delta_M

Hmm, that gives |T2|/|T0| = infinity, not sqrt(2).

So T0 and T2 must refer to something else in the claim.
Let me re-read: "the ratio |T2|^2/|T0|^2 = 2 identically".

Maybe T0 and T2 are the projections onto the TWO Jacobian eigenmodes?
Or the projections onto J=0 and J=2 spherical harmonics on S^2?
""")

# Project the two eigenvectors onto isotropic vs anisotropic directions
evals_s, evecs_s = np.linalg.eig(J_full)
phys = [(abs(ev), ev, vec) for ev, vec in zip(evals_s, evecs_s.T) if 1e-8 < abs(ev) < 1-1e-8]
phys.sort(key=lambda x: -x[0])

print("The two physical eigenvectors of J:")
for i, (mag, ev, vec) in enumerate(phys[:2]):
    v = vec.real
    print(f"\n  mode {i} (lambda={ev.real:.6f}):")
    print(f"  raw = {v}")
    # Decompose into isotropic (s1+s2+s3)/sqrt(3) and anisotropic
    iso = np.array([1,1,1,0])/np.sqrt(3)
    ani1 = np.array([2,-1,-1,0])/np.sqrt(6)
    ani2 = np.array([0,1,-1,0])/np.sqrt(2)
    thet = np.array([0,0,0,1])

    proj_iso = np.dot(v, iso)
    proj_ani1 = np.dot(v, ani1)
    proj_ani2 = np.dot(v, ani2)
    proj_thet = np.dot(v, thet)
    print(f"  iso (0++) component: {proj_iso:.6f}")
    print(f"  ani1 (2++ sym) component: {proj_ani1:.6f}")
    print(f"  ani2 (2++ asym) component: {proj_ani2:.6f}")
    print(f"  theta component: {proj_thet:.6f}")
    print(f"  |ani|^2/|iso|^2 = {(proj_ani1**2+proj_ani2**2)/(proj_iso**2+1e-30):.6f}")
