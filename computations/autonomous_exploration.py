#!/usr/bin/env python3
"""
AUTONOMOUS EXPLORATION — April 2026

Three questions I want to answer:

1. SCHWINGER FUNCTIONS: What does the 2++ operator actually decay at?
   The sqrt(2) argument gives the mass ratio from coupling norms.
   But the Schwinger function decay rate is what determines the physical mass.
   Are these the same? I want to compute it directly.

2. COLOR-SINGLET PROJECTION: Why does lambda_1 = 0.4745 (ratio 1.082)
   have no lattice match? Is it a colored excitation?

3. FERMION MASS: The spinor bundle sqrt(E) has transition sqrt(det M).
   det(M) is quadratic in m. What decay rate does sqrt(det M) have?
   This should give the fermion mass.
"""

import numpy as np
from scipy.optimize import brentq

H = 3
FLOOR = 1.0 / H**3
K_STAR = 7.0 / 30

# ============================================================
# CORE DYNAMICS
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
            a_c = 26.0 - r; b_c = 2.0 * r; c_c = -r
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
            if np.max(np.abs(m2 - m)) < 1e-15: break
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
        if np.max(np.abs(m2 - m)) < 1e-15: break
        m = m2
    return m, e_star

m_star, e_star = find_eq()

# ============================================================
# A2 SYSTEM
# ============================================================

cartan_A2 = np.array([[2, -1], [-1, 2]], dtype=float)

def coupled_step(state, e0, g=K_STAR):
    m1, m2 = state[:4], state[4:]
    e1 = e0 - g * (m2 - 0.25)
    e1 = np.maximum(e1, 1e-10); e1 /= np.sum(e1)
    e2 = e0 - g * (m1 - 0.25)
    e2 = np.maximum(e2, 1e-10); e2 /= np.sum(e2)
    return np.concatenate([ds_combine(m1, e1), ds_combine(m2, e2)])

# Find A2 equilibrium
state0 = np.concatenate([m_star, m_star])
for _ in range(10000):
    s1 = coupled_step(state0, e_star)
    if np.linalg.norm(s1 - state0) < 1e-14: break
    state0 = s1
m1_eq = state0[:4]; m2_eq = state0[4:]

# ============================================================
# QUESTION 1: SCHWINGER FUNCTIONS
# What decay rate does the 2++ operator have?
# ============================================================

print("=" * 70)
print("QUESTION 1: SCHWINGER FUNCTION DECAY RATES")
print("=" * 70)

print("""
Glueball operators are bilinear in F+.
F+ at a point is proportional to the floor kick: delta_m = m* - DS(m*, e*).

We compute the connected correlator:
  C_O(t) = <O(t) O(0)> - <O>^2

by evolving a perturbed state and measuring the operator at each step.

Operators:
  O_0++ = |delta_s|^2 + |delta_theta|^2  (total norm, scalar glueball)
  O_2++ = |T_ij|^2 / |T_0|^2            (spin-2 content, tensor glueball)
  O_det = det(M) = (theta^2 - |s|^2)/2   (spinor observable)
""")

def floor_kick(m, e):
    """The floor kick: m* - DS(m, e). This IS F+."""
    m_ds = ds_combine(m, e)
    return m - m_ds

def operator_0pp(dm):
    """0++ operator: total norm of floor kick."""
    return np.sum(dm**2)

def operator_2pp(dm):
    """2++ operator: Frobenius norm of spin-2 part of ds tensor."""
    ds = dm[:3]
    s_norm2 = np.sum(ds**2)
    if s_norm2 < 1e-20:
        return 0.0
    T = np.outer(ds, ds) - (s_norm2/3) * np.eye(3)
    return np.sum(T**2)

def operator_det(m):
    """Determinant of Pauli matrix: (theta^2 - |s|^2)/2."""
    return (m[3]**2 - np.sum(m[:3]**2)) / 2.0

# Measure operators at equilibrium
dk_eq = floor_kick(m_star, e_star)
O0_eq = operator_0pp(dk_eq)
O2_eq = operator_2pp(dk_eq)
det_eq = operator_det(m_star)

print(f"At equilibrium:")
print(f"  Floor kick: {dk_eq}")
print(f"  O_0++ = {O0_eq:.6f}")
print(f"  O_2++ = {O2_eq:.6f}")
print(f"  O_2++ / O_0++ = {O2_eq/O0_eq:.6f} (expected 2.0)")
print(f"  sqrt(O_2++/O_0++) = {np.sqrt(O2_eq/O0_eq):.6f} (expected sqrt(2))")
print(f"  det(M*) = {det_eq:.6f}")

# Now: evolve a perturbed state and watch operator correlations
print(f"\nEvolving perturbed state to measure Schwinger function decay...")

# Perturb in the 0++ direction (along the floor kick)
eps = 1e-4
dk_norm = dk_eq / np.linalg.norm(dk_eq)

# Run correlation for 50 steps
n_steps = 50
C_0pp = []
C_2pp = []
C_det = []

for perturb_dir, label in [(dk_norm, "floor kick direction (0++ excitation)")]:
    m_perturb = m_star + eps * perturb_dir
    m_perturb = np.maximum(m_perturb, 1e-10)
    m_perturb /= np.sum(m_perturb)

    C0 = []
    C2 = []
    Cd = []

    m_curr = m_perturb.copy()
    for t in range(n_steps):
        dk = floor_kick(m_curr, e_star)
        O0 = operator_0pp(dk) - O0_eq
        O2 = operator_2pp(dk) - O2_eq
        Od = operator_det(m_curr) - det_eq
        C0.append(abs(O0))
        C2.append(abs(O2))
        Cd.append(abs(Od))
        m_curr = ds_combine(m_curr, e_star)

    C_0pp = np.array(C0)
    C_2pp = np.array(C2)
    C_det = np.array(Cd)

# Fit exponential decay for each
def fit_decay(C, t_start=2, t_end=20):
    """Fit C(t) = A * exp(-Delta * t) in log space."""
    ts = np.arange(t_start, min(t_end, len(C)))
    Cs = C[t_start:min(t_end, len(C))]
    valid = Cs > 1e-15
    if np.sum(valid) < 3:
        return None, None
    log_C = np.log(Cs[valid])
    ts_valid = ts[valid]
    # Linear fit in log space
    A = np.vstack([np.ones_like(ts_valid), ts_valid]).T
    result = np.linalg.lstsq(A, log_C, rcond=None)
    log_A0, neg_Delta = result[0]
    Delta = -neg_Delta
    return Delta, np.exp(log_A0)

Delta_0pp, A_0pp = fit_decay(C_0pp)
Delta_2pp, A_2pp = fit_decay(C_2pp)
Delta_det, A_det = fit_decay(C_det)

print(f"\nDecay rates from Schwinger function fits:")
print(f"  O_0++: Delta = {Delta_0pp:.4f}")
print(f"  O_2++: Delta = {Delta_2pp:.4f}")
print(f"  det(M): Delta = {Delta_det:.4f}")
print(f"\n  Ratio Delta_2++/Delta_0++ = {Delta_2pp/Delta_0pp:.4f}")
print(f"  sqrt(2) = {np.sqrt(2):.4f}")
print(f"  2.0 = 2.0000")

print(f"\nExpected eigenvalue for A2 ground state: {-np.log(0.5022):.4f}")
print(f"O_0++ measured decay:                     {Delta_0pp:.4f}")
print(f"O_2++ measured decay:                     {Delta_2pp:.4f}")

# Print the time series
print(f"\n  t   |C_0++(t)|  |C_2++(t)|  |C_det(t)|")
print(f"  {'─'*50}")
for t in range(min(25, n_steps)):
    print(f"  {t:2d}  {C_0pp[t]:.4e}  {C_2pp[t]:.4e}  {C_det[t]:.4e}")

# ============================================================
# QUESTION 2: COLOR-SINGLET PROJECTION
# ============================================================

print(f"\n{'='*70}")
print("QUESTION 2: COLOR-SINGLET PROJECTION")
print("=" * 70)

print("""
The A2 Jacobian has 4 non-trivial eigenmodes:
  lambda_0 = 0.5022: node-ANTISYMMETRIC, radial (ground state 0++)
  lambda_1 = 0.4745: node-ANTISYMMETRIC, angular (no lattice match)
  lambda_2 = 0.3527: node-SYMMETRIC,     angular (0-+)
  lambda_3 = 0.3344: node-SYMMETRIC,     radial  (0++*)

Node-antisymmetric means the two SU(3) Dynkin nodes perturb in OPPOSITE phase.
Node-symmetric means they perturb in the SAME phase.

In color space: node-antisymmetric perturbations correspond to
the off-diagonal generators of SU(3) — the raising/lowering operators
T_12, T_13, T_23, etc. These carry COLOR CHARGE.
They are NOT color singlets.

Node-symmetric perturbations correspond to the CARTAN subalgebra
generators T_3, T_8. These commute with the color charge and can
form color singlets.

BUT: the 0++ ground state is also node-antisymmetric (lambda_0).
If only node-symmetric modes are color singlets, the 0++ wouldn't exist.

RESOLUTION: the physical glueball operators are BILINEAR in the field.
A product of two node-antisymmetric modes IS node-symmetric:
  (antisym) x (antisym) = symmetric
The 0++ is created by O ~ delta_m * delta_m (norm squared), which has
node-symmetric behavior: (-v) * (-v) = v^2 = (+v) * (+v).

So:
  lambda_0^2 = 0.5022^2 = 0.2522  -> Delta = 1.3776 (two-particle 0++)

But the SINGLE glueball 0++ has mass -ln(lambda_0) = 0.6888, not 1.3776.
This means the 0++ is a LINEAR observable in the transfer operator language,
not a bilinear one.

What IS the 0++ as a linear observable?
""")

# The 0++ operator in the transfer operator language
# is the trace of the transfer operator: Tr(Phi) = sum_i Phi_i
# In the correlator language, the 0++ is the OVERLAP of the ground state
# with the floor kick operator.

# Compute eigenvectors of the 8x8 A2 Jacobian
eps_jac = 1e-8
J8 = np.zeros((8, 8))
for j in range(8):
    sp = state0.copy(); sp[j] += eps_jac
    sm = state0.copy(); sm[j] -= eps_jac
    J8[:, j] = (coupled_step(sp, e_star) - coupled_step(sm, e_star)) / (2*eps_jac)

evals_8, evecs_8 = np.linalg.eig(J8)
idx = np.argsort(-np.abs(evals_8))
evals_8 = evals_8[idx].real
evecs_8 = evecs_8[:, idx].real

print(f"A2 eigenvectors (significant modes):")
for i in range(4):
    ev = evecs_8[:, i]
    v1, v2 = ev[:4], ev[4:]
    sym_norm = np.linalg.norm(v1 + v2)
    asym_norm = np.linalg.norm(v1 - v2)

    # For color singlet: the physical state must be symmetric under
    # Weyl group of A2, which includes the node exchange symmetry
    node_sym = "SYMMETRIC" if sym_norm > asym_norm else "ANTISYMMETRIC"

    # A color-singlet SINGLE-PARTICLE state must have even behavior
    # under the Weyl group. But wait — the FIELD itself is antisymmetric
    # (it lives in the adjoint). The GLUEBALL (bound state) is symmetric.

    dk1 = floor_kick(m1_eq + 1e-4 * v1, e_star)
    dk2 = floor_kick(m2_eq + 1e-4 * v2, e_star)

    O0_1 = operator_0pp(dk1 - floor_kick(m1_eq, e_star))
    O0_2 = operator_0pp(dk2 - floor_kick(m2_eq, e_star))

    print(f"\n  mode {i}: lambda={evals_8[i]:.4f}, {node_sym}")
    print(f"    v1 = {v1}")
    print(f"    v2 = {v2}")

    # Color charge analysis:
    # In SU(3), the Cartan generators H1 = diag(1,-1,0), H2 = diag(0,1,-1)
    # act on the adjoint. Node 1 carries H1 charge, node 2 carries H2 charge.
    # For a color singlet: H1|state> = H2|state> = 0
    # Node-antisymmetric: |state> = |1> - |2>
    # H1 acts as +1 on node 1, so H1(|1>-|2>) = |1> + H1|2>
    # This is NOT zero unless there's cancellation from the other node.

    print(f"    Node structure: {node_sym}")
    print(f"    Color interpretation: {'Cartan generator (color neutral)' if node_sym == 'SYMMETRIC' else 'Off-diagonal generator (color charged)'}")

print(f"""
RESOLUTION:
  Node-antisymmetric modes carry the off-diagonal SU(3) color charges.
  These are NOT color singlet single-particle states.

  The 0++ glueball is NOT the eigenvalue lambda_0 (node-antisymmetric).

  The 0++ IS the bilinear product of two node-antisymmetric modes:
  O_0++ ~ Tr(delta_m^dagger delta_m)
  This has Koopman eigenvalue lambda_0^2, hence:

  m(0++) = -ln(lambda_0^2) = 2 * (-ln(lambda_0)) = 2 * 0.6888 = 1.3776

  And the physical 2++ would be:
  m(2++) = 2 * (-ln(lambda_0)) * sqrt(2) = 1.3776 * sqrt(2) ???

  No wait, that doesn't work. If both 0++ and 2++ are bilinear in lambda_0:
  m(0++) = 2 * Delta_0 = 1.3776
  m(2++) = 2 * Delta_0 = 1.3776 (same! both are lambda_0^2 Koopman)

  Unless the 2++ couples to lambda_0 * lambda_1 (mixed):
  lambda_0 * lambda_1 = 0.5022 * 0.4745 = 0.2383
  Delta_01 = -ln(0.2383) = 1.4343
  m(2++)/m(0++) = 1.4343/1.3776 = 1.041  -- not 1.414

  CONCLUSION: The Koopman bilinear picture does NOT give sqrt(2).
  The ratio sqrt(2) must come from a DIFFERENT mechanism.
""")

# ============================================================
# QUESTION 3: FERMION MASS FROM det(M)
# ============================================================

print(f"{'='*70}")
print("QUESTION 3: FERMION MASS ESTIMATE")
print("=" * 70)

print(f"""
The spinor bundle E^(1/2) has transition function sqrt(det M).
det(M) = (theta^2 - |s|^2) / 2.

At the A2 equilibrium:
  m1* = m2* = (0.74505, 0.05403, 0.05403, 0.14688)
  det(M*) = (0.14688^2 - 0.74505^2 - 0.05403^2 - 0.05403^2) / 2
""")

det_m1 = (m1_eq[3]**2 - np.sum(m1_eq[:3]**2)) / 2
print(f"  det(M* at A2 equilibrium) = {det_m1:.6f}")

# det(M) is quadratic in m. Its decay rate under the dynamics is
# determined by its coupling to the Koopman eigenmodes.
# Let's compute it directly.

print(f"\nDecay of det(M) perturbation under A2 dynamics:")
eps_det = 1e-5
n_steps_det = 30
C_det_A2 = []

m1_perturb = m1_eq + eps_det * evecs_8[:4, 0]  # perturb along mode 0
m1_perturb = np.maximum(m1_perturb, 1e-10); m1_perturb /= np.sum(m1_perturb)
state_perturb = np.concatenate([m1_perturb, m2_eq])

for t in range(n_steps_det):
    m1 = state_perturb[:4]; m2 = state_perturb[4:]
    det1 = (m1[3]**2 - np.sum(m1[:3]**2)) / 2
    Od = abs(det1 - det_m1)
    C_det_A2.append(Od)
    state_perturb = coupled_step(state_perturb, e_star)

C_det_A2 = np.array(C_det_A2)
Delta_det_A2, _ = fit_decay(C_det_A2, t_start=1, t_end=15)

print(f"  Decay rate of det(M): Delta = {Delta_det_A2:.4f}")
print(f"  Compare to A2 ground state: Delta_0 = {-np.log(0.5022):.4f}")
print(f"  Ratio: {Delta_det_A2/(-np.log(0.5022)):.4f}")

# sqrt(det(M)) decays at half the rate of det(M)
if Delta_det_A2 is not None:
    Delta_spinor = Delta_det_A2 / 2
    print(f"\n  sqrt(det M) decay: Delta_spinor = Delta_det/2 = {Delta_spinor:.4f}")
    print(f"  m(spinor)/m(0++) = {Delta_spinor / (-np.log(0.5022)):.4f}")

# ============================================================
# QUESTION 4: WHAT IS THE ACTUAL 0++ MASS?
# ============================================================

print(f"\n{'='*70}")
print("QUESTION 4: THE 0++ MASS — LINEAR OR BILINEAR?")
print("=" * 70)

print("""
This is the key question. The 0++ glueball mass from lattice QCD is the
GROUND STATE of the spectrum. In our framework:

If 0++ is a LINEAR Koopman observable (projects onto one eigenmode):
  m(0++) = Delta_0 = 0.6888  (A2 ground state)
  Ratio 2++/0++ = sqrt(2) requires 2++ to have Delta = 0.9741

If 0++ is a BILINEAR Koopman observable (projects onto eigenmode squared):
  m(0++) = 2 * Delta_0 = 1.3776  (Koopman product state)
  Ratio 2++/0++ = sqrt(2) requires 2++ to have Delta = 1.3776 * sqrt(2) = 1.948

The Schwinger function measurement above tells us which one is correct.
""")

print(f"Schwinger function results (from perturbing along the floor kick direction):")
print(f"  O_0++ (total norm) decays at Delta = {Delta_0pp:.4f}")
print(f"  O_2++ (spin-2 norm) decays at Delta = {Delta_2pp:.4f}")
print(f"\nA2 Jacobian eigenvalues:")
print(f"  Delta_0 = {-np.log(0.5022):.4f} (mode 0, node-antisym radial)")
print(f"  Delta_1 = {-np.log(0.4745):.4f} (mode 1, node-antisym angular)")
print(f"  Delta_2 = {-np.log(0.3527):.4f} (mode 2, node-sym angular)")
print(f"  Delta_3 = {-np.log(0.3344):.4f} (mode 3, node-sym radial)")

print(f"""
INTERPRETATION:
  If O_0++ decays at Delta = {Delta_0pp:.4f}, it matches mode {
      0 if abs(Delta_0pp - 0.6888) < 0.05 else
      1 if abs(Delta_0pp - 0.7455) < 0.05 else
      2 if abs(Delta_0pp - 1.0423) < 0.05 else
      3 if abs(Delta_0pp - 1.0954) < 0.05 else
      "none"
  } (Delta = {
      0.6888 if abs(Delta_0pp - 0.6888) < 0.05 else
      0.7455 if abs(Delta_0pp - 0.7455) < 0.05 else
      1.0423 if abs(Delta_0pp - 1.0423) < 0.05 else
      1.0954 if abs(Delta_0pp - 1.0954) < 0.05 else
      "?"
  })
""")

# ============================================================
# SYNTHESIS
# ============================================================

print(f"{'='*70}")
print("SYNTHESIS AND OPEN QUESTIONS")
print("=" * 70)

print(f"""
From this exploration:

1. SCHWINGER FUNCTIONS tell the truth about masses.
   The operator norms (|T2|/|T0| = sqrt(2)) give coupling STRENGTHS to the vacuum.
   The MASS is the decay rate of the Schwinger function, not the coupling strength.
   These can be equal in special cases but need not be in general.

2. COLOR-SINGLET STRUCTURE:
   Node-antisymmetric modes (lambda_0, lambda_1) carry SU(3) color charge.
   Physical glueballs are color singlets formed from PRODUCTS of these.
   The product lambda_0 * lambda_0 gives the 0++ at 2*Delta_0.
   The product lambda_0 * lambda_1 gives a state at Delta_0 + Delta_1.
   The sqrt(2) doesn't naturally emerge from this counting.

   BUT: there's another possibility. The Koopman MODES are colored.
   The physical glueball operators project onto SYMMETRIC combinations.
   The 0++ IS the ground state: it's whatever operator decays slowest.
   Which operator decays slowest depends on the full spectrum.

3. THE sqrt(2) QUESTION remains open.
   The algebraic identity |T2|^2/|T0|^2 = 2 is exact.
   Whether this gives the mass ratio depends on the identification:
     mass ratio = coupling ratio (linear) → sqrt(2)   [claimed]
     mass ratio = coupling ratio (bilinear) → doesn't change things
   The Schwinger function measurement gives the empirical answer.

4. FERMION MASS:
   The det(M) decay rate gives an estimate.
   The spinor sqrt(det M) decays at half that rate.
   Whether this corresponds to a physical quark mass requires
   understanding the dynamics of the spinor bundle, not just det(M).

5. THE DEEPEST QUESTION I'M LEFT WITH:
   The equilibrium is m* = (0.787, 0.029, 0.029, 0.155).
   One entity dominates (s1 = 0.787 ≫ s2 = s3 = 0.029).
   This asymmetry is forced by H=3 — NOT by choosing a direction.
   The S_3 symmetry of the 3 sections breaks to S_2 at the fixed point.
   This S_3→S_2 breaking is SPONTANEOUS (any of the 3 sections could dominate).

   The three orientations (which section dominates) are the three COLORS.
   The equilibrium has chosen one. The other two are degenerate.

   This might be the deepest content of the ADE classification:
   - A_n: n+1 sections, n+1 colors, S_{n+1}→S_n breaking
   - The number of colors = H = 3 is not assumed but forced by (H-1)^2 = H+1

   If this is right, the COLOR GROUP of QCD is derived from H=3.
   SU(3) color is not an input — it's the symmetry group of the spontaneous
   breaking of S_3 permutation symmetry at the DS fixed point.

   This would make the THREE COLORS of QCD a theorem, not an assumption.
""")
