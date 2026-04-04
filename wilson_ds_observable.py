"""
Wilson loop as a DS observable — sub-problem (i).

The audit's most damning point: W(C) was never defined in DS terms.

In standard gauge theory:
  W(C) = Tr P exp(ig ∮_C A·dx)
  On the lattice: W(C) = Tr ∏_{links} U_link

In the DS framework:
  The gauge connection A is encoded in the holomorphic bundle E → CP³.
  The Ward correspondence reconstructs A from the transition functions.
  At the equilibrium, dW is an isomorphism (Theorem ward_spectral).

So the Wilson loop in DS terms is:
  1. At each link (x, x+μ), the connection A_μ(x) is reconstructed
     from the bundle transition function via Ward.
  2. The transition function IS M(x)^{-1} M(x+μ) — but NOT as a
     pure gauge (which telescopes). Rather, it's the HOLOMORPHIC
     part of the transition function, with the anti-holomorphic
     part (the floor correction) producing the curvature.
  3. The Wilson loop is the trace of the path-ordered product of
     these connection elements.

Wait. The problem is subtler. Let me think about what the connection
actually IS in the DS framework.

The bundle E → CP³ has transition function M(Z) at each point Z ∈ CP³.
The connection on E is:
  a^{1,0} = M^{-1} ∂M     (holomorphic part → F^-)
  a^{0,1} = M^{-1} ∂̄M    (anti-holomorphic part → F^+)

The FULL connection is a = a^{1,0} + a^{0,1}.
At the equilibrium: a^{0,1} = 0 (floor inactive), so a = a^{1,0}.
Away from equilibrium: both parts contribute.

The Wilson loop uses the FULL connection:
  W(C) = Tr P exp(∮_C a)

But this is on CP³, not on S⁴. The spacetime Wilson loop on S⁴
is obtained by restricting to a contour C ⊂ S⁴ and lifting it to
CP³ via the twistor fibration.

For each link (x, x+μ) on S⁴, the parallel transport in E is:
  g(x, x+μ) = P exp(∫_x^{x+μ} a)

In the DS framework, this parallel transport IS the DS propagation:
combining m(x) with the evidence from the gauge field along the link.

Let me make this precise.
"""

import numpy as np
from numpy.linalg import norm, det, inv

def ds_combine(m, e):
    th_m, s1m, s2m, s3m = m
    th_e, s1e, s2e, s3e = e
    K = s1m*s2e + s1m*s3e + s2m*s1e + s2m*s3e + s3m*s1e + s3m*s2e
    if abs(K) >= 1.0:
        return m.copy()
    fac = 1.0 / (1.0 - K)
    return np.array([
        fac * th_m * th_e,
        fac * (s1m*s1e + s1m*th_e + th_m*s1e),
        fac * (s2m*s2e + s2m*th_e + th_m*s2e),
        fac * (s3m*s3e + s3m*th_e + th_m*s3e),
    ])

def born_floor(m, H=3):
    th = m[0]
    s = m[1:]
    s_sq = np.dot(s, s)
    total_sq = th**2 + s_sq
    if total_sq < 1e-30:
        return np.array([1.0, 0.0, 0.0, 0.0])
    born = th**2 / total_sq
    if born < 1.0/H**3:
        th_new = np.sqrt(s_sq / (H**3 - 1))
        m = np.array([th_new, m[1], m[2], m[3]])
    total = np.sum(np.abs(m))
    if total > 0:
        m = m / total
    return m

def ds_step(m, e):
    return born_floor(ds_combine(m, e))

def pauli_matrix(m):
    th, s1, s2, s3 = m
    return np.array([
        [th + s3, s1 - 1j*s2],
        [s1 + 1j*s2, th - s3]
    ], dtype=complex) / np.sqrt(2)

# ============================================================
# THE DS PARALLEL TRANSPORTER
# ============================================================
print("=" * 70)
print("THE DS PARALLEL TRANSPORTER")
print("=" * 70)

# The parallel transport from site x to site y in the DS framework
# is NOT M(x)^{-1} M(y) (pure gauge, trivial holonomy).
#
# It IS the DS combination map viewed as a linear operator on M_2(C).
#
# When site x combines with evidence from site y:
#   m'(x) = DS(m(x), m(y))
#
# In the Pauli matrix representation:
#   M'(x) = (1/(1-K)) * [{M(x), M(y)}/√2 + D - (s·e)I] / √2
#   (from the DS-Matrix Decomposition, Theorem dsmatrix)
#
# This is a LINEAR map on M(x) for fixed M(y):
#   M'(x) = T_{xy} · M(x)
# where T_{xy} is a linear operator on 2×2 matrices.
#
# T_{xy} IS the parallel transporter. It's not a group element —
# it's a linear map on the fibre. In gauge theory terms, it's
# the holonomy of the connection along the link (x,y).
#
# The Wilson loop is: start with a test matrix, apply T_{12},
# then T_{23}, ..., then T_{n1}, and take the trace.

def ds_transporter(m_evidence):
    """The DS parallel transporter for fixed evidence m_e.

    Returns a function that maps M (2×2) to M' (2×2).
    This is the linearisation of the DS step in the Pauli basis.
    """
    e = m_evidence
    th_e, s1e, s2e, s3e = e

    def transport(M_in, m_state):
        # DS combine m_state with evidence e
        m_out = ds_step(m_state, e)
        M_out = pauli_matrix(m_out)
        return M_out

    return transport

def ds_wilson_loop(lattice, contour):
    """Compute Wilson loop by DS-propagating a test state around contour.

    lattice: array of mass functions at each site
    contour: list of site indices forming a closed loop

    The Wilson loop = Tr(M_final) / Tr(M_initial)
    where M_initial is the Pauli matrix at the starting site
    and M_final is the result of DS-propagating around the loop.
    """
    if len(contour) < 2:
        return 1.0

    # Start with the mass function at the first site
    m = lattice[contour[0]].copy()
    M_initial = pauli_matrix(m)

    # Propagate around the loop: at each step, use the NEXT site as evidence
    for i in range(len(contour) - 1):
        evidence = lattice[contour[i+1]]
        m = ds_step(m, evidence)

    # Close the loop: use the first site as evidence one more time
    # (the contour includes the return to start)
    M_final = pauli_matrix(m)

    # The Wilson loop is how much the state changed
    # W = Tr(M_initial^{-1} M_final) / 2
    if abs(det(M_initial)) > 1e-10:
        holonomy = inv(M_initial) @ M_final
        W = np.trace(holonomy) / 2
        return np.real(W)
    else:
        return 1.0

# ============================================================
# TEST ON A 2D LATTICE
# ============================================================
print("\n" + "=" * 70)
print("WILSON LOOP ON 2D LATTICE (DS PROPAGATION)")
print("=" * 70)

N = 16
np.random.seed(42)

# Random configuration on B
lat = np.zeros((N, N, 4))
for i in range(N):
    for j in range(N):
        lat[i, j] = born_floor(np.random.dirichlet([2, 2, 2, 2]))

def make_contour(x, y, R, T, N):
    """Rectangle contour starting at (x,y), size R×T."""
    contour = []
    # Right
    for dx in range(R+1):
        contour.append(((x+dx)%N, y))
    # Up
    for dy in range(1, T+1):
        contour.append(((x+R)%N, (y+dy)%N))
    # Left
    for dx in range(1, R+1):
        contour.append(((x+R-dx)%N, (y+T)%N))
    # Down (back to start)
    for dy in range(1, T+1):
        contour.append((x, (y+T-dy)%N))
    return contour

# Flatten lattice for easy indexing
lat_flat = {}
for i in range(N):
    for j in range(N):
        lat_flat[(i,j)] = lat[i,j]

def wilson_on_lattice(lat_flat, x, y, R, T, N):
    contour = make_contour(x, y, R, T, N)
    # Convert to flat indexing
    m = lat_flat[contour[0]].copy()
    M_init = pauli_matrix(m)

    for i in range(1, len(contour)):
        evidence = lat_flat[contour[i]]
        m = ds_step(m, evidence)

    M_final = pauli_matrix(m)

    if abs(det(M_init)) > 1e-10:
        holonomy = inv(M_init) @ M_final
        return np.real(np.trace(holonomy)) / 2
    return 1.0

# Measure Wilson loops at various sizes
n_meas = 100

print(f"\n{'R×T':>6} {'Area':>5} {'Perim':>5} {'<W>':>10} {'<|W|>':>10} "
      f"{'-ln<|W|>':>10} {'-ln/A':>8}")

results = []
for R in range(1, 7):
    for T in [R]:  # Square loops for cleaner signal
        Ws = []
        for _ in range(n_meas):
            x, y = np.random.randint(0, N, 2)
            W = wilson_on_lattice(lat_flat, x, y, R, T, N)
            Ws.append(W)

        W_mean = np.mean(Ws)
        W_abs = np.mean(np.abs(Ws))
        area = R * T
        perim = 2*(R+T)
        logW = -np.log(W_abs) if W_abs > 1e-10 else float('inf')
        results.append((R, T, area, perim, W_mean, W_abs, logW))
        print(f"{R}×{T:>2} {area:5d} {perim:5d} {W_mean:10.6f} {W_abs:10.6f} "
              f"{logW:10.4f} {logW/area:8.4f}")

# Also non-square
print("\nNon-square loops:")
for R, T in [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]:
    Ws = []
    for _ in range(n_meas):
        x, y = np.random.randint(0, N, 2)
        W = wilson_on_lattice(lat_flat, x, y, R, T, N)
        Ws.append(W)

    W_abs = np.mean(np.abs(Ws))
    area = R * T
    perim = 2*(R+T)
    logW = -np.log(W_abs) if W_abs > 1e-10 else float('inf')
    results.append((R, T, area, perim, np.mean(Ws), W_abs, logW))
    print(f"{R}×{T:>2} {area:5d} {perim:5d} {np.mean(Ws):10.6f} {W_abs:10.6f} "
          f"{logW:10.4f} {logW/area:8.4f}")

# ============================================================
# FIT: area law vs perimeter law
# ============================================================
print("\n" + "=" * 70)
print("AREA LAW VS PERIMETER LAW")
print("=" * 70)

areas = np.array([r[2] for r in results])
perims = np.array([r[3] for r in results])
logWs = np.array([r[6] for r in results])

# Filter out any infinities
mask = np.isfinite(logWs) & (logWs > 0)
if mask.sum() > 2:
    a_fit = areas[mask]
    p_fit = perims[mask]
    l_fit = logWs[mask]

    # Area fit
    A_mat = np.column_stack([a_fit, np.ones(len(a_fit))])
    c_area = np.linalg.lstsq(A_mat, l_fit, rcond=None)[0]
    pred_a = A_mat @ c_area
    SS_res_a = np.sum((l_fit - pred_a)**2)
    SS_tot = np.sum((l_fit - np.mean(l_fit))**2)
    R2_a = 1 - SS_res_a/SS_tot if SS_tot > 0 else 0

    # Perimeter fit
    P_mat = np.column_stack([p_fit, np.ones(len(p_fit))])
    c_perim = np.linalg.lstsq(P_mat, l_fit, rcond=None)[0]
    pred_p = P_mat @ c_perim
    SS_res_p = np.sum((l_fit - pred_p)**2)
    R2_p = 1 - SS_res_p/SS_tot if SS_tot > 0 else 0

    # Combined
    B_mat = np.column_stack([a_fit, p_fit, np.ones(len(a_fit))])
    c_both = np.linalg.lstsq(B_mat, l_fit, rcond=None)[0]
    pred_b = B_mat @ c_both
    SS_res_b = np.sum((l_fit - pred_b)**2)
    R2_b = 1 - SS_res_b/SS_tot if SS_tot > 0 else 0

    print(f"Area law:      -ln<|W|> = {c_area[0]:.6f}·A + {c_area[1]:.6f},  R² = {R2_a:.4f}")
    print(f"Perimeter law: -ln<|W|> = {c_perim[0]:.6f}·P + {c_perim[1]:.6f},  R² = {R2_p:.4f}")
    print(f"Combined:      -ln<|W|> = {c_both[0]:.6f}·A + {c_both[1]:.6f}·P + {c_both[2]:.6f},  R² = {R2_b:.4f}")

    if R2_a > R2_p:
        print(f"\n*** AREA LAW WINS ***")
    else:
        print(f"\n*** PERIMETER LAW WINS ***")

    print(f"\nMeasured σ = {c_area[0]:.6f}")
    print(f"Predicted σ = -ln(23/30) = {-np.log(23/30):.6f}")
else:
    print("Not enough data points for fitting")
