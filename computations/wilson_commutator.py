"""
Wilson loop from the commutator [M,E] — the pure rotation content.

The DS-Matrix Decomposition (Theorem dsmatrix):
  √2 M''_pre = {M,E} + D - (s·e)I

The commutator [M,E] = i(s×e)·σ is ABSENT from the DS output
but computable from the inputs. It is the su(2)-valued curvature.

The Wilson loop is the path-ordered product of exp([M,E]) around
a closed contour — the holonomy of the rotation content only,
with the stretching (anticommutator) removed.
"""

import numpy as np
from numpy.linalg import norm, det, eigvals
from scipy.linalg import expm

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

def pauli_matrix(m):
    th, s1, s2, s3 = m
    return np.array([
        [th + s3, s1 - 1j*s2],
        [s1 + 1j*s2, th - s3]
    ], dtype=complex) / np.sqrt(2)

def commutator_ME(m, e):
    """Compute [M, E] = i(s×e)·σ, the su(2) curvature content.

    Returns a 2×2 traceless anti-Hermitian matrix (Lie algebra element).
    """
    s = m[1:]  # s₁, s₂, s₃
    se = e[1:]  # e₁, e₂, e₃

    # Cross product s × e
    cross = np.array([
        s[1]*se[2] - s[2]*se[1],
        s[2]*se[0] - s[0]*se[2],
        s[0]*se[1] - s[1]*se[0],
    ])

    # [M, E] = i(s×e)·σ / 2  (the 1/2 from M = (θI+s·σ)/√2 for both)
    # Actually: M = (θI + s·σ)/√2, E = (φI + se·σ)/√2
    # ME = (1/2)(θφI + θ·se·σ + φ·s·σ + (s·σ)(se·σ))
    # (s·σ)(se·σ) = (s·se)I + i(s×se)·σ
    # [M,E] = ME - EM = i(s×se)·σ (the symmetric parts cancel)
    # With the √2 normalizations: [M,E] = i(s×se)·σ / 2... let me just compute it.

    M = pauli_matrix(m)
    E = pauli_matrix(e)
    comm = M @ E - E @ M

    return comm

def link_holonomy(m_site, m_neighbor):
    """The SU(2) holonomy element for the link between two sites.

    Extracts the commutator [M, E], which is in su(2) (the Lie algebra),
    and exponentiates to get the group element.
    """
    comm = commutator_ME(m_site, m_neighbor)

    # comm is anti-Hermitian (su(2) element): comm = iH where H is Hermitian
    # exp(comm) ∈ SU(2) (unitary, det=1)
    U = expm(comm)
    return U

def wilson_loop_commutator(lat_flat, contour):
    """Wilson loop from path-ordered product of link holonomies.

    Each link contributes exp([M_site, M_neighbor]) ∈ SU(2).
    The Wilson loop is Tr(∏ U_link) / 2.
    """
    W = np.eye(2, dtype=complex)

    for i in range(len(contour) - 1):
        site = contour[i]
        neighbor = contour[i + 1]
        U = link_holonomy(lat_flat[site], lat_flat[neighbor])
        W = W @ U

    # Close the loop
    U_close = link_holonomy(lat_flat[contour[-1]], lat_flat[contour[0]])
    W = W @ U_close

    return np.real(np.trace(W)) / 2

# ============================================================
# TEST
# ============================================================
print("=" * 70)
print("WILSON LOOP FROM COMMUTATOR [M,E]")
print("=" * 70)

N = 16
np.random.seed(42)

lat = {}
for i in range(N):
    for j in range(N):
        lat[(i,j)] = born_floor(np.random.dirichlet([2,2,2,2]))

def make_contour(x, y, R, T, N):
    contour = []
    for dx in range(R+1):
        contour.append(((x+dx)%N, y))
    for dy in range(1, T+1):
        contour.append(((x+R)%N, (y+dy)%N))
    for dx in range(1, R+1):
        contour.append(((x+R-dx)%N, (y+T)%N))
    for dy in range(1, T+1):
        contour.append((x, (y+T-dy)%N))
    return contour

# Verify link holonomy is SU(2)
m1 = born_floor(np.random.dirichlet([2,2,2,2]))
m2 = born_floor(np.random.dirichlet([2,2,2,2]))
U = link_holonomy(m1, m2)
print(f"Link holonomy check:")
print(f"  det(U) = {det(U):.8f} (should be ~1)")
print(f"  U·U† = {np.max(np.abs(U @ U.conj().T - np.eye(2))):.2e} (should be ~0)")
print(f"  Tr(U) = {np.trace(U):.6f}")

comm = commutator_ME(m1, m2)
print(f"\n  [M,E] trace = {np.trace(comm):.2e} (should be 0 — traceless)")
print(f"  [M,E] + [M,E]† = {np.max(np.abs(comm + comm.conj().T)):.2e} (should be 0 — anti-Hermitian)")

# Measure Wilson loops
n_meas = 200
print(f"\n{'R×T':>6} {'Area':>5} {'Perim':>5} {'<W>':>10} {'<|W|>':>10} "
      f"{'-ln<|W|>':>10} {'-ln/A':>8}")

results = []
for R in range(1, 8):
    for T in range(R, 8):
        Ws = []
        for _ in range(n_meas):
            x, y = np.random.randint(0, N, 2)
            contour = make_contour(x, y, R, T, N)
            W = wilson_loop_commutator(lat, contour)
            Ws.append(W)

        W_mean = np.mean(Ws)
        W_abs = np.mean(np.abs(Ws))
        area = R * T
        perim = 2*(R+T)

        if W_abs > 1e-10 and W_abs < 10:
            logW = -np.log(W_abs)
        else:
            logW = float('nan')

        results.append((R, T, area, perim, W_mean, W_abs, logW))

        if T == R or T <= 3:
            print(f"{R}×{T:>2} {area:5d} {perim:5d} {W_mean:10.6f} {W_abs:10.6f} "
                  f"{logW:10.4f} {logW/area if area > 0 and not np.isnan(logW) else float('nan'):8.4f}")

# ============================================================
# FIT
# ============================================================
print("\n" + "=" * 70)
print("AREA LAW VS PERIMETER LAW")
print("=" * 70)

areas = np.array([r[2] for r in results])
perims = np.array([r[3] for r in results])
logWs = np.array([r[6] for r in results])

mask = np.isfinite(logWs) & (logWs > 0)
if mask.sum() > 3:
    a_f, p_f, l_f = areas[mask], perims[mask], logWs[mask]

    A_m = np.column_stack([a_f, np.ones(len(a_f))])
    c_a = np.linalg.lstsq(A_m, l_f, rcond=None)[0]
    ss_r_a = np.sum((l_f - A_m @ c_a)**2)
    ss_t = np.sum((l_f - np.mean(l_f))**2)
    R2_a = 1 - ss_r_a/ss_t if ss_t > 0 else 0

    P_m = np.column_stack([p_f, np.ones(len(p_f))])
    c_p = np.linalg.lstsq(P_m, l_f, rcond=None)[0]
    ss_r_p = np.sum((l_f - P_m @ c_p)**2)
    R2_p = 1 - ss_r_p/ss_t if ss_t > 0 else 0

    B_m = np.column_stack([a_f, p_f, np.ones(len(a_f))])
    c_b = np.linalg.lstsq(B_m, l_f, rcond=None)[0]
    ss_r_b = np.sum((l_f - B_m @ c_b)**2)
    R2_b = 1 - ss_r_b/ss_t if ss_t > 0 else 0

    print(f"Area fit:      -ln<|W|> = {c_a[0]:.6f}·A + {c_a[1]:.6f},  R² = {R2_a:.4f}")
    print(f"Perimeter fit: -ln<|W|> = {c_p[0]:.6f}·P + {c_p[1]:.6f},  R² = {R2_p:.4f}")
    print(f"Combined:      -ln<|W|> = {c_b[0]:.6f}·A + {c_b[1]:.6f}·P + {c_b[2]:.6f},  R² = {R2_b:.4f}")

    print(f"\nString tension σ = {c_a[0]:.6f}")
    print(f"Predicted σ = -ln(23/30) = {-np.log(23/30):.6f}")

    if R2_a > R2_p:
        print(f"\n*** AREA LAW WINS: R²(area)={R2_a:.4f} > R²(perim)={R2_p:.4f} ***")
    else:
        print(f"\n*** PERIMETER LAW WINS: R²(perim)={R2_p:.4f} > R²(area)={R2_a:.4f} ***")
else:
    print("Not enough valid data points")
