"""
Physical observables along the K*=7/30 moduli curve.

For each point on the K*=7/30 curve (parametrised by φ), compute:
  1. Transfer operator eigenvalues λ₀, λ₁, λ₂
  2. Penrose residue ||ρ₋₁|| and |F⁺|²
  3. SO(4) decomposition of d̄Φ
  4. Commutator content C·det(M*)²
  5. Eigenvalue splitting
  6. Non-abelian fraction of the kink

Everything computed at the K*=7/30 equilibrium for each evidence distribution.
"""

import numpy as np
from scipy.optimize import brentq

H = 3
FLOOR = 1.0 / H**3

# Pauli matrices
sig = [np.array([[0,1],[1,0]], dtype=complex),
       np.array([[0,-1j],[1j,0]], dtype=complex),
       np.array([[1,0],[0,-1]], dtype=complex)]
I2 = np.eye(2, dtype=complex)

def to_M(m):
    return (m[3]*I2 + m[0]*sig[0] + m[1]*sig[1] + m[2]*sig[2]) / np.sqrt(2)

def enforce_floor(m):
    s, theta = m[:3], m[3]
    ssq = sum(si**2 for si in s)
    if theta**2 / (ssq + theta**2) >= FLOOR:
        return m.copy()
    ss = sum(s)
    if ss < 1e-15: return m.copy()
    r = ssq / ss**2
    disc = (2*r)**2 + 4*(26-r)*r
    t = (-2*r + np.sqrt(disc)) / (2*(26-r))
    alpha = (1-t) / ss
    return np.array([s[0]*alpha, s[1]*alpha, s[2]*alpha, t])

def ds_step(m, e):
    s, theta = m[:3], m[3]
    se, phi = e[:3], e[3]
    s_new = s*se + s*phi + theta*se
    theta_new = theta*phi
    K = sum(s[i]*se[j] for i in range(3) for j in range(3) if i != j)
    d = 1.0 - K
    if abs(d) < 1e-15: return m.copy(), K
    m_out = np.zeros(4)
    m_out[:3] = s_new/d
    m_out[3] = theta_new/d
    m_out = m_out / np.sum(m_out)
    return enforce_floor(m_out), K

def find_fp(e, n=1000):
    m = np.array([0.7, 0.05, 0.05, 0.2])
    for _ in range(n):
        m2, K = ds_step(m, e)
        if np.max(np.abs(m2-m)) < 1e-14: break
        m = m2
    return m, sum(m[i]*e[j] for i in range(3) for j in range(3) if i != j)

def K_at_pdom_phi(p_dom, phi):
    p_w = (1-p_dom)/2
    e = np.array([p_dom*(1-phi), p_w*(1-phi), p_w*(1-phi), phi])
    e = e / np.sum(e)
    m, K = find_fp(e)
    return K

def wirtinger_jacobian(m, e, eps=1e-8):
    """Anti-holomorphic Jacobian d̄Φ via Wirtinger derivatives.
    d/dz̄ = (1/2)(d/dx + i d/dy) for z = x + iy."""
    m_c = m.astype(complex)
    e_c = e.astype(complex)

    def Phi(m_in):
        out, _ = ds_step(np.real(m_in), np.real(e_c))
        return out

    J_anti = np.zeros((4, 4), dtype=complex)
    m0 = Phi(m_c)

    for j in range(4):
        # d/dx_j
        m_px = m_c.copy(); m_px[j] += eps
        m_mx = m_c.copy(); m_mx[j] -= eps
        dfdx = (Phi(m_px) - Phi(m_mx)) / (2*eps)

        # d/dy_j (imaginary perturbation)
        m_py = m_c.copy(); m_py[j] += 1j*eps
        m_my = m_c.copy(); m_my[j] -= 1j*eps
        dfdy = (Phi(m_py) - Phi(m_my)) / (2*eps)

        # Wirtinger: d/dz̄ = (1/2)(d/dx + i d/dy)
        J_anti[:, j] = 0.5 * (dfdx + 1j * dfdy)

    return J_anti

def so4_decomposition(J_anti):
    """Decompose d̄Φ into SO(4) irreps: (1,1) + (3,1)⊕(1,3) + (3,3)."""
    J = np.real(J_anti)  # Take real part for SO(4) decomposition
    n = J.shape[0]
    if n != 4:
        return 0, 0, 0

    # (1,1): scalar = trace/4
    scalar = np.trace(J) / 4 * np.eye(4)

    # Antisymmetric: (3,1)⊕(1,3)
    antisym = 0.5 * (J - J.T)

    # Symmetric traceless: (3,3)
    sym_traceless = 0.5 * (J + J.T) - scalar

    norm_total = np.linalg.norm(J, 'fro')**2
    if norm_total < 1e-30:
        return 0, 0, 0

    f_scalar = np.linalg.norm(scalar, 'fro')**2 / norm_total
    f_gauge = np.linalg.norm(antisym, 'fro')**2 / norm_total
    f_grav = np.linalg.norm(sym_traceless, 'fro')**2 / norm_total

    return f_scalar, f_gauge, f_grav

# ============================================================
# Compute observables along the K*=7/30 curve
# ============================================================
print("=" * 100)
print("PHYSICAL OBSERVABLES ALONG THE K*=7/30 MODULI CURVE")
print("=" * 100)

results = []

phi_values = np.linspace(0.005, 0.60, 60)

for phi in phi_values:
    try:
        def obj(p): return K_at_pdom_phi(p, phi) - 7.0/30
        p_lo, p_hi = 0.35, 0.99
        if obj(p_lo) * obj(p_hi) > 0: continue
        p_dom = brentq(obj, p_lo, p_hi, xtol=1e-12)

        p_w = (1-p_dom)/2
        e = np.array([p_dom*(1-phi), p_w*(1-phi), p_w*(1-phi), phi])
        e = e / np.sum(e)
        m, K = find_fp(e)
        if abs(K - 7/30) > 1e-6: continue

        # 1. Transfer operator eigenvalues (3x3 projected Jacobian)
        eps = 1e-8
        m0, _ = ds_step(m, e)
        J3 = np.zeros((3, 3))
        for j in range(3):
            idx = [0, 1, 3][j]
            m_p = m.copy(); m_p[idx] += eps; m_p[3 if idx != 3 else 0] -= eps
            mp_out, _ = ds_step(m_p, e)
            for i in range(3):
                J3[i, j] = (mp_out[[0,1,3][i]] - m0[[0,1,3][i]]) / eps
        eigs_3 = np.sort(np.abs(np.linalg.eigvals(J3)))[::-1]
        lam0, lam1, lam2 = eigs_3[0], eigs_3[1], eigs_3[2]
        gap = -np.log(lam0) if 0 < lam0 < 1 else 0
        splitting = abs(lam0-lam1)/lam0 if lam0 > 0 else 0

        # 2. Anti-holomorphic Jacobian and Penrose residue
        J_anti = wirtinger_jacobian(m, e)
        dbar_norm = np.linalg.norm(J_anti, 'fro')

        # Singular values of d̄Φ (rank structure)
        sv = np.linalg.svd(np.real(J_anti), compute_uv=False)
        rank1_ratio = sv[1]/sv[0] if sv[0] > 1e-10 else 0

        # Penrose residue: d̄M = (∂M/∂m)·d̄Φ
        # ∂M/∂m maps 4-vector to 2x2 matrix via Pauli embedding
        dbar_m = np.real(J_anti)  # Use real part
        # The image in M₂(C): δM = (δθ·I + δs₁·σ₁ + δs₂·σ₂ + δs₃·σ₃)/√2
        # For each column of d̄Φ, compute the 2x2 matrix
        M_star = to_M(m)
        M_inv = np.linalg.inv(M_star) if abs(np.linalg.det(M_star)) > 1e-10 else None

        rho_norm = 0
        F_plus_sq = 0
        comm_norm = 0

        if M_inv is not None:
            # Compute the (0,1)-form a* = M*⁻¹ · Pauli(d̄Φ)
            # For each anti-holomorphic direction α:
            a_star = []
            for alpha in range(4):
                delta_m = dbar_m[:, alpha]
                delta_M = (delta_m[3]*I2 + delta_m[0]*sig[0] +
                          delta_m[1]*sig[1] + delta_m[2]*sig[2]) / np.sqrt(2)
                a_alpha = M_inv @ delta_M
                a_star.append(a_alpha)

            # Penrose residue ~ ||a*|| (leading order)
            a_norms = [np.linalg.norm(a, 'fro') for a in a_star]
            rho_norm = np.sqrt(sum(n**2 for n in a_norms))
            F_plus_sq = rho_norm**2

            # Commutator content: ||[a*_α, a*_β]||
            comm_total = 0
            for i in range(4):
                for j in range(i+1, 4):
                    comm = a_star[i] @ a_star[j] - a_star[j] @ a_star[i]
                    comm_total += np.linalg.norm(comm, 'fro')**2
            comm_norm = np.sqrt(comm_total)

        # 3. SO(4) decomposition
        f_scalar, f_gauge, f_grav = so4_decomposition(J_anti)

        # 4. Commutator content C·det(M*)²
        det_M = np.linalg.det(M_star)
        # C·det² from the paper's formula
        n_born = H**3 - 1  # = 26
        C_det_sq = 4*(3*n_born**2 + 2*n_born + 3) / (n_born - 1)**2  # = 8332/625
        # But this is the ALGEBRAIC identity — same everywhere at Born floor

        # Cross product |s×e|
        sx = np.array([m[1]*e[2]-m[2]*e[1], m[2]*e[0]-m[0]*e[2], m[0]*e[1]-m[1]*e[0]])
        cross_sq = np.sum(sx**2)

        results.append({
            'phi': phi, 'p_dom': p_dom, 'K': K,
            'lam0': lam0, 'lam1': lam1, 'lam2': lam2,
            'gap': gap, 'splitting': splitting,
            'dbar_norm': dbar_norm, 'rank1_ratio': rank1_ratio,
            'rho_norm': rho_norm, 'F_plus_sq': F_plus_sq,
            'comm_norm': comm_norm,
            'f_scalar': f_scalar, 'f_gauge': f_gauge, 'f_grav': f_grav,
            'det_M': abs(det_M), 'cross_sq': cross_sq,
            'theta': m[3], 'm': m, 'e': e,
        })
    except Exception as ex:
        continue

print(f"\nComputed observables at {len(results)} points on the K*=7/30 curve\n")

# ============================================================
# Display results
# ============================================================
print(f"{'phi':>6s} {'p_dom':>6s} {'lam0':>8s} {'lam1':>8s} {'Delta':>6s} "
      f"{'split%':>6s} {'||dbar||':>8s} {'rk1_rat':>8s} "
      f"{'||rho||':>8s} {'|F+|^2':>8s} {'[a,a]':>8s} "
      f"{'scalar':>6s} {'gauge':>6s} {'grav':>6s} {'|sxe|^2':>8s}")
print("-" * 140)

for r in results[::max(1, len(results)//30)]:
    print(f"{r['phi']:6.3f} {r['p_dom']:6.3f} {r['lam0']:8.4f} {r['lam1']:8.4f} "
          f"{r['gap']:6.3f} {r['splitting']*100:6.2f} {r['dbar_norm']:8.4f} "
          f"{r['rank1_ratio']:8.5f} {r['rho_norm']:8.4f} {r['F_plus_sq']:8.4f} "
          f"{r['comm_norm']:8.4f} {r['f_scalar']:6.3f} {r['f_gauge']:6.3f} "
          f"{r['f_grav']:6.3f} {r['cross_sq']:8.6f}")

# ============================================================
# Invariance analysis
# ============================================================
print("\n" + "=" * 100)
print("WHAT IS CONSTANT vs WHAT VARIES along the K*=7/30 curve")
print("=" * 100)

quantities = {
    'K*': [r['K'] for r in results],
    'λ₀': [r['lam0'] for r in results],
    'λ₁': [r['lam1'] for r in results],
    'Δ': [r['gap'] for r in results],
    'splitting': [r['splitting'] for r in results],
    '||d̄Φ||': [r['dbar_norm'] for r in results],
    'rank-1 ratio': [r['rank1_ratio'] for r in results],
    '||ρ₋₁||': [r['rho_norm'] for r in results],
    '|F⁺|²': [r['F_plus_sq'] for r in results],
    '||[a,a]||': [r['comm_norm'] for r in results],
    'scalar (1,1)': [r['f_scalar'] for r in results],
    'gauge (3,1)⊕(1,3)': [r['f_gauge'] for r in results],
    'graviton (3,3)': [r['f_grav'] for r in results],
    '|det(M*)|': [r['det_M'] for r in results],
    '|s×e|²': [r['cross_sq'] for r in results],
    'θ*': [r['theta'] for r in results],
    'σ = -ln(23/30)': [0.26593 for r in results],  # constant by construction
    'C·det²=8332/625': [8332/625 for r in results],  # algebraic identity
}

print(f"\n{'Quantity':>25s} {'Mean':>12s} {'Std':>12s} {'CV':>8s} {'Min':>10s} {'Max':>10s}  Status")
print("-" * 105)

for name, vals in quantities.items():
    vals = np.array([v for v in vals if np.isfinite(v) and v != 0])
    if len(vals) == 0: continue
    mn, sd = np.mean(vals), np.std(vals)
    cv = sd/abs(mn) if abs(mn) > 1e-10 else float('inf')
    status = "*** UNIVERSAL ***" if cv < 0.005 else ("~ universal" if cv < 0.03 else "GAUGE-SPECIFIC")
    print(f"{name:>25s} {mn:12.6f} {sd:12.6f} {cv:8.4f} {vals.min():10.6f} {vals.max():10.6f}  {status}")

# ============================================================
# Key relationships
# ============================================================
print("\n" + "=" * 100)
print("KEY RELATIONSHIPS along the curve")
print("=" * 100)

phis = np.array([r['phi'] for r in results])
lam0s = np.array([r['lam0'] for r in results])
rhos = np.array([r['rho_norm'] for r in results])
comms = np.array([r['comm_norm'] for r in results])
f_gravs = np.array([r['f_grav'] for r in results])
f_gauges = np.array([r['f_gauge'] for r in results])
dbars = np.array([r['dbar_norm'] for r in results])
gaps = np.array([r['gap'] for r in results])
crosses = np.array([r['cross_sq'] for r in results])

# λ₀ vs φ
c = np.polyfit(phis, lam0s, 1)
print(f"\nλ₀ ≈ {c[0]:.4f}·φ + {c[1]:.4f} (linear fit)")
c2 = np.polyfit(phis, lam0s, 2)
print(f"λ₀ ≈ {c2[0]:.4f}·φ² + {c2[1]:.4f}·φ + {c2[2]:.4f} (quadratic fit)")

# ||ρ₋₁|| vs φ
if np.any(rhos > 0):
    mask = rhos > 0
    c_rho = np.polyfit(phis[mask], rhos[mask], 1)
    print(f"\n||ρ₋₁|| ≈ {c_rho[0]:.4f}·φ + {c_rho[1]:.4f}")

# Graviton fraction vs φ
if np.any(f_gravs > 0):
    mask = f_gravs > 0
    c_grav = np.polyfit(phis[mask], f_gravs[mask], 1)
    print(f"\nGraviton fraction ≈ {c_grav[0]:.4f}·φ + {c_grav[1]:.4f}")

# |s×e|² vs λ₀
if np.any(crosses > 0) and np.any(lam0s > 0):
    mask = (crosses > 0) & (lam0s > 0)
    c_cross = np.polyfit(crosses[mask], lam0s[mask], 1)
    print(f"\nλ₀ ≈ {c_cross[0]:.4f}·|s×e|² + {c_cross[1]:.4f}")

# ||d̄Φ|| vs φ
if np.any(dbars > 0):
    print(f"\n||d̄Φ|| range: [{dbars.min():.4f}, {dbars.max():.4f}]")
    c_db = np.polyfit(phis, dbars, 1)
    print(f"||d̄Φ|| ≈ {c_db[0]:.4f}·φ + {c_db[1]:.4f}")
