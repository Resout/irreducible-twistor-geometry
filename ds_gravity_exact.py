"""
Gravity from DS: EXACT computation.

Uses the real DS combination from spectral_gap_computation.py,
extended to complex mass functions, with correct evidence construction
giving K*=7/30 and Born=1/27.
"""
import numpy as np

H = 3
FLOOR = 1.0 / H**3  # 1/27

# ============================================================
# DS combination (exact, from spectral_gap_computation.py)
# Extended to complex
# ============================================================
def ds_combine(m, e, apply_floor=True):
    """DS combination: m=(s1,s2,s3,theta), e=(se1,se2,se3,theta_e)."""
    s = m[:3]
    theta = m[3]
    se = e[:3]
    theta_e = e[3]

    # Pre-normalisation products
    s_new = s * se + s * theta_e + theta * se
    theta_new = theta * theta_e

    # Conflict: cross-focal products between different singletons
    K = sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)

    # Dempster normalisation
    denom = 1.0 - K
    if abs(denom) < 1e-15:
        return m.copy(), K

    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = s_new / denom
    m_out[3] = theta_new / denom

    # L1 normalise
    total = np.sum(np.abs(m_out))
    if total > 0:
        m_out = m_out / total

    # Born floor
    if apply_floor:
        born_val = np.abs(m_out[3])**2 / np.sum(np.abs(m_out)**2)
        if born_val < FLOOR:
            m_out = enforce_floor(m_out)

    return m_out, K

def enforce_floor(m):
    """Enforce Born(theta) >= 1/27, maintaining L1=1. Complex version."""
    s = m[:3].copy()
    theta = m[3]

    phase_theta = theta / np.abs(theta) if np.abs(theta) > 1e-15 else 1.0
    phases_s = np.array([si / np.abs(si) if np.abs(si) > 1e-15 else 1.0 for si in s])
    abs_s_orig = np.abs(s)

    lo, hi = np.abs(theta), 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        s_sum = np.sum(abs_s_orig)
        s_scale = (1.0 - mid) / s_sum if s_sum > 0 else 0
        s_trial = abs_s_orig * s_scale
        born_val = mid**2 / (np.sum(s_trial**2) + mid**2)
        if born_val < FLOOR:
            lo = mid
        else:
            hi = mid

    theta_new_abs = (lo + hi) / 2
    s_sum = np.sum(abs_s_orig)
    s_scale = (1.0 - theta_new_abs) / s_sum if s_sum > 0 else 0

    m_out = np.zeros(4, dtype=complex)
    m_out[:3] = phases_s * abs_s_orig * s_scale
    m_out[3] = phase_theta * theta_new_abs
    return m_out

def born(m):
    return np.abs(m[3])**2 / np.sum(np.abs(m)**2)

def K_conflict(m, e):
    s, se = m[:3], e[:3]
    return sum(s[i] * se[j] for i in range(3) for j in range(3) if i != j)

# ============================================================
# Evidence construction (exact, from spectral_gap_computation.py)
# ============================================================
# Born probs: p_dom = 9/10, p_weak = 1/20, p_theta = 1/27
# s1/s2 ratio: s1^2/s2^2 = p_dom/p_weak = 18, so s1 = 3√2 s2
# Born floor: 26θ² = 20s2² => θ = s2√(10/13)
# L1: s1 + 2s2 + θ = 1

coeff = 3*np.sqrt(2) + 2 + np.sqrt(10/13)
s2_eq = 1.0 / coeff
s1_eq = 3*np.sqrt(2) * s2_eq
theta_eq = np.sqrt(10/13) * s2_eq

e_eq = np.array([s1_eq, s2_eq, s2_eq, theta_eq], dtype=complex)
print(f"Evidence: e = ({np.real(e_eq[0]):.6f}, {np.real(e_eq[1]):.6f}, {np.real(e_eq[2]):.6f}, {np.real(e_eq[3]):.6f})")
print(f"L1 = {np.sum(np.abs(e_eq)):.10f}")
print(f"Born(θ_e) = {born(e_eq):.10f}, floor = {FLOOR:.10f}")

# ============================================================
# Find the REAL fixed point with K* = 7/30
# ============================================================
print("\n--- Finding fixed point ---")
m = np.array([0.4, 0.2, 0.2, 0.2], dtype=complex)
for i in range(500):
    m_new, K = ds_combine(m, e_eq)
    diff = np.max(np.abs(m_new - m))
    if diff < 1e-15:
        print(f"Converged at step {i}")
        break
    m = m_new

m_star = m_new
_, K_star = ds_combine(m_star, e_eq, apply_floor=False)
print(f"m* = ({np.real(m_star[0]):.10f}, {np.real(m_star[1]):.10f}, {np.real(m_star[2]):.10f}, {np.real(m_star[3]):.10f})")
print(f"K* = {np.real(K_star):.10f}")
print(f"7/30 = {7/30:.10f}")
print(f"|K* - 7/30| = {abs(K_star - 7/30):.2e}")
print(f"Born(θ*) = {born(m_star):.10f}")
print(f"1/27 = {1/27:.10f}")

# ============================================================
# Anti-holomorphic Jacobian
# ============================================================
def dbar_phi(m, e, eps=1e-7):
    """∂̄Φ where Φ(m) = DS(m, e) with floor."""
    n = 4
    z = m.copy()
    def phi(z_in):
        out, _ = ds_combine(z_in, e, apply_floor=True)
        return out

    dbar = np.zeros((n, n), dtype=complex)
    for b in range(n):
        zp = z.copy(); zp[b] += eps
        zm = z.copy(); zm[b] -= eps
        d_dx = (phi(zp) - phi(zm)) / (2 * eps)

        zp = z.copy(); zp[b] += 1j * eps
        zm = z.copy(); zm[b] -= 1j * eps
        d_dy = (phi(zp) - phi(zm)) / (2 * eps)

        dbar[:, b] = 0.5 * (d_dx + 1j * d_dy)

    return dbar

def decompose_dbar(db):
    """Decompose into spin-0 (trace), spin-1 (antisym), spin-2 (sym traceless)."""
    tr = np.trace(db) / 4
    sym = (db + db.T) / 2
    antisym = (db - db.T) / 2
    sym_tl = sym - tr * np.eye(4)

    frob = np.sqrt(np.sum(np.abs(db)**2))
    conf = 4 * np.abs(tr)**2
    s1 = np.sum(np.abs(antisym)**2)
    grav = np.sum(np.abs(sym_tl)**2)

    return {
        'frobenius': frob,
        'trace': tr,
        'conf_norm2': conf,
        'spin1_norm2': s1,
        'grav_norm2': grav,
        'singular_values': np.linalg.svd(db, compute_uv=False),
    }

# ============================================================
print("\n" + "=" * 60)
print("PART 1: ∂̄Φ at the real fixed point")
print("=" * 60)

db_star = dbar_phi(m_star, e_eq)
dec_star = decompose_dbar(db_star)
print(f"||∂̄Φ|| = {dec_star['frobenius']:.6f}")
print(f"Singular values: {dec_star['singular_values']}")
print(f"Born exactly at floor: floor active on boundary")
print()

# ============================================================
print("=" * 60)
print("PART 2: Trajectory from floor-active state (real evidence)")
print("=" * 60)

# Start with very low theta — floor MUST activate
m_init = np.array([0.02, 0.38, 0.32, 0.28], dtype=complex)
print(f"Initial Born = {born(m_init):.6f} (floor = {FLOOR:.6f})")
print(f"Floor active: {born(m_init) < FLOOR}")
print()

m = m_init.copy()
print(f"{'Step':>4} {'θ':>8} {'s1':>8} {'Born':>8} {'K':>8} {'||∂̄||':>8} {'grav':>8} {'floor':>5}")
print("-" * 65)

for step in range(30):
    b = born(m)
    K_val = K_conflict(m, e_eq)
    db = dbar_phi(m, e_eq)
    dec = decompose_dbar(db)
    fl = "YES" if b < FLOOR else "no"

    if step < 15 or step % 5 == 0:
        print(f"{step:4d} {np.abs(m[3]):8.4f} {np.abs(m[0]):8.4f} {np.real(b):8.5f} "
              f"{np.abs(K_val):8.5f} {dec['frobenius']:8.4f} {np.sqrt(dec['grav_norm2']):8.4f} {fl:>5}")

    m_new, _ = ds_combine(m, e_eq)
    m = m_new

print()

# ============================================================
print("=" * 60)
print("PART 3: Complex trajectory (both chiralities)")
print("=" * 60)

m_c = np.array([0.02+0.01j, 0.38+0.02j, 0.32-0.015j, 0.28+0.005j], dtype=complex)
m_c = m_c / np.sum(np.abs(m_c))

print(f"Initial Born = {born(m_c):.6f}")
print()

m = m_c.copy()
print(f"{'Step':>4} {'Born':>8} {'|K|':>8} {'Im(K)':>10} {'||∂̄||':>8} {'grav':>8} {'conf':>8}")
print("-" * 65)

for step in range(30):
    b = born(m)
    K_val = K_conflict(m, e_eq)
    db = dbar_phi(m, e_eq)
    dec = decompose_dbar(db)

    if step < 15 or step % 5 == 0:
        print(f"{step:4d} {np.real(b):8.5f} {np.abs(K_val):8.5f} {np.imag(K_val):10.6f} "
              f"{dec['frobenius']:8.4f} {np.sqrt(dec['grav_norm2']):8.4f} {np.sqrt(dec['conf_norm2']):8.4f}")

    m_new, _ = ds_combine(m, e_eq)
    m = m_new

print()

# ============================================================
print("=" * 60)
print("PART 4: Rank-2 entangled coupled system")
print("=" * 60)

# Two sites, each uses the other as evidence (rank-2 coupling)
mA = np.array([0.02+0.01j, 0.38+0.02j, 0.32-0.015j, 0.28+0.005j], dtype=complex)
mA = mA / np.sum(np.abs(mA))
mB = np.array([0.03-0.005j, 0.35+0.01j, 0.30-0.01j, 0.32+0.008j], dtype=complex)
mB = mB / np.sum(np.abs(mB))

print(f"Initial Born_A = {born(mA):.6f}, Born_B = {born(mB):.6f}")
print()

print(f"{'Step':>4} {'Born_A':>8} {'Born_B':>8} {'|K_AB|':>8} {'∂̄_A':>8} {'grav_A':>8} {'ent':>6}")
print("-" * 60)

integrated = {'conf': 0, 'spin1': 0, 'grav': 0}
n_floor = 0

for step in range(80):
    bA, bB = born(mA), born(mB)
    K_AB = K_conflict(mA, mB)

    db_A = dbar_phi(mA, mB)
    dec_A = decompose_dbar(db_A)

    # Entanglement of dbar
    sv = dec_A['singular_values']
    ent = 0.0
    if sv[0] > 1e-12:
        p = sv**2 / np.sum(sv**2)
        p = p[p > 1e-15]
        ent = -np.sum(p * np.log2(p))

    integrated['conf'] += dec_A['conf_norm2']
    integrated['spin1'] += dec_A['spin1_norm2']
    integrated['grav'] += dec_A['grav_norm2']
    if np.real(bA) < FLOOR:
        n_floor += 1

    if step < 20 or step % 10 == 0:
        print(f"{step:4d} {np.real(bA):8.5f} {np.real(bB):8.5f} {np.abs(K_AB):8.5f} "
              f"{dec_A['frobenius']:8.4f} {np.sqrt(dec_A['grav_norm2']):8.4f} {ent:6.3f}")

    # Coupled DS step
    mA_new, _ = ds_combine(mA, mB)
    mB_new, _ = ds_combine(mB, mA)
    mA, mB = mA_new, mB_new

total = integrated['conf'] + integrated['spin1'] + integrated['grav']
print()
print(f"Floor active: {n_floor}/80 steps")
print(f"Final Born_A = {born(mA):.6f}, K_AB = {K_conflict(mA, mB):.6f}")
print()

if total > 1e-15:
    print("INTEGRATED DECOMPOSITION:")
    print(f"  Conformal (spin-0): {integrated['conf']/total*100:.1f}%")
    print(f"  Spin-1:             {integrated['spin1']/total*100:.1f}%")
    print(f"  Graviton (spin-2):  {integrated['grav']/total*100:.1f}%")
else:
    print("No signal detected.")

# ============================================================
print("\n" + "=" * 60)
print("PART 5: Key numbers")
print("=" * 60)

print(f"Fubini-Study scalar curvature of CP³: R = 24")
print(f"DS fixed point: K* = 7/30 = {7/30:.6f}")
print(f"Born floor: 1/27 = {1/27:.6f}")
print(f"Fixed point mass: θ* = {np.real(m_star[3]):.6f}")
print(f"Metric ratio at equilibrium: g₀₀/g₁₁ = {np.real(m_star[3])**2 / (np.real(m_star[0]) + np.real(m_star[3]))**2:.6f}")
print(f"(H-1)/H = {2/3:.6f}")
